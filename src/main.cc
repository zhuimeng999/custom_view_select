#include <iostream>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/program_options.hpp>
#include <opencv2/core.hpp>
#include <unordered_map>
#include "ColmapSparseInfo.hpp"
#include "MVSNetParser.hpp"

namespace po = boost::program_options;

double angle_sigma1 = 1;
double angle_sigma2 = 10;
double angle_theta = 5;
int num_view = 10;

std::string in_image_dir;
std::string sparse_dir;
std::string output_dir;

void parse_commandline(int argc, char *argv[])
{
  try
  {
    po::options_description desc{"Options"};
    desc.add_options()("help,h", "Help screen")
        ("in_image_dir", po::value<std::string>(), "in_image_dir")
        ("sparse_dir", po::value<std::string>(), "sparse_dir")
        ("output_dir", po::value<std::string>(), "output_dir")
        ("sigma1", po::value<double>(), "sigma1")
        ("sigma2", po::value<double>(), "sigma2")
        ("theta", po::value<double>(), "theta")
        ("num_view", po::value<int>(), "num_view");

    po::positional_options_description pos_desc;
    pos_desc.add("in_image_dir", 1).add("sparse_dir", 1).add("output_dir", 1);

    po::command_line_parser parser{argc, argv};
//    parser.options(desc).positional(pos_desc).allow_unregistered();
    parser.options(desc).positional(pos_desc);
    po::parsed_options parsed_options = parser.run();

    po::variables_map vm;
    store(parsed_options, vm);
    notify(vm);

    if (vm.count("help")){
      std::cout << desc << '\n';
      exit(EXIT_FAILURE);
    }

    if (vm.count("in_image_dir") != 1){
      BOOST_LOG_TRIVIAL(error) << "you must provide output_dir";
      exit(EXIT_FAILURE);
    } else {
      in_image_dir = vm["in_image_dir"].as<std::string>();
    }
    if (vm.count("sparse_dir") != 1){
      BOOST_LOG_TRIVIAL(error) << "you must provide sparse_dir";
      exit(EXIT_FAILURE);
    } else {
      sparse_dir = vm["sparse_dir"].as<std::string>();
    }
    if (vm.count("output_dir") != 1){
      BOOST_LOG_TRIVIAL(error) << "you must provide output_dir";
      exit(EXIT_FAILURE);
    } else {
      output_dir = vm["output_dir"].as<std::string>();
    }

    if (vm.count("sigma1")){
      angle_sigma1 = vm["sigma1"].as<double>();
    }
    if (vm.count("sigma2")){
      angle_sigma2 = vm["sigma2"].as<double>();
    }
    if (vm.count("theta")){
      angle_theta = vm["theta"].as<double>();
    }
    if (vm.count("num_view")){
      num_view = vm["num_view"].as<double>();
    }
  }
  catch (const po::error &ex)
  {
    std::cerr << ex.what() << '\n';
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char *argv[]) {
  boost::log::core::get()->set_filter(
          boost::log::trivial::severity >= boost::log::trivial::info
  );

  parse_commandline(argc, argv);

  ColmapSparseInfo csi;
  BOOST_LOG_TRIVIAL(info) << "reading sparse reconstruction result ...";
  csi.Read(sparse_dir);
  auto num_image = csi.images_.size();

  std::vector<ColmapSparseInfo::image_t> index2imageid;
  std::unordered_map<ColmapSparseInfo::image_t, uint64_t> imageid2index;
  std::vector<ImageInfoAdditional> imageInfoAdditionals;
  index2imageid.reserve(num_image);
  imageid2index.reserve(num_image);
  imageInfoAdditionals.resize(num_image);

  BOOST_LOG_TRIVIAL(info) << "building continue index ...";
  std::vector<const ColmapSparseInfo::Image *> image_vecotor;
  image_vecotor.reserve(num_image);
  for (const auto &it: csi.images_) {
    image_vecotor.push_back(&it.second);
  }
  std::sort(image_vecotor.begin(), image_vecotor.end(),
            [](const ColmapSparseInfo::Image *x, const ColmapSparseInfo::Image *y){ return x->image_id < y->image_id;});
  for(const auto it: image_vecotor){
    imageid2index.emplace(it->image_id, index2imageid.size());
    index2imageid.push_back(it->image_id);
  }

  BOOST_LOG_TRIVIAL(info) << "extract image information ...";
#pragma omp parallel for default(none) shared(index2imageid, csi, imageInfoAdditionals, num_image)
  for (auto i = 0; i < num_image; i++) {
    const auto &image = csi.images_[index2imageid[i]];
    Eigen::Matrix4d extr = Eigen::Matrix4d::Identity();
    Eigen::Matrix3d intr = Eigen::Matrix3d::Identity();

    const auto &R = Eigen::Quaterniond(image.Qvec.x(), image.Qvec.y(), image.Qvec.z(), image.Qvec.w()).toRotationMatrix();
    extr.block(0, 0, 3, 3) = R;
    extr.block(0, 3, 3, 1) = image.Tvec;

    std::vector<double> depths;
    for (auto it: image.point3D_ids) {
      if (it != ColmapSparseInfo::kInvalidPoint3DId) {
        Eigen::Vector3d &&proj = R * csi.points3D_[it].XYZ + image.Tvec;
        depths.push_back(proj.z());
      }
    }
    std::sort(depths.begin(), depths.end());
    auto &additional = imageInfoAdditionals[i];
    additional.name = image.name;
    additional.extr = extr;
    additional.intr = csi.cameras_[image.camera_id].GetK();
    additional.centor = -R.transpose() * image.Tvec;
    additional.direction = R.col(2);
    additional.depth_min = depths[static_cast<int>(depths.size() * 0.01)];
    additional.depth_max = depths[static_cast<int>(depths.size() * 0.99)];
  }

  BOOST_LOG_TRIVIAL(info) << "prepare score matrix and vectorlize point map ...";
  std::vector<std::vector<double>> score_matrix(num_image);
  for (int i = 0; i < num_image; i++) {
    score_matrix[i].resize(num_image, 0);
  }

  std::vector<const ColmapSparseInfo::Point3D *> points3D_vector;
  points3D_vector.reserve(csi.points3D_.size());
  for (const auto &it: csi.points3D_) {
    points3D_vector.push_back(&it.second);
  }

  BOOST_LOG_TRIVIAL(info) << "building score matrix ...";
#pragma omp parallel for default(shared)
  for (uint64_t i = 0; i < points3D_vector.size(); i++) {
    const auto &tracks = points3D_vector[i]->track;
    const auto &pos = points3D_vector[i]->XYZ;
    for (int j = 0; j < tracks.size(); j++) {
      auto image_index1 = imageid2index.at(tracks[j].first);
      const auto &c1 = imageInfoAdditionals[image_index1].centor;
      for (int k = j + 1; k < tracks.size(); k++) {
        auto image_index2 = imageid2index.at(tracks[k].first);
        const auto &c2 = imageInfoAdditionals[image_index2].centor;

        if (image_index1 == image_index2) {
          auto image = csi.images_[tracks[k].first];

          if ((image.point3D_ids.at(tracks[j].second) != points3D_vector[i]->point3D_id) or
              (image.point3D_ids.at(tracks[j].second) != image.point3D_ids.at(tracks[k].second))) {
            BOOST_LOG_TRIVIAL(error) << "consistant check failed";
            exit(EXIT_FAILURE);
          }
          if (tracks[j].second != tracks[k].second) {
            const auto &&displace = image.points2D[tracks[j].second] - image.points2D[tracks[k].second];
            BOOST_LOG_TRIVIAL(debug) << "point have two project point in same image: displace in image "
                                     << displace.transpose() << " meam reprojection error "
                                     << points3D_vector[i]->error;
          }

          continue;
        }

        const auto &&v1 = c1 - pos;
        const auto &&v2 = c2 - pos;
        double angle = (180.0f / M_PIf64) * acos(v1.dot(v2) / (v1.norm() * v2.norm()));
        double sigma = (angle < angle_theta) ? angle_sigma1 : angle_sigma2;
        double kernel = (angle - angle_theta) / sigma;
        double score = exp(-kernel * kernel / 2);
#pragma omp critical
        {
          score = score_matrix[image_index1][image_index2] + score;
          score_matrix[image_index1][image_index2] = score;
          score_matrix[image_index2][image_index1] = score;
        }
      }
    }
  }

  std::vector<std::vector<int>> score_sorted(num_image);

#pragma omp parallel for default(none) shared(score_sorted, score_matrix, num_image)
  for (auto i = 0; i < num_image; i++) {
    score_sorted[i].resize(num_image);
    const auto &scores = score_matrix[i];
    for (int j = 0; j < num_image; j++) {
      score_sorted[i][j] = j;
    }
    std::sort(score_sorted[i].begin(), score_sorted[i].end(),
              [&scores](int x, int y) { return scores[x] > scores[y]; });
  }

  BOOST_LOG_TRIVIAL(info) << "write output stream";
  MVSNetParser pairWriter;
  pairWriter.imagePairs_.resize(num_image);
  for(int i = 0; i < num_image; i++){
    for(int j = 0; j < num_view; j++){
      pairWriter.imagePairs_[i].emplace_back(score_sorted[i][j], score_matrix[i][score_sorted[i][j]]);
    }
  }
  pairWriter.WritePair(output_dir + "/pair.txt");
  pairWriter.WriteCams(output_dir + "/cams", imageInfoAdditionals, 32);
  pairWriter.WriteImages(output_dir + "/images", in_image_dir, imageInfoAdditionals);
  BOOST_LOG_TRIVIAL(info) << "done ...";
//  MVSNetParser pairReader;
//  pairReader.ReadPair("/home/lucius/data/workspace/pair.txt");
//  for(int i = 0; i < pairWriter.imagePairs_.size(); i++){
//    for(int j = 0; j < num_view; j++){
//      if(pairWriter.imagePairs_[i][j].first != pairReader.imagePairs_[i][j].first){
//        BOOST_LOG_TRIVIAL(error) << "mismatch " << i << " vs" << j;
//        BOOST_LOG_TRIVIAL(error) << "    mismatch " << pairWriter.imagePairs_[i][j].first << " vs" << pairReader.imagePairs_[i][j].first;
//        BOOST_LOG_TRIVIAL(error) << "    mismatch " << pairWriter.imagePairs_[i][j].second << " vs" << pairReader.imagePairs_[i][j].second;
//      }
//    }
//  }
  return 0;
}
