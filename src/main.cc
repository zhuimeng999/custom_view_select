#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <fstream>

#include "ColmapSparseInfo.hpp"
#include "MVSNetViewSelect.hpp"
#include "ColmapViewSelect.hpp"
#include "Options.hpp"

int main(int argc, char *argv[]) {
  boost::log::core::get()->set_filter(
      boost::log::trivial::severity >= boost::log::trivial::info
  );

  parse_commandline(argc, argv);

  ColmapSparseInfo csi;
  BOOST_LOG_TRIVIAL(info) << "reading sparse reconstruction result ...";
  csi.Read(options.sparse_dir);
  BOOST_LOG_TRIVIAL(info) << "building continue index ...";
  csi.BuildContinueIndex();

  BOOST_LOG_TRIVIAL(info) << "extract image information ...";
  csi.ComputeExtraInfo();

  const auto num_image = csi.index2imageid.size();
  BOOST_LOG_TRIVIAL(info) << "building score matrix ...";
  std::vector<std::vector<double>> score_matrix(num_image);
  for (int i = 0; i < num_image; i++) {
    score_matrix[i].resize(num_image, 0);
  }

  if (options.alg == "mvsnet") {
    mvsnet_view_select(csi, score_matrix);
  } else {
    colmap_view_select(csi, score_matrix);
  }


  std::vector<std::vector<int>> score_sorted(num_image);

#pragma omp parallel for default(none) shared(score_sorted, score_matrix)
  for (auto i = 0; i < score_sorted.size(); i++) {
    score_sorted[i].resize(score_sorted.size());
    const auto &scores = score_matrix[i];
    for (int j = 0; j < score_sorted.size(); j++) {
      score_sorted[i][j] = j;
    }
    std::sort(score_sorted[i].begin(), score_sorted[i].end(),
              [&scores](int x, int y) { return scores[x] > scores[y]; });
  }

  for(auto i = 0; i < csi.index2imageid.size(); i++){
    auto image_id = csi.index2imageid[i];
    const auto &image = csi.images_.at(image_id);
    const auto &camera = csi.cameras_.at(image.camera_id);
    std::vector<std::vector<Eigen::Vector3d>> ref_pos;
    ref_pos.resize(camera.height);
    for(auto h = 0; h < camera.height; h++){
      const auto K_inv = image.intr.inverse();
      ref_pos[h].resize(camera.width);
      for(auto w = 0; w < camera.width; w++){
        ref_pos[h][w] = K_inv*Eigen::Vector3d(w + 0.5, h + 0.5, 1.0);
        ref_pos[h][w] = ref_pos[h][w].normalized();
      }
    }
    for(auto j = 0; j < options.num_view; j++){
      auto image_id_y = csi.index2imageid[score_matrix[i][j]];
      const auto &image_y = csi.images_.at(image_id_y);
      const auto &camera_y = csi.cameras_.at(image_y.camera_id);
      const auto relative_R = image_y.R*image.R.transpose();
      const auto relative_T = image_y.Tvec - relative_R*image.Tvec;
      const auto centor = relative_T.hnormalized();

      const auto src_K_inv = image_y.intr.inverse();
      const auto pos_min = (src_K_inv*Eigen::Vector3d(0, 0, 1.)).hnormalized();
      const auto pos_max = (src_K_inv*Eigen::Vector3d(camera_y.width, camera_y.height, 1.)).hnormalized();
      for(auto h = 0; h < camera.height; h++){
        for(auto w = 0; w < camera.width; w++) {
          const auto src_pos = relative_R*ref_pos[h][w];
          const Eigen::Vector2d src_proj_pos = src_pos.hnormalized();
          Eigen::Vector2d plane_direction = centor - src_proj_pos;
          const auto distance = plane_direction.norm();
          plane_direction = plane_direction/distance;
          const auto alpha = src_pos.z()/relative_T.z();
          {
            Eigen::Vector2d range_min = (pos_min - src_proj_pos).array()/plane_direction.array();
            Eigen::Vector2d range_max = (pos_max - src_proj_pos).array()/plane_direction.array();
            if(plane_direction.x() > 0){
              if(plane_direction.y() > 0){
                Eigen::Vector2d(range_min.maxCoeff(), range_max.minCoeff());
              } else {
                Eigen::Vector2d(std::max(range_min.maxCoeff(), range_max.maxCoeff()), std::numeric_limits<double>::infinity());
              }
            } else {
              if(plane_direction.y() > 0){
                Eigen::Vector2d(std::numeric_limits<double>::infinity(), std::min(range_min.minCoeff(), range_max.minCoeff()));
              } else {
                Eigen::Vector2d(range_max.maxCoeff(),range_min.minCoeff());
              }
            }
          }
        }
      }
    }
  }

  BOOST_LOG_TRIVIAL(info) << "write output stream";
  if(options.selection_only){
    std::ofstream out(options.output_dir + "/patch-match.cfg", std::ios::trunc);
    if (!out.is_open()) {
      BOOST_LOG_TRIVIAL(error) << "can not open output file";
      exit(EXIT_FAILURE);
    }
    for (int i = 0; i < score_sorted.size(); i++) {
      const auto &ref_img_name = csi.images_[csi.index2imageid[i]].name;
      out << ref_img_name << std::endl;
      if(score_matrix[i][score_sorted[i][0]] < 1e-8){
        BOOST_LOG_TRIVIAL(error) << "no proper source view selected for image " << ref_img_name;
        exit(EXIT_FAILURE);
      }
      out << csi.images_[csi.index2imageid[score_sorted[i][0]]].name;
      for(int j = 1; j < options.num_view; j++){
        if(score_matrix[i][score_sorted[i][j]] > 1e-8){
          out << ", " << csi.images_[csi.index2imageid[score_sorted[i][j]]].name;
        }
      }
      out << std::endl;
    }
    out.close();
  } else {
    MVSNetParser pairWriter;
    pairWriter.imagePairs_.resize(num_image);
    for (int i = 0; i < num_image; i++) {
      for (int j = 0; j < options.num_view; j++) {
        pairWriter.imagePairs_[i].emplace_back(score_sorted[i][j], score_matrix[i][score_sorted[i][j]]);
      }
    }
    pairWriter.WritePair(options.output_dir + "/pair.txt");
    pairWriter.WriteCams(options.output_dir + "/cams", csi, options.max_d);
    pairWriter.WriteSparseDepth(options.output_dir + "/sparse_depth", csi, options.max_d);
    pairWriter.WriteImages(options.output_dir + "/images", options.in_image_dir, csi);
  }

  BOOST_LOG_TRIVIAL(info) << "done ...";

  return 0;
}
