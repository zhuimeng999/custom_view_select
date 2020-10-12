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
