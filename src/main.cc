#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

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
  MVSNetParser pairWriter;
  pairWriter.imagePairs_.resize(num_image);
  for (int i = 0; i < num_image; i++) {
    for (int j = 0; j < options.num_view; j++) {
      pairWriter.imagePairs_[i].emplace_back(score_sorted[i][j], score_matrix[i][score_sorted[i][j]]);
    }
  }
  pairWriter.WritePair(options.output_dir + "/pair.txt");
  pairWriter.WriteCams(options.output_dir + "/cams", csi, 32);
  pairWriter.WriteImages(options.output_dir + "/images", options.in_image_dir, csi);
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
