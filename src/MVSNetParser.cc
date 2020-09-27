//
// Created by lucius on 9/24/20.
//

#include <boost/log/trivial.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <opencv2/imgcodecs.hpp>
#include <fstream>
#include "MVSNetParser.hpp"
#include "Options.hpp"

void MVSNetParser::ReadPair(const std::string &path) {
  std::ifstream in(path);
  if (!in.is_open()) {
    BOOST_LOG_TRIVIAL(error) << "can not open output file";
    exit(EXIT_FAILURE);
  }
  uint32_t total_view_num = 0;
  in >> total_view_num;

  imagePairs_.resize(total_view_num);
  uint32_t ref_id;
  uint32_t src_num;
  uint32_t src_id;
  double src_score;
  for (auto i = 0; i < total_view_num; i++) {
    in >> ref_id >> src_num;
    if (ref_id != i) {
      BOOST_LOG_TRIVIAL(error) << "unmatch index in inpute file expect" << i << " got " << ref_id;
      exit(EXIT_FAILURE);
    }
    imagePairs_[ref_id].reserve(src_num);
    for (auto j = 0; j < src_num; j++) {
      in >> src_id >> src_score;
      imagePairs_[ref_id].emplace_back(src_id, src_score);
    }
  }
  std::string dummy;
  in >> dummy;
  if ((dummy != "") or !in.eof()) {
    BOOST_LOG_TRIVIAL(error) << "file size is too big";
    exit(EXIT_FAILURE);
  }
  in.close();
}

void MVSNetParser::WritePair(const std::string &path) {
  std::ofstream out(path, std::ios::trunc);
  if (!out.is_open()) {
    BOOST_LOG_TRIVIAL(error) << "can not open output file";
    exit(EXIT_FAILURE);
  }
  out.precision(5);
  out << imagePairs_.size() << std::endl;
  for (int i = 0; i < imagePairs_.size(); i++) {
    out << i << std::endl;
    out << imagePairs_[i].size();
    for (int j = 0; j < imagePairs_[i].size(); j++) {
      out << " " << imagePairs_[i][j].first << " " << imagePairs_[i][j].second;
    }
    out << std::endl;
  }
  out.close();
}

void MVSNetParser::WriteCams(const std::string &cam_dir, const ColmapSparseInfo &csi, int max_d) {
  boost::filesystem::create_directory(cam_dir);
  for (int index = 0; index < csi.index2imageid.size(); index++) {
    auto i = csi.index2imageid[index];
    const auto &image = csi.images_.at(i);
    std::ofstream out(cam_dir + boost::str(boost::format("/%08d_cam.txt") % index), std::ios::trunc);
    out.precision(17);
    out << "extrinsic" << std::endl;
    for (int j = 0; j < 4; j++) {
      out << image.extr(j, 0) << " " << image.extr(j, 1)
          << " " << image.extr(j, 2) << " " << image.extr(j, 3) << std::endl;
    }
    out << std::endl << "intrinsic" << std::endl;
    for (int j = 0; j < 3; j++) {
      out << image.intr(j, 0) << " " << image.intr(j, 1)
          << " " << image.intr(j, 2) << std::endl;
    }
    out << std::endl;
    double internal = (image.depth_max - image.depth_min) / (max_d - 1);
    internal =  internal/ options.interval_scale;
    out << image.depth_min << " " << internal << " " << max_d
        << " " << image.depth_max << std::endl;
  }
}

void MVSNetParser::WriteImages(const std::string &out_image_dir, const std::string &in_image_dir,
                               const ColmapSparseInfo &csi) {
  boost::filesystem::create_directory(out_image_dir);
  for (int index = 0; index < csi.index2imageid.size(); index++) {
    auto i = csi.index2imageid[index];
    const auto &image = csi.images_.at(i);
    auto in_img_path = boost::filesystem::path(in_image_dir + "/" + image.name);
    auto out_img_path = out_image_dir + "/" + str(boost::format("%08d.jpg") % index);
    if (boost::filesystem::extension(in_img_path) == ".jpg") {
      boost::filesystem::rename(in_img_path, out_img_path);
    } else {
      auto img = cv::imread(in_img_path.string(), cv::IMREAD_UNCHANGED);
      cv::imwrite(out_img_path, img);
    }
  }
}
