//
// Created by lucius on 9/24/20.
//

#ifndef CUSTOM_VIEW_SELECT_MVSNETPARSER_HPP
#define CUSTOM_VIEW_SELECT_MVSNETPARSER_HPP

#include <Eigen/Eigen>
#include <vector>
#include "ColmapSparseInfo.hpp"

struct ImageInfoAdditional {
  Eigen::Matrix4d extr;
  Eigen::Matrix3d intr;
  Eigen::Vector3d centor;
  Eigen::Vector3d direction;
  double depth_min;
  double depth_max;
  std::string name;
};

class MVSNetParser {
public:
  void ReadPair(const std::string &path);
  void WritePair(const std::string &path);

  void WriteCams(const std::string &cam_dir, const std::vector<ImageInfoAdditional> &imageInfoAdditionals, int max_d = 128);

  void WriteImages(const std::string &out_image_dir, const std::string &in_image_dir, const std::vector<ImageInfoAdditional> &imageInfoAdditionals);

  std::vector<std::vector<std::pair<uint32_t, double>>> imagePairs_;
};


#endif //CUSTOM_VIEW_SELECT_MVSNETPARSER_HPP
