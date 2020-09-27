//
// Created by lucius on 9/26/20.
//

#ifndef CUSTOM_VIEW_SELECT_OPTIONS_HPP
#define CUSTOM_VIEW_SELECT_OPTIONS_HPP

#include <string>

struct Options {
  std::string in_image_dir;
  std::string sparse_dir;
  std::string output_dir;

  bool selection_only = false;
  std::string alg = "mvsnet";
  /*common options */
  uint64_t num_view = 10;

  /*mvsnet */
  double angle_sigma1 = 1;
  double angle_sigma2 = 10;
  double angle_theta = 5;
  double interval_scale = 1;
  uint32_t max_d = 128;

  /*colmap*/
  double kTriangulationAnglePercentile = 75;
  double min_triangulation_angle = 1.0f;   // Minimum triangulation angle in degrees.
};

extern Options options;

void parse_commandline(int argc, char *argv[]);

#endif //CUSTOM_VIEW_SELECT_OPTIONS_HPP
