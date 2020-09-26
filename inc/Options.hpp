//
// Created by lucius on 9/26/20.
//

#ifndef CUSTOM_VIEW_SELECT_OPTIONS_HPP
#define CUSTOM_VIEW_SELECT_OPTIONS_HPP

#include <string>

struct Options {
  double angle_sigma1 = 1;
  double angle_sigma2 = 10;
  double angle_theta = 5;
  int num_view = 10;

  std::string in_image_dir;
  std::string sparse_dir;
  std::string output_dir;
};

extern Options options;
void parse_commandline(int argc, char *argv[]);

#endif //CUSTOM_VIEW_SELECT_OPTIONS_HPP
