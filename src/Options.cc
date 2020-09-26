//
// Created by lucius on 9/26/20.
//

#include <boost/program_options.hpp>
#include <boost/log/trivial.hpp>
#include "Options.hpp"
namespace po = boost::program_options;

Options options;

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
      BOOST_LOG_TRIVIAL(error) << desc << '\n';
      exit(EXIT_FAILURE);
    }

    if (vm.count("in_image_dir") != 1){
      BOOST_LOG_TRIVIAL(error) << "you must provide output_dir";
      exit(EXIT_FAILURE);
    } else {
      options.in_image_dir = vm["in_image_dir"].as<std::string>();
    }
    if (vm.count("sparse_dir") != 1){
      BOOST_LOG_TRIVIAL(error) << "you must provide sparse_dir";
      exit(EXIT_FAILURE);
    } else {
      options.sparse_dir = vm["sparse_dir"].as<std::string>();
    }
    if (vm.count("output_dir") != 1){
      BOOST_LOG_TRIVIAL(error) << "you must provide output_dir";
      exit(EXIT_FAILURE);
    } else {
      options.output_dir = vm["output_dir"].as<std::string>();
    }

    if (vm.count("sigma1")){
      options.angle_sigma1 = vm["sigma1"].as<double>();
    }
    if (vm.count("sigma2")){
      options.angle_sigma2 = vm["sigma2"].as<double>();
    }
    if (vm.count("theta")){
      options.angle_theta = vm["theta"].as<double>();
    }
    if (vm.count("num_view")){
      options.num_view = vm["num_view"].as<double>();
    }
  }
  catch (const po::error &ex)
  {
    BOOST_LOG_TRIVIAL(error) << ex.what() << '\n';
    exit(EXIT_FAILURE);
  }
}
