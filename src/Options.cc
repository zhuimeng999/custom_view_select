//
// Created by lucius on 9/26/20.
//

#include <boost/program_options.hpp>
#include <boost/log/trivial.hpp>
#include "Options.hpp"

namespace po = boost::program_options;

Options options;

void parse_commandline(int argc, char *argv[]) {
  try {
    po::options_description desc{"Options"};
    desc.add_options()("help,h", "Help screen")
        ("in_image_dir", po::value<std::string>(&options.in_image_dir), "in_image_dir")
        ("sparse_dir", po::value<std::string>(&options.sparse_dir), "sparse_dir")
        ("output_dir", po::value<std::string>(&options.output_dir), "output_dir")
        ("alg", po::value<std::string>(&options.alg), "alg{mvsnet, colmap}")
        ("selection_only", po::bool_switch(&options.selection_only), "selection_only")
        ("num_view", po::value<uint64_t>(&options.num_view), "num_view")
        ("sigma1", po::value<double>(&options.angle_sigma1), "sigma1")
        ("sigma2", po::value<double>(&options.angle_sigma2), "sigma2")
        ("theta", po::value<double>(&options.angle_theta), "theta")
        ("max_d", po::value<uint32_t>(&options.max_d), "max_d")
        ("interval_scale", po::value<double>(&options.interval_scale), "interval_scale")
        ("kTriangulationAnglePercentile", po::value<double>(&options.kTriangulationAnglePercentile), "kTriangulationAnglePercentile")
        ("min_triangulation_angle", po::value<double>(&options.min_triangulation_angle), "min_triangulation_angle");

    po::positional_options_description pos_desc;
    pos_desc.add("in_image_dir", 1).add("sparse_dir", 1).add("output_dir", 1);

    po::command_line_parser parser{argc, argv};
//    parser.options(desc).positional(pos_desc).allow_unregistered();
    parser.options(desc).positional(pos_desc);
    po::parsed_options parsed_options = parser.run();

    po::variables_map vm;
    store(parsed_options, vm);
    notify(vm);

    if (vm.count("help")) {
      BOOST_LOG_TRIVIAL(error) << desc << '\n';
      exit(EXIT_FAILURE);
    }

    if (options.in_image_dir.empty() or options.sparse_dir.empty() or options.output_dir.empty()) {
      BOOST_LOG_TRIVIAL(error) << "you must provide [in_image_dir sparse_dir output_dir]";
      exit(EXIT_FAILURE);
    }

    if ((options.alg != "mvsnet") and (options.alg != "colmap")) {
      BOOST_LOG_TRIVIAL(error) << "only {mvsnet, colmap} are support, got " << options.alg;
      exit(EXIT_FAILURE);
    }
    if ((options.kTriangulationAnglePercentile < 0) or (options.kTriangulationAnglePercentile > 100)) {
      BOOST_LOG_TRIVIAL(error) << "kTriangulationAnglePercentile must between [0, 100], got "
                               << options.kTriangulationAnglePercentile;
      exit(EXIT_FAILURE);
    }
  }
  catch (const po::error &ex) {
    BOOST_LOG_TRIVIAL(error) << ex.what() << '\n';
    exit(EXIT_FAILURE);
  }
}
