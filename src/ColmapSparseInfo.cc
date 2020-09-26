//
// Created by lucius on 9/24/20.
//
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/log/trivial.hpp>
#include <boost/filesystem.hpp>
#include "ColmapSparseInfo.hpp"

namespace fs = boost::filesystem;


template<typename T>
T ReverseBytes(const T &data) {
  T data_reversed = data;
  std::reverse(reinterpret_cast<char *>(&data_reversed),
               reinterpret_cast<char *>(&data_reversed) + sizeof(T));
  return data_reversed;
}

inline bool IsLittleEndian() {
#ifdef BOOST_BIG_ENDIAN
  return false;
#else
  return true;
#endif
}

inline bool IsBigEndian() {
#ifdef BOOST_BIG_ENDIAN
  return true;
#else
  return false;
#endif
}

template<typename T>
T LittleEndianToNative(const T x) {
  if (IsLittleEndian()) {
    return x;
  } else {
    return ReverseBytes(x);
  }
}

template<typename T>
T BigEndianToNative(const T x) {
  if (IsBigEndian()) {
    return x;
  } else {
    return ReverseBytes(x);
  }
}

template<typename T>
T NativeToLittleEndian(const T x) {
  if (IsLittleEndian()) {
    return x;
  } else {
    return ReverseBytes(x);
  }
}

template<typename T>
T NativeToBigEndian(const T x) {
  if (IsBigEndian()) {
    return x;
  } else {
    return ReverseBytes(x);
  }
}

template<typename T>
T ReadBinaryLittleEndian(std::istream *stream) {
  T data_little_endian;
  stream->read(reinterpret_cast<char *>(&data_little_endian), sizeof(T));
  return LittleEndianToNative(data_little_endian);
}

template<typename T>
void ReadBinaryLittleEndian(std::istream *stream, std::vector<T> *data) {
  for (size_t i = 0; i < data->size(); ++i) {
    (*data)[i] = ReadBinaryLittleEndian<T>(stream);
  }
}

template<typename T>
void WriteBinaryLittleEndian(std::ostream *stream, const T &data) {
  const T data_little_endian = NativeToLittleEndian(data);
  stream->write(reinterpret_cast<const char *>(&data_little_endian), sizeof(T));
}

template<typename T>
void WriteBinaryLittleEndian(std::ostream *stream, const std::vector<T> &data) {
  for (const auto &elem : data) {
    WriteBinaryLittleEndian<T>(stream, elem);
  }
}

ColmapSparseInfo::ColmapSparseInfo() {
  union {
    uint32_t i;
    char c[4];
  } bint = {0x01020304};

  big = (bint.c[0] == 1);
  if (big) {
    BOOST_LOG_TRIVIAL(fatal) << "we only support little endian";
    exit(EXIT_FAILURE);
  }
}

void ColmapSparseInfo::Read(const std::string &path) {
  if (fs::is_regular_file(path + "/cameras.bin") &&
      fs::is_regular_file(path + "/images.bin") &&
      fs::is_regular_file(path + "/points3D.bin")) {
    ReadBinary(path);
  } else if (fs::is_regular_file(path + "/cameras.txt") &&
             fs::is_regular_file(path + "/images.txt") &&
             fs::is_regular_file(path + "/points3D.txt")) {
    ReadText(path);
  } else {
    BOOST_LOG_TRIVIAL(error) << "cameras, images, points3D files do not exist at " << path;
    exit(EXIT_FAILURE);
  }
}

const std::vector<std::pair<std::string, int>> CAMERA_INFOS = {
    {"SIMPLE_PINHOLE", 3},  // 0
    {"PINHOLE",        4},         // 1
    {"SIMPLE_RADIAL",  4},   // 2
    {"RADIAL",         5},          // 3
    {"OPENCV",         8},          // 4
    {"OPENCV_FISHEYE", 8},  // 5
    {"FULL_OPENCV",    12}     //6
};

int ColmapSparseInfo::GetModelId(const std::string &model_name) {
  for (int i = 0; i < CAMERA_INFOS.size(); i++) {
    if (CAMERA_INFOS[i].first == model_name) {
      return i;
    }
  }
  BOOST_LOG_TRIVIAL(error) << "unkown camera model name " << model_name;
  exit(EXIT_FAILURE);
}

std::string ColmapSparseInfo::GetModelName(const int model_id) {
  return CAMERA_INFOS.at(model_id).first;
}

Eigen::Matrix3d ColmapSparseInfo::Camera::GetK() const {
  Eigen::Matrix3d ret = Eigen::Matrix3d::Identity();
  switch (model_id) {
    case 0:
      ret(0, 0) = ret(1, 1) = params[0];
      ret(0, 2) = params[1];
      ret(1, 2) = params[2];
      break;
    case 1:
      ret(0, 0) = params[0];
      ret(1, 1) = params[1];
      ret(0, 2) = params[2];
      ret(1, 2) = params[3];
      break;
    case 2:
      ret(0, 0) = ret(1, 1) = params[0];
      ret(0, 2) = params[1];
      ret(1, 2) = params[2];
//      if(params[3] != 0){
//        BOOST_LOG_TRIVIAL(error) << "camera  SIMPLE_RADIAL got k = " << params[3] << " not 0";
//        exit(EXIT_FAILURE);
//      }
      break;
    case 3:
      ret(0, 0) = ret(1, 1) = params[0];
      ret(0, 2) = params[1];
      ret(1, 2) = params[2];
//      if((params[3] != 0) || params[4] != 0){
//        BOOST_LOG_TRIVIAL(error) << "camera  RADIAL got k1 = " << params[3] << " k2 = " << params[4] << " not 0";
//        exit(EXIT_FAILURE);
//      }
      break;
    default:
      BOOST_LOG_TRIVIAL(error) << "unsupport camera model " << model_id;
      exit(EXIT_FAILURE);
  }
  return ret;
}

const ColmapSparseInfo::point3D_t ColmapSparseInfo::kInvalidPoint3DId = std::numeric_limits<ColmapSparseInfo::point3D_t>::max();

void ColmapSparseInfo::BuildContinueIndex() {
  index2imageid.reserve(images_.size());
  imageid2index.reserve(images_.size());

  std::vector<const ColmapSparseInfo::Image *> image_vecotor;
  image_vecotor.reserve(images_.size());
  for (const auto &it: images_) {
    image_vecotor.push_back(&it.second);
  }
  std::sort(image_vecotor.begin(), image_vecotor.end(),
            [](const ColmapSparseInfo::Image *x, const ColmapSparseInfo::Image *y) {
                return x->image_id < y->image_id;
            });
  for (const auto it: image_vecotor) {
    imageid2index.emplace(it->image_id, index2imageid.size());
    index2imageid.push_back(it->image_id);
  }
}

void ColmapSparseInfo::ComputeExtraInfo() {
#pragma omp parallel for default(none)
  for (auto i = 0; i < index2imageid.size(); i++) {
    auto &image = images_[index2imageid[i]];
    Eigen::Matrix4d extr = Eigen::Matrix4d::Identity();
    Eigen::Matrix3d intr = Eigen::Matrix3d::Identity();

    const auto &R = Eigen::Quaterniond(image.Qvec.x(), image.Qvec.y(), image.Qvec.z(),
                                       image.Qvec.w()).toRotationMatrix();
    extr.block(0, 0, 3, 3) = R;
    extr.block(0, 3, 3, 1) = image.Tvec;

    std::vector<double> depths;
    for (auto it: image.point3D_ids) {
      if (it != ColmapSparseInfo::kInvalidPoint3DId) {
        depths.push_back(R.row(2).dot(points3D_[it].XYZ) + image.Tvec.z());
      }
    }
    std::sort(depths.begin(), depths.end());
    image.extr = extr;
    image.intr = cameras_[image.camera_id].GetK();
    image.centor = -R.transpose() * image.Tvec;
    image.direction = R.col(2);
    image.depth_min = depths[static_cast<int>(depths.size() * 0.01)];
    image.depth_max = depths[static_cast<int>(depths.size() * 0.99)];
  }
}

void ColmapSparseInfo::ReadText(const std::string &path) {
  ReadCamerasText(path + "/cameras.txt");
  ReadImagesText(path + "/images.txt");
  ReadPoints3DText(path + "/points3D.txt");
}

void ColmapSparseInfo::ReadBinary(const std::string &path) {
  ReadCamerasBinary(path + "/cameras.bin");
  ReadImagesBinary(path + "/images.bin");
  ReadPoints3DBinary(path + "/points3D.bin");
}

void ColmapSparseInfo::WriteText(const std::string &path) const {
  WriteCamerasText(path + "/cameras.txt");
  WriteImagesText(path + "/images.txt");
  WritePoints3DText(path + "/points3D.txt");
}

void ColmapSparseInfo::WriteBinary(const std::string &path) const {
  WriteCamerasBinary(path + "/cameras.bin");
  WriteImagesBinary(path + "/images.bin");
  WritePoints3DBinary(path + "/points3D.bin");
}

void ColmapSparseInfo::ReadCamerasText(const std::string &path) {
  cameras_.clear();

  std::ifstream file(path);
  if (!file.is_open()) {
    BOOST_LOG_TRIVIAL(error) << "can not open file" << path;
    exit(EXIT_FAILURE);
  }

  std::string line;
  std::string item;

  while (std::getline(file, line)) {
    boost::trim(line);

    if (line.empty() || line[0] == '#') {
      continue;
    }

    std::stringstream line_stream(line);

    Camera camera;

    // ID
    std::getline(line_stream, item, ' ');
    camera.camera_id = std::stoul(item);

    // MODEL
    std::getline(line_stream, item, ' ');
    camera.model_id = GetModelId(item);

    // WIDTH
    std::getline(line_stream, item, ' ');
    camera.width = std::stoll(item);

    // HEIGHT
    std::getline(line_stream, item, ' ');
    camera.height = std::stoll(item);

    // PARAMS
    camera.params.clear();
    while (!line_stream.eof()) {
      std::getline(line_stream, item, ' ');
      camera.params.push_back(std::stold(item));
    }

    if (CAMERA_INFOS.at(camera.model_id).second != camera.params.size()) {
      BOOST_LOG_TRIVIAL(fatal) << "camera model params mismatch: expect " << CAMERA_INFOS.at(camera.model_id).second
                               << " Got " << camera.params.size();
      exit(EXIT_FAILURE);
    }

    cameras_.emplace(camera.camera_id, camera);
  }
}

void ColmapSparseInfo::ReadImagesText(const std::string &path) {
  images_.clear();

  std::ifstream file(path);
  if (!file.is_open()) {
    BOOST_LOG_TRIVIAL(error) << "can not open file" << path;
    exit(EXIT_FAILURE);
  }

  std::string line;
  std::string item;

  while (std::getline(file, line)) {
    boost::trim(line);

    if (line.empty() || line[0] == '#') {
      continue;
    }

    std::stringstream line_stream1(line);

    // ID
    std::getline(line_stream1, item, ' ');
    const image_t image_id = std::stoul(item);

    class Image image;
    image.image_id = image_id;

//    image.SetRegistered(true);
//    reg_image_ids_.push_back(image_id);

    // QVEC (qw, qx, qy, qz)
    std::getline(line_stream1, item, ' ');
    image.Qvec(0) = std::stold(item);

    std::getline(line_stream1, item, ' ');
    image.Qvec(1) = std::stold(item);

    std::getline(line_stream1, item, ' ');
    image.Qvec(2) = std::stold(item);

    std::getline(line_stream1, item, ' ');
    image.Qvec(3) = std::stold(item);

    image.NormalizeQvec();


    // TVEC
    std::getline(line_stream1, item, ' ');
    image.Tvec(0) = std::stold(item);

    std::getline(line_stream1, item, ' ');
    image.Tvec(1) = std::stold(item);

    std::getline(line_stream1, item, ' ');
    image.Tvec(2) = std::stold(item);

    // CAMERA_ID
    std::getline(line_stream1, item, ' ');
    image.camera_id = std::stoul(item);

    // NAME
    std::getline(line_stream1, item, ' ');
    image.name = item;

    // POINTS2D
    if (!std::getline(file, line)) {
      break;
    }

    boost::trim(line);
    std::stringstream line_stream2(line);

//    std::vector<Eigen::Vector2d> points2D;
//    std::vector<point3D_t> point3D_ids;

    if (!line.empty()) {
      while (!line_stream2.eof()) {
        Eigen::Vector2d point;

        std::getline(line_stream2, item, ' ');
        point.x() = std::stold(item);

        std::getline(line_stream2, item, ' ');
        point.y() = std::stold(item);

        image.points2D.push_back(point);

        std::getline(line_stream2, item, ' ');
        if (item == "-1") {
          image.point3D_ids.push_back(kInvalidPoint3DId);
        } else {
          image.point3D_ids.push_back(std::stoll(item));
        }
      }
    }

//    image.SetUp(Camera(image.CameraId()));
//    image.SetPoints2D(points2D);

//    for (point2D_t point2D_idx = 0; point2D_idx < image.NumPoints2D();
//         ++point2D_idx) {
//      if (point3D_ids[point2D_idx] != kInvalidPoint3DId) {
//        image.SetPoint3DForPoint2D(point2D_idx, point3D_ids[point2D_idx]);
//      }
//    }

    images_.emplace(image.image_id, image);
  }
}

void ColmapSparseInfo::ReadPoints3DText(const std::string &path) {
  points3D_.clear();

  std::ifstream file(path);
  if (!file.is_open()) {
    BOOST_LOG_TRIVIAL(error) << "can not open file" << path;
    exit(EXIT_FAILURE);
  }

  std::string line;
  std::string item;

  while (std::getline(file, line)) {
    boost::trim(line);

    if (line.empty() || line[0] == '#') {
      continue;
    }

    std::stringstream line_stream(line);

    // ID
    std::getline(line_stream, item, ' ');
    const point3D_t point3D_id = std::stoll(item);

    // Make sure, that we can add new 3D points after reading 3D points
    // without overwriting existing 3D points.
//    num_added_points3D_ = std::max(num_added_points3D_, point3D_id);

    class Point3D point3D;
    point3D.point3D_id = point3D_id;
    // XYZ
    std::getline(line_stream, item, ' ');
    point3D.XYZ(0) = std::stold(item);

    std::getline(line_stream, item, ' ');
    point3D.XYZ(1) = std::stold(item);

    std::getline(line_stream, item, ' ');
    point3D.XYZ(2) = std::stold(item);

    // Color
    std::getline(line_stream, item, ' ');
    point3D.Color(0) = static_cast<uint8_t>(std::stoi(item));

    std::getline(line_stream, item, ' ');
    point3D.Color(1) = static_cast<uint8_t>(std::stoi(item));

    std::getline(line_stream, item, ' ');
    point3D.Color(2) = static_cast<uint8_t>(std::stoi(item));

    // ERROR
    std::getline(line_stream, item, ' ');
    point3D.error = std::stod(item);

    // TRACK
    while (!line_stream.eof()) {
      std::getline(line_stream, item, ' ');
      boost::trim(item);
      if (item.empty()) {
        break;
      }
      image_t image_id = std::stoul(item);

      std::getline(line_stream, item, ' ');
      uint32_t point2D_idx = std::stoul(item);

      point3D.track.emplace_back(image_id, point2D_idx);
    }

    point3D.track.shrink_to_fit();

    points3D_.emplace(point3D_id, point3D);
  }
}

void ColmapSparseInfo::ReadCamerasBinary(const std::string &path) {
  std::ifstream file(path, std::ios::binary);
  if (!file.is_open()) {
    BOOST_LOG_TRIVIAL(error) << "can not open file" << path;
    exit(EXIT_FAILURE);
  }

  const size_t num_cameras = ReadBinaryLittleEndian<uint64_t>(&file);
  for (size_t i = 0; i < num_cameras; ++i) {
    class Camera camera;
    camera.camera_id = ReadBinaryLittleEndian<camera_t>(&file);
    camera.model_id = ReadBinaryLittleEndian<int>(&file);
    camera.width = ReadBinaryLittleEndian<uint64_t>(&file);
    camera.height = ReadBinaryLittleEndian<uint64_t>(&file);
    camera.params.resize(CAMERA_INFOS[camera.model_id].second);
    ReadBinaryLittleEndian<double>(&file, &camera.params);
    if (CAMERA_INFOS.at(camera.model_id).second != camera.params.size()) {
      BOOST_LOG_TRIVIAL(fatal) << "camera model params mismatch: expect " << CAMERA_INFOS.at(camera.model_id).second
                               << " Got " << camera.params.size();
      exit(EXIT_FAILURE);
    }
    cameras_.emplace(camera.camera_id, camera);
  }
  if (file.peek() != EOF) {
    auto current = file.tellg();
    file.seekg(0, file.end);
    auto size = file.tellg();
    BOOST_LOG_TRIVIAL(error) << "unexpect file size " << size << ", only used " << current;
    exit(EXIT_FAILURE);
  }
}

void ColmapSparseInfo::ReadImagesBinary(const std::string &path) {
  std::ifstream file(path, std::ios::binary);
  if (!file.is_open()) {
    BOOST_LOG_TRIVIAL(error) << "can not open file" << path;
    exit(EXIT_FAILURE);
  }

  const size_t num_reg_images = ReadBinaryLittleEndian<uint64_t>(&file);
  for (size_t i = 0; i < num_reg_images; ++i) {
    class Image image;

    image.image_id = ReadBinaryLittleEndian<image_t>(&file);

    image.Qvec(0) = ReadBinaryLittleEndian<double>(&file);
    image.Qvec(1) = ReadBinaryLittleEndian<double>(&file);
    image.Qvec(2) = ReadBinaryLittleEndian<double>(&file);
    image.Qvec(3) = ReadBinaryLittleEndian<double>(&file);
    image.NormalizeQvec();

    image.Tvec(0) = ReadBinaryLittleEndian<double>(&file);
    image.Tvec(1) = ReadBinaryLittleEndian<double>(&file);
    image.Tvec(2) = ReadBinaryLittleEndian<double>(&file);

    image.camera_id = ReadBinaryLittleEndian<camera_t>(&file);

    char name_char;
    do {
      file.read(&name_char, 1);
      if (name_char != '\0') {
        image.name += name_char;
      }
    } while (name_char != '\0');

    const size_t num_points2D = ReadBinaryLittleEndian<uint64_t>(&file);

//    std::vector<Eigen::Vector2d> points2D;
    image.points2D.reserve(num_points2D);
//    std::vector<point3D_t> point3D_ids;
    image.point3D_ids.reserve(num_points2D);
    for (size_t j = 0; j < num_points2D; ++j) {
      const double x = ReadBinaryLittleEndian<double>(&file);
      const double y = ReadBinaryLittleEndian<double>(&file);
      image.points2D.emplace_back(x, y);
      image.point3D_ids.push_back(ReadBinaryLittleEndian<point3D_t>(&file));
    }

//    image.SetUp(Camera(image.CameraId()));
//    image.SetPoints2D(points2D);
//
//    for (point2D_t point2D_idx = 0; point2D_idx < image.NumPoints2D();
//         ++point2D_idx) {
//      if (point3D_ids[point2D_idx] != kInvalidPoint3DId) {
//        image.SetPoint3DForPoint2D(point2D_idx, point3D_ids[point2D_idx]);
//      }
//    }
//
//    image.SetRegistered(true);
//    reg_image_ids_.push_back(image.ImageId());

    images_.emplace(image.image_id, image);
  }
}

void ColmapSparseInfo::ReadPoints3DBinary(const std::string &path) {
  std::ifstream file(path, std::ios::binary);
  if (!file.is_open()) {
    BOOST_LOG_TRIVIAL(error) << "can not open file" << path;
    exit(EXIT_FAILURE);
  }

  const size_t num_points3D = ReadBinaryLittleEndian<uint64_t>(&file);
  for (size_t i = 0; i < num_points3D; ++i) {
    class Point3D point3D;

    const point3D_t point3D_id = ReadBinaryLittleEndian<point3D_t>(&file);
    point3D.point3D_id = point3D_id;
//    num_added_points3D_ = std::max(num_added_points3D_, point3D_id);

    point3D.XYZ(0) = ReadBinaryLittleEndian<double>(&file);
    point3D.XYZ(1) = ReadBinaryLittleEndian<double>(&file);
    point3D.XYZ(2) = ReadBinaryLittleEndian<double>(&file);
    point3D.Color(0) = ReadBinaryLittleEndian<uint8_t>(&file);
    point3D.Color(1) = ReadBinaryLittleEndian<uint8_t>(&file);
    point3D.Color(2) = ReadBinaryLittleEndian<uint8_t>(&file);
    point3D.error = ReadBinaryLittleEndian<double>(&file);

    const size_t track_length = ReadBinaryLittleEndian<uint64_t>(&file);
    for (size_t j = 0; j < track_length; ++j) {
      const image_t image_id = ReadBinaryLittleEndian<image_t>(&file);
      const point2D_t point2D_idx = ReadBinaryLittleEndian<point2D_t>(&file);
      point3D.track.emplace_back(image_id, point2D_idx);
    }
    point3D.track.shrink_to_fit();

    points3D_.emplace(point3D_id, point3D);
  }
}

void ColmapSparseInfo::WriteCamerasText(const std::string &path) const {
  std::ofstream file(path, std::ios::trunc);
  if (!file.is_open()) {
    BOOST_LOG_TRIVIAL(error) << "can not open file" << path;
    exit(EXIT_FAILURE);
  }

  // Ensure that we don't loose any precision by storing in text.
  file.precision(17);

  file << "# Camera list with one line of data per camera:" << std::endl;
  file << "#   CAMERA_ID, MODEL, WIDTH, HEIGHT, PARAMS[]" << std::endl;
  file << "# Number of cameras: " << cameras_.size() << std::endl;

  for (const auto &camera : cameras_) {
    std::ostringstream line;

    line << camera.first << " ";
    line << GetModelName(camera.second.camera_id) << " ";
    line << camera.second.width << " ";
    line << camera.second.height << " ";

    for (const double param : camera.second.params) {
      line << param << " ";
    }

    std::string line_string = line.str();
    line_string = line_string.substr(0, line_string.size() - 1);

    file << line_string << std::endl;
  }
}

void ColmapSparseInfo::WriteImagesText(const std::string &path) const {
  std::ofstream file(path, std::ios::trunc);
  if (!file.is_open()) {
    BOOST_LOG_TRIVIAL(error) << "can not open file" << path;
    exit(EXIT_FAILURE);
  }

  // Ensure that we don't loose any precision by storing in text.
  file.precision(17);

  file << "# Image list with two lines of data per image:" << std::endl;
  file << "#   IMAGE_ID, QW, QX, QY, QZ, TX, TY, TZ, CAMERA_ID, "
          "NAME"
       << std::endl;
  file << "#   POINTS2D[] as (X, Y, POINT3D_ID)" << std::endl;
  file << "# Number of images: " << images_.size()
       << ", mean observations per image: "
       << 0 << std::endl;

  for (const auto &image : images_) {
//    if (!image.second.IsRegistered()) {
//      continue;
//    }

    std::ostringstream line;
    std::string line_string;

    line << image.first << " ";

    // QVEC (qw, qx, qy, qz)
    const Eigen::Vector4d normalized_qvec = image.second.Qvec / image.second.Qvec.norm();
    line << normalized_qvec(0) << " ";
    line << normalized_qvec(1) << " ";
    line << normalized_qvec(2) << " ";
    line << normalized_qvec(3) << " ";

    // TVEC
    line << image.second.Tvec(0) << " ";
    line << image.second.Tvec(1) << " ";
    line << image.second.Tvec(2) << " ";

    line << image.second.camera_id << " ";

    line << image.second.name;

    file << line.str() << std::endl;

    line.str("");
    line.clear();

    for (auto i = 0; i < image.second.points2D.size(); i++) {
      const Eigen::Vector2d &point2D = image.second.points2D[i];
      const point3D_t point3D_id = image.second.point3D_ids[i];
      line << point2D.x() << " ";
      line << point2D.y() << " ";
      if (point3D_id != kInvalidPoint3DId) {
        line << point3D_id << " ";
      } else {
        line << -1 << " ";
      }
    }
    line_string = line.str();
    line_string = line_string.substr(0, line_string.size() - 1);
    file << line_string << std::endl;
  }
}

void ColmapSparseInfo::WritePoints3DText(const std::string &path) const {
  std::ofstream file(path, std::ios::trunc);
  if (!file.is_open()) {
    BOOST_LOG_TRIVIAL(error) << "can not open file" << path;
    exit(EXIT_FAILURE);
  }

  // Ensure that we don't loose any precision by storing in text.
  file.precision(17);

  file << "# 3D point list with one line of data per point:" << std::endl;
  file << "#   POINT3D_ID, X, Y, Z, R, G, B, ERROR, "
          "TRACK[] as (IMAGE_ID, POINT2D_IDX)"
       << std::endl;
  file << "# Number of points: " << points3D_.size()
       << ", mean track length: " << 0 << std::endl;

  for (const auto &point3D : points3D_) {
    file << point3D.first << " ";
    file << point3D.second.XYZ(0) << " ";
    file << point3D.second.XYZ(1) << " ";
    file << point3D.second.XYZ(2) << " ";
    file << static_cast<int>(point3D.second.Color(0)) << " ";
    file << static_cast<int>(point3D.second.Color(1)) << " ";
    file << static_cast<int>(point3D.second.Color(2)) << " ";
    file << point3D.second.error << " ";

    std::ostringstream line;

    for (const auto &track_el : point3D.second.track) {
      line << track_el.first << " ";
      line << track_el.second << " ";
    }

    std::string line_string = line.str();
    line_string = line_string.substr(0, line_string.size() - 1);

    file << line_string << std::endl;
  }
}

void ColmapSparseInfo::WriteCamerasBinary(const std::string &path) const {
  std::ofstream file(path, std::ios::trunc | std::ios::binary);
  if (!file.is_open()) {
    BOOST_LOG_TRIVIAL(error) << "can not open file" << path;
    exit(EXIT_FAILURE);
  }

  WriteBinaryLittleEndian<uint64_t>(&file, cameras_.size());

  for (const auto &camera : cameras_) {
    WriteBinaryLittleEndian<camera_t>(&file, camera.first);
    WriteBinaryLittleEndian<int>(&file, camera.second.model_id);
    WriteBinaryLittleEndian<uint64_t>(&file, camera.second.width);
    WriteBinaryLittleEndian<uint64_t>(&file, camera.second.height);
    for (const double param : camera.second.params) {
      WriteBinaryLittleEndian<double>(&file, param);
    }
  }
}

void ColmapSparseInfo::WriteImagesBinary(const std::string &path) const {
  std::ofstream file(path, std::ios::trunc | std::ios::binary);
  if (!file.is_open()) {
    BOOST_LOG_TRIVIAL(error) << "can not open file" << path;
    exit(EXIT_FAILURE);
  }

  WriteBinaryLittleEndian<uint64_t>(&file, images_.size());

  for (const auto &image : images_) {
//    if (!image.second.IsRegistered()) {
//      continue;
//    }

    WriteBinaryLittleEndian<image_t>(&file, image.first);

    const Eigen::Vector4d normalized_qvec = image.second.Qvec.normalized();
//            NormalizeQuaternion(image.second.Qvec());
    WriteBinaryLittleEndian<double>(&file, normalized_qvec(0));
    WriteBinaryLittleEndian<double>(&file, normalized_qvec(1));
    WriteBinaryLittleEndian<double>(&file, normalized_qvec(2));
    WriteBinaryLittleEndian<double>(&file, normalized_qvec(3));

    WriteBinaryLittleEndian<double>(&file, image.second.Tvec(0));
    WriteBinaryLittleEndian<double>(&file, image.second.Tvec(1));
    WriteBinaryLittleEndian<double>(&file, image.second.Tvec(2));

    WriteBinaryLittleEndian<camera_t>(&file, image.second.camera_id);

    const std::string name = image.second.name + '\0';
    file.write(name.c_str(), name.size());

    WriteBinaryLittleEndian<uint64_t>(&file, image.second.points2D.size());
    for (auto i = 0; i < image.second.points2D.size(); i++) {
      const Eigen::Vector2d &point2D = image.second.points2D[i];
      WriteBinaryLittleEndian<double>(&file, point2D.x());
      WriteBinaryLittleEndian<double>(&file, point2D.y());
      WriteBinaryLittleEndian<point3D_t>(&file, image.second.point3D_ids[i]);
    }
  }
}

void ColmapSparseInfo::WritePoints3DBinary(const std::string &path) const {
  std::ofstream file(path, std::ios::trunc | std::ios::binary);
  if (!file.is_open()) {
    BOOST_LOG_TRIVIAL(error) << "can not open file" << path;
    exit(EXIT_FAILURE);
  }

  WriteBinaryLittleEndian<uint64_t>(&file, points3D_.size());

  for (const auto &point3D : points3D_) {
    WriteBinaryLittleEndian<point3D_t>(&file, point3D.first);
    WriteBinaryLittleEndian<double>(&file, point3D.second.XYZ(0));
    WriteBinaryLittleEndian<double>(&file, point3D.second.XYZ(1));
    WriteBinaryLittleEndian<double>(&file, point3D.second.XYZ(2));
    WriteBinaryLittleEndian<uint8_t>(&file, point3D.second.Color(0));
    WriteBinaryLittleEndian<uint8_t>(&file, point3D.second.Color(1));
    WriteBinaryLittleEndian<uint8_t>(&file, point3D.second.Color(2));
    WriteBinaryLittleEndian<double>(&file, point3D.second.error);

    WriteBinaryLittleEndian<uint64_t>(&file, point3D.second.track.size());
    for (const auto &track_el : point3D.second.track) {
      WriteBinaryLittleEndian<image_t>(&file, track_el.first);
      WriteBinaryLittleEndian<point2D_t>(&file, track_el.second);
    }
  }
}
