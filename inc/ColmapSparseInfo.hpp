//
// Created by lucius on 9/24/20.
//

#ifndef CUSTOM_VIEW_SELECT_COLMAPSPARSEINFO_HPP
#define CUSTOM_VIEW_SELECT_COLMAPSPARSEINFO_HPP

#include <vector>
#include <unordered_map>
#include <Eigen/Eigen>

class ColmapSparseInfo {
public:

  typedef uint32_t camera_t;
  typedef uint32_t image_t;
  typedef uint32_t point2D_t;
  typedef uint64_t point3D_t;
  static const point3D_t kInvalidPoint3DId;

  struct Camera {
    camera_t camera_id;
    int model_id;
    uint64_t width;
    uint64_t height;
    std::vector<double> params;

    Eigen::Matrix3d GetK() const;
  };

  struct Image {
    image_t image_id;
    camera_t camera_id;
    Eigen::Vector4d Qvec;
    Eigen::Vector3d Tvec;
    std::string name;
    std::vector<Eigen::Vector2d> points2D;
    std::vector<point3D_t> point3D_ids;

    /* extra info */
    Eigen::Matrix4d extr;
    Eigen::Matrix3d R;
    Eigen::Matrix3d intr;
    Eigen::Vector3d centor;
    Eigen::Vector3d direction;
    double depth_min;
    double depth_max;

    void NormalizeQvec() {
      auto qn = Qvec.norm();
      if (qn == 0.0) {
        Qvec = Eigen::Vector4d(1.0, 0, 0, 0);
      } else {
        Qvec = Qvec / Qvec.norm();
      }
    }
  };

  struct Point3D {
    Eigen::Vector3d XYZ;
    Eigen::Matrix<uint8_t, 3, 1> Color;
    double error;
    std::vector<std::pair<image_t, point2D_t>> track;
    point3D_t point3D_id;
  };

  ColmapSparseInfo();

  void BuildContinueIndex();

  void ComputeExtraInfo();

  void Read(const std::string &path);

  void ReadText(const std::string &path);

  void ReadBinary(const std::string &path);

  static int GetModelId(const std::string &model_name);

  static std::string GetModelName(const int model_id);

  std::unordered_map<camera_t, Camera> cameras_;
  std::unordered_map<image_t, Image> images_;
  std::unordered_map<point3D_t, Point3D> points3D_;

  /*image index mapper*/
  std::vector<image_t> index2imageid;
  std::unordered_map<image_t, uint32_t> imageid2index;
private:

  void WriteText(const std::string &path) const;

  void WriteBinary(const std::string &path) const;

  void ReadCamerasText(const std::string &path);

  void ReadImagesText(const std::string &path);

  void ReadPoints3DText(const std::string &path);

  void ReadCamerasBinary(const std::string &path);

  void ReadImagesBinary(const std::string &path);

  void ReadPoints3DBinary(const std::string &path);

  void WriteCamerasText(const std::string &path) const;

  void WriteImagesText(const std::string &path) const;

  void WritePoints3DText(const std::string &path) const;

  void WriteCamerasBinary(const std::string &path) const;

  void WriteImagesBinary(const std::string &path) const;

  void WritePoints3DBinary(const std::string &path) const;

  bool big;
};


#endif //CUSTOM_VIEW_SELECT_COLMAPSPARSEINFO_HPP
