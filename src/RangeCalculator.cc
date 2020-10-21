//
// Created by lucius on 10/20/20.
//

#include <boost/log/trivial.hpp>
#include <Eigen/Eigen>
#include "RangeCalculator.hpp"
#include "Options.hpp"

static inline void solve3dxy(const double pxy, const double pz, const double nxy, const double nz,
                           const double val_min, const double val_max, Eigen::Vector2d &range)
{
  auto denominator1 = nxy - nz*val_min;
  auto numerator1 = pz*val_min - pxy;
  const auto denominator2 = nxy - nz*val_max;
  const auto numerator2 = pz*val_max - pxy;

  if(denominator1 == 0.){
    assert(denominator2 != 0);
    assert(numerator1 != 0);
    auto val = numerator2/denominator2;
    /* check size of median value between non trival value*/
    if((numerator1 > 0) ^ (denominator2 > 0)){
      range.x() = -std::numeric_limits<double>::infinity();
      range.y() = val;
    } else {
      range.x() = val;
      range.y() = std::numeric_limits<double>::infinity();
    }
  } else if(denominator2 == 0.){
    assert(numerator2 != 0);
    auto val = numerator1/denominator1;
    /* check size of median value between non trival value*/
    if((numerator2 > 0) ^ (denominator1 > 0)){
      range.x() = -std::numeric_limits<double>::infinity();
      range.y() = val;
    } else {
      range.x() = val;
      range.y() = std::numeric_limits<double>::infinity();
    }
  } else if((denominator1 > 0.) ^ (denominator2 > 0.)){
    auto val = numerator1/denominator1;
    auto val2 = numerator2/denominator2;
    /* we only keep value > 0, -infinite is droped */
    if(val2 < val){
      range.x() = val;
    } else {
      range.x() = val2;
    }
    range.y() = std::numeric_limits<double>::infinity();
  } else {
    auto val = numerator1/denominator1;
    auto val2 = numerator2/denominator2;
    if(val2 < val){
      range.x() = val2;
      range.y() = val;
    } else {
      range.x() = val;
      range.y() = val2;
    }
  }
}

static inline void solve2d(const Eigen::Vector2d &projection, const Eigen::Vector2d &direction,
                      const double wmin, const double wmax, const double hmin, const double hmax, Eigen::Vector2d &range)
{
  if(direction.x() == 0.){
    BOOST_LOG_TRIVIAL(info) << " unhandled branch " << direction;
    exit(EXIT_FAILURE);
  } else if(direction.y() == 0.){
    BOOST_LOG_TRIVIAL(info) << " unhandled branch " << direction;
    exit(EXIT_FAILURE);
  }

  auto range_wmin = (wmin - projection.x())/direction.x();
  auto range_wmax = (wmax - projection.x())/direction.x();
  auto range_hmin = (hmin - projection.y())/direction.y();
  auto range_hmax = (hmax - projection.y())/direction.y();

  if(direction.x() > 0.){
    if(direction.y() > 0.){
      range.x() = std::max(range_wmin, range_hmin);
      range.y() = std::min(range_wmax, range_hmax);
    } else {
      range.x() = std::max(range_wmin, range_hmax);
      range.y() = std::min(range_wmax, range_hmin);
    }
  } else {
    if(direction.y() > 0.){
      range.x() = std::max(range_wmax, range_hmin);
      range.y() = std::min(range_wmin, range_hmax);
    } else {
      range.x() = std::max(range_wmax, range_hmax);
      range.y() = std::min(range_wmin, range_hmin);
    }
  }
}

void compute_range(const ColmapSparseInfo &csi, std::vector<std::vector<int>> &score_sorted){
  for(auto i = 0; i < csi.index2imageid.size(); i++){
    auto image_id = csi.index2imageid[i];
    const auto &image = csi.images_.at(image_id);
    const auto &camera = csi.cameras_.at(image.camera_id);
    std::vector<std::vector<Eigen::Vector3d>> ref_pos;
    ref_pos.resize(camera.height);
    for(auto h = 0; h < camera.height; h++){
      const auto K_inv = image.intr.inverse();
      ref_pos[h].resize(camera.width);
      for(auto w = 0; w < camera.width; w++){
        ref_pos[h][w] = K_inv*Eigen::Vector3d(w + 0.5, h + 0.5, 1.0);
        assert(ref_pos[h][w].z() == 1.);
        ref_pos[h][w] = ref_pos[h][w].normalized();
      }
    }
    for(auto j = 0; j < options.num_view; j++){
      auto image_id_y = csi.index2imageid[score_sorted[i][j]];
      const auto &image_y = csi.images_.at(image_id_y);
      const auto &camera_y = csi.cameras_.at(image_y.camera_id);
      const Eigen::Matrix3d relative_R = image_y.R*image.R.transpose();
      const Eigen::Vector3d relative_T = image_y.Tvec - relative_R*image.Tvec;
      /* here we do not check whether relative_T.z() == 0., we check later */
      const Eigen::Vector2d center = relative_T.hnormalized();

      const Eigen::Matrix3d src_K_inv = image_y.intr.inverse();
      const Eigen::Vector3d pos_min = src_K_inv*Eigen::Vector3d(0, 0, 1.);
      const Eigen::Vector3d pos_max = src_K_inv*Eigen::Vector3d(camera_y.width, camera_y.height, 1.);
      BOOST_LOG_TRIVIAL(info) << "processing " <<  image_y.name;

      std::vector<std::vector<Eigen::Vector2d>> solutions(camera.height);
      for(auto h = 0; h < camera.height; h++){
        solutions[h].resize(camera.width);
        for(auto w = 0; w < camera.width; w++){
          Eigen::Vector3d src_pos = relative_R*ref_pos[h][w];
          Eigen::Vector2d range1, range2;
          solve3dxy(relative_T.x(), relative_T.z(), src_pos.x(), src_pos.z(), pos_min.x(), pos_max.x(), range1);
          solve3dxy(relative_T.y(), relative_T.z(), src_pos.y(), src_pos.z(), pos_min.y(), pos_max.y(), range2);
          double min_val = std::max(range1.x(), range2.x());
          double max_val = std::min(range1.y(), range2.y());
          if(src_pos.z() > 0.){
            min_val = std::max(min_val, -relative_T.z()/src_pos.z());
          } else {
            max_val = std::min(max_val, -relative_T.z()/src_pos.z());
          }
          min_val = std::max(min_val, 0.);
          assert(src_pos.z() != 0);
          auto projection = src_pos.hnormalized();
          Eigen::Vector2d direction;
          if(relative_T.z() == 0.){
            direction.x() = relative_T.x();
            direction.y() = relative_T.y();
          } else if(relative_T.z() > 0.) {   /* projection direction along depth radial, determined by the sign of z */
            direction = center - projection;
          } else { /* projection direction along depth radial, determined by the sign of z */
            direction = projection - center;
          }
          auto distance = direction.norm();
          direction = direction/distance;

          Eigen::Vector2d range;
          solve2d(projection, direction, pos_min.x(), pos_max.x(), pos_min.y(), pos_max.y(), range);
          if(range.x() < 0.){
            range.x() = 0.;
          }
          if(relative_T.z() > 0){ /* we can not cross the centor of camera when the direction point to centor */
            range.y() = std::min(range.y(), distance);
          }

          if(min_val >= max_val){
            assert((range.x() >= range.y()) or (range.y() <= 0));
          } else {
            assert(range.x() <= range.y());
            auto disparity_min = distance/std::abs(1. + min_val*src_pos.z()/relative_T.z());
            auto disparity_max = distance/std::abs(1. + max_val*src_pos.z()/relative_T.z());
            if(disparity_min > disparity_max){
              std::swap(disparity_min, disparity_max);
            }
            auto err1 = disparity_min - range.x();
            auto err2 = disparity_max - range.y();
            if((abs(err1/distance) > 1e-8) or (abs(err2/distance) > 1e-8)){
              BOOST_LOG_TRIVIAL(info) << pos_min.transpose() << " # " << pos_max.transpose();
              BOOST_LOG_TRIVIAL(info) << relative_T.transpose() << " # " << src_pos.transpose();
              BOOST_LOG_TRIVIAL(info) << range.transpose() << " # " << disparity_min << " " << disparity_max << " # " << min_val << " " << max_val;
              exit(EXIT_FAILURE);
            }
          }
        }
      }
    }
  }
}
