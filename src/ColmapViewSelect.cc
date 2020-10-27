//
// Created by lucius on 9/26/20.
//
#include <boost/log/trivial.hpp>
#include "ColmapViewSelect.hpp"
#include "Options.hpp"


template<typename T>
inline T Percentile(const std::vector<T> &elems, const double p) {
  if (elems.empty()) {
    BOOST_LOG_TRIVIAL(error) << "elems size must greater than zero";
    exit(EXIT_FAILURE);
  }

  const int idx = static_cast<int>(std::round(p / 100 * (elems.size() - 1)));
  const size_t percentile_idx =
      std::max(0, std::min(static_cast<int>(elems.size() - 1), idx));

  std::vector<T> ordered_elems = elems;
  std::nth_element(ordered_elems.begin(),
                   ordered_elems.begin() + percentile_idx, ordered_elems.end());

  return ordered_elems.at(percentile_idx);
}

inline double CalculateTriangulationAngle(const Eigen::Vector3d &proj_center1,
                                          const Eigen::Vector3d &proj_center2,
                                          const Eigen::Vector3d &point3D) {
  const double baseline_length_squared =
      (proj_center1 - proj_center2).squaredNorm();

  const double ray_length_squared1 = (point3D - proj_center1).squaredNorm();
  const double ray_length_squared2 = (point3D - proj_center2).squaredNorm();

  // Using "law of cosines" to compute the enclosing angle between rays.
  const double denominator =
      2.0 * std::sqrt(ray_length_squared1 * ray_length_squared2);
  if (denominator == 0.0) {
    return 0.0;
  }
  const double nominator =
      ray_length_squared1 + ray_length_squared2 - baseline_length_squared;
  const double angle = std::abs(std::acos(nominator / denominator));

  // Triangulation is unstable for acute angles (far away points) and
  // obtuse angles (close points), so always compute the minimum angle
  // between the two intersecting rays.
  return std::min(angle, M_PI - angle);
}

void colmap_view_select(const ColmapSparseInfo &csi, std::vector<std::vector<double>> &score_matrix) {
  const auto num_image = csi.index2imageid.size();

  // Use maximum number of overlapping images as source images. Overlapping
  // will be sorted based on the number of shared points to the reference
  // image and the top ranked images are selected. Note that images are only
  // selected if some points have a sufficient triangulation angle.

  std::vector<const ColmapSparseInfo::Point3D *> points3D_vector;
  points3D_vector.reserve(csi.points3D_.size());
  for (const auto &it: csi.points3D_) {
    points3D_vector.push_back(&it.second);
  }

  std::vector<std::vector<std::pair<uint64_t, std::vector<double>>>> shared_points(csi.index2imageid.size());
  for (int i = 0; i < num_image; i++) {
    shared_points[i].resize(num_image, {0, {}});
  }
#pragma omp parallel for default(none) shared(points3D_vector, csi, options, shared_points, boost::log::keywords::severity)
  for (uint64_t i = 0; i < points3D_vector.size(); i++) {
    const auto &tracks = points3D_vector[i]->track;
    const auto &pos = points3D_vector[i]->XYZ;
    for (int j = 0; j < tracks.size(); j++) {
      const auto image_index1 = csi.imageid2index.at(tracks[j].first);
      const auto &c1 = csi.images_.at(tracks[j].first).centor;
      for (int k = 0; k < j; k++) {
        const auto image_index2 = csi.imageid2index.at(tracks[k].first);
        const auto &c2 = csi.images_.at(tracks[k].first).centor;

        if (image_index1 == image_index2) {
          auto image = csi.images_.at(tracks[k].first);

          if ((image.point3D_ids.at(tracks[j].second) != points3D_vector[i]->point3D_id) or
              (image.point3D_ids.at(tracks[j].second) != image.point3D_ids.at(tracks[k].second))) {
            BOOST_LOG_TRIVIAL(error) << "consistant check failed";
            exit(EXIT_FAILURE);
          }
          if (tracks[j].second != tracks[k].second) {
            const auto &&displace = image.points2D[tracks[j].second] - image.points2D[tracks[k].second];
            BOOST_LOG_TRIVIAL(debug) << "point have two project point in same image: displace in image "
                                     << displace.transpose() << " meam reprojection error "
                                     << points3D_vector[i]->error;
          }

          continue;
        }
        const auto angle = (180.0 / M_PIf64) * CalculateTriangulationAngle(c1, c2, pos);

#pragma omp critical
        {
          if (image_index1 > image_index2) {
            shared_points[image_index1][image_index2].first += 1;
            shared_points[image_index1][image_index2].second.push_back(angle);
          } else {
            shared_points[image_index2][image_index1].first += 1;
            shared_points[image_index2][image_index1].second.push_back(angle);
          }
        }
      }
    }
  }

#pragma omp parallel for default(none) schedule(dynamic) shared(score_matrix, shared_points, options)
  for (auto i = 0; i < shared_points.size(); i++) {
    for (auto k = 0; k < i; k++) {
      const auto &pair_info = shared_points[i][k];
      if (pair_info.second.size() > 0) {
        const auto angle = Percentile(pair_info.second, options.kTriangulationAnglePercentile);
        if (angle >= options.min_triangulation_angle) {
          score_matrix[i][k] = pair_info.first;
          score_matrix[k][i] = pair_info.first;
        }
      }
    }
  }
}
