//
// Created by lucius on 9/26/20.
//

#include <boost/log/trivial.hpp>
#include "MVSNetViewSelect.hpp"

#include "Options.hpp"


void mvsnet_view_select(const ColmapSparseInfo &csi, std::vector<std::vector<double>> &score_matrix) {
  std::vector<const ColmapSparseInfo::Point3D *> points3D_vector;
  points3D_vector.reserve(csi.points3D_.size());
  for (const auto &it: csi.points3D_) {
    points3D_vector.push_back(&it.second);
  }

#pragma omp parallel for default(none) shared(points3D_vector, csi, options, score_matrix, boost::log::keywords::severity)
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
                                     << displace.transpose() << " mean reprojection error "
                                     << points3D_vector[i]->error;
          }

          continue;
        }

        const auto &&v1 = c1 - pos;
        const auto &&v2 = c2 - pos;
        double angle = (180.0f / M_PIf64) * acos(v1.dot(v2) / std::sqrt(v1.squaredNorm() * v2.squaredNorm()));
        double sigma = (angle < options.angle_theta) ? options.angle_sigma1 : options.angle_sigma2;
        double kernel = (angle - options.angle_theta) / sigma;
        double score = exp(-kernel * kernel / 2);
#pragma omp critical
        {
          score = score_matrix[image_index1][image_index2] + score;
          score_matrix[image_index1][image_index2] = score;
          score_matrix[image_index2][image_index1] = score;
        }
//#pragma omp atomic
//        score_matrix[image_index1][image_index2] += score;
//#pragma omp atomic
//        score_matrix[image_index2][image_index1] += score;
      }
    }
  }
}
