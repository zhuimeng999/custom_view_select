//
// Created by lucius on 9/26/20.
//

#ifndef CUSTOM_VIEW_SELECT_MVSNET_VIEW_SELECT_H
#define CUSTOM_VIEW_SELECT_MVSNET_VIEW_SELECT_H

#include "ColmapSparseInfo.hpp"
#include "MVSNetParser.hpp"

void mvsnet_view_select(const ColmapSparseInfo &csi, std::vector<std::vector<double>> &score_matrix);

#endif //CUSTOM_VIEW_SELECT_MVSNET_VIEW_SELECT_H
