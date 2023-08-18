#pragma once

#include "igl/decimate_callback_types.h"
#include <Eigen/Core>

void qslim_aligned_callbacks(
    Eigen::MatrixXi &E,
    std::vector<std::tuple<Eigen::MatrixXd, Eigen::RowVectorXd, double>>
        &quadrics,
    int &v1, int &v2, Eigen::MatrixXd &PD1, Eigen::MatrixXd &PD2,
    Eigen::VectorXd &PV1, Eigen::VectorXd &PV2,
    igl::decimate_cost_and_placement_callback &cost_and_placement,
    igl::decimate_pre_collapse_callback &pre_collapse,
    igl::decimate_post_collapse_callback &post_collapse);
