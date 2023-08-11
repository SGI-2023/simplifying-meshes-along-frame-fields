#ifndef FRAME_FIELD_H
#define FRAME_FIELD_H

#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/principal_curvature.h>
#include <Eigen/Core>

#include "utils.h"

void get_frame_field(const Eigen::MatrixXd& V,
										 const Eigen::MatrixXi& F,
										 Eigen::MatrixXd& PD1,
										 Eigen::MatrixXd& PD2,
										 Eigen::MatrixXd& PV1,
										 Eigen::MatrixXd& PV2);

void frame_field_alignment_data(const Eigen::MatrixXd& V,
																const Eigen::MatrixXi& F,
																Eigen::VectorXd& A);

#endif
