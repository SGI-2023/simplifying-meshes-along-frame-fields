#ifndef UTILS_H
#define UTILS_H

#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/decimate.h>
#include <Eigen/Core>
#include <igl/circulation.h>

#include <iostream>
#include <string>
#include <cmath>

void view_mesh(const Eigen::MatrixXd& V,
							 const Eigen::MatrixXi& F);

void normal_vector(const Eigen::RowVector3d& v,
									 const Eigen::RowVector3d& w,
									 Eigen::RowVector3d& r);

void projection_onto_plane(const Eigen::RowVector3d& v,
													 const Eigen::RowVector3d& n,
													 Eigen::RowVector3d& r);

void bisector_vector(const Eigen::RowVector3d& v,
										 const Eigen::RowVector3d& w,
										 Eigen::RowVector3d& r);

double angle_between_vectors(const Eigen::RowVector3d& v,
														 const Eigen::RowVector3d w);

double alignment_function(const Eigen::RowVector3d& v,
													const Eigen::RowVector3d& w,
													const Eigen::RowVector3d& p_edge);
void compute_before(const int e,
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi & E,
	const Eigen::MatrixXi & F,
	const Eigen::VectorXi & EMAP,
	const Eigen::MatrixXi & EF,
	const Eigen::MatrixXi & EI, 
	const Eigen::MatrixXi & v,//frame field 1s
	const Eigen::MatrixXi & w,//frame field 2s
	int & valenceSum, 
  	double & alignment);
#endif
