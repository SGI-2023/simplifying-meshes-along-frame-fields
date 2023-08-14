#pragma once

#include <Eigen/Core>
#include <iostream>
#include <string>
#include <cmath>

using namespace Eigen;

void normal_vector(const RowVector3d &v,
				   const RowVector3d &w,
				   RowVector3d &r);

void projection_onto_plane(const RowVector3d &v,
						   const RowVector3d &n,
						   RowVector3d &r);

void bisector_vector(const RowVector3d &v,
					 const RowVector3d &w,
					 RowVector3d &r);

double angle_between_vectors(const RowVector3d &v,
							 const RowVector3d w);

double alignment_function(const RowVector3d &v,
						  const RowVector3d &w,
						  const RowVector3d &p_edge);
