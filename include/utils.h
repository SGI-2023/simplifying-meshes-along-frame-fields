#ifndef UTILS_H
#define UTILS_H

#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/decimate.h>
#include <Eigen/Core>

#include <iostream>
#include <string>

void view_mesh(const Eigen::MatrixXd& V,
							 const Eigen::MatrixXi& F);

#endif
