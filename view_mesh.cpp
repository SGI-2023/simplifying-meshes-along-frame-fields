#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/decimate.h>
#include <Eigen/Core>

#include <iostream>
#include <string>

#include "utils.h"

int main(int argc, char* argv[])
{
	std::string filename = argc < 3 ? "../data/spot.obj" : argv[2];

	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	igl::read_triangle_mesh(filename, V, F);

	view_mesh(V, F);

	return 0;
}
