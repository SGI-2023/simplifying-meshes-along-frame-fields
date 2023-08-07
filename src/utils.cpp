#include "utils.h"

void view_mesh(const Eigen::MatrixXd& V,
							 const Eigen::MatrixXi& F)
{
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V, F);
	viewer.launch();
}

