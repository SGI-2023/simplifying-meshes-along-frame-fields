#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/decimate.h>
#include <Eigen/Core>

#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
	std::string filename = "../data/spot.obj";

	Eigen::MatrixXd V, U;
	Eigen::MatrixXi F, G;
	Eigen::VectorXi J, I;

	igl::read_triangle_mesh(filename, V, F);

	int num_faces = (F.rows())/4;
	std::cout << "Decimating into 1/4 of faces." << std::endl;
	igl::decimate(V, F, num_faces, U, G, J, I);

	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(U, G);
	viewer.launch();

	return 0;
}
