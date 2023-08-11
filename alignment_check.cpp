#include <igl/read_triangle_mesh.h>
#include <igl/avg_edge_length.h>
#include <igl/opengl/glfw/Viewer.h>
#include <string>

#include "Frame_field.h"

int main(int argc, char* argv[])
{
	std::string filename = argc < 2 ? "../data/spot.obj" : argv[1];

	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	igl::read_triangle_mesh(filename, V, F);

	Eigen::MatrixXd PD1, PD2;
	Eigen::MatrixXd PV1, PV2;

	get_frame_field(V, F, PD1, PD2, PV1, PV2);

	Eigen::VectorXd A;
	frame_field_alignment_data(V, F, A);
	std::cout << "A = \n" << A << std::endl;

	const double avg = igl::avg_edge_length(V, F);
	const double factor = avg * 0.1;
	const Eigen::RowVector3d red(0.8, 0.2, 0.2), blue(0.2, 0.2, 0.8);

	igl::opengl::glfw::Viewer vw;
	
	vw.callback_key_pressed = 
		[&](igl::opengl::glfw::Viewer &, unsigned char key, int)->bool
	{
		switch(key)
		{
			case ' ':
				vw.data().add_edges(
					V.array() + (PD1.array() * PV1.replicate(1, 3).array() * factor),
					V.array() - (PD1.array() * PV1.replicate(1, 3).array() * factor),
					red
				);

				vw.data().add_edges(
					V.array() + (PD2.array() * PV2.replicate(1, 3).array() * factor),
					V.array() - (PD2.array() * PV2.replicate(1, 3).array() * factor),
					blue
				);
				return true;
			case 'R':
			case 'r':
				vw.data().clear();
				vw.data().set_mesh(V, F);
				vw.core().align_camera_center(V, F);
				vw.data().set_data(A);
		}
		return false;
	};
	
	vw.data().set_mesh(V, F);
	vw.core().align_camera_center(V, F);
	vw.data().set_data(A);
	vw.launch();
}
