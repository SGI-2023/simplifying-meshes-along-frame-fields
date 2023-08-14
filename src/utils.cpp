#include "utils.h"

void view_mesh(const Eigen::MatrixXd& V,
							 const Eigen::MatrixXi& F)
{
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V, F);
	viewer.launch();
}

void normal_vector(const Eigen::RowVector3d& v,
									 const Eigen::RowVector3d& w,
									 Eigen::RowVector3d& r)
{
	r = v.cross(w);
}

void projection_onto_plane(const Eigen::RowVector3d& v,
													 const Eigen::RowVector3d& n,
													 Eigen::RowVector3d& r)
{
	r = v - (v.dot(n)/n.norm()) * n;
}

void bisector_vector(const Eigen::RowVector3d& v,
										 const Eigen::RowVector3d& w,
										 Eigen::RowVector3d& r)
{
	r = (1/v.norm()) * v + (1/w.norm()) * w;
}

double angle_between_vectors(const Eigen::RowVector3d& v,
														 const Eigen::RowVector3d w)
{
	return std::acos( v.dot(w)/(v.norm() * w.norm()) );
}

double alignment_function(const Eigen::RowVector3d& v,
													const Eigen::RowVector3d& w,
													const Eigen::RowVector3d& p_edge)
{
	double alpha = angle_between_vectors(v, w)/2;
	double beta = M_PI - alpha;
	Eigen::RowVector3d u = Eigen::RowVector3d::Zero();
	bisector_vector(v, w, u);

	// Projecting p_edge onto the the plane defined by v and w.
	Eigen::RowVector3d plane_normal;
	normal_vector(v, w, plane_normal);
	Eigen::RowVector3d p;
	projection_onto_plane(p_edge, plane_normal, p);

	// Angle between p and bisector u on the plane defined by v and w.
	double theta = angle_between_vectors(p, u);

	double f = 1;
	if(theta < alpha)
		f = std::pow(std::cos((M_PI/(2*alpha))*theta), 2);
	else if(theta < alpha+2*beta)
		f = std::pow(std::cos(M_PI/(2*beta) * (theta-M_PI/2)), 2);
	else if(theta < M_PI)
		f = std::pow(std::cos(M_PI/(2*alpha) * (theta-M_PI)), 2);
	
	return f;
}
