#include "utils.h"

void normal_vector(const RowVector3d &v,
				   const RowVector3d &w,
				   RowVector3d &r)
{
	r = v.cross(w);
}

void projection_onto_plane(const RowVector3d &v,
						   const RowVector3d &n,
						   RowVector3d &r)
{
	r = v - (v.dot(n) / n.norm()) * n;
}

void bisector_vector(const RowVector3d &v,
					 const RowVector3d &w,
					 RowVector3d &r)
{
	r = (1 / v.norm()) * v + (1 / w.norm()) * w;
}

double angle_between_vectors(const RowVector3d &v,
							 const RowVector3d w)
{
	return std::acos(v.dot(w) / (v.norm() * w.norm()));
}

double alignment_function(const RowVector3d &v,
						  const RowVector3d &w,
						  const RowVector3d &p_edge)
{
	double alpha = angle_between_vectors(v, w) / 2;
	double beta = M_PI - alpha;
	RowVector3d u = RowVector3d::Zero();
	bisector_vector(v, w, u);

	// Projecting p_edge onto the the plane defined by v and w.
	RowVector3d plane_normal;
	normal_vector(v, w, plane_normal);
	RowVector3d p;
	projection_onto_plane(p_edge, plane_normal, p);

	// Angle between p and bisector u on the plane defined by v and w.
	double theta = angle_between_vectors(p, u);

	double f = 1;
	if (theta < alpha)
		f = std::pow(std::cos((M_PI / (2 * alpha)) * theta), 2);
	else if (theta < alpha + 2 * beta)
		f = std::pow(std::cos(M_PI / (2 * beta) * (theta - M_PI / 2)), 2);
	else if (theta < M_PI)
		f = std::pow(std::cos(M_PI / (2 * alpha) * (theta - M_PI)), 2);

	return 1-f;
}
// double alignment_function(const RowVector3d &v,
// 						  const RowVector3d &w,
// 						  const RowVector3d &p_edge, 
// 						  bool & v_align)
// {
// 	RowVector3d plane_normal;
// 	normal_vector(v, w, plane_normal);
// 	RowVector3d p;
// 	projection_onto_plane(p_edge, plane_normal, p);
// 	double theta1 = angle_between_vectors(p, v);
// 	double theta2 = angle_between_vectors(p, w);
// 	v_align = theta1 < theta2;
// 	return theta1 < theta2 ? theta1 : theta2;
// }
double alignment_function(const RowVector3d &v1,
						  const RowVector3d &w1,
						  const RowVector3d &v2,
						  const RowVector3d &w2,
						  const RowVector3d &p_edge){
	if (p_edge[0] == 0 && p_edge[1] == 0 && p_edge[2]==0) {
		return 0;
	}
	RowVector3d plane_normal1;
	RowVector3d plane_normal2;
	double thetav1;
	double thetaw1;
	double thetav2;
	double thetaw2;
	RowVector3d p1;
	RowVector3d p2;
	double pi =  3.14159;
    // v1 = PD1.row(E.col(0)(e));
    // w1 = PD2.row(E.col(0)(e));
    // RowVector3d p_edge = p1-p2;
    normal_vector(v1, w1, plane_normal1);
    projection_onto_plane(p_edge, plane_normal1, p1);
    thetav1 = angle_between_vectors(p1, v1);
    thetaw1 = angle_between_vectors(p1, w1);
    thetav1 = thetav1 > pi/2 ? pi - thetav1 : thetav1;
    thetaw1 = thetaw1 > pi/2 ? pi - thetaw1 : thetaw1;  
    // v2 = PD1.row(E.col(1)(e));
    // w2 = PD2.row(E.col(1)(e));
    normal_vector(v2, w2, plane_normal2);
    projection_onto_plane(p_edge, plane_normal2, p2);
    thetav2 = angle_between_vectors(p2, v2);
    thetaw2 = angle_between_vectors(p2, w2);
    thetav2 = thetav2 > pi/2 ? pi - thetav2 : thetav2;
    thetaw2 = thetaw2 > pi/2 ? pi - thetaw2 : thetaw2; 

	// std::cout << "p_edge: " << p_edge[0] << std::endl;

	// std::cout << "v1: " << v1[0] << std::endl;
	// std::cout << "w1: " << w1[0] << std::endl;
	// std::cout << "v2: " << v2[0] << std::endl;
	// std::cout << "w2: " << w2[0] << std::endl;

	// std::cout << "plane_normal1: " << plane_normal1[0] << std::endl;
	// std::cout << "plane_normal2: " << plane_normal2[0] << std::endl;
	// std::cout << "p1: " << p1[0] << std::endl;
	// std::cout << "p2: " << p2[0] << std::endl;

	// std::cout << "thetav1: " << thetav1 << std::endl;
	// std::cout << "thetaw1: " << thetaw1 << std::endl;
	// std::cout << "thetav2: " << thetav2 << std::endl;
	// std::cout << "thetaw2: " << thetaw2 << std::endl;
	// std::exit(1);
    double val = (thetav1 + thetav2)/2 < (thetaw1 + thetaw2)/2 ? (thetav1 + thetav2)/2 : (thetaw1 + thetaw2)/2;
	if (val != val) {
			std::cout << "p_edge: " << p_edge[0] << std::endl;

	std::cout << "v1: " << v1[0] << std::endl;
	std::cout << "w1: " << w1[0] << std::endl;
	std::cout << "v2: " << v2[0] << std::endl;
	std::cout << "w2: " << w2[0] << std::endl;

	std::cout << "plane_normal1: " << plane_normal1[0] << std::endl;
	std::cout << "plane_normal2: " << plane_normal2[0] << std::endl;
	std::cout << "p1: " << p1[0] << std::endl;
	std::cout << "p2: " << p2[0] << std::endl;

		std::cout << "thetav1: " << thetav1 << std::endl;
	std::cout << "thetaw1: " << thetaw1 << std::endl;
	std::cout << "thetav2: " << thetav2 << std::endl;
	std::cout << "thetaw2: " << thetaw2 << std::endl;
		std::exit(1);
	}
	return val;
}
