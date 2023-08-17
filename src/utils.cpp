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
