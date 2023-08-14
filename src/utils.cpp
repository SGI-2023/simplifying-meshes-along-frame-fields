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
void compute_before(const int e,
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi & E,
	const Eigen::MatrixXi & F,
	const Eigen::VectorXi & EMAP,
	const Eigen::MatrixXi & EF,
	const Eigen::MatrixXi & EI, 
	const Eigen::MatrixXd & v,//frame field 1s
	const Eigen::MatrixXd & w,//frame field 2s
	int & valenceSum, 
  	double & alignment){
		std::vector<int> Nv;
  		std::vector<int> Nf;
		igl::circulation(e, true, F,EMAP,EF,EI, Nv,Nf);
		valenceSum = Nv.size();
		int alignmentTotal = 0;
		for (int i = 0; i < valenceSum; i++) {
			alignmentTotal += alignment_function(v.row(Nv[i]), w.row(Nv[i]),V.row(Nv[i])-V.row(E.col(0)(e)));
		}
		igl::circulation(e, false, F,EMAP,EF,EI, Nv,Nf);
		valenceSum += Nv.size();
		for (int i = 1; i < valenceSum; i++) {
			
			alignmentTotal += alignment_function(v.row(Nv[i]), w.row(Nv[i]),V.row(Nv[i])-V.row(E.col(1)(e)));
		}
		alignment = alignmentTotal/(valenceSum-1);
  }
// cost function, thinking about helper function
	// needs: p1, p2, valence p1, valence p2,
	// 		avged alignment for surroundings pre-collapse, 
	// 		avged alignment for surroundings post-collapse.--needs new point.
// double avg_alignment(v1, list_of_edges_indexed_into_vertex_list, &vertex_list, &frame_field_list)
//		takes in a list of edges (see below), calls alignment function on all of them and avgs/normalizes.
//		returns this value. 
// void get_edges(vertex_list, vertex_index_from, vertex_index_exclude, &edge_list)
// 		returns list of edges (indexed into vertices) emanating from vertex_index_from
//		except for the one that goes to vertex_index_exclude. 
//		If vertex_index_exclude is not present, then just return all of the edges
//		returns them in &edge_list
// void get_all_edges(vertex_list, edge_index, &edge_list)
//		use the above function to get all edges emanating from a pair of vertices/an edge
//		return the edges in &edge_list