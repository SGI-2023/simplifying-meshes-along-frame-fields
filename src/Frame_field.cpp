#include "Frame_field.h"


void get_frame_field(const Eigen::MatrixXd& V,
										 const Eigen::MatrixXi& F,
										 Eigen::MatrixXd& PD1,
										 Eigen::MatrixXd& PD2,
										 Eigen::MatrixXd& PV1,
										 Eigen::MatrixXd& PV2)
{
	Eigen::SparseMatrix<double> L, M, Minv;
	Eigen::MatrixXd HN;

	igl::cotmatrix(V, F, L);
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
	igl::invert_diag(M, Minv);

	// Laplace-Beltrami of position.
	HN = -Minv * (L * V);

	// Exact magnitude as mean curvature
	Eigen::VectorXd H = HN.rowwise().norm();

	igl::principal_curvature(V, F, PD1, PD2, PV1, PV2);
}

void frame_field_alignment_data(const Eigen::MatrixXd& V,
																const Eigen::MatrixXi& F,
																Eigen::VectorXd& A)
{
	Eigen::MatrixXd PD1, PD2;
	Eigen::MatrixXd PV1, PV2;

	get_frame_field(V, F, PD1, PD2, PV1, PV2);

	// Initialize resulting data.
	A = Eigen::VectorXd(V.rows());
	for(int i = 0; i < V.rows(); ++i)
		A[i] = 1;

	for(int i = 0; i < F.rows(); ++i){
		Eigen::RowVector3i face = F.row(i);
		
		for(int j = 0; j < 3; ++j){
			int adj1 = face[(i+1)%3];
			int adj2 = face[(i+2)%3];
			Eigen::RowVector3d P = V.row(adj1);
			Eigen::RowVector3d Q = V.row(adj2);
			Eigen::RowVector3d p = Q - P;

			Eigen::RowVector3d v1 = PD1.row(adj1);
			Eigen::RowVector3d w1 = PD2.row(adj1);
			Eigen::RowVector3d v2 = PD1.row(adj2);
			Eigen::RowVector3d w2 = PD2.row(adj2);

			double f1 = alignment_function(v1, w1, p);
			double f2 = alignment_function(v2, w2, p);
		
			A[adj1] = f1 < A[adj1] ? f1 : A[adj1];
			A[adj2] = f2 < A[adj2] ? f2 : A[adj2];
		}
	}
}
