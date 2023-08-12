// Source:
// https://github.com/libigl/libigl/blob/main/tutorial/203_CurvatureDirections/main.cpp

#include <igl/avg_edge_length.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/parula.h>
#include <igl/per_corner_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/principal_curvature.h>
#include <igl/read_triangle_mesh.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

int main(int argc, char *argv[]) {
  using namespace Eigen;
  std::string filename = "../data/spot.obj";
  if (argc > 1) {
    filename = argv[1];
  }
  // Load a mesh in OFF format
  igl::read_triangle_mesh(filename, V, F);

  // Alternative discrete mean curvature
  MatrixXd HN;
  SparseMatrix<double> L, M, Minv;
  igl::cotmatrix(V, F, L);
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
  igl::invert_diag(M, Minv);
  // Laplace-Beltrami of position
  HN = -Minv * (L * V);
  // Extract magnitude as mean curvature
  VectorXd H = HN.rowwise().norm();

  // Compute curvature directions via quadric fitting
  MatrixXd PD1, PD2; // V * 3      Where V is the number of the Vertices
  VectorXd PV1, PV2; // V          Where V is the number of the Vertices
  igl::principal_curvature(V, F, PD1, PD2, PV1, PV2);
  // mean curvature
  H = 0.5 * (PV1 + PV2);

  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);

  viewer.data().set_data(H);

  // Average edge length for sizing
  const double avg = igl::avg_edge_length(V, F);
  const double factor = avg * 0.1;

  // Draw a red segment parallel to the maximal curvature direction
  const RowVector3d red(0.8, 0.2, 0.2), blue(0.2, 0.2, 0.8);
  viewer.data().add_edges(
      V + (PD1.array().colwise() * PV1.array()).matrix() * factor,
      V - (PD1.array().colwise() * PV1.array()).matrix() * factor, red);

  // Draw a blue segment parallel to the minimal curvature direction
  viewer.data().add_edges(
      V + (PD2.array().colwise() * PV2.array()).matrix() * factor,
      V - (PD2.array().colwise() * PV2.array()).matrix() * factor, blue);

  // Hide wireframe
  viewer.data().show_lines = false;

  viewer.launch();
}
