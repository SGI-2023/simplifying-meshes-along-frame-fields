#include "../include/window.h"

void Window::visualizeFrameFields(RowVector3d maxCurvatureColor, RowVector3d minCurvatureColor)
{
  mesh.get_frame_field();
  // Average edge length for sizing
  const double avg = igl::avg_edge_length(mesh.V, mesh.F);
  const double factor = avg * 0.1;

  // Draw a red segment parallel to the maximal curvature direction
  viewer.data().add_edges(
      mesh.V.array() + (mesh.PD1.array() * mesh.PV1.replicate(1, 3).array() * factor),
      mesh.V.array() - (mesh.PD1.array() * mesh.PV1.replicate(1, 3).array() * factor),
      maxCurvatureColor);

  // Draw a blue segment parallel to the minimal curvature direction
  viewer.data().add_edges(
      mesh.V.array() + (mesh.PD2.array() * mesh.PV2.replicate(1, 3).array() * factor),
      mesh.V.array() - (mesh.PD2.array() * mesh.PV2.replicate(1, 3).array() * factor),
      minCurvatureColor);
}

void Window::visualizeMetricOnEdges()
{
  // Get edges indexed into vertex ids
  MatrixXi E;
  igl::edges(mesh.F, E);
  //std::cout << mesh.F.rows() << " " << E.rows() << " " << E.cols() << std::endl;
  // Z is a dummy measure that should calculate the length of each edge
  VectorXd Z;
  Z.resize(E.rows());
  for (int e = 0; e < E.rows(); e++)
  {
    Z(e) = (mesh.V.row(E.col(0)(e)) - mesh.V.row(E.col(1)(e))).norm();
  }
  // Generate colormap based on Z
  MatrixXd C;
  igl::colormap(igl::COLOR_MAP_TYPE_VIRIDIS, Z, true, C);
  viewer.data().show_lines = false;
  viewer.data().set_edges(mesh.V,E,C);
}

void Window::initialize_callbacks()
{
  const auto &key_down =
      [&](igl::opengl::glfw::Viewer &viewer, unsigned char key, int mod) -> bool
  {
    switch (key)
    {
    case 'C':
    case 'c':
      visualizeFrameFields(RowVector3d(0.8, 0.2, 0.2), RowVector3d(0.2, 0.2, 0.8));
      break;
    case 'V':
    case 'v':
      mesh.decimate();

      viewer.data().clear();
      viewer.data().set_mesh(mesh.V, mesh.F);
      viewer.data().set_face_based(true);
      break;
    case 'R':
    case 'r':
      mesh.frame_field_alignment_data();

      viewer.data().clear();
      viewer.data().set_mesh(mesh.V, mesh.F);
      viewer.core().align_camera_center(mesh.V, mesh.F);
      viewer.data().set_data(mesh.A);
    case 'K':
    case 'k':
      visualizeMetricOnEdges();
    default:
      return false;
    }
    return true;
  };

  viewer.callback_key_down = key_down;
};