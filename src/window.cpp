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
    default:
      return false;
    }
    return true;
  };

  viewer.callback_key_down = key_down;
};