#include <igl/circulation.h>
#include <igl/collapse_edge.h>
#include <igl/edge_flaps.h>
#include <igl/resolve_duplicated_faces.h>
#include <igl/decimate.h>
#include <igl/invert_diag.h>
#include <igl/principal_curvature.h>
#include <igl/invert_diag.h>
#include <igl/avg_edge_length.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/parallel_for.h>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <iostream>
#include <set>

using namespace Eigen;
using namespace std;
using namespace igl;

class InteractiveVisualizer {
  igl::opengl::glfw::Viewer viewer;

  bool toggle_frame_fields = false;
  MatrixXd V, OV;
  MatrixXi F, OF;

  // Prepare array-based edge data structures and priority queue
  VectorXi EMAP;
  MatrixXi E, EF, EI;
  igl::min_heap<std::tuple<double, int, int>> Q;
  Eigen::VectorXi EQ;
  // If an edge were collapsed, we'd collapse it to these points:
  MatrixXd C;
  int num_collapsed;

  void visualizeFrameFields(RowVector3d maxCurvatureColor, RowVector3d minCurvatureColor)
  {
    //std::cout << "framefields:" << V.rows() << " " << F.rows() << std::endl;

    // Compute curvature directions via quadric fitting
    MatrixXd PD1, PD2; // V * 3      Where V is the number of the Vertices
    VectorXd PV1, PV2; // V          Where V is the number of the Vertices
    igl::principal_curvature(V, F, PD1, PD2, PV1, PV2);
    // mean curvature
    // H = 0.5 * (PV1 + PV2);

    // viewer.data().set_data(H);

    // Average edge length for sizing
    const double avg = igl::avg_edge_length(V, F);
    const double factor = avg * 0.1;

    // Draw a red segment parallel to the maximal curvature direction
    viewer.data().add_edges(
      V.array() + (PD1.array() * PV1.replicate(1, 3).array() * factor),
      V.array() - (PD1.array() * PV1.replicate(1, 3).array() * factor),
      maxCurvatureColor);

    // Draw a blue segment parallel to the minimal curvature direction
    viewer.data().add_edges(
      V.array() + (PD2.array() * PV2.replicate(1, 3).array() * factor),
      V.array() - (PD2.array() * PV2.replicate(1, 3).array() * factor),
      minCurvatureColor);
  }

  void reset()
  {
    F = OF;
    V = OV;
    edge_flaps(F, E, EMAP, EF, EI);
    C.resize(E.rows(), V.cols());
    VectorXd costs(E.rows());
    // https://stackoverflow.com/questions/2852140/priority-queue-clear-method
    // Q.clear();
    Q = {};
    EQ = Eigen::VectorXi::Zero(E.rows());
    {
      Eigen::VectorXd costs(E.rows());
      igl::parallel_for(
        E.rows(), [&](const int e)
      {
        double cost = e;
        RowVectorXd p(1, 3);
        shortest_edge_and_midpoint(e, V, F, E, EMAP, EF, EI, cost, p);
        C.row(e) = p;
        costs(e) = cost; },
        10000);
      for (int e = 0; e < E.rows(); e++)
      {
        Q.emplace(costs(e), e, 0);
      }
    }

    num_collapsed = 0;
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
    viewer.data().set_face_based(true);
  };

  void decimation()
  {
    if (!Q.empty())
    {
      bool something_collapsed = false;
      // collapse edge
      const int max_iter = std::ceil(0.01 * Q.size());
      for (int j = 0; j < max_iter; j++)
      {
        if (!collapse_edge(shortest_edge_and_midpoint, V, F, E, EMAP, EF, EI, Q, EQ, C))
        {
          break;
        }
        //resolve_duplicated_faces(F);
        something_collapsed = true;
        num_collapsed++;
      }

      if (something_collapsed)
      {
        //std::cout << V << " " << F << std::endl;
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        viewer.data().set_face_based(true);
      }
    }
  };

  const auto& key_down =
    [](igl::opengl::glfw::Viewer& viewer, unsigned char key, int mod) -> bool
  {
    switch (key)
    {
    case 'C':
    case 'c':
      if (toggle_frame_fields)
      {
        viewer.data().clear_edges();
      }
      else
      {
        visualizeFrameFields(RowVector3d(0.8, 0.2, 0.2), RowVector3d(0.2, 0.2, 0.8));
      }
      toggle_frame_fields = !toggle_frame_fields;
      break;
    case 'D':
    case 'd':
      decimation();
      break;
    case 'R':
    case 'r':
      reset();
      break;
    default:
      return false;
    }
    return true;
  };
};



//int main(int argc, char *argv[])
//{
//  cout << "  'r'  reset." << endl;
//  // Load a closed manifold mesh
//  string filename = "../data/spot.obj";
//  if (argc >= 2)
//  {
//    filename = argv[1];
//  }
//  read_triangle_mesh(filename, OV, OF);
//
//  reset();
//  viewer.core().background_color.setConstant(1);
//  viewer.callback_key_down = key_down;
//  return viewer.launch();
//}