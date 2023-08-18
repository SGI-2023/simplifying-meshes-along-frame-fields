#include "qslim_aligned.h"
#include "per_vertex_point_to_plane_quadrics_aligned.h"
#include "qslim_aligned_callbacks.h"

#include "igl/connect_boundary_to_infinity.h"
#include "igl/decimate.h"
#include "igl/decimate_callback_types.h"
#include "igl/edge_flaps.h"
#include "igl/is_edge_manifold.h"
#include "igl/max_faces_stopping_condition.h"
#include "igl/per_vertex_point_to_plane_quadrics.h"
#include "igl/principal_curvature.h"
#include "igl/remove_unreferenced.h"

bool qslim_aligned(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                   const size_t max_m, Eigen::MatrixXd &U, Eigen::MatrixXi &G,
                   Eigen::VectorXi &J, Eigen::VectorXi &I) {
  // Original number of faces
  const int orig_m = F.rows();
  // Tracking number of faces
  int m = F.rows();
  typedef Eigen::MatrixXd DerivedV;
  typedef Eigen::MatrixXi DerivedF;
  DerivedV VO;
  DerivedF FO;
  igl::connect_boundary_to_infinity(V, F, VO, FO);
  // decimate will not work correctly on non-edge-manifold meshes. By extension
  // this includes meshes with non-manifold vertices on the boundary since these
  // will create a non-manifold edge when connected to infinity.
  if (!igl::is_edge_manifold(FO)) {
    return false;
  }
  Eigen::VectorXi EMAP;
  Eigen::MatrixXi E, EF, EI;
  Eigen::MatrixXd PD1, PD2;
  Eigen::VectorXd PV1, PV2;
  igl::principal_curvature(V, F, PD1, PD2, PV1, PV2);
  igl::edge_flaps(FO, E, EMAP, EF, EI);
  // Quadrics per vertex
  typedef std::tuple<Eigen::MatrixXd, Eigen::RowVectorXd, double> Quadric;
  std::vector<Quadric> quadrics;
  per_vertex_point_to_plane_quadrics_aligned(VO, FO, EMAP, EF, EI, PD1, PD2,
                                             PV1, PV2, quadrics);
  // State variables keeping track of edge we just collapsed
  int v1 = -1;
  int v2 = -1;
  // Callbacks for computing and updating metric
  igl::decimate_cost_and_placement_callback cost_and_placement;
  igl::decimate_pre_collapse_callback pre_collapse;
  igl::decimate_post_collapse_callback post_collapse;
  qslim_aligned_callbacks(E, quadrics, v1, v2, PD1, PD2, PV1, PV2,
                          cost_and_placement, pre_collapse, post_collapse);
  // Call to greedy decimator
  bool ret =
      igl::decimate(VO, FO, cost_and_placement,
                    igl::max_faces_stopping_condition(m, orig_m, max_m),
                    pre_collapse, post_collapse, E, EMAP, EF, EI, U, G, J, I);
  // Remove phony boundary faces and clean up
  const Eigen::Array<bool, Eigen::Dynamic, 1> keep = (J.array() < orig_m);
  igl::slice_mask(Eigen::MatrixXi(G), keep, 1, G);
  igl::slice_mask(Eigen::VectorXi(J), keep, 1, J);
  Eigen::VectorXi _1, I2;
  igl::remove_unreferenced(Eigen::MatrixXd(U), Eigen::MatrixXi(G), U, G, _1,
                           I2);
  igl::slice(Eigen::VectorXi(I), I2, 1, I);

  return ret;
}
