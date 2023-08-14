#include "igl/circulation.h"
#include "igl/decimate_callback_types.h"
#include "igl/shortest_edge_and_midpoint.h"
#include <Eigen/src/Core/Matrix.h>

void optimal_collapse_edge_callbacks(
    Eigen::MatrixXi &E,
    /* std::vector<std::tuple<Eigen::MatrixXd, Eigen::RowVectorXd, double>>
     */
    /*     &quadrics, */
    Eigen::MatrixXd &PD1, Eigen::MatrixXd &PD2, Eigen::VectorXd &PV1,
    Eigen::VectorXd &PV2, int &v1, int &v2,
    igl::decimate_cost_and_placement_callback &cost_and_placement,
    igl::decimate_pre_collapse_callback &pre_collapse,
    igl::decimate_post_collapse_callback &post_collapse) {
  /* typedef std::tuple<Eigen::MatrixXd, Eigen::RowVectorXd, double> Quadric; */
  cost_and_placement = igl::shortest_edge_and_midpoint;
  /* [&PD1, &PD2, &PV1, &PV2, &v1, &v2]( */
  //     const int e, const Eigen::MatrixXd &V, const Eigen::MatrixXi & /*F*/,
  //     const Eigen::MatrixXi &E, const Eigen::VectorXi & /*EMAP*/,
  //     const Eigen::MatrixXi & /*EF*/, const Eigen::MatrixXi & /*EI*/,
  /*     double &cost, Eigen::RowVectorXd &p) { */
  /* // Combined quadric */
  /* Quadric quadric_p; */
  /* quadric_p = quadrics[E(e, 0)] + quadrics[E(e, 1)]; */
  /* // Quadric: p'Ap + 2b'p + c */
  /* // optimal point: Ap = -b, or rather because we have row vectors:
   * pA=-b */
  /* const auto &A = std::get<0>(quadric_p); */
  /* const auto &b = std::get<1>(quadric_p); */
  /* const auto &c = std::get<2>(quadric_p); */
  /* p = -b * A.inverse(); */
  /* cost = p.dot(p * A) + 2 * p.dot(b) + c; */
  /* // Force infs and nans to infinity */
  /* if (std::isinf(cost) || cost != cost) { */
  /*   cost = std::numeric_limits<double>::infinity(); */
  /*   // Prevent NaNs. Actually NaNs might be useful for debugging. */
  /*   p.setConstant(0); */
  /* } */
  /* }; */
  // Remember endpoints
  pre_collapse =
      [&v1, &v2](const Eigen::MatrixXd &,                             /*V*/
                 const Eigen::MatrixXi &,                             /*F*/
                 const Eigen::MatrixXi &E, const Eigen::VectorXi &,   /*EMAP*/
                 const Eigen::MatrixXi &,                             /*EF*/
                 const Eigen::MatrixXi &,                             /*EI*/
                 const igl::min_heap<std::tuple<double, int, int>> &, /*Q*/
                 const Eigen::VectorXi &,                             /*EQ*/
                 const Eigen::MatrixXd &,                             /*C*/
                 const int e) -> bool {
    v1 = E(e, 0);
    v2 = E(e, 1);
    return true;
  };
  // update quadric
  post_collapse = [&v1, &v2, &PD1, &PD2, &PV1, &PV2](
                      const Eigen::MatrixXd &,     /*V*/
                      const Eigen::MatrixXi &F,    /*F*/
                      const Eigen::MatrixXi &,     /*E*/
                      const Eigen::VectorXi &EMAP, /*EMAP*/
                      const Eigen::MatrixXi &EF,   /*EF*/
                      const Eigen::MatrixXi &EI,   /*EI*/
                      const igl::min_heap<std::tuple<double, int, int>> &, /*Q*/
                      const Eigen::VectorXi &, /*EQ*/
                      const Eigen::MatrixXd &, /*C*/
                      const int e,             /*e*/
                      const int,               /*e1*/
                      const int,               /*e2*/
                      const int,               /*f1*/
                      const int,               /*f2*/
                      const bool collapsed) -> void {
    std::vector<int> neighbor_vertices, neighbor_faces;
    igl::circulation(e, false, F, EMAP, EF, EI, neighbor_vertices,
                     neighbor_faces);
    Eigen::RowVectorXd avg_neighbor_PD1(neighbor_vertices.size()),
        avg_neighbor_PD2(neighbor_vertices.size());
    double avg_neighbor_PV1 = 0;
    double avg_neighbor_PV2 = 0;

    for (int i = 0; i < neighbor_vertices.size(); i++) {
      for (int j = 0; j < PD1.cols(); j++) {
        avg_neighbor_PD1(j) += PD1(neighbor_vertices[i], j);
        avg_neighbor_PD2(j) += PD2(neighbor_vertices[i], j);
      }
      avg_neighbor_PV1 += PV1(neighbor_vertices[i]);
      avg_neighbor_PV2 += PV2(neighbor_vertices[i]);
    }

    avg_neighbor_PD1 /= neighbor_vertices.size();
    avg_neighbor_PD2 /= neighbor_vertices.size();
    avg_neighbor_PV1 /= neighbor_vertices.size();
    avg_neighbor_PV2 /= neighbor_vertices.size();

    if (v1 < v2) {
      for (int j = 0; j < PD1.cols(); j++) {
        PD1(v1, j) = avg_neighbor_PD1(j);
        PD2(v1, j) = avg_neighbor_PD2(j);
      }
    } else {
      for (int j = 0; j < PD1.cols(); j++) {
        PD1(v2, j) = avg_neighbor_PD1(j);
        PD2(v2, j) = avg_neighbor_PD2(j);
      }
    }
    PV1(v1 < v2 ? v1 : v2) = avg_neighbor_PV1;
    PV2(v1 < v2 ? v1 : v2) = avg_neighbor_PV2;
  };
}
