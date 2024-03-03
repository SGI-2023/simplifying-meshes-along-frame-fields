#include "qslim_aligned_callbacks.h"
#include "igl/decimate_callback_types.h"
#include "igl/quadric_binary_plus_operator.h"
#include "utils.h"
#include <Eigen/LU>

void qslim_aligned_callbacks(
    Eigen::MatrixXi &E,
    std::vector<std::tuple<Eigen::MatrixXd, Eigen::RowVectorXd, double>>
        &quadrics,
    int &v1, int &v2, Eigen::MatrixXd &PD1, Eigen::MatrixXd &PD2,
    Eigen::VectorXd &PV1, Eigen::VectorXd &PV2,
    igl::decimate_cost_and_placement_callback &cost_and_placement,
    igl::decimate_pre_collapse_callback &pre_collapse,
    igl::decimate_post_collapse_callback &post_collapse) {
  typedef std::tuple<Eigen::MatrixXd, Eigen::RowVectorXd, double> Quadric;
  cost_and_placement =
      [&quadrics, &PD1, &PD2](
          const int e, const Eigen::MatrixXd &V, const Eigen::MatrixXi & /*F*/,
          const Eigen::MatrixXi &E, const Eigen::VectorXi & /*EMAP*/,
          const Eigen::MatrixXi & /*EF*/, const Eigen::MatrixXi & /*EI*/,
          double &cost, Eigen::RowVectorXd &p) {
        // Combined quadric
        Quadric quadric_p;
        quadric_p = igl::operator+(quadrics[E(e, 0)], quadrics[E(e, 1)]);
        // Quadric: p'Ap + 2b'p + c
        // optimal point: Ap = -b, or rather because we have row vectors: pA=-b
        const auto &A = std::get<0>(quadric_p);
        const auto &b = std::get<1>(quadric_p);
        const auto &c = std::get<2>(quadric_p);
        p = -b * A.inverse();
        cost = (p.dot(p * A) + 2 * p.dot(b) + c) *
               (alignment_function(PD1.row(E(e, 0)), PD2.row(E(e, 0)),
                                   V.row(E(e, 1)) - V.row(E(e, 0))) +
                1);
        // Force infs and nans to infinity
        if (std::isinf(cost) || cost != cost) {
          cost = std::numeric_limits<double>::infinity();
          // Prevent NaNs. Actually NaNs might be useful for debugging.
          p.setConstant(0);
        }
      };
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
  post_collapse = [&v1, &v2, &PD1, &PD2, &PV1, &PV2, &quadrics](
                      const Eigen::MatrixXd &, /*V*/
                      const Eigen::MatrixXi &, /*F*/
                      const Eigen::MatrixXi &, /*E*/
                      const Eigen::VectorXi &, /*EMAP*/
                      const Eigen::MatrixXi &, /*EF*/
                      const Eigen::MatrixXi &, /*EI*/
                      const igl::min_heap<std::tuple<double, int, int>> &, /*Q*/
                      const Eigen::VectorXi &, /*EQ*/
                      const Eigen::MatrixXd &, /*C*/
                      const int,               /*e*/
                      const int,               /*e1*/
                      const int,               /*e2*/
                      const int,               /*f1*/
                      const int,               /*f2*/
                      const bool collapsed) -> void {
    if (collapsed) {
      quadrics[v1 < v2 ? v1 : v2] = igl::operator+(quadrics[v1], quadrics[v2]);

      // update frame field by averaging that of vertices of collapsed edge
      PD1.row(v1 < v2 ? v1 : v2) = 0.5 * (PD1.row(v1) + PD1.row(v2));
      PV1.row(v1 < v2 ? v1 : v2) = 0.5 * (PV1.row(v1) + PV1.row(v2));
      PD2.row(v1 < v2 ? v1 : v2) = 0.5 * (PD2.row(v1) + PD2.row(v2));
      PV2.row(v1 < v2 ? v1 : v2) = 0.5 * (PV2.row(v1) + PV2.row(v2));

      // Normalize
      PD1.row(v1 < v2 ? v1 : v2).normalize();
      PD2.row(v1 < v2 ? v1 : v2).normalize();
    }
  };
}
