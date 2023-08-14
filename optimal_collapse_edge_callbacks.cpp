#include "igl/decimate_callback_types.h"

std::vector<int> V2Fe;

const igl::decimate_pre_collapse_callback custom_pre_collapse_callback =
    [&V2Fe](const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
            const Eigen::MatrixXi &E, const Eigen::VectorXi &EMAP,
            const Eigen::MatrixXi &EF, const Eigen::MatrixXi &EI,
            const igl::min_heap<std::tuple<double, int, int>> &Q,
            const Eigen::VectorXi &EQ, const Eigen::MatrixXd &C,
            const int e) -> bool {
  V2Fe = igl::circulation(e, true, EMAP, EF, EI);
  return true;
};

const igl::decimate_post_collapse_callback custom_post_collapse_callback =
    [&PD1, &PD2, &PV1, &PV2,
     &V2Fe](const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
            const Eigen::MatrixXi &E, const Eigen::VectorXi &EMAP,
            const Eigen::MatrixXi &EF, const Eigen::MatrixXi &EI,
            const igl::min_heap<std::tuple<double, int, int>> &Q,
            const Eigen::VectorXi &EQ, const Eigen::MatrixXd &C, const int e,
            const int e1, const int e2, const int f1, const int f2,
            const bool collapsed) -> bool {
  if (!collapsed)
    return false;

  int v1 = E(e, 0);
  int v2 = E(e, 1);

  PD1.row(v1) = 0.5 * (PD1.row(v1) + PD1.row(v2));
  PV1.row(v1) = 0.5 * (PV1.row(v1) + PV1.row(v2));
  PD2.row(v1) = 0.5 * (PD2.row(v1) + PD2.row(v2));
  PV2.row(v1) = 0.5 * (PV2.row(v1) + PV2.row(v2));

  // 4 normalized vectors: PD1.row(v1), PD2.row(v1), - PD1.row(v1), -
  // PD2.row(v1)
  std::vector<Eigen::VectorXd> Dir = {
      PD1.row(v1).normalized(), PD2.row(v1).normalized(),
      -PD1.row(v1).normalized(), -PD2.row(v1).normalized()};

  // store max alignment -> pair (max_aligment, (edge_index, if_flipped))
  std::vector<std::pair<double, std::pair<int, bool>>> max_alignment = {
      std::make_pair(-1, std::make_pair(-1, false)),
      std::make_pair(-1, std::make_pair(-1, false)),
      std::make_pair(-1, std::make_pair(-1, false)),
      std::make_pair(-1, std::make_pair(-1, false))};

  //// first get the edges of each face in V2Fe
  for (int f : V2Fe) {
    if (f == f1 || f == f2 || f >= F.rows())
      continue;

    for (int i = 0; i < 3; i++) {
      int ei = EMAP(f * 3 + i); // Need to Make sure this indexing is right
      if (E(ei, 0) == v1 ||
          E(ei, 1) == v1) // the edge is already there, no flip to generate it
      {
        // just making vs always v1 and ve the other vertex
        int vs = v1;
        int ve = E(ei, 1);
        if (E(ei, 1) == v1)
          ve = E(ei, 0);

        Eigen::Vector3d vs_e_vec = V.row(ve) - V.row(vs);
        vs_e_vec.normalize();

        for (int j = 0; j < 4; j++) {
          double alignment = vs_e_vec.dot(Dir[j]);
          if (alignment > max_alignment[j].first) {
            max_alignment[j].first = alignment;
            max_alignment[j].second.first = ei;
            max_alignment[j].second.second = false;
          }
        }
      } else // the edge is not there, flip to generate it
      {
        int vs = v1;
        // ve would be the on the other side of the edge, first get the face:
        int f_opposite = EF(ei, 0);
        if (f_opposite == f)
          f_opposite = EF(ei, 0);

        if (f_opposite >= F.rows())
          continue;

        // now get the vertex that is not E(ei, 0) or E(ei, 1)
        int ve = F(f_opposite, 0);
        if (ve == E(ei, 0) || ve == E(ei, 1))
          ve = F(f_opposite, 1);
        if (ve == E(ei, 0) || ve == E(ei, 1))
          ve = F(f_opposite, 2);

        Eigen::Vector3d vs_e_vec = V.row(ve) - V.row(vs);
        vs_e_vec.normalize();

        for (int j = 0; j < 4; j++) {
          double alignment = vs_e_vec.dot(Dir[j]);
          if (alignment > max_alignment[j].first) {
            max_alignment[j].first = alignment;
            max_alignment[j].second.first = ei;
            max_alignment[j].second.second = true;
          }
        }
      }
    }
  }
  return true;
};

igl::decimate_stopping_condition_callback stopping_condition_callback;
