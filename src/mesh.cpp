#include "../include/mesh.h"

void Mesh::initialize_decimation_callbacks()
{

  custom_pre_collapse_callback = [&](
                                     const MatrixXd &V,
                                     const MatrixXi &F,
                                     const MatrixXi &E,
                                     const VectorXi &EMAP,
                                     const MatrixXi &EF,
                                     const MatrixXi &EI,
                                     const igl::min_heap<std::tuple<double, int, int>> &Q,
                                     const VectorXi &EQ,
                                     const MatrixXd &C,
                                     const int e) -> bool
  {
    V2Fe = igl::circulation(e, true, EMAP, EF, EI);
    return true;
  };

  custom_post_collapse_callback = [&](
                                      const MatrixXd &V,
                                      const MatrixXi &F,
                                      const MatrixXi &E,
                                      const VectorXi &EMAP,
                                      const MatrixXi &EF,
                                      const MatrixXi &EI,
                                      const igl::min_heap<std::tuple<double, int, int>> &Q,
                                      const VectorXi &EQ,
                                      const MatrixXd &C,
                                      const int e,
                                      const int e1,
                                      const int e2,
                                      const int f1,
                                      const int f2,
                                      const bool collapsed) -> bool
  {
    if (!collapsed)
      return false;

    int v1 = E(e, 0);
    int v2 = E(e, 1);

    PD1.row(v1) = 0.5 * (PD1.row(v1) + PD1.row(v2));
    PV1.row(v1) = 0.5 * (PV1.row(v1) + PV1.row(v2));
    PD2.row(v1) = 0.5 * (PD2.row(v1) + PD2.row(v2));
    PV2.row(v1) = 0.5 * (PV2.row(v1) + PV2.row(v2));

    // 4 normalized vectors: PD1.row(v1), PD2.row(v1), - PD1.row(v1), - PD2.row(v1)
    std::vector<VectorXd> Dir =
        {
            PD1.row(v1).normalized(),
            PD2.row(v1).normalized(),
            -PD1.row(v1).normalized(),
            -PD2.row(v1).normalized()};

    // store max alignment -> pair (max_aligment, (edge_index, if_flipped))
    std::vector<std::pair<double, std::pair<int, bool>>> max_alignment = {
        std::make_pair(-1, std::make_pair(-1, false)),
        std::make_pair(-1, std::make_pair(-1, false)),
        std::make_pair(-1, std::make_pair(-1, false)),
        std::make_pair(-1, std::make_pair(-1, false))};

    //// first get the edges of each face in V2Fe
    for (int f : V2Fe)
    {
      if (f == f1 || f == f2 || f >= F.rows())
        continue;

      for (int i = 0; i < 3; i++)
      {
        int ei = EMAP(f * 3 + i);             // Need to Make sure this indexing is right
        if (E(ei, 0) == v1 || E(ei, 1) == v1) // the edge is already there, no flip to generate it
        {
          // just making vs always v1 and ve the other vertex
          int vs = v1;
          int ve = E(ei, 1);
          if (E(ei, 1) == v1)
            ve = E(ei, 0);

          Vector3d vs_e_vec = V.row(ve) - V.row(vs);
          vs_e_vec.normalize();

          for (int j = 0; j < 4; j++)
          {
            double alignment = vs_e_vec.dot(Dir[j]);
            if (alignment > max_alignment[j].first)
            {
              max_alignment[j].first = alignment;
              max_alignment[j].second.first = ei;
              max_alignment[j].second.second = false;
            }
          }
        }
        else // the edge is not there, flip to generate it
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

          Vector3d vs_e_vec = V.row(ve) - V.row(vs);
          vs_e_vec.normalize();

          for (int j = 0; j < 4; j++)
          {
            double alignment = vs_e_vec.dot(Dir[j]);
            if (alignment > max_alignment[j].first)
            {
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
}

void Mesh::decimate()
{
  std::cout << "Decimating into 3/4 of faces." << std::endl;

  std::cout << "Face and Vert Size: " << F.rows() << " " << V.rows() << std::endl;
  std::cout << "PD PV Size: " << PD1.rows() << " " << PD2.rows() << " " << PV1.rows() << " " << PV2.rows() << std::endl;

  int curr_nF, orig_nF, target_nF;
  curr_nF = F.rows();
  orig_nF = F.rows();
  target_nF = 3 * F.rows() / 4;

  igl::max_faces_stopping_condition(curr_nF, orig_nF, target_nF, stopping_condition_callback);

  VectorXi J, I;
  igl::decimate(
      V,
      F,
      igl::shortest_edge_and_midpoint,
      stopping_condition_callback,
      custom_pre_collapse_callback,
      custom_post_collapse_callback,
      V,
      F,
      J,
      I);
}

// Below 2 functions came from Frame_field.cpp
void Mesh::get_frame_field()
{
  igl::principal_curvature(V, F, PD1, PD2, PV1, PV2);
}

void Mesh::frame_field_alignment_data()
{
  get_frame_field();

  // Initialize resulting data.
  A = VectorXd(V.rows());
  for (int i = 0; i < V.rows(); ++i)
    A[i] = 1;

  for (int i = 0; i < F.rows(); ++i)
  {
    RowVector3i face = F.row(i);

    for (int j = 0; j < 3; ++j)
    {
      int adj1 = face[(i + 1) % 3];
      int adj2 = face[(i + 2) % 3];
      RowVector3d P = V.row(adj1);
      RowVector3d Q = V.row(adj2);
      RowVector3d p = Q - P;

      RowVector3d v1 = PD1.row(adj1);
      RowVector3d w1 = PD2.row(adj1);
      RowVector3d v2 = PD1.row(adj2);
      RowVector3d w2 = PD2.row(adj2);

      double f1 = alignment_function(v1, w1, p);
      double f2 = alignment_function(v2, w2, p);

      A[adj1] = f1 < A[adj1] ? f1 : A[adj1];
      A[adj2] = f2 < A[adj2] ? f2 : A[adj2];
    }
  }
}
