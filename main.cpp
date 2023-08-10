#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/decimate.h>
#include <igl/principal_curvature.h>
#include <igl/max_faces_stopping_condition.h>
#include <igl/decimate_trivial_callbacks.h>
#include <igl/flip_edge.h>
#include <igl/circulation.h>
#include <Eigen/Core>
//#include "collapsed_dec.h"
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
  std::string filename = "../data/spot.obj";

  Eigen::MatrixXd V, U;
  Eigen::MatrixXi F, G;
  Eigen::VectorXi J, I;

  igl::read_triangle_mesh(filename, V, F);

  // Compute curvature directions via quadric fitting
  Eigen::MatrixXd PD1, PD2; // V * 3      Where V is the number of the Vertices
  Eigen::VectorXd PV1, PV2; // V          Where V is the number of the Vertices
  igl::principal_curvature(V, F, PD1, PD2, PV1, PV2);

  /// In case we needed a more customized cost & placement call back, we can uncomment this

  //const igl::decimate_cost_and_placement_callback custom_cost_and_placement_callback = [=](
  //  const int e,
  //  const Eigen::MatrixXd& V,
  //  const Eigen::MatrixXi& F,
  //  const Eigen::MatrixXi& E,
  //  const Eigen::VectorXi& EMAP,
  //  const Eigen::MatrixXi& EF,
  //  const Eigen::MatrixXi& EI,
  //  double& cost,
  //  Eigen::RowVectorXd& p
  //  ) -> void
  //{
  //  ;
  //};

  /// In case we needed a more customized stopping condition, we can uncomment this

  //const igl::decimate_stopping_condition_callback custom_stopping_condition_callback = [=](
  //  const Eigen::MatrixXd& V,
  //  const Eigen::MatrixXi& F,
  //  const Eigen::MatrixXi& E,
  //  const Eigen::VectorXi& EMAP,
  //  const Eigen::MatrixXi& EF,
  //  const Eigen::MatrixXi& EI,
  //  const igl::min_heap<std::tuple<double, int, int>>& Q,
  //  const Eigen::VectorXi& EQ,
  //  const Eigen::MatrixXd& C,
  //  const int e,
  //  const int e1,
  //  const int e2,
  //  const int f1,
  //  const int f2) -> bool
  //{
  //  return Q.size() < 8000; // TODO: clean up
  //};

  std::vector<int> V2Fe;

  const igl::decimate_pre_collapse_callback custom_pre_collapse_callback = [&V2Fe](
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXi& E,
    const Eigen::VectorXi& EMAP,
    const Eigen::MatrixXi& EF,
    const Eigen::MatrixXi& EI,
    const igl::min_heap<std::tuple<double, int, int>>& Q,
    const Eigen::VectorXi& EQ,
    const Eigen::MatrixXd& C,
    const int e) -> bool
  {
    V2Fe = igl::circulation(e, true, EMAP, EF, EI);
    return true;
  };

  const igl::decimate_post_collapse_callback custom_post_collapse_callback = [&PD1, &PD2, &PV1, &PV2, &V2Fe](
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXi& E,
    const Eigen::VectorXi& EMAP,
    const Eigen::MatrixXi& EF,
    const Eigen::MatrixXi& EI,
    const igl::min_heap<std::tuple<double, int, int>>& Q,
    const Eigen::VectorXi& EQ,
    const Eigen::MatrixXd& C,
    const int e,
    const int e1,
    const int e2,
    const int f1,
    const int f2,
    const bool collapsed) -> bool
  {
    //   EMAP #F*3 list of indices into E, mapping each directed edge to unique
    //     edge in E
    //   EF  #E by 2 list of edge flaps, EF(e,0)=f means e=(i-->j) is the edge of
    //     F(f,:) opposite the vth corner, where EI(e,0)=v. Similarly EF(e,1) "
    //     e=(j->i)
    //   EI  #E by 2 list of edge flap corners (see above).
    //   E   #E by 2 list of edge indices, each row containing the indices of the
    //       two vertices that the edge connects

    if (!collapsed)
      return false;

    // Assumptions (to be cleaned/handled later)
    // - v1 is the vertex the edge is collapsed into
    // - all are manifold
    // - no boundary
    // EMAP (f, i) = EMAP (f * 3 + i); 

    int v1 = E(e, 0);
    int v2 = E(e, 1);

    PD1.row(v1) = 0.5 * (PD1.row(v1) + PD1.row(v2));
    PV1.row(v1) = 0.5 * (PV1.row(v1) + PV1.row(v2));
    PD2.row(v1) = 0.5 * (PD2.row(v1) + PD2.row(v2));
    PV2.row(v1) = 0.5 * (PV2.row(v1) + PV2.row(v2));

    // 4 normalized vectors: PD1.row(v1), PD2.row(v1), - PD1.row(v1), - PD2.row(v1)
    std::vector<Eigen::VectorXd> Dir =
    {
        PD1.row(v1).normalized(),
        PD2.row(v1).normalized(),
        -PD1.row(v1).normalized(),
        -PD2.row(v1).normalized() 
    };

    // store max alignment -> pair (max_aligment, (edge_index, if_flipped))
    std::vector<std::pair<double, std::pair<int, bool>>> max_alignment = {
        std::make_pair(-1, std::make_pair(-1, false)),
        std::make_pair(-1, std::make_pair(-1, false)),
        std::make_pair(-1, std::make_pair(-1, false)),
        std::make_pair(-1, std::make_pair(-1, false)) };

    //// first get the edges of each face in V2Fe
    for (int f : V2Fe)
    {
      if (f == f1 || f == f2 || f >= F.rows())
        continue;

      for (int i = 0; i < 3; i++)
      {
        int ei = EMAP(f * 3 + i); // Need to Make sure this indexing is right
        if (E(ei, 0) == v1 || E(ei, 1) == v1) // the edge is already there, no flip to generate it
        {
          // just making vs always v1 and ve the other vertex
          int vs = v1;
          int ve = E(ei, 1);
          if (E(ei, 1) == v1)
            ve = E(ei, 0);

          Eigen::Vector3d vs_e_vec = V.row(ve) - V.row(vs);
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

          Eigen::Vector3d vs_e_vec = V.row(ve) - V.row(vs);
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

    //std::cout << "Direction 1: " << max_alignment[0].first << ", "
    //  << max_alignment[0].second.first << ", " << max_alignment[0].second.second << std::endl;
    //std::cout << "Direction 2: " << max_alignment[1].first << ", "
    //  << max_alignment[1].second.first << ", " << max_alignment[1].second.second << std::endl;
    //std::cout << "Direction 3: " << max_alignment[2].first << ", "
    //  << max_alignment[2].second.first << ", " << max_alignment[2].second.second << std::endl;
    //std::cout << "Direction 4: " << max_alignment[3].first << ", "
    //  << max_alignment[3].second.first << ", " << max_alignment[3].second.second << std::endl;

    // typedef Eigen::MatrixXi::Scalar Index;
    // typedef Eigen::Matrix<Index,Eigen::Dynamic,2> MatrixX2I;
    // MatrixX2I En,uEn;
    // Eigen::VectorXi EMAPn;
    // std::vector<std::vector<Index> > uE2En;
    // igl::unique_edge_map(F, En, uEn, EMAPn, uE2En);


    // // now we have the best alignment for each direction, we can flip the edges accordingly
    // for (int j = 0; j < 4; j++)
    // {
    // 	if (max_alignment[j].second.first != -1)
    // 	{
    // 		int ei = max_alignment[j].second.first;
    // 		bool flip = max_alignment[j].second.second;
    // 		if (flip)
    // 		{
    // 			// use the edge flip function
    // 			igl::flip_edge(F, En, uEn, EMAPn, uE2En,  ei);
    // 		}
    // 	}
    // }

    return true;
  };

  std::cout << "Decimating into 1/10 of faces." << std::endl;

  int curr_nF = F.rows();
  int orig_nF = F.rows();
  int target_nF = F.rows() / 10;

  std::cout << F.rows() << std::endl;

  igl::decimate_stopping_condition_callback stopping_condition_callback;

  igl::max_faces_stopping_condition(curr_nF, orig_nF, target_nF, stopping_condition_callback);

  // igl::decimate(V, F, num_faces, U, G, J, I);
  igl::decimate(
    V,
    F,
    igl::shortest_edge_and_midpoint,
    stopping_condition_callback,
    custom_pre_collapse_callback,
    custom_post_collapse_callback,
    U,
    G,
    J,
    I);

  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(U, G);
  viewer.launch();

  return 0;
}
