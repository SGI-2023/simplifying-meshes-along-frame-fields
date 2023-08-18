#pragma once

#include <Eigen/Core>
#include <tuple>
#include <vector>

// Compute quadrics per vertex of a "closed" triangle mesh (V,F). Rather than
// follow the qslim paper, this implements the lesser-known _follow up_
// "Simplifying Surfaces with Color and Texture using Quadric Error Metrics".
// This allows V to be n-dimensional (where the extra coordiantes store
// texture UVs, color RGBs, etc.
//
// Inputs:
//   V  #V by n list of vertex positions. Assumes that vertices with
//     infinite coordinates are "points at infinity" being used to close up
//     boundary edges with faces. This allows special subspace quadrice for
//     boundary edges: There should never be more than one "point at
//     infinity" in a single triangle.
//   F  #F by 3 list of triangle indices into V
//   E  #E by 2 list of edge indices into V.
//   EMAP #F*3 list of indices into E, mapping each directed edge to unique
//     unique edge in E
//   EF  #E by 2 list of edge flaps, EF(e,0)=f means e=(i-->j) is the edge of
//     F(f,:) opposite the vth corner, where EI(e,0)=v. Similarly EF(e,1) "
//     e=(j->i)
//   EI  #E by 2 list of edge flap corners (see above).
// Outputs:
//   quadrics  #V list of quadrics, where a quadric is a tuple {A,b,c} such
//     that the quadratic energy of moving this vertex to position x is
//     given by x'Ax - 2b + c
//
void per_vertex_point_to_plane_quadrics_aligned(
    const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
    const Eigen::MatrixXi &EMAP, const Eigen::MatrixXi &EF,
    const Eigen::MatrixXi &EI, const Eigen::MatrixXd &PD1,
    const Eigen::MatrixXd &PD2, const Eigen::VectorXd &PV1,
    const Eigen::VectorXd &PV2,
    std::vector<std::tuple<Eigen::MatrixXd, Eigen::RowVectorXd, double>>
        &quadrics);
