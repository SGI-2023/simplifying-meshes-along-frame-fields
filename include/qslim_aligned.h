#pragma once

#include <Eigen/Core>

/// Decimate (simplify) a triangle mesh in nD according to the paper
/// "Simplifying Surfaces with Color and Texture using Quadric Error Metrics"
/// by [Garland and Heckbert, 1987] (technically a followup to qslim). The
/// mesh can have open boundaries but should be edge-manifold.
///
/// @param[in] V  #V by dim list of vertex positions. Assumes that vertices with
///     infinite coordinates are "points at infinity" being used to close up
///     boundary edges with faces. This allows special subspace quadrice for
///     boundary edges: There should never be more than one "point at
///     infinity" in a single triangle.
/// @param[in] F  #F by 3 list of triangle indices into V
/// @param[in] max_m  desired number of output faces
/// @param[out] U  #U by dim list of output vertex posistions (can be same ref
/// as V)
/// @param[out] G  #G by 3 list of output face indices into U (can be same ref
/// as F)
/// @param[out] J  #G list of indices into F of birth face
/// @param[out] I  #U list of indices into V of birth vertices
bool qslim_aligned(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                   const size_t max_m, Eigen::MatrixXd &U, Eigen::MatrixXi &G,
                   Eigen::VectorXi &J, Eigen::VectorXi &I);
