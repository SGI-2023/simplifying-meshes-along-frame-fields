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
    // The circulation functino only returns faces around one vertex of the edge (circulates the edge ccw or cw)
    // Thus, we get both and then insert them into a set to get a unique list of faces
    V2Fe.clear();
    std::vector<int> V2FeVec = igl::circulation(e, true, EMAP, EF, EI);
    V2Fe.insert(V2FeVec.begin(), V2FeVec.end());
    V2FeVec = igl::circulation(e, false, EMAP, EF, EI);
    V2Fe.insert(V2FeVec.begin(), V2FeVec.end());

    v1 = E(e, 0);
    v2 = E(e, 1);

    return true;
  };

  // You can change this cost and placement functions with yours
  cost_and_placement_callback = [&](
                                    const int e,
                                    const Eigen::MatrixXd &V,
                                    const Eigen::MatrixXi &F,
                                    const Eigen::MatrixXi &E,
                                    const Eigen::VectorXi &EMAP,
                                    const Eigen::MatrixXi &EF,
                                    const Eigen::MatrixXi &EI,
                                    double &cost,
                                    Eigen::RowVectorXd &p)
  {
    // cost = frame_field_alignment_data(e, E); -> another cost function based on only frame-field alignment
    p = 0.5 * (V.row(E(e, 0)) + V.row(E(e, 1)));
    qcoarsen_based_cost(e, V, E, F, EMAP, EF, EI, p, cost); 
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

    // update the frame field by averaging the frame field of the two vertices that are collapsed
    PD1.row(v2) = PD1.row(v1) = 0.5 * (PD1.row(v1) + PD1.row(v2));
    PV1.row(v2) = PV1.row(v1) = 0.5 * (PV1.row(v1) + PV1.row(v2));
    PD2.row(v2) = PD2.row(v1) = 0.5 * (PD2.row(v1) + PD2.row(v2));
    PV2.row(v2) = PV2.row(v1) = 0.5 * (PV2.row(v1) + PV2.row(v2));
    // Note, we are updating both as it seems there's no convention for which we know which vertex is the one kept


    // dir should be normalized
    PD1.row(v1).normalize();
    PD2.row(v1).normalize();
    PD1.row(v2).normalize();
    PD2.row(v2).normalize();

    // 4 normalized vectors: PD1.row(v1), PD2.row(v1), - PD1.row(v1), - PD2.row(v1)
    std::vector<VectorXd> Dir =
        {
            PD1.row(v1),
            PD2.row(v1),
            -PD1.row(v1),
            -PD2.row(v1)};

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

      int vs = 0;
      if (F(f, 0) == v1 || F(f, 1) == v1 || F(f, 2) == v1)
        vs = v1;
      else if (F(f, 0) == v2 || F(f, 1) == v2 || F(f, 2) == v2)
        vs = v2;
      else
        continue;

      for (int i = 0; i < 3; i++)
      {
        int ei = EMAP(i * F.rows() + f);             // Need to Make sure this indexing is right

        // the edge is already there, no flip to generate it
        if (E(ei, 0) == vs || E(ei, 1) == vs)
        {
          // just making vs always v1 and ve the other vertex
          int ve = E(ei, 1);
          if (E(ei, 1) == vs)
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
      
    /*std::cout << "direction 1: " << max_alignment[0].first << ", "
      << max_alignment[0].second.first << ", " << max_alignment[0].second.second << std::endl;
    std::cout << "direction 2: " << max_alignment[1].first << ", "
      << max_alignment[1].second.first << ", " << max_alignment[1].second.second << std::endl;
    std::cout << "direction 3: " << max_alignment[2].first << ", "
      << max_alignment[2].second.first << ", " << max_alignment[2].second.second << std::endl;
    std::cout << "direction 4: " << max_alignment[3].first << ", "
      << max_alignment[3].second.first << ", " << max_alignment[3].second.second << std::endl;*/

    //typedef Eigen::Index int;
    typedef Eigen::Matrix<int,Eigen::Dynamic,2> MatrixX2I;
    MatrixX2I En,uEn;
    Eigen::VectorXi EMAPn;
    std::vector<std::vector<int> > uE2En;
    Eigen::Matrix<int, Eigen::Dynamic, 3> Fn(F);
    igl::unique_edge_map(F, En, uEn, EMAPn, uE2En);

    // now we have the best alignment for each direction, we can flip the edges accordingly
    /// The code below crashes after some iterations (index out of range), not yet sure why tho.
    //for (int j = 0; j < 4; j++)
    //{
    //  if (max_alignment[j].second.first != -1)
    //  {
    //  	int ei = max_alignment[j].second.first;
    //  	bool flip = max_alignment[j].second.second;
    //  	if (flip)
    //  	{
    //  	  // use the edge flip function
    //  	  igl::flip_edge(Fn, En, uEn, EMAPn, uE2En,  ei);
    //  	}
    //  }
    //}

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
  target_nF = F.rows() * 3/4;

  igl::max_faces_stopping_condition(curr_nF, orig_nF, target_nF, stopping_condition_callback);

  VectorXi J, I;
  igl::decimate(
      V,
      F,
      cost_and_placement_callback,
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

// vertex-based
void Mesh::frame_field_alignment_data()
{
  get_frame_field();

  // Initialize resulting data.
  A = VectorXd(V.rows());
  for (int i = 0; i < V.rows(); ++i)
    A[i] = 0;

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

      A[adj1] = f1 > A[adj1] ? f1 : A[adj1];
      A[adj2] = f2 > A[adj2] ? f2 : A[adj2];
    }
  }
}

// edge-based
double Mesh::frame_field_alignment_data(const int e, const MatrixXi &E)
{
  // Initialize resulting data for edges.
  // Assuming the edge matrix is called 'edges' with dimensions E x 2.
  //std::cout << E.rows() << std::endl;

  int adj1 = E(e, 0);
  int adj2 = E(e, 1);

  RowVector3d P = V.row(adj1);
  RowVector3d Q = V.row(adj2);
  RowVector3d p = Q - P;

  RowVector3d v1 = PD1.row(adj1);
  RowVector3d w1 = PD2.row(adj1);
  RowVector3d v2 = PD1.row(adj2);
  RowVector3d w2 = PD2.row(adj2);

  double f1 = alignment_function(v1, w1, p);
  double f2 = alignment_function(v2, w2, p);

  return std::max(f1, f2);
}

// qcoarsen_based_cost and its helper functions 

void Mesh::compute_before(const int e,
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi & E,
	const Eigen::MatrixXi & F,
	const Eigen::VectorXi & EMAP,
	const Eigen::MatrixXi & EF,
	const Eigen::MatrixXi & EI,
	int & valenceSum,
  	double & alignment,
	std::vector<int> & neighbors)
{
		std::vector<int> Nv;
  	std::vector<int> Nf;
		igl::circulation(e, true, F,EMAP,EF,EI, Nv,Nf);
		valenceSum = Nv.size();

		double alignmentTotal = 0;
		for (int i = 0; i < valenceSum; i++) {
			if (i != 0) {
				neighbors.push_back(Nv[i]);
			}
			alignmentTotal += alignment_function(PD1.row(Nv[i]), PD2.row(Nv[i]),V.row(Nv[i])-V.row(E.col(0)(e)));
		}
		igl::circulation(e, false, F,EMAP,EF,EI, Nv,Nf);

		for (int i = 1; i < Nv.size(); i++) {
			neighbors.push_back(Nv[i]);
			alignmentTotal += alignment_function(PD1.row(Nv[i]), PD2.row(Nv[i]),V.row(Nv[i])-V.row(E.col(1)(e)));
    }
    alignment = alignmentTotal/(valenceSum+Nv.size()-1);
    valenceSum = valenceSum == 6 ? 1 : 0;
		valenceSum += Nv.size() == 6 ? 1 : 0;
}
void Mesh::compute_after(
	const Eigen::MatrixXd& V,
	const std::vector<int> & neighbors,//as generated by compute_before
  const Eigen::RowVector3d & p, //new point
	int & valenceSum,
  double & alignment)
{
		double alignmentTotal = 0;
		for (int i = 0; i < neighbors.size(); i++) {
			alignmentTotal += alignment_function(PD1.row(neighbors[i]), PD2.row(neighbors[i]),V.row(neighbors[i])-p);
    }
		alignment = alignmentTotal/neighbors.size();
		valenceSum = neighbors.size() == 6 ? 1 : 0;
}
void Mesh::qcoarsen_based_cost(
  const int e,
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi & E,
	const Eigen::MatrixXi & F,
	const Eigen::VectorXi & EMAP,
	const Eigen::MatrixXi & EF,
	const Eigen::MatrixXi & EI,
  const Eigen::RowVector3d & p,
  double & cost
	)
{
  int idealCountB, idealCountA;
  double alignB, alignA;
  std::vector<int> neighbors;
  compute_before(e,V,E,F,EMAP,EF,EI,idealCountB,alignB,neighbors);
  compute_after(V, neighbors, p, idealCountA, alignA);
  double alpha = 0.01;
  double total_ideal = V.rows();//for now, as a placeholder
  double vB_A = (V.cols()+1)*total_ideal/(V.cols()*(total_ideal+idealCountA-idealCountB));
  double dist = (V.row(E(e,0))-V.row(E(e,1))).norm();
  double weight = 1;//placeholder
  cost = (vB_A + alpha)*(vB_A + alpha)*(alignB/alignA+alpha)*(alignB/alignA+alpha)*(dist*weight+alpha)*(dist*weight+alpha);
}
