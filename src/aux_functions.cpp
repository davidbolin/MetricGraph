#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
//

//' @name assemble_fem
//' @title Construction of FEM matrices
//' @description Function used to construct FEM matrices on metric graphs.
//' @param E [nx2 matrix] Matrix of edges
//' @param h_e [n vector] Vector of h's
//' @param nV [int] Number of vertices
//' @noRd
//'
// [[Rcpp::export]]

Rcpp::List assemble_fem(Eigen::MatrixXd E, Eigen::VectorXd h_e, int nV, bool petrov){

    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> trp_C, trp_G, trp_B;
    int i;
    int v1,v2;
    int nE  = E.rows();
    Eigen::SparseMatrix<double> C(nV,nV), G(nV,nV), B(nV,nV);

    std::vector<Trip> trp_Cpet, trp_Gpet;
    Eigen::SparseMatrix<double> Cpet(nV,nE), Gpet(nV,nE);

    for(i=0; i<nE; i++){
        v1 = E(i, 0)-1;
        v2 = E(i, 1)-1;

        // Assembling C
        trp_C.push_back(Trip(v1,v1,h_e(i)/3));
        trp_C.push_back(Trip(v2,v2,h_e(i)/3));
        trp_C.push_back(Trip(v1,v2,h_e(i)/6));
        trp_C.push_back(Trip(v2,v1,h_e(i)/6));

        // Assembling G
        trp_G.push_back(Trip(v1,v1,1/h_e(i)));
        trp_G.push_back(Trip(v2,v2,1/h_e(i)));
        trp_G.push_back(Trip(v1,v2,-1/h_e(i)));
        trp_G.push_back(Trip(v2,v1,-1/h_e(i)));

        // Assembling B
        trp_B.push_back(Trip(v1,v1,-0.5));
        trp_B.push_back(Trip(v1,v2,-0.5));
        trp_B.push_back(Trip(v2,v2,0.5));
        trp_B.push_back(Trip(v2,v1,0.5));

        if(petrov){
          // Assembling Cpet
          trp_Cpet.push_back(Trip(v1,i,h_e(i)/2));
          trp_Cpet.push_back(Trip(v2,i,h_e(i)/2));

          // Assembling Gpet
          trp_Gpet.push_back(Trip(v1,i,-1));
          trp_Gpet.push_back(Trip(v2,i,1));
        }
    }

    C.setFromTriplets(trp_C.begin(), trp_C.end());
    G.setFromTriplets(trp_G.begin(), trp_G.end());
    B.setFromTriplets(trp_B.begin(), trp_B.end());

    Rcpp::List out;
    out["C"] = C;
    out["G"] = G;
    out["B"] = B;

    if(petrov){
      Cpet.setFromTriplets(trp_Cpet.begin(), trp_Cpet.end());
      Gpet.setFromTriplets(trp_Gpet.begin(), trp_Gpet.end());
      out["Cpet"] = Cpet;
      out["Gpet"] = Gpet;
    }

    return(out);
}

// Obtain the coordinates of the projection of a point in a line. line = p0 + vt, point p
//' @name proj_vec
//' @noRd
Eigen::Vector2d proj_vec(Eigen::VectorXd p0, Eigen::VectorXd v, Eigen::VectorXd p){
  return p0 + (p-p0).dot(v) * v/(v.dot(v));
}

// Obtain the coordinates of the point in the line that is closest to the point. line = p0 + vt, point p
//' @name proj_vec2
//' @noRd
Eigen::Vector2d proj_vec2(Eigen::VectorXd p0, Eigen::VectorXd v, Eigen::VectorXd p){
   Eigen::VectorXd proj = p0 + (p-p0).dot(v) * v/(v.dot(v));
   if((p-p0).dot(v)/(v.dot(v)) > 1){
        proj = p0+v;
   } else if((p-p0).dot(v)/(v.dot(v)) < 0){
        proj = p0;
   }
   return proj;
}

// Obtain the distance along the line between the projected point
//' @name proj_vec
//' @noRd
double proj_vec_dist(Eigen::MatrixXd line, Eigen::VectorXd point){
  Eigen::VectorXd p0 = line.row(0);
  Eigen::VectorXd v = line.row(1) - line.row(0);
  Eigen::VectorXd proj = proj_vec(p0, v, point);
  if((point-p0).dot(v)/(v.dot(v)) > 1){
    proj = line.row(1);
  } else if((point-p0).dot(v)/(v.dot(v)) < 0){
    proj = line.row(0);
  }
  return (proj-p0).norm();
}

// Obtain the distance of a point along a linestring.
//' @name proj_vec_line
//' @noRd

double proj_vec_line(Eigen::MatrixXd line, Eigen::VectorXd point, int normalized = 0){
  Eigen::VectorXd dist_vec(line.rows());
  int i, min_index;
  double min_dist = -1.0;
  dist_vec(0) = 0;
  for(i=0; i<line.rows()-1; i++){
    Eigen::VectorXd p0 = line.row(i);
    Eigen::VectorXd v = line.row(i+1) - line.row(i);
    Eigen::VectorXd proj_temp = proj_vec2(p0, v, point);
    if(min_dist < 0){
        min_dist = (point - proj_temp).norm();
        min_index = i;
    }
    double dist_temp = (point - proj_temp).norm();
    if(dist_temp < min_dist){
        min_index = i;
        min_dist = dist_temp;
    }
    dist_vec(i+1) = dist_vec(i) + v.norm();
  }

  Eigen::MatrixXd line_temp = line.block(min_index,0,2,2);
  double proj_dist = proj_vec_dist(line_temp, point);
  double dist_return = dist_vec(min_index) + proj_dist;
  if(dist_return > dist_vec(min_index+1)){
    dist_return = dist_vec(min_index+1);
  }
  if(dist_return < 0){
    dist_return = 0;
  }
  if(normalized != 0){
    dist_return = dist_return/dist_vec(line.rows()-1);
  }
  return(dist_return);
}

//' @name projectVecLine
//' @title Projects SpatialPoints into SpatialLines
//' @description Obtain the coordinates of the projection of points into lines.
//' @param lines [nx2 matrix] Matrix of the points of the lines
//' @param points [nx2 matrix] Matrix of the points
//' @param normalized [int] 0 means not normalized, 1 means normalized
//' @noRd
//'
// [[Rcpp::export]]

Eigen::VectorXd projectVecLine(Eigen::MatrixXd lines, Eigen::MatrixXd points, int normalized = 0){
    int size_return = points.rows();
    int i;
    Eigen::VectorXd out_vec(size_return);
    for(i = 0; i < size_return; i++){
        Eigen::VectorXd point = points.row(i);
        out_vec(i) = proj_vec_line(lines, point, normalized);
    }
    return(out_vec);
}

//' @name interpolate2
//' @title Finds the point with respect to a distance along the line
//' @description Finds the point with respect to a distance along the line
//' @param lines [nx2 matrix] Matrix of the points of the lines
//' @param pos [k vector] vector of positions.
//' @param normalized [int] 0 means not normalized, 1 means normalized
//' @noRd
//'
// [[Rcpp::export]]

Rcpp::List interpolate2_aux(Eigen::MatrixXd lines, Eigen::VectorXd pos, int normalized = 0){
    int size_return = pos.size();
    int i,j;
    Eigen::MatrixXd out_mat(size_return,2);
    Eigen::VectorXd dist_vec(lines.rows());
    Eigen::VectorXd idx_pos(pos.size());
    dist_vec(0) = 0;
    for(i=0; i<lines.rows()-1; i++){
        Eigen::VectorXd p0 = lines.row(i);
        Eigen::VectorXd v = lines.row(i+1) - lines.row(i);
        dist_vec(i+1) = dist_vec(i) + v.norm();
    }
    dist_vec = dist_vec/dist_vec(lines.rows()-1);
    Eigen::VectorXd pos_rel;
    if(normalized != 0){
        pos_rel = pos;
    } else{
        pos_rel = pos/dist_vec(lines.rows()-1);
    }

    for(i=0; i< pos.size(); i++){
        int tmp_ind = -1;
        if(pos_rel(i) < 0){
            pos_rel(i) = 0;
        } else if(pos_rel(i)>1){
            pos_rel(i) = 1;
        }
        for(j=0; j<dist_vec.size()-1; j++){
            if(pos_rel(i) >= dist_vec(j) && pos_rel(i) <= dist_vec(j+1)){
                tmp_ind = j;
            }
        }

        double dist_pos = (pos_rel(i) - dist_vec(tmp_ind))/(dist_vec(tmp_ind+1)-dist_vec(tmp_ind));

        out_mat.row(i) = lines.row(tmp_ind) + (lines.row(tmp_ind+1) - lines.row(tmp_ind))*dist_pos;
        idx_pos(i) = tmp_ind+1;
    }

    return  Rcpp::List::create(
      Rcpp::Named("coords")        = out_mat,
      Rcpp::Named("idx") = idx_pos
    );
}


//' @name compute_length
//' @title Compute length
//' @description Computes the length of a piecewise-linear function whose coordinates are given in a matrix.
//' @param coords nx2 matrix Matrix of the points of the lines
//' @noRd
//'
// [[Rcpp::export]]
double compute_length(Eigen::MatrixXd coords) {

    double arclength = 0;

    int i;

    for(i = 0 ; i < coords.rows()-1; i++){
         Eigen::VectorXd v = coords.row(i+1) - coords.row(i);
         arclength = arclength + v.norm();
     }

    return(arclength);
}


// [[Rcpp::export]]
List generate_mesh(int n_edges, NumericVector edge_lengths, double h, IntegerVector n_e, 
                   IntegerMatrix E, IntegerVector ind, bool continuous) {
  std::vector<int> PtE_edge;
  std::vector<double> PtE_pos;
  std::vector<double> h_e;
  std::vector<int> E_start, E_end;
  int current_max_index = ind.size(); // Start from the existing number of vertices

  for (int i = 0; i < n_edges; ++i) {
    if (n_e[i] > 0) {
      // Generate d.e sequence
      std::vector<double> d_e(n_e[i]);
      for (int j = 0; j < n_e[i]; ++j) {
        d_e[j] = (j + 1) / static_cast<double>(n_e[i] + 1);
        PtE_edge.push_back(i + 1); // R uses 1-based indexing
        PtE_pos.push_back(d_e[j]);
      }

      // Compute h_e based on the first d_e value
      double segment_length = edge_lengths[i] * d_e[0];
      for (int j = 0; j < n_e[i] + 1; ++j) {
        h_e.push_back(segment_length);
      }

      // Create internal vertex indices
      std::vector<int> V_int(n_e[i]);
      for (int j = 0; j < n_e[i]; ++j) {
        V_int[j] = current_max_index + 1 + j;
      }
      current_max_index += n_e[i];

      // Construct edges: connect original start, internal vertices, and original end
      E_start.push_back(E(i, 0)); // Original start vertex
      E_end.push_back(V_int[0]);  // Connect to first internal vertex

      for (int j = 0; j < n_e[i] - 1; ++j) {
        E_start.push_back(V_int[j]);
        E_end.push_back(V_int[j + 1]);
      }

      E_start.push_back(V_int[n_e[i] - 1]);
      E_end.push_back(E(i, 1)); // Connect last internal vertex to original end
    } else {
      // No internal vertices: keep the original edge
      E_start.push_back(E(i, 0));
      E_end.push_back(E(i, 1));
      h_e.push_back(edge_lengths[i]);
    }
  }

  // Return the results as a List
  return List::create(
    Named("PtE_edge") = wrap(PtE_edge),
    Named("PtE_pos") = wrap(PtE_pos),
    Named("h_e") = wrap(h_e),
    Named("E_start") = wrap(E_start),
    Named("E_end") = wrap(E_end)
  );
}