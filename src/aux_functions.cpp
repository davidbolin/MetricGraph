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

Rcpp::List assemble_fem(Eigen::MatrixXd E, Eigen::VectorXd h_e, int nV){

    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> trp_C, trp_G, trp_B;
    int i;
    int v1,v2;

    Eigen::SparseMatrix<double> C(nV,nV), G(nV,nV), B(nV,nV);

    for(i=0; i<E.rows(); i++){
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
    }

    C.setFromTriplets(trp_C.begin(), trp_C.end());
    G.setFromTriplets(trp_G.begin(), trp_G.end());          
    B.setFromTriplets(trp_B.begin(), trp_B.end());      

    Rcpp::List out;
    out["C"] = C;
    out["G"] = G;
    out["B"] = B;
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

Eigen::MatrixXd interpolate2_aux(Eigen::MatrixXd lines, Eigen::VectorXd pos, int normalized = 0){
    int size_return = pos.size();
    int i,j;
    Eigen::MatrixXd out_mat(size_return,2);
    Eigen::VectorXd dist_vec(lines.rows());
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
    }

    return(out_mat);
}
