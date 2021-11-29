#pragma once

#include <Eigen/Sparse>
#include <igl/edge_lengths.h>
#include <igl/face_areas.h>
#include <igl/dihedral_angles.h>
#include <igl/volume.h>
#include "dihedral_sine.hpp"
#include <math.h>


// TODO: HW3
// Assignment 3, Part 3.2.
/* Implement your code here. */
// Implement the function to compute the cotangent laplacian matrix L 
// V is the vertex matrix of shape (n, 3), each row is the position of a vertex in the mesh
// F is the element index matrix of shape (m, 4), each row is the vertex indices of a tetrahedron
// L is the output cotangent laplacian matrix of shape (n, n), and it's a sparse matrix.
// Hints:
	// 1. For each tetrahedron, loop over each of its edge,
	//    consider which part of the L matrix this edge in this tetrahedron contributes to
	// 2. compute the cos and sin of the dihedral angle by the law of diehedral angles http://mathworld.wolfram.com/Tetrahedron.html
	//	  specifically, compute the sin and cos of dihedral angles from the edge lengths, face areas and tet volume
	// 3. build the triplets <row, col, value> in IJV
using Vector3 = Eigen::Matrix<double, 3, 1>;
//fill in IJV
void check( 
    const Eigen::MatrixXd& V, 
    int v1,
    int v2,
    double theta,
    std::vector<Eigen::Triplet<double> >& IJV,
    int v3,
    int v4)
{
    Vector3 le = Vector3(V(v3,0)-V(v4,0),V(v3,1)-V(v4,1),V(v3,2)-V(v4,2));
    double value = sqrt(le.dot(le))*cos(theta)/(6*sin(theta));
    IJV.push_back(Eigen::Triplet<double>(v1,v2,value));
    IJV.push_back(Eigen::Triplet<double>(v2,v1,value));
}

void cotangent_laplacian(
	const Eigen::MatrixXd& V, 
	const Eigen::MatrixXi& F, 
	Eigen::SparseMatrix<double>& L) 
{
	L.resize(V.rows(), V.rows());
	std::vector<Eigen::Triplet<double> > IJV;
	IJV.clear();
    Eigen::Matrix<double,Eigen::Dynamic,6> theta;
    Eigen::Matrix<double,Eigen::Dynamic,6> cos_theta;
    igl::dihedral_angles(V,F,theta,cos_theta);//compute theta
    //sstd::cout<<theta;
    for(int i = 0;i<cos_theta.rows();i++){
	    check(V,F(i,1),F(i,2),theta(i,0),IJV,F(i,0),F(i,3));
        check(V,F(i,0),F(i,2),theta(i,1),IJV,F(i,1),F(i,3));
        check(V,F(i,0),F(i,1),theta(i,2),IJV,F(i,2),F(i,3));
        check(V,F(i,0),F(i,3),theta(i,3),IJV,F(i,1),F(i,2));
        check(V,F(i,1),F(i,3),theta(i,4),IJV,F(i,0),F(i,2));
        check(V,F(i,2),F(i,3),theta(i,5),IJV,F(i,0),F(i,1));
    }
	/* Implement your code here. */
    for(int i = 0;i<V.rows();i++){
        double value = 0;
        for(int j = 0;j<IJV.size();j++){
            if (IJV[j].row()==i){
                value-=IJV[j].value();
            }
        }
        IJV.push_back(Eigen::Triplet<double>(i,i,value));
    }
    // for(int i = 0;i<IJV.size();i++){
    //     std::cout<<IJV[i].row()<<" "<<IJV[i].col()<<" "<<IJV[i].value()<<"\n";
    // }
	// Set From Triplets Sums all Triplets with the same indices
	L.setFromTriplets(IJV.begin(), IJV.end());
}