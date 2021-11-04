#include <Eigen/Core>
#include <limits>
#include <math.h>

// TODO: HW3
// Assignment 3, Part 2.1. 
// Fill the values inside linear weight matrix W of shape (n, m)
// V is a vertices matrix in shape (n, 3), each row is the position of a vertex in mesh
// C is a control points matrix in shape (m, 3), each row is the position of a control point
void ComputeLinearSkinningWeights(
	const Eigen::MatrixXd& V, 
	const Eigen::MatrixXd& C, 
	Eigen::MatrixXd& W) 
{
	W.resize(V.rows(), C.rows());
	for(int i = 0;i<V.rows();i++){
	    for(int j = 0;j<C.rows();j++){
	        W(i,j) = 1/sqrt((V(i,0)-C(j,0))*(V(i,0)-C(j,0))+(V(i,1)-C(j,1))*(V(i,1)-C(j,1))+(V(i,2)-C(j,2))*(V(i,2)-C(j,2)));
        }
    }
}