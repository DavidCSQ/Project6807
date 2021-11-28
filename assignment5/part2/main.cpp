// include FEM
#include "poly_mesh.hpp"
#include "linear_material.hpp"
#include "deformable_body.hpp"
#include "tetrahedral_mesh.hpp"
#include "tet_deformable_body.hpp"
#include "typedefs.hpp"
// include Mesh
#include "voxelizer.hpp"
#include "marching_cube.hpp"
// include Geometry
#include "GeometryExploration.hpp"
// include Eigen
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
// include std
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <set>

// compute the performance metrics for a given model stored in stl_name
// input: 
//  stl_name: the filename for the model
// output:
//  compliance: a scalar value for the compliance (L = F'U = U'*K*U)
//  num_voxels: a scalar value for the number of solid voxels
void get_K_and_fe(const materials::TetrahedralMesh<double> tet_mesh, // the tet mesh with rest configuration
    const materials::TetDeformableBody<double> tet_def_body, 
    const Eigen::MatrixXd vertices,                                  // the vertices coordinates in current configuration (size 3 * N)
    const std::vector<bool>& index_mask,                             // the index mask, which obtained from set_boundary_conditions
    Eigen::SparseMatrix<double>& K,                                  // output: the stiffness matrix
    Eigen::VectorXd& force_elastic                                   // output: the elastic force in current configuration
    ) {
    
    const size_t rows = tet_mesh.vertex().rows();
    const size_t cols = tet_mesh.vertex().cols();
    const size_t total_dim = rows * cols;

    Eigen::SparseMatrix<double> K_full = tet_def_body.ComputeStiffnessMatrix(vertices);

    Eigen::SparseMatrix<double> regularizer(K_full.rows(), K_full.cols());
    regularizer.setIdentity();

    // add a regularization term on the diagonal
    K_full += regularizer * 1e-4;

    Eigen::Matrix3Xd force_elastic_full = tet_def_body.ComputeElasticForce(vertices);

    // trim K_full and force_elastic_full with vertex_mask
    std::vector<int> index_mapping(total_dim);
    int dofs = 0;
    for (size_t i = 0;i < total_dim;i++)
        if (index_mask[i])
            index_mapping[i] = dofs ++;

    force_elastic = Eigen::VectorXd::Zero(dofs);
    std::vector<Eigen::Triplet<double> > triplet_list;
    triplet_list.clear();
    for (int k = 0;k < K_full.outerSize();k++)
        for (Eigen::SparseMatrix<double>::InnerIterator it(K_full, k); it; ++it) {
            if (index_mask[it.row()] && index_mask[it.col()]) {
                triplet_list.push_back(Eigen::Triplet<double>(index_mapping[it.row()], index_mapping[it.col()], it.value()));
            }
        }
    K = Eigen::SparseMatrix<double>(dofs, dofs);
    K.setFromTriplets(triplet_list.begin(), triplet_list.end());

    for (size_t i = 0;i < total_dim;i++)
        if (index_mask[i])
            force_elastic[index_mapping[i]] = force_elastic_full(i % 3, i / 3);
}
void write_facet(std::ofstream& file, Eigen::Vector3d p0, Eigen::Vector3d p1, Eigen::Vector3d p2) {
    file << "facet normal 0 0 0\n";
    file << "outer loop\n";
    file << "vertex " + std::to_string(p0.x()) + " " + std::to_string(p0.y()) + " " + std::to_string(p0.z()) + "\n";
    file << "vertex " + std::to_string(p1.x()) + " " + std::to_string(p1.y()) + " " + std::to_string(p1.z()) + "\n";
    file << "vertex " + std::to_string(p2.x()) + " " + std::to_string(p2.y()) + " " + std::to_string(p2.z()) + "\n";
    file << "endloop\n";
    file << "endfacet\n";
}
void write_tet_mesh(const std::string& filename, const Eigen::Matrix<double, 3, Eigen::Dynamic>& vertices, const Eigen::Matrix<int, 4, Eigen::Dynamic>& elements) {

    std::ofstream file;
    file.open(filename, std::ios::out);

    file << "solid tet_mesh \n";

    for (size_t i = 0; i < elements.cols(); ++i) {
        write_facet(file, vertices.col(elements(0, i)), vertices.col(elements(1, i)), vertices.col(elements(2, i)));
        write_facet(file, vertices.col(elements(1, i)), vertices.col(elements(3, i)), vertices.col(elements(2, i)));
        write_facet(file, vertices.col(elements(1, i)), vertices.col(elements(0, i)), vertices.col(elements(3, i)));
        write_facet(file, vertices.col(elements(0, i)), vertices.col(elements(2, i)), vertices.col(elements(3, i)));
    }

    file << "endsolid";

    file.close();
}
void SolvePerformance(std::string stl_name, double& compliance, double& num_voxels) {
    // TODO: HW5
    // part 2.2 implement the pipeline from design space to the performance space
    // compute the compliance and the num_voxels for the model
    const int dim = 3;

    // using this linear material
    materials::LinearElasticityMaterial<dim, double> linear_elasticity_material(10000000, 0.45);
    double dx = 0.25; // spacing parameter(need to be fixed)
    double density = 1.0;

    // step 1: voxelize the mesh
    mesh::Voxelizer<double> voxelizer(stl_name, dx);
    voxelizer.AdvancedVoxelization();

    // step 2: convert the voxelization results to tetrahedral mesh
    materials::Matrix3X<double> V; // Vertex Positions
	materials::Matrix4Xi<double> T; // Tetrahedral Elements
    std::vector<int> force;
    force.clear();
    std::vector<int> fixed;
    fixed.clear();
    
    //std::cout<<hexm.vertex().size();
    materials::Matrix3X<double> vertex;
    materials::Matrix8Xi<int> element;
    bool a = voxelizer.ConvertToHexMeshm(vertex,element,num_voxels);
    //std::cout<<vertex.size();
    V.resize(3,vertex.cols());
    T.resize(4,element.cols()*5);
    for(int i = 0;i<vertex.cols();i++){
        V(0,i)=vertex(0,i);
        V(1,i)=vertex(1,i);
        V(2,i)=vertex(2,i);
    }
    for(int i = 0;i<element.cols();i++){
        T(0,i*5)=element(0,i);
        T(1,i*5)=element(1,i);
        T(2,i*5)=element(2,i);
        T(3,i*5)=element(4,i);
        T(0,i*5+1)=element(1,i);
        T(1,i*5+1)=element(2,i);
        T(2,i*5+1)=element(3,i);
        T(3,i*5+1)=element(7,i);
        T(0,i*5+2)=element(2,i);
        T(1,i*5+2)=element(4,i);
        T(2,i*5+2)=element(6,i);
        T(3,i*5+2)=element(7,i);
        T(0,i*5+3)=element(1,i);
        T(1,i*5+3)=element(4,i);
        T(2,i*5+3)=element(5,i);
        T(3,i*5+3)=element(7,i);
        T(0,i*5+4)=element(1,i);
        T(1,i*5+4)=element(2,i);
        T(2,i*5+4)=element(4,i);
        T(3,i*5+4)=element(7,i);
    }
    //std::cout<<vertex.size();
    materials::TetrahedralMesh<double> tet_mesh = materials::TetrahedralMesh<double>(V,T);
    // step 3: construct boundary conditions for fem
    const size_t rows = tet_mesh.vertex().rows();
    const size_t cols = tet_mesh.vertex().cols();
    //std::cout<<rows<<" "<<cols;
    const size_t total_dim = rows * cols;
    std::vector<bool> index_mask;
    Eigen::VectorXd f_ext;
    index_mask.resize(total_dim);
    for (size_t i = 0;i < total_dim;i++)
        index_mask[i] = true;

    // fix the vertices on the left face and remove them from variable list
    //std::cout<<"test"<<voxelizer.pmax()(2)<<dx<<'\n'<<tet_mesh.vertex().row(2);
    for (size_t i = 0; i < tet_mesh.vertex().cols(); ++i) {
        if ((tet_mesh.vertex()(0, i) <= voxelizer.pmin()(0)+dx)||(tet_mesh.vertex()(0, i) >= voxelizer.pmax()(0)-2*dx)) {
            for (size_t j = 0;j < 3;j++) {
                index_mask[i * 3 + j] = false;//std::cout<<"true";
            }
        }
    }

    // apply external forces on each nodes of bottom right region
    Eigen::VectorXd f_ext_full = Eigen::VectorXd::Zero(total_dim);
    for (size_t i = 0; i < tet_mesh.vertex().cols(); ++i) {
        if (tet_mesh.vertex()(2, i) >= voxelizer.pmax()(2)-2*dx) {
            f_ext_full(3 * i + 2) = -5000.0;//std::cout<<"true";
        }
    }
    //std::cout<<vertex.size();
    // trim f_ext_full
    std::vector<int> index_mapping(total_dim);
    int dofs = 0;
    for (size_t i = 0;i < total_dim;i++)
        if (index_mask[i])
            index_mapping[i] = dofs ++;
    f_ext = Eigen::VectorXd(dofs); // trimed external forces
    for (size_t i = 0;i < total_dim; i++)
        if (index_mask[i])
            f_ext[index_mapping[i]] = f_ext_full[i];
    // step 4: compute compliance by the fem implemented in assignment 4, please use the linear model we pre-defined
    materials::TetDeformableBody<double> tet_def_body(linear_elasticity_material, tet_mesh.vertex(), density, tet_mesh);

    Eigen::SparseMatrix<double> K;
    Eigen::VectorXd fe;
    //std::cout<<vertex.size();
    get_K_and_fe(tet_mesh, tet_def_body, tet_mesh.vertex(), index_mask, K, fe);
    //std::cout<<vertex.size();
    // solve K*u = f_ext by eigen cg solver
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper> solver;
    solver.setMaxIterations(1000000); //just a large number

    std::cout << "Time to solve" << std::endl;
    Eigen::VectorXd u = solver.compute(K).solve(f_ext);
    std::cout<<u.size()<<f_ext.size();
    int cnt = 0;
    compliance = 0;
    size_t num_vertices;
    num_vertices = tet_mesh.vertex().cols();
    for (size_t i = 0;i < num_vertices; i++) 
        for (size_t j = 0;j < 3;j++) {
            if (index_mask[i * 3 + j]) {
                compliance += u(cnt)*f_ext(cnt);
                cnt+=1;
            }
    }
    Eigen::Matrix<double, 3, Eigen::Dynamic> x = tet_mesh.vertex();
    cnt = 0;
    for (size_t i = 0;i < num_vertices; i++) 
        for (size_t j = 0;j < 3;j++) {
            if (index_mask[i * 3 + j]) {
                x(j, i) += u(cnt++);
            }
    }
    write_tet_mesh(std::string(PROJECT_SOURCE_DIR) + "/data/assignment5/mydeformed_linear.stl", x, tet_mesh.element());

}
    
    // step 5: compute num_voxels

int main(int argc, char *argv[])
{ 
    std::cerr << "Welcome to Assignment 5" << std::endl;
    
    int N = 200000;
    std::vector<Eigen::Vector2d> p1_input; 

    p1_input.clear();

    // fix random seed to get repeatable results
    srand (1);

    for (int i = 0; i < N; i++) {
        p1_input.push_back(Eigen::Vector2d(i + 1, N - i) + Eigen::Vector2d::Random() * 2.0);
    }

    // test 2.1: fast 2d Pareto front
    std::vector<Eigen::Vector2d> p1_result = geometry::ParetoFront2D(p1_input);
    // print 2.1 results
    std::ofstream file1;
    file1.open(PROJECT_SOURCE_DIR"/data/assignment5/q1_result.txt");
    file1 << "print 2.1 test result" << std::endl;
    file1 << "Totol number of points: " << p1_result.size() << std::endl;
    for (int i = 0; i < p1_result.size(); i++) {
        file1 << "P" << i << ": " << std::endl;
        file1 << p1_result[i] << std::endl;
    }
    file1.close();

    // 2.2: implement the pipeline from design space to the performance space
    std::vector<Eigen::Vector2d> sample_performance;
    
    double compliance, num_voxels;
    
    // TODO: HW5
    // part 2.2 mapping from design space to performance space
    // Once you debug the bridge example correct, comment the code before dash line
    // and uncomment the code after dash line to run the test on 121 bridges
    //SolvePerformance(PROJECT_SOURCE_DIR"/data/assignment5/bridge.stl", compliance, num_voxels);
    //SolvePerformance(PROJECT_SOURCE_DIR"/data/assignment5/CSG/assn5_meshes/bridge_r_40_o_-25.stl", compliance, num_voxels);
    //std::cout << "result " << compliance << " " << num_voxels << std::endl;
    // -----------------------------------------------------------------------
     std::string base(PROJECT_SOURCE_DIR"/data/assignment5/CSG/assn5_meshes/bridge");
     int radius_start = 30;
     int radius_end   = 40;
     int offset_start = -30;
     int offset_end   = -20;

     std::ofstream file2;
     file2.open("q2_result.txt");
     file2 << "print 2.2 test result" << std::endl << std::endl;

     int count = 0;
     for (int r = radius_start; r <= radius_end; r++) {
         for (int o = offset_start; o <= offset_end; o++) {
             std::string mesh_name = base + "_r_" + std::to_string(r) + "_o_" + std::to_string(o) + ".stl";
             std::cout << "bridge_r_" + std::to_string(r) + "_o_" + std::to_string(o) << std::endl;
             SolvePerformance(mesh_name, compliance, num_voxels);
             std::cerr << "compliance = " << compliance << ", num_voxels = " << num_voxels << std::endl;
             sample_performance.push_back(Eigen::Vector2d(compliance, num_voxels));
             file2 << "bridge_r_" + std::to_string(r) + "_o_" + std::to_string(o) + ".stl" << std::endl;
             file2 << "Compliance: " +  std::to_string(compliance) << std::endl;
             file2 << "Total mass: " + std::to_string(num_voxels) << std::endl << std::endl;
             count++;
         }
     }
     file2.close();

     std::vector<Eigen::Vector2d> p2_result = geometry::ParetoFront2D(sample_performance);
     // print Q3 results
     file2.open("q2_pareto_front_result.txt");
     file2 << "print 2.2 pareto front result" << std::endl;
     file2 << "Totol number of points: " << p2_result.size() << std::endl;
     for (int i = 0; i < p2_result.size(); i++) {
         file2 << "P" << i << ": " << std::endl;
         file2 << p2_result[i] << std::endl;
     }
     file2.close();

    return 0;
}
