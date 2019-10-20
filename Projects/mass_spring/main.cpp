#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <random>
#include <chrono>

#include "Cloth.h"
#include "SimulationDriver.h"
#include "TetMesh.h"

int main(int argc, char* argv[])
{
    using T = float;
    constexpr int dim = 3;
    // using TV = Eigen::Matrix<T,dim,1>;

    SimulationDriver<T,dim> driver;

    // set up mass spring system
    T youngs_modulus = 2;
    T damping_coeff = 2; // 0

    // TESTING TETMESH STUFF
    // VTKParser<T, dim> parser;
    // vector<TV> test_x;
    // int test_N;
    // vector<Eigen::Matrix<int,2,1>> test_segments;
    // parser.buildSegments("../../torus.vtk", test_N, test_x, test_segments);
    
    
    // for (int i=0; i < test_x.size(); i++) {
    //     cout << "Vertex "  << i << " is:" << test_x[i](0) << ", " << test_x[i](1) << ", " << test_x[i](2) << "\n";
    // }    

    // for (int i=0; i < test_segments.size(); i++) {
    //     cout << "Segment "  << i << " is:" << test_segments[i](0) << ", " << test_segments[i](1) << "\n";
    // } 
    // END TESTING TETMESH STUFF


    SegmentMesh<T, dim> mesh = SegmentMesh<T, dim>("../../torus.vtk");

    // for (int i=0; i < mesh.m.size(); i++) {
    //     cout << "Mass "  << i << " is:" << mesh.m[i] << "\n";
    // }

    
    // for (int i=0; i < mesh.segments.size(); i++) {
    //     cout << "Segment "  << i << " is: " << mesh.segments[i](0) << ", " << mesh.segments[i](1) << "\n";
    // }  
    

    // for (int i=0; i < mesh.x.size(); i++) {
    //     cout << "Pos " << i << " is: " << mesh.x[i](0) << " " << mesh.x[i](1) << " " << mesh.x[i](2) << "\n";
    // }
    // Cloth<T, dim> cloth;

    // simulate
    driver.ms.segments = mesh.segments;
    driver.ms.m = mesh.m;
    driver.ms.v = mesh.v;
    driver.ms.x = mesh.x;
    driver.ms.youngs_modulus = youngs_modulus;
    driver.ms.damping_coeff = damping_coeff;
    driver.ms.node_is_fixed = mesh.node_is_fixed;
    driver.ms.rest_length = mesh.rest_length;

    driver.run(50);

    return 0;
}
