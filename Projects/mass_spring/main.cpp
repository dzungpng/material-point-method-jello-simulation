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

    SimulationDriver<T,dim> driver;

    // set up mass spring system
    T youngs_modulus = 2000;
    T damping_coeff = 2; // 0

    SegmentMesh<T, dim> mesh = SegmentMesh<T, dim>("../../torus.vtk");

    // simulate
    driver.ms.segments = mesh.segments;
    driver.ms.m = mesh.m;
    driver.ms.v = mesh.v;
    driver.ms.x = mesh.x;
    driver.ms.youngs_modulus = youngs_modulus;
    driver.ms.damping_coeff = damping_coeff;
    driver.ms.node_is_fixed = mesh.node_is_fixed;
    driver.ms.rest_length = mesh.rest_length;

    driver.run(120);

    return 0;
}
