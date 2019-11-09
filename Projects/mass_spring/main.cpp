#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <random>
#include <chrono>

#include "SimulationDriver.h"
#include "Geometry/Grid.h"
#include "Sampling.h"

int main(int argc, char* argv[])
{
    using T = float;
    constexpr int dim = 3; 
    using TV = Eigen::Matrix<T,dim,1>;

    // Set up particle system
    SimulationDriver<T,dim> driver;

    // Grid parameters
    CartesianGrid<T, dim> grid;
    TV res = (grid.maxCorner-grid.minCorner);
    res(0) = res(0) / grid.dx + 1;
    res(1) = res(1) / grid.dx + 1;
    res(2) = res(2) / grid.dx + 1;

    int numSamples = 1000;
    // Sample particles
    std::vector<TV> samples = Sampling<T, dim>::stratefiedSampling(numSamples);

    // Set up particle attributes
    std::vector<T> mp;
    std::vector<TV> xp;
    std::vector<TV> vp;

    int N = 512;
    int N_points = N*N*N;
    T dx = (T)1/(N-1);
    mp.resize(N_points);
    xp.resize(N_points);
    vp.resize(N_points);
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            int id = i*N+j;
            mp[id] = (T)1/N_points;
            xp[id](0) = (i-1)*dx;
            xp[id](2) = (j-1)*dx;
            xp[id](1) = 1;
            vp[id] = TV::Zero();
        }
    }

    // simulate
    driver.run(1);


    return 0;
}
