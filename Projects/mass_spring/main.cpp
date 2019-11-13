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
#include "Sampling.h"


int main(int argc, char* argv[])
{
    using T = float;
    constexpr int dim = 3; 
    using TV = Eigen::Matrix<T,dim,1>;

    // Set up particle system----------------------------------------------
    SimulationDriver<T,dim> driver;

    // Set up particles----------------------------------------------
    int N = 8;
    int Np = N*N*N;
    // Distance between per particle
    T dx = (T)1/N;

    // Grid parameters----------------------------------------------
    CartesianGrid<T, dim> grid;
    Eigen::Matrix<int, dim,1> res;

    // 5 x 5 x 5 grid 
    grid.cellWidth = (T)1/4;
    res(0) = 1/grid.cellWidth + 1;
    res(1) = 1/grid.cellWidth + 1;
    res(2) = 1/grid.cellWidth + 1;
    grid.res = res;
    grid.nCells = res(0)*res(1)*res(2);

    // Sample particles----------------------------------------------
    std::vector<T> mp;
    std::vector<TV> vp;
    std::vector<TV> xp;
    mp.resize(Np);
    vp.resize(Np);
    xp.resize(Np);

    pcg32 rng;
    T uniform_mass = (T)1/Np;

    for(int x = 0; x < N; x++) {
        for(int y = 0; y < N; y++) {
            for(int z = 0; z < N; z++) {
                int i = x;
                int j = y * (N-1);
                int k = z * (N-1) * (N-1);
                int index = i + j + k;
                mp[index] = uniform_mass;
                xp[index](0) = (x + rng.nextFloat())*dx;
                xp[index](1) = (y + rng.nextFloat())*dx;
                xp[index](2) = (z + rng.nextFloat())*dx;
                vp[index] = TV::Zero();
            }
        }
    }


    // Set up particle attributes
    T E = 1e4;
    T nu = .3;
    T mu = E/((T)2 * ((T)1 + nu));
    T lambda = E * nu / (((T)1 + nu)*((T)1 - (T)2 * nu));
    T rho = (T)1000;
    // T Vp0 = (dx*dx*dx)/6;

    driver.grid = grid;
    driver.ms.m = mp;
    driver.ms.x = xp;
    driver.ms.v = vp;
    driver.run(2);

    return 0;
    
}
