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

    // Sample particles----------------------------------------------
    int N = 16;
    int Np = N*N*N;
    // Distance between per particle
    T dx = 0.02;
    //T dx = 1;

    // Grid parameters----------------------------------------------
    CartesianGrid<T, dim> grid;
    TV res; // grid dimensions 
    res(0) = N/4; // X
    res(1) = N/4; // Y 
    res(2) = N/4; // Z
    grid.res = res;
    grid.cellWidth = (T)N/(N/4);
    std::cout << grid.cellWidth << "\n";
    grid.nCells = res(0)*res(1)*res(2);

    std::vector<T> mp;
    std::vector<TV> vp;
    std::vector<TV> xp;
    mp.resize(Np);
    vp.resize(Np);
    xp.resize(Np);

    pcg32 rng;
    T uniform_mass = 1/(T)Np;

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
    T mu = E/(2 * (1 + nu));
    T lambda = E * nu / ((1 + nu)*(1-2*nu));
    T rho = 1000.f;
    T Vp0 = (dx*dx*dx)/6;

    driver.grid = grid;
    driver.ms.m = mp;
    driver.ms.x = xp;
    driver.ms.v = vp;
    driver.run(2);

    return 0;
    
}
