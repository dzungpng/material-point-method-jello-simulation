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
    int N = 16;
    int Np = N*N*N;
    // Distance between per particle
    T dx = (T)1/(T)N;

    // Grid parameters----------------------------------------------
    CartesianGrid<T, dim> grid;
    Eigen::Matrix<int, dim,1> res;

    // 5 x 5 x 5 grid 
    grid.cellWidth = /*(T)1/(T)4*/0.02;
    res(0) = (T)1/grid.cellWidth + 1;
    res(1) = (T)1/grid.cellWidth + 1;
    res(2) = (T)1/grid.cellWidth + 1;
    grid.res = res;
    grid.nCells = res(0)*res(1)*res(2);


    // Sample particles----------------------------------------------
    std::vector<T> mp;
    std::vector<TV> vp;
    std::vector<TV> xp_og;
    std::vector<TV> Fp;
    // mp.resize(Np);
    // vp.resize(Np);
    // Fp.resize(Np);
    xp_og.resize(Np);

    // Set up particle attributes
    T E = 1e4;
    T nu = .3;
    T mu = E/((T)2 * ((T)1 + nu));
    T lambda = E * nu / (((T)1 + nu)*((T)1 - (T)2 * nu));
    T rho = (T)1000;
    T Vp0 = (T)grid.cellWidth*grid.cellWidth*grid.cellWidth/8;
    pcg32 rng;
    T uniform_mass = Vp0*rho;
    // T uniform_mass = (T)1;

    for(int x = 0; x < N; x++) {
        for(int y = 0; y < N; y++) {
            for(int z = 0; z < N; z++) {
                int i = x;
                int j = y * N;
                int k = z * N * N;
                int index = i + j + k;
                // mp[index] = uniform_mass;
                xp_og[index](0) = (x + rng.nextFloat())*dx;
                xp_og[index](1) = (y + rng.nextFloat())*dx;
                xp_og[index](2) = (z + rng.nextFloat())*dx;
                // vp[index] = TV::Zero();
                // Fp[index] = TV::Zero();
            }
        }
    }

    // Resampling to fit inside a box
    std::vector<TV> xp_new;
    TV minCorner(TV::Ones()*0.2);
    TV maxCorner(TV::Ones()*0.8);
    Sampling<T, dim>::selectInBox(xp_og, xp_new, minCorner, maxCorner);
    Fp.resize(xp_new.size());
    mp.resize(xp_new.size());
    vp.resize(xp_new.size());

    for(int i = 0; i < xp_new.size(); i++) {
        mp[i] = uniform_mass;
        vp[i] = TV::Zero();
        Fp[i] = TV::Zero();
    }

    driver.grid = grid;
    driver.ms.m = mp;
    driver.ms.x = xp_new;
    driver.ms.v = vp;
    driver.ms.f = Fp;
    
    driver.run(2);

    return 0;
    
}
