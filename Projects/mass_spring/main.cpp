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
    using Mat = Eigen::Matrix<T, dim, dim>;

    // Set up particle system----------------------------------------------
    SimulationDriver<T,dim> driver;

    // Set up particles----------------------------------------------
    int N = 50;
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
    //Forces
    std::vector<Mat> Fp;
    std::vector<Mat> Fe;
    std::vector<Mat> F;

    std::vector<T> Vp0;

    // Set up particle attributes
    T E = (T)1.4e5; // initial young modulus
    T nu = (T) 0.2; // Poisson ratio
    T mu = E/((T)2 * ((T)1 + nu));
    T lambda = E * nu / (((T)1 + nu)*((T)1 - (T)2 * nu));
    T rho = (T)4 * pow(10, 2); // initial density

    pcg32 rng;
    T vp0 = grid.cellWidth*grid.cellWidth*grid.cellWidth/(T)8;
    T uniform_mass = vp0*rho;
    driver.ms.lambda = lambda;
    driver.ms.mu = mu;

    // Snow particle attributes
    T theta_c = (T) 1.9 * pow(10, -2); // critical compression
    T theta_t = (T) 7.5 * pow(10, -3); // critical stretch
    T epsilon = (T)10; // hardening coeff
    
    driver.ms.theta_c = theta_c;
    driver.ms.theta_t = theta_t;
    driver.ms.epsilon = epsilon;

    // Sampling particle positions before limiting to a boundary
    xp_og.resize(Np);
    for(int x = 0; x < N; x++) {
        for(int y = 0; y < N; y++) {
            for(int z = 0; z < N; z++) {
                int i = x;
                int j = y * N;
                int k = z * N * N;
                int index = i + j + k;
                xp_og[index](0) = (x + rng.nextFloat())*dx;
                xp_og[index](1) = (y + rng.nextFloat())*dx;
                xp_og[index](2) = (z + rng.nextFloat())*dx;
            }
        }
    }

    // Resampling to fit inside a box
    std::vector<TV> xp_new;
    TV minCorner(TV::Ones()*0.3);
    TV maxCorner(TV::Ones()*0.7);
    Sampling<T, dim>::selectInBox(xp_og, xp_new, minCorner, maxCorner);

    Fp.resize(xp_new.size());
    Fe.resize(xp_new.size());
    F.resize(xp_new.size());

    mp.resize(xp_new.size());
    vp.resize(xp_new.size());
    // Material space volume of each particle inside a 3d grid
    Vp0.resize(xp_new.size(), vp0);

    for(int i = 0; i < xp_new.size(); i++) {
        mp[i] = uniform_mass;
        vp[i] = TV::Zero();
        Fp[i] = Mat::Identity();
        Fe[i] = Mat::Identity();
        F[i] = Mat::Identity();
    }

    driver.grid = grid;
    driver.ms.m = mp;
    driver.ms.x = xp_new;
    driver.ms.v = vp;
    driver.ms.Fp = Fp;
    driver.ms.Fe = Fe;
    driver.ms.F = F;
    driver.ms.Vp0 = Vp0;

    driver.run(40);

    return 0;
    
}
