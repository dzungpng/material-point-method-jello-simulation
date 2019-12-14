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
#include "Sampling/Sampling.h"
#include "Geometry/SampledMesh.h"
#include "Geometry/PolyParser.h"

int main(int argc, char* argv[])
{
    using T = float;
    constexpr int dim = 3; 
    using TV = Eigen::Matrix<T,dim,1>;
    using Mat = Eigen::Matrix<T, dim, dim>;

    // Set up particle system----------------------------------------------
    SimulationDriver<T,dim> driver;

    // Set up particles----------------------------------------------
    int N = 100;
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
    std::vector<Mat> Fp;
    std::vector<T> Vp0;

    // Set up particle attributes
    T E = 1e4; // 1e5 good for sphere
    T nu = (T)0.3; // previously 0.3
    T mu = E/((T)2 * ((T)1 + nu));
    T lambda = E * nu / (((T)1 + nu)*((T)1 - (T)2 * nu));
    T rho = (T)1000;
    pcg32 rng;
    T vp0 = grid.cellWidth*grid.cellWidth*grid.cellWidth/(T)8;
    T uniform_mass = vp0*rho;
    driver.ms.lambda = lambda;
    driver.ms.mu = mu;

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

    // BOX
    TV minCorner(TV(0.4, 0.7, 0.4));
    TV maxCorner(TV(0.6, 1, 0.6));
    Sampling<T, dim>::selectInBox(xp_og, xp_new, minCorner, maxCorner);

    // SPHERE
    // Sampling<T, dim>::sampleSphere(xp_new, Np, 0.1, 0.3);

    // OBJ
    // PolyParser<T, dim> parser = PolyParser<T, dim>("mesh/torus-small.poly");
    // xp_new = parser.x;

    Fp.resize(xp_new.size());
    mp.resize(xp_new.size());
    vp.resize(xp_new.size());
    // Material space volume of each particle inside a 3d grid
    Vp0.resize(xp_new.size(), vp0);

    for(int i = 0; i < xp_new.size(); i++) {
        vp[i] = TV::Zero();
        mp[i] = uniform_mass;
        Fp[i] = Mat::Identity();
    }

    driver.grid = grid;
    driver.ms.m = mp;
    driver.ms.x = xp_new;
    driver.ms.v = vp;
    driver.ms.Fp = Fp;
    driver.ms.Vp0 = Vp0;

    // Adding geometry to sampled geometry so we can continually add to the scene later
    SampledMesh<T, dim> sampledMesh(mp, xp_new, vp, Fp, Vp0);
    driver.sampledMesh = &sampledMesh;

    driver.run(100);

    return 0;
    
}
