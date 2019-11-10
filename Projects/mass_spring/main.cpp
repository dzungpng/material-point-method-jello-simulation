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
    grid.res = res;
 
    int numGridCells = grid.res(0)*grid.res(1)*grid.res(2);
    std::cout << "Number of grid cells: " << numGridCells << "\n";

    int numSamples = 100;
    // Sample particles
    std::vector<TV> samples = Sampling<T, dim>::stratefiedSampling(numSamples);
    std::cout << "Number of samples: " << samples.size() << "\n";

    // Set up particle attributes
    T E = 1e4;
    T nu = .3;
    T mu = E/(2 * (1 + nu));
    T lambda = E * nu / ((1 + nu)*(1-2*nu));
    T rho = 1000.f;
    int Np = samples.size();
    T Vp0 = (grid.dx*grid.dx*grid.dx)/6;
    
    std::vector<T> mp;
    std::vector<TV> xp;
    std::vector<TV> vp;
    std::vector<TV> fp;

    mp.resize(Np);
    xp.resize(Np);
    vp.resize(Np);
    fp.resize(Np);
    for(int i=0; i<Np; i++) {
        mp[i] = Vp0;
        xp[i] = TV::Zero();
        vp[i] = TV::Zero();
        fp[i] = TV::Zero();
    }

    // ParticleSimulationSystem<T, dim> ps;
    T dt = (T)0.001;
    TV gravity;
    gravity.setZero();
    gravity(1) = -9.8;

    int max_frame = 100;
    for(int frame=1; frame<max_frame; frame++) {
        std::cout << "Frame " << frame << std::endl;

        int N_substeps = (int)(((T)1/24)/dt);
        for (int step = 1; step <= N_substeps; step++) {
            //std::cout << "Step " << step << std::endl;
            
            // Clear and reset grid
            grid.clear();
            grid.mg.resize(numGridCells);
            grid.vg.resize(numGridCells);
            grid.force.resize(numGridCells);
            for(int cellIdx = 0; cellIdx < numGridCells; cellIdx++) {
                grid.mg[cellIdx] = 0.f;
                grid.vg[cellIdx] = TV::Zero();
                grid.force[cellIdx] = TV::Zero();
            }
            
            // Particle 2 Grid
            // computeParticleMomentum
            TV result;
            for(int pIdx = 0; pIdx < Np; pIdx++) {
                for(int d = 0; d < dim; d++) {
                    result(d) = result(d) + mp[pIdx] * vp[pIdx](d);
                }
            }
            // std::cout << "Part momentum before p2g: " << result(0) << ", " << result(1) << ", " << result(2) << "\n";

            // transferP2G
            // At this point mg and cg should be resized and initialized
            for(int pIdx = 0; pIdx < Np; pIdx++) {
                TV x_p = samples[pIdx];

                // Converting x_p to index space (aka dx between each grid cell = 1)
                // ------------------------------------------------
                TV x_p_index_space;
                x_p_index_space(0) = x_p(0)/grid.dx;
                x_p_index_space(1) = x_p(1)/grid.dx;
                x_p_index_space(2) = x_p(2)/grid.dx;

                // ComputeWeights1D for each x, y, and z
                // -------------------------------------------------
                // X
                TV w1; 
                T base_node1;
                Sampling<T, dim>::computeWeights1D(x_p_index_space(0), base_node1, w1);
                // Y
                TV w2;
                T base_node2;
                Sampling<T, dim>::computeWeights1D(x_p_index_space(1), base_node2, w2);
                // Z
                TV w3;
                T base_node3;
                Sampling<T, dim>::computeWeights1D(x_p_index_space(0), base_node3, w3);

                std::cout << "w1: " << w1(0) << ", " << w1(1) << ", " << w1(2) << "\n";
                std::cout << "w2: " << w2(0) << ", " << w2(1) << ", " << w2(2) << "\n";
                std::cout << "w3: " << w3(0) << ", " << w3(1) << ", " << w3(2) << "\n";
                std::cout << "\n";
            }
        }
    }

    return 0;
}
