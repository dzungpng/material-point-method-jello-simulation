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

    // T mg[(int)grid.res(0)][(int)grid.res(1)][(int)grid.res(2)];
    // TV vg[(int)grid.res(0)][(int)grid.res(1)][(int)grid.res(2)];
    // TV force[(int)grid.res(0)][(int)grid.res(1)][(int)grid.res(2)];
 
    // int numGridCells = grid.res(0)*grid.res(1)*grid.res(2);
    // std::cout << "Number of grid cells: " << numGridCells << "\n";

    // Sample particles----------------------------------------------
    int N = 28;
    int Np = N*N*N;
    // Distance between per particle
    T dx = (T)1/(N-1);

    // Grid parameters----------------------------------------------
    CartesianGrid<T, dim> grid;
    TV res;
    res(0) = N/2;
    res(1) = N/2;
    res(2) = N/2;
    grid.res = res;
    grid.nCells = res(0)*res(1)*res(2);
    grid.cellWidth = (T)1/N;

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

    // for(int frame=1; frame<max_frame; frame++) {

    //     int N_substeps = (int)(((T)1/24)/dt);
    //     for (int step = 1; step <= N_substeps; step++) {
    //         //std::cout << "Step " << step << std::endl;
            
    //         // Clear and reset grid
    //         memset(mg, 0, (sizeof(mg)));
    //         memset(vg, 0, (sizeof(vg)));
    //         memset(force, 0, (sizeof(force)));
            
    //         // Particle 2 Grid
    //         // computeParticleMomentum
    //         TV result;
    //         for(int pIdx = 0; pIdx < Np; pIdx++) {
    //             for(int d = 0; d < dim; d++) {
    //                 result(d) = result(d) + mp[pIdx] * vp[pIdx](d);
    //             }
    //         }
    //         // std::cout << "Part momentum before p2g: " << result(0) << ", " << result(1) << ", " << result(2) << "\n";

    //         // transferP2G
    //         // At this point mg and cg should be resized and initialized
    //         for(int pIdx = 0; pIdx < Np; pIdx++) {
    //             TV x_p = samples[pIdx];

    //             // Converting x_p to index space (aka dx between each grid cell = 1)
    //             // ------------------------------------------------
    //             TV x_p_index_space;
    //             x_p_index_space(0) = x_p(0)/grid.dx;
    //             x_p_index_space(1) = x_p(1)/grid.dx;
    //             x_p_index_space(2) = x_p(2)/grid.dx;

    //             // ComputeWeights1D for each x, y, and z
    //             // -------------------------------------------------
    //             // X
    //             TV w1; 
    //             T base_node1;
    //             Sampling<T, dim>::computeWeights1D(x_p_index_space(0), base_node1, w1);
    //             // Y
    //             TV w2;
    //             T base_node2;
    //             Sampling<T, dim>::computeWeights1D(x_p_index_space(1), base_node2, w2);
    //             // Z
    //             TV w3;
    //             T base_node3;
    //             Sampling<T, dim>::computeWeights1D(x_p_index_space(0), base_node3, w3);

    //             // std::cout << "w1: " << w1(0) << ", " << w1(1) << ", " << w1(2) << "\n";
    //             // std::cout << "w2: " << w2(0) << ", " << w2(1) << ", " << w2(2) << "\n";
    //             // std::cout << "w3: " << w3(0) << ", " << w3(1) << ", " << w3(2) << "\n";
    //             // std::cout << "\n";
    //             int node_num = 0;
    //             for(int i1=0; i1 < dim; i1++) {
    //                 float w_i1 = w1(i1);
    //                 int node_i1 = base_node1 + (i1 - 1);

    //                 for(int i2=0; i2 < dim; i2++) {
    //                     T w_i1i2 = w_i1 * w2(i2);
    //                     int node_i2 = base_node2 + (i2 - 1);

    //                     for(int i3=0; i3 < dim; i3++) {
    //                         T w_i1i2i3 = w_i1i2 * w3(i3);
    //                         int node_i3 = base_node3 + (i3 - 1);

    //                         // splat mass
    //                         std::cout << node_i1 << ", " << node_i2 << ", " << node_i3 << "\n";
    //                         // if(node_i1 )
    //                         // mg[node_i1][node_i2][node_i3] = mg[node_i1][node_i2][node_i3] + mp[pIdx] * w_i1i2i3; 
                            
    //                         //
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    return 0;
    
}
