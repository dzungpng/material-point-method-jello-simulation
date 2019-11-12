#pragma once
#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>

#include <sys/stat.h>
#include <iostream>
#include "ParticleSystem.h"
#include "Geometry/Grid.h"
#include "Sampling.h"

template<class T, int dim>
class SimulationDriver{
public:
    using TV = Eigen::Matrix<T,dim,1>;
    using SpMat = Eigen::SparseMatrix<T>;
    using Vec = Eigen::Matrix<T,Eigen::Dynamic,1>;

    MassSpringSystem<T,dim> ms;
    CartesianGrid<T, dim> grid;

    T dt;
    TV gravity;
    T ground;
    T collision_stiffness;
    T EPSILON = 1e-8;

    SimulationDriver()
    : dt((T)0.00001) // 0.0015 for implicit
    {
        gravity.setZero();
        gravity(1) = -9.8;
        collision_stiffness = 5e3;

        //sphere = Sphere(collision_stiffness, TV::Ones()*0.4, 0.25);
        ground = 0.1;
    }

    void run(const int max_frame)
    {
        for(int frame=1; frame<max_frame; frame++) {
            std::cout << "Frame " << frame << std::endl;

            int N_substeps = (int)(((T)1/24)/dt);
            for (int step = 1; step <= N_substeps; step++) {
                // std::cout << "Step " << step << std::endl;
                transferParticleToGrid();
                return;
            }
            mkdir("output/", 0777);
            std::string filename = "output/" + std::to_string(frame) + ".poly";
            ms.dumpPoly(filename);
            std::cout << std::endl;
        }
    }


    TV computeParticleMomentum(const std::vector<T> mp, const std::vector<TV> vp) {
        TV result = TV::Zero();
        for(int pIdx = 0; pIdx < mp.size(); pIdx++) {
            for(int d = 0; d < dim; d++) {
                result(d) += mp[pIdx] * vp[pIdx](d);
            }
        }
        // std::cout << result(0) << ", " << result(1) << ", " << result(2) << "\n";
        return result;
    }

    void transferP2G() 
    {
        int Np = ms.x.size();
        for(int p = 0; p < Np; p++) {
            TV X = ms.x[p];
            TV X_index_space;
            X_index_space(0) = X(0)/grid.cellWidth;
            X_index_space(1) = X(1)/grid.cellWidth;
            X_index_space(2) = X(2)/grid.cellWidth;
            
            // X
            TV w1; 
            T base_node1;
            Sampling<T, dim>::computeWeights1D(X_index_space(0), base_node1, w1);
            // Y
            TV w2;
            T base_node2;
            Sampling<T, dim>::computeWeights1D(X_index_space(1), base_node2, w2);
            // Z
            TV w3;
            T base_node3;
            Sampling<T, dim>::computeWeights1D(X_index_space(2), base_node3, w3);

            for(int i1=0; i1 <= dim; i1++) {
                float w_i1 = w1(i1);
                int node_i1 = base_node1 + i1;

                for(int i2=0; i2 < dim; i2++) {
                    T w_i1i2 = w_i1 * w2(i2);
                    int node_i2 = base_node2 + i2;

                    for(int i3=0; i3 < dim; i3++) {
                        T w_i1i2i3 = w_i1i2 * w3(i3);
                        int node_i3 = base_node3 + i3;
                        
                        int g_idx = node_i1 + (grid.res(0)-1) * node_i2 + (grid.res(1)-1) * (grid.res(2)-1) * node_i3;
                        
                        // splat mass 
                        std::cout << g_idx << "\n";
                         
                        grid.mg[g_idx] += ms.m[p] * w_i1i2i3;

                        // Splat momentum
                        for(int d = 0; d < dim; d++) {
                            grid.vg[g_idx](d) += (w_i1i2i3 * ms.m[p]) * ms.v[p](d);
                        }
                    }
                }
            }
        }

        for(int i = 0; i < grid.mg.size(); i++) {
            if(grid.mg[i] > EPSILON) {
                grid.active_nodes[i] = 1;
                grid.vg[i](0) /= grid.mg[i];
                grid.vg[i](1) /= grid.mg[i];
                grid.vg[i](2) /= grid.mg[i];
            } else {
                grid.vg[i] = TV::Zero();
            }
        }
    }

    void transferParticleToGrid()
    {
        grid.clear();
        TV Lg = computeParticleMomentum(ms.m, ms.v);
        transferP2G();
    }
};
