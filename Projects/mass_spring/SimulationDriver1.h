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

    SimulationDriver()
    : dt((T)1e-3) // 0.0015 for implicit
    {
        gravity.setZero();
        gravity(1) = -9.8;
        collision_stiffness = 5e3;

        //sphere = Sphere(collision_stiffness, TV::Ones()*0.4, 0.25);
        //ground = 0.1;
    }

    TV computeParticleMomentum() {
        TV result = TV::Zero();
        for(int pIdx = 0; pIdx < ms.m.size(); pIdx++) {
            for(int d = 0; d < dim; d++) {
                result(d) += ms.m[pIdx] * ms.v[pIdx](d);
            }
        }
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

            // if(X_index_space(0) > grid.res(0) || X_index_space(1) > grid.res(1) || X_index_space(2) > grid.res(2)) {
            //     std::cout << X(0) << ", " << X(1) << ", " << X(2) << "\n";
            //     std::cout << X_index_space(0) << ", " << X_index_space(1) << ", " << X_index_space(2) << "\n";
            //     std::cout << grid.cellWidth << "\n";
            // }
            
            // X
            TV w1 = TV::Zero(); 
            T base_node1 = 0;
            Sampling<T, dim>::computeWeights1D(X_index_space(0), base_node1, w1);
            // Y
            TV w2 = TV::Zero();
            T base_node2 = 0;
            Sampling<T, dim>::computeWeights1D(X_index_space(1), base_node2, w2);
            // Z
            TV w3 = TV::Zero();
            T base_node3 = 0;
            Sampling<T, dim>::computeWeights1D(X_index_space(2), base_node3, w3);

            for(int i1=0; i1 < dim; i1++) {
                float w_i1 = w1(i1);
                int node_i1 = base_node1 + i1 - 1;

                for(int i2=0; i2 < dim; i2++) {
                    T w_i1i2 = w_i1 * w2(i2);
                    int node_i2 = base_node2 + i2 - 1;

                    for(int i3=0; i3 < dim; i3++) {
                        T w_i1i2i3 = w_i1i2 * w3(i3);
                        int node_i3 = base_node3 + i3 - 1;
                        
                        int g_idx = node_i1 + grid.res(0) * node_i2 + grid.res(1) * grid.res(2) * node_i3;

                        // if (g_idx > 125) {
                        //     std::cout << node_i1 << ", " << node_i2 << ", " << node_i3 << "\n";
                        // }

                        //splat mass                            
                        grid.mg[g_idx] += (ms.m[p] * w_i1i2i3);

                        // // Splat momentum
                        for(int d = 0; d < dim; d++) {
                            grid.vgn[g_idx](d) += (w_i1i2i3 * ms.m[p] * ms.v[p](d));
                        }
                    }
                }
            }
        }

        for(int i = 0; i < grid.mg.size(); i++) {
            if(grid.mg[i] != 0) {
                grid.active_nodes.push_back(i);
                grid.vgn[i](0) /= grid.mg[i];
                grid.vgn[i](1) /= grid.mg[i];
                grid.vgn[i](2) /= grid.mg[i];
            } else {
                grid.vgn[i] = TV::Zero();
            }
        }
    }

    TV computeGridMomentum(int useVg) {
        TV result = TV::Zero();
        if(useVg == 1) {
            for(int i = 0; i < grid.mg.size(); i++) {
                for(int d = 0; d < dim; d++) {
                    result(d) += grid.mg[i] * grid.vg[i](d);
                }
            }
        }
        else 
        {
            for(int i = 0; i < grid.mg.size(); i++) {
                for(int d = 0; d < dim; d++) {
                    result(d) += grid.mg[i] * grid.vgn[i](d);
                }
            }
        }
        return result;
    }

    void addGravity() {
        for(int i = 0; i < grid.active_nodes.size(); i++) {
            int idx = grid.active_nodes[i];
            for(int d = 0; d < dim; d++) {
                grid.force[idx](d) += (grid.mg[idx] * gravity(d));
            }
        }
    }

    void addElasticity() {
        // *** TODO ***
    }

    void updateGridVelocity() {
        for(int i = 0; i < grid.active_nodes.size(); i++) {
            int idx = grid.active_nodes[i];
            for(int d=0; d < dim; d++) {
                grid.vg[idx](d) = grid.vgn[idx](d) + dt*grid.force[idx](d)/grid.mg[idx];
            }
        }
    }

    /**
     * Set domain boundary velocities
     */ 
    void setBoundaryVelocities(const int thickness) {
        int N = grid.res(0); // would need to change if decide not to do box
        int N_reverse = N - thickness;

        // X Direction
        for(int x = 0; x < thickness; x++) {
            for(int y = 0; y < N; y++) {
                for(int z = 0; z < N; z++) {
                    int idx = x + y * N + z * N * N;
                    grid.vg[idx] = TV::Zero();
                }
            }
        }

        for(int x = N_reverse; x < N; x++) {
            for(int y = 0; y < N; y++) {
                for(int z = 0; z < N; z++) {
                    int idx = x + y * N + z * N * N;
                    grid.vg[idx] = TV::Zero();
                }
            }
        }

        // Y Direction
        for(int x = 0; x < N; x++) {
            for(int y = 0; y < thickness; y++) {
                for(int z = 0; z < N; z++) {
                    int idx = x + y * N + z * N * N;
                    grid.vg[idx] = TV::Zero();
                }
            }
        }

        for(int x = 0; x < N; x++) {
            for(int y = N_reverse; y < N; y++) {
                for(int z = 0; z < N; z++) {
                    int idx = x + y * N + z * N * N;
                    grid.vg[idx] = TV::Zero();
                }
            }
        }


        // Z Direction
        for(int x = 0; x < N; x++) {
            for(int y = 0; y < N; y++) {
                for(int z = 0; z < thickness; z++) {
                    int idx = x + y * N + z * N * N;
                    grid.vg[idx] = TV::Zero();
                }
            }
        }

        for(int x = 0; x < N; x++) {
            for(int y = 0; y < N; y++) {
                for(int z = N_reverse; z < N; z++) {
                    int idx = x + y * N + z * N * N;
                    grid.vg[idx] = TV::Zero();
                }
            }
        }

    }

    /**
     * Evolve deformation gradient
     */ 
    void evolveF() {
        // *** TODO ***
    }

    void transferG2P(const T flip) {
        int Np = ms.x.size();

        for(int p = 0; p < Np; p++) {
            TV X = ms.x[p];
            TV X_index_space;
            X_index_space(0) = X(0)/grid.cellWidth;
            X_index_space(1) = X(1)/grid.cellWidth;
            X_index_space(2) = X(2)/grid.cellWidth;

            // X
            TV w1 = TV::Zero(); 
            T base_node1 = 0;
            Sampling<T, dim>::computeWeights1D(X_index_space(0), base_node1, w1);
            // Y
            TV w2 = TV::Zero();
            T base_node2 = 0;
            Sampling<T, dim>::computeWeights1D(X_index_space(1), base_node2, w2);
            // Z
            TV w3 = TV::Zero();
            T base_node3 = 0;
            Sampling<T, dim>::computeWeights1D(X_index_space(2), base_node3, w3);

            TV v_pic = TV::Zero();
            TV v_flip = ms.v[p];

            for(int x = 0; x < dim; x++) {
                T wx = w1(x);
                int node_x = base_node1 + x - 1;

                for (int y = 0; y < dim; y++) {
                    T wy = wx * w2(y);
                    int node_y = base_node2 + y - 1;

                    for (int z = 0; z < dim; z++) {
                        T wz = wy * w3(z);
                        int node_z = base_node3 + z - 1;
                        
                        int g_idx = node_x + grid.res(0) * node_y + grid.res(1) * grid.res(2) * node_z;
                        


                        for (int d = 0; d < dim; d++) {
                            v_pic(d) = v_pic(d) + wz * grid.vg[g_idx](d);
                            v_flip(d) = v_flip(d) + (wz * (grid.vg[g_idx](d) - grid.vgn[g_idx](d)));
                        }
                    }
                }
            }

            for(int d = 0; d < dim; d++) {
                ms.v[p](d) = ((T)1 - flip) * v_pic(d) + flip * v_flip(d);
                ms.x[p](d) += (dt * v_pic(d));
            }
        }
    }


    void transferParticleToGrid()
    {
        // Initialize grid data to zero
        grid.clear();

        // Particle 2 Grid
        //TV Lp = computeParticleMomentum();
        // std::cout << Lp(0) << ", " << Lp(1) << ", " << Lp(2) << "\n";
        transferP2G();

        //TV Lg0 = computeGridMomentum(0);
        // std::cout << "Grid momentum after p2g: " << Lg0(0) << ", " << Lg0(1) << ", " << Lg0(2) << "\n";

        // Compute force
        addGravity();

        // *** UNCOMMENT WHEN DONE ****
        // addElasticity();

        // Update grid velocity
        updateGridVelocity();

        // Boundary conditions
        // setBoundaryVelocities(1);

        // // Transfer Grid to Particle (including particle)
        // //TV Lg1 = computeGridMomentum(1);
        // // std::cout << "Grid momentum before g2p: " << Lg1(0) << ", " << Lg1(1) << ", " << Lg1(2) << "\n";
        // // *** UNCOMMENT WHEN DONE ****
        // // evolveF();

        transferG2P((T)0.99); // bigger faster
    }

    void run(const int max_frame)
    {
        for(int frame=1; frame<max_frame; frame++) {
            std::cout << "Frame " << frame << std::endl;

            int N_substeps = (int)(((T)1/24)/dt);
            for (int step = 1; step <= N_substeps; step++) {
                // std::cout << "Step " << step << std::endl;
                transferParticleToGrid();
            }
            mkdir("output/", 0777);
            std::string filename = "output/" + std::to_string(frame) + ".poly";
            ms.dumpPoly(filename);
            std::cout << std::endl;
        }
    }
};
