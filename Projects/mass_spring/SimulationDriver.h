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
    : dt((T)0.00001) // 0.0015 for implicit
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
            
            // X
            TV w1 = TV::Zero(); 
            T base_node1 = (T)0;;
            Sampling<T, dim>::computeWeights1D(X_index_space(0), base_node1, w1);
            // Y
            TV w2 = TV::Zero();
            T base_node2 = (T)0;
            Sampling<T, dim>::computeWeights1D(X_index_space(1), base_node2, w2);
            // Z
            TV w3 = TV::Zero();
            T base_node3 = (T)0;
            Sampling<T, dim>::computeWeights1D(X_index_space(2), base_node3, w3);

            for(int i1=0; i1 <= dim; i1++) {
                float w_i1 = w1(i1);
                int node_i1 = base_node1 + i1 - 1;

                for(int i2=0; i2 < dim; i2++) {
                    T w_i1i2 = w_i1 * w2(i2);
                    int node_i2 = base_node2 + i2 - 1;

                    for(int i3=0; i3 < dim; i3++) {
                        T w_i1i2i3 = w_i1i2 * w3(i3);
                        int node_i3 = base_node3 + i3 - 1;
                        
                        int g_idx = node_i1 + (grid.res(0)-1) * node_i2 + (grid.res(1)-1) * (grid.res(2)-1) * node_i3;

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
            if(grid.mg[i] > (T)0) {
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
            for(int d = 0; d < dim; d++) {
                grid.force[i](d) += (grid.mg[i] * gravity(d));
            }
        }
    }

    void addElasticity() {
        // *** TODO ***
    }

    void updateGridVelocity() {
        for(int i = 0; i < grid.active_nodes.size(); i++) {
            for(int d=0; d < dim; d++) {
                grid.vg[i](d) = grid.vgn[i](d) + dt*grid.force[i](d)/grid.mg[i];
            }
        }
    }

    /**
     * Set domain boundary velocities
     */ 
    void setBoundaryVelocities(const int thickness) {
        int N = grid.vg.size();
        int N_reverse = (N-1) - thickness;

        // X Direction
        for(int i = 0; i < thickness; i++) {
            for(int y = 0; y < N; y++) {
                for(int z = 0; z < N; z++) {
                    int idx = i + y * (N-1) + z * (N-1) * (N-1);
                    grid.vg[idx] = (T)0;
                }
            }
        }

        for(int i = N_reverse; i < N; i--) {
            for(int y = 0; y < N; y++) {
                for(int z = 0; z < N; z++) {
                    int idx = i + y * (N-1) + z * (N-1) * (N-1);
                    grid.vg[idx] = (T)0;
                }
            }
        }

        // *** TO DO ***
        // Finish the rest of the directions

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
            TV X_index_space = X;
            X_index_space(0) /= grid.cellWidth;
            X_index_space(1) /= grid.cellWidth;
            X_index_space(2) /= grid.cellWidth;

            // X
            TV w1 = TV::Zero(); 
            T base_node1 = (T)0;;
            Sampling<T, dim>::computeWeights1D(X_index_space(0), base_node1, w1);
            // Y
            TV w2 = TV::Zero();
            T base_node2 = (T)0;
            Sampling<T, dim>::computeWeights1D(X_index_space(1), base_node2, w2);
            // Z
            TV w3 = TV::Zero();
            T base_node3 = (T)0;
            Sampling<T, dim>::computeWeights1D(X_index_space(2), base_node3, w3);

            TV v_pic = TV::Zero();
            TV v_flip = ms.v[p];

            for(int x = 0; x < dim; x++) {
                T wx = w1(x);
                T node_x = base_node1 + (x - 1);

                for (int y = 0; y < dim; y++) {
                    T wy = wx * w2(y);
                    T node_y = base_node2 + (y - 1);

                    for (int z = 0; z < dim; z++) {
                        T wz = wy * w3(z);
                        T node_z = base_node3 + (z - 1);
                        
                        int g_idx = node_x + (grid.res(0)-1) * node_y + (grid.res(1)-1) * (grid.res(2)-1) * node_z;

                        for (int d = 0; d < dim; d++) {
                            v_pic(d) += (wz * grid.vg[g_idx](d));
                            v_flip(d) += (wz * (grid.vg[g_idx](d) - grid.vgn[g_idx](d)));
                        }
                    }
                }
            }

            for(int d = 0; d < dim; d++) {
                ms.v[p](d) = ((T)1 - flip) * v_pic(d) + flip * v_flip(d);
                ms.x[p](d) += (dt*v_pic(d));
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
        // *** UNCOMMENT WHEN DONE ****
        // setBoundaryVelocities(3);

        // Transfer Grid to Particle (including particle)
        //TV Lg1 = computeGridMomentum(1);
        // std::cout << "Grid momentum before g2p: " << Lg1(0) << ", " << Lg1(1) << ", " << Lg1(2) << "\n";
        // *** UNCOMMENT WHEN DONE ****
        // evolveF();

        transferG2P((T)0.95);
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
