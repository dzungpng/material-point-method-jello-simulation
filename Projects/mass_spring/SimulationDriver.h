#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>

#include <sys/stat.h>
#include <iostream>
#include "MassSpringSystem.h"
#include "ParticleSimulationSystem.h"

template<class T, int dim>
class SimulationDriver{
public:
    using TV = Eigen::Matrix<T,dim,1>;
    using SpMat = Eigen::SparseMatrix<T>;
    using Vec = Eigen::Matrix<T,Eigen::Dynamic,1>;

    // MassSpringSystem<T,dim> ms;
    ParticleSimulationSystem<T, dim> ps;

    T dt;
    TV gravity;


    SimulationDriver()
      // : dt((T)0.00001) 
      : dt((T)0.0015)  // 150 times bigger dt than explicit. 
      //We can't go arbitrarily large because we are still doing approximations 
      //to the non-linear problem using taylor expansion.
    {
        gravity.setZero();
        gravity(1) = -9.8;

    }

    void run(const int max_frame)
    {
        for(int frame=1; frame<max_frame; frame++) {
            std::cout << "Frame " << frame << std::endl;

            int N_substeps = (int)(((T)1/24)/dt);
            for (int step = 1; step <= N_substeps; step++) {
                //std::cout << "Step " << step << std::endl;
                // advanceOneStepExplicitIntegration();
		        // advanceOneStepImplicitIntegration();
                ps.computeParticleMomentum();
            }
            mkdir("output/", 0777);
            std::string filename = "output/" + std::to_string(frame) + ".poly";
            ps.dumpPoly(filename);
            std::cout << std::endl;
        }
    }
};
