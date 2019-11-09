#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <fstream>

#include "Grid.h"

template<class T, int dim>
class ParticleSimulationSystem{
public:
    using TV = Eigen::Matrix<T,dim,1>;

    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> v;

    CartesianGrid<3> grid;

    
    ParticleSimulationSystem() {}

    void evaluateParticlePositions() {
        grid.clear();

    }
};