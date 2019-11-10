#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <fstream>


template<class T, int dim>
class ParticleSimulationSystem{
public:
    using TV = Eigen::Matrix<T,dim,1>;

    std::vector<T> mp;
    std::vector<TV> xp;
    std::vector<TV> vp;
    
    ParticleSimulationSystem() {}

    void computeParticleMomentum() {
        // TODO
        return;
    }

    void dumpPoly(std::string filename)
    {
        std::ofstream fs;
        fs.open(filename);
        fs << "POINTS\n";
        int count = 0;
        for (auto X : xp) {
            fs << ++count << ":";
            for (int i = 0; i < dim; i++)
                fs << " " << X(i);
            if (dim == 2)
                fs << " 0";
            fs << "\n";
        }
        fs << "POLYS\n";
        count = 0;
        fs << "END\n";
        fs.close();
    }
};