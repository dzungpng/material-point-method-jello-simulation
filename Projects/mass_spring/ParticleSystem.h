#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <fstream>

#include "Geometry/SampledMesh.h"

template<class T, int dim>
class MPM{
public:
    using TV = Eigen::Matrix<T,dim,1>;
    using Mat = Eigen::Matrix<T, dim, dim>;

    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> v;
    std::vector<Mat> Fp;
    std::vector<T> Vp0;

    T lambda;
    T mu;
    T theta_t;
    T theta_c;
    T epsilon;

    MPM() {}
    /**
     * dumping poly files. 
     * @param  {std::string} filename : 
     */
    void dumpPoly(std::string filename)
    {
        std::ofstream fs;
        fs.open(filename);
        fs << "POINTS\n";
        int count = 0;
        for (auto X : x) {
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

    void addGeometry(const SampledMesh<T, dim> *mesh) 
    /**
     * This function is called inside of run every x seconds to generate new geometry (if wanted)
     * @param mesh: 
     */
    {
        m.insert(m.end(), mesh->m.begin(), mesh->m.end());
        x.insert(x.end(), mesh->x.begin(), mesh->x.end());
        v.insert(v.end(), mesh->v.begin(), mesh->v.end());
        Fp.insert(Fp.end(), mesh->Fp.begin(), mesh->Fp.end());
        Vp0.insert(Vp0.end(), mesh->Vp0.begin(), mesh->Vp0.end());
    }    
    
};
