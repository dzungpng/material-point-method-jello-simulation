#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <fstream>

template<class T, int dim>
class MPM{
public:
    using TV = Eigen::Matrix<T,dim,1>;
    using Mat = Eigen::Matrix<T, dim, dim>;

    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> v;

    std::vector<Mat> Fp;
    std::vector<Mat> Fe;
    std::vector<Mat> F;

    std::vector<T> Vp0;

    T lambda;
    T mu;
    T zeta;
    
    MPM() {}

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
    
};
