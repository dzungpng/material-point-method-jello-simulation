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
    T theta_t;
    T theta_c;
    T epsilon;

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

    void addGeometry(std::vector<T> m_new,
                     std::vector<TV> x_new,
                     std::vector<TV> v_new,
                     std::vector<Mat> Fp_new,
                     std::vector<Mat> Fe_new,
                     std::vector<Mat> F_new,
                     std::vector<T> Vp0_new,) 
    {
        m.insert(m.end(), m_new.begin(), m_new.end());
        x.insert(x.end(), x_new.begin(), x_new.end());
        v.insert(v.end(), v_new.begin(), v_new.end());
        F.insert(F.end(), F_new.begin(), F_new.end());
        Fp.insert(Fp.end(), Fp_new.begin(), Fp_new.end());
        Fe.insert(Fe.end(), Fe_new.begin(), Fe_new.end());
        Vp0.insert(Vp0.end(), Vp0_new.begin(), Vp0_new.end());
    }    
    
};
