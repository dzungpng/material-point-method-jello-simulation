#ifndef SAMPLEDMESH_H
#define SAMPLEDMESH_H

#include <Eigen/Core>
#include <Eigen/Dense>

// Geometry for storing sampled particles
template<class T, int dim>
class SampledMesh {
public:
    using TV = Eigen::Matrix<T,dim,1>;
    using Mat = Eigen::Matrix<T, dim, dim>;

    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> v;
    std::vector<Mat> Fp;
    std::vector<T> Vp0;

    SampledMesh(const std::vector<T> m, const std::vector<TV> x, const std::vector<TV> v, const std::vector<Mat> Fp, const std::vector<T> Vp0) : m(m), x(x), v(v), Fp(Fp), Vp0(Vp0) {}
};

#endif // SAMPLEDMESH_H
