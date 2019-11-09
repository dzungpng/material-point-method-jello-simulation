#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>
#include <sys/stat.h>
#include <iostream>


template<class T, int dim>
class Point {
public:
    T mass;
    Eigen::Matrix<T,dim,1> position;
    Eigen::Matrix<T,dim,1> velocity;
};