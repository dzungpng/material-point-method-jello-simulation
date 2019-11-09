#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>
#include <sys/stat.h>
#include <iostream>

// template<class T, int dim>
// class GridNode {
// public:
//     // MEMBERS
//     T mass;
//     using TV = Eigen::Matrix<T,dim,1>;

//     TV position;
//     TV velocity; 

//     // CONSTRUCTORS
//     GridNode(T mass, TV position, TV velocity): mass(mass), position(position), velocity(velocity) {} 
//     GridNode() : mass(0), position(TV(0, 0, 0), TV(0, 0, 0)) {}
// };

template<class T, int dim>
class CartesianGrid {
public:
    using TV = Eigen::Matrix<T,dim,1>;

    // MEMBERS
    TV minCorner;
    TV maxCorner;
    T dx = 0.02;

    CartesianGrid() : minCorner(TV::Zero()), maxCorner(TV::Ones()) {}
    CartesianGrid(TV minCorner,TV maxCorner) : minCorner(minCorner), maxCorner(maxCorner) {}

    void clear() {
        return;
    }
};