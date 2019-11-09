#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>
#include <sys/stat.h>
#include <iostream>

template<class T, int dim>
class GridNode {
public:
    // MEMBERS
    T mass;
    Eigen::Matrix<T,dim,1> position;
    Eigen::Matrix<T,dim,1> velocity; 

    // CONSTRUCTORS
    GridNode(mass, position, velocity): mass(mass), position(potition), velocity(velocity) {} 
    GridNode() : mass(0), position(Eigen::Matrix<T,dim,1>(0, 0, 0), Eigen::Matrix<T,dim,1>(0, 0, 0)) {}
};

template<int dim>
class CartesianGrid {
public:
    // MEMBERS
    Eigen::Matrix<GridNode, dim, dim> grid;
    CartesianGrid() {}

    void clear() {
        this = CartesianGrid<dim>();
    }
};