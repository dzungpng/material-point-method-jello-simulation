#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>
#include <sys/stat.h>
#include <iostream>
#include <map> 

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

    // cell width
    T dx = 0.02;

    // std::vector<T> mg;
    // std::vector<TV> vg;
    std::map<TV, T> mg; // mapping grid coords to mass
    std::map<TV, TV> vg; // mapping grid coords to velocity
    std::vector<TV> force; // mapping grid coords to force
 
    // grid dimensions
    TV res;
    
    CartesianGrid() : minCorner(TV::Zero()), maxCorner(TV::Ones()) {}
    CartesianGrid(TV minCorner,TV maxCorner) : minCorner(minCorner), maxCorner(maxCorner) {}

    void clear() {
        mg.clear();
        vg.clear();
        force.clear();
    }

};