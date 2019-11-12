#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>
#include <sys/stat.h>
#include <iostream>
#include <map>

template<class T, int dim>
class CartesianGrid {
public:
    using TV = Eigen::Matrix<T,dim,1>;

    // MEMBERS
    TV minCorner;
    TV maxCorner;

    // cell width
    T cellWidth; // cell width

    std::vector<T> mg;
    std::vector<TV> vg;
    std::vector<TV> xg;
    std::vector<TV> force; // mapping grid coords to force
    std::vector<T> active_nodes;

    TV res; // resolution of the grid
    T nCells; // number of cells in the grid
    
    CartesianGrid() : minCorner(TV::Zero()), maxCorner(TV::Ones()) {}
    CartesianGrid(TV minCorner,TV maxCorner) : minCorner(minCorner), maxCorner(maxCorner) {}

    void clear() {
        mg.clear();
        vg.clear();
        xg.clear();
        force.clear();
        active_nodes.clear();
        mg.resize(nCells);
        vg.resize(nCells);
        xg.resize(nCells);
        force.resize(nCells);
        active_nodes.resize(nCells);
    }
};
