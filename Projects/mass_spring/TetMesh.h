#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>     /* strtof */
#include <cstring>
using namespace std;

/*
 * Class that parses a vtk file and construct a mapping representing segments
 * of a segment mesh.
 */
template<class T, int dim>
class VTKParser {
public:
    using TV = Eigen::Matrix<T,dim,1>;

    string first_numberstring(std::string const & str)
    {
        size_t const n = str.find_first_of("0123456789");
        if (n != string::npos)
        {
            size_t const m = str.find_first_not_of("0123456789", n);
            return str.substr(n, m != string::npos ? m-n : m);
        }
        return string();
    }

    void parseXYZPosition(string str, TV& position) 
    { 
        char cstr[str.size() + 1];
        strcpy(cstr, str.c_str());
        char* xEnd;
        char* yEnd;
        float x, y, z;
        x = strtof (cstr, &xEnd);
        y = strtof (xEnd, &yEnd);
        z = strtof (yEnd, NULL);
        position(0) = x;
        position(1) = y;
        position(2) = z;
    } 

    std::vector<Eigen::Matrix<int,2,1>> buildSegments(string vtkFilePath,
                                                      int& N,
                                                      vector<TV>& x){
        ifstream file(vtkFilePath);
        string line;
        if (file.is_open()) 
        {
            while(getline(file, line))
            {
                if (line.rfind("POINTS", 0) == 0) {
                    // Adding the positions of each point into x vector
                    string position;
                    N = int(first_numberstring(line));
                    x.resize(N*N);
                    int id = 0;
                    while(getline(file, position) && !position.empty()) {
                        parseXYZPosition(position, x[id]);
                        id+=1;
                    }
                }
            }
            file.close();
        }
        else cout << "Unable to open file.";
        cout << "Total vertices read: " << N << '\n';
        return std::vector<Eigen::Matrix<int,2,1>>(); 
    }
}; 

/* 
 * Class that contains a hashmap representing the ordering
 * of vertices parsed from a TetMesh. 
*/
template<class T, int dim>
class SegmentMesh {
public:
    using TV = Eigen::Matrix<T,dim,1>;

    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> v;

    std::vector<Eigen::Matrix<int,2,1>> segments;
    std::vector<T> rest_length;

    int N; // Number of vertices.
};
