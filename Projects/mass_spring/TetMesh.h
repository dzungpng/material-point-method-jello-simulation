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

    void parseXYZPosition(const string str, TV& position) 
    { 
        char cstr[str.size() + 1];
        strcpy(cstr, str.c_str());
        char* xEnd;
        char* yEnd;
        float x = strtof (cstr, &xEnd);
        float y = strtof (xEnd, &yEnd);
        float z = strtof (yEnd, NULL);
        position(0) = x;
        position(1) = y;
        position(2) = z;
    } 

    void parseSegments(const string str, vector<Eigen::Matrix<int,2,1>>& segments) 
    { 
        char cstr[str.size() + 1];
        strcpy(cstr, str.c_str());
        char* End;
        char* aEnd;
        char* bEnd;
        char* cEnd;
        strtof (cstr, &End);
        float a = strtof (End, &aEnd);
        float b = strtof (aEnd, &bEnd);
        float c = strtof (bEnd, &cEnd);
        float d = strtof (cEnd, NULL);

        Eigen::Matrix<int,2,1> seg1;
        Eigen::Matrix<int,2,1> seg2;
        Eigen::Matrix<int,2,1> seg3;

        seg1 << b,a;
        seg2 << c,a;
        seg3 << d,a;

        segments.push_back(seg1);
        segments.push_back(seg2);
        segments.push_back(seg3);
    } 

    void buildSegments(const string vtkFilePath,
                       int& N, vector<TV>& x,
                       vector<Eigen::Matrix<int,2,1>>& segments){
        ifstream file(vtkFilePath);
        string line;
        if (file.is_open()) 
        {
            while(getline(file, line))
            {
                if (line.rfind("POINTS", 0) == 0) {
                    // Adding the pox[i](1) << ", " << test_x[i](2) <<sitions of each point into x vector
                    string position;
                    // Getting number of vertices
                    string::size_type s; 
                    N = stoi(first_numberstring(line),&s);
                    x.resize(N);
                    int id = 0;
                    while(getline(file, position) && !position.empty()) {
                        parseXYZPosition(position, x[id]);
                        id++;
                    }
                } 
                if (line.rfind("CELLS", 0) == 0) {
                    string tet;
                    string::size_type s;
                    while(getline(file, tet) && !tet.empty()) {
                        parseSegments(tet, segments);
                    }
                }
            }
            file.close();
        }
        else cout << "Unable to open file.";
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

    // MEMBERS 
    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> v;

    std::vector<Eigen::Matrix<int,2,1>> segments;
    std::vector<T> rest_length;

    int N; // Number of vertices.

    // CONSTRUCTORS
    SegmentMesh(string vtkFilePath) {
        VTKParser<T, dim> parser;
        parser.buildSegments("../../torus.vtk", N, x, segments);
    }
};
