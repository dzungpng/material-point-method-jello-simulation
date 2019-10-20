#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>     /* strtof */
#include <cstring>
#include <unordered_set>
using namespace std;


struct hash_pair { 
    template <class T1, class T2> 
    size_t operator()(const pair<T1, T2>& p) const
    { 
        auto hash1 = hash<T1>{}(p.first); 
        auto hash2 = hash<T2>{}(p.second); 
        return hash1 ^ hash2; 
    } 
}; 

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

    void parseSegments(const string str,
                       unordered_set<pair<int, int>, hash_pair>& segment_set,
                       vector<TV>& x,
                       const T density) 
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

        pair<int, int> seg1;
        pair<int, int> seg2;
        pair<int, int> seg3;
        pair<int, int> seg4;
        pair<int, int> seg5;
        pair<int, int> seg6;

        // seg1 << a,b;
        // seg2 << a,c;
        // seg3 << a,d;
        // seg4 << b,c;
        // seg5 << b,d;
        // seg6 << c,d;

        seg1.first = a;
        seg1.second = b;

        seg2.first = a;
        seg2.second = c;
        
        seg3.first = a;
        seg3.second = d;

        seg4.first = b;
        seg4.second = c;

        seg5.first = b;
        seg5.second = d;

        seg6.first = c;
        seg6.second = d;

        segment_set.insert(seg1);
        segment_set.insert(seg2);
        segment_set.insert(seg3);
        segment_set.insert(seg4);
        segment_set.insert(seg5);
        segment_set.insert(seg6);

        // T rest_length_ab = (x[b] - x[a]).norm();
        // T rest_length_ac = (x[c] - x[a]).norm();
        // T rest_length_ad = (x[d] - x[a]).norm();
        // T rest_length_bc = (x[b] - x[c]).norm();
        // T rest_length_bd = (x[b] - x[d]).norm();
        // T rest_length_cd = (x[d] - x[c]).norm();

        // rest_length.push_back(rest_length_ab);
        // rest_length.push_back(rest_length_ac);
        // rest_length.push_back(rest_length_ad);
        // rest_length.push_back(rest_length_bc);
        // rest_length.push_back(rest_length_bd);
        // rest_length.push_back(rest_length_cd);

        // m[a] += rest_length_ab*density/2;
        // m[a] += rest_length_ac*density/2;
        // m[a] += rest_length_ad*density/2;

        // m[b] += rest_length_ab*density/2;
        // m[b] += rest_length_bc*density/2;
        // m[b] += rest_length_bd*density/2;

        // m[c] += rest_length_ac*density/2;
        // m[c] += rest_length_bc*density/2;
        // m[c] += rest_length_cd*density/2;

        // m[d] += rest_length_ad*density/2;
        // m[d] += rest_length_bd*density/2;
        // m[d] += rest_length_cd*density/2;

        // segments.push_back(seg1);
        // segments.push_back(seg2);
        // segments.push_back(seg3);
        // segments.push_back(seg4);
        // segments.push_back(seg5);
        // segments.push_back(seg6);
    } 

    void buildSegments(const string vtkFilePath,
                       int& N,
                       vector<TV>& x,
                       unordered_set<pair<int, int>, hash_pair>& segment_set,
                       const T density,
                       vector<TV>& v){
        ifstream file(vtkFilePath);
        string line;
        if (file.is_open()) 
        {
            while(getline(file, line))
            {
                // position vector
                if (line.rfind("POINTS", 0) == 0) {
                    string position;
                    // Getting number of vertices
                    string::size_type s; 
                    N = stoi(first_numberstring(line),&s);
                    x.resize(N);
                    v.resize(N);
                    int id = 0;
                    while(getline(file, position) && !position.empty()) {
                        parseXYZPosition(position, x[id]);
                        v[id] = TV::Zero();
                        id++;
                    }
                } 

                // segment vector
                if (line.rfind("CELLS", 0) == 0) {
                    string tet;
                    while(getline(file, tet) && !tet.empty()) {
                        parseSegments(
                            tet,
                            segment_set,
                            x,
                            density
                        );
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
    vector<T> m; // built in constructor
    vector<TV> x; // built in parser
    vector<TV> v; // built in sim

    vector<Eigen::Matrix<int,2,1>> segments; // built in constructor
    vector<T> rest_length; // built in parser
    vector<bool> node_is_fixed;

    unordered_set<pair<int, int>, hash_pair> segment_set; 

    int N; // Number of vertices.
    T density = 2.f; // 1.4
    VTKParser<T, dim> parser;
    
    // CONSTRUCTORS
    SegmentMesh(string vtkFilePath) {
        parser.buildSegments(
            vtkFilePath,
            N,
            x,
            segment_set,
            density,
            v
        );
        node_is_fixed.resize(N, false);

        m.resize(N);

        unordered_set<pair<int, int>, hash_pair>::iterator itr; 
        for (itr = segment_set.begin(); itr != segment_set.end(); itr++) {
            int idx1 = (*itr).first;
            int idx2 = (*itr).second;
            T len = (x[idx2]-x[idx1]).norm();
            rest_length.push_back(len);

            Eigen::Matrix<int,2,1> segment;
            segment << idx1, idx2;
            segments.push_back(segment);

            m[idx1] = len*density/2;
            m[idx2] = len*density/2;
        }
    }
};
