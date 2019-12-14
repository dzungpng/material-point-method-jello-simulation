#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h> /* strtof */
#include <cstring>
using namespace std;

// Class to parse a poly file and extract position information from it
template<class T, int dim>
class PolyParser {
public:
    using TV = Eigen::Matrix<T,dim,1>;
    std::vector<TV> x;

    void parseXYZPosition(const string str, std::vector<TV> &x) 
    { 
        char cstr[str.size() + 1];
        strcpy(cstr, str.c_str());
        char* xEnd;
        char* yEnd;
        float px = strtof (cstr, &xEnd);
        float py = strtof (xEnd, &yEnd);
        float pz = strtof (yEnd, NULL);
        x.push_back(TV(px, py, pz));
    } 

    void parsePoints(string filePath) {
        ifstream file(filePath);
        string line;
        if (file.is_open()) 
        {
            while(getline(file, line))
            {
                // position vector
                if (line.rfind("POLYS", 0) != 0) {
                    string position = line.substr(line.find(":") + 1);
                    parseXYZPosition(position, x);
                } else {
                    break;
                }
            }
            file.close();
        }
        else cout << "Unable to open file.";
    }

    PolyParser(string filePath) {
        parsePoints(filePath);
    }

};