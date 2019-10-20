#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>
#include <sys/stat.h>
#include <iostream>

template<class T, int dim>
class Cloth {
public:
    // using TV = Eigen::Matrix<T,dim,1>;

    std::vector<T> m;
    std::vector<Eigen::Matrix<T,dim,1>> x;
    std::vector<Eigen::Matrix<T,dim,1>> v;
    std::vector<bool> node_is_fixed;

    std::vector<Eigen::Matrix<int,2,1>> segments;
    std::vector<T> rest_length;
    int N = 64;
    int N_points = N*N;
    T dx = (T)1/(N-1);

    Cloth() {
        m.resize(N_points);
        x.resize(N_points);
        v.resize(N_points);

        for(int i=0; i<N; i++){
            for(int j=0; j<N; j++){
                int id = i*N+j;
                m[id] = (T)1/N_points;
                x[id](0) = (i-1)*dx;
                x[id](2) = (j-1)*dx;
                x[id](1) = 1;
                v[id] = Eigen::Matrix<T,dim,1>::Zero();
            }
        }
        node_is_fixed.resize(N_points,false);

        // structure
        for(int i=0; i<N-1; i++){
            for(int j=0; j<N; j++){
                Eigen::Matrix<int,2,1> seg;
                int p = i*N+j, q=(i+1)*N+j;
                seg << p,q;
                segments.push_back(seg);
                rest_length.push_back((x[p]-x[q]).norm());
            }
        }
        for(int i=0; i<N; i++){
            for(int j=0; j<N-1; j++){
                Eigen::Matrix<int,2,1> seg;
                int p = i*N+j, q=i*N+j+1;
                seg << p,q;
                segments.push_back(seg);
                rest_length.push_back((x[p]-x[q]).norm());
            }
        }

        // shear
        for(int i=0; i<N-1; i++){
            for(int j=0; j<N-1; j++){
                Eigen::Matrix<int,2,1> seg;
                int p = i*N+j, q=(i+1)*N+j+1;
                seg << p,q;
                segments.push_back(seg);
                rest_length.push_back((x[p]-x[q]).norm());
            }
        }
        for(int i=0; i<N-1; i++){
            for(int j=0; j<N-1; j++){
                Eigen::Matrix<int,2,1> seg;
                int p = (i+1)*N+j, q=i*N+j+1;
                seg << p,q;
                segments.push_back(seg);
                rest_length.push_back((x[p]-x[q]).norm());
            }
        }

        //bending
        for(int i=0; i<N-2; i++){
            for(int j=0; j<N; j++){
                Eigen::Matrix<int,2,1> seg;
                int p = i*N+j, q=(i+2)*N+j;
                seg << p,q;
                segments.push_back(seg);
                rest_length.push_back((x[p]-x[q]).norm());
            }
        }
        for(int i=0; i<N; i++){
            for(int j=0; j<N-2; j++){
                Eigen::Matrix<int,2,1> seg;
                int p = i*N+j, q=i*N+j+2;
                seg << p,q;
                segments.push_back(seg);
                rest_length.push_back((x[p]-x[q]).norm());
            }
        }
    }
};