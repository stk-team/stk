/*** Created by Gurprit Singh ***/
#include "read-pointset.h"
#include <iostream>
#include <fstream>
#include <limits>

template <typename T>
void read_pointsetnD(std::string filename, std::vector<T> &vec, int dims){
    std::ifstream file(filename.c_str());
    if(dims == 1){
        if (file) {
           T x=0;

            while (file >> x){
                vec.push_back(x);
            }
        }
    }
    else if(dims==2){
        if (file) {
            T x[2];
            x[0]=0,x[1]=0;
            while (file >> x[0] >> x[1]){
                vec.push_back(x[0]);
                vec.push_back(x[1]);
            }
        }
    }
    else if(dims==3){
        if (file) {
            T x[3];
            x[0]=0,x[1]=0,x[2]=0;
            while (file >> x[0] >> x[1] >> x[2]){
                vec.push_back(x[0]);
                vec.push_back(x[1]);
                vec.push_back(x[2]);
            }
        }
    }
    else if(dims==4){
        if (file) {
            T x[4];
            x[0]=0,x[1]=0,x[2]=0,x[3]=0;
            while (file >> x[0] >> x[1] >> x[2] >> x[3]){
                vec.push_back(x[0]);
                vec.push_back(x[1]);
                vec.push_back(x[2]);
                vec.push_back(x[3]);
            }
        }
    }
    else if(dims==5){
        if (file) {
            T x[5];
            x[0]=0,x[1]=0,x[2]=0,x[3]=0,x[4]=0;
            while (file >> x[0] >> x[1] >> x[2] >> x[3] >> x[4]){
                vec.push_back(x[0]);
                vec.push_back(x[1]);
                vec.push_back(x[2]);
                vec.push_back(x[3]);
                vec.push_back(x[4]);
            }
        }
    }
    else{
        std::cerr << " File doesn't exist! " << filename << std::endl;
        exit(-1);
    }
    //std::cerr <<"Number of points in the data: " << vec.size() / dims << std::endl;
}

template void read_pointsetnD(std::string filename, std::vector<float>& vec, int dim);
template void read_pointsetnD(std::string filename, std::vector<double>& vec, int dim);
template void read_pointsetnD(std::string filename, std::vector<unsigned long>& vec, int dim);



