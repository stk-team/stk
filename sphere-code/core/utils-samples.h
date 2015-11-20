#ifndef UTILSSAMPLES_H
#define UTILSSAMPLES_H

#include <vector>
#include <iostream>

inline void write_samples(std::vector<double> &data, int dim){

    int tempN = data.size() / double(dim);

    std::cerr << "Numpts: " << tempN << std::endl;

    for(int i=0; i < tempN; i++){
        for(int k=0; k < dim; k++){
            std::cout << data[dim*i+k] <<" ";
        }
        std::cout << std::endl;
    }
}

#endif // UTILSSAMPLES_H
