#include "sampler-sphere.h"
#include "./../core/utils-sphere.h"
#include "./../domaintransform/domaintransform.h"

#include <cmath>
#include <iostream>
#include <set>
#include <unordered_set>

//Slow Dart -throwing version of Poisson Disk Sampling.
//For Fast version, look the implementation after this one!
std::vector<double> spherical_dart_throwing_samples(const int &N){

    double radius = 2.793/sqrt(N); //The magic number 2.763
    std::vector<double> buffer;
    double sx,sy,sz;
    uniformSampleSphere(drand48(), drand48(), &sx,&sy,&sz);
    buffer.push_back(sx);
    buffer.push_back(sy);
    buffer.push_back(sz);
    for(int i = 0; i < N*500; i++){
        double sx=0,sy=0,sz=0;
        uniformSampleSphere(drand48(), drand48(), &sx,&sy,&sz);
        bool conflict = false;
        for(int j = 0; j < buffer.size(); j+=3){
            double dx = sx - buffer[j+0];
            double dy = sy - buffer[j+1];
            double dz = sz - buffer[j+2];

            double euclid = sqrt(dx*dx+dy*dy+dz*dz);
            if(euclid < radius){
                conflict = true;
                break;
            }
        }
        if(!conflict){
            buffer.push_back(sx);
            buffer.push_back(sy);
            buffer.push_back(sz);
        }
        if(buffer.size() == (3*N))
            break;
    }
    //std::cerr << "# darts: " << buffer.size() / 3.0 << std::endl;

    int tempN = N;
    if(buffer.size() != (3*N)){
        tempN = buffer.size()/3.0;
        //sphereSamples.resize(tempN);
       //std::cerr << "# samples generated: " << tempN << std::endl;
    }
    //vec3 initval(0,0,0);
    //arr<vec3> sphereSamples(tempN, initval);
    std::vector<double> dartSamples;
    for(int i=0; i < tempN; i++){
        dartSamples.push_back(buffer[3*i+0]);
        dartSamples.push_back(buffer[3*i+1]);
        dartSamples.push_back(buffer[3*i+2]);
        //std::cout << buffer[3*i] << " " << buffer[3*i+1] << " " << buffer[3*i+2] << " " << 1 << std::endl;
    }
    return dartSamples;
}





// Fast Poisson Disk Samplings
void griditize(cell* grid, std::vector<unsigned int> &validcells, const double* sampleSpace, int sampleSpaceSize,
                       int gridSize, double grid_dx, int xmin){
    std::unordered_set<unsigned int> temp;
     double inv_grid_dx = 1.0/grid_dx;

    for(int i = 0; i < sampleSpaceSize; i++){

        int col =  (sampleSpace[3*i+0] - xmin) * inv_grid_dx;
        int row =  (sampleSpace[3*i+1] - xmin) * inv_grid_dx;
        int slice =  (sampleSpace[3*i+2] - xmin) * inv_grid_dx;

        if(row > gridSize-1 || col > gridSize-1 || slice > gridSize-1){
            std::cerr << row <<" " << col <<" " << slice << " : " << "OUT OF BOUND !!!" << std::endl;
            exit(-2);
        }
        int index = (row*gridSize+col) + slice*gridSize*gridSize;

        grid[index].c = col;
        grid[index].r = row;
        grid[index].s = slice;
        grid[index].cellid = index; //This should be randomly selected, for now it is useless
        grid[index].samples.push_back(sampleSpace[3*i+0]);
        grid[index].samples.push_back(sampleSpace[3*i+1]);
        grid[index].samples.push_back(sampleSpace[3*i+2]);
        grid[index].valid = true;

        if(temp.insert(index).second){
            validcells.push_back(index);
        }
    }
}

std::vector<double> spherical_poisson_disk_samples(const int &N, const double* sampleSpace, int sampleSpaceSize){

    //double tradius = 0.086;
    double tradius = 2.893/sqrt(N);
    double xmax = 1.0, xmin = -1.0;
    double grid_dx = (0.999)*tradius/std::sqrt(3);
    unsigned int gridSize = (unsigned int)std::ceil( (xmax-xmin) / grid_dx);
    std::vector<double> darts;
    //std::cerr << "gridSize: " << gridSize << std::endl;
    //exit(-2);
    //######################################################################

    //Create a 3D grid of cell size r/sqrt(3);
    cell *grid = new cell[gridSize*gridSize*gridSize] ();
    std::vector<unsigned int> validcells;

    griditize(grid, validcells, sampleSpace, sampleSpaceSize, gridSize, grid_dx, xmin);

    //######################################################################

    //Pre-compute random ordering for the selection of cells
    std::unordered_set <unsigned int> randomIndices;
    std::vector<unsigned int> indices;
    int b = validcells.size();
    int a = 0;
    for(;;){
        if(indices.size() == validcells.size())
            break;
        int p = (b-a)*drand48() + a;
        if(randomIndices.insert(p).second)
            indices.push_back(p);
    }

    //######################################################################
    //Go through all valid indices in pre-computed randomized order.
    for(int ii=0; ii<indices.size(); ii++){
        //fprintf(stderr,"\r index (%d) ",ii);

        //std::cout << ii  <<" " << validcells[indices[ii]] << " "<< validcells[ii] << std::endl;
        int index = validcells[indices[ii]];
        int row = grid[index].r;
        int col = grid[index].c;
        int slice = grid[index].s;
        //std::cout << "#### " << col <<" " << row <<" " << slice << std::endl;
        if(index != grid[index].cellid){
            std::cerr << "Not the right grid cell !!!" << std::endl;
            exit(-2);
        }

        //##################################################################
        ///For the first sample, choose any random sample from any valid cell,
        ///Store that sample and remove all other samples from that grid cell.
        if(ii == 0){
            int isamp = (grid[index].samples.size()/3.0) * drand48();
            darts.push_back(grid[index].samples[3*isamp+0]);
            darts.push_back(grid[index].samples[3*isamp+1]);
            darts.push_back(grid[index].samples[3*isamp+2]);

            grid[index].visited = true;
            grid[index].samples.clear();
            grid[index].samples.push_back(darts[0]);
            grid[index].samples.push_back(darts[1]);
            grid[index].samples.push_back(darts[2]);
        }
        //######################################################################
        else{  // go through each sample and check whether is it far enough !
            ///Check if the next sample is far away from already visited valid grid cells
            /// If yes, then just take that sample as a valid sample.
            bool is_nbh_visited = false;
            for(int rr = row-2; rr <= row+2; rr++){
                for(int cc = col-2; cc <= col+2; cc++){
                    for(int ss = slice-2; ss <= slice+2; ss++){
                        if(rr < 0 || rr > gridSize-1)
                            continue;
                        if(cc < 0 || cc > gridSize-1)
                            continue;
                        if(ss < 0 || ss >  gridSize-1)
                            continue;

                        int tindex = (rr*gridSize+cc) + (ss*gridSize*gridSize);

                        if(grid[tindex].visited && grid[tindex].valid){
                            is_nbh_visited = true;
                            break;
                        }
                    }
                }
            }
            if(!is_nbh_visited){
                int isamp = (grid[index].samples.size()/3.0) * drand48();
                double x = grid[index].samples[3*isamp+0];
                double y = grid[index].samples[3*isamp+1];
                double z = grid[index].samples[3*isamp+2];
                darts.push_back(x);
                darts.push_back(y);
                darts.push_back(z);

                grid[index].visited = true;
                grid[index].samples.clear();
                grid[index].samples.push_back(x);
                grid[index].samples.push_back(y);
                grid[index].samples.push_back(z);
            } /*******************/
            else if(is_nbh_visited){   //If atleast one of the visited cell is the nbh of this new valid cell
                std::vector<double> valid_nbhs;
                for(int rr = row-2; rr <= row+2; rr++){
                    for(int cc = col-2; cc <= col+2; cc++){
                        for(int ss = slice-2; ss <= slice+2; ss++){
                            if(rr < 0 || rr > gridSize-1)
                                continue;
                            if(cc < 0 || cc > gridSize-1)
                                continue;
                            if(ss < 0 || ss >  gridSize-1)
                                continue;

                            int tindex = (rr*gridSize+cc) + (ss*gridSize*gridSize);

                            if(grid[tindex].valid && grid[tindex].visited && tindex != index ){
                                for(int i = 0; i < grid[tindex].samples.size(); i++)
                                    valid_nbhs.push_back(grid[tindex].samples[i]);
                            }
                        }
                    }
                }
                //Go through all the valid neighbors
                int numTrialSamples = 30;
                if(numTrialSamples > (grid[index].samples.size()/3.0))
                    numTrialSamples = grid[index].samples.size()/3.0;
                bool found_sample = false;
                std::unordered_set<unsigned int> tSamples;
                for(;;){
                    int nspp = grid[index].samples.size()/3.0;
                    unsigned int isamp = nspp * drand48();
                    tSamples.insert(isamp);
                    if(tSamples.size() == numTrialSamples)
                        break;
                }
                for(auto t = tSamples.begin(); t != tSamples.end(); ++t){
                    int irand = *t;
                    bool conflict = false;
                    for(int j = 0; j < valid_nbhs.size(); j+=3){
                        double dx = grid[index].samples[3*irand+0] - valid_nbhs[j+0];
                        double dy = grid[index].samples[3*irand+1] - valid_nbhs[j+1];
                        double dz = grid[index].samples[3*irand+2] - valid_nbhs[j+2];

                        double dist = sqrt(dx*dx + dy*dy + dz*dz);
                        if(dist < tradius){
                            conflict = true;
                            break;
                        }
                    }
                    if(!conflict){
                        found_sample = true;
                        double x = grid[index].samples[3*irand+0];
                        double y = grid[index].samples[3*irand+1];
                        double z = grid[index].samples[3*irand+2];
                        darts.push_back(x);
                        darts.push_back(y);
                        darts.push_back(z);

                        grid[index].visited = true;
                        grid[index].samples.clear();
                        grid[index].samples.push_back(x);
                        grid[index].samples.push_back(y);
                        grid[index].samples.push_back(z);
                        break;
                    }
                }
                if(!found_sample){
                    grid[index].samples.clear();
                    grid[index].valid = false;
                    grid[index].visited = true;
                }
            }
            /**********************/
        }
    }

    int numDarts = darts.size()/3.0;
    //vec3 initval(0,0,0);
    //arr<vec3> poissonSamples(numDarts, initval);
    std::vector<double> poissonSamples;
    for(int i=0; i < numDarts; i++){
        poissonSamples.push_back(darts[3*i+0]);
        poissonSamples.push_back(darts[3*i+1]);
        poissonSamples.push_back(darts[3*i+2]);
    }

    delete [] grid;
    return poissonSamples;
}


