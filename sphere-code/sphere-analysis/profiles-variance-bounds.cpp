#include "sphere-analysis.h"

double constantProfile(int l, double a){
    return a;
}

double linearProfile(int l, double limit, double alpha, int n){
    if(l < (alpha*sqrt(n)))
        return limit * l/(alpha*sqrt(n));
    else
        return limit;
}

double quadraticProfile(int l, double limit, double alpha, int n){
    if(l < (alpha*sqrt(n)))
        return limit * l * l/(alpha*alpha*n);
    else
        return limit;
}

double linearProfileLower(int l, double slope, double limit, double scale){
    if(scale *l*slope > limit)
        return limit;
    else
        return scale * l*slope;
}

