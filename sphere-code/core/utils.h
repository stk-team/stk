#ifndef GARAGE_CORE_UTILS_H
#define GARAGE_CORE_UTILS_H

#if defined(__linux__)
#define GARAGE_IS_LINUX
#elif defined(__APPLE__)
  #define GARAGE_IS_APPLE
  #if !(defined(__i386__) || defined(__amd64__))
  #define GARAGE_IS_APPLE_PPC
  #else
  #define GARAGE_IS_APPLE_X86
  #endif
#endif
#include <math.h>       /* floor */
#include "mymacros.h"
#include <iostream>
#include <vector>
#include <stdlib.h>


void create_folders(std::string homedir, std::string &data, std::string &images, std::string &graphs);
void read_fibonacci_table(std::vector<unsigned long> &Table, const std::string &filename);

// Global Inline Functions
inline float Lerp(float t, float v1, float v2) {
    return (1.f - t) * v1 + t * v2;
}


inline float Clamp(float val, float low, float high) {
    if (val < low) return low;
    else if (val > high) return high;
    else return val;
}


inline int Clamp(int val, int low, int high) {
    if (val < low) return low;
    else if (val > high) return high;
    else return val;
}


inline int Mod(int a, int b) {
    int n = int(a/b);
    a -= n*b;
    if (a < 0) a += b;
    return a;
}

inline int Floor2Int(float val) {
    return (int)floorf(val);
}


inline int Round2Int(float val) {
    return Floor2Int(val + 0.5f);
}


inline int Float2Int(float val) {
    return (int)val;
}


inline int Ceil2Int(float val) {
    return (int)ceilf(val);
}

inline bool IsPowerOfTwo(int a){
    if(a==0 || a < 0)
        return false;
    else if(a%2 != 0)
        return false;
    else if(a/2 == 1)
        return true;
    else{
        a = a >> 1;
        IsPowerOfTwo(a);
    }
}

inline int sgn(double x){
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}

inline float radicalInverse(int n, int base){
    float val = 0;
    float invBase = 1. / float(base), invBi = invBase;
    while(n > 0){
        int d_i = n % base;
        val += d_i * invBi;
        n *= invBase;
        invBi *= invBase;
    }
    return val;
}

inline void addingOffsetForSinusoids(float* data_ptr, int width, int height){
    For(r, width){
        For(c, height){
            float fcol = data_ptr[r*width+c];
            fcol = (1.f + fcol)/2.0f;
            float clamped_col = Clamp(fcol, 0.0f, 1.0f);
            data_ptr[r*width+c] = clamped_col;
        }
    }
}

inline int getnthbit(const int &N, int nthbit){
    int mask =  1 << nthbit;
    int masked_n = N & mask;
    int thebit = masked_n >> nthbit;
    return thebit;
}

inline int getnbits(const int &N){
    int temp = N;
    int nbits=0;
    do{
        temp = temp >> 1;
        int mask =  1 << nbits;
        int masked_n = N & mask;

        if(masked_n != 0){
            std::cerr << "N is not of pow 2 !!!" << std::endl;
            exit(-2);
        }
        nbits++;
    }while(temp != 1);
    return nbits;
}

inline void paddedzerosN(std::string &s1, int fileNumber){

    if(fileNumber > 100000){
        std::cerr << "fileNumber is too big !!!" << std::endl;
        exit(-2);
    }

    do{
        s1 = "0"+s1;
    }while(s1.length() != 6);
}

#endif
