/* An ANN class. */


#ifndef ANN_H
#define ANN_H

#include <string>
#include <mkl.h>
#include "MVOperations.h"
#include "VectorOperations.h"
#include "SDNN_data.h"
#include <cmath>
#include <algorithm>
#include <fstream>


constexpr int s = 7;

class ANN {

double *W1;
double *B1;
double *W2;
double *B2;
double *W3;
double *B3;
double *W4;
double *B4;

const double alpha;

public :

    ANN(std::string W1_filename, std::string W2_filename, 
        std::string W3_filename, std::string W4_filename, 
        std::string B1_filename, std::string B2_filename, 
        std::string B3_filename, std::string B4_filename, double _alpha);   

    ANN(double _alpha) : alpha{_alpha}{set_NN_weights();}

    void set_NN_weights(std::string W1_file, std::string B1_file, 
        std::string W2_file, std::string B2_file, std::string W3_file, 
        std::string B3_file, std::string W4_file, std::string B4_file);

    void set_NN_weights();

    ANN(const ANN & ann) :  alpha{ann.alpha} {set_NN_weights();}

    ANN & operator= (const ANN &ann);

    virtual ~ANN() {free_mem();}

    double getAlpha() const {return alpha;}

    void fwdClassif(int *tau, const double *regStencils, int N0);

private :

    inline void elu(const int N, double* input, const double alpha) const;

    void softmax(int N, double* input) const;

    void readMatrix(double* A, std::string filename_A);

    void free_mem();

};


#endif