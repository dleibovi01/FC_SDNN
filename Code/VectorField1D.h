/* 1D vector field */

#ifndef VECTORFIELD1D_H
#define VECTORFIELD1D_H

#include <vector>
#include <iostream>
#include <algorithm>
// #include "SpatDiffScheme.h"
#include "mkl.h"



class VectorField1D{

std::vector<double*>  flow;
int unknowns;
int length;

public:
    //VectorField1D(Node1D* nodes);
    VectorField1D();
    
    VectorField1D(const std::vector<double*> &field, int n);
    
    // Copy constructor
    VectorField1D(const VectorField1D &field);

    // Move constructor
    VectorField1D(VectorField1D &&field);

    VectorField1D(int _unknowns, int _length);

    // Copy assignment
    VectorField1D & operator= (const VectorField1D &v);

    // Move assignment
    VectorField1D & operator= (VectorField1D &&field);

    VectorField1D & operator+= (const VectorField1D &v);

    VectorField1D & operator-= (const VectorField1D &v);

    VectorField1D & operator*= (const VectorField1D &v);

    VectorField1D & operator*= (double *v);

    VectorField1D & operator*= (double alpha);

    VectorField1D & operator/= (const VectorField1D &v);

    VectorField1D & operator>>= (const VectorField1D &v);

    VectorField1D & operator/= (double alpha){return (*this *= (1.0 / alpha));} 


    ~VectorField1D();
    
    int getUnknowns() const {return unknowns;}

    int getLength() const {return length;}

    double* getField(int n) const {return flow[n];}

    void setField(int n, int i, const double* flow_values)
    {
        std::copy(flow_values, flow_values + i, flow[n]);
    }

    void setFieldValue(int i, int j, double flow_value)
        {(flow[i])[j] = flow_value;}

    double getFieldValue(int i, int j) const {return (flow[i])[j];}        

    void setUnknows(int u){unknowns = u;}

    void setLength(int l){length = l;}

    std::vector<double *> getFlow() const{return flow;}

    void circshift (int i);   

    VectorField1D extract(std::vector<int> indices) const;

    void setFields(const VectorField1D &v,  std::vector<int> indices);

    void linComb(std::vector<double> coeffs, 
        std::vector<std::vector<int> > extractions, std::vector<int> target);

};

inline VectorField1D operator+ (const VectorField1D &a, const VectorField1D &b)
    {return VectorField1D{a}+=b;}

inline VectorField1D operator- (const VectorField1D &a, const VectorField1D &b)
    {return VectorField1D{a}-=b;}

inline VectorField1D operator* (const VectorField1D &a, const VectorField1D &b)
    {return VectorField1D{a}*=b;}    

inline VectorField1D operator* (const VectorField1D &a, double alpha)
    {return VectorField1D{a}*=alpha;}     

inline VectorField1D operator* (double alpha, const VectorField1D &a)
    {return VectorField1D{a}*alpha;}

inline VectorField1D operator/ (const VectorField1D &a, double alpha)
    {return VectorField1D{a}/=alpha;}    

VectorField1D circshift(const VectorField1D &v, int j);

// template<typename Sp_Diff>
// void diffLinComb(VectorField1D* v, 
//     const std::vector<std::vector<int> > &extractions, 
//     const std::vector<int> &target, const std::vector<double> & coeffs,
//     const VectorField1D &flux, const Sp_Diff &sp)
// {
//     // extractions[ncoeffs, unknowns]
//     // target [unknowns]
//     int ncoeffs = coeffs.size();
//     int N = v->getLength();
//     double data[v->getLength()];
//     int unknowns = v->getUnknowns();
//     for(int i = 0; i < length; i++)
//     {
//         data[i] = 0.0;
//     }
//     for(int i = 0; i < extractions.size(); i++) 
//     {
//         sp.diff(flux->getField(i), data);
//         cblas_dscal(N, coeffs[ncoeffs - 1], data, 1);
//         for(int j = 0; j < ncoeffs - 1; j++)
//         {
//             cblas_daxpy(N, coeffs[j], v->getField(extractions[j][i]), 1, data, 1);
//         }
//         v->setField(target[i], N, data);
//     }
// };
      

#endif 