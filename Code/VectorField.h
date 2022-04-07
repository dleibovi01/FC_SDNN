/* Vector field */

#ifndef VECTORFIELD_H
#define VECTORFIELD_H

#include <vector>
#include <iostream>
#include <algorithm>
#include "mkl.h"



class VectorField{

std::vector<double*>  flow;
int unknowns;
int length;

public:
    //VectorField(Node1D* nodes);
    VectorField();
    
    VectorField(const std::vector<double*> &field, int n);
    
    // Copy constructor
    VectorField(const VectorField &field);

    // Move constructor
    VectorField(VectorField &&field);

    VectorField(int _unknowns, int _length);

    // Copy assignment
    VectorField & operator= (const VectorField &v);

    // Move assignment
    VectorField & operator= (VectorField &&field);

    VectorField & operator+= (const VectorField &v);

    VectorField & operator+= (const double a);

    VectorField & operator-= (const VectorField &v);

    VectorField & operator*= (const VectorField &v);

    VectorField & operator*= (double *v);

    VectorField & operator*= (double alpha);

    VectorField & operator/= (const VectorField &v);

    VectorField & operator>>= (const VectorField &v);

    VectorField & operator/= (double alpha){return (*this *= (1.0 / alpha));} 

    VectorField & sqr();

    VectorField & abs();


    ~VectorField();
    
    int getUnknowns() const {return unknowns;}

    int getLength() const {return length;}

    double* getField(int n) const {return flow[n];}

    void setField(int n, int i, const double* flow_values)
    {
        std::copy(flow_values, flow_values + i, flow[n]);
    }

    void setFieldValue(int i, int j, double flow_value)
        {(flow[i])[j] = flow_value;}

    void setFlow(const std::vector<double*> flow_ext);

    double getFieldValue(int i, int j) const {return (flow[i])[j];}        

    void setUnknows(int u){unknowns = u;}

    void setLength(int l){length = l;}

    std::vector<double *> getFlow() const{return flow;}

    void circshift (int i);   

    VectorField extract(std::vector<int> indices) const;

    void setFields(const VectorField &v,  std::vector<int> indices);

    void linComb(std::vector<double> coeffs, 
        std::vector<std::vector<int> > extractions, std::vector<int> target);

};

inline VectorField operator+ (const VectorField &a, const VectorField &b)
    {return VectorField{a}+=b;}

inline VectorField operator+ (const VectorField &a, const double b)
    {return VectorField{a}+=b;}

inline VectorField operator+ (const double b, const VectorField &a)
    {return VectorField{a}+=b;}

inline VectorField operator- (const VectorField &a, const VectorField &b)
    {return VectorField{a}-=b;}

inline VectorField operator* (const VectorField &a, const VectorField &b)
    {return VectorField{a}*=b;}    

inline VectorField operator* (const VectorField &a, double alpha)
    {return VectorField{a}*=alpha;}     

inline VectorField operator* (double alpha, const VectorField &a)
    {return VectorField{a}*alpha;}

inline VectorField operator/ (const VectorField &a, double alpha)
    {return VectorField{a}/=alpha;}  

inline VectorField operator/ (const VectorField &a, const VectorField &b)
    {return VectorField{a}/=b;} 

inline VectorField sqr(const VectorField &a) 
    {return VectorField{a}.sqr();}

inline VectorField abs(const VectorField &a) 
    {return VectorField{a}.abs();}

VectorField circshift(const VectorField &v, int j);


template<typename VectorField>
void linComb(VectorField* v, 
    const std::vector<std::vector<int> > &extractions, 
    const std::vector<int> &target, const std::vector<double> & coeffs,
    const VectorField &flux)
{
    int ncoeffs = coeffs.size();
    int N = v->getLength();
    double data[v->getLength()];
    int unknowns = v->getUnknowns();
    for(int i = 0; i < N; i++)
    {
        data[i] = 0.0;
    }
    for(int i = 0; i < extractions[0].size(); i++) 
    {
        std::copy(flux.getField(i), flux.getField(i) + N, data);
        cblas_dscal(N, coeffs[ncoeffs - 1], data, 1);
        for(int j = 0; j < ncoeffs - 1; j++)
        {
            cblas_daxpy(N, coeffs[j], v->getField(extractions[j][i]), 1, data, 1);
        }
        v->setField(target[i], N, data);
    }
};

template<typename VectorField>
void linComb2(const VectorField &v1, const VectorField &v2, VectorField* v3, 
    double alpha, double beta)
{
    int N = v1.getLength();
    int unknowns = v1.getUnknowns();
    double * a1;
    double * a2;
    double * a3;
    for(int i = 0; i < unknowns; i++) 
    {
        a1 = v1.getField(i);
        a2 = v2.getField(i);
        a3 = v3->getField(i);
        for(int j = 0; j < N; j++)
        {
            a3[j] = alpha*a1[j] + beta*a2[j];
        }
    }
};


template<typename VectorField>
void linComb3(const VectorField &v1, const VectorField &v2, 
    const VectorField & v3, VectorField* v4, double alpha, double beta,
    double gamma)
{
    int N = v1.getLength();
    int unknowns = v1.getUnknowns();
    double * a1;
    double * a2;
    double * a3;
    double * a4;
    for(int i = 0; i < unknowns; i++) 
    {
        a1 = v1.getField(i);
        a2 = v2.getField(i);
        a3 = v3.getField(i);
        a4 = v4->getField(i);
        for(int j = 0; j < N; j++)
        {
            a4[j] = alpha*a1[j] + beta*a2[j] + gamma*a3[j];
        }
    }
};


template<typename VectorField>
void linComb4(const VectorField &v1, const VectorField &v2, 
    const VectorField & v3, const VectorField & v4, VectorField* v5,
    double alpha, double beta, double gamma, double delta)
{
    int N = v1.getLength();
    int unknowns = v1.getUnknowns();
    double * a1;
    double * a2;
    double * a3;
    double * a4;
    double * a5;
    for(int i = 0; i < unknowns; i++) 
    {
        a1 = v1.getField(i);
        a2 = v2.getField(i);
        a3 = v3.getField(i);
        a4 = v4.getField(i);
        a5 = v5->getField(i);
        for(int j = 0; j < N; j++)
        {
            a5[j] = alpha*a1[j] + beta*a2[j] + gamma*a3[j] + delta*a4[j];
        }
    }
};
      

#endif 