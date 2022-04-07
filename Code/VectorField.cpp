/* 1D vector field */


#include "VectorField.h"
#include "printing.h"
#include <typeinfo>


VectorField::VectorField()
{
    unknowns = 0;
    length = 0;
    flow.push_back(new double[0]);
}

VectorField::VectorField(const std::vector<double*> &field, int n)
{
    unknowns = field.size();
    length = n;
    for(int i = 0; i < unknowns; i++)
    {
        flow.push_back(new double[length]);
        setField(i, length, field[i]);
    }
}

VectorField::VectorField(const VectorField &field)
{
    unknowns = field.unknowns;
    length = field.length;
    for (int i = 0; i < unknowns; i++)
    {
        flow.push_back(new double[length]);
        setField(i, length, (field.flow)[i]);
    }
}

VectorField::VectorField(VectorField &&field)
{
    // std::cout <<"calling move construction" << std::endl;
    unknowns = field.unknowns;
    length = field.length;
    flow = field.flow;
    // while(!field.flow.empty())
    // {
    //     delete[] field.flow.back();
    //     field.flow.pop_back();
    // } 
    for(int i = 0; i < unknowns; i++)
    {
        field.flow[i] = nullptr;
    }
}


VectorField::VectorField(int _unknowns, int _length)
{
    unknowns = _unknowns;
    length = _length;
    for(int i = 0; i < unknowns; i++)
    {
        flow.push_back(new double[length]);
    }
}

VectorField & VectorField::operator= (const VectorField &v)
{
    if(this == &v)
    {
        return *this;
    }
    else
    {
        if((length == v.length) && (unknowns == v.unknowns))
        {
            for(int i = 0; i < unknowns; i++)
            {
                setField(i, length, v.getField(i));
            }
        }
        else
        {
            while(!flow.empty())
            {
                delete[] flow.back();
                unknowns--;
                flow.pop_back();
            } 

            unknowns = v.unknowns;
            length = v.length;
            for(int i = 0; i < unknowns; i++)
            {
                flow.push_back(new double[length]);
                setField(i, length, v.getField(i));
            }
        }

    }
    return *this;  
}

VectorField & VectorField::operator= (VectorField &&field)
{
    // std::cout <<"calling move assignment" << std::endl;
    if(this == &field)
    {
        return *this;
    }
    unknowns = field.unknowns;
    length = field.length;  
    while(!flow.empty())
    {
        delete[] flow.back();
        flow.pop_back();
    }   
    flow = field.flow;
    for(int i = 0; i < unknowns; i++)
    {
        field.flow[i] = nullptr;
    }    
    return *this;
}

void VectorField::setFlow(const std::vector<double*> flow_ext)
{
    for(int i = 0; i < unknowns; i++)
    {
        setField(i, length, flow_ext[i]);
    }
}


VectorField & VectorField::operator+= (const VectorField &v)
{
    double alpha = 1.0;
    int incx = 1;
    int incy = 1;
    double *y;
    for(int i = 0; i < unknowns; i++)
    {
        y = getField(i);
        cblas_daxpy(length, alpha, v.getField(i), incx, y, incy);
        setField(i, length, y);
    }
    return *this;
}

VectorField & VectorField::operator+= (const double a)
{
    double *y;
    for(int i = 0; i < unknowns; i++)
    {
        // y = flow[i];
        y = getField(i);
        for(int j = 0; j < length; j++)
        {
            y[j] += a;
        }
    }
    return *this;
}

VectorField & VectorField::operator-= (const VectorField &v)
{
    double alpha = - 1.0;
    int incx = 1;
    int incy = 1;
    double *y;
    for(int i = 0; i < unknowns; i++)
    {
        y = getField(i);
        cblas_daxpy(length, alpha, v.getField(i), incx, y, incy);
        setField(i, length, y);
    }
    return *this;
}

VectorField & VectorField::operator*= (const VectorField &v)
{
    double *y;
    for(int i = 0; i < unknowns; i++)
    {
        y = getField(i);
        vdMul(length, v.getField(i), y, y);
        setField(i, length, y);
    }
    return *this;
}



VectorField & VectorField::operator*= (double* v)
{
    double *y;
    for(int i = 0; i < unknowns; i++)
    {
        y = getField(i);
        vdMul(length, v, y, y);
        setField(i, length, y);
    }
    return *this;
}





VectorField & VectorField::operator*= (double alpha)
{
    double *y;
    int incx = 1;
    int incy = 1;
    for(int i = 0; i < unknowns; i++)
    {
        y = getField(i);
        cblas_dscal(length, alpha, y, incx);
        setField(i, length, y);
    }
    return *this;
}

VectorField & VectorField::sqr()
{
    double *y;
    for(int i = 0; i < unknowns; i++)
    {
        vdMul(length, getField(i), getField(i), getField(i));
    }
    return *this;
}


VectorField & VectorField::operator/= (const VectorField &v)
{
    // double *y;
    for(int i = 0; i < unknowns; i++)
    {
        vdDiv(length, getField(i), v.getField(i), getField(i));
    }
    return *this;
}


VectorField & VectorField::operator>>= (const VectorField &v)
{
    double *y;
    for(int i = 0; i < unknowns; i++)
    {
        y = getField(i);
        std::rotate(y, y + 1, y + length);
        setField(i, length, y);
    }
    return *this;
}

VectorField & VectorField::abs ()
{
    for(int i = 0; i < unknowns; i++)
    {
        vdAbs(length, getField(i), getField(i));
    }
    return *this;
}


void VectorField::circshift (int j)
{
    double *y;
    std::cout <<"circshift" << std::endl; 
    for(int i = 0; i < unknowns; i++)
    {
        y = getField(i);
        if(j > 0)
        {
            std::rotate(y, y + j, y + length);
        }
        else if ((j < 0) && (j > -length))
        {
            std::rotate(y, y + length + j, y + length);
        }
        setField(i, length, y);
    }
}


VectorField::~VectorField()
{
    while(!flow.empty())
    {
        delete[] flow.back();
        unknowns--;
        flow.pop_back();
    }
}

VectorField VectorField::extract(std::vector<int> indices) const
{
    VectorField w{indices.size(), length};
    for(int k = 0; k < indices.size(); k++)
    {
        w.setField(k, length, flow[indices[k]]);
    } 
    return w;
}

void VectorField::setFields(const VectorField &v,  std::vector<int> indices)
{
    for(int i = 0; i < indices.size(); i++)
    {
        setField(indices[i], length, v.getField(i));
    }
}

void VectorField::linComb(std::vector<double> coeffs, 
    std::vector<std::vector<int> > extractions, std::vector<int> target)
{
    // extractions[ncoeffs, unknowns]
    // target [unknowns]
    int ncoeffs = coeffs.size();
    double data[length];
    for(int i = 0; i < length; i++)
    {
        data[i] = 0.0;
    }
    for(int i = 0; i < extractions.size(); i++) // NOT SURE HERE
    {
        for(int j = 0; j < ncoeffs; j++)
        {
            cblas_daxpy(length, coeffs[j], flow[extractions[j][i]], 1, data, 1);
        }
        setField(target[i], length, data);
    }
}



VectorField circshift (const VectorField & v, int j)
{
    VectorField u{v};
    u.circshift(j);
    return u;
}


