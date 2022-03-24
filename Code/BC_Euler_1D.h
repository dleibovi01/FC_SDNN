/* Boundary conditions for 1D Linear advection equations*/

#ifndef BC_Euler_1D_H
#define BC_Euler_1D_H


#include "BC.h"
#include <vector>


class BC_Euler_Sod : public BC
{
public:

    BC_Euler_Sod() : BC{3} {};

    std::vector<double> enforceBC(const Node1D* node, const double t)
    {
        std::vector<double> bdry_values;
        double x = node->getPos();
        if(x < 0.5)
        {
            bdry_values.push_back(1.0);
            bdry_values.push_back(0.0);
            bdry_values.push_back(2.5);
        }  
        else if(x > 0.5)
        {
            bdry_values.push_back(0.125);
            bdry_values.push_back(0.0);
            bdry_values.push_back(0.25);
        }  
        return bdry_values;
    }

    std::vector<double> getBC_L(const double t)
    {
        std::vector<double> bdry_values;
        bdry_values.push_back(1.0);
        bdry_values.push_back(0.0);
        bdry_values.push_back(2.5);  
        return bdry_values;     
    }

    std::vector<double> getBC_R(const double t)
    {
        std::vector<double> bdry_values;
        bdry_values.push_back(0.125);
        bdry_values.push_back(0.0);
        bdry_values.push_back(0.25);   
        return bdry_values;    
    }
};




#endif