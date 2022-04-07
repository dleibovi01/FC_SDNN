/* Boundary conditions for 1D Linear advection equations*/

#ifndef BC_LA_1D_H
#define BC_LA_1D_H


#include "BC.h"
#include <vector>


class BC_LA_One_Wave : public BC
{

public:

    BC_LA_One_Wave() : BC{1} {};

    std::vector<double> enforceBC(const Node* node, const double t) const
    {
        std::vector<double> bdry_values;
        // double x = node->getPos();
        double x = node->getPosition()[0];
        if(x == 0)
        {
            bdry_values.push_back(0.0);
        }    
        return bdry_values;
    }

    std::vector<double> getBC_L(const double t) const
    {
        std::vector<double> bdry_values;
        bdry_values.push_back(0.0); 
        return bdry_values;     
    }

    std::vector<double> getBC_R(const double t) const
    {
        std::vector<double> bdry_values;
        bdry_values.push_back(0.0); 
        return bdry_values;    
    }

    std::vector<int> getEnforceableBdries(const Node* node, const double t)
    {
        std::vector<int> v;
        v.push_back(0);
        return v;

    }
};




#endif