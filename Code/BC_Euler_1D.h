/* Boundary conditions for 1D Euler equations*/

#ifndef BC_Euler_1D_H
#define BC_Euler_1D_H


#include "BC.h"
#include <vector>


class BC_Euler_Sod : public BC
{

    std::vector<double> bdry_values_l;
    std::vector<double> bdry_values_r;
    std::vector<int> enforceable_bdries;

public:

    BC_Euler_Sod() : BC{3}
    {
        bdry_values_l.push_back(1.0);
        bdry_values_l.push_back(0.0);
        bdry_values_l.push_back(2.5); 
        bdry_values_r.push_back(0.125);
        bdry_values_r.push_back(0.0);
        bdry_values_r.push_back(0.25);   
        enforceable_bdries.push_back(0);   
        enforceable_bdries.push_back(1);  
        enforceable_bdries.push_back(2);   
    }


    std::vector<double> enforceBC(const Node* node, const double t) 
    {
        std::vector<double> bdry_values;
        // double x = node->getPos();
        double x = node->getPosition()[0];
        if(x == 0.0)
        {
            return bdry_values_l;
        }  
        else if(x == 1.0)
        {
            return bdry_values_r;
        }  
        return bdry_values;
    }

    std::vector<double> getBC_L(const double t) const {return bdry_values_l;}

    std::vector<double> getBC_R(const double t) const {return bdry_values_r;}

    std::vector<int> getEnforceableBdries(const Node* node,
        const double t) const {return enforceable_bdries;}
};

class BC_Euler_Sod_NN : public BC
{

    std::vector<double> bdry_values_l;
    std::vector<double> bdry_values_r;
    std::vector<int> enforceable_bdries;
public:

    BC_Euler_Sod_NN() : BC{3}
    {
        bdry_values_l.push_back(0.0);
        bdry_values_l.push_back(0.0);
        bdry_values_l.push_back(0.0); 
        bdry_values_r.push_back(0.0);
        bdry_values_r.push_back(0.0);
        bdry_values_r.push_back(0.0);       
    }


    std::vector<double> enforceBC(const Node* node, const double t) 
    {
        std::vector<double> bdry_values;
        // double x = node->getPos();
        double x = node->getPosition()[0];
        if(x == 0.0)
        {
            return bdry_values_l;
        }  
        else if(x == 1.0)
        {
            return bdry_values_r;
        }  
        return bdry_values;
    }

    std::vector<double> getBC_L(const double t) const {return bdry_values_l;}

    std::vector<double> getBC_R(const double t) const {return bdry_values_r;}

    std::vector<int> getEnforceableBdries(const Node* node, const double t)
        const {return enforceable_bdries;}

};


class BC_Euler_Sod_ND : public BC
{

    std::vector<double> bdry_values_l;
    std::vector<double> bdry_values_r;
    std::vector<int> enforceable_bdries_l;
    std::vector<int> enforceable_bdries_r;
public:

    BC_Euler_Sod_ND() : BC{3}
    {
        bdry_values_l.push_back(0.0);
        bdry_values_l.push_back(0.0);
        bdry_values_l.push_back(0.0); 
        bdry_values_r.push_back(0.125);
        bdry_values_r.push_back(0.0);
        bdry_values_r.push_back(0.25); 
        enforceable_bdries_r.push_back(0);   
        enforceable_bdries_r.push_back(1);  
        enforceable_bdries_r.push_back(2);        
    }


    std::vector<double> enforceBC(const Node* node, const double t) 
    {
        std::vector<double> bdry_values;
        // double x = node->getPos();
        double x = node->getPosition()[0];
        if(x == 0.0)
        {
            return bdry_values_l;
        }  
        else if(x == 1.0)
        {
            return bdry_values_r;
        }  
        return bdry_values;
    }

    std::vector<double> getBC_L(const double t) const {return bdry_values_l;}

    std::vector<double> getBC_R(const double t) const {return bdry_values_r;}

    std::vector<int> getEnforceableBdries(const Node* node, const double t)
        const 
    {
        std::vector<int> enforceable_bdries;
        // double x = node->getPos();
        double x = node->getPosition()[0];
        if(x == 0.0)
        {
            return enforceable_bdries_l;
        }  
        else if(x == 1.0)
        {
            return enforceable_bdries_r;
        }
        return enforceable_bdries;
    }

};

#endif