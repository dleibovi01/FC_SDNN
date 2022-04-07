/* Boundary condition */

#ifndef BC_H
#define BC_H

#include "Patch1D.h"
#include "Patch1DUniform.h"
#include "Node.h"
#include "Mesh.h"
#include <vector>
// #include "BC_LA_1D.h"


// class BC{

// private: 
//     class function{  

//         std::string problem;
//         public:
//             function(std::string _problem){problem = _problem;}
//             std::vector<double> operator()(double x, double t)
//             {
//                 std::string problem0 = "";
//                 std::string problem1 = "LA_Gaussian";
//                 std::string problem2 = "step";
//                 std::string problem3 = "LA_test";
//                 std::string problem4 = "three_waves";
//                 std::string problem5 = "one_wave";
//                 std::string problem6 = "smooth_LA";
//                 std::string problem7 = "Euler1D_Sod";
//                 std::vector<double> y;
//                 if(problem.compare(problem0) == 0)
//                 {
//                     if(x == 0)
//                     {
//                         y.push_back(0.0);
//                     }                    
//                 }
//                if(problem.compare(problem3) == 0)
//                 {
//                     if(x == 0)
//                     {
//                         y.push_back(0.0);
//                     }                    
//                 }
//                 else if(problem.compare(problem4) == 0)
//                 {
//                     if(x == 0)
//                     {
//                         y.push_back(0.0);
//                     }                  
//                 }
//                 else if(problem.compare(problem5) == 0)
//                 {
//                     if(x == 0)
//                     {
//                         y.push_back(0.0);
//                     }                  
//                 }    
//                 else if(problem.compare(problem6) == 0)
//                 {
//                     y.push_back(0.0);
//                 }           
//                 else if(problem.compare(problem7) == 0)
//                 {
//                    if(x < 0.5)
//                     {
//                         y.push_back(1.0);
//                         y.push_back(0.0);
//                         y.push_back(2.5);
//                     }  
//                    else if(x > 0.5)
//                     {
//                         y.push_back(0.125);
//                         y.push_back(0.0);
//                         y.push_back(0.25);
//                     }  
//                 }                         
//                 return y;
//             }
//     };
//     std::string problem;
//     int unknowns;

// public:
//     BC(std::string _problem, int _unknowns) : problem{_problem}, 
//         unknowns{_unknowns} {};

//     template<typename Patch>
//     void setBC(Patch* patch, double t)
//     {
//         setBC(patch, t, 0);
//     }

//     template<typename Patch>
//     void setBC(Patch* patch, double t, int stage)
//     {
//         function f{problem};
//         std::vector<int> bdry_nodes = patch->getPhysBdryNodes(); 
//         std::vector<double> bdry_values;
//         for(int i = 0; i < bdry_nodes.size(); i++)
//         {
//             bdry_values = f(patch->getNode(bdry_nodes[i])->getPos(), t);
//             for(int j = 0; j < unknowns; j++)
//             {
//                 patch->setFlowValue(stage*unknowns + j, bdry_nodes[i],
//                     bdry_values[j]);
//             }
//         }
//     }

//     template <typename Mesh>
//     void operator() (Mesh* mesh, double t)
//     {
//         function f{problem};
//         for(int i = 0; i < (mesh->getPatches()).size(); i++)
//         {
//             setBC(mesh->getPatches()[i], t);
//         }
//     }    
// };

// #endif

class BC{

int unknowns;

public:

    BC(int _unknowns) : unknowns{_unknowns} {};

    template<typename Patch>
    void setBC(Patch* patch, double t)
    {
        setBC(patch, t, 0);
    }

    template<typename Patch>
    void setBC(Patch* patch, double t, int stage)
    {
        std::vector<int> bdry_nodes = patch->getPhysBdryNodes(); 
        std::vector<double> bdry_values;
        std::vector<int> enforceable_bdries;
        for(int i = 0; i < bdry_nodes.size(); i++)
        {
            bdry_values = enforceBC(patch->getNode(bdry_nodes[i]), t);
            enforceable_bdries = 
                getEnforceableBdries(patch->getNode(bdry_nodes[i]), t);
            for(int j = 0; j < enforceable_bdries.size(); j++)
            {
                patch->setFlowValue(stage*unknowns + enforceable_bdries[j],
                    bdry_nodes[i], bdry_values[j]);
            }
        }
    }

    template <typename Mesh>
    void operator() (Mesh* mesh, double t)
    {
        for(int i = 0; i < (mesh->getPatches()).size(); i++)
        {
            setBC(mesh->getPatches()[i], t);
        }
    }    

protected:

    virtual std::vector<double> enforceBC(const Node* node, const double t)
    {
        std::vector<double> v;
        return v;
    }

    virtual std::vector<int> getEnforceableBdries(const Node* node,
        const double t) const
    {
        std::vector<int> v;
        return v;        
    }

};





#endif