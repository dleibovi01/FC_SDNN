/* Boundary condition */

#ifndef BC_H
#define BC_H

#include "Patch1D.h"
#include "Patch1DUniform.h"
#include "Mesh.h"


class BC{

private: 
    class function{  

        std::string problem;
        public:
            function(std::string _problem){problem = _problem;}
            double operator()(double x, double t)
            {
                std::string problem0 = "";
                std::string problem1 = "LA_Gaussian";
                std::string problem2 = "step";
                std::string problem3 = "LA_test";
                std::string problem4 = "three_waves";
                std::string problem5 = "one_wave";
                std::string problem6 = "smooth_LA";
                double y = 0.0;
                if(problem.compare(problem0) == 0)
                {
                    if(x == 0)
                    {
                        y = 0.0;
                    }                    
                }
                else if(problem.compare(problem4) == 0)
                {
                    if(x == 0)
                    {
                        y = 0.0;
                    }                  
                }
                else if(problem.compare(problem5) == 0)
                {
                    if(x == 0)
                    {
                        y = 0.0;
                    }                  
                }    
                else if(problem.compare(problem6) == 0)
                {
                    y = 0.0;
                }                               
                return y;
            }
    };
    std::string problem;

public:
    BC(std::string _problem){problem = _problem;};

    template<typename Patch>
    void setBC(Patch* patch, double t)
    {
        function f{problem};
        std::vector<int> bdry_nodes = patch->getPhysBdryNodes(); 
        for(int i = 0; i < bdry_nodes.size(); i++)
        {
            // patch->getNode(bdry_nodes[i])->setValue(0, 
            //     f(patch->getNode(bdry_nodes[i])->getPos(), t));
            patch->setFlowValue(0, bdry_nodes[i], 
                f(patch->getNode(bdry_nodes[i])->getPos(), t));
            
        }
        // patch->NodesToVectorField();
    }

    template <typename Mesh>
    void operator() (Mesh* mesh, double t)
    {
        function f{problem};
        for(int i = 0; i < (mesh->getPatches()).size(); i++)
        {
            setBC(mesh->getPatches()[i], t);
        }
    }    
};

#endif