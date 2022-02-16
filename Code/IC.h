/* Initial condition */

#ifndef IC_H
#define IC_H

#include <string>
#include <math.h>
#include <iostream>
#include <vector>
#include "Patch1D.h"
#include "Patch1DUniform.h"
#include "Mesh.h"



/* An initial condition functor */
class IC{

private: 
    class function{  

        std::string problem;
        public:
            function(std::string _problem){problem = _problem;}
            std::vector<double> operator()(double x)
            {
                std::string problem0 = "";
                std::string problem1 = "LA_Gaussian";
                std::string problem2 = "step";
                std::string problem3 = "LA_test";
                std::string problem4 = "three_waves";
                std::string problem5 = "one_wave";
                std::string problem6 = "smooth_LA";
                std::string problem7 = "Euler1D_Sod";
                std::vector<double> y;
                if(problem.compare(problem0) == 0)
                {
                    y.push_back(0.0);
                }
                else if(problem.compare(problem1) == 0)
                {
                    double sigma = 0.1;
                    double pi = 3.1415926535897932384686;
                    y.push_back(1.0 / (sqrt(2.0 * pi) * sigma)*
                        exp(- 0.5 * (x - 0.5) * (x - 0.5) / sigma / sigma));
                }
                else if(problem.compare(problem2) == 0)
                {
                    if(x < 0.5)
                    {
                        y.push_back(1.0);
                    }
                    else
                    {
                        y.push_back(0.0);
                    }
                }
                else if(problem.compare(problem3) == 0)
                {
                    y.push_back(exp(- 160 * (x - 0.5) * (x - 0.5)));
                }   
                else if(problem.compare(problem4) == 0)
                {
                    if(0.2 < x && x <= 0.3)
                    {
                        y.push_back(10.0*(x - 0.2));
                    }
                    else if(0.3 < x && x <= 0.4)
                    {
                        y.push_back(10.0*(0.4 - x));
                    }
                    else if(0.6 < x && x <= 0.8)
                    {
                        y.push_back(1.0);
                    }
                    else if(1 < x && x <= 1.2)
                    {
                        y.push_back(100.0*(x - 1.0)*(1.2 - x));
                    }
                    else
                    {
                        y.push_back(0.);
                    }    
                }   
                else if(problem.compare(problem5) == 0)
                {
                    if(0.2 < x && x <= 0.3)
                    {
                        y.push_back(1.0);
                    }
                    else
                    {
                        y.push_back(0.);
                    }    
                }      
                else if(problem.compare(problem6) == 0)
                {
                    y.push_back(std::exp(- 160 * (x - 0.5) * (x - 0.5)));
                }     
                else if(problem.compare(problem7) == 0)
                {
                    if(0 <= x && x <= 0.5)
                    {
                        y.push_back(1.0);
                        y.push_back(0.0);
                        y.push_back(2.5);
                    }
                    else
                    {
                        y.push_back(0.125);
                        y.push_back(0.0);
                        y.push_back(0.25);
                    }    
                }                              
                return y;
            }
    };
    std::string problem;
    int unknowns;

public:
    IC(std::string _problem, int _unknowns) : problem{_problem}, 
        unknowns{_unknowns} {};
    template<typename Patch>
    void setIC(Patch* patch)
    {
        function f{problem};
        std::vector<double> init_values;
        for(int i = 0; i < patch->getNnodes(); i++)
        {
            init_values = f(patch->getNode(i)->getPos());
            for(int j = 0; j < unknowns; j++)
            {
                patch->getNode(i)->setValue(j, init_values[j]);
            }
        }
        patch->NodesToVectorField();
    }

    template <typename Mesh>
    void operator() (Mesh* mesh)
    {
        function f{problem};
        for(int i = 0; i < (mesh->getPatches()).size(); i++)
        {
            setIC(mesh->getPatches()[i]);
        }
    }    
};




#endif