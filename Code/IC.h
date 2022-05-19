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
            std::vector<double> operator()(std::vector<double> x)
            {
                std::string problem0 = "";
                std::string problem1 = "LA_Gaussian";
                std::string problem2 = "step";
                std::string problem3 = "LA_test";
                std::string problem4 = "three_waves";
                std::string problem5 = "one_wave";
                std::string problem6 = "smooth_LA";
                std::string problem7 = "Euler1D_Sod";
                std::string problem8 = "2D_test";
                std::string problem9 = "Burgers_Square";
                std::string problem10 = "Euler2D_Riemann4";
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
                        exp(- 0.5 * (x[0] - 0.5) * (x[0] - 0.5) / sigma / sigma));
                }
                else if(problem.compare(problem2) == 0)
                {
                    if(x[0] < 0.5)
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
                    y.push_back(exp(- 160 * (x[0] - 0.5) * (x[0] - 0.5)));
                }   
                else if(problem.compare(problem4) == 0)
                {
                    if(0.2 < x[0] && x[0] <= 0.3)
                    {
                        y.push_back(10.0*(x[0] - 0.2));
                    }
                    else if(0.3 < x[0] && x[0] <= 0.4)
                    {
                        y.push_back(10.0*(0.4 - x[0]));
                    }
                    else if(0.6 < x[0] && x[0] <= 0.8)
                    {
                        y.push_back(1.0);
                    }
                    else if(1 < x[0] && x[0] <= 1.2)
                    {
                        y.push_back(100.0*(x[0] - 1.0)*(1.2 - x[0]));
                    }
                    else
                    {
                        y.push_back(0.);
                    }    
                }   
                else if(problem.compare(problem5) == 0)
                {
                    if(0.2 < x[0] && x[0] <= 0.3)
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
                    y.push_back(std::exp(- 160 * (x[0] - 0.5) * (x[0] - 0.5)));
                }     
                else if(problem.compare(problem7) == 0)
                {
                    if(0 <= x[0] && x[0] <= 0.5)
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
                else if(problem.compare(problem8) == 0)
                {
                    y.push_back(std::exp(std::cos(x[0]) * std::sin(x[1])));
                }      
                else if(problem.compare(problem9) == 0)
                {
                    if((x[0] >= 0.5) && (x[1] >= 0.5))
                        y.push_back(-1.0);
                    if((x[0] < 0.5) && (x[1] >= 0.5))
                        y.push_back(-0.2);
                    if((x[0] < 0.5) && (x[1] < 0.5))
                        y.push_back(0.5);
                    if((x[0] >= 0.5) && (x[1] < 0.5))
                        y.push_back(0.8);                                             
                }     
                else if(problem.compare(problem10) == 0)
                {
                    double gamma = 7.0/5.0;

                    double rho_1 = 1.1;
                    double vx_1 = 0;
                    double vy_1 = 0;
                    double p_1 = 1.1;
                    double E_1 = p_1/(gamma - 1.0)
                        + 0.5*rho_1*(vx_1*vx_1 + vy_1*vy_1);

                    double rho_2 = 0.5065;
                    double vx_2 = 0.8939;
                    double vy_2 = 0;
                    double p_2 = 0.35;
                    double E_2 = p_2/(gamma - 1.0)
                        + 0.5*rho_2*(vx_2*vx_2 + vy_2*vy_2);

                    double rho_3 = 1.1;
                    double vx_3 = 0.8939;
                    double vy_3 = 0.8939;
                    double p_3 = 1.1;
                    double E_3 = p_3/(gamma - 1.0)
                         + 0.5*rho_3*(vx_3*vx_3 + vy_3*vy_3);

                    double rho_4 = 0.5065;
                    double vx_4 = 0;
                    double vy_4 = 0.8939;
                    double p_4 = 0.35;
                    double E_4 = p_4/(gamma - 1)
                         + 0.5*rho_4*(vx_4*vx_4 + vy_4*vy_4);


                    if((x[0] > 0.6) && (x[1] > 0.6))
                    {
                        y.push_back(rho_1);
                        y.push_back(rho_1*vx_1);
                        y.push_back(rho_1*vy_1);
                        y.push_back(E_1);
                    }
                    if((x[0] <= 0.6) && (x[1] > 0.6))
                    {
                        y.push_back(rho_2);
                        y.push_back(rho_2*vx_2);
                        y.push_back(rho_2*vy_2);
                        y.push_back(E_2);
                    }
                    if((x[0] <= 0.6) && (x[1] <= 0.6))
                    {
                        y.push_back(rho_3);
                        y.push_back(rho_3*vx_3);
                        y.push_back(rho_3*vy_3);
                        y.push_back(E_3);
                    }
                    if((x[0] > 0.6) && (x[1] <= 0.6))
                    {
                        y.push_back(rho_4);
                        y.push_back(rho_4*vx_4);
                        y.push_back(rho_4*vy_4);
                        y.push_back(E_4); 
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
            init_values = f(patch->getNode(i)->getPosition());
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