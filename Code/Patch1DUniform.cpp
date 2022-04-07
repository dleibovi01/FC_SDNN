/* Uniform 1D patch. */

#include "Patch1DUniform.h"

Patch1DUniform::Patch1DUniform(int N, int _unknowns, double _a, double _b, 
    int lb, int rb, int intrbl, int intrbr)
{
    a = _a;
    b = _b;
    Nnodes = N;
    double h = (b - a)/double(N - 1);
    std::vector<double> position;
    position.push_back(0.0);
    for(int i = 0; i < N; i++)
    {
        nodes.push_back(new Node(_unknowns));
        nodes[i]->setIndex(i);
        position[0] = a + i*h;
        nodes[i]->setPosition(position);
    }
    for(int i = 0; i < lb; i++)
    {
        phys_bdry_nodes.push_back(i);
    }
    for(int i = 0; i < rb; i++)
    {
        phys_bdry_nodes.push_back(N - rb + i);
    }
    intra_patch_nodes_l = intrbl;
    intra_patch_nodes_r = intrbr;
    v = VectorField{_unknowns, N};
}
