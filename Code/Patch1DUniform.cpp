/* Uniform 1D patch. */

#include "Patch1DUniform.h"


Patch1DUniform::Patch1DUniform(int N, int _unknowns, double _a, double _b, 
bool lb, bool rb, int intrbl, int intrbr)
{
    a = _a;
    b = _b;
    Nnodes = N;
    double h = (b - a)/double(N - 1);
    for(int i = 0; i < N; i++)
    {
        nodes.push_back(new Node1D(_unknowns));
        nodes[i]->setIndex(i);
        nodes[i]->setPos(a + i*h);
    }
    if(lb)
    {
        phys_bdry_nodes.push_back(0);
    }
    if(rb)
    {
        phys_bdry_nodes.push_back(N - 1);
    }
    intra_patch_nodes_l = intrbl;
    intra_patch_nodes_r = intrbr;
    //v.setLength(N);
    //v.setUnknows(_unknowns);
    v = VectorField1D{_unknowns, N};
}
