/* 1D Mesh */

#include "Mesh.h"



// Mesh1DUniform::Mesh1DUniform(double a, double b, int n_patches, int patchsize, 
//     int _overlap, int intrb, int unknowns, bool l_b, bool r_b)
Mesh1DUniform::Mesh1DUniform(double a, double b, int n_patches, int patchsize, 
    int _overlap, int intrb, int unknowns, int l_b, int r_b)
{
    double L;
    double h;
    int N;
    // bool lb;
    // bool rb;
    int lb;
    int rb;    
    int intrbl;
    int intrbr;
    overlap = _overlap;
    N = n_patches*patchsize - (n_patches - 1)*overlap;
    L = b - a;
    h = L/(double (N - 1));
    for(int i = 0; i < n_patches; i++)
    {
        if (i == 0)
        {
            lb = l_b;
            intrbl = 0;
        }
        else
        {
            lb = 0;
            intrbl = intrb;
        }
        
        if(i == n_patches - 1)
        {          
            rb = r_b;
            intrbr = 0;
        }
        else
        {
            rb = 0;
            intrbr = intrb; 
        }
        patches.push_back(new Patch1DUniform(patchsize, unknowns, 
            a + double(i*(patchsize - overlap))*h, 
            a + double(i*(patchsize - overlap) + (patchsize - 1))*h , lb, rb, 
            intrbl, intrbr));
    }
}
