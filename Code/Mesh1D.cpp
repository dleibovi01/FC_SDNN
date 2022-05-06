/* 1D Mesh */

#include "Mesh1D.h"
#include <numeric>



Mesh1DUniform::Mesh1DUniform(double a, double b, int n_patches, int patchsize, 
    int _overlap, int intrb, int unknowns, int l_b, int r_b)
{
    double L;
    double h;
    int N;
    int lb;
    int rb;    
    int intrbl;
    int intrbr;
    overlap = _overlap;
    N = n_patches*patchsize - (n_patches - 1)*overlap;
    L = b - a;
    h = L / (double (N - 1));
    std::vector<int> inner_nodes_onepatch(patchsize);
    std::vector<int> inner_nodes_leftpatch(patchsize - intrb);
    std::vector<int> inner_nodes_innerpatch(patchsize - 2*intrb);
    std::vector<int> inner_nodes_rightpatch(patchsize - intrb);
    std::iota(inner_nodes_onepatch.data(), inner_nodes_onepatch.data() +
        patchsize, 0);
    std::iota(inner_nodes_leftpatch.data(), inner_nodes_leftpatch.data() +
        patchsize - intrb, 0);
    std::iota(inner_nodes_innerpatch.data(), inner_nodes_innerpatch.data() +
        patchsize - 2*intrb, intrb);
    std::iota(inner_nodes_rightpatch.data(), inner_nodes_rightpatch.data() +
        patchsize - intrb, intrb);
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
        patches.push_back(new Patch1D(patchsize, unknowns, 
            a + double(i*(patchsize - overlap))*h, 
            a + double(i*(patchsize - overlap) + (patchsize - 1))*h , lb, rb, 
                intrbl, intrbr));
        if(n_patches == 1)
        {
            patches[i]->setInnerNodes(inner_nodes_onepatch);
        }
        else
        {
            if (i == 0)
            {
                patches[i]->setInnerNodes(inner_nodes_leftpatch);
            }
            else if(i == n_patches - 1)
            {
                patches[i]->setInnerNodes(inner_nodes_rightpatch);
            }
            else
            {
                patches[i]->setInnerNodes(inner_nodes_innerpatch);
            }
        }
    }
}


void Mesh1DUniform::setIntraPatchBC(int unknowns, int stage)
{
    int N;
    int bdry_elems_l;
    int bdry_elems_r;
    N = patches.size();
    if (N > 1)
    {
        for(int i = 0; i < patches.size(); i++)
        {
            if(i == 0)
            {
                patches[i]->setInnerBdry(patches[i + 1], overlap, unknowns, 
                    stage, false);
            }
            else if(i == patches.size() - 1)
            {
                patches[i]->setInnerBdry(patches[i - 1], overlap, unknowns, 
                    stage, true);
            }
            else
            {
                patches[i]->setInnerBdry(patches[i - 1], patches[i + 1], 
                    overlap, overlap, unknowns, stage);                   
            }              
        }
    }
}


