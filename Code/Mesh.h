/* 1D Mesh */

#ifndef MESH_H
#define MESH_H

#include "Patch1D.h"
#include "Patch1DUniform.h"

template <typename Patch>
class Mesh{

protected:

std::vector<Patch*> patches;

public:

    Mesh() {}
    Mesh(const std::vector<Patch*> p)
    {
        int n = p.size();
        for(int i = 0; i < n; i++)
        {
            patches.push_back(p[i]);
        }
    }
    /* Copy constructor */
    Mesh(const Mesh &mesh) 
    {
        std::vector<Patch*> p = mesh.getPatches();
        int n = p.size();
        for(int i = 0; i < n; i++)
        {
            patches.push_back(p[i]);
        }        
    }
    /* Copy assignment */
    Mesh & operator=(const Mesh &mesh)
    {
        pacthes.clear();
        std::vector<Patch*> p = mesh.getPatches();
        int n = p.size();
        for(int i = 0; i < n; i++)
        {
            patches.push_back(p[i]);
        }   
        return *this;       
    }
    void setPatch(Patch* patch, int i) {patches[i] = patch;}
    void setPatches(const std::vector<Patch*> p) {patches = p;}
    std::vector<Patch*> getPatches() const {return patches;}
    std::vector<Patch*>* getPatchesPtr() {return &patches;}
};


class Mesh1DUniform : public Mesh<Patch1DUniform>
{
public:
    int overlap;
    /* Constructor */
   
    Mesh1DUniform(double a, double b, int n_patches, int patchsize, 
        int _overlap, int intrb, int unknowns, int l_b, int r_b);  

    /* Copy constructor */

    int getOverlap() const {return overlap;}

    double getH() const {return (patches[0]->getH());}

    void setIntraPatchBC(int unknowns, int stage)
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
};
#endif
