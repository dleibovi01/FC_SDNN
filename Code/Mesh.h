/* 1D Mesh */

#ifndef MESH_H
#define MESH_H

#include "Patch1D.h"
#include "Patch2D.h"
// #include "Patch1DUniform.h"

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

    // ~Mesh()
    // {
    //     while(!patches.empty())
    //     {
    //         delete patches.back();
    //         patches.pop_back();
    //     }
    // }

    
    void setPatch(Patch* patch, int i) {patches[i] = patch;}
    void setPatches(const std::vector<Patch*> p) {patches = p;}
    std::vector<Patch*> getPatches() const {return patches;}
    std::vector<Patch*>* getPatchesPtr() {return &patches;}
};

#endif
