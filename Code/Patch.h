/* Patch */

#ifndef PATCH_H
#define PATCH_H

#include "Node.h"
#include "VectorField.h"
#include <iostream>

class Patch{

protected:

VectorField v;
std::vector<Node*> nodes;
int Nnodes;
std::vector<int> phys_bdry_nodes;

public:

    const VectorField & getFlow() const {return v;}
    VectorField & getFlowRef() {return v;}
    VectorField* getFlowPtr() {return &v;}
    std::vector<Node*> getNodes() const {return nodes;}
    Node* getNode(int i) const {return nodes[i];}
    int getNnodes() const {return Nnodes;}
    int getUnknowns() const {return v.getUnknowns();}
    const std::vector<int> & getPhysBdryNodes() const {return phys_bdry_nodes;}
    void setField(const VectorField &_flow){v = _flow;}
    void setFlowValue(int i, int j, double d){v.setFieldValue(i, j, d);}

    void NodesToVectorField();
    void VectorFieldToNodes();

};


#endif