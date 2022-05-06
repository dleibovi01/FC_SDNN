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

// vector of indices of the nodes that are in the interior of the patch, i.e.
// not in the fringe. Needed for determining the SDNN viscosity assignments.
std::vector<int> inner_nodes;


public:

    Patch(){};
    Patch(const Patch &patch);
    ~Patch();
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
    const std::vector<int> & getInnerNodes() const {return inner_nodes;}
    void addInnerNode(int i) {inner_nodes.push_back(i);}
    void setInnerNodes(const std::vector<int> &in) {inner_nodes = in;}

    void NodesToVectorField();
    void VectorFieldToNodes();

};


#endif