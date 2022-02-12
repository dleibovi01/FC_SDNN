/* 1D node */

#ifndef NODE1D_H
#define NODE1D_H

#include <iostream>
#include <vector>
#include <cmath>

class Node1D{

/* The abscissa of the node. */
double pos;
/* The index of the node. */
int index;
/* The number of flow variables. */
int unknowns;
/* The flow data at the location of the node. */
std::vector<double> values;

public:

    // Node1D();

    Node1D(int _unknowns);
    
    Node1D(double pos, int n, int N, std::vector<double> values);

    //~Node1D();

    double getPos() const {return pos;}

    int getIndex() const {return index;}

    int getUnknowns() const {return unknowns;}

    std::vector<double>  getValues() const {return values;}

    double getValue(int i) const {return values[i];}

    void setValues(const std::vector<double> input_values);

    void setValue(int u, double value) {values[u] = value;}

    void setIndex(int i) {index = i;}

    void setUnknowns(int _unknowns);

    void setPos(double _pos) {pos = _pos;}

    double getDist(const Node1D & node) const {return std::abs(pos - node.pos);}



};




#endif 