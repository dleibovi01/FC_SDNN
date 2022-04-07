/* Node */

#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <vector>
#include <cmath>

class Node{

/* The abscissa of the node. */
// double pos;
std::vector<double> position;
/* The index of the node. */
int index;
/* The number of flow variables. */
int unknowns;
/* The flow data at the location of the node. */
std::vector<double> values;

public:

    // Node();

    Node(int _unknowns);
    
    Node(const std::vector<double> & _pos, int _n, int _N, 
        const std::vector<double> & _values) : position{_pos}, index{_n}, 
        unknowns{_N}, values{_values} {};

    std::vector<double> getPosition() const {return position;}

    int getIndex() const {return index;}

    int getUnknowns() const {return unknowns;}

    std::vector<double>  getValues() const {return values;}

    double getValue(int i) const {return values[i];}

    void setValues(const std::vector<double> input_values);

    void setValue(int u, double value) {values[u] = value;}

    void setIndex(int i) {index = i;}

    void setUnknowns(int _unknowns);

    void setPosition(const std::vector<double> & _pos) {position = _pos;}

    double getDist(const Node & node) const;


};




#endif 