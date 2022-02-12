

#include "Node1D.h"


// Node1D::Node1D()
// {
//     pos = 0.0;
//     unknowns = 1;
//     index = -1;
//     values = new double[1];
// }

Node1D::Node1D(int _unknowns)
{
    pos = 0.0;
    unknowns = _unknowns;
    index = -1;
    //values = new double[unknowns];
    for(int i = 0; i < unknowns; i++)
    {
        values.push_back(0.0);
    }
}

Node1D::Node1D(double _pos, int _n, int _N, std::vector<double> _values)
{
    pos = _pos;
    index = _n;
    unknowns = _N;
    // values = new double[unknowns];
    // std::copy(_values, _values + unknowns, values);
    for(int i = 0; i < unknowns; i++)
    {
        values.push_back(_values[i]);
    }
}

// Node1D::~Node1D()
// {
//     delete[] values;
// }

void Node1D::setValues(const std::vector<double> input_values)
{
   //std::copy(input_values, input_values + unknowns, values);
   for(int i = 0; i < unknowns; i++)
   {
        values[i] = input_values[i];
   }
}

void Node1D::setUnknowns(int _unknowns)
{
    if(unknowns == 0)
    {
        unknowns = _unknowns;
        //values = new double[unknowns];
        for(int i = 0; i < unknowns; i++)
        {
            values.push_back(0.0);
        }
    }
}