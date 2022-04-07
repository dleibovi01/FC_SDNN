

#include "Node.h"


Node::Node(int _unknowns)
{
    unknowns = _unknowns;
    index = -1;
    //values = new double[unknowns];
    for(int i = 0; i < unknowns; i++)
    {
        values.push_back(0.0);
    }
}

void Node::setValues(const std::vector<double> input_values)
{
   for(int i = 0; i < unknowns; i++)
   {
        values[i] = input_values[i];
   }
}

void Node::setUnknowns(int _unknowns)
{
    if(unknowns == 0)
    {
        unknowns = _unknowns;
        for(int i = 0; i < unknowns; i++)
        {
            values.push_back(0.0);
        }
    }
}

double Node::getDist(const Node & node) const
{
    double distance = 0.0;
    std::vector<double> position2 = node.position;
    for(int i = 0; i < position.size(); i++)
    {
        distance += (position[i] - position2[i]) * (position[i] - position2[i]);
    }
    distance = std::sqrt(distance);
    return distance;
}