
#ifndef PRINTING_H
#define PRINTING_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <string>
#include "VectorField1D.h"
#include "Node1D.h"
#include "Patch1D.h"
#include "Mesh.h"

void Print_Mat(const double * A, int nrows, int ncols);
void Print_Mat(const double * A, int nrows, int ncols, bool flag);
void Print_Mat(const std::complex<double> * A, int nrows, int ncols);
void Print_Mat(const int* A, int nrows, int ncols);

void Print_VectorField1D(const VectorField1D &field);
void Print_VectorField1D(const VectorField1D &field, int precision);
void Print_VectorField1D(const VectorField1D &field, bool print, int precision);
void Print_VectorField1D(const VectorField1D &field, bool print);
void Print_Node1D(const Node1D &node);
void Print_SparseMatrix_csr(int rows, int cols, int* rows_start, int* rows_end,
    int* col_indx, double* values);

template<class Patch>
void Print_Patch1D(const Patch &patch)
{
    int unknowns;
    std::cout << patch.getNnodes() << " nodes" << std::endl;
    for(int i = 0; i < patch.getNnodes(); i++)
    {
        std::cout << "Position = " << patch.getNode(i)->getPos() << " ";
        unknowns = patch.getNode(i)->getUnknowns();
        for(int j = 0; j < unknowns; j++)
        {
            std::cout << "Value " << j << " = " << (patch.getNode(i)->getValues())[j] << "  ";
        }
        std::cout << std::endl;
    }
}

template<class Patch>
void Print_Mesh1D(const Mesh<Patch> &mesh)
{
    std::vector<Patch*> patches = mesh.getPatches();
    for(int i = 0; i < patches.size(); i++)
    {
        std::cout << "Patch " << i + 1 << std::endl;
        Print_Patch1D(*patches[i]);
        std::cout << std::endl;
        std::cout << std::endl;

    }
}

template<class Patch>
void Print_Mesh1D(const Mesh<Patch> &mesh, int unknowns, int intrb,
    std::string filename)
{
    int N;
    std::vector<Patch*> patches = mesh.getPatches();
    std::ofstream myfile (filename);
    if (myfile.is_open())
    {
        myfile << patches.size();
        myfile << "\n";   
        myfile << unknowns;
        myfile << "\n";   
        myfile << intrb;
        myfile << "\n";            
        for(int i = 0; i < patches.size(); i++)
        {
            auto patch = patches[i];
            auto v = patches[i]->getFlow();
            N = v.getLength();
            myfile << N;
            myfile << "\n";             
            std::cout << "Patch " << i + 1 << std::endl;
            for(int j = 0; j < N; j++)
            {
                myfile << std::fixed << std::setprecision(17) << patch->getNode(j)->getPos();
                myfile << "\n";
            }
            for(int j = 0; j < unknowns; j++)
            {
                for(int k = 0; k < N; k++)
                {
                    // myfile << v.getFieldValue(j, k);
                    myfile << std::fixed << std::setprecision(17) << v.getFieldValue(j, k);
                    myfile << "\n";
                }
            }     
        }
        myfile.close();
    }
    else std::cout << "Unable to open file";    
}



template<class Patch>
void Print_Mesh1D(const Mesh<Patch> &mesh, int unknowns, int intrb, int vector,
    std::string filename)
{
    int N;
    unknowns = 1;
    std::vector<Patch*> patches = mesh.getPatches();
    std::ofstream myfile (filename);
    if (myfile.is_open())
    {
        myfile << patches.size();
        myfile << "\n";   
        myfile << unknowns;
        myfile << "\n";   
        myfile << intrb;
        myfile << "\n";            
        for(int i = 0; i < patches.size(); i++)
        {
            auto patch = patches[i];
            auto v = patches[i]->getFlow();
            N = v.getLength();
            myfile << N;
            myfile << "\n";             
            // std::cout << "Patch " << i + 1 << std::endl;
            for(int j = 0; j < N; j++)
            {
                myfile << std::fixed << std::setprecision(17) << patch->getNode(j)->getPos();
                myfile << "\n";
            }
            for(int j = 0; j < unknowns; j++)
            {
                for(int k = 0; k < N; k++)
                {
                    // myfile << v.getFieldValue(j, k);
                    myfile << std::fixed << std::setprecision(17) << v.getFieldValue(vector + j, k);
                    myfile << "\n";
                }
            }     
        }
        myfile.close();
    }
    else std::cout << "Unable to open file";    
}





#endif 