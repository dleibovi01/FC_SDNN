
#include "printing.h"




void Print_Mat(const int * A, int nrows, int ncols)
{
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            std::cout << A[j*nrows + i] << " ";
        }
        std::cout << std::endl;
    }   
}


void Print_Mat(const double * A, int nrows, int ncols)
{
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            std::cout << A[j*nrows + i] << " ";
        }
        std::cout << std::endl;
    }   
}

void Print_Mat(const double * A, int nrows, int ncols, bool flag)
{
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            std::cout << A[j + i*ncols] << " ";
        }
        std::cout << std::endl;
    }   
}

void Print_Mat(const std::complex<double> * A, int nrows, int ncols)
{
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            std::cout << A[j*nrows + i].real() << " + " << A[j*nrows + i].imag() << " I ";
        }
        std::cout << std::endl;
    }   
}

void Print_VectorField(const VectorField &field)
{
    int N = field.getUnknowns();
    int length = field.getLength();
    std::vector<double*> flow = field.getFlow();
    double data [length*N] ; 
    for (int i = 0; i < N; i++)
    {
        for(int j = 0; j < length; j++)
        {
            data[i*length + j] = (flow[i])[j];
        }    
    }
    for(int i = 0; i < length; i++)
    {
        for(int j = 0; j < N; j++)
        {
            std::cout << data[j*length + i] << "  ";
        }
        std::cout << std::endl; 
    }
}

void Print_VectorField(const VectorField &field, int precision)
{
    int N = field.getUnknowns();
    int length = field.getLength();
    std::cout.precision(precision);
    std::vector<double*> flow = field.getFlow();
    double data [length*N] ; 
    for (int i = 0; i < N; i++)
    {
        for(int j = 0; j < length; j++)
        {
            data[i*length + j] = (flow[i])[j];
        }    
    }
    for(int i = 0; i < length; i++)
    {
        for(int j = 0; j < N; j++)
        {
            std::cout << data[j*length + i] << "  ";
        }
        std::cout << std::endl; 
    }
}

void Print_VectorField(const VectorField &field, bool print)
{
    if(print)
    {
        VectorField u{field.getUnknowns() + 1, field.getLength()};
        double temp[field.getLength()];
        for(int i = 0; i < field.getLength(); i++)
        {
            temp[i] = double(i + 1);
        }
        u.setField(0, field.getLength(), temp);
        for(int i = 0; i < field.getUnknowns(); i++)
        {
            u.setField(i+1, field.getLength(), field.getField(i));
        }
        Print_VectorField(u);
    }      
    else
    {
        Print_VectorField(field);
    } 
}


void Print_VectorField(const VectorField &field, bool print, int precision)
{
    if(print)
    {
        VectorField u{field.getUnknowns() + 1, field.getLength()};
        double temp[field.getLength()];
        for(int i = 0; i < field.getLength(); i++)
        {
            temp[i] = double(i + 1);
        }
        u.setField(0, field.getLength(), temp);
        for(int i = 0; i < field.getUnknowns(); i++)
        {
            u.setField(i+1, field.getLength(), field.getField(i));
        }
        Print_VectorField(u, precision);
    }      
    else
    {
        Print_VectorField(field, precision);
    } 
}

void Print_Node(const Node &node)
{
    std::cout << "Node position: " << node.getPosition()[0] << std::endl;
    std::cout << "Node index: " << node.getIndex() << std::endl;
    std::cout << "Unknowns: " << node.getUnknowns() << std::endl;
    std::cout << "Stored values:" << std::endl;
    //Print_Mat(node.getValues(), node.getUnknowns(), 1);
    double a[node.getUnknowns()];
    for(int i = 0; i < node.getUnknowns(); i++)
    {
        a[i] = node.getValues()[i];
    } 
    Print_Mat(a, node.getUnknowns(), 1);
}

void Print_SparseMatrix_csr(int rows, int cols, int* rows_start, int* rows_end,
    int* col_indx, double* values)
{
    double entries[rows*cols];
    int current_col = 0;
    int val_index = 0;
    for(int i = 0; i < rows*cols; i++)
    {
        entries[i] = 0.0;
    }
    for(int i = 0; i < rows; i++)
    {
        current_col = rows_start[i];
        for(int j = rows_start[i]; j < rows_end[i]; j++)
        {
            entries[rows*col_indx[j] + i] = values[j];
        }
    }
    Print_Mat(entries, rows, cols);

}

void Print_Patch2D(const Patch2D &patch)
{
    int unknowns = patch.getNode(0)->getUnknowns();
    // if(patch.getNnodes() > 0)
    //     unknowns = patch.getNode(0)->getUnknowns();
    int patchsize_x = patch.getNx();
    int patchsize_y = patch.getNy();
    for(int u = 0; u < unknowns; u++)
    {
        std::cout << "Unknown " << u << std::endl;
        for(int i = 0; i < patchsize_y; i++)
        {
            for(int j = 0; j < patchsize_x; j++)
            {
                std::cout << patch.getNode(patchsize_y*j + i)->getValues()[u] 
                    << " ";
            }
            std::cout << std::endl;
        }
    }
}


void Print_Mesh2DUniform(const Mesh2DUniform &mesh)
{
    int npatches_x = mesh.getNpatches_x();
    int npatches_y = mesh.getNpatches_y();
    int npatches = npatches_x * npatches_y;
    auto patches = mesh.getPatches();
    for(int j = 0; j < npatches_x; j++)
    {
        for(int i = 0; i < npatches_y; i++)
        {
            std::cout << "PATCH " << npatches_y*j + i << ". position : (" << 
                i << ", " << j << ")" << std::endl;
            Print_Patch2D(*patches[npatches_y*j + i]);
            std::cout << std::endl;
            std::cout << std::endl;
        }
    }   
}

