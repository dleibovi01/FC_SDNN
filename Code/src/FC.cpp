/* FC related functions */

#include "FC.h"


void read_FC_Data(double *A, double *Q, int d, int C, std::string filename_A, std::string filename_Q)
{
    int i = 0;
    double data;
    std::ifstream Adata (filename_A);
    if(Adata.is_open())
    {
        while(Adata >> data)
        {
            A[i] = data;
            i = i + 1;       
        }
    }

    i = 0;
    std::ifstream Qdata (filename_Q);
    if(Qdata.is_open())
    {
        while(Qdata >> data)
        {
            Q[i] = data;
            i = i + 1;       
        }
    }


}