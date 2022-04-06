// Sparse Matrix class file

#ifndef SPMATRIX_CSR_H
#define SPMATRIX_CSR_H

#include <vector>
#include "mkl_spblas.h"
#include "printing.h"

class SpMatrix_csr {

sparse_index_base_t indexing;
int rows;
int cols;
int* rows_start;
int* rows_end;
int* col_indx;
double* values;
sparse_matrix_t* Matrix_handle;
struct matrix_descr descr;

public:
    
    SpMatrix_csr(const sparse_index_base_t &_indexing, const int &_rows, 
        const int &_cols, const int* _rows_start, const int* _rows_end, 
        const int* _col_indx, const double* _values);

    SpMatrix_csr(const SpMatrix_csr &spmat);

    SpMatrix_csr & operator=(const SpMatrix_csr &spmat);

    ~SpMatrix_csr();

    int getRows() const {return rows;}
    int getCols() const {return cols;}
    int* getRowsStart() const {return rows_start;}
    int* getRowsEnd() const {return rows_end;}
    int* getColIndx() const {return col_indx;}
    double* getValues() const {return values;}
    sparse_matrix_t* getHandle() const {return Matrix_handle;}
    matrix_descr getDescr() const {return descr;}
    void optimize() {int status = mkl_sparse_optimize(*Matrix_handle);}
    void setMVHint();

    void rowSum(double* sum) const;
    void MV(const double* x, double* y) const;

    // CSR sparse matrix dense vector multiplication
    void SMatVect(const double* x, double* y) const;
    void MV(double alpha, const double* x, double beta, double* y) const;
    void dotMV(const double* x);
    // void operator*(double* x, double* y) const;

    void print() const;

    

    
};




#endif 