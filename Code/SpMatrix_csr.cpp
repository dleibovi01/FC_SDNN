/* Sparse matrix CSR */

#include "SpMatrix_csr.h"


SpMatrix_csr::SpMatrix_csr(const sparse_index_base_t &_indexing, const int &_rows, 
        const int &_cols, const int* _rows_start, const int* _rows_end, 
        const int* _col_indx, const double* _values) : indexing{_indexing}, 
        rows{_rows}, cols{_cols}
{
    int status;
    int nnz = _rows_end[rows - 1] - 1;
    // std::cout << "nnz = " << nnz << std::endl;
    rows_start = new int[rows + 1];
    rows_end = new int[rows];
    col_indx = new int[nnz];
    values = new double[nnz];
    Matrix_handle = new sparse_matrix_t;

    std::copy(_rows_start, _rows_start + rows, rows_start);
    std::copy(_rows_end, _rows_end + rows, rows_end);
    rows_start[rows] = rows_end[rows-1];

    std::copy(_col_indx, _col_indx + nnz, col_indx);
    std::copy(_values, _values + nnz, values);

    status = mkl_sparse_d_create_csr(Matrix_handle, indexing, rows,
        cols, rows_start, rows_end, col_indx, values);
    descr.type = SPARSE_MATRIX_TYPE_GENERAL;

}

SpMatrix_csr::SpMatrix_csr(const SpMatrix_csr &spmat) : indexing{spmat.indexing}, 
        rows{spmat.rows}, cols{spmat.cols}
{
    int status;
    int nnz = spmat.rows_end[rows - 1] - 1;  
    rows_start = new int[rows + 1];
    rows_end = new int[rows];
    col_indx = new int[nnz];
    values = new double[nnz];
    Matrix_handle = new sparse_matrix_t;

    std::copy(spmat.rows_start, spmat.rows_start + rows, rows_start);
    std::copy(spmat.rows_end, spmat.rows_end + rows, rows_end);
    rows_start[rows] = rows_end[rows-1];

    std::copy(spmat.col_indx, spmat.col_indx + nnz, col_indx);
    std::copy(spmat.values, spmat.values + nnz, values);

    status = mkl_sparse_d_create_csr(Matrix_handle, indexing, rows,
        cols, rows_start, rows_end, col_indx, values);  
    descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    // descr.diag = SPARSE_DIAG_NON_UNIT;
}

SpMatrix_csr & SpMatrix_csr::operator=(const SpMatrix_csr &spmat)
{
    if(this == &spmat)
    {
        return *this;
    }
    else
    {
        delete [] rows_start;
        delete [] rows_end;
        delete [] col_indx;
        delete [] values;
        delete Matrix_handle;

        indexing = spmat.indexing;
        rows = spmat.rows;
        cols = spmat.cols;
        int status;
        // int nnz = spmat.rows_end[rows - 1];  
        int nnz = spmat.rows_end[rows - 1] - 1; 
        rows_start = new int[rows + 1];
        rows_end = new int[rows];
        col_indx = new int[nnz];
        values = new double[nnz];
        Matrix_handle = new sparse_matrix_t;

        std::copy(spmat.rows_start, spmat.rows_start + rows, rows_start);
        std::copy(spmat.rows_end, spmat.rows_end + rows, rows_end);
        rows_start[rows] = rows_end[rows-1];

        std::copy(spmat.col_indx, spmat.col_indx + nnz, col_indx);
        std::copy(spmat.values, spmat.values + nnz, values);

        status = mkl_sparse_d_create_csr(Matrix_handle, indexing, rows,
            cols, rows_start, rows_end, col_indx, values); 
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        // descr.diag = SPARSE_DIAG_NON_UNIT;
        return *this;         
    }
}

SpMatrix_csr::~SpMatrix_csr()
{
    // std::cout << "values " << std::endl;
    // Print_Mat(values, rows_end[rows - 1] - 1, 1);
    // std::cout << std::endl;
    delete [] rows_start;
    delete [] rows_end;
    delete [] col_indx;
    delete [] values;
    delete Matrix_handle;
}

void SpMatrix_csr::rowSum(double* sum) const
{
    std::vector<double> ones(cols, 1.0);
    MV(ones.data(), sum);
}

inline void SpMatrix_csr::MV(const double* x, double* y) const
{
    MV(1.0, x, 0.0, y);
}

void SpMatrix_csr::MV(double alpha, const double* x, double beta, double* y)
    const
{
    int status;
    status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, 
        *Matrix_handle, descr, x, beta, y);
}

void SpMatrix_csr::SMatVect(const double* x, double* y) const
{
    for(int i = 0; i < rows; i++)
    {
        // y[i] = 0.0;
        for(int j = rows_start[i] - 1; j < rows_start[i + 1] - 1; j++)
        {
            y[i] += values[j] * x[col_indx[j] - 1];
            // index++;
        }
    }
}

void SpMatrix_csr::dotMV(const double* x)
{
    int status;
    for(int i = 0; i < rows; i++)
    {
        for(int j = rows_start[i] - 1; j < rows_end[i] - 1; j++)
        {
            values[j] = values[j] * x[i];
        }
    }
    status = mkl_sparse_d_create_csr(Matrix_handle, indexing, rows,
        cols, rows_start, rows_end, col_indx, values); 
}

void SpMatrix_csr::print() const
{
    double entries[rows*cols];
    for(int i = 0; i < rows*cols; i++)
    {
        entries[i] = 0.0;
    }    
    for(int i = 0; i < rows; i++)
    {
        for(int j = rows_start[i] - 1; j < rows_end[i] - 1; j++)
        {
            entries[(col_indx[j] - 1)*rows + i] = values[j];
        }
    }
    Print_Mat(entries, rows, cols);
}

void SpMatrix_csr::setMVHint()
{
    int status;
    int calls = 1000000;
    status = mkl_sparse_set_mv_hint(*Matrix_handle,
        SPARSE_OPERATION_NON_TRANSPOSE, descr, calls);
}