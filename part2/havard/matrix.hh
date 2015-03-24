// An extremely simple matrix template
//
// Written by Arne Morten Kvarving. NO rights reserved.
//
// This code is in the public domain.

#pragma once

#include <iostream>
#include <vector>
#include <mpi.h>
#include <omp.h>
//#include <unistd.h> // sleep
using namespace std;

//! \brief Enum for data layout in memory
enum DataOrder {
  RowMajor, //!< Row major - C ordering
  ColMajor  //!< Column major - Fortran ordering
};

//! \brief Matrix template
//! \param scalar Type of matrix elements
//! \param order Data layout in memory
template<typename scalar, DataOrder order=RowMajor>
class Matrix {
 protected:
     size_t m_rows; 		//!< Number of rows in matrix
     size_t m_cols; 		//!< Number of columns in matrix
     std::vector<scalar> m_data; //!< The actual data
 public:
     //! \brief Constructor
     Matrix(): m_rows(0), m_cols(0), m_data(0) {}

     //! \brief Constructor
     //! \param rows Number of rows in matrix
     //! \param cols Number of columns in matrix
     Matrix(size_t rows, size_t cols) : m_rows(rows), m_cols(cols)
    {
        m_data.resize(m_rows*m_cols);
    }

     //! \brief Returns a reference to a matrix element
     //! \brief i Row number
     //! \brief j Column number
     //! \details This hides the data layout in memory from user.
     //!          Indices are zero based.
     inline scalar& operator()(size_t i, size_t j)
     {
         if (order == RowMajor)
             return m_data[i*m_cols+j];
         else
             return m_data[j*m_rows+i];
     }

     //! \brief Print matrix to a stream
     //! \param precision Precision for the scalars
     //! \param out The stream to print to. Defaults to std::cout
     void print(std::streamsize precision=0, std::ostream& out=std::cout)
     {
         std::streamsize oldprec;
         if (precision > 0) {
             oldprec = out.precision();
             out.precision(precision);
         }

         for (size_t i=0;i<m_rows;++i) {
             for (size_t j=0;j<m_cols;++j)
                 out << (*this)(i,j) << " ";
             out << std::endl;
         }
         if (precision > 0)
             out.precision(oldprec);
     }

     scalar getRows() {
         return m_rows;
     }
     scalar getCols() {
         return m_cols;
     }
};


//! \brief Matrix template for distributed with MPI
//! \param scalar Type of matrix elements
//! \param order Data layout in memory
template<typename scalar, DataOrder order=RowMajor>
class MatrixMPI: public Matrix<scalar, order>{
 private:
     // Variables
     int m_commSize;                    //!< Size of communicator
     int m_commRank;                    //!< Rank in communicator
     MPI_Comm &m_comm;                  //!< Communicator
     std::vector<size_t> m_cols_all;    //!< Vector containing all local columns
     //-----------------------------------------------------------
     // Functions
     //! \brief Copy blocks from temp to this, and order by col instead of row
     //! \param temp Temporary MatrixMPI receive buffer
     //! \param senddisps Send displacements array (for MPI_Alltoallv)
     void copyBlocksByRowToCol(MatrixMPI<double, RowMajor> &temp, vector<int> senddisps){
         double *tempptr = temp.colFront(0);
         size_t tempit = 0;
         size_t rowStart, rowNextStart;
         // iterate over block k
         for (size_t k = 0; k < m_commSize; ++k) {
             rowStart      = senddisps[k];
             rowNextStart  = rowStart + m_cols_all[k];
             for (size_t c = 0; c < this->m_cols; ++c) {
                 for (size_t r = rowStart; r < rowNextStart; ++r) {
                     (*this)(r,c) = tempptr[tempit++];
                 }
             }
         }
     }

 public:
     //! \brief Constructor
     //! \param glob_dim Global dim of quadratic matrix
     //! \param comm MPI communicator
     MatrixMPI(size_t glob_dim, MPI_Comm &comm): Matrix<scalar, order>(), m_comm(comm) {
     //MatrixMPI(size_t glob_dim, MPI_Comm &comm) {
         //m_comm = comm;
         MPI_Comm_size(comm, &m_commSize);
         MPI_Comm_rank(comm, &m_commRank);
         size_t cols = glob_dim/m_commSize;  //local cols
         size_t rows = glob_dim;
         m_cols_all.resize(m_commSize, cols);

         // Add an remaining columns to rank 0, 1, ...
         if ((glob_dim%m_commSize) > (size_t) m_commRank)
             cols++;
         for (size_t i = 0; i < glob_dim%m_commSize; ++i) {
             m_cols_all[i]++;
         }
         this->m_rows = rows;
         this->m_cols = cols;
         this->m_data.resize(this->m_rows * this->m_cols);
     }

     std::vector<size_t> getCols_all() {
         return m_cols_all;
     }

     //! \brief returns pointer to first element in local column
     //! \param j Local column number
     scalar* colFront(size_t col) {
         if (order == RowMajor) {
             if (col == 0)
                 return &this->m_data[0];
             std::cout << "No implementation for RowMajor other than col = 0" << std::endl;
             exit(1);
         } 
         else
             return &this->m_data[col*this->m_rows];
     }

     //! \brief Transpose matrix, but double memory usage
     void transpose() {
         // Create data type consisting of one local row
         MPI_Datatype locrow;
         create_types(&locrow, this->m_cols, 1, this->m_rows, sizeof(double));

         // Create types to MPI_Alltoallv
         vector<int> sendcnts(m_commSize);
         for (size_t i = 0; i < m_commSize; ++i) {
             // Need to convert size_t to int
             sendcnts[i] = m_cols_all[i];
         }
         vector<int> senddisps(m_commSize, 0);
         for (size_t i = 1; i < m_commSize; ++i) {
             senddisps[i] = senddisps[i-1] + sendcnts[i-1];
         }

         // Receive data of type double 
         vector<int> receivecnts(m_commSize);
         for (size_t i = 0; i < m_commSize; ++i) {
             receivecnts[i] = sendcnts[i]*this->m_cols;
         }
         vector<int> receivedisps(m_commSize, 0);
         for (size_t i = 1; i < m_commSize; ++i) {
             receivedisps[i] = receivedisps[i-1] + receivecnts[i-1];
         }

         MatrixMPI<double, RowMajor> temp(this->m_rows, m_comm);
         //cout << "Ready to send" << endl;
         MPI_Alltoallv(&this->m_data[0], &sendcnts.front(), 
                 &senddisps.front(), locrow, temp.colFront(0), 
                 &receivecnts.front(), &receivedisps.front(), MPI_DOUBLE, m_comm);

         // Reorder blocks to get transpose
         this->copyBlocksByRowToCol(temp, senddisps);
         //cout << "Done reordering" << endl;

         MPI_Type_free(&locrow);
     }

};

// Create and commit MPI datatypes
void create_types(MPI_Datatype* newtype, int block_count, int block_length, int stride, MPI_Aint extent){
    MPI_Datatype tmp;
    //Create new datatype
    MPI_Type_vector(block_count, block_length, stride, MPI_DOUBLE, &tmp);
    // Set extent (used for displacement)
    MPI_Type_create_resized(tmp, 0, extent, newtype);
    MPI_Type_free(&tmp);
    //Commit the above
    MPI_Type_commit(newtype);
}


double WallTime () {
    return omp_get_wtime();
}


