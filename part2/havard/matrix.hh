// An extremely simple matrix template
//
// Written by Arne Morten Kvarving. NO rights reserved.
//
// This code is in the public domain.

#pragma once

#include <iostream>
#include <vector>
#include <mpi.h>

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

     scalar getrows() {
         return m_rows;
     }

     //// finds next index for "transpose"
     //inline void nextIxTrans(size_t i, size_t j, size_t &ii, size_t &jj) {
         //ii = m_cols*j + i/m_cols;
         //jj = i % m_cols;
     //}


     //void mkTransposed() {
         //if (order == RowMajor)
             //exit(1);

         //size_t ii, jj, kk, r0, c0, r1, c1, r2, c2;
         //ii = 0;
         //scalar tmp1, tmp2;
         //for (size_t i = m_cols; i > 1; --i) {
             //jj = ii-1;
             //ii++;
             //for (size_t j = 0; j < i; ++j) {
                 //jj++;
                 //kk = ii;
                 //for (size_t k = 0; k < (i-1); ++k) {
                     //// index = (kk + jj*m_cols, ii-1)
                     //r0 = kk + jj*m_cols;
                     //c0 = ii - 1;
                     //nextIxTrans(r0, c0, r1, c1);
                     //nextIxTrans(r1, c1, r2, c2);
                     //tmp1 = (*this)(r1, c1);
                     //tmp2 = (*this)(r2, c2);
                     //(*this)(r1, c1) = (*this)(r0, c0);
                     //(*this)(r2, c2) = tmp1;
                     //(*this)(r0, c0) = tmp2;
                     //kk++;
                     ////cout << r0 << ", " << c0 << endl;
                 //}
             //}
         //}
     //}
};


//! \brief Matrix template for distributed with MPI
//! \param scalar Type of matrix elements
//! \param order Data layout in memory
template<typename scalar, DataOrder order=RowMajor>
class MatrixMPI: public Matrix<scalar, order>{
 private:
     int m_commSize;
     int m_commRank;
     MPI_Comm &m_comm;
 public:
     //! \brief Constructor
     //! \param rows Number of rows in matrix
     //! \param cols Number of columns in matrix
     MatrixMPI(size_t rows, size_t cols, MPI_Comm &comm): Matrix<scalar, order>(rows, cols), m_comm(comm) {
     //MatrixMPI(size_t rows, size_t cols): Matrix<scalar, order>(rows, cols)  {
         //MPI_Comm_size(MPI_COMM_WORLD, &m_commSize);
         //MPI_Comm_rank(MPI_COMM_WORLD, &m_commRank);
         MPI_Comm_size(m_comm, &m_commSize);
         MPI_Comm_rank(m_comm, &m_commRank);
     }

     //! \brief returns pointer to first element in local column
     //! \param j Local column number
     scalar* colFront(size_t j) {
         if (order == RowMajor) {
             std::cout << "No implementation for RowMajor" << std::endl;
             exit(1);
         } 
         else
             return &this->m_data[j*this->m_rows];
     }

     //! \brief transpose sub matrix m_cols x m_cols
     //! \param rowStart Upper left element in submatrix
     void transposeSub(size_t rowStart) {
         size_t r, r2, c2;
         for (size_t c = 0; c < this->m_cols; ++c) {
             r = rowStart + c + 1;
             r2 = c + rowStart;
             for (size_t j = c+1; j < this->m_cols; ++j) {
                 c2 = r-rowStart;
                 swap(r, c, r2, c2); 
                 r++;
             }
         }
     }

     //! \brief swap elemts in matrix
     //! \param r1 Row 1 index
     //! \param c1 Col 1 index
     //! \param r2 Row 2 index
     //! \param c2 Col 2 index
     inline void swap(size_t r1, size_t c1, size_t r2, size_t c2) {
         scalar tmp = (*this)(r1, c1);
         (*this)(r1, c1) = (*this)(r2, c2);
         (*this)(r2, c2) = tmp;
     }

     //! \brief Does transpose using MPI
     void transpose() {
                     
         // New datatype
         //MPI_Datatype MPI_scalar;
         //MPI_Type_create_resized(MPI_INT, 0, sizeof(scalar), &MPI_scalar);
         //MPI_Type_commit(&MPI_scalar);

         //Sending to get transpose
         for (size_t i = 0; i < this->m_cols; ++i) {
             //MPI_Alltoall(this->colFront(i), this->m_cols, MPI_DOUBLE, 
                     //this->colFront(i), this->m_cols, MPI_DOUBLE, 
                     //MPI_COMM_WORLD);
             //MPI_Alltoall(this->colFront(i), this->m_cols, MPI_scalar, 
                     //this->colFront(i), this->m_cols, MPI_scalar, 
                     //m_comm);
             MPI_Alltoall(this->colFront(i), this->m_cols, MPI_DOUBLE, 
                     this->colFront(i), this->m_cols, MPI_DOUBLE, 
                     m_comm);
         }
         // Restructuring to get transposed
         size_t rowStart = 0;
         for (size_t i = 0; i < m_commSize; ++i) { 
             this->transposeSub(rowStart);
             rowStart += this->m_cols;
         }
     }


};

// Create and commit MPI datatypes
void create_types(MPI_Datatype* newtype, int block_count, int block_length, int stride){
    //For coloumns that neighbour other processes
    MPI_Type_vector(block_count, block_length, stride, MPI_UNSIGNED_CHAR, newtype);
    //Commit the above
    MPI_Type_commit(newtype);
}






