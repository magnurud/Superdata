// An extremely simple matrix template
//
// Written by Arne Morten Kvarving. NO rights reserved.
//
// This code is in the public domain.

#pragma once

#include <iostream>
#include <vector>
#include <iostream>
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
 private:
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

     // swap elemts
     inline void swap(size_t r1, size_t c1, size_t r2, size_t c2) {
         scalar tmp = (*this)(r1, c1);
         (*this)(r1, c1) = (*this)(r2, c2);
         (*this)(r2, c2) = tmp;
     }

     // finds next index for "transpose"
     inline void nextIxTrans(size_t i, size_t j, size_t &ii, size_t &jj) {
         ii = m_cols*j + i/m_cols;
         jj = i % m_cols;
     }

     void mkTransposed() {
         if (order == RowMajor)
             exit(1);

         size_t ii, jj, kk, r0, c0, r1, c1, r2, c2;
         ii = 0;
         scalar tmp1, tmp2;
         for (size_t i = m_cols; i > 1; --i) {
             jj = ii-1;
             ii++;
             for (size_t j = 0; j < i; ++j) {
                 jj++;
                 kk = ii;
                 for (size_t k = 0; k < (i-1); ++k) {
                     // index = (kk + jj*m_cols, ii-1)
                     r0 = kk + jj*m_cols;
                     c0 = ii - 1;
                     nextIxTrans(r0, c0, r1, c1);
                     nextIxTrans(r1, c1, r2, c2);
                     tmp1 = (*this)(r1, c1);
                     tmp2 = (*this)(r2, c2);
                     (*this)(r1, c1) = (*this)(r0, c0);
                     (*this)(r2, c2) = tmp1;
                     (*this)(r0, c0) = tmp2;
                     kk++;
                     cout << r0 << ", " << c0 << endl;
                 }
             }
         }
     }
};
