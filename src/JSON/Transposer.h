//-*-C++-*-

#ifndef PSIMAG_DCA_Transposer_H
#define PSIMAG_DCA_Transposer_H

#include <new>
#include <stdexcept>
#include <vector>
#include <cstddef>
#include <iostream>
#include <iomanip>

//#include "PSIMAGAssert.h"
#include "BLAS.h"
#include "LAPACK.h"
#include "OperationClosure.h"
#include "MatrixLike.h"
//#include "Tag.h"

namespace psimag {

  template<typename T> class Matrix;
  
  //======================================================================

  template<typename MatrixLikeType>
  class Transposer {
  public:
    
    typedef Transposer<MatrixLikeType>          ThisType;
    typedef typename MatrixLikeType::value_type value_type;
    typedef value_type                          FieldType;
    
    MatrixLikeType& mat;
    
    Transposer(MatrixLikeType& m):
      mat(m)
    {}
    
    MatrixLikeType& transpose() {
      return mat;
    }

    const MatrixLikeType& transpose() const {
      return mat;
    }

    //======================================================================

    ThisType& operator = (ThisType& other) {
      mat = other.mat;
    }

    //======================================================================

    void toJSN(std::ostream& os, std::string title, int width) const {
      MatrixLike::toJSN((*this),os,title,width);
    }

    //======================================================================

    // ----- MatrixLike
    
    size_t n_col() const { return mat.n_row(); }
    size_t n_row() const { return mat.n_col(); }
    
    value_type& operator () (size_t rowIndex, size_t colIndex) {
      // Assuming mat provides it's own asserts!
      return mat(colIndex,rowIndex);
    }
    

    const value_type& operator () (size_t rowIndex, size_t colIndex) const {
      // Assuming mat provides it's own asserts!
      return mat(colIndex,rowIndex);
    }

    void resize(size_t nRow, size_t nCol) {
      mat.resize(nCol,nRow);
    }

    // ----- VectorLike 

    // these probably should return elements is a different order ?? &*&*

    size_t size()  const { return mat.size(); }
    
    const value_type& operator[] (size_t componentIndex) const {
      return mat[componentIndex];
    }
    
    value_type& operator[] (size_t componentIndex)  {
      return mat[componentIndex];
    }
  };
  
  //====================================================================== This one owns it's Matrix

  template<typename MatrixLikeType>
  class Transposed {
  public:
    
    typedef typename MatrixLikeType::value_type value_type;
    typedef value_type                          FieldType;
    
    MatrixLikeType mat;
    
    Transposed(size_t nrows, size_t nCols):
      mat(nrows,nCols)
    {}
    
    Transposed(size_t nrows):
      mat(nrows)
    {}
    
    MatrixLikeType& transpose() {
      return mat;
    }

    const MatrixLikeType& transpose() const {
      return mat;
    }

    void eraseRowAndColumn(size_t r, size_t c){
      mat.eraseRowAndColumn(c,r);
    }

    // ----- MatrixLike
    
    size_t n_col() const { return mat.n_row(); }
    size_t n_row() const { return mat.n_col(); }
    
    value_type& operator () (size_t rowIndex, size_t colIndex) {
      // Assuming mat provides it's own asserts!
      return mat(colIndex,rowIndex);
    }
    

    value_type operator () (size_t rowIndex, size_t colIndex) const {
      // Assuming mat provides it's own asserts!
      return mat(colIndex,rowIndex);
    }

    void resize(size_t nRow, size_t nCol) {
      mat.resize(nCol,nRow);
    }

    // ----- VectorLike

    size_t size()  const { return mat.size(); }
    
    const value_type& operator[] (size_t componentIndex) const {
      return mat[componentIndex];
    }
    
    value_type& operator[] (size_t componentIndex)  {
      return mat[componentIndex];
    }
  };
  
} /* namespace psimag */


#endif 
