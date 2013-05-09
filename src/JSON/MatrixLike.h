//-*-C++-*-
// Author: Michael S. Summers (ORNL)
//
#ifndef PSIMAG_MatrixLike_H
#define PSIMAG_MatrixLike_H

#include <new>
#include <stdexcept>
#include <vector>
#include <cstddef>
#include <iostream>
#include <iomanip>
#include "String.h"

#include "PSIMAGAssert.h"
#include "BLAS.h"
#include "LAPACK.h"
//#include "Tag.h"

namespace psimag {

  
  /*! \brief Template functions for types that provide member functions:
    
      - n_row, 
      - n_col, and 
      - operator()
      
      and support:     
      
      typedef typename MatrixLikeType::value_type FieldType;
      
      \author Mike Summers
  */
  namespace MatrixLike {
    
    template<typename VectorLikeType>
    class DoubleVectorWrap {

    public:

      typedef typename PsimagLite::Vector<VectorLikeType>::Type        DoubleVectorType;
      typedef typename DoubleVectorType::value_type InnerVectorType;
      typedef typename InnerVectorType::value_type  value_type;
      typedef typename InnerVectorType::value_type  FieldType;

      const DoubleVectorType& doubleVector;
      
      DoubleVectorWrap(DoubleVectorType& dv):
	doubleVector(dv)
      {}

      size_t n_row() const { return doubleVector.size(); }
      size_t n_col() const { return doubleVector[0].size(); }
      
      value_type operator() (size_t i, size_t j) {
	return doubleVector[i][j];
      }

      const value_type operator() (size_t i, size_t j) const {
	return doubleVector[i][j];
      }
      
      size_t size() const {return n_row()*n_col();}
      
    };

    template<typename MatrixLikeType>
    void size(const MatrixLikeType& matrix)  {
      return matrix.n_row() * matrix.n_col();
    }
    
    template<typename MatrixLikeType>
    void printRow(const MatrixLikeType& matrix,
		  size_t row,
		  std::ostream& os,
		  int width=13) {
      os << "[";
      for(size_t j=0; j<matrix.n_col(); j++) 
	os << " " << std::setw(width) << matrix(row,j) << ",";
      os << "]";
    }
    
    template<typename MatrixLikeType>
    void print(const MatrixLikeType& matrix,
	       std::ostream& os,
	       int width=13) {
      for(size_t i=0; i<matrix.n_row(); i++) {
	for(size_t j=0; j<matrix.n_col; j++) 
	  os << " " << std::setw(width) << matrix(i,j);
	os << "\n";
      }
    }
    
    template<typename MatrixLikeType>
    void printList(const MatrixLikeType& matrix,
		   std::ostream& os,
		   int           width  = 13,
		   PsimagLite::String   offset = "") {
      os << "[";
      
      size_t lastRow = matrix.n_row()-1;
      for(size_t row=0; row<matrix.n_row(); row++) {
	if(row!= 0)
	  os << offset;
	printRow(matrix,row,os,width);
	if(row != lastRow)
	  os << ",\n";
	else
	  os << "";
      }
      os << "]";
    }
    
    template<typename MatrixLikeType>
    void printArray(const MatrixLikeType& matrix,
		    std::ostream& os,
		    int width=13) {
      os << "array(";
      printList(matrix,os,width,"       ");
      os << ")";
    }
    
    /** \ingroup ostream
     *
     * Output stream operator for SuperCrystal
     */
    /*template<typename MatrixLikeType>
    Tag toXML(const MatrixLikeType& matrix,
	      String name="MatrixLike") {
      Tag tag(name);
      tag["rows"] = matrix.n_row();
      tag["cols"] = matrix.n_col();
      tag.content << std::endl;
      tag.content << std::endl;
      print(matrix,tag.content);
      return tag;
    }*/
    
    // ======================================================================  
    
    template<typename MatrixLikeType>
    void toJSN(const MatrixLikeType& matrix,
	       std::ostream& os, 
	       PsimagLite::String title="MatrixLike",
	       int width=13) {
      
      os.precision(width);
      os << std::fixed;
      os << "{ \'tile\': \'"    << title        << "\', \n"
	 << " \'rows\':"    << matrix.n_row() << ", \n"
	 << " \'cols\':"    << matrix.n_col() << ", \n"
	 << " \'data\': ";
      printArray(matrix,os,width);
      os << "\n}\n";
    }
    
    // ======================================================================  
    
    template<typename MatrixLikeType>
    void writeVTK(const MatrixLikeType& matrix,
		  PsimagLite::String fileName,
		  PsimagLite::String version    = "0.1",
		  PsimagLite::String byte_order = "BigEndian") {
      
      throw PsimagLite::LogicError("================================================= Later we will write a VTK file for fileName << \n");    
    }
    

    //======================================================================

    template<typename VectorLikeType,
	     typename MatrixLikeType> 
    VectorLikeType& setVectorFromRow(const MatrixLikeType& matrix,
				     size_t row, 
				     VectorLikeType& result)  {
      for (size_t col=0; col < matrix.n_col(); col++) {
	result[col] = matrix(row,col);
      }
      return result;
    }
    
    template<typename VectorLikeType,
	     typename MatrixLikeType> 
    VectorLikeType& setVectorFromCol(const MatrixLikeType& matrix,
				 size_t col, 
				 VectorLikeType& result)  {
      for (size_t row=0; row < matrix.n_row(); row++) {
	result[row] = matrix(row,col);
      }
      return result;
    }
    
    template<typename VectorLikeType1, 
	     typename VectorLikeType2>
    typename VectorLikeType2::value_type scalarProduct(const VectorLikeType1& v1,
						       const VectorLikeType2& v2) {
      typename VectorLikeType2::value_type result = v1[0] * v2[0];
      for(size_t i=1; i < v2.size(); i++)
	result += v1[i] * v2[i];
      return result;
    }

    template<typename MatrixLikeType> 
    inline
    MatrixLikeType& increment(MatrixLikeType& matrix,
			      typename MatrixLikeType::value_type val) {
      for (size_t row=0; row < matrix.n_row(); row++) 
	for (size_t col=0; col < matrix.n_col(); col++) 
	  matrix(row,col) += val;
      return matrix;
    }
    
    template<typename MatrixLikeType> 
    inline
    MatrixLikeType& decrement(MatrixLikeType& matrix,
			      typename MatrixLikeType::value_type val) {
      for (size_t row=0; row < matrix.n_row(); row++) 
	for (size_t col=0; col < matrix.n_col(); col++) 
	  matrix(row,col) -= val;
      return matrix;
    }
    
    template<typename MatrixLikeType> 
    inline
    MatrixLikeType& times(MatrixLikeType& matrix,
			  typename MatrixLikeType::value_type val) {
      for (size_t row=0; row < matrix.n_row(); row++) 
	for (size_t col=0; col < matrix.n_col(); col++) 
	  matrix(row,col) *= val;
      return matrix;
    }
    
    template<typename MatrixLikeType> 
    inline
    MatrixLikeType& divide(MatrixLikeType& matrix,
			   typename MatrixLikeType::value_type val) {
      for (size_t row=0; row < matrix.n_row(); row++) 
	for (size_t col=0; col < matrix.n_col(); col++) 
	  matrix(row,col) /= val;
      return matrix;
    }
    
    template<typename MatrixLikeType> 
    inline
    MatrixLikeType& squareElements(MatrixLikeType& matrix) {
      for (size_t row=0; row < matrix.n_row(); row++) 
	for (size_t col=0; col < matrix.n_col(); col++) 
	  matrix(row,col) *= matrix(row,col);
      return matrix;
    }
    
    template<typename MatrixLikeType1,
	     typename MatrixLikeType2> 
    inline
    MatrixLikeType2& copy(const MatrixLikeType1& matrix1,
			  MatrixLikeType2& matrix2) {
      for (size_t row=0; row < matrix2.n_row(); row++) 
	for (size_t col=0; col < matrix2.n_col(); col++) 
	  matrix2(row,col) = matrix1(row,col);
      return matrix2;
    }

    template<typename MatrixLikeType1,
	     typename MatrixLikeType2> 
    inline
    bool equals(const MatrixLikeType1& matrix1,
		const MatrixLikeType2& matrix2) {
      for (size_t row=0; row < matrix1.n_row(); row++) 
	for (size_t col=0; col < matrix1.n_col(); col++) 
	  if (matrix2(row,col) != matrix1(row,col)) 
	    return false;
      return true;
    }
    
    template<typename MatrixLikeType>
    class INVERT {
    public:

      typedef typename MatrixLikeType::value_type T;

      int operator () (const MatrixLikeType& A, PsimagLite::Matrix<T>& INV) {
	
	ASSERT(A.n_row()==A.n_col(),
	       std::range_error("A.n_row() != A.n_col() in INVERT for Matrix<T>"));
	ASSERT(A.n_row()==INV.n_row(),
	       std::range_error("A.n_row() != INV.n_row() in INVERT for Matrix<T>"));
	ASSERT(A.n_col()==INV.n_col(),
	       std::range_error("A.n_col() != INV.n_col() in INVERT for Matrix<T>"));
	
	copy(A,INV);
	int N  (static_cast<int>(A.n_row()));
	int LDA(N);  // Give this more though later &*&*

	pivots.resize(N);

	int info;

	// GET LU Decomposition which will be stored in INV
	LAPACK::GETRF(N,N,&INV(0,0),LDA,&pivots[0],info);
      
	if (info != 0) {
	  PsimagLite::OstringStream message;
	  message << "INVERT(A,B): GETRF failed! \n";
	  if (info < 0) 
	    message << "The " << -info << "th argument had an illegal value!\n";
	  else
	    message << "U("<<info<<","<<info<<") was exactly zero. \n";
	  throw PsimagLite::LogicError(message.str());
	}

	int    lwork(-1);
	double blocksize;

	// Get the optimal block size
	LAPACK::GETRI(N,&INV(0,0),LDA,&pivots[0],blocksize,-1,info);
      
	lwork  = static_cast<int>(blocksize)*N;
	PsimagLite::Vector<double>::Type work(lwork);
      
	LAPACK::GETRI(N,&INV(0,0),LDA,&pivots[0],&work[0],lwork,info);
 
	if (info != 0) {
	  PsimagLite::OstringStream message;
	  message << "INVERT(A,B): GETRI failed! \n";
	  if (info < 0) 
	    message << "The " << -info << "th argument had an illegal value!\n";
	  else
	    message << "U("<<info<<","<<info<<") was exactly zero. \n";
	  throw PsimagLite::LogicError(message.str());
	}
	
	return info;
      }
    private:
      PsimagLite::Vector<int>::Type pivots;
    };

  } /* namespace MatrixLike */
  
} /* namespace psimag */


#endif 
