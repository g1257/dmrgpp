//-*-C++-*-

/** \ingroup Data_Structures */
/*@{*/

/*! \file IndexedMatrix.h  
 *
 */

#ifndef PSIMAG_IndexedMatrix_H
#define PSIMAG_IndexedMatrix_H

namespace psimag {
  
  //====================================================================== 

  template<typename MatrixLikeType,
	   typename IndexType>
  class IndexedMatrix {
  public:
    
    typedef IndexedMatrix<MatrixLikeType,IndexType>       ThisType;
    typedef typename MatrixLikeType::value_type value_type;
    typedef typename MatrixLikeType::value_type FieldType;
    
    const MatrixLikeType&       mat;
    IndexType                   rowIndices;
    IndexType                   colIndices;
    
    IndexedMatrix(const MatrixLikeType& m, 
		  IndexType   rIndices,
		  IndexType   cIndices):
      mat(m),
      rowIndices(rIndices),
      colIndices(cIndices)
    {}
    
    // ----- MatrixLike
    
    inline size_t n_col() const { return colIndices.size(); }
    inline size_t n_row() const { return rowIndices.size(); }
    
    inline
    value_type& operator () (size_t rowIndex, size_t colIndex) {
      ASSERT(rowIndex < rowIndices.size(),
	     std::range_error("IndexedMatrix(size_t rowIndex[bad], size_t colIndex)"));
      ASSERT(colIndex < colIndices.size(),
	     std::range_error("IndexedMatrix(size_t rowIndex[bad], size_t colIndex)"));
      // Assuming mat provides it's own asserts!
      return mat(rowIndices[rowIndex],colIndices[colIndex]);
    }
    
    inline
    const value_type& operator () (size_t rowIndex, size_t colIndex) const {
      ASSERT(rowIndex < rowIndices.size(),
	     std::range_error("IndexedMatrix(size_t rowIndex[bad], size_t colIndex)"));
      ASSERT(colIndex < colIndices.size(),
	     std::range_error("IndexedMatrix(size_t rowIndex[bad], size_t colIndex)"));
      // Assuming mat provides it's own asserts!
      return mat(rowIndices[rowIndex],colIndices[colIndex]);
    }
    
  };
  
} // end namespace psimag


/*@}*/
#endif
