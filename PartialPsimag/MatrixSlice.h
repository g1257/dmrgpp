//-*-C++-*-

/** \ingroup Data_Structures */
/*@{*/

/*! \file MatrixSlice.h  
 *
 */

#ifndef PSIMAG_MatrixSlice_H
#define PSIMAG_MatrixSlice_H

namespace psimag {
  
  //====================================================================== Row Slice
  
  template<typename MatrixLikeType> class TransposedRowSlice;
  
  template<typename MatrixLikeType>
  class RowSlice {
  public:
    
    typedef RowSlice<MatrixLikeType>            ThisType;
    typedef typename MatrixLikeType::value_type value_type;
    typedef typename MatrixLikeType::value_type FieldType;
    
    MatrixLikeType&       mat;
    const size_t          sliceRow;
    
    RowSlice(MatrixLikeType& m, size_t rowIndex):
      mat(m),
      sliceRow(rowIndex)
    {}
    
    // ----- MatrixLike
    
    inline size_t n_col() const { return mat.n_col(); }
    inline size_t n_row() const { return 1;           }
    
    inline
    value_type& operator () (size_t rowIndex, size_t colIndex) {
      ASSERT(rowIndex==0,
	     std::range_error("MatrixRowSlice(size_t rowIndex[bad], size_t colIndex)"));
      // Assuming mat provides it's own asserts!
      return mat(sliceRow,colIndex);
    }
    
    inline
    const value_type& operator () (size_t rowIndex, size_t colIndex) const {
      ASSERT(rowIndex==0,
	     std::range_error("MatrixRowSlice(size_t rowIndex[bad], size_t colIndex)"));
      // Assuming mat provides it's own asserts!
      return mat(sliceRow,colIndex);
    }

    TransposedRowSlice<MatrixLikeType> transpose() {
      return TransposedRowSlice<MatrixLikeType>(mat,sliceRow);
    }
    
    // ----- VectorLike

    template<typename VectorLike>
    ThisType& operator = (const VectorLike& v) {
      for(size_t i=0; i<v.size(); i++) 
	mat(sliceRow,i) = v[i];
      return *this;
    }

    inline size_t size()  const { return mat.n_col(); }
    
    inline
    const value_type& operator[] (size_t componentIndex) const {
      return mat(sliceRow,componentIndex);
    }
    
    inline
    value_type operator[] (size_t componentIndex)  {
      return mat(sliceRow,componentIndex);
    }
  };
  
  //====================================================================== Transposed Row Slice

  template<typename MatrixLikeType>
  class TransposedRowSlice:
    public RowSlice<MatrixLikeType>
  {
  public:
    
    typedef TransposedRowSlice<MatrixLikeType>  ThisType;
    typedef RowSlice<MatrixLikeType>            BaseType;
    typedef typename MatrixLikeType::value_type value_type;
    typedef typename MatrixLikeType::value_type FieldType;
    
    TransposedRowSlice(const MatrixLikeType& m, size_t rowIndex):
      BaseType(m,rowIndex)
    {}
    
    // ----- MatrixLike
    
    inline size_t n_col() const { return 1; }
    inline size_t n_row() const { return this->mat.n_col();           }
    
    RowSlice<MatrixLikeType> transpose() {
      return RowSlice<MatrixLikeType>(this->mat,this->sliceRow);
    }

    inline
    value_type& operator () (size_t rowIndex, size_t colIndex) {
      ASSERT(colIndex==0,
	     std::range_error("TransposedRowSlice(size_t rowInde, size_t colIndexx[bad])"));
      // Assuming mat provides it's own asserts!
      return this->mat(this->sliceRow,rowIndex);
    }
    
    inline
    const value_type& operator () (size_t rowIndex, size_t colIndex) const {
      ASSERT(colIndex==0,
	     std::range_error("TransposedRowSlice(size_t rowIndex, size_t colIndex[bad])"));
      // Assuming mat provides it's own asserts!
      return this->mat(this->sliceRow,rowIndex);
    }
    
    // ----- VectorLike (inherited)
  };
  

  //====================================================================== Col Slice
  
  template<typename MatrixLikeType> class TransposedColSlice;

  template<typename MatrixLikeType>
  class ColSlice {
  public:
    
    typedef ColSlice<MatrixLikeType>            ThisType;
    typedef typename MatrixLikeType::value_type value_type;
    typedef typename MatrixLikeType::value_type FieldType;
    
    MatrixLikeType&       mat;
    const size_t          sliceCol;
    
    ColSlice(MatrixLikeType& m, size_t colIndex):
      mat(m),
      sliceCol(colIndex)
    {}
    
    // ----- MatrixLike
    
    size_t n_col() const { return 1;           }
    size_t n_row() const { return mat.n_row(); }
    
    TransposedColSlice<MatrixLikeType> transpose() {
      return TransposedColSlice<MatrixLikeType>(mat,sliceCol);
    }

    value_type& operator () (size_t rowIndex, size_t colIndex) {
      ASSERT(colIndex==0,
	     std::range_error("ColSlice::(size_t rowIndex, size_t colIndex[bad])"));
      return mat(rowIndex,sliceCol);
    }
    
    inline
    const value_type& operator () (size_t rowIndex, size_t colIndex) const {
      ASSERT(colIndex==0,
	     std::range_error("ColSlice(size_t rowIndex, size_t colIndex[bad])"));
      // Assuming mat provides it's own asserts!
      return mat(rowIndex,sliceCol);
    }
    
    // ----- VectorLike
    
    size_t size()  const { return mat.n_row(); }
    
    inline
    const value_type& operator[] (size_t componentIndex) const {
      return mat(componentIndex,sliceCol);
    }
    
    inline
    value_type& operator[] (size_t componentIndex)  {
      return mat(componentIndex,sliceCol);
    }
  };
  
  //====================================================================== Transposed Col Slice
  
  template<typename MatrixLikeType>
  class TransposedColSlice {
  public:
    
    typedef TransposedRowSlice<MatrixLikeType>  ThisType;
    typedef ColSlice<MatrixLikeType>            BaseType;
    typedef typename MatrixLikeType::value_type value_type;
    typedef typename MatrixLikeType::value_type FieldType;

    TransposedColSlice(const MatrixLikeType& m, size_t colIndex):
      BaseType(m,colIndex)
    {}
    
    // ----- MatrixLike
    
    size_t n_col() const { return this-> mat.n_row();  }
    size_t n_row() const { return 1;                   }
    
    ColSlice<MatrixLikeType> transpose() {
      return ColSlice<MatrixLikeType>(this->mat,this->sliceCol);
    }

    value_type& operator () (size_t rowIndex, size_t colIndex) {
      ASSERT(rowIndex==0,
	     std::range_error("TransposedColSlice(size_t rowIndex[bad], size_t colIndex)"));
      return this->mat(colIndex,this->sliceCol);
    }
    
    inline
    const value_type& operator () (size_t rowIndex, size_t colIndex) const {
      ASSERT(rowIndex==0,
	     std::range_error("TransposedColSlice(size_t rowIndex[bad], size_t colIndex)"));
      // Assuming mat provides it's own asserts!
      return mat(colIndex,this->sliceCol);
    }
    
    // ----- VectorLike (inherited)
  };
  
  //======================================================================
  //======================================================================


  template<typename MatrixLikeType>
  std::ostream& operator << (std::ostream& os, const RowSlice<MatrixLikeType> ms) {
    os << "RowSlice[";
    for(size_t i=0; i< ms.size(); i++)
      os << ms[i] << " ";
    os << "]";
    return os;
  }
  template<typename MatrixLikeType>
  std::ostream& operator << (std::ostream& os, const TransposedRowSlice<MatrixLikeType> ms) {
    os << "TransposedRowSlice[";
    for(size_t i=0; i< ms.size(); i++)
      os << ms[i] << " ";
    os << "]";
    return os;
  }
  template<typename MatrixLikeType>
  std::ostream& operator << (std::ostream& os, const ColSlice<MatrixLikeType> ms) {
    os << "ColSlice[";
    for(size_t i=0; i< ms.size(); i++)
      os << ms[i] << " ";
    os << "]";
    return os;
  }
  template<typename MatrixLikeType>
  std::ostream& operator << (std::ostream& os, const TransposedColSlice<MatrixLikeType> ms) {
    os << "TransposedColSlice[";
    for(size_t i=0; i< ms.size(); i++)
      os << ms[i] << " ";
    os << "]";
    return os;
  }
  


} // end namespace PSIMAG


/*@}*/
#endif
