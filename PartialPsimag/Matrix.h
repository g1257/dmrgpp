//-*-C++-*-

#ifndef PSIMAG_Matrix_H
#define PSIMAG_Matrix_H
#include <cstring>
#include <new>
#include <stdexcept>
#include <vector>
#include <cstddef>
#include <iostream>
#include <iomanip>
#include <typeinfo>
#include "PSIMAGAssert.h"
#include "BLAS.h"
#include "LAPACK.h"
#include "Tag.h"
#include "MatrixSlice.h"
#include "IndexedMatrix.h"

namespace psimag {

  template<typename ScalarType>
  class MatrixScalarPrinter {
  public:
    static
    void print(std::ostream& os, const ScalarType& scalar, int width) {
      os << " " << std::setw(width) << scalar;
    };
  };

  template<typename ScalarType>
  class MatrixScalarPrinter<std::complex<ScalarType> > {
  public:
    static
    void print(std::ostream& os, const std::complex<ScalarType>& scalar, int width) {
      os << " [" << std::setw(width) << scalar.real() << ", " << std::setw(width) << scalar.imag() << "]";
    };
  };

  /** \brief A flexible-size matrix class
   *
   * \author Thomas Schulthess
   */
  template <class T>
  class Matrix {
    //    template <class S> 
    //    friend std::istream& operator >> (std::istream& is, Matrix<S>& A);
    //    template <class S>
    //    friend std::ostream& operator << (std::ostream& os, Matrix<S>& A);
  public:
    typedef Matrix<T>    ThisType;
    typedef size_t       size_type;
    typedef T            value_type;
    typedef T*           pointer;
    typedef T*           iterator;
    typedef const T*     const_iterator;
    typedef T&           reference;
    typedef const T&     const_reference;

    Matrix() : nRow(0), nCol(0), lDim(0), owner(true), data(0),col(0) {};

    Matrix(size_type nRows,size_type nCols,size_type ldim=0,const T& val=T()) {
      if(ldim<nRows) 
	ldim=nRows;
      if(bool(nRows*nCols)) 
	init(true,nRows,nCols,ldim);
      else { 
	nRow=0; nCol=0; lDim=0; owner=true; data=0; col=0; 
      } 
      if(val==T(0))
	memset(data,0,sizeof(T)*lDim*nCol);
      else
	for(size_type i=0;i<nCol*lDim;i++) 
	  data[i]=val;
    }

    /** matrix that does not own data (user is responsible for proper dim. */
    Matrix(size_type nRows, size_type nCols, T* dat, size_type ldim=0) 
      : data(dat) {
      if(ldim<nRows) ldim=nRows;
      if(bool(nRows*nCols)) init(false,nRows,nCols,ldim);
      else { nRow=0; nCol=0; lDim=0; owner=false; col=0; }
    }

    /** makes a deep copy of mat */
    Matrix(const Matrix<T> & mat) 
      : nRow(mat.nRow), nCol(mat.nCol), lDim(mat.nRow),owner(true) {
      if(bool(nRow*nCol)) {
	init(true,nRow,nCol,lDim);
	if(mat.lDim==lDim)
	  memcpy(data,mat.data,sizeof(T)*nCol*lDim);
	else
	  for(size_type j=0;j<nCol;j++) for(size_type i=0;i<nRow;i++) 
					  col[j][i] = mat.col[j][i];
      }
      else { nRow=0; nCol=0; lDim=0; owner=true; data=0; col=0; } 
    }

    /*! Makes a deep copy of something Matrix-like. 
      
      A type is matrix like if it has the memebr functions:
      - n_row
      - n_col
      - operator () (i,j)
      
      \author MSS 
    */
    template<typename MatrixLikeType>
    Matrix(const MatrixLikeType& mat):
      nRow(mat.n_row()), 
      nCol(mat.n_col()), 
      lDim(mat.n_row()),
      owner(true) 
    {
      if(bool(nRow*nCol)) {
	init(true,nRow,nCol,lDim);
	for(size_type j=0;j<nCol;j++) 
	  for(size_type i=0;i<nRow;i++) 
	    col[j][i] = mat(i,j);
      }
      else { 
	nRow=0; nCol=0; lDim=0; owner=true; data=0; col=0; 
      } 
    }

    ~Matrix() { kill(); }
    
    iterator       begin() { return data; }
    const_iterator begin() const { return data; }
    iterator       end()   { return data+(nCol-1)*lDim+nRow; }
    const_iterator end() const { return data+(nCol-1)*lDim+nRow; }
    T&             first() { return data[0]; }
    const T&       first() const {return data[0]; }
    T&             last() { return col[nCol-1][nRow-1]; }
    const T&       last() const { return col[nCol-1][nRow-1]; }
    
    size_type size() const { return nRow*nCol; }
    size_type max_size() const { return lDim*nCol; }

    void resize(size_type m,size_type n,size_type ldim=0) {
      if(ldim<m) ldim=m;
      if(ldim*n <= nCol*lDim) {
	if(bool(ldim*n) && data) {
	  //if(n>nCol) { delete [] col; col = new (T*)[n]; } // not ISO compliant?
	  if(n>nCol) { delete [] col; col = new T*[n]; }
	  nCol = n; nRow = m; lDim = ldim;
	  T* tmp = data;
	  for(size_type i=0; i<nCol; i++) { col[i] = tmp; tmp += ldim; }
	}
	else {
	  if(col) delete [] col;
	  nRow=0; nCol=0; lDim=0;
	}
      }
      else {
	kill(); init(owner,m,n,ldim);
      }
    }
    size_type n_row() const { return nRow; }
    size_type n_col() const { return nCol; }
    size_type l_dim() const { return lDim; }
    
    T& operator() (size_type i, size_type j) {
      ASSERT(i<nRow,
	     std::range_error("i>=n_row in T& Matrix<T>::operator()(size_type i,size_type j"));
      ASSERT(j<nCol,
	     std::range_error("j>=n_col in T& Matrix<T>::operator()(size_type i,size_type j"));
      return col[j][i];
    }
    
    const T& operator() (size_type i, size_type j) const {
      ASSERT(i<nRow,
	     std::range_error("i>=n_row in const T& Matrix<T>::operator()(size_type i,size_type j"));
      ASSERT(j<nCol,
	     std::range_error("j>=n_col in const T& Matrix<T>::operator()(size_type i,size_type j"));
      return col[j][i];
    }
  
    bool operator == (const Matrix<T>&);
    bool operator != (const Matrix<T>&); 

    void print(std::ostream& os) const {
      for (size_t o = 0; o < nRow; o++) {
	for (size_t p = 0; p < nCol; p++) 
	  os << (*this)(o,p) << " ";
	os << std::endl;
      }
    }
    
    void printMemoryOrder(std::ostream& os) const {
      for (size_t i = 0; i < size(); i++) {
	os << data[i] << " ";
	os << std::endl;
      }
    }
    
    // ======================================================================  MSS

	void toJSN(std::ostream& os, std::string title="Matrix", int width=13, bool showIndices= false,bool pythonFormat=true) const 
	{
      
		os.precision(13);
		os 	<< std::fixed;
		os 	<< "{ \"tile\": \""    << title        << "\", \n"
			<< " \"rows\":"    << nRow << ", \n"
			<< " \"cols\":"    << nCol << ", \n"
			<< " \"data\": ";
		if (pythonFormat) os<<"array(";
		os<<"[\n";
      
		for(size_t i=0; i<nRow; i++) {
	//os << "(" << i << "," << domain1[i] << "), "
			if (showIndices) os << std::setw(3) << i << ",";
			os << " [";

			for(size_t j=0; j<nCol; j++) {
				const T& val = (*this)(i,j);    
				MatrixScalarPrinter<T>::print(os,val, width);
				//os << " " << std::setw(width) << (*this)(i,j);
				if (j<nCol-1) os << ",";
			}
			os << "]";
			if (i<nRow-1) os<<",";
			os << "\n";
		}
		os << "]";
		if (pythonFormat) os<<")";
		os<<"\n}\n";
	}
    
    // ======================================================================  MSS

    void writeVTK(std::string fileName,
		  std::string version    = "0.1",
		  std::string byte_order = "BigEndian") const{
      
      std::ostringstream msg;
      msg << "================================================= Later we will write a VTK file for " << fileName << "\n";    
      throw std::logic_error(msg.str());
    }

    /* operator= added by G.A. */
    Matrix<T> &operator=(const Matrix<T>& mat)
    {
      if (this == &mat) return *this;   // Gracefully handle self assignment[12.1]
      kill();
      nRow=mat.n_row(), nCol=mat.n_col(), lDim=mat.n_row(),owner=true;
      if(bool(nRow*nCol)) {
	init(true,nRow,nCol,lDim);
	if(mat.lDim==lDim)
	  memcpy(data,mat.data,sizeof(T)*nCol*lDim);
	else
	  for(size_type j=0;j<nCol;j++) for(size_type i=0;i<nRow;i++)
					  col[j][i] = mat.col[j][i];
      }
      else { nRow=0; nCol=0; lDim=0; owner=true; data=0; col=0; }
      return *this;
    }

    /* operator[] added by M.S. */
    template<typename VectorType> 
    VectorType& getRow(size_t row, VectorType& result) const {
      for (size_t col=0; col < (*this).n_col(); col++) {
	result[col] = (*this)(row,col);
      }
      return result;
    }
    
    /* operator[] added by M.S. */
    RowSlice<ThisType> operator [] (size_t row) {
      return RowSlice<ThisType>(*this,row);
    }
    
    /* operator[] added by M.S. */
    const RowSlice<const ThisType> operator [] (size_t row) const {
      return RowSlice<const ThisType>(*this,row);
    }
    
    RowSlice<ThisType> getRow(size_t rowIndex)  {
      return RowSlice<ThisType>(*this,rowIndex);
    }
    
    RowSlice<ThisType> getCol(size_t colIndex)  {
      return ColSlice<ThisType>(*this,colIndex);
    }
    
    IndexedMatrix<ThisType,const std::vector<size_t> > 
    operator () (const std::vector<size_t> rowIndexes,
		 const std::vector<size_t> colIndexes) const {
      return IndexedMatrix<ThisType, const std::vector<size_t> >(*this,rowIndexes,colIndexes);
    }
    
    /* operator added by M.S. */
    inline
    ThisType& operator = (const T& val)
    {
      for(size_type i=0;i<nCol*lDim;i++) 
	data[i]=val;
      return *this;
    }

   /* operator added by M.S. */
    inline
    ThisType& operator += (const ThisType& other)
    {
      // domain checking ??? 
      for(size_type i=0;i<nCol*lDim;i++) 
	data[i] += other.data[i];
      return *this;
    }

    /* added by M.S. */
    inline
    ThisType& squareElements()
    {
      for(size_type i=0;i<nCol*lDim;i++) 
	data[i] = data[i] * data[i];
      return *this;
    }

    /* added by M.S. */
    inline
    ThisType& operator /= (const T& val)
    {
      for(size_type i=0;i<nCol*lDim;i++) 
	data[i] /= val;
      return *this;
    }
	
	    /* added by G.A. */
    inline
    ThisType& operator *= (const T& val)
    {
      for(size_type i=0;i<nCol*lDim;i++)
        data[i] *= val;
      return *this;
    }

    /* added by M.S. */
    inline
    T& atOffset(size_t offset)
    {
      return data[offset];
    }

    /* operator added by M.S. */
    inline
    const T& atOffset(size_t offset) const 
    {
      return data[offset];
    }

    /* operator added by M.S. */
    template<typename Formatter>
    inline
    void load(std::string directory, int number)
    {
      // from the filename from directory and number
      std::ostringstream fileNameBuff;
      fileNameBuff << "./" << directory << "/" << number;
     
     typename Formatter::In formatter(fileNameBuff.str());

      // read the size
      int irows=formatter.irows,icols=formatter.icols;
       this->resize(irows,icols);
       
       formatter.read((char *)data);

    }

    /* operator added by M.S. */
    template<typename Formatter>
    inline
    void save(std::string directory, int number)
    {
      // from the filename from directory and number
      std::ostringstream fileNameBuff;
      fileNameBuff << "./" << directory << "/" << number;
	
     typename Formatter::Out formatter(fileNameBuff.str(),nRow,nCol);
     
     formatter.write((char *)data);
     
    }

  private:
    //Matrix<T>& operator = (const Matrix<T>&);
    size_type nRow,nCol,lDim;
    bool owner;
    T* data;
    T** col;
    void init(bool own,size_type m,size_type n, size_type ld) 
    {
      // init assumes that m*n is non zero and that ld>=m
      ASSERT(bool(m*n)&&bool(ld>=m),
	     std::runtime_error("n*m=0 || ld>=m in private Matrix<T>::init()"));
      owner = own;
      if(own) {
	try {
	  data = new T[m*n];
	}
	catch(std::bad_alloc& e) {
	  std::ostringstream msg;
	  msg << e.what() << "\n"
	      << "ERROR in Matrix<T>: Bad_alloc! failed to allocate memory for data\n"
	      << "new " << typeid(T).name() << "[" << m << "*" << n << "=" << n*m << "];\n";
	  throw std::runtime_error(msg.str());
	}
      }
      else {
	if(!data) {
	  throw std::runtime_error("ERROR in Matrix<T>: data array is empty");
	}
      }
      try {
	// col = new (T*)[n];  // not ISO compliant?
	col = new T*[n];
      }
      catch(std::bad_alloc& e) {
	std::ostringstream msg;
	msg << e.what() << "\n"
	    << "Bad_alloc ERROR in Matrix<T>: failed to allocate memory for columns\n";
	throw std::runtime_error(msg.str());
      }
      nRow = m;
      nCol = n;
      lDim = ld;
      if(data && col) {
	T* tmp = data;
	for(size_type i=0; i<nCol; i++) { col[i] = tmp; tmp += lDim; }
      }
    }
    void kill () {
      if(owner && data) delete [] data;
      if(col) delete [] col;
    }

  };

  template <class T>
  std::istream& operator >> (std::istream& is, Matrix<T>& A) {
    typedef typename Matrix<T>::size_type size_type;
    size_type m,n,l;
    is >> m >> n >> l;
    if(is) {
      A.resize(m,n,l);
      for (size_type j=0; j<A.n_col(); j++) for (size_type i=0; i<A.n_row(); i++) {
	  is >> A(i,j) >> std::endl;
	}
    }
    if(!is) {
      throw std::range_error("ERROR istream& operator >> (std::istream&, Matrix<T>&): read past end stream");
    }
    return is;
  }

	template<class T>
	std::ostream &operator<<(std::ostream &os,Matrix<T> const &A)
	{
		size_t i,j;
		os<<A.n_row()<<" "<<A.n_col()<<"\n";
		for (i=0;i<A.n_row();i++) {
			for (j=0;j<A.n_col();j++) os<<A(i,j)<<" ";
			os<<"\n";
		}
		
		return os;
	}

  template<typename T>
  bool Matrix<T>::operator == (const Matrix<T>& m)
  {
    if( m.n_row()!=this->n_row() ) return false;
    if( m.n_col()!=this->n_col() ) return false;
    const_iterator titr = this->begin(); 
    for(const_iterator mitr=m.begin(); mitr!=m.end(); mitr++)
      {
	if( *mitr!=*titr ) return false;
	titr++;
      }
    return true;
  }

  template<typename T>
  bool Matrix<T>::operator != (const Matrix<T>& m)
  {
    if( m.n_row()!=this->n_row() ) return true;
    if( m.n_col()!=this->n_col() ) return true;
    const_iterator titr = this->begin(); 
    for(const_iterator mitr=m.begin(); mitr!=m.end(); mitr++)
      {
	if( *mitr!=*titr ) return true;;
	titr++;
      }
    return false;;
  }

  // C = alpha * op(A) * op(B) + beta * C
  template<class T> inline
  void GEMM(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C, 
	    const T& alpha=1, const T& beta=0, char opA = 'N', char opB = 'N') 
  {
    ASSERT(A.n_row() == C.n_row(),
	   std::range_error("A.n_row() != C.n_row() in GEMM for Matrix<T>"));
    ASSERT(A.n_col() == C.n_col(),
	   std::range_error("A.n_col() != C.n_col() in GEMM for Matrix<T>"));
    ASSERT(A.n_col() == B.n_row(),
	   std::range_error("A.n_col() != B.n_row() in GEMM for Matrix<T>"));
    BLAS::GEMM(char(opA),
               char(opB),
               static_cast<int>(A.n_row()),
               static_cast<int>(B.n_col()),
               static_cast<int>(A.n_col()),
               alpha,
	       &A(0,0),
               static_cast<int>(A.l_dim()),
               &B(0,0),
               static_cast<int>(B.l_dim()),
               beta,
               &C(0,0),
               static_cast<int>(C.l_dim()));
  }

  template<class T>
  class GESV {
  public:
    int operator () (Matrix<T>& A, Matrix<T>& B) {
      ASSERT(A.n_row()==A.n_col(),
	     std::range_error("A.n_row() != A.n_col() in GESV for Matrix<T>"));
      ASSERT(A.n_row()==B.n_row(),
	     std::range_error("A.n_row() != B.n_row() in GESV for Matrix<T>"));
      pivots.resize(A.n_row());
      int info;
      LAPACK::GESV(static_cast<int>(A.n_row()),
		   static_cast<int>(B.n_row()),
		   &A(0,0),
		   static_cast<int>(A.l_dim()),
		   &pivots[0],
		   &B(0,0),
		   static_cast<int>(B.l_dim()),
		   info);
      return info;
    }
  private:
    std::vector<int> pivots;
  };

//#ifdef PSIMAG_MATRIX_INVERT
  template<class T>
  class INVERT {
  public:
    int operator () (const Matrix<T>& A, Matrix<T>& INV) {

      ASSERT(A.n_row()==A.n_col(),
	     std::range_error("A.n_row() != A.n_col() in INVERT for Matrix<T>"));
      ASSERT(A.n_row()==INV.n_row(),
	     std::range_error("A.n_row() != INV.n_row() in INVERT for Matrix<T>"));
      ASSERT(A.n_col()==INV.n_col(),
	     std::range_error("A.n_col() != INV.n_col() in INVERT for Matrix<T>"));

      INV = A;
      int N  (static_cast<int>(A.n_row()));
      int LDA(static_cast<int>(A.l_dim()));

      pivots.resize(A.n_row());

      int info;

      // GET LU Decomposition which will be stored in INV
      LAPACK::DGETRF(N,N,&INV(0,0),LDA,&pivots[0],info);
      
      if (info != 0) {
	std::ostringstream message;
	message << "INVERT(A,B): GETRF failed! \n";
	if (info < 0) 
	  message << "The " << -info << "th argument had an illegal value!\n";
	else
	  message << "U("<<info<<","<<info<<") was exactly zero. \n";
	throw std::logic_error(message.str());
      }

      int    lwork(-1);
      double blocksize;

      // Get the optimal block size
      LAPACK::DGETRI(N,&INV(0,0),LDA,&pivots[0],blocksize,-1,info);
      
      lwork  = static_cast<int>(blocksize)*N;

      std::vector<double> work(lwork);
      
      LAPACK::DGETRI(N,&INV(0,0),LDA,&pivots[0],&work[0],lwork,info);
 
      if (info != 0) {
	std::ostringstream message;
	message << "INVERT(A,B): GETRI failed! \n";
	if (info < 0) 
	  message << "The " << -info << "th argument had an illegal value!\n";
	else
	  message << "U("<<info<<","<<info<<") was exactly zero. \n";
	throw std::logic_error(message.str());
      }
      
      return info;
    }
  private:
    std::vector<int> pivots;
  };
//#endif
  
  /** \ingroup ostream
   *
   * Output stream operator for SuperCrystal
   */
  template<typename Field>
  Tag toXML(const Matrix<Field>& matrix,
	    std::string name="Matrix") {
    Tag tag(name);
    tag["rows"] = matrix.n_row();
    tag["cols"] = matrix.n_col();
    tag.content << std::endl;
    tag.content << std::endl;
    matrix.print(tag.content);

    return tag;
  }
  
} /* namespace psimag */


#endif /* PSIMAG_Matrix_H */
