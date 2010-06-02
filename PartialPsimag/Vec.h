//-*-C++-*-

#ifndef PSIMAG_Vec_H
#define PSIMAG_Vec_H

#include <iostream>
#include <cmath>
#include <stdexcept>

#include<cstddef>
#include "PSIMAGAssert.h"
#include "MatTraits.h"

/** \file Vec.h
 *  \author Thomas C. Schulthess
 *
 * \note The use of Traits and metaprogramming is only started. Much
 *       of the functionality here is already implemented in the Mat
 *       directory (or should be). Mike
 *
 **/

namespace psimag {
  
  template <typename T, size_t DIM, typename VecTraits> class Vec;
 
}

// These cause problems later on
//  template <class T, size_t DIM>
//  std::istream& operator >> (std::istream& is, psimag::Vec<T,DIM>& v);
//  
//  template <class T, size_t DIM>
//  std::ostream& operator << (std::ostream& os, const psimag::Vec<T,DIM> &v);

namespace psimag {
  
  // --------------------------------------------------------------------------

  /*   These methods are overloaded to provide proper behavior for scalars.
   *   Having these methods defined allows many algorithms that work for 
   *   Heisenberg models to also work for Ising models.
   */

  /**  L2Norm for scalar returns the absolute value of the scalar.
    *  This function is defined to allow generic algorithms that only need
    *  to know the length of a vector work equally well with scalar
    *  quantities. This is frequently needed for algorithms that are 
    *  equally valid for Heisenberg and Ising models.
    *  There is a good implementation question here: are these specializations
    *  better replaced by a true templated function that returns the absolute
    *  value, so that only the Vec forms are specializations?
    *  @author G. Brown, ORNL, September 2004
    */
  template<typename T> T L2Norm(T a) { return std::abs(a); }
#if 0
  template<> float L2Norm<float>(float a)          { return std::abs(a); }
  template<> double L2Norm<double>(double a)       { return std::abs(a); }
  template<> short L2Norm<short>(int a)            { return std::abs(a); }
  template<> int L2Norm<int>(int a)                { reutrn std::abs(a); }
  template<> long L2Norm<long>(int a)              { return std::abs(a); }
  template<> unsigned int L2Norm<unsigned int>(int a}        { return a; } 
  template<> unsigned short L2Norm<unsigned short>(int a)    { return a; }
  template<> unsigned long L2Norm<unsigned long>(int a)      { return a; }
#endif

  // --------------------------------------------------------------------------

  /** \class Vec
   *  \brief Cartesian vector with fixed dimensions. 
   *
   *  It is implemented as a template to allow different types of elements, and
   *  the lower dimensions are implemented through partial specialization
   *  for better efficiency.
   *  \ingroup VecGroup
   *  @author T. Schulthess, G. Brown, M. Summers
   *  @date October 2000
   */
  template < typename T, 
	     size_t   DIM, 
	     typename TRAITS=ColMajorTraits<T,DIM,1> > 
  class Vec {
    // Friends:

    template <typename S, size_t D> friend
    std::istream& operator >> (std::istream& is, Vec<S,D>& v);

    ///
    template <typename S, size_t D> friend
    std::ostream& operator << (std::ostream& os, const Vec<S,D> &v);
    
    
  public:

    // Note Traits and their use are  only partially integrated at this time.
    typedef TRAITS Traits;

    ///
    /** Default constructor, initializes all components to 0.
     */
    Vec() 
    { 
      for(size_t i=0;i<DIM;i++) 
	d[i]=0; 
    }
    
    /** Construct a vector with all components set to value of a.
     */
    explicit Vec(const T& a) 
    { 
      for(size_t i=0;i<DIM;i++) 
	d[i]=a; 
    }
    
    /** Construct a vector with components \f$v_i = a[ {\rm stride} * i]\f$.
     */
    Vec(const T* const a, size_t stride=1)
    { 
      for(size_t i=0;i<DIM;i++) 
	d[i]=a[i*stride]; 
    }
    
    /** This copy constructor makes deep copies. 
     */
    Vec(const Vec<T,DIM>& v)
    { 
      for(size_t i=0;i<DIM;i++) 
	d[i]= v.d[i]; 
    }
    
    /** Create a null default constructor for safety. 
     */
    ~Vec() 
    {}
 
    /** Gives the size of the vector 
    */
    static size_t size () { return DIM; }
    static const size_t static_size=DIM;

    /** Return constant reference to i<SUP>th</SUP> component, note zero i=0,1, ,DIM-1. 
     */
    const T& operator [] (size_t i) const 
    { 
      ASSERT(i<DIM,
	     std::range_error("i>=DIM in Vec<T,DIM>::operator[](size_t i)")); 
      return d[i]; 
    }
    
    /** Return reference to i<SUP>th</SUP> component,
     */
    T& operator [] (size_t i) 
    { 
      ASSERT(i<DIM,
	     std::range_error("i>=DIM in Vec<T,DIM>::operator[](size_t i)"));
      return d[i]; 
    }
    
    /** Assign the value of a to each component.
     */
    Vec<T,DIM>& operator = (const Vec<T,DIM>& v) 
    { 
      for(size_t i=0;i<DIM;i++) 
	d[i] = v.d[i]; 
      return *this; 
    }
    
    /** Assign the value of a to each component.
     */
    Vec<T,DIM>& operator = (const T& a) 
    { 
      for(size_t i=0;i<DIM;i++) 
	d[i] = a; 
      return *this; 
    }

    /** Increment component-wise by vector v.
     */
    Vec<T,DIM>& operator += (const Vec<T,DIM>& v)
    { 
      for(int i=0; i<DIM; i++) d[i] += v.d[i];
      return *this;
    }

    /** Decrement component-wise by vector v.
     */
    Vec<T,DIM>& operator -= (const Vec<T,DIM>& v)
    { 
      for(int i=0; i<DIM; i++) d[i] -= v.d[i];
      return *this;
    }

    /** Scale vector by scalar a.
     */
    Vec<T,DIM>& operator *= (const T& a) 
    {
      for(int i=0; i<DIM; i++) d[i] *= a;
      return *this;
    }

    /** Scale vector by scalar 1/a.
     */
    Vec<T,DIM>& operator /= (const T& a) 
    {
      for(int i=0; i<DIM; i++) d[i] /= a;
      return *this;
    }

 public:

    /// The iterator is just a pointer
    typedef T* iterator;
    
    /// Constant iterator 
    typedef T const* const_iterator;

    /// Iterator directed at first component of vector
    iterator begin() { return &(d[0]); }

    /// Constant iterator directed at first component of vector
    const_iterator begin() const { return &(d[0]); }

    /// Iterator directed at one past the last component of vector
    iterator end() { return &(d[DIM]); }

    /// Constant iterator directed at one past the last component of the vector
    const_iterator end() const { return &(d[DIM]); }


  private:
    
    T d[DIM];
    
  };
  
  ///
  template <typename T, size_t DIM>
  inline std::istream& operator >> (std::istream& is, Vec<T,DIM>& v)
  { 
    for(size_t i=0;i<DIM;i++) 
      is >> v.d[i]; 
    return is; 
  }  
  
  ///
  template <typename T, size_t DIM>
  inline std::ostream& operator << (std::ostream& os, const Vec<T,DIM> &v)
  { 
    os << v.d[0]; 
    for(size_t i=1;i<DIM;i++) 
      os << "\t" << v.d[i]; 
    return os; 
  }
  
  /** L1Norm \f$= \sum_{j=0}^{DIM-1} \|x_i\|\f$.
   */
  template < typename T, size_t DIM > 
  T L1Norm(const Vec <T,DIM> &a )
  { 
    T val=0;
    for(size_t i=0; i<DIM; i++)
      val += std::abs(a[i]);
    return val;
  }
  
  /** L2Norm \f$= sqrt{(\sum_{j=0}^{DIM-1} x_i^2)}\f$.
   */
  template < typename T, size_t DIM > 
  T L2Norm(const Vec <T,DIM> &a)
  {
    T val=0;
    for(size_t i=0; i<DIM; i++)
      val += a[i]*a[i];
    return sqrt(val);
  }
  
  /** Returns scalar product of a and b. 
   */
  template < typename T, size_t DIM >
  T operator * (const Vec<T,DIM>& a, 
		const Vec<T,DIM>& b)
  { 
    T val=0;
    for(size_t i=0; i<DIM; i++)
      val += a[i]*b[i];
    return val;
  }
  
  /** Left-multiplication of vector by scalar a.
   *
   *  /bug Corrected initialization of temporary Vec from val = a to val = b.
   */
  template <typename T, size_t DIM>
  Vec<T,DIM> operator * (const T& a, 
			  const Vec<T,DIM>& b)
  { 
    Vec<T,DIM> val=b;
    for(size_t i=0; i<DIM; i++)
      val[i] *= a;
    return val;
  }

  /** Right-multiplication of vector by scalar a.
   */
  template < typename T, size_t DIM >
  Vec<T,DIM> operator * (const Vec<T,DIM>& b, 
			  const T& a)
  { 
    Vec<T,DIM> val=b;
    for(size_t i=0; i<DIM; i++)
      val[i] *= a;
    return val;
  }

  /** Division of a vector by a scalar a.
   */
  template < typename T, size_t DIM, typename Traits>
  Vec<T,DIM, Traits> operator / (const Vec<T,DIM,Traits>& b, 
				 const T& a)
  { 
    Vec<T,DIM,Traits> val=b;
    for(size_t i=0; i<DIM; i++)
      val[i] /= a;
    return val;
  }

  /** Vector addition.
   */
  template < typename T, size_t DIM >
  Vec<T,DIM> operator + (const Vec<T,DIM>& a, 
			  const Vec<T,DIM>& b)
  { 
    Vec<T,DIM> val;
    for(size_t i=0; i<DIM; i++) 
      val[i] = a[i] + b[i];
    return val;
  }

  /** Vector subtraction.
   */
  template < typename T, size_t DIM >
  Vec<T,DIM> operator - (const Vec<T,DIM>& a, 
			  const Vec<T,DIM>& b)
  { 
    Vec<T,DIM> val;
    for(size_t i=0; i<DIM; i++) 
      val[i] = a[i] - b[i];
    return val;
  }

  /** Negation.
   */
  template < typename T, size_t DIM >
  Vec<T,DIM> operator - (const Vec<T,DIM>& a)
  { 
    Vec<T,DIM> val;
    for(size_t i=0; i<DIM; i++) 
      val[i] = -a[i];
    return val;
  }

  /** Equivalence operator
   */
  template<typename T, size_t DIM>
  bool operator == (const Vec<T,DIM>& a, const Vec<T,DIM>& b)
  {
     //
     // old code, semantics easy to see
     // bool equal = true;
     // for(size_t i=0; equal && i<DIM; i++) if(a[i]!=b[i]) equal=false;
     // return equal;
     //
     // new code, HKL thinks its faster. Depends on optimizer
     for(size_t i=0; i<DIM; i++) if(a[i]!=b[i]) return false;
     return true;
  }

  /** Nonequivalence operator
   */
  template<typename T, size_t DIM>
  bool operator != (const Vec<T,DIM>& a, const Vec<T,DIM>& b)
  {
    //bool nequal = false;
    //for(size_t i=0; !nequal && i<DIM; i++) if(a[i]!=b[i]) nequal=true;
    //return nequal;
    for(size_t i=0; i<DIM; i++) if(a[i]!=b[i]) return true; 
    return false;
  }
 
  /** Ordering operator
   */
  template<typename T, size_t DIM>
  bool operator < ( const Vec<T,DIM>& a, const Vec<T,DIM>& b)
  {
      // bool lessthan = false;
      // bool equal = true;
      // for(size_t i=0; equal && !lessthan && i<DIM; i++) {
      //    if(a[i]<b[i]) lessthan=true;
      //    else if(a[i]>b[i]) equal = false;
      // }
      // return lessthan;
      for(size_t i=0; i<DIM; i++)
      {
        if(a[i]<b[i]) return true;
        if(a[i]>b[i]) return false;
      }
      return false;
  }

  // --------------------------------------------------------------------------
  
  /**  \brief Specialized Cartesian vector in 3D.
   *
   *   All methods of Vec<T,DIM> apply to Vec<T,3>.
   *
   *   \ingroup VecGroup
   *   @author T. Schulthess and G. Brown
   *   @date October 2000
   */
  template < typename T, typename TRAITS > 
  class Vec< T, 3, TRAITS > {

  public:
  
    // Note Traits and their use are  only partially integrated at this time.
    typedef TRAITS Traits;

    ///
    template < typename S> friend 
    std::istream& operator >> (std::istream& is, Vec<S,3>& v);  
    ///
    template < typename S> friend 
    std::ostream& operator << (std::ostream& os, const Vec<S,3> &v);

    /** Default constructor, set all components to zero.
     */
    Vec()
    {
      d[0]=d[1]=d[2]=0;
    }

    /** Construct a vector with all components set to value of a.
     */
    explicit Vec(const T& a) 
    {
      d[0]=d[1]=d[2]=a;
    }

    /** Construct a vector with elements set to (x,y,z)
     */
    Vec(const T& x,const T& y,const T& z) 
    {
      d[0]=x; d[1]=y; d[2]=z;
    }

    /** Construct a vector with components \f$v_i = a[ {\rm stride} * i]\f$.
     */
    Vec(const T * const a,size_t stride=1) 
    { 
      d[0]=a[0]; d[1]=a[stride]; d[2]=a[2*stride]; 
    }

    /** This copy constructor makes deep copies. 
     */
    Vec(const Vec<T,3>& v) 
    { 
      d[0]=v.d[0]; d[1]=v.d[1]; d[2]=v.d[2]; 
    }

    /** Create a null default constructor for safety. 
     */
    ~Vec() 
    {}
 
    /** Gives the size of the vector 
    */
    static size_t size () { return 3; }
    static const size_t static_size=3;

    /** Return constant reference to i<SUP>th</SUP> component, note i=0,1,3. 
     */
    const T& operator [] (size_t i) const 
    { 
      ASSERT(i<3,
	     std::range_error("i>=3 in Vec<T,3>::operator[](size_t i)"));
      return d[i]; 
    }

    /** Return reference to i<SUP>th</SUP> component, note i=0,1,3. 
     */
    T& operator [] (size_t i) 
    { 
      ASSERT(i<3,
	     std::range_error("i>=3 in Vec<T,3>::operator[](size_t i)")); 
      return d[i]; 
    }

    /** This assignment operator makes deep copy.
     */
    Vec<T,3>& operator = (const Vec<T,3>& v)
    { 
      d[0]=v.d[0]; d[1]=v.d[1]; d[2]=v.d[2]; 
      return *this; 
    }
  
    /** Assign the value of a to each component.
     */
    Vec<T,3>& operator = (const T& a) 
    {
      d[0]=d[1]=d[2]=a; 
      return *this;
    }

    /** Increment component-wise by vector v.
     */
    Vec<T,3>& operator += (const Vec<T,3>& v)
    { 
      d[0]+=v.d[0]; d[1]+=v.d[1]; d[2]+=v.d[2]; 
      return *this;
    }

    /** Decrement component-wise by vector v.
     */
    Vec<T,3>& operator -= (const Vec<T,3>& v)
    { 
      d[0]-=v.d[0]; d[1]-=v.d[1]; d[2]-=v.d[2]; 
      return *this;
    }

    /** Scale vector by scalar a.
     */
    Vec<T,3>& operator *= (const T& a) 
    {
      d[0]*=a; d[1]*=a; d[2]*=a; 
      return *this;
    }

 public:

    /// The iterator is just a pointer
    typedef T* iterator;

    /// Constant iterator
    typedef T const* const_iterator;

    /// Iterator directed at first component of vector
    iterator begin() { return &(d[0]); }

    /// Constant iterator directed at first component of vector
    const_iterator begin() const { return &(d[0]); }

    /// Iterator directed at one past the last component of vector
    iterator end() { return &(d[3]); }

    /// Constant iterator directed at one past the last component of the vector
    const_iterator end() const { return &(d[3]); }

  private:

    T d[3];

  };


  ///
  template < typename T > 
  std::istream& operator >> (std::istream& is, Vec<T,3>& v)
  { 
    for(size_t i=0;i<3;i++) 
      is >> v.d[i]; 
    return is; 
  }  


  ///
  template < typename T > 
  std::ostream& operator << (std::ostream& os, const Vec<T,3> &v)
  { 
    os << v.d[0]; 
    for(size_t i=1;i<3;i++) 
      os << "\t" << v.d[i]; 
    return os; 
  }


  /** L1Norm \f$= \sum_{j=0}^{DIM-1} \|x_i\|\f$.
   */
  template < typename T > T L1Norm(const Vec<T,3>& a)
  { 
    return std::abs(a[0])+std::abs(a[1])+std::abs(a[2]); 
  }

  /** L2Norm \f$= sqrt{\sum_{j=0}^{DIM-1} x_i^2}\f$.
   */
  template < typename T > T L2Norm(const Vec<T,3>& a)
  { 
    return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]); 
  }

  /** Returns scalar product of a and b. 
   */
  template < typename T >
  T operator * (const Vec<T,3>& a, const Vec<T,3>& b)
  { 
    return T(a[0]*b[0]+a[1]*b[1]+a[2]*b[2]); 
  }

  /** Returns vector product of a and b 
   */
  template < typename T >
  Vec<T,3> operator % (const Vec<T,3>& a, const Vec<T,3>& b)
  { 
    return Vec<T,3> (a[1] * b[2] - a[2] * b[1],
		     a[2] * b[0] - a[0] * b[2],
		     a[0] * b[1] - a[1] * b[0]);
  }

  /** Left-multiplication of vector by scalar a.
   */
  template < typename T >
  Vec<T,3> operator * (const T& a, const Vec<T,3>& b)
  { 
    return Vec<T,3> (a * b[0], a * b[1], a * b[2]); 
  }

  /** Right-multiplication of vector by scalar a.
   */
  template < typename T >
  Vec<T,3> operator * (const Vec<T,3>& b, const T& a)
  {
    return Vec<T,3> (a * b[0], a * b[1], a * b[2]); 
  }

  /** Division of vector by scalar a.
   */
  template < typename T, typename Traits>
  Vec<T,3,Traits> operator / (const Vec<T,3,Traits>& b, const T& a)
  {
    return Vec<T,3,Traits> (b[0]/a, b[1]/a, b[2]/a); 
  }

  /** Vector addition.
   */
  template < typename T >
  Vec<T,3> operator + (const Vec<T,3>& a, const Vec<T,3>& b)
  { 
    return Vec<T,3> (a[0]+b[0],a[1]+b[1],a[2]+b[2]); 
  }

  /** Vector subtraction.
   */
  template < typename T >
  Vec<T,3> operator - (const Vec<T,3>& a, const Vec<T,3>& b)
  { 
    return Vec<T,3> (a[0]-b[0],a[1]-b[1],a[2]-b[2]); 
  }

  /** Negation.
   */
  template < typename T >
  Vec<T,3> operator - (const Vec<T,3>& a)
  { 
    return Vec<T,3> (-a[0],-a[1],-a[2]); 
  }

  // --------------------------------------------------------------------------

  /**  \brief Specialized to Cartesian vector in 2D.
   *
   *   All methods of Vec<T,DIM> apply to Vec<T,2>.
   *
   *   \ingroup VecGroup
   *   @author T. Schulthess and G. Brown
   *   @date  October 2000
   */
  template < typename T , typename TRAITS > class Vec< T, 2, TRAITS> {

  public:

    // Note Traits and their use are  only partially integrated at this time.
    typedef TRAITS Traits;

    ///
    template <typename S> friend
    std::istream& operator >> (std::istream& is, Vec<S,2>& v);
    ///
    template <typename S> friend
    std::ostream& operator << (std::ostream& os, const Vec<S,2> &v);
  
    /** Default constructor, initializes all components to 0. 
     */
    Vec() 
    {
      d[0]=d[1]=0;
    }

    /** Sets all components to value of a, no implicite conversion. 
     */
    explicit Vec(const T& a) 
    {
      d[0]=d[1]=a;
    }

    /** Initialize 2D vector with (x,y). 
     */
    Vec(const T& x,const T& y) 
    {
      d[0]=x; d[1]=y;
    }

    /** Construct Cartesian vector \f$\{a[0],a[stride]\}\f$. 
     */
    Vec(const T * const a,size_t stride=1) 
    {  
      d[0]=a[0]; d[1]=a[stride]; 
    }
 
    /** Copy constructor makes deep copy.
     */
    Vec(const Vec<T,2>& v) 
    { 
      d[0]=v.d[0]; d[1]=v.d[1]; 
    }

    /** Destructor that does nothing.
     */
    ~Vec() 
    {}
 
    /** Gives the size of the vector 
    */
    static size_t size () { return 2; }
    static const size_t static_size=2;

    /** Return constant reference to i<SUP>th</SUP> component, note that i=0,1. 
     */
    const T& operator [] (size_t i) const 
    { 
      ASSERT(i<2,
	     std::range_error("i>=2 in Vec<T,2>::operator[](size_t i)"));
      return d[i]; 
    }

    /** Return reference to i<SUP>th</SUP> component, note that i=0,1.
     */
    T& operator [] (size_t i) 
    { 
      ASSERT(i<2,
	     std::range_error("i>=2 in Vec<T,2>::operator[](size_t i)"));
      return d[i]; 
    }

    /** Assignment operator makes deep copy.
     */
    Vec<T,2>& operator = (const Vec<T,2>& v)
    { 
      d[0]=v.d[0]; d[1]=v.d[1]; 
      return *this; 
    }

    /** Assign value of a to all fields of vector.
     */
    Vec<T,2>& operator = (const T& a) 
    {
      d[0]=d[1]=a;
      return *this;
    }

    /** Componentwise increment by another vector.
     */
    template<size_t OtherDIM>
    Vec<T,2>& operator += (const Vec<T,OtherDIM>& v)
    { 
      d[0]+=v.d[0]; d[1]+=v.d[1]; 
      return *this; 
    }

    /** Componentwise decrement by another vector.
     */
    Vec<T,2>& operator -= (const Vec<T,2>& v)
    { 
      d[0]-=v.d[0]; d[1]-=v.d[1]; 
      return *this; 
    }

    /** Scale vector by a.
     */
    Vec<T,2>& operator *= (const T& a) 
    {
      d[0]*=a; d[1]*=a; 
      return *this;
    }

 public:

    /// The iterator is just a pointer
    typedef T* iterator;

    /// Constant iterator
    typedef T const* const_iterator;

    /// Iterator directed at first component of vector
    iterator begin() { return &(d[0]); }

    /// Constant iterator directed at first component of vector
    const_iterator begin() const { return &(d[0]); }

    /// Iterator directed at one past the last component of vector
    iterator end() { return &(d[2]); }

    /// Constant iterator directed at one past the last component of the vector
    const_iterator end() const { return &(d[2]); }

  private:

    T d[2];

  };


  ///
  template < typename T > 
  std::istream& operator >> (std::istream& is, Vec<T,2>& v)
  { 
    for(size_t i=0;i<2;i++) 
      is >> v.d[i]; 
    return is; 
  }  


  ///
  template < typename T > 
  std::ostream& operator << (std::ostream& os, const Vec<T,2> &v)
  { 
    os << v.d[0]; 
    for(size_t i=1;i<2;i++) 
      os << "\t" << v.d[i]; 
    return os; 
  }

  /** L1Norm \f$= \|a_0\| + \|a_1\|\f$.
   */
  template < typename T > T L1Norm(const Vec<T,2>& a)
  { 
    return std::abs(a[0])+std::abs(a[1]); 
  }

  /** L2Norm \f$= sqrt{\sum_{j=0}^{DIM-1} x_i^2}\f$.
   */
  template < typename T > T L2Norm(const Vec<T,2>& a)
  { 
    return sqrt(a[0]*a[0]+a[1]*a[1]); 
  }

  /** Returns scalar product of a and b. 
   */
  template < typename T >
  T operator * (const Vec<T,2>& a, const Vec<T,2>& b)
  { 
    return T(a[0]*b[0]+a[1]*b[1]); 
  }

  /** Returns vector product of a and b.
   */
  template < typename T >
  T operator % (const Vec<T,2>& a, const Vec<T,2>& b)
  { 
    return T(a[0] * b[1] - a[1] * b[0]);
  }

  /** Scalar left-multiplication.
   */
  template < typename T >
  Vec<T,2> operator * (const T& a, const Vec<T,2>& b)
  { 
    return Vec<T,2> (a * b[0], a * b[1]); 
  }

  /** Scalar right-multiplication.
   */
  template < typename T >
  Vec<T,2> operator * (const Vec<T,2>& b, const T& a)
  { 
    return Vec<T,2> (a * b[0], a * b[1]); 
  }

  /** Scalar division.
   */
  template < typename T, typename Traits>
  Vec<T,2,Traits> operator / (const Vec<T,2,Traits>& b, const T& a)
  { 
    return Vec<T,2,Traits> (b[0]/a, b[1]/a); 
  }


  /** Vector addition.
   */
  template <typename T, typename Traits>
  Vec<T,2> operator + (const Vec<T,2,Traits>& a, const Vec<T,2,Traits>& b)
  { 
    return Vec<T,2> (a[0]+b[0],a[1]+b[1]); 
  }

  /** Vector subtraction.
   */
  template < typename T >
  Vec<T,2> operator - (const Vec<T,2>& a, const Vec<T,2>& b)
  { 
    return Vec<T,2> (a[0]-b[0],a[1]-b[1]); 
  }

  /** Negation.
   */
  template < typename T >
  Vec<T,2> operator - (const Vec<T,2>& a)
  { 
    return Vec<T,2> (-a[0],-a[1]); 
  }

} /* namspace psimag */

#endif /* PSIMAG_Vec_H */
