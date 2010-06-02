//-*-C++-*-

#ifndef PSIMAG_VectorLike_H
#define PSIMAG_VectorLike_H

#include <new>
#include <stdexcept>
#include <vector>
#include <cstddef>
#include <iostream>
#include <iomanip>

#include "PSIMAGAssert.h"
#include "BLAS.h"
#include "LAPACK.h"
#include "Tag.h"

namespace psimag {
  
  /*! \brief Template functions for types that provide member functions:
    
    - size, and 
    - operator[]
      
    and support:     
      
    typedef typename VectorLikeType::value_type FieldType;
      
    \author Mike Summers
  */
  namespace VectorLike {
    
    template<typename VectorLikeType>
    std::string toString(const VectorLikeType& vector,
			 int width=13) {
      std::ostringstream buff;
      printList(vector,buff);
      return buff.str();
    }
    
    template<typename VectorLikeType>
    void print(const VectorLikeType& vector,
	       std::ostream& os,
	       int width=13) {
      for(size_t i=0; i<vector.size(); i++) 
	os << " " << std::setw(width) << vector[i];
      os << "\n";
    }
    
    template<typename VectorLikeType>
    void printList(const VectorLikeType& vector,
		   std::ostream& os,
		   int width=13) {
      os << "[ ";
      
      for(size_t i=0; i<vector.size(); i++) 
	os  << std::setw(width) << vector[i] << ", ";
      
      os << "]";
    }
    
    template<typename VectorLikeType>
    void printArray(const VectorLikeType& vector,
		    std::ostream& os,
		    int width=13) {
      os << "array(";
      printList(vector,os,width);
      os << ")";
    }
    
    /** \ingroup ostream
     *
     * Output stream operator for SuperCrystal
     */
    template<typename VectorLikeType>
    Tag toXML(const VectorLikeType& vector,
	      std::string name="VectorLike") {
      Tag tag(name);
      tag["size"] = vector.size();
      tag.content << std::endl;
      print(vector,tag.content);
      return tag;
    }
    
    // ======================================================================  
    
    template<typename VectorLikeType>
    void toJSN(const VectorLikeType& vector,
	       std::ostream& os, 
	       std::string title="VectorLike",
	       int width=13) {
      
      os.precision(width);
      os << std::fixed;
      os << "{ \'tile\': \'"    << title        << "\', \n"
	 << " \'size\':"    << vector.size() << ", \n"
	 << " \'data\': ";
      printArray(vector,os,width);
      os << "\n}\n";
    }
    
    // ======================================================================  
    
    template<typename VectorLikeType>
    void writeVTK(const VectorLikeType& vector,
		  std::string fileName,
		  std::string version    = "0.1",
		  std::string byte_order = "BigEndian") {
      
      throw std::logic_error("writeVTK ===== Later we will write a VTK file for \n");    
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

    template<typename VectorLikeType> 
    inline
    VectorLikeType& increment(VectorLikeType& vector,
			      typename VectorLikeType::value_type val) {
      for (size_t i=0; i < vector.size(); i++) 
	vector[i] += val;
      return vector;
    }
    
    template<typename VectorLikeType> 
    inline
    VectorLikeType& decrement(VectorLikeType& vector,
			      typename VectorLikeType::value_type val) {
      for (size_t i=0; i < vector.size(); i++) 
	vector[i] -= val;
      return vector;
    }
    
    template<typename VectorLikeType> 
    inline
    VectorLikeType& times(VectorLikeType& vector,
			  typename VectorLikeType::value_type val) {
      for (size_t i=0; i < vector.size(); i++) 
	vector[i] *= val;
      return vector;
    }
    
    template<typename VectorLikeType> 
    inline
    typename VectorLikeType::value_type squaredLength(const VectorLikeType& vector) {
      typename VectorLikeType::value_type result=0;
      for (size_t i=0; i < vector.size(); i++) 
	result += vector[i] * vector[i];
      return result;
    }
    
    template<typename VectorLikeType> 
    inline
    typename VectorLikeType::value_type length(const VectorLikeType& vector) {
      return sqrt(squaredLength(vector));
    }
    
    template<typename VectorLikeType> 
    inline
    VectorLikeType& divide(VectorLikeType& vector,
			   typename VectorLikeType::value_type val) {
      for (size_t i=0; i < vector.size(); i++) 
	vector[i] /= val;
      return vector;
    }
    
    template<typename VectorLikeType> 
    inline
    VectorLikeType& squareElements(VectorLikeType& vector) {
      for (size_t i=0; i < vector.size(); i++) 
	vector[i] *= vector[i];
      return vector;
    }
    
    template<typename VectorLikeType1,
	     typename VectorLikeType2> 
    inline
    std::vector<typename VectorLikeType1::value_type>
    difference(const VectorLikeType1& vector1,
	       const VectorLikeType2& vector2) {
      std::vector<typename VectorLikeType1::value_type> result(vector1.size());
      for (size_t i=0; i < vector1.size(); i++) 
	result[i] = vector1[i] - vector2[i];
      return result;
    }
    
    template<typename VectorLikeType1,
	     typename VectorLikeType2> 
    inline
    std::vector<typename VectorLikeType1::value_type>
    sum(const VectorLikeType1& vector1,
	const VectorLikeType2& vector2) {
      std::vector<typename VectorLikeType1::value_type> result(vector1.size());
      for (size_t i=0; i < vector1.size(); i++) 
	result[i] = vector1[i] + vector2[i];
      return result;
    }
    
    template<typename VectorLikeType1,
	     typename VectorLikeType2> 
    inline
    void accumulate(VectorLikeType1& vector1,
		    const VectorLikeType2& vector2) {
      for (size_t i=0; i < vector1.size(); i++) 
	vector1[i] += vector2[i];
    }
    
    template<typename VectorLikeType1,
	     typename VectorLikeType2> 
    inline
    VectorLikeType2& copy(const VectorLikeType1& vector,
			  VectorLikeType2& targetVector) {
      for (size_t i=0; i < targetVector.size(); i++) 
	targetVector[i] = vector[i];
      return targetVector;
    }
    
    template<typename VectorLikeType1,
	     typename VectorLikeType2> 
    inline
    bool equals(const VectorLikeType1& vector1,
		const VectorLikeType2& vector2) {
      for (size_t i=0; i < vector1.size(); i++) 
	if (vector2[i] != vector1[i]) 
	  return false;
      return true;
    }
    
  } /* namespace VectorLike  */
  
} /* namespace psimag */


#endif 
