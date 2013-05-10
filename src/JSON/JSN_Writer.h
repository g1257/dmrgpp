//-*-C++-*-

/** \ingroup DCA */
/*@{*/

/*! \file JSN_Writer.h
Author: Michael S. Summers

*/

#ifndef dca_JSN_Writer_H
#define dca_JSN_Writer_H

#include "VectorLike.h"
#include "String.h"

namespace dca {

  using namespace psimag::VectorLike;

  template<typename T> class PrintArrayWidth         {public: enum {value=7 }; };
  template<>           class PrintArrayWidth<double> {public: enum {value=13}; };

  template<template<typename,typename> class MapType=std::map>
  class JSN {
  public:

    typedef enum {First,Middle,Last} Position;
    typedef JSN<MapType> ThisType;

    //======================================================================

    JSN(size_t offset_=0):
      offset          (offset_),
      precis          (10),
      printedFirstLine(false)
    {}
    
    //======================================================================

    class KeyReference {
    public:
      JSN& writer;
      String key;

      KeyReference(String k, JSN& w):
	writer(w),
	key   (k)
      {}

      template<typename T>
      void operator= (const T& val) {
	writer.add(key,val);
      }

      operator JSN& () {
	return writer;
      }

      KeyReference operator[](String subKey) { return writer.getComponentWriter(key)[subKey]; }
    };

    //======================================================================

    KeyReference operator[](String key) { return KeyReference(key,*this); }

    //======================================================================

    ThisType& getComponentWriter(String key) { return writers[key]; }

    //======================================================================

    mutable int  offset;
    int  precis;
    mutable bool printedFirstLine;

    MapType<String, bool>                          bools;
    MapType<String, String>                   strings;
    MapType<String, double>                        numbers;
    MapType<String, typename Vector<double>::Type >          vectors;
    MapType<String, typename Vector<int>::Type >             intvectors;
    MapType<String, typename Vector<size_t>::Type >          sizetvectors;
    MapType<String, typename Vector<String>::Type >     stringLists;
    MapType<String, const psimag::Matrix<int>*>    intMatrices;
    MapType<String, const psimag::Matrix<double>*> dblMatrices;
    MapType<String, const psimag::Matrix<std::complex<double> >*> cMatrices;

    MapType<String, ThisType>                   writers;
    
    //======================================================================

    void add(String key, const bool&          b)   { bools  [key] = b;    }
    void add(String key, const char*        str)   { strings[key] = str;  }
    void add(String key, const String& str)   { strings[key] = str;  }
    void add(String key, const size_t& val)              { numbers[key] = val;  }
    void add(String key, const int&    val)              { numbers[key] = val;  }
    void add(String key, const double& val)              { numbers[key] = val;  }

    void add(String key, const typename Vector<double>::Type&    vals) { vectors[key]     = vals; }
    void add(String key, const typename Vector<int>::Type&       vals) { intvectors[key]  = vals; }
    void add(String key, const typename Vector<size_t>::Type&    vals) { sizetvectors[key]= vals; }
    void add(String key, const psimag::Matrix<int>&    mat ) { intMatrices[key] = &mat; }
    void add(String key, const psimag::Matrix<double>& mat ) { dblMatrices[key] = &mat; }

    //======================================================================

    template<class T>
    void add(String key, const T& obj){ 
      obj.toJSN(getComponentWriter(key));
    }

    //======================================================================

    template<typename T>
    void add(String key, const Map<String, T>::Type& map) {
      typedef typename Map<String, T>::Type::const_iterator itr;
      ThisType& writer = writers[key]; 
      for(itr i=map.begin(); i!= map.end(); i++) 
	writer.add(i->first,i->second);
    }
    
    //======================================================================

    template<typename ValType>
    static
    int maxKeyWidth( const MapType<String, ValType>& map) {
      size_t result = 0;
      typedef typename MapType<String, ValType>::const_iterator itr;
      for(itr i=map.begin(); i!= map.end(); i++) 
	if( i->first.length() > result)
	  result = i->first.length();
      return result+2;
    }

    //======================================================================

    int maxKeyWidth() const {
      int result = 0;
      int maxKeyLen = 0; 

      maxKeyLen = maxKeyWidth(strings);
      if (maxKeyLen > result) result = maxKeyLen;

      maxKeyLen = maxKeyWidth(bools);
      if (maxKeyLen > result) result = maxKeyLen;

      maxKeyLen = maxKeyWidth(numbers);
      if (maxKeyLen > result) result = maxKeyLen;

      maxKeyLen = maxKeyWidth(vectors);
      if (maxKeyLen > result) result = maxKeyLen;

      maxKeyLen = maxKeyWidth(intvectors);
      if (maxKeyLen > result) result = maxKeyLen;

      maxKeyLen = maxKeyWidth(sizetvectors);
      if (maxKeyLen > result) result = maxKeyLen;

      maxKeyLen = maxKeyWidth(stringLists);
      if (maxKeyLen > result) result = maxKeyLen;

      maxKeyLen = maxKeyWidth(writers);
      if (maxKeyLen > result) result = maxKeyLen;

      maxKeyLen = maxKeyWidth(intMatrices);
      if (maxKeyLen > result) result = maxKeyLen;

      maxKeyLen = maxKeyWidth(dblMatrices);
      if (maxKeyLen > result) result = maxKeyLen;

      maxKeyLen = maxKeyWidth(cMatrices);
      if (maxKeyLen > result) result = maxKeyLen;

      return result;
    }

    //======================================================================

    template<typename ValType>
    void printLines(std::ostream& os,
		    const MapType<String, ValType>& map,
		    int offset, 
		    int keyWidth,
		    Position position=Middle) const {
      
      typedef typename MapType<String, ValType>::const_iterator itr;
      
      for(itr i=map.begin(); i!= map.end(); i++) {
	printLine(os,i->first,i->second,offset,keyWidth,printedFirstLine);
	printedFirstLine = true;
      }
    }

    //======================================================================
    
    template<typename ValType>
    static
    void printLine(std::ostream& os,
		   const String& key,
		   const ValType&     val,
		   int  offset, 
		   int  keyWidth,
		   bool printedFirstLine_,
		   bool printVal=true)  {

      if (printedFirstLine_) 
	os << ",\n" << String(offset,' ') ;
      
      os << std::setw(keyWidth) << std::left << quoted(key);
      os << " : ";
      if (printVal)
	valString(os,offset,keyWidth,val) ;
    }

    //======================================================================
    //======================================================================

    template<typename T>
    static 
    void valString(std::ostream& os,int offset, int keyWidth,
		   const T& val) {
      os << val;
    }

    template<typename T>
    static 
    void valString(std::ostream& os,int offset, int keyWidth,
		   const psimag::Matrix<T>* val) {

      const psimag::Matrix<T>& matrix(*val);

      int kWidth    = 6;
      int newOffset = offset + keyWidth + 4; 

      os << "{";
      printLine(os,"rows", matrix.n_row(), newOffset, kWidth, false);
      printLine(os,"cols", matrix.n_col(), newOffset, kWidth, true );
      printLine(os,"data", 0,              newOffset, kWidth, true, false);
      psimag::MatrixLike::printArray(matrix,os,PrintArrayWidth<T>::value,newOffset+8);
      os << "}";

    }

    static
    void valString(std::ostream& os, int offset, int keyWidth,
		   const String& val) {
      os << quoted(val);
    }

    static
    void valString(std::ostream& os,int offset, int keyWidth,
		   const ThisType& val) {
      val.print(os, offset + keyWidth + 4);
    }

    template<typename T>
    static
    void valString(std::ostream& os, int offset, int keyWidth,
		   const typename Vector<T>::Type& val) {
      size_t width = maxElementStringWidth(val);
      printList(val,os,',',width);
    }

    //======================================================================

    void print(std::ostream& os,int off=1,bool printLastBrace=true) const {

      offset           = off;
      printedFirstLine = false;

      int keyWidth   = maxKeyWidth();

      os << "{";
      printLines(os,strings,     offset, keyWidth);
      printLines(os,bools,       offset, keyWidth);
      printLines(os,numbers,     offset, keyWidth);
      printLines(os,vectors,     offset, keyWidth);
      printLines(os,sizetvectors,offset, keyWidth);
      printLines(os,intvectors,  offset, keyWidth);
      printLines(os,stringLists, offset, keyWidth);
      printLines(os,intMatrices, offset, keyWidth);
      printLines(os,dblMatrices, offset, keyWidth);
      printLines(os,cMatrices,   offset, keyWidth);
      printLines(os,writers,     offset, keyWidth);
      if (printLastBrace)
	os << "}";
    }

    //======================================================================
      
    static
    String toString(const Map<String, String>::Type& map) {

      PsimagLite::OstringStream buff;
      buff << "{";
      Map<String,String>::Type::const_iterator itr;
      for (itr = map.begin(); itr != map.end(); itr++) {
	if (itr != map.begin())
	  buff << ",";
	buff << " \"" << itr->first << "\": \"" << itr->second << "\"";
      }
      buff << "}";
      return buff.str();
    }

    template<typename T>
    static
    String toString(const Map<String, T>::Type& map) {
      typedef Map<String,T>::Type MapType;
      typedef typename MapType::const_iterator CITR;

      PsimagLite::OstringStream buff;
      buff << "{";
      for (CITR itr = map.begin(); itr != map.end(); itr++) {
	if (itr != map.begin())
	  buff << ",";
	buff << " \"" << itr->first << "\": " << itr->second;
      }
      buff << "}";
      return buff.str();
    }

    template<typename T>
    String toString(const Map<String, typename Vector<T> >::Type::Type& map) {

      typedef Map<String,T>::Type MapType;
      typedef typename MapType::const_iterator CITR;

      PsimagLite::OstringStream buff;
      buff << "{";

      for (CITR itr = map.begin(); itr != map.end(); itr++) {
	if (itr != map.begin())
	  buff << ",";
	buff << " \"" << itr->first << "\": " <<  psimag::VectorLike::toString(itr->second);
      }
      buff << "}";
      return buff.str();
    }

    template<class T>
    inline
    void printVector(typename Vector<T>::Type vec,
		     String title,
		     std::ostream& os)  {
      os.precision(precis);
      os << std::fixed; // scientific;
      os << "{ "
	 << " \'tile\': \'" << title   << "\', \n"
	 << " \'type\': \'Function1D', \n"
	 << " \'size\': " << vec.size() << ", \n";
      os << " \'data\': ";
      os << "array([[";
      for(size_t j=0; j<vec.size(); j++) 
	os << " " << std::setw(precis) << vec[j] << ",\n";
      os << "]])\n";
      
      os << "}";
    }
    //====================================================================== 
    static 
    String quoted(String str) {
      PsimagLite::OstringStream result;
      result << "\"" << str << "\"";
      return result.str();
    }
    
    //======================================================================
   
  };
  
} // end namespace DCA


/*@}*/
#endif
