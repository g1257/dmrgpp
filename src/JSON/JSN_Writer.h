//-*-C++-*-

/** \ingroup DCA */
/*@{*/

/*! \file JSN_Writer.h

*/

#ifndef dca_JSN_Writer_H
#define dca_JSN_Writer_H

#include "VectorLike.h"

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
      std::string key;

      KeyReference(std::string k, JSN& w):
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

      KeyReference operator[](std::string subKey) { return writer.getComponentWriter(key)[subKey]; }
    };

    //======================================================================

    KeyReference operator[](std::string key) { return KeyReference(key,*this); }

    //======================================================================

    ThisType& getComponentWriter(std::string key) { return writers[key]; }

    //======================================================================

    mutable int  offset;
    int  precis;
    mutable bool printedFirstLine;

    MapType<std::string, bool>                          bools;
    MapType<std::string, std::string>                   strings;
    MapType<std::string, double>                        numbers;
    MapType<std::string, std::vector<double> >          vectors;
    MapType<std::string, std::vector<int> >             intvectors;
    MapType<std::string, std::vector<size_t> >          sizetvectors;
    MapType<std::string, std::vector<std::string> >     stringLists;
    MapType<std::string, const psimag::Matrix<int>*>    intMatrices;
    MapType<std::string, const psimag::Matrix<double>*> dblMatrices;
    MapType<std::string, const psimag::Matrix<std::complex<double> >*> cMatrices;

    MapType<std::string, ThisType>                   writers;
    
    //======================================================================

    void add(std::string key, const bool&          b)   { bools  [key] = b;    }
    void add(std::string key, const char*        str)   { strings[key] = str;  }
    void add(std::string key, const std::string& str)   { strings[key] = str;  }
    void add(std::string key, const size_t& val)              { numbers[key] = val;  }
    void add(std::string key, const int&    val)              { numbers[key] = val;  }
    void add(std::string key, const double& val)              { numbers[key] = val;  }

    void add(std::string key, const std::vector<double>&    vals) { vectors[key]     = vals; }
    void add(std::string key, const std::vector<int>&       vals) { intvectors[key]  = vals; }
    void add(std::string key, const std::vector<size_t>&    vals) { sizetvectors[key]= vals; }
    void add(std::string key, const psimag::Matrix<int>&    mat ) { intMatrices[key] = &mat; }
    void add(std::string key, const psimag::Matrix<double>& mat ) { dblMatrices[key] = &mat; }

    //======================================================================

    template<class T>
    void add(std::string key, const T& obj){ 
      obj.toJSN(getComponentWriter(key));
    }

    //======================================================================

    template<typename T>
    void add(std::string key, const std::map<std::string, T>& map) {
      typedef typename std::map<std::string, T>::const_iterator itr;
      ThisType& writer = writers[key]; 
      for(itr i=map.begin(); i!= map.end(); i++) 
	writer.add(i->first,i->second);
    }
    
    //======================================================================

    template<typename ValType>
    static
    int maxKeyWidth( const MapType<std::string, ValType>& map) {
      size_t result = 0;
      typedef typename MapType<std::string, ValType>::const_iterator itr;
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
		    const MapType<std::string, ValType>& map,
		    int offset, 
		    int keyWidth,
		    Position position=Middle) const {
      
      typedef typename MapType<std::string, ValType>::const_iterator itr;
      
      for(itr i=map.begin(); i!= map.end(); i++) {
	printLine(os,i->first,i->second,offset,keyWidth,printedFirstLine);
	printedFirstLine = true;
      }
    }

    //======================================================================
    
    template<typename ValType>
    static
    void printLine(std::ostream& os,
		   const std::string& key,
		   const ValType&     val,
		   int  offset, 
		   int  keyWidth,
		   bool printedFirstLine_,
		   bool printVal=true)  {

      if (printedFirstLine_) 
	os << ",\n" << std::string(offset,' ') ;
      
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
		   const std::string& val) {
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
		   const std::vector<T>& val) {
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
    std::string toString(const std::map<std::string, std::string>& map) {

      std::ostringstream buff;
      buff << "{";
      std::map<std::string,std::string>::const_iterator itr;
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
    std::string toString(const std::map<std::string, T>& map) {
      typedef std::map<std::string,T> MapType;
      typedef typename MapType::const_iterator CITR;

      std::ostringstream buff;
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
    std::string toString(const std::map<std::string, std::vector<T> >& map) {

      typedef std::map<std::string,T> MapType;
      typedef typename MapType::const_iterator CITR;

      std::ostringstream buff;
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
    void printVector(std::vector<T> vec,
		     std::string title,
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
    std::string quoted(std::string str) {
      std::ostringstream result;
      result << "\"" << str << "\"";
      return result.str();
    }
    
    //======================================================================
   
  };
  
} // end namespace DCA


/*@}*/
#endif
