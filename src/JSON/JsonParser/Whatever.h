//-*-C++-*-

/** \ingroup JsonParser */
/*@{*/

/*! \file JsonParser.h  
 *
 *  
 *
 */

#ifndef  JsonParser_Whatever_H
#define  JsonParser_Whatever_H
#include <complex>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include "CharacterMapper.h"

namespace JsonParser {
  
  template<typename T> class Type; //       { public: static Whatever::WhateverType to(); };
      
  class Whatever {
  public:
    
    typedef enum {
      WHATEVER_MAT,
      WHATEVER_MAP,
      WHATEVER_VECTOR,
      WHATEVER_MATRIX,
      WHATEVER_STRING,
      WHATEVER_INTEGER,
      WHATEVER_DOUBLE,
      WHATEVER_BOOL,
      WHATEVER_NULL,
      WHATEVER_UNKNOWN
    } WhateverType;
  
    typedef std::map<std::wstring,Whatever> WhateverMap;
    typedef std::vector<Whatever>           WhateverVector;

    typedef std::map<std::wstring,const JsonParser::Whatever*> FlatMapType;

    
    //====================================================================== WhateverType
  
    static std::wstring typeName(WhateverType t) {
      switch (t) {
      case WHATEVER_NULL:
	return L"WHATEVER_NULL";
      case WHATEVER_MAT:
	return L"WHATEVER_MAT";
      case WHATEVER_MAP:
	return L"WHATEVER_MAP";
      case WHATEVER_VECTOR:
	return L"WHATEVER_VECTOR";
      case WHATEVER_MATRIX:
	return L"WHATEVER_MATRIX";
      case WHATEVER_STRING:
	return L"WHATEVER_STRING";
      case WHATEVER_INTEGER:
	return L"WHATEVER_INTEGER";
      case WHATEVER_DOUBLE:
	return L"WHATEVER_DOUBLE";
      case WHATEVER_BOOL:
	return L"WHATEVER_BOOL";
      case WHATEVER_UNKNOWN:
	return L"WHATEVER_UNKNOWN";
      default:
	throw std::logic_error("Whatever::typeName given wrong type");
      }
    }
  
    //======================================================================

    static std::string ntypeName(WhateverType t) {
      std::wstring wname(typeName(t));
      return std::string(wname.begin(), wname.end());
    }

    //======================================================================

    WhateverType   type;
    Whatever*      parent;
    std::wstring   valueString;
    WhateverMap    whateverMap;
    WhateverVector whateverVector;
    bool           whateverBool;
    int            whateverInteger;
    double         whateverDouble;
    std::wstring   myKey;
    int            myIndex;
    std::string    filename;
    int            startPos;
    int            endPos;

    //======================================================================

    Whatever():
      type(WHATEVER_UNKNOWN),
      parent(0),
      valueString(),
      whateverMap(),
      whateverVector(),
      whateverBool(true),
      whateverInteger(0),
      whateverDouble(0),
      myKey(L"?"),
      myIndex(-1),
      startPos(0),
      endPos(0)
    {}

    //======================================================================

    Whatever(WhateverType t):
      type(t),
      parent(0),
      valueString((t = WHATEVER_NULL)? L"NULL" : L""),
      whateverMap(),
      whateverVector(),
      whateverBool(true),
      whateverInteger(0),
      whateverDouble(0),
      myKey(L"?"),
      myIndex(-1),
      startPos(0),
      endPos(0)
    {}

    //======================================================================

    static 
    const Whatever& null() {
      static Whatever result(WHATEVER_NULL);
      //      result.setNull();
      return result;
    }

    //======================================================================

    std::string sname() const {
      std::wstring wname(name());
      return std::string(wname.begin(),wname.end());
    }

    //======================================================================

    std::wstring name() const {
      std::wostringstream nam;
      collectName(nam);
      nam << "{" << typeName(type) << "}";
      return nam.str();
    }

    //======================================================================

    void collectName(std::wostringstream& nam) const {
      if(parent == 0) {
	nam << "root";
	return;
      }
      parent->collectName(nam);
      switch (parent->type) {
      case WHATEVER_MAP:
	nam << L"[" << myKey << L"]";
	break;
      case WHATEVER_VECTOR:
	nam << L"[" << myIndex << L"]";
	break;
      default:
	nam << L"[ ?? not vector or map! ]";
      }
    }

    //======================================================================

    void assertOkWhateverType(WhateverType t, std::string location) const {
      if (type != t) {
	std::ostringstream msg;
	std::wstring tName(typeName(t));
	std::wstring t2Name(typeName(type));
	msg << "Error in '" << location << "' '" << sname() << "'\n"
	    << " assertOkWhateverType(" << std::string(tName.begin(),tName.end()) << ") actual type is " <<  std::string(t2Name.begin(),t2Name.end())  << "\n";
	throw std::logic_error(msg.str());
      }
    }
      
    //======================================================================

    void assertWhateverMap(std::string moreMsg="") const {
      if (type != WHATEVER_MAP) {
	std::wstring t2Name(typeName(type));
	std::ostringstream msg;
	msg << sname();
	msg << " is the wrong type should be WHATEVER_MAP actual type is " <<  std::string(t2Name.begin(),t2Name.end())  << "\n";
	msg << moreMsg;
	throw std::logic_error(msg.str());
      }
    }

    //======================================================================
    //
    void flattenInto(FlatMapType& flatMap, const std::wstring prefix=L"") const {
      
      this->assertWhateverMap  ("In flattenInto");

      WhateverMap::const_iterator itr = whateverMap.begin();

      for (; itr != whateverMap.end(); itr++) {

	std::wostringstream newKey;
	if (prefix.length() == 0)
	  newKey << itr->first;
	else
	  newKey << prefix << L":" << itr->first;

	flatMap[newKey.str()] = &(itr->second); 

	//	std::wcout << newKey.str() << " " << flatMap[newKey.str()] << " : " << typeName(flatMap[newKey.str()]->type) << L"\n";
	
	if (itr->second.type == WHATEVER_MAP) 
	  itr->second.flattenInto(flatMap, newKey.str());
      }
    }
    
    //======================================================================

    void assertKey(std::wstring key) const {

      this->assertWhateverMap("In assertKey");

      if (whateverMap.count(key) !=  1) {
	std::ostringstream msg;
	msg << sname() << "[" << std::string(key.begin(),key.end()) << "] key not found!\n";
	throw std::logic_error(msg.str());
      }
    }
    
    //====================================================================== map[]
    
    Whatever& operator[] (const std::wstring key) {
      assertKey(key);
      return whateverMap[key];
    }
    
    //====================================================================== map[]
    
    const Whatever& operator[] (const std::wstring key) const {
      assertKey(key);
      WhateverMap& wm = const_cast<WhateverMap&>(whateverMap);
      return wm[key];
    }
    
    //====================================================================== map[]
    
    Whatever& operator[] (const std::string key) {
      std::wstring wKey(key.begin(),key.end());
      assertKey(wKey);
      return whateverMap[wKey];
    }
    
    //====================================================================== map[]
    
    const Whatever& operator[] (const std::string key) const {
      std::wstring wKey(key.begin(),key.end());
      assertKey(wKey);
      WhateverMap& wm = const_cast<WhateverMap&>(whateverMap);
      return wm[wKey];
    }
    
    //======================================================================
    
    size_t size() const {
      
      assertOkWhateverType(WHATEVER_VECTOR, "vector.size()");
      return whateverVector.size();
    }
      
    //====================================================================== count
    
    int count(std::string key) const {
      return count(std::wstring(key.begin(),key.end()));
    }
      
    //====================================================================== count
    
    int count(std::wstring key) const {
      
      std::ostringstream msg;
      msg << "count(" << std::string(key.begin(),key.end())  << ")";
      assertOkWhateverType(WHATEVER_MAP, msg.str());
      
      return whateverMap.count(key);
    }
      
    //====================================================================== Vector[]
    
    Whatever& operator[] (size_t index) {
      
      std::ostringstream msg;
      msg << "Whatever[" << index << "]";
      assertOkWhateverType(WHATEVER_VECTOR, msg.str());
      
      if (!(index <= whateverVector.size()) ) {
	std::ostringstream msg;
	msg << "There was no value for, Whatever[" << index << "] \n";
	throw std::range_error(msg.str());
      }

      return whateverVector[index];
    }

    //====================================================================== Vector[]
    
    const Whatever& operator[] (size_t index) const {
      
      std::ostringstream msg;
      msg << "Whatever[" << index << "]";
      assertOkWhateverType(WHATEVER_VECTOR, msg.str());
      
      if (!(index <= whateverVector.size()) ) {
	std::ostringstream msg;
	msg << "There was no value for, Whatever[" << index << "] \n";
	throw std::range_error(msg.str());
      }

      return whateverVector[index];
    }

    //====================================================================== Vector[]
    
    Whatever& back() {
      
      std::ostringstream msg;
      msg << "Whatever.back()";
      assertOkWhateverType(WHATEVER_VECTOR, msg.str());
      
      return whateverVector.back();
    }

    //====================================================================== Vector[]

    template<typename T>
    Whatever& push_back(T& value) {
      
      std::ostringstream msg;
      msg << "Whatever.push_back()";
      assertOkWhateverType(WHATEVER_VECTOR, msg.str());
      
      whateverVector.push_back(Whatever());
      whateverVector.back() = value;
      return *this;
    }

    Whatever& push_back() {
      
      std::ostringstream msg;
      msg << "Whatever.push_back()";
      assertOkWhateverType(WHATEVER_VECTOR, msg.str());
      
      whateverVector.push_back(Whatever());
      return *this;
    }

    Whatever& push_back_null() {
      
      std::ostringstream msg;
      msg << "Whatever.push_back_null()";
      assertOkWhateverType(WHATEVER_VECTOR, msg.str());
      
      whateverVector.push_back(Whatever());
      whateverVector.back().setNull();
      return *this;
    }

    //====================================================================== get

    template<typename T>
    void get (T& value) {
      std::istringstream valueStream(valueString);
      valueStream >> value;
    }

    //====================================================================== get

    void setNull () {
      type = WHATEVER_NULL;
      valueString = L"NULL";
    }
    
    //====================================================================== set

    //     template<typename T>
    //     void set (T& value) {
    //       type = TYPE<T>::to();
    //       std::wostringstream valueStream;
    //       valueStream << value;
    //       valueString = valueStream.str();
    //     }
    
    void set(const std::wstring& value) {
      type = WHATEVER_STRING;
      valueString = value;
    }

    void set(bool value) {
      type = WHATEVER_BOOL;
      whateverBool = value;
    }

    void set(int value) {
      type = WHATEVER_INTEGER;
      whateverInteger = value;
    }

    void set(double value) {
      type = WHATEVER_DOUBLE;
      whateverDouble = value;
    }

    //     void set(WhateverType t, std::wstring& value) {
    //       type = t;
    //       valueString = value;
    //     }

    //====================================================================== assign

    template<typename T>
    Whatever& operator = (const T& value) {
      set(value);
      return *this;
    }

    //======================================================================

  };

  template<typename T> class TYPE       { public: static Whatever::WhateverType to() {return Whatever::WHATEVER_UNKNOWN;  } };

  template<> class TYPE<std::string>    { public: static Whatever::WhateverType to() {return Whatever::WHATEVER_STRING;   } };
  template<> class TYPE<int>            { public: static Whatever::WhateverType to() {return Whatever::WHATEVER_INTEGER;  } };
  template<> class TYPE<double>         { public: static Whatever::WhateverType to() {return Whatever::WHATEVER_DOUBLE;   } };
  template<> class TYPE<bool>           { public: static Whatever::WhateverType to() {return Whatever::WHATEVER_BOOL;     } };
  template<> class TYPE<std::map<std::wstring,Whatever> >  { public: Whatever::WhateverType to() {return Whatever::WHATEVER_MAP;    } };
  template<> class TYPE<std::vector<Whatever> >            { public: Whatever::WhateverType to() {return Whatever::WHATEVER_VECTOR; } };
  
  //======================================================================
  
  int& operator <= (int& lhs, const Whatever& w) {
    switch (w.type) {
    case Whatever::WHATEVER_INTEGER: {
      lhs = w.whateverInteger;
      return lhs;
    }
    case Whatever::WHATEVER_MAT: 
    case Whatever::WHATEVER_VECTOR:
    case Whatever::WHATEVER_MATRIX:
    case Whatever::WHATEVER_STRING:
    case Whatever::WHATEVER_DOUBLE:
    case Whatever::WHATEVER_BOOL:
    case Whatever::WHATEVER_NULL:
    case Whatever::WHATEVER_UNKNOWN:
    default: {
      std::ostringstream msg;
      msg << "int d <= " << w.sname() << " produced a type error!\n";
      msg << " trying to assign a " << Whatever::ntypeName(w.type) << " to a int!\n";
      throw std::logic_error(msg.str());
    }
    }
  }
  
  template<typename T>
  std::complex<T>& operator <= (std::complex<T>& lhs, const Whatever& w) {
    
    switch (w.type) {
    case Whatever::WHATEVER_VECTOR: {
      T v0;
      v0 <= w.whateverVector[0];
      T v1;
      v1 <= w.whateverVector[1];
      std::complex<T> result(v0,v1);
      lhs = result;
      return lhs;
    }
    case Whatever::WHATEVER_MAT: 
    case Whatever::WHATEVER_MATRIX:
    case Whatever::WHATEVER_STRING:
    case Whatever::WHATEVER_INTEGER:
    case Whatever::WHATEVER_DOUBLE:
    case Whatever::WHATEVER_BOOL:
    case Whatever::WHATEVER_NULL:
    case Whatever::WHATEVER_UNKNOWN:
    default: {
      std::ostringstream msg;
      msg << "std::complex<T>  <= " << w.sname() << " produced a type error!\n";
      msg << " trying to assign a " << Whatever::ntypeName(w.type) << " to a std::complex<T>!\n";
      throw std::logic_error(msg.str());
    }
    }
  }

  std::string& operator <= (std::string& lhs, const Whatever& w) {
    switch (w.type) {
    case Whatever::WHATEVER_STRING:
      lhs = std::string(w.valueString.begin(), w.valueString.end());
      return lhs;
    case Whatever::WHATEVER_MAT: 
    case Whatever::WHATEVER_MATRIX:
    case Whatever::WHATEVER_INTEGER:
    case Whatever::WHATEVER_DOUBLE:
    case Whatever::WHATEVER_BOOL:
    case Whatever::WHATEVER_NULL:
    case Whatever::WHATEVER_UNKNOWN:
    default: {
      std::ostringstream msg;
      msg << "std::string d <= " << w.sname() << " produced a type error!\n";
      msg << " trying to assign a " << Whatever::ntypeName(w.type) << " to a std::string!\n";
      throw std::logic_error(msg.str());
    }
    }
  }

  float& operator <= (float& lhs, const Whatever& w) {
    switch (w.type) {
    case Whatever::WHATEVER_INTEGER:
      lhs = static_cast<float>(w.whateverInteger);
      return lhs;
    case Whatever::WHATEVER_DOUBLE:
      lhs = static_cast<float>(w.whateverDouble);
      return lhs;
    case Whatever::WHATEVER_MAT: 
    case Whatever::WHATEVER_MATRIX:
    case Whatever::WHATEVER_STRING:
    case Whatever::WHATEVER_BOOL:
    case Whatever::WHATEVER_NULL:
    case Whatever::WHATEVER_UNKNOWN:
    default: {
      std::ostringstream msg;
      msg << "float d <= " << w.sname() << " produced a type error!\n";
      msg << " trying to assign a " << Whatever::ntypeName(w.type) << " to a float!\n";
      throw std::logic_error(msg.str());
    }
    }
  }

  bool& operator <= (bool& lhs, const Whatever& w) {
    switch (w.type) {
    case Whatever::WHATEVER_BOOL:
      lhs = static_cast<bool>(w.whateverBool);
      return lhs;
    case Whatever::WHATEVER_INTEGER:
      lhs = static_cast<bool>(w.whateverInteger);
      return lhs;
    case Whatever::WHATEVER_MAT: 
    case Whatever::WHATEVER_MATRIX:
    case Whatever::WHATEVER_STRING:
    case Whatever::WHATEVER_DOUBLE:
    case Whatever::WHATEVER_NULL:
    case Whatever::WHATEVER_UNKNOWN:
    default: {
      std::ostringstream msg;
      msg << "float bool <= " << w.sname() << " produced a type error!\n";
      msg << " trying to assign a " << Whatever::ntypeName(w.type) << " to a bool!\n";
      throw std::logic_error(msg.str());
    }
    }
  }

  double& operator <= (double& lhs, const Whatever& w) {

    switch (w.type) {
    case Whatever::WHATEVER_INTEGER:
      lhs = static_cast<double>(w.whateverInteger);
      return lhs;
    case Whatever::WHATEVER_DOUBLE:
      lhs = w.whateverDouble;
      return lhs;;
    case Whatever::WHATEVER_MAT: 
    case Whatever::WHATEVER_MATRIX:
    case Whatever::WHATEVER_STRING:
    case Whatever::WHATEVER_BOOL:
    case Whatever::WHATEVER_NULL:
    case Whatever::WHATEVER_UNKNOWN:
    default: {
      std::ostringstream msg;
      msg << "double d <= " << w.sname() << " produced a type error!\n";
      msg << " trying to assign a " << Whatever::ntypeName(w.type) << " to a double!\n";
      throw std::logic_error(msg.str());
    }
    }
  }

  template<typename T>
  std::vector<T>& operator <= (std::vector<T>& lhs, const Whatever& w) {

    switch (w.type) {
    case Whatever::WHATEVER_VECTOR:
      lhs.resize(w.whateverVector.size(),0);
      for (size_t i=0; i< lhs.size(); i++)
	lhs[i] <= w.whateverVector[i];
      return lhs;
    case Whatever::WHATEVER_MAT: 
    case Whatever::WHATEVER_MATRIX:
    case Whatever::WHATEVER_STRING:
    case Whatever::WHATEVER_INTEGER:
    case Whatever::WHATEVER_DOUBLE:
    case Whatever::WHATEVER_BOOL:
    case Whatever::WHATEVER_NULL:
    case Whatever::WHATEVER_UNKNOWN:
    default: {
      std::ostringstream msg;
      msg << "std::vector<t> <= " << w.sname() << " produced a type error!\n";
      msg << " trying to assign a " << Whatever::ntypeName(w.type) << " to a std::vector<t>!\n";
      throw std::logic_error(msg.str());
    }
    }
  }

  template<typename T>
  std::map<std::string,T>& operator <= (std::map<std::string,T>& lhs, const Whatever& w) {

    switch (w.type) {
    case Whatever::WHATEVER_MAP: {
      
      typedef std::map<std::wstring,Whatever> WhateverMapType;
      lhs.clear();
    
      WhateverMapType::const_iterator itr = w.whateverMap.begin();
      for (; itr != w.whateverMap.end(); itr++) {
	//	const std::pair<std::wstring,Whatever>& pair(*itr);
	std::string key(itr->first.begin(), itr->first.end());
	lhs[key] <= itr->second;
      }
      return lhs;
    }
    case Whatever::WHATEVER_MAT: 
    case Whatever::WHATEVER_VECTOR:
    case Whatever::WHATEVER_MATRIX:
    case Whatever::WHATEVER_STRING:
    case Whatever::WHATEVER_INTEGER:
    case Whatever::WHATEVER_DOUBLE:
    case Whatever::WHATEVER_BOOL:
    case Whatever::WHATEVER_NULL:
    case Whatever::WHATEVER_UNKNOWN:
    default: {
      std::ostringstream msg;
      msg << "std::map<t,t2> <= " << w.sname() << " produced a type error!\n";
      msg << " trying to assign a " << Whatever::ntypeName(w.type) << " to a std::map!\n";
      throw std::logic_error(msg.str());
    }
    }
  }

  //======================================================================

  template<typename StreamType>
  StreamType& operator << (StreamType& os, const Whatever& w) {

    typedef Whatever::WhateverMap        WhateverMap;
    typedef WhateverMap::const_iterator  WhateverMapItr;
    typedef std::wstring                 StringType;

    switch (w.type) {
    case Whatever::WHATEVER_MAT: {
      StringType wfilename(w.filename.begin(),w.filename.end());
      os << "{ 'fileName': '" << wfilename << "'"
	 << ", 'startPos': "  << w.startPos 
	 << ", 'endPos': "    << w.endPos << "}";
      break;
    }
    case Whatever::WHATEVER_MAP:
      os << "{";
      for (WhateverMapItr itr = w.whateverMap.begin();
	   itr != w.whateverMap.end();
	   itr++) {
	const StringType& wkey((*itr).first);
	os << "\"" << wkey << "\": " << (*itr).second << ",";
      }
      os << "}";
      break;
    case Whatever::WHATEVER_VECTOR:
      os << "[";
      for (size_t i=0; i<w.whateverVector.size();i++)
	os << w.whateverVector[i] << ",";
      os << "]";
      break;
    case Whatever::WHATEVER_MATRIX:
      os << "WHATEVER_MATRIX";
      break;
    case Whatever::WHATEVER_STRING:
      os << "\"" << w.valueString << "\"";
      break;
    case Whatever::WHATEVER_INTEGER:
      os << w.whateverInteger;
      break;
    case Whatever::WHATEVER_DOUBLE:
      os << w.whateverDouble;
      break;
    case Whatever::WHATEVER_BOOL:
      os << w.whateverBool;
      break;
    case Whatever::WHATEVER_UNKNOWN:
      os <<"WHATEVER_UNKNOWN";
      break;
    default:
      throw std::logic_error("Whatever::typeName given wrong type");
    }
    return os;
  }
    
  //======================================================================
  
} // end namespace JsonParser


  /*@}*/
#endif

