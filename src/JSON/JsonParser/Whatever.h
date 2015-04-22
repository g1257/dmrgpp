//-*-C++-*-
// Author: Michael S. Summers (ORNL)
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
#include "Map.h"
#include "Vector.h"
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

    typedef PsimagLite::Map<std::wstring,Whatever>::Type WhateverMap;
    typedef PsimagLite::Vector<Whatever>::Type           WhateverVector;

    typedef PsimagLite::Map<std::wstring,const JsonParser::Whatever*>::Type FlatMapType;


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
	throw PsimagLite::LogicError("Whatever::typeName given wrong type");
      }
    }

    //======================================================================

    static PsimagLite::String ntypeName(WhateverType t) {
      std::wstring wname(typeName(t));
      return PsimagLite::String(wname.begin(), wname.end());
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
    PsimagLite::String    filename;
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

    PsimagLite::String sname() const {
      std::wstring wname(name());
      return PsimagLite::String(wname.begin(),wname.end());
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

    void assertOkWhateverType(WhateverType t, PsimagLite::String location) const {
      if (type != t) {
	PsimagLite::OstringStream msg;
	std::wstring tName(typeName(t));
	std::wstring t2Name(typeName(type));
	msg << "Error in '" << location << "' '" << sname() << "'\n"
	    << " assertOkWhateverType(" << PsimagLite::String(tName.begin(),tName.end()) << ") actual type is " <<  PsimagLite::String(t2Name.begin(),t2Name.end())  << "\n";
	throw PsimagLite::LogicError(msg.str());
      }
    }

    //======================================================================

    void assertWhateverMap(PsimagLite::String moreMsg="") const {
      if (type != WHATEVER_MAP) {
	std::wstring t2Name(typeName(type));
	PsimagLite::OstringStream msg;
	msg << sname();
	msg << " is the wrong type should be WHATEVER_MAP actual type is " <<  PsimagLite::String(t2Name.begin(),t2Name.end())  << "\n";
	msg << moreMsg;
	throw PsimagLite::LogicError(msg.str());
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
	PsimagLite::OstringStream msg;
	msg << sname() << "[" << PsimagLite::String(key.begin(),key.end()) << "] key not found!\n";
	throw PsimagLite::LogicError(msg.str());
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

    Whatever& operator[] (const PsimagLite::String key) {
      std::wstring wKey(key.begin(),key.end());
      assertKey(wKey);
      return whateverMap[wKey];
    }

    //====================================================================== map[]

    const Whatever& operator[] (const PsimagLite::String key) const {
      std::wstring wKey(key.begin(),key.end());
      assertKey(wKey);
      WhateverMap& wm = const_cast<WhateverMap&>(whateverMap);
      return wm[wKey];
    }

    //======================================================================

    SizeType size() const {

      assertOkWhateverType(WHATEVER_VECTOR, "vector.size()");
      return whateverVector.size();
    }

    //====================================================================== count

    int count(PsimagLite::String key) const {
      return count(std::wstring(key.begin(),key.end()));
    }

    //====================================================================== count

    int count(std::wstring key) const {

      PsimagLite::OstringStream msg;
      msg << "count(" << PsimagLite::String(key.begin(),key.end())  << ")";
      assertOkWhateverType(WHATEVER_MAP, msg.str());

      return whateverMap.count(key);
    }

    //====================================================================== Vector[]

    Whatever& operator[] (SizeType index) {

      PsimagLite::OstringStream msg;
      msg << "Whatever[" << index << "]";
      assertOkWhateverType(WHATEVER_VECTOR, msg.str());

      if (!(index <= whateverVector.size()) ) {
	PsimagLite::OstringStream msg;
	msg << "There was no value for, Whatever[" << index << "] \n";
	throw PsimagLite::RangeError(msg.str());
      }

      return whateverVector[index];
    }

    //====================================================================== Vector[]

    const Whatever& operator[] (SizeType index) const {

      PsimagLite::OstringStream msg;
      msg << "Whatever[" << index << "]";
      assertOkWhateverType(WHATEVER_VECTOR, msg.str());

      if (!(index <= whateverVector.size()) ) {
	PsimagLite::OstringStream msg;
	msg << "There was no value for, Whatever[" << index << "] \n";
	throw PsimagLite::RangeError(msg.str());
      }

      return whateverVector[index];
    }

    //====================================================================== Vector[]

    Whatever& back() {

      PsimagLite::OstringStream msg;
      msg << "Whatever.back()";
      assertOkWhateverType(WHATEVER_VECTOR, msg.str());

      return whateverVector.back();
    }

    //====================================================================== Vector[]

    template<typename T>
    Whatever& push_back(T& value) {

      PsimagLite::OstringStream msg;
      msg << "Whatever.push_back()";
      assertOkWhateverType(WHATEVER_VECTOR, msg.str());

      whateverVector.push_back(Whatever());
      whateverVector.back() = value;
      return *this;
    }

    Whatever& push_back() {

      PsimagLite::OstringStream msg;
      msg << "Whatever.push_back()";
      assertOkWhateverType(WHATEVER_VECTOR, msg.str());

      whateverVector.push_back(Whatever());
      return *this;
    }

    Whatever& push_back_null() {

      PsimagLite::OstringStream msg;
      msg << "Whatever.push_back_null()";
      assertOkWhateverType(WHATEVER_VECTOR, msg.str());

      whateverVector.push_back(Whatever());
      whateverVector.back().setNull();
      return *this;
    }

    //====================================================================== get

    template<typename T>
    void get (T& value) {
      PsimagLite::IstringStream valueStream(valueString);
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

//     operator bool()         {  bool        result <= *this; return result;  }
//     operator int()          {  int         result <= *this; return result;  }
//     operator double()       {  double      result <= *this; return result;  }
//     operator float()        {  double      result <= *this; return static_cast<float>(result);  }
//     operator String()  {  String result <= *this; return result;  }

  };

  template<typename T> class TYPE       { public: static Whatever::WhateverType to() {return Whatever::WHATEVER_UNKNOWN;  } };

  template<> class TYPE<PsimagLite::String>    { public: static Whatever::WhateverType to() {return Whatever::WHATEVER_STRING;   } };
  template<> class TYPE<int>            { public: static Whatever::WhateverType to() {return Whatever::WHATEVER_INTEGER;  } };
  template<> class TYPE<double>         { public: static Whatever::WhateverType to() {return Whatever::WHATEVER_DOUBLE;   } };
  template<> class TYPE<bool>           { public: static Whatever::WhateverType to() {return Whatever::WHATEVER_BOOL;     } };
  template<> class TYPE<PsimagLite::Map<std::wstring,Whatever>::Type >  { public: Whatever::WhateverType to() {return Whatever::WHATEVER_MAP;    } };
  template<> class TYPE<PsimagLite::Vector<Whatever>::Type >            { public: Whatever::WhateverType to() {return Whatever::WHATEVER_VECTOR; } };

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
      PsimagLite::OstringStream msg;
      msg << "int d <= " << w.sname() << " produced a type error!\n";
      msg << " trying to assign a " << Whatever::ntypeName(w.type) << " to a int!\n";
      throw PsimagLite::LogicError(msg.str());
    }
    }
  }

  SizeType& operator <=(SizeType& lhs, const Whatever& w) {
	  int x = 0;
	  x <= w;
	  if (x<0) {
		  PsimagLite::OstringStream msg;
		  msg << "Expecting SizeType got negative int: "<<x<<"\n";
		  throw PsimagLite::LogicError(msg.str());
	  }
	  lhs = x;
	  return lhs;
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
      PsimagLite::OstringStream msg;
      msg << "std::complex<T>  <= " << w.sname() << " produced a type error!\n";
      msg << " trying to assign a " << Whatever::ntypeName(w.type) << " to a std::complex<T>!\n";
      throw PsimagLite::LogicError(msg.str());
    }
    }
  }

  PsimagLite::String& operator <= (PsimagLite::String& lhs, const Whatever& w) {
    switch (w.type) {
    case Whatever::WHATEVER_STRING:
      lhs = PsimagLite::String(w.valueString.begin(), w.valueString.end());
      return lhs;
    case Whatever::WHATEVER_MAT:
    case Whatever::WHATEVER_MATRIX:
    case Whatever::WHATEVER_INTEGER:
    case Whatever::WHATEVER_DOUBLE:
    case Whatever::WHATEVER_BOOL:
    case Whatever::WHATEVER_NULL:
    case Whatever::WHATEVER_UNKNOWN:
    default: {
      PsimagLite::OstringStream msg;
      msg << "String d <= " << w.sname() << " produced a type error!\n";
      msg << " trying to assign a " << Whatever::ntypeName(w.type) << " to a String!\n";
      throw PsimagLite::LogicError(msg.str());
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
      PsimagLite::OstringStream msg;
      msg << "float d <= " << w.sname() << " produced a type error!\n";
      msg << " trying to assign a " << Whatever::ntypeName(w.type) << " to a float!\n";
      throw PsimagLite::LogicError(msg.str());
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
      PsimagLite::OstringStream msg;
      msg << "float bool <= " << w.sname() << " produced a type error!\n";
      msg << " trying to assign a " << Whatever::ntypeName(w.type) << " to a bool!\n";
      throw PsimagLite::LogicError(msg.str());
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
      PsimagLite::OstringStream msg;
      msg << "double d <= " << w.sname() << " produced a type error!\n";
      msg << " trying to assign a " << Whatever::ntypeName(w.type) << " to a double!\n";
      throw PsimagLite::LogicError(msg.str());
    }
    }
  }

  template<typename T>
  typename PsimagLite::Vector<T>::Type& operator <= (typename PsimagLite::Vector<T>::Type& lhs, const Whatever& w) {

    switch (w.type) {
    case Whatever::WHATEVER_VECTOR:
      lhs.resize(w.whateverVector.size(),0);
      for (SizeType i=0; i< lhs.size(); i++)
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
      PsimagLite::OstringStream msg;
      msg << "typename Vector<t>::Type <= " << w.sname() << " produced a type error!\n";
      msg << " trying to assign a " << Whatever::ntypeName(w.type) << " to a std::vector<t>!\n";
      throw PsimagLite::LogicError(msg.str());
    }
    }
  }

  template<typename T>
  typename PsimagLite::Map<PsimagLite::String,T>::Type& operator <= (typename PsimagLite::Map<PsimagLite::String,T>::Type& lhs, const Whatever& w) {

    switch (w.type) {
    case Whatever::WHATEVER_MAP: {

      typedef PsimagLite::Map<std::wstring,Whatever>::Type WhateverMapType;
      lhs.clear();

      WhateverMapType::const_iterator itr = w.whateverMap.begin();
      for (; itr != w.whateverMap.end(); itr++) {
	//	const std::pair<std::wstring,Whatever>& pair(*itr);
	PsimagLite::String key(itr->first.begin(), itr->first.end());
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
      PsimagLite::OstringStream msg;
      msg << "Map<t,t2>::Type <= " << w.sname() << " produced a type error!\n";
      msg << " trying to assign a " << Whatever::ntypeName(w.type) << " to a std::map!\n";
      throw PsimagLite::LogicError(msg.str());
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
      for (SizeType i=0; i<w.whateverVector.size();i++)
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
      throw PsimagLite::LogicError("Whatever::typeName given wrong type");
    }
    return os;
  }

  //======================================================================

} // end namespace JsonParser


  /*@}*/
#endif

