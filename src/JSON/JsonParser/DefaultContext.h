//-*-C++-*-
// Author: Michael S. Summers (ORNL)
/** \ingroup JsonParser */
/*@{*/

/*! \file DefaultContext.h  
 *
 *  
 *
 */

#ifndef  JsonParser_DefaultContext_H
#define  JsonParser_DefaultContext_H

#include "Whatever.h"
#include <iostream>
#include <sstream>
#include "assert.h"

namespace JsonParser {
  
  //======================================================================

  class DefaultContext {
  public:

    //======================================================================

    Whatever                  result;
    std::vector<Whatever*>    stack;
    std::wstring              key;
    bool                      trace;

    //======================================================================
    
    DefaultContext():
      result(),
      stack(1,&result),
      key(L""),
      trace(false)
    {}

    //======================================================================

    Whatever& currentObject() {
      return *stack.back();
    }
    
    //======================================================================

    //* Object and Array Referenced objects are created when needed on LHS of assignment! */
    Whatever& referencedObject() {
      
      Whatever&   current(currentObject());
      switch (current.type) {
      case Whatever::WHATEVER_MAP: {
	Whatever& result(current.whateverMap[key]);
	result.myKey  = key;
	result.parent = &current;
	if (trace) std::wcout << L"   referencedMap => '" <<  result.name() << L"' !\n";
	return result;
      }
      case Whatever::WHATEVER_VECTOR:  {
	size_t idx = current.size();
	current.push_back();
	Whatever& result(current.back());
	result.myIndex = idx;
	result.parent = &current;
	if (trace) std::wcout << L"   referencedVector => '" <<  result.name() << L"' !\n";
	return result;
      }
      default:
	if (trace) std::wcout << L"   referencedObject => '" <<  current.name() << L"' !\n";
	return current;
      }
    }
    
    //======================================================================

    template<int TypeId>
    void endObjectOrArray() {
      if (trace) std::wcout << L"   DefaultContext is ending object '" 
			    <<  currentObject().name() 
			    << L"' by popping the object stack!\n";
      if (stack.size() != 0) {
	assert(currentObject().type == TypeId);
	stack.pop_back();
	if (trace) {
	  if(stack.size() > 0) {
	    std::wcout << L"   current object is now '" 
		       <<  currentObject().name() <<  L"\n";
	  }
	  else
	    std::wcout << L"   Parsing completed! \n "; 
	}
      }
      else {
	std::wcout << L"   DefaultContext is ending an object that does not exist!."
		   << L"   Is there an extra '}'?\n"
		   << L"   The last result was : " << result.name() << L"\n";
      }
    }
    
    //======================================================================

    template<typename T>
    void set(const T& v) {
      Whatever& refObj(referencedObject());
      refObj = v;
      if (trace) std::wcout << L"  " <<  refObj.name() << L" = " << v << "\n";
    }

    //======================================================================

    template<Whatever::WhateverType TypeId>
    void beginObjectOrArray() {

      Whatever& refObj = referencedObject(); // Generally creates the object or array 
                                             // (except for the first time when refObject == result )
      refObj.type = TypeId;                  // Sets the type
      if (&refObj != &result)    	     // result is already on the stack, 
                                             // so we dont' need to push it on.
	stack.push_back(&refObj);            // Make the referenced object the current object
      
      if (trace) std::wcout << L" Set the type of " 
			    << refObj.name() 
			    << L" to " 
			    << Whatever::typeName(TypeId) 
			    << " and make it the current object\n";
    }

    //======================================================================
    
    std::string CurrentContext() {
      std::wstring result(referencedObject().name());
      return std::string(result.begin(),result.end());
    }
    
    //======================================================================
    
    void End(size_t numChars, size_t numLines) {
      std::cout << "Parsing completed! read " << numChars << " characters and " << numLines << " lines.\n";
    }
    
    //======================================================================
    
    void ObjectEnd() {
      endObjectOrArray<Whatever::WHATEVER_MAP>();
    }
    
    //======================================================================
    
    void MatrixEnd(size_t charNum) {
      currentObject().endPos = charNum;
      endObjectOrArray<Whatever::WHATEVER_MAT>();
    }
    
    //======================================================================
    
    void ArrayEnd() {
      endObjectOrArray<Whatever::WHATEVER_VECTOR>();
    }
    
    //======================================================================
    
    void ArrayBegin() {
      beginObjectOrArray<Whatever::WHATEVER_VECTOR>();
    }

    //======================================================================
    
    void ObjectBegin() {
      beginObjectOrArray<Whatever::WHATEVER_MAP>();
    }

    //======================================================================
    
    void MatrixBegin(std::string filename, size_t charNum) {
      beginObjectOrArray<Whatever::WHATEVER_MAT>();
      currentObject().filename = filename;
      currentObject().startPos = charNum;
    }

    //======================================================================
    
    void Integer(ParseBuffer& s) {
      int i;
      s >> i;
      set(i);
    }
    
    //======================================================================
    
    void Float(ParseBuffer& s) {
      double d;
      s >> d;
      set(d);
    }

    //======================================================================
    
    void Null( ) {
      set(Whatever::null());
    }
    
    //======================================================================
    
    void True( ) {
      set(true);
    }
    
    //======================================================================
    
    void False() {
      set(false);
    }
    
    //======================================================================
    
    void String(ParseBuffer& s) {
      set(s.str());
    }
    
    //======================================================================
    
    void Key(const std::wstring& s) {
      key = s;
      if (trace) std::wcout << L"   key = '" <<  key << L"'\n";
    }
    
    //======================================================================
    
    void Max(std::string s) {
      throw std::logic_error("Max! What is this?");
    }
    
    //======================================================================
    
    void Comment(std::wstring s) {
    }
    
    //======================================================================
    
    void None() {
      throw std::logic_error("None! What is this?");
    }

  
  };

  //======================================================================

  std::wostream& operator << (std::wostream& os, const DefaultContext& ctx) {

    os << ctx.result;

    return os;
  }

} // end namespace JsonParser


/*@}*/
#endif
