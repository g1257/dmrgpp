//-*-C++-*-
// Author: Michael S. Summers (ORNL)
/** \ingroup JsonParser */
/*@{*/

/*! \file ParserBuffer.h  
 *
 *  
 *
 */

#ifndef  JsonParser_ParseBuffer_H
#define  JsonParser_ParseBuffer_H

#include <iostream>
#include <sstream>
#include <iomanip>
#include "TypesMixin.h"
#include <stdexcept>

namespace JsonParser {
  
  class ParseBuffer: 
    public TypesMixin
  //    public std::wstringstream
  {

  public:

    //    std::wstringstream& superStream;
    std::vector<wchar_t> theCharacters;
    bool                 trace;

    //======================================================================

    ParseBuffer():
      TypesMixin(),
      theCharacters(),
      trace(false)
    {}

    void clear() {
      theCharacters.clear();
      if (trace) std::wcout << L"   ParseBuffer: clear()\n";
    }
    
    void put(wchar_t wc) {
      theCharacters.push_back(wc);
    }

//     template<typename ReturnType>
//     ReturnType get() {
//       ReturnType value;
//       superStream >> value;
//       return value;
//     }

    std::wstring str() {
      return std::wstring(theCharacters.begin(), theCharacters.end());
    }
    
  };

  template<typename T>
  ParseBuffer& operator >> (ParseBuffer& buffer, T& value) {
    std::wistringstream is(buffer.str());
    is >> value;
    return buffer;
  }


} // end namespace JsonParser


/*@}*/
#endif
