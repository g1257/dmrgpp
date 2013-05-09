//-*-C++-*-
// Author: Michael S. Summers (ORNL)
/** \ingroup JsonParser */
/*@{*/

/*! \file ModesMixin.h  
 *
 *  
 *
 */

#ifndef  JsonParser_ModesMixin_H
#define  JsonParser_ModesMixin_H

#include <vector>
#include "String.h"
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace JsonParser {
  
  class ModesMixin {
  public:

    /* These modes can be pushed on the stack. */
    typedef enum modes {
      MODE_ARRAY = 1, 
      MODE_DONE = 2,  
      MODE_KEY = 3,   
      MODE_OBJECT = 4
    } ModeType;

    PsimagLite::Vector<ModeType>::Type stack;
    
    ModesMixin():
      stack(1,MODE_DONE)
    {}
    
    //======================================================================

    PsimagLite::String modeName(ModeType m) {
      switch (m) {
      case MODE_ARRAY:
	return "MODE_ARRAY";
      case MODE_DONE:
	return "MODE_DONE";
      case MODE_KEY:
	return "MODE_KEY";
      case MODE_OBJECT:
	return "MODE_OBJECT";
      default:
	return "MODE_UNKNOWN";
      }
    }

    //======================================================================

    void push(ModeType mode) {
      stack.push_back(mode);
    }
    
    void pop(ModeType expectedMode) {
      if (stack.size() == 0) {
	PsimagLite::OstringStream msg;
	msg << "JsonParser.pop was expecting mode " << modeName(expectedMode) << " to be on the back of the stack. \n"
	    << "However the stack was empty!\n";
	throw std::logic_error(msg.str());
      }
      if (expectedMode != stack.back()) {
	PsimagLite::OstringStream msg;
	msg << "JsonParser.pop was expecting mode " << modeName(expectedMode) << " to be on the back of the stack. \n"
	    << "However the back of the stack contained " << modeName(stack.back()) << "\n";
	throw std::logic_error(msg.str());
      }
      stack.pop_back();
    }
    
    const ModeType& currentMode() {
      return stack.back();
    }
  };

} // end namespace JsonParser


/*@}*/
#endif

