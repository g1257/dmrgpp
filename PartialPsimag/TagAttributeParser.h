//-*-C++-*-

#ifndef TAG_Attribute_Parser_H
#define TAG_Attribute_Parser_H

/** \ingroup ostream **/
/*@{*/

/** \file TagAttributeParser.h
 *  Contains the class definition for XML tag attributes.
 */

#include <cstddef>
#include <limits>
#include <list>
#include <stdexcept>
#include <vector>
#include <ostream>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <functional>

#include "PSIMAGAssert.h"

namespace psimag {

  //======================================================================

  /**
   * \brief A class representing the attributes of a Tag.
   */

  using namespace std;
  
  class TagAttributeParser {
    
  private:
    
    string::size_type len;
    string::size_type start;
    string::size_type nextPos;
    
    const string     input;
    const string     delimiters;
    const string     ignore;
    const string     quotes;
    const string     slashes;

    string     token;
    
    map<string,string>& attributes;

  public:
    
    
    TagAttributeParser(const string& in, map<string,string>& attributes_):
      len(in.length()),
      start(0),
      input(in),
      delimiters("="),
      ignore(" 	\n"),
      quotes("'"),
      slashes("\\"),
      attributes(attributes_)
    {
      while (getKeyValuePair());
    }

  private:

    string::size_type remainingChars() {
      return len - start;
    }

    bool getWhiteSpace() {

      string::size_type wsEndPos = input.find_first_not_of(ignore, start);
      if (wsEndPos == string::npos) {
	start = len;
	return true;
      }
      if (start == wsEndPos) 
	return false;
      start = wsEndPos;
      return true;
    }

    bool getQuotedToken() {
      
      getWhiteSpace();
      
      if (remainingChars() < 2)  return false;
      
      if (input.find_first_of(quotes,start) != start) return false;

      char              qchar = input[start];
      string::size_type cursor = start+1;
      string::size_type nextQuotePos;

      while(true) {
	nextQuotePos = input.find_first_of(qchar,cursor);
	if (nextQuotePos == string::npos) return false;
	if (input[nextQuotePos-1] != '\\')  break;
	cursor = nextQuotePos+1;
      }

      token = input.substr(start+1, nextQuotePos-(start+1));
      start = nextQuotePos +1;
      return true;
    }

    bool getEqual() {
      
      getWhiteSpace();

      if (input.find_first_of('=',start) != start)
	return false;
      
      start = start + 1;
      return true;
    }

    bool getBasicToken() {
      
      getWhiteSpace();

      string::size_type pos = input.find_first_of(ignore+delimiters,start);
      if (pos == string::npos) {
	token = input.substr(start, remainingChars());
	start = len;
	return true;
      }

      token = input.substr(start, pos-start);
      start = pos;
      return true;
    }

    bool getQuotedOrBasicToken() {
      if(getQuotedToken()) return true;
      return getBasicToken();
    }

    bool getKeyValuePair() {
      
      if (remainingChars() == 0) return false;

      if (!getBasicToken()) return false;

      string key = token;

      if (!getEqual() || !getQuotedOrBasicToken()) return false;

      attributes[key] = token;
      return true;
    }

    bool ignoreTillSeperator() {

      if (remainingChars() == 0) return false;

      string::size_type delimPos = input.find_first_of(delimiters, start);
      if (delimPos == string::npos) {
	start = len;
	return true;
      }
      if (start == delimPos) return false;
      start = delimPos + 1;
      return true;
    }

  };

}  

#endif

/*@}*/

       
