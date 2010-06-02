//-*-C++-*-

#ifndef TAG_Attributes_H
#define TAG_Attributes_H

/** \ingroup ostream **/
/*@{*/

/** \file TagAttributes.h
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
#include "TagAttributeParser.h"

namespace psimag {

  //======================================================================

  class TagAttributes;

  /**
   * \brief A class representing the attributes of a Tag.
   */
  class TagAttributesClosure {
  public:
    std::map<std::string, std::string>& tagAttributes;
    std::string                         key;

    typedef std::string ValueType;

    TagAttributesClosure( const std::string                    _key,
			  std::map<std::string, std::string>&  _tagAttributes):
      tagAttributes(_tagAttributes),
      key(_key)
    {}
    
    template<typename InType>
    TagAttributesClosure& operator = (const InType& obj) {
      std::ostringstream objBuffer;
      objBuffer << obj;
      std::string objStr = objBuffer.str();
      tagAttributes[key] = objStr;
      return (*this);
    }
    
    template<typename InType1, typename InType2>
    TagAttributesClosure& operator = (const std::pair<InType1,InType2>& obj) {
      std::ostringstream objBuffer;
      objBuffer << "(" << obj.first << "," << obj.second << ")" << std::endl;
      std::string objStr = objBuffer.str();
      tagAttributes[key] = objStr;

      return (*this);
    }
    
    template<typename InType>
    TagAttributesClosure& operator = (const std::vector<InType>& vec) {
      std::ostringstream objBuffer;
      for (size_t i=0; i< vec.size(); i++)
	objBuffer << vec[i] << " ";
      std::string objStr = objBuffer.str();
      tagAttributes[key] = objStr;

      return (*this);
    }

    /** Provide a conversion operator to make this object look like a string.  */
    operator ValueType() const {			
      return tagAttributes[key];
    } 

  };
  

  /**
   * \brief A class representing the attributes of a Tag.
   */
  class TagAttributes  {
  public:

    typedef std::map<std::string, std::string>::const_iterator const_iterator;

    std::map<std::string, std::string> attrs;

    size_t size() const { return attrs.size(); }

    const_iterator begin() const { return attrs.begin(); }
    const_iterator end()   const { return attrs.end();   }

    void update(const std::string attrString) {
      TagAttributeParser(attrString, attrs);
    }

    TagAttributesClosure operator [] (std::string key) {
      TagAttributesClosure closure(key,attrs);
      return closure;
    }
    
  };
  
  //======================================================================
/* ATTENTION/FIXME: Header files should not contain implementations
        (because it causes multiple definitions when using more than one compilation unit)
        quick fix: I added the keyword inline (Gonzalo) */

  /**
   * \brief Attribute writer.
   */
  inline std::ostream& operator << (std::ostream& os, const TagAttributes& attr) {

    typedef TagAttributes::const_iterator AttributeItr;

    static const size_t maxAttr = 4;

    if (attr.size() <= maxAttr) {
      AttributeItr        a;
      for (a=attr.begin(); a!=attr.end(); a++) 
	os << " " << a->first << "=\"" << a->second << "\" ";
      return os;
    }
    
    bool         first = true;
    AttributeItr        a;
    for (a=attr.begin(); a!=attr.end(); a++) {
      if (first)
	first = false;
      else
	os << std::endl;
      os << "      " << a->first << "=\"" << a->second << "\" ";
    }
    return os;
  }

  //======================================================================

}  

#endif

/*@}*/
