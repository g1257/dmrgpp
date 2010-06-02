//-*-C++-*-

#ifndef TAG_H
#define TAG_H

/** \ingroup ostream **/
/*@{*/

/** \file Tag.h
 *  Contains the class definition for XML Tag and related objects.
 */

#include <cstddef>
#include <limits>
#include <list>
#include <stdexcept>
#include <vector>
#include <ostream>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <algorithm>
#include <functional>

#include "PSIMAGAssert.h"
#include "TagAttributes.h"

namespace psimag {

  //======================================================================

  class Tag: public TagAttributes {
  public:

    class Elements: public std::vector<Tag> {
    public:
      typedef std::vector<Tag> SuperType;
      
      template<typename TagCollectionType>
      void add(const TagCollectionType& tags) {
	this->insert(this->end(),tags.begin(),tags.end());
      }
      
      void add(const Tag& tag) {
	this->push_back(tag);
      }
    };
    
    
    std::string         name;
    Elements            elements;
    std::ostringstream  content;
    
    //------------------------------------------------------------ Default Constructor

    Tag():                                    name("span")                       {}

    //------------------------------------------------------------ One argument constructors

    Tag(const char* name_):                         name(name_)                         {}
    Tag(const std::string name_):                   name(name_)                         {}

    Tag(const Tag& tag): 
      TagAttributes(tag), 
      name(tag.name), 
      elements(tag.elements) 
    { 
      content.str(tag.content.str());  
    }

    //------------------------------------------------------------ Two argument constructors

    Tag(char* name_, const Elements&      elements_  ): name(name_), elements(elements_)        {}
    Tag(char* name_, const TagAttributes& attributes ): TagAttributes(attributes), name(name_)  {}
    Tag(char* name_, const Tag&        tag        ):    name(name_)                             { elements.add(tag); }
    template<typename ContentType>
    Tag(char* name_, const ContentType& cont      ):    name(name_)                             { content << " " << cont << " "; }

    //------------------------------------------------------------ Three argument constructors

    Tag(char* name_, const TagAttributes& attributes, const Elements& elements_): 
      TagAttributes(attributes), name(name_), elements(elements_) {}

    //------------------------------------------------------------

  /**
   * \brief Tag Assignment. 
   */
    Tag& operator= (const Tag& tag) {
      name                            = tag.name;
      elements                        = tag.elements;
      static_cast<TagAttributes>(*this)  = static_cast<TagAttributes>(tag);
      content.str(tag.content.str());
      return (*this);
    }

  /**
   * \brief Tag content Assignment. 
   */
    template<typename Type>
    Tag& operator << (const Type& something) {
      content << something;
      return *this;
    }

  /**
   * \brief Tag content Assignment. 
   */
    template<typename Type>
    Tag& operator << (const std::vector<Type>& something) {
      for(size_t i=0; i<something.size(); i++) 
	content << something[i] << " ";
      content << std::endl;
      return *this;
    }

  /**
   * \brief Add a collection of tags to this tags elements.
   */
    template< template<typename> class TagCollectionType >
    void add(const TagCollectionType<Tag>& tags) {
      elements.insert(elements.end(),tags.begin(),tags.end());
    }

  /**
   * \brief Add a tag to this tags elements.
   */
    void add(const Tag& tag) {
      elements.push_back(tag);
    }
  };

  //======================================================================
/* ATTENTION/FIXME: Header files should not contain implementations
        (because it causes multiple definitions when using more than one compilation unit)
        quick fix: I added the keyword inline (Gonzalo) */

  inline std::ostream& operator << (std::ostream& os, const Tag& t) {
    
    typedef TagAttributes::const_iterator AttributeItr;

    static const size_t contentBreakPoint = 100;
    std::string content = t.content.str();

    os << "<" 
       << t.name 
       << static_cast<TagAttributes>(t) ;

    if (t.elements.size() == 0 && content.size() == 0 ) {
      os << "/>" ;
      return os;
    }

    os << ">" ;

    if (t.elements.size() > 0 && content.size() > contentBreakPoint) 
      os << std::endl;

    for (size_t i=0; i < t.elements.size(); i++)
      os << std::endl << t.elements[i] ;

    if (t.elements.size() > 0)
      os << std::endl;

    os << " " << content << " ";
    
    if (content.size() >  contentBreakPoint) 
      os << std::endl;

    os << "</" << t.name << ">" ;
    return os;
  }

  //======================================================================

  template<typename Field>
  std::string toString(Field value) {
    std::ostringstream buff;
    buff << value ;
    return buff.str();
  }
  
  
  //======================================================================

}  

#endif

/*@}*/
