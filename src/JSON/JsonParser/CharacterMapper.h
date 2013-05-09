//-*-C++-*-
// Author: Michael S. Summers (ORNL)
/** \ingroup JsonParser */
/*@{*/

/*! \file CharacterMapper.h  
 *
 *  
 *
 */

#ifndef  JsonParser_CharacterMapper_H
#define  JsonParser_CharacterMapper_H

#include <wchar.h>
#include "String.h"
#include <stdexcept>

namespace JsonParser {

  class CharacterMapper {

  public:

    /*
      Characters are mapped into these 31 character classes. This allows for
      a significant reduction in the size of the state transition table.
    */
    typedef enum classes {
      C_SPACE,  /* space */
      C_WHITE,  /* other whitespace */
      C_LPARN,  /* (  */
      C_RPARN,  /* ) */
      C_LCURB,  /* {  */
      C_RCURB,  /* } */
      C_LSQRB,  /* [ */
      C_RSQRB,  /* ] */
      C_COLON,  /* : */
      C_COMMA,  /* , */
      C_QUOTE,  /* " */
      C_BACKS,  /* \ */
      C_SLASH,  /* / */
      C_PLUS,   /* + */
      C_MINUS,  /* - */
      C_POINT,  /* . */
      C_ZERO ,  /* 0 */
      C_DIGIT,  /* 123456789 */
      C_LOW_A,  /* a */
      C_LOW_B,  /* b */
      C_LOW_C,  /* c */
      C_LOW_D,  /* d */
      C_LOW_E,  /* e */
      C_LOW_F,  /* f */
      C_LOW_L,  /* l */
      C_LOW_N,  /* n */
      C_LOW_R,  /* r */
      C_LOW_S,  /* s */
      C_LOW_T,  /* t */
      C_LOW_U,  /* u */
      C_LOW_Y,  /* y */
      C_ABCDF,  /* ABCDF */
      C_E,      /* E */
      C_ETC,    /* everything else */
      C_STAR,   /* * */   
      C_EOF,    /* end of file */
      C_ERR,    /* error */
      NR_CLASSES
    } CharacterClass;

     static PsimagLite::String clsName(CharacterClass cls) {
      switch (cls) {
      case C_SPACE: return "C_SPACE space ";
      case C_WHITE: return "C_WHITE other whitespace ";
      case C_LPARN: return "C_LPARN ( ";
      case C_RPARN: return "C_RPARN ) ";
      case C_LCURB: return "C_LCURB { ";
      case C_RCURB: return "C_RCURB } ";
      case C_LSQRB: return "C_LSQRB [ ";
      case C_RSQRB: return "C_RSQRB ] ";
      case C_COLON: return "C_COLON : ";
      case C_COMMA: return "C_COMMA , ";
      case C_QUOTE: return "C_QUOTE \\\" ";
      case C_BACKS: return "C_BACKS \\ ";
      case C_SLASH: return "C_SLASH / ";
      case C_PLUS:  return "C_PLUS + ";
      case C_MINUS: return "C_MINUS - ";
      case C_POINT: return "C_POINT . ";
      case C_ZERO:  return "C_ZERO 0 ";
      case C_DIGIT: return "C_DIGIT 123456789 ";
      case C_LOW_A: return "C_LOW_A a ";
      case C_LOW_B: return "C_LOW_B b ";
      case C_LOW_C: return "C_LOW_C c ";
      case C_LOW_D: return "C_LOW_D d ";
      case C_LOW_E: return "C_LOW_E e ";
      case C_LOW_F: return "C_LOW_F f ";
      case C_LOW_L: return "C_LOW_L l ";
      case C_LOW_N: return "C_LOW_N n ";
      case C_LOW_R: return "C_LOW_R r ";
      case C_LOW_S: return "C_LOW_S s ";
      case C_LOW_T: return "C_LOW_T t ";
      case C_LOW_U: return "C_LOW_U u ";
      case C_LOW_Y: return "C_LOW_Y y ";
      case C_ABCDF: return "C_ABCDF ABCDF ";
      case C_E:     return "C_E E ";
      case C_ETC:   return "C_ETC everything else ";
      case C_STAR:  return "C_STAR * ";
      case C_EOF:   return "C_ERR error ";
      case C_ERR:   return "C_ERR error ";
      case NR_CLASSES: return "NP_CLASSES ";
      default:    return "Unkown CharacterClass ";
      }
  }

    static CharacterClass wcharToClass(wchar_t widec) {

      /* Determine the character's class. */
      if (widec >= 128) return C_ETC;

      static CharacterClass ascii_class[128] = {
	/*
	  This array maps the 128 ASCII characters into character classes.
	  The remaining Unicode characters should be mapped to C_ETC.
	  Non-whitespace control characters are errors.
	*/
	C_ERR,     C_ERR,     C_ERR,     C_ERR,     C_ERR,     C_ERR,     C_ERR,    C_ERR,
	C_ERR,     C_WHITE, C_WHITE,     C_ERR,     C_ERR,     C_WHITE,   C_ERR,    C_ERR,
	C_ERR,     C_ERR,     C_ERR,     C_ERR,     C_ERR,     C_ERR,     C_ERR,    C_ERR,
	C_ERR,     C_ERR,     C_ERR,     C_ERR,     C_ERR,     C_ERR,     C_ERR,    C_ERR,
	
	C_SPACE, C_ETC,   C_QUOTE, C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,
	C_LPARN, C_RPARN, C_STAR,   C_PLUS,  C_COMMA, C_MINUS, C_POINT, C_SLASH,
	C_ZERO,  C_DIGIT, C_DIGIT, C_DIGIT, C_DIGIT, C_DIGIT, C_DIGIT, C_DIGIT,
	C_DIGIT, C_DIGIT, C_COLON, C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,
	
	C_ETC,   C_ABCDF, C_ABCDF, C_ABCDF, C_ABCDF, C_E,     C_ABCDF, C_ETC,
	C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,
	C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,
	C_ETC,   C_ETC,   C_ETC,   C_LSQRB, C_BACKS, C_RSQRB, C_ETC,   C_ETC,
	
	C_ETC,   C_LOW_A, C_LOW_B, C_LOW_C, C_LOW_D, C_LOW_E, C_LOW_F, C_ETC,
	C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_LOW_L, C_ETC,   C_LOW_N, C_ETC,
	C_ETC,   C_ETC,   C_LOW_R, C_LOW_S, C_LOW_T, C_LOW_U, C_ETC,   C_ETC,
	C_ETC,   C_LOW_Y, C_ETC,   C_LCURB, C_ETC,   C_RCURB, C_ETC,   C_ETC
      };
      
      CharacterClass& result = ascii_class[wctob(widec)];

      if (result == C_ERR) {
	std::wostringstream msg;
	msg << L"CharacterMapper::wcharToClass(" << widec << L") resulted in an error!\n";
	std::wstring m(msg.str());
	throw std::logic_error(PsimagLite::String(m.begin(), m.end()));
      }
      
      return result;
    }

    static bool isWhiteSpace(CharacterClass nextClass) {
      if (nextClass == C_SPACE) return true;
      if (nextClass == C_WHITE) return true;
      return false;
    }
    
    static bool foo() { return false;}

    static bool is_legal_white_space(const std::wstring& s) {

      CharacterClass char_class;
      wchar_t c;
      
      if (s.size() == 0) 
	return false;
      
      for (size_t i=0; i< s.size(); i++) {
	c = s[i];
	if (c < 0 || c >= 128) 
	  return false;
      }
     
      char_class = wcharToClass(c);
      
      if (char_class != C_SPACE && 
	  char_class != C_WHITE ) 
	return false;

      return true;
    }

    //======================================================================
    
    static wchar_t getEscapedCharacter(std::wistream& inputStream) {
	
      wchar_t secondChar = inputStream.get();
      
      switch(secondChar) {
      case L'b':
	return L'\b';
      case L'f':
	return L'\f';
      case L'n':
	return L'\n';
      case L'r':
	return L'\r';
      case L't':
	return L'\t';
      case L'"':
	return L'"';
      case L'\\':
	return L'\\';
      case L'/':
	return L'/';
      // 	case 'u':
      // 	  buffer.push_back( '\\');
      // 	  buffer.push_back( 'u');
      default:
	throw std::logic_error("JsonParser: Encountered an escapped character that was not recognized.");
      }
    }
      
    //======================================================================

    static void skipWhiteAnd(std::wstring chars, std::wistream& inputStream) {
      
      while(true) {

	wchar_t        nextChar;
	nextChar = inputStream.get();
	if (chars.find_first_of(nextChar) != std::wstring::npos  ||
	    isWhiteSpace(wcharToClass(nextChar)) )
	  continue;
	
	inputStream.unget();
	return;
      }
    }

    //======================================================================
    
  };

} // end namespace JsonParser


/*@}*/
#endif
