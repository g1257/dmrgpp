//-*-C++-*-
// Author: Michael S. Summers (ORNL)
/** \ingroup JsonParser */
/*@{*/

/*! \file JsonParser.h  
 *
 *  
 *
 */

#ifndef  JsonParser_H
#define  JsonParser_H

#include "wchar.h"
#include "StatesMixin.h"
#include "ActionsMixin.h"
#include "ModesMixin.h"
#include "TypesMixin.h"
#include "CharacterMapper.h"
#include "ParseBuffer.h"
#include "DefaultContext.h"
#include "AugmentedStateTransitionTable.h"

namespace JsonParser {
  
  template<typename ContextType=DefaultContext>
  class JsonParser: 
    public CharacterMapper,
    public StatesMixin,
    public ActionsMixin,
    public ModesMixin,
    public TypesMixin
  {
  public:

    //======================================================================

    ContextType           ctx;
    ParseBuffer           buffer;
    StateType             state;
    JsonType              lastItemType;    

    //======================================================================

    bool                 escaped; 

    bool                 allow_comments; 
    StateType            before_comment_state;
    bool                 comment; 
    size_t               comment_begin_offset;

    //======================================================================

    size_t               numChar;
    size_t               numLines;

    bool                 trace;

    std::string          filename;

    //======================================================================

    JsonParser():
      ctx                  (),
      buffer               (),
      state                (GO),
      lastItemType         (JSON_T_NONE),
      escaped              (false),
      allow_comments       (true),
      before_comment_state (GO),
      comment              (false),
      comment_begin_offset (0),
      numChar              (0),
      numLines             (0),
      trace                (false),
      filename             ("")
    {}

    //======================================================================
    
    std::pair<wchar_t,CharacterClass> getNextCharacterAndClass(std::wistream& inputStream) {
      
      std::pair<wchar_t,CharacterClass> result;
      CharacterClass& nextClass(result.second);
      wchar_t&        nextChar (result.first);
      
      nextClass = C_WHITE;
      
      // Arrange to skip white space!
      while (isWhiteSpace(nextClass)) { 
	
	if (inputStream.eof()) {
	  nextChar = EOF;
	  nextClass = C_EOF;
	  return result;
	}

	nextChar = inputStream.get();

	if (nextChar == static_cast<wchar_t>(WEOF)) {
	  nextChar = EOF;
	  nextClass = C_EOF;
	  return result;
	}

	numChar++;

	if (nextChar == L'\n')
	  numLines++;
	
	if (nextChar == L'\\' and state == ST)  {
	  nextChar =  getEscapedCharacter(inputStream);
	  numChar++;
	}

	nextClass = wcharToClass(nextChar);

	// If we are a looking for a string or a comment don't repeat
	if (state == ST || comment) return result; 
      }
      
      if (nextChar < 0) {
	std::ostringstream msg;
	msg << "JsonParser::getNextCharacterandClass produced a unrecognized next_char " << nextChar << " |-> " << ((char) wctob(nextChar)) << "\n";
	msg << " In file '" << filename << "' \n";
	msg << "     around line          : " << numLines << "\n";
	msg << "     at character position: " << numChar << "\n";
	throw std::range_error(msg.str());
      }
      
      return result;
    }

    //======================================================================

    void setCommentState() {
      
      switch (currentMode()) {
	  
      case MODE_ARRAY:
      case MODE_OBJECT:   
	switch(state) {
	case VA:    /* value */
	case AR:    /* array */
	  before_comment_state = state;
	  break;
	default:
	  before_comment_state = GO;
	  break;
	}
	break;
      default:
	before_comment_state = state;
	break;
      }
      comment = true;
    }
  
    void unloadBuffer() {
      switch (lastItemType) {
      case JSON_T_STRING:
	ctx.String(buffer);	
	break;
      default:
	break;
      }
      lastItemType = JSON_T_NONE;
      buffer.clear();
    }

    //======================================================================
    
    void performAction(wchar_t nextChar, CharacterClass nextClass, StateType state, ActionType action) {

      if(trace) 
	std::cout << "-------------------------------------------------- performAction( '"<< ((char) wctob(nextChar)) << "' , (" << actionName(action) << ")  )\n";

      //if (trace) std::wcout << L"   buffer contents: '" << buffer.str() << L"' \n";

      switch (action) {
	
      case Consume:
	//if (trace) std::cout << "   consuming " << nextChar <<  " \n";
	return;
	
      case BeginObject:
	if (trace) std::cout << "  Starting a new object. \n";
	ctx.ObjectBegin();
	buffer.clear();
	return;
        
      case BeginArray:
	if (trace) std::cout << "   Starting a new array. \n";
	ctx.ArrayBegin();
	buffer.clear();
	return;
        
      case BeginMatrix:
	if (trace) std::cout << "   Starting a new matrix at " << filename << ":" <<numChar << ". \n";
	ctx.MatrixBegin(filename,numChar);
	buffer.clear();
	return;
        
      case EndObject:
	if (trace) std::cout << "  Ending an object. \n";
	unloadBuffer();
	ctx.ObjectEnd();
	return;
        
      case EndArray:
	if (trace) std::cout << "   Ending an new array. \n";
	unloadBuffer();
	ctx.ArrayEnd();
	return;
        
      case EndMatrix:
	if (trace) std::cout << "   Ending a new matrix at " << numChar << ". \n";
	ctx.MatrixEnd(numChar);
	buffer.clear();
	return;
        
      case DoNext: // Encountered a comma assign last value

	unloadBuffer();
	if (trace) std::wcout << L"   Encountered a comma assign last value to array or object nextChar =  " << nextChar << " \n";
	return; 
 
      case RecordChar:
	buffer.put(nextChar);
	return;

      case RecordString:
	lastItemType = JSON_T_STRING;
	return;

      case RecordKey:
	lastItemType = JSON_T_NONE;
	ctx.Key(buffer.str());
	buffer.clear();
	return;

      case RecordFloat:
	lastItemType = JSON_T_FLOAT;
	ctx.Float(buffer);
	buffer.clear();
	return;

      case RecordInteger:
	lastItemType = JSON_T_INTEGER;
	ctx.Integer(buffer);
	buffer.clear();
	return;

      case RecordTrue:
	lastItemType = JSON_T_TRUE;
	ctx.True();
	return;

      case RecordFalse:
	lastItemType = JSON_T_FALSE;
	ctx.False();
	return;

      case RecordNull:
	lastItemType = JSON_T_NULL;
	ctx.Null();
	return;

      case Abort: {
	std::ostringstream msg;
	msg << "JsonParser::performAction was sent abort from AugmentedStateTranslationTable\n"
	    << "  nextChar   = '" << ((char) wctob(nextChar)) << "'\n" 
	    << "  nextClass  = "  << clsName(nextClass) << "\n"
	    << "  state      = "  << name(state) << "\n"
	    << "  action     = "  << actionName(action) << "\n"
	    << "  character# = " << numChar << "\n"
	    << "  line#      = " << numLines << "\n"
	    << "  context    = " << ctx.CurrentContext();
	throw std::logic_error(msg.str());
      }
      case EndFile: {
	ctx.End(numChar,numLines);
	return;
      }
      default: {
	std::ostringstream msg;
	msg << "JsonParser::performAction was passed a unkown action " << actionName(action) << "\n";
	throw std::logic_error(msg.str());
      }
      }
    }
  
    //======================================================================
    
    /*    */
    bool parseChar(std::wistream& inputStream) {
      
      typedef AugmentedStateTranslationTable::Pair Pair;
      
      std::pair<wchar_t,CharacterClass> next(getNextCharacterAndClass(inputStream));
      wchar_t&        nextChar (next.first);
      CharacterClass& nextClass(next.second);
	
      Pair pair = AugmentedStateTranslationTable::getStateActionPair(state,nextClass);
      
      PsimagLite::Vector<ActionType>::Type& actions(pair.actions);
      

      if(trace) std::cout << "actions.size() = " << actions.size() ;

      for (size_t i=0; i< actions.size(); i++) {
	
	performAction(nextChar, nextClass, state, actions[i]);

	if (actions[i] == EndFile) return false;
      }
      state = pair.newState;
      
      bool doMore = !inputStream.eof();
      return doMore;
    }

    
  };

} // end namespace JsonParser




  /*@}*/
#endif
