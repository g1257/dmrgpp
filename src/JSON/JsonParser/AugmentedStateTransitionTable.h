//-*-C++-*-
// Author: Michael S. Summers (ORNL)
/** \ingroup JsonParser */
/*@{*/

/*! \file AugmentedStateTranslationTable.h  
 *
 *  
 *
 */

#ifndef  JsonParser_AugmentedStateTranslationTable_H
#define  JsonParser_AugmentedStateTranslationTable_H

#include "StatesMixin.h"
#include "ActionsMixin.h"
#include "CharacterMapper.h"

namespace JsonParser {

  class AugmentedStateTranslationTable: 
    public CharacterMapper, 
    public StatesMixin,
    public ActionsMixin
  {
  public:
    
    class Pair {
    public:
      StateType               newState;
      PsimagLite::Vector<ActionType>::Type actions;
      
      Pair(const Pair& other):
	newState(other.newState),
	actions(other.actions)
      {}

      Pair(StateType s, ActionType a):
	newState(s),
	actions(1,a)
      {}

      Pair(StateType s, ActionType a, ActionType b):
	newState(s),
	actions(2,a)
      {
	actions[1] = b;
      }
    };

    static Pair getStateActionPair(StateType state, CharacterClass cls) {

      if (!(state < NR_STATES) || !(cls < NR_CLASSES) ) {
	std::ostringstream msg;
	msg << "StateTranslationTable ("<< state << " < " << NR_STATES << "," << cls << " < " << NR_CLASSES << ")\n";
	msg << "StateTranslationTable ("<< name(state) << "," << clsName(cls) << ")\n";
	throw std::logic_error(msg.str());
      }
      
      switch (state) {
	
      case GO:  // Looking for a start
	switch (cls) {
	case C_SPACE: return Pair(GO,Consume);
	case C_WHITE: return Pair(GO,Consume);
	case C_LCURB: return Pair(BS,BeginObject);
	case C_LSQRB: return Pair(VA,BeginArray);
	default:      return Pair(END,Abort);
	}

      case VA:  // Looking for a value
	switch (cls) {
	case C_SPACE: return Pair(VA,Consume);
	case C_WHITE: return Pair(VA,Consume);
	case C_LCURB: return Pair(BS,BeginObject);
	case C_LSQRB: return Pair(VA,BeginArray);
	case C_QUOTE: return Pair(ST,Consume);
	case C_LOW_A: return Pair(A,Consume);
	case C_MINUS:
	case C_PLUS:
	case C_ZERO:
	case C_DIGIT:
	  return Pair(IT,RecordChar);
	case C_POINT:
	  return Pair(FR,RecordChar);
	case C_LOW_T: return Pair(T,Consume);
	case C_LOW_F: return Pair(F,Consume);
	case C_LOW_N: return Pair(N,Consume);
	default:      return Pair(END,Abort);
	}

      case EV:  // End Value, Looking for a next
	switch (cls) {
	case C_SPACE: return Pair(EV,Consume);
	case C_WHITE: return Pair(EV,Consume);
	case C_RCURB: return Pair(EV,EndObject);
	case C_RSQRB: return Pair(EV,EndArray);
	case C_COMMA: return Pair(VA,DoNext);
	case C_COLON: return Pair(VA,RecordKey);
	case C_EOF:   return Pair(EV,EndFile);
	default:      return Pair(END,Abort);
	}

      case BS: // Looking for start of string
	switch (cls) {
	case C_SPACE: return Pair(BS,Consume);
	case C_WHITE: return Pair(BS,Consume);
	case C_QUOTE: return Pair(ST,Consume);
	default:      return Pair(END,Abort);
	}
      case ST: // Reading String
	switch (cls) {
	case C_QUOTE: return Pair(EV,RecordString);
	default:      return Pair(ST,RecordChar);
	}

      case IT: // Integer
	switch (cls) {
	case C_ZERO:
	case C_DIGIT:
	  return Pair(IT,RecordChar);
	case C_POINT:
	  return Pair(FR,RecordChar);
	case C_E:
	  return Pair(EX,RecordChar);
	case C_LOW_E:
	  return Pair(EX,RecordChar);
	case C_RCURB: return Pair(EV,RecordInteger,EndObject);
	case C_RSQRB: return Pair(EV,RecordInteger,EndArray);
	case C_COMMA: return Pair(VA,RecordInteger,DoNext);
	default:      
	  return Pair(END,Abort); 
	}
	
      case FR: //   Problem with accepting things to short!?
	switch (cls) {
	case C_ZERO:
	case C_DIGIT:
	  return Pair(FR,RecordChar);
	case C_POINT:
	  return Pair(END,Abort);
	case C_E:
	  return Pair(EX,RecordChar);
	case C_LOW_E:
	  return Pair(EX,RecordChar);
	case C_RCURB: return Pair(EV,RecordFloat,EndObject);
	case C_RSQRB: return Pair(EV,RecordFloat,EndArray);
	case C_COMMA: return Pair(VA,RecordFloat,DoNext);
	default: 
	  return Pair(END,Abort); 
	}
	
      case EX: // True
	switch (cls) {
	case C_MINUS:
	case C_PLUS:
	case C_ZERO:
	case C_DIGIT:
	  return Pair(EX2,RecordChar);
	default:      
	  return Pair(END,Abort);
	}
	
      case EX2: // True
	switch (cls) {
	case C_ZERO:
	case C_DIGIT:
	  return Pair(EX2,RecordChar);
	case C_RCURB: return Pair(EV,RecordFloat,EndObject);
	case C_RSQRB: return Pair(EV,RecordFloat,EndArray);
	case C_COMMA: return Pair(VA,RecordFloat,DoNext);
	default: 
	  return Pair(END,Abort);
	}
	
      case T: // True
	switch (cls) {
	case C_LOW_R: return Pair(TR,Consume);
	default:      return Pair(END,Abort);
	}
      case TR: // True
	switch (cls) {
	case C_LOW_U: return Pair(TRU,Consume);
	default:      return Pair(END,Abort);
	}
      case TRU: // True
	switch (cls) {
	case C_LOW_E: return Pair(EV,RecordTrue);
	default:      return Pair(END,Abort);
	}
	
      case N: // Null
	switch (cls) {
	case C_LOW_U: return Pair(NU,Consume);
	default:      return Pair(END,Abort);
	}
      case NU: // Null
	switch (cls) {
	case C_LOW_L: return Pair(NUL,Consume);
	default:      return Pair(END,Abort);
	}
      case NUL: // Null
	switch (cls) {
	case C_LOW_L: return Pair(EV,RecordNull);
	default:      return Pair(END,Abort);
	}
	
      case F: // False
	switch (cls) {
	case C_LOW_A: return Pair(FA,Consume);
	default:      return Pair(END,Abort);
	}
      case FA: // False
	switch (cls) {
	case C_LOW_L: return Pair(FAL,Consume);
	default:      return Pair(END,Abort);
	}
      case FAL: // False
	switch (cls) {
	case C_LOW_S: return Pair(FALS,Consume);
	default:      return Pair(END,Abort);
	}
      case FALS: // False
	switch (cls) {
	case C_LOW_E: return Pair(EV,RecordFalse);
	default:      return Pair(END,Abort);
	}

      case A: // False
	switch (cls) {
	case C_LOW_R: return Pair(AR,Consume);
	default:      return Pair(END,Abort);
	}
      case AR: // False
	switch (cls) {
	case C_LOW_R: return Pair(ARR,Consume);
	default:      return Pair(END,Abort);
	}
      case ARR: // False
	switch (cls) {
	case C_LOW_A: return Pair(ARRA,Consume);
	default:      return Pair(END,Abort);
	}
      case ARRA: // False
	switch (cls) {
	case C_LOW_Y: return Pair(ARRAY,Consume);
	default:      return Pair(END,Abort);
	}
      case ARRAY: // False
	switch (cls) {
	case C_LPARN: return Pair(MS,Consume,BeginMatrix);
	default:      return Pair(END,Abort);
	}

      case MS: // False
	switch (cls) {
	case C_RPARN: return Pair(EV,EndMatrix);
	default:      return Pair(MS,Consume);
	}

      default:
	std::ostringstream msg;
	msg << "AugmentedStateTranslationTable::getStateActionPair was passed a unkown state " << name(state) << "\n";
	throw std::logic_error(msg.str());
      }
    }
  };

} // end namespace JsonParser


/*@}*/
#endif
