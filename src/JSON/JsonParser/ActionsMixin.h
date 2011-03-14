//-*-C++-*-

/** \ingroup JsonParser */
/*@{*/

/*! \file ActionsMixin.h  
 *
 *  
 *
 */

#ifndef  JsonParser_ActionsMixin_H
#define  JsonParser_ActionsMixin_H

#include <string>

namespace JsonParser {
  
  class ActionsMixin {
    
  public:

    /* The state and action codes.  */
    typedef enum actions {
      Consume,       /* Consume the character    */
      BeginObject,   /* BeginObject              */
      BeginArray,    /* BeginArray               */
      BeginMatrix,   /* BeginMatrix              */
      EndObject,     /* endObject              */
      EndArray,      /* EndArray               */
      EndMatrix,      /* EndMatrix               */
      DoNext,        /* DoNext                 */
      RecordKey,     /* RecordKey                */
      RecordChar,    /* RecordChar               */
      RecordString,  /* RecordString             */
      RecordInteger, /* RecordInteger            */
      RecordFloat,   /* RecordFloat            */
      RecordTrue,    /* RecordTrue            */
      RecordFalse,   /* RecordFalse            */
      RecordNull,    /* RecordNull            */
      EndFile,       /* EndFile               */
      Abort,         /* Abort Parsing            */
      NR_ACTIONS
    } ActionType;
    
    static std::string  actionName(ActionType action) {
      
      switch (action) {
      case Consume: return " Consume the character    ";
      case BeginObject: return " BeginObject              ";
      case BeginArray: return " BeginArray               ";
      case BeginMatrix: return " BeginMatrix               ";
      case EndObject: return " endObject              ";
      case EndArray: return " EndArray               ";
      case EndMatrix: return " EndMatrix               ";
      case DoNext: return " DoNext                 ";
      case RecordKey: return " RecordKey                ";
      case RecordChar: return " RecordChar               ";
      case RecordString: return " RecordString             ";
      case RecordInteger: return " RecordInteger            ";
      case RecordFloat: return " RecordFloat            ";
      case RecordTrue: return " RecordTrue            ";
      case RecordFalse: return " RecordFalse            ";
      case RecordNull: return " RecordNull            ";
      case EndFile: return " EndFile            ";
      case Abort: return " Abort Parsing            ";
      default: return " Unkown action code";
      }
    }

  };

} // end namespace JsonParser


/*@}*/
#endif
