//-*-C++-*-
// Author: Michael S. Summers (ORNL)
/** \ingroup JsonParser */
/*@{*/

/*! \file StatesMixin.h  
 *
 *  
 *
 */

#ifndef  JsonParser_StatesMixin_H
#define  JsonParser_StatesMixin_H

#include "String.h"

namespace JsonParser {

  class StatesMixin {

  public:

    /* The state and action codes.  */
    typedef enum states {
      GO,  /* start    */
      VA,  /* looking for a value    */
      EV,  /* Looking after a value  */
      BS,  /* Looking for the begining of a string  */
      ST,  /* Looking for string characters         */
      IT,  /* Looking for an integer  */
      FR,  /* Looking for the integer after the point */
      EX,  /* Looking for the exponent */
      EX2, /* Looking for the integer part of the exponent */
      T,   /* Looking for _rue (true)   */
      TR,  /* Looking for __ue (true)        */
      TRU, /* Looking for ___e (true)       */
      F,   /* Looking for _alse (false)   */
      FA,  /* Looking for __lse (false)   */
      FAL, /* Looking for ___se (false)   */
      FALS,/* Looking for ____e (false)   */
      N,   /* Looking for _ull (null)     */
      NU,  /* Looking for __ll (null)     */
      NUL, /* Looking for ___l (null)     */
      A,   /* Looking for _rray (array)     */
      AR,  /* Looking for __ray (array)     */
      ARR, /* Looking for ___ay (array)     */
      ARRA, /* Looking for ___y (array)     */
      ARRAY, /* Looking for ____( (array)     */
      MS,    /* Skipping matrix chars     */
      END,  /* END of Parsing               */
      NR_STATES
    } StateType;
    
    static PsimagLite::String  name(StateType state) {

      switch (state) {
      case GO: return "GO: start    ";
      case VA: return "VA: looking for a value    ";
      case EV: return "EV: Looking after a value  ";
      case BS: return "BS: Looking for the begining of a string  ";
      case ST: return "ST: Looking for string characters         ";
      case IT: return "IT: Looking for an integer  ";
      case FR: return "FR: Looking for the integer after the point ";
      case EX: return "EX: Looking for the exponent ";
      case EX2: return "EX2: Looking for the integer part of the exponent ";
      case T: return "T: Looking for _rue (true)   ";
      case TR: return "TR: Looking for __ue (true)        ";
      case TRU: return "TRU: Looking for ___e (true)       ";
      case F: return "F: Looking for _alse (false)   ";
      case FA: return "FA: Looking for __lse (false)   ";
      case FAL: return "FAL: Looking for ___se (false)   ";
      case FALS: return "FALS: Looking for ____e (false)   ";
      case N: return "N: Looking for _ull (null)     ";
      case NU: return "NU: Looking for __ll (null)     ";
      case NUL: return "NUL: Looking for ___l (null)     ";
      case A: return "A: Looking for _rray (array)     ";
      case AR: return "AR: Looking for __ray (array)     ";
      case ARR: return "ARR: Looking for ___ay (array)     ";
      case ARRA: return "ARRA: Looking for ___y (array)     ";
      case ARRAY: return "ARRAY: Looking for ____( (array)     ";
      case MS: return "skipping matrix chars     ";
      case END: return "END: End of Parsing                ";
      default: return " Unkown state code";
      }
    }

  };

} // end namespace JsonParser


/*@}*/
#endif
