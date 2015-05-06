//-*-C++-*-
// Author: Michael S. Summers (ORNL)
/** \ingroup JsonParser */
/*@{*/

/*! \file StatesAndActionsMixin.h
 *
 *
 *
 */

#ifndef  JsonParser_StatesAndActionsMixin_H
#define  JsonParser_StatesAndActionsMixin_H

#include "AllocatorCpu.h"

namespace JsonParser {

  class StatesAndActionsMixin {

  public:

    /* The state and action codes.  */
    typedef enum statesAndActions {
      GO,  /* start    */
      OK,  /* ok       */
      OB,  /* object   */
      KE,  /* key      */
      CO,  /* colon    */
      VA,  /* value    */
      AR,  /* array    */
      ST,  /* string   */
      ES,  /* escape   */
      U1,  /* u1       */
      U2,  /* u2       */
      U3,  /* u3       */
      U4,  /* u4       */
      MI,  /* minus    */
      ZE,  /* zero     */
      IT,  /* integer  */
      FR,  /* fraction */
      E1,  /* e        */
      E2,  /* ex       */
      E3,  /* exp      */
      T1,  /* tr       */
      T2,  /* tru      */
      T3,  /* true     */
      F1,  /* fa       */
      F2,  /* fal      */
      F3,  /* fals     */
      F4,  /* false    */
      N1,  /* nu       */
      N2,  /* nul      */
      N3,  /* null     */
      C1,  /* /        */
      C2,  /* / *     */
      C3,  /* *        */
      FX,  /* *.* *eE* */
      D1,  /* second UTF-16 character decoding started by \ */
      D2,  /* second UTF-16 character proceeded by u */
      NR_STATES,
      __ = -1,  /* Error */
      CN = -2,  /* Colon */
      CM = -3,  /* Comma */
      SE = -4,  /* STRING END */
      BA = -5,  /* Begin Array */
      BO = -6,  /* Begin Object */
      EA = -7,  /* End Array */
      EE = -8,  /* End Empty Object */
      EO = -9,  /* End Object */
      CB = -10, /* comment begin */
      CE = -11, /* comment end */
      FA = -12, /* false */
      TR = -13, /* false */
      NU = -14, /* null */
      DE = -15, /* double detected by exponent e E */
      DF = -16, /* double detected by fraction . */
      SB = -17, /* string begin */
      MX = -18, /* integer detected by minus */
      ZX = -19, /* integer detected by zero */
      IX = -20, /* integer detected by 1-9 */
      EX = -21, /* next char is escaped */
      UC = -22, /* Unicode character read */
    } StatesAndActionsType;

    static bool isState(StatesAndActionsType stateOrAction) {
      if (stateOrAction < NR_STATES)
	return true;
      return false;
    }

    static String  name(StatesAndActionsType stateOrAction) {

      switch (stateOrAction) {
      case GO: return "GO: start    ";
      case OK: return "OK: ok       ";
      case OB: return "OB: object   ";
      case KE: return "KE: key      ";
      case CO: return "CO: colon    ";
      case VA: return "VA: value    ";
      case AR: return "AR: array    ";
      case ST: return "ST: string   ";
      case ES: return "ES: escape   ";
      case U1: return ": u1       ";
      case U2: return ": u2       ";
      case U3: return ": u3       ";
      case U4: return ": u4       ";
      case MI: return "MI: minus    ";
      case ZE: return "ZE: zero     ";
      case IT: return "IT: integer  ";
      case FR: return "FR: fraction ";
      case E1: return ": e        ";
      case E2: return ": ex       ";
      case E3: return ": exp      ";
      case T1: return ": tr       ";
      case T2: return ": tru      ";
      case T3: return ": true     ";
      case F1: return ": fa       ";
      case F2: return ": fal      ";
      case F3: return ": fals     ";
      case F4: return ": false    ";
      case N1: return ": nu       ";
      case N2: return ": nul      ";
      case N3: return ": null     ";
      case C1: return ": /        ";
      case C2: return ": / *      ";
      case C3: return ": *         ";
      case FX: return "FX: *.* *eE* ";
      case D1: return ": second UTF-16 character decoding started by \\ ";
      case D2: return ": second UTF-16 character proceeded by u ";
      case NR_STATES: return "NR_STATES:number of states";
      case __: return "__: Error ";
      case CN: return "CN: Colon Action";
      case CM: return "CM: Comma ";
      case SE: return "SE: STRING END ";
      case BA: return "BA: Begin Array ";
      case BO: return "BO: Begin Object ";
      case EA: return "EA: End Array ";
      case EE: return "EE: End Empty Object ";
      case EO: return "EO: End Object ";
      case CB: return "CB: comment begin ";
      case CE: return "CE: comment end ";
      case FA: return "FA: false ";
      case TR: return "TR: false ";
      case NU: return "NU: null ";
      case DE: return "DE: double detected by exponent e E ";
      case DF: return "DF: double detected by fraction . ";
      case SB: return "SB: string begin ";
      case MX: return "MX: integer detected by minus ";
      case ZX: return "ZX: integer detected by zero ";
      case IX: return "IX: integer detected by 1-9 ";
      case EX: return "EX: next char is escaped ";
      case UC: return "UC: Unicode character read ";
      default: return " Unkown code";
      }
    }

  };

} // end namespace JsonParser


/*@}*/
#endif
