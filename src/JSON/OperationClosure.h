//-*- mode: C++; -*-
// Author: Michael S. Summers (ORNL)
//
/** \ingroup PSIMAG */
/*@{*/

/*! \file OperationClosure.h  
 *
 */

#ifndef PSIMAG_Operation_Closures_H
#define PSIMAG_Operation_Closures_H


namespace psimag {
  
  class OP {
  public:
    typedef enum{FourierTransform,Integrate,PLUS,MINUS,TIMES,DIVIDE,INV,NORM} Type;
  };

  //====================================================================== 
  
  template<typename Operand1Type,
	   int      Operator,
	   typename Operand2Type>
  class OperationClosure {
  public:

    const Operand1Type& lhs;
    const Operand2Type& rhs;
    
    OperationClosure(const Operand1Type& l,
		     const Operand2Type& r):
      lhs(l),
      rhs(r)
    {}
  };      

  //====================================================================== 
  
  template<int      Operator,
	   typename OperandType>
  class UnaryOperationClosure {
  public:
    
    const OperandType& operand;
    SizeType iterationCount;
    PsimagLite::String msg;

    UnaryOperationClosure(const OperandType& theOperand):
      operand(theOperand),
      iterationCount(0),
      msg("")
    {}
  };     

  //====================================================================== OP::TIMES Compound Closures

  template<typename T,
	   typename Operand1Type,
	   int      Operator,
	   typename Operand2Type>
  OperationClosure<T,OP::TIMES,OperationClosure<Operand1Type,Operator,Operand2Type> > operator * (T& lhs, 
											OperationClosure<Operand1Type,Operator,Operand2Type> rhs) {
    return OperationClosure<T,OP::TIMES, OperationClosure<Operand1Type,Operator,Operand2Type> >(lhs, rhs);
  }

  //======================================================================
  
  template<typename Operand1Type,
	   int      Operator,
	   typename Operand2Type,
	   typename T>
  OperationClosure<OperationClosure<Operand1Type,Operator,Operand2Type>,OP::TIMES,T > operator * (OperationClosure<Operand1Type,Operator,Operand2Type> lhs, 
											T& rhs) {
    return OperationClosure<OperationClosure<Operand1Type,Operator,Operand2Type>,OP::TIMES, T>(lhs, rhs);
  }

  //====================================================================== OP::PLUS Compound Closures

  template<typename Operand11Type,
	   int      Operator1,
	   typename Operand12Type,
	   typename Operand21Type,
	   int      Operator2,
	   typename Operand22Type>
  OperationClosure<OperationClosure<Operand11Type,Operator1,Operand12Type>,
		   OP::PLUS,
		   OperationClosure<Operand21Type,Operator2,Operand22Type> > operator + (OperationClosure<Operand11Type,Operator1,Operand12Type> lhs, 
											 OperationClosure<Operand21Type,Operator2,Operand22Type> rhs) {
    return OperationClosure<OperationClosure<Operand11Type,Operator1,Operand12Type>,OP::PLUS,OperationClosure<Operand21Type,Operator2,Operand22Type> >(lhs, rhs);
  }

} // end namespace DCA


/*@}*/
#endif
