//-*-C++-*-
// Author: Michael S. Summers (ORNL)
//
#ifndef PSIMAG_Assert_H
#define PSIMAG_Assert_H

#include<iostream>
#include<exception>
#include <sstream>

namespace psimag { 

// ================================== COMPILE TIME ASSERTION ================
/** Compile time assertion methods taken from: 
 *   A. Alexandrescu, Modern C++ Design
 */
template<bool> struct CompileTimeError;
template<> struct CompileTimeError<true> { };
#define STATIC_CHECK(expr)  {CompileTimeError<(expr) != 0>(); }

#  if 0
   this method would give better error messages but seems not to work
   with g++
   
   template<bool> struct CompileTimeChecker
   {
   CompileTimeChecker(...) {};
   };
   template<> struct CompileTimeChecker<false> { };
   #define STATIC_CHECK(expr, msg) \
   {\
   class ERROR_##msg {};\
   (void)sizeof(CompileTimeChecker<(expr)>(ERROR_##msg()));\
   }
#  endif   


// ================================== RUN TIME ASSERTION ====================
#ifndef NDEBUG
  
/** \define Runtime Assertion Macro
 *  \author Hwee Kuan Lee
 */
#define ASSERT( a, e )  psimag::_Assert(__FILE__, __LINE__, a, e )


/** If assertion A fails exception E is thrown. 
 *  Code generated only if NDEBUG not defined
 *  \author Hwee Kuan Lee and Thomas Schulthess
 */
template < class A, class E > 
inline 
void _Assert(const char * filename, int lineno, A assertion, E except)
{
  if (!assertion) 
    {
      PsimagLite::OstringStream msg;
      msg << "ASSERTION FAILED AT:"<<filename<<":"<<lineno<<"\n";
      msg << except.what() << std::endl;
      throw E(msg.str());
    }
}

template < class A >
inline
void _Assert(const char * filename, int lineno, 
             A assertion, std::exception except=std::exception())
{
  if (!assertion) 
    {
      PsimagLite::OstringStream msg;
      msg << "ASSERTION FAILED AT:"<<filename<<":"<<lineno<<"\n";
      throw std::exception(msg.str());
    }
}

#endif /* ifndef NDEBUG */


#ifdef NDEBUG

// Does nothing if NDEBUG is defined
#define ASSERT(expr,except)

#endif // NDEBUG

}      // namespace psimag 


#endif // PSIMAG_Assert_H
