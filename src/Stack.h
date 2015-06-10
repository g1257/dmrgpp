/*
Copyright (c) 2009 , UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."

*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************


*/
/** \ingroup PsimagLite */
/*@{*/

/*! \file Stack.h
 *
 *  std::stack companions
 */

#ifndef PSIMAGLITE_STACK_H_
#define PSIMAGLITE_STACK_H_
#include <stack>
#include "AllocatorCpu.h"
#include "Vector.h"

namespace PsimagLite {

template<typename T>
class Stack {

	typedef std::deque<T,typename Allocator<T>::Type> DequeType_;

public:

	typedef std::stack<T,DequeType_> Type;
}; // class Stack

template<typename T>
class IsStackLike {
public:
	enum {True = false};
};

template<typename T>
class IsStackLike<std::stack<T,std::deque<T,typename Allocator<T>::Type> > > {
public:
	enum {True = true};
};

template<typename StackType>
typename EnableIf<IsStackLike<StackType>::True,std::ostream>::Type&
operator<<(std::ostream& os,const StackType& st)
{
	StackType st2 = st;
	os<<st2.size()<<"\n";
	while(!st2.empty()) {
		typename StackType::value_type x = st2.top();
		os<<x<<"\n";
		st2.pop();
	}
	return os;

}

template<typename StackType>
typename EnableIf<IsStackLike<StackType>::True,std::istream>::Type&
operator>>(std::istream& is,StackType& x)
{
	typedef typename StackType::value_type ValueType;
	typename Vector<ValueType>::Type tmpVec;
	is>>tmpVec;
	for (int i=tmpVec.size()-1;i>=0;i--) {
		x.push(tmpVec[i]);
	}
	return is;
}

} // namespace PsimagLite

/*@}*/
#endif // PSIMAGLITE_STACK_H_
