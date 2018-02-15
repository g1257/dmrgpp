/*
Copyright (c) 2009-2018, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 2.]
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

/*! \file IoNg.h
 *
 *  This class handles Input/Output for PsimagLite
 */

#ifndef PSI_IO_NG_H
#define PSI_IO_NG_H

#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "Stack.h"
#include "Map.h"

namespace PsimagLite {

class IoNg {

public:

	class Out {

	public:

		Out()  { throw RuntimeError("IoNg:: not implemented\n"); }

		Out(std::ostream& os)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		Out(const String& fn)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		const String& filename() const { throw RuntimeError("IoNg:: not implemented\n"); }

		void open(String const &fn,
		          std::ios_base::openmode mode)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		void close()
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		void printline(const String &s)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		void printline(OstringStream &s)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		template<typename X>
		void printVector(X const &x,
		                 String const &label,
		                 typename EnableIf<IsVectorLike<X>::True, int>::Type = 0)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		template<class T>
		void print(const String& label, const T& something)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		void print(const String& something)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		template<typename SomePrintableType>
		void print(const SomePrintableType& something)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		template<typename X>
		void printMatrix(const X& mat,
		                 String const &s,
		                 typename EnableIf<IsMatrixLike<X>::True, int>::Type = 0)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		int rank() { throw RuntimeError("IoNg:: not implemented\n"); }

		void flush()
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		void setPrecision(SizeType x)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		template<typename ActionType>
		void action(ActionType& a)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		template<typename T>
		friend Out& operator<<(Out& io, const T& t)
		{ throw RuntimeError("IoNg:: not implemented\n"); }
	};

	class In {

	public:

		typedef int long LongIntegerType;
		static const LongIntegerType LAST_INSTANCE=-1;
		typedef unsigned int long LongSizeType;

		In() { throw RuntimeError("IoNg:: not implemented\n"); }

		In(String const &fn)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		void open(String const &fn)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		void close()
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		template<typename X>
		SizeType readline(X &x,const String &s,LongIntegerType level=0)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		template<typename X>
		typename EnableIf<IsVectorLike<X>::True,std::pair<String,SizeType> >::Type
		read(X &x,
		     String const &s,
		     LongIntegerType level=0,
		     bool beQuiet = false)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		template<typename X>
		typename EnableIf<IsStackLike<X>::True,std::pair<String,SizeType> >::Type
		read(X &x,
		     String const &s,
		     LongIntegerType level=0,
		     bool beQuiet = false)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		//! Assumes something of the form
		//! label[key]=value
		template<typename MapType>
		typename EnableIf<IsMapLike<MapType>::True,void>::Type
		read(MapType& x,
		     String const &s,
		     LongIntegerType level=0)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		template<typename X>
		std::pair<String,SizeType> readKnownSize(X &x,
		                                         String const &s,
		                                         LongIntegerType level=0)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		std::pair<String,SizeType> advance(String const &s,
		                                   LongIntegerType level=0,
		                                   bool beQuiet=false)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		void readFullLine(String& temp)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		LongSizeType advanceToLine(LongSizeType line)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		SizeType count(const String& s)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		template<typename X,template<typename> class SomeType>
		void readSparseVector(SomeType<X> &x,
		                      String const &s,
		                      LongIntegerType level=0)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		template<typename X>
		void readMatrix(X &mat,
		                String const &s,
		                LongIntegerType level= 0,
		                typename EnableIf<IsMatrixLike<X>::True, int>::Type  = 0)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		template<
		        typename FieldType,
		        template <typename> class SparseMatrixTemplate,
		        template<typename,template<typename> class>
		        class X>
		void readMatrix(X<FieldType,SparseMatrixTemplate>& op,
		                const String& s,
		                LongIntegerType level=0)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		void rewind()
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		void move(int x)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		bool eof() const { throw RuntimeError("IoNg:: not implemented\n"); }

		const char* filename() const
		{ throw RuntimeError("IoNg:: not implemented\n"); }


		template<typename T>
		friend void operator>>(In& io, T& t)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

	};
}; //class IoNg

template<>
struct IsInputLike<IoNg::In> {
	enum {True = true};
};

template<>
struct IsOutputLike<IoNg::Out> {
	enum {True = true};
};

} // namespace PsimagLite

/*@}*/
#endif

