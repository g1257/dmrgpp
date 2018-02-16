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
#include "H5Cpp.h"

namespace PsimagLite {

template<typename T>
struct ToH5 {};

template<>
struct ToH5<double> {
	static const H5::PredType type;
	static const H5T_class_t super;
};

class IoNg {

public:

	class Out {

	public:

		typedef std::vector<String> VectorStringType;

		Out() : hdf5File_(0), groupDef_(0) {}

		Out(const String& fn)
		    : hdf5File_(new H5::H5File(fn, H5F_ACC_TRUNC)),
		      groupDef_(new H5::Group(hdf5File_->createGroup("/Def")))
		{}

		~Out()
		{
			delete groupDef_; // should I close something first? FIXME
			groupDef_ = 0;
			delete hdf5File_; // should I close something first? FIXME
			hdf5File_ = 0;
		}

		const String& filename() const
		{
			// find member that returns filename FIXME
			// assert(hdf5File_);
			// return hdf5File_.filename();
			throw RuntimeError("IoNg:: not implemented\n");
		}

		void open(String const &fn,
		          std::ios_base::openmode mode)
		{
			if (hdf5File_) delete hdf5File_;

			// deal with mode
			hdf5File_ = new H5::H5File(fn, H5F_ACC_TRUNC);

			labels_.clear();

			throw RuntimeError("IoNg:: open cannot handle mode yet\n");
		}

		void close()
		{
			delete groupDef_; // should I close something first? FIXME
			groupDef_ = 0;
			delete hdf5File_; // should I close something first? FIXME
			hdf5File_ = 0;
			labels_.clear();
		}

		void printline(const String &s)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		void printline(OstringStream &s)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		template<typename T>
		void write(const std::vector<T>& v,
		           const String& label)
		{
			assert(hdf5File_);

			SizeType count = findCount(label);
			hsize_t dims[1];
			dims[0] = v.size();
			H5::DataSpace *dataspace = new H5::DataSpace(1, dims); // create new dspace
			H5::DSetCreatPropList dsCreatPlist; // What properties here? FIXME
			String name = "/Def/" + label + ttos(count);
			internalWrite<T>(&(v[0]), name, *dataspace, dsCreatPlist);
			delete dataspace;
			labels_.push_back(label);
		}

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

	private:

		template<typename T>
		void internalWrite(const void *ptr,
		                   String name,
		                   H5::DataSpace& dataspace,
		                   H5::DSetCreatPropList& dsCreatPlist)
		{

			H5::DataSet* dataset = new H5::DataSet(hdf5File_->createDataSet(name,
			                                                                ToH5<T>::type,
			                                                                dataspace,
			                                                                dsCreatPlist));
			dataset->write(ptr, ToH5<T>::type);
			delete dataset;
		}

		SizeType findCount(const String& label) const
		{
			SizeType c = 0;
			SizeType total = labels_.size();
			for (SizeType i = 0; i < total; ++i) {
				if (std::find(labels_.begin(), labels_.end(), label) == labels_.end())
					continue;
				++c;
			}

			return c;
		}

		H5::H5File* hdf5File_;
		H5::Group* groupDef_;
		VectorStringType labels_;
	};

	class In {

	public:

		typedef int long LongIntegerType;
		static const LongIntegerType LAST_INSTANCE=-1;
		typedef unsigned int long LongSizeType;

		In() : hdf5File_(0), groupDef_(0) {}

		In(String const &fn)
		    : hdf5File_(new H5::H5File(fn, H5F_ACC_RDONLY)),
		      groupDef_(new H5::Group(hdf5File_->openGroup("Def")))
		{}

		~In()
		{
			delete groupDef_; // should I close something first? FIXME
			groupDef_ = 0;
			delete hdf5File_; // should I close something first? FIXME
			hdf5File_ = 0;
		}

		void open(String const &fn)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		void close()
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		template<typename X>
		SizeType readline(X &x,const String &s,LongIntegerType level=0)
		{ throw RuntimeError("IoNg:: not implemented\n"); }

		template<typename T>
		void read(std::vector<T>& x,
		          String const &s,
		          LongIntegerType level = 0,
		          bool beQuiet = false)
		{
			assert(hdf5File_);
			assert(groupDef_);

			// proper failure reporting is needed here
			String name = "/Def/" + s + ttos(level);
			H5::DataSet dataset = H5::DataSet(groupDef_->openDataSet(name));
			H5::DataSpace dataspace = dataset.getSpace();
			int rank = dataspace.getSimpleExtentNdims();
			if (rank != 1)
				throw RuntimeError("Reading " + s + " is not a vector\n");

			hsize_t dimsOut[1];
			dataspace.getSimpleExtentDims(dimsOut, NULL);
			x.resize(dimsOut[0]);

			internalRead<T>(&(x[0]), s, dataset);
		}

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

	private:

		template<typename T>
		void internalRead(void* ptr, String label, H5::DataSet& dataset) const
		{
			H5T_class_t typeClass = dataset.getTypeClass();
			if (typeClass != ToH5<T>::super)
				throw RuntimeError("Reading " + label + " has incorrect type\n");
			// H5::FloatType ft = dataset.getFloatType(); // <-- check correct subtype FIXME

			dataset.read(ptr, ToH5<T>::type);
		}

		H5::H5File* hdf5File_;
		H5::Group* groupDef_;
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

