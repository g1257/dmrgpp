/*
Copyright (c) 2009-2015, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 3.0]
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
/** \ingroup DMRG */
/*@{*/

/*! \file VectorWithOffset.h
 *
 *  A class to represent a vector like this 000000 XXXXXXXX 0000000000000
 *  offset_ is where the first X (non-zero element) is.
 *  data_ contains the nonzero part.
 *  sizE_ is the size of the vector
 */
#ifndef VECTOR_WITH_OFFSET_H
#define VECTOR_WITH_OFFSET_H
#include "Vector.h"
#include "ProgressIndicator.h"

namespace Dmrg {
template<typename FieldType>
struct VectorWithOffset {
public:
	typedef FieldType value_type;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef typename PsimagLite::Vector<FieldType>::Type VectorType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	static const FieldType zero_;

	VectorWithOffset()
	    : progress_("VectorWithOffset"), size_(0), offset_(0), m_(0)
	{}

	template<typename SomeBasisType>
	VectorWithOffset(const VectorSizeType& weights,
	                 const SomeBasisType& someBasis)
	    : progress_("VectorWithOffset"),size_(someBasis.size())
	{
		bool found = false;
		for (SizeType i=0;i<weights.size();i++) {
			if (weights[i]>0) {
				if (found) {
					PsimagLite::String msg("FATAL: VectorWithOffset::");
					msg += " more than one non-zero sector found. ";
					msg += " Maybe you should be using VectorWithOffsets instead?\n";
					throw PsimagLite::RuntimeError(msg);
				}

				data_.resize(weights[i]);
				offset_ = someBasis.partition(i);
				m_ = i;
				found = true;
			}
		}
	}

	void resize(SizeType x)
	{
		size_ = x;
		data_.clear();
		offset_=0;
		m_=0;
	}

	template<typename SomeBasisType>
	void set(const typename PsimagLite::Vector<VectorType>::Type& v,
	         const SomeBasisType& someBasis)
	{
		bool found = false;
		size_ = someBasis.size();
		for (SizeType i=0;i<v.size();i++) {
			if (v[i].size()>0) {
				if (found) {
					PsimagLite::String msg("FATAL: VectorWithOffset::");
					msg += " more than one non-zero sector found. ";
					msg += " Maybe you should be using VectorWithOffsets instead?\n";
					throw PsimagLite::RuntimeError(msg);
				}

				data_ = v[i];
				offset_ = someBasis.partition(i);
				m_ = i;
				found = true;
			}
		}

		if (!found) throw PsimagLite::RuntimeError("Set failed\n");
	}

	template<typename SomeBasisType>
	void fromFull(const VectorType& v,const SomeBasisType& someBasis)
	{
		size_ = someBasis.size();
		try {
			m_ = findPartition(v,someBasis);
			offset_ = someBasis.partition(m_);
			SizeType total = someBasis.partition(m_+1) - offset_;
			data_.resize(total);
			for (SizeType i=0;i<total;i++) data_[i] = v[i+offset_];
		} catch (std::exception& e) {
			std::cout<<e.what();
			m_=0;
			offset_=0;
			data_.resize(0);
		}
	}

	SizeType sectors() const { return (size_ == 0) ? 0 : 1; }

	SizeType sector(SizeType) const { return m_; }

	SizeType offset(SizeType) const { return offset_; }

	SizeType effectiveSize(SizeType) const { return data_.size(); }

	void setDataInSector(const VectorType& v,SizeType)
	{
		data_=v;
	}

	void extract(VectorType& v, SizeType = 0) const
	{
		v=data_;
	}

	template<typename SparseVectorType>
	void toSparse(SparseVectorType& sv) const
	{
		sv.resize(size_);
		for (SizeType i=0;i<data_.size();i++)
			sv[i+offset_] = data_[i];
	}

	template<typename IoOutputter>
	void save(IoOutputter& io,const PsimagLite::String& label) const
	{
		io.printline(label);
		PsimagLite::String s="#size="+ttos(size_);
		io.printline(s);
		s="#offset="+ttos(offset_);
		io.printline(s);
		s="#m="+ttos(m_);
		io.printline(s);
		io.printVector(data_,"#data");
	}

	template<typename IoInputter>
	void load(IoInputter& io,const PsimagLite::String& label,SizeType counter=0)
	{
		io.advance(label,counter);
		int x = 0;
		io.readline(x,"#size=");
		if (x<0)
			throw PsimagLite::RuntimeError("VectorWithOffset::load(...): size<0\n");
		size_ = x;
		io.readline(x,"#offset=");
		if (x<0)
			throw PsimagLite::RuntimeError("VectorWithOffset::load(...): offset<0\n");
		offset_ = x;
		io.readline(x,"#m=");
		if (x<0)
			throw PsimagLite::RuntimeError("VectorWithOffset::load(...): m<0\n");
		m_ = x;
		io.read(data_,"#data");
	}

	template<typename IoInputter>
	void loadOneSector(IoInputter& io,
	                   const PsimagLite::String& label,
	                   SizeType counter=0)
	{
		load(io,label,counter);
	}

	template<typename SomeBasisType>
	void populateSectors(const SomeBasisType&)
	{
		throw PsimagLite::RuntimeError("VectorWithOffset cannot populateSectors\n");
	}

	template<typename SomeBasisType>
	void populateFromQns(const VectorSizeType& qns,
	                     const SomeBasisType& someBasis)
	{
		if (qns.size() == 0) return;
		if (qns.size() != 1) {
			PsimagLite::String msg("FATAL: VectorWithOffset::populateFromQns: ");
			throw PsimagLite::RuntimeError(msg + "more than one sector found\n");
		}

		size_ = someBasis.size();

		SizeType q = qns[0];
		m_ = findPartitionWithThisQn(q,someBasis);
		offset_ = someBasis.partition(m_);
		SizeType total = someBasis.partition(m_+1) - offset_;
		VectorType tmpV(total,0);
		data_ = tmpV;

		PsimagLite::OstringStream msg;
		msg<<"populateFromQns succeeded";
		progress_.printline(msg,std::cout);
	}

	void collapseSectors() {}

	SizeType size() const { return size_; }

	SizeType effectiveSize() const { return data_.size(); }

	SizeType offset() const { return offset_; }

	VectorWithOffset<FieldType> operator+=(const VectorWithOffset<FieldType>& v)
	{
		if (size_ == 0 && offset_ == 0 && m_ == 0) {
			data_ = v.data_;
			size_ = v.size_;
			offset_ = v.offset_;
			m_ = v.m_;
			return *this;
		}

		if (size_ != v.size_ || offset_ != v.offset_ || m_ != v.m_)
			throw PsimagLite::RuntimeError("VectorWithOffset::operator+=\n");

		data_ += v.data_;

		return *this;
	}

	const FieldType& slowAccess(SizeType i) const
	{
		assert(i>=offset_ && i<(offset_+data_.size()));
		return data_[i-offset_];
	}

	FieldType& slowAccess(SizeType i)
	{
		assert(i >= offset_);
		assert(i-offset_ < data_.size());
		return data_[i-offset_];
	}

	const FieldType& fastAccess(SizeType,SizeType j) const
	{
		assert(j<data_.size());
		return data_[j];
	}

	FieldType& fastAccess(SizeType,SizeType j)
	{
		assert(j<data_.size());
		return data_[j];
	}

	int index2Sector(SizeType i) const
	{
		return ((i < offset_) || (i >= (offset_+data_.size()))) ? (-1) : (0);
	}

	template<typename FieldType2>
	friend FieldType2 norm(const Dmrg::VectorWithOffset<FieldType2>& v);

	template<typename FieldType2>
	friend FieldType2 norm(const Dmrg::VectorWithOffset<std::complex<FieldType2> >& v);

	template<typename FieldType2>
	friend FieldType2 operator*(const Dmrg::VectorWithOffset<FieldType2>& v1,
	                            const Dmrg::VectorWithOffset<FieldType2>& v2);

	template<typename FieldType3,typename FieldType2>
	friend VectorWithOffset<FieldType2> operator*(const FieldType3& value,
	                                              const VectorWithOffset<FieldType2>& v);

private:

	template<typename SomeBasisType>
	SizeType findPartitionWithThisQn(SizeType qn,
	                                 const SomeBasisType& someBasis) const
	{
		SizeType np = someBasis.partition()-1;
		for (SizeType i=0;i<np;i++) {
			SizeType state = someBasis.partition(i);
			if (SizeType(someBasis.qn(state))==qn) return i;
		}

		throw PsimagLite::RuntimeError("VectorWithOffset: findPartitionWithThisQn\n");
	}

	template<typename SomeBasisType>
	SizeType findPartition(const VectorType& v,const SomeBasisType& someBasis) const
	{
		bool found = false;
		SizeType p = 0;
		for (SizeType i=0;i<someBasis.partition()-1;i++) {
			if (nonZeroPartition(v,someBasis,i)) {
				if (found) {
					PsimagLite::String msg("FATAL: VectorWithOffset::");
					msg += " More than one partition found\n";
					throw PsimagLite::RuntimeError(msg);
				}

				found = true;
				p = i;
			}
		}
		if (!found) {
			PsimagLite::String msg("VectorWithOffset::");
			msg += " No partition found\n";
			throw PsimagLite::RuntimeError(msg);
		}

		return p;
	}

	template<typename SomeBasisType>
	bool nonZeroPartition(const VectorType& v,
	                      const SomeBasisType& someBasis,SizeType i) const
	{
		typename VectorType::value_type zero = 0;
		for (SizeType j=someBasis.partition(i);j<someBasis.partition(i+1);j++) {
			if (v[j]!=zero) return true;
		}
		return false;
	}

	PsimagLite::ProgressIndicator progress_;
	SizeType size_;
	VectorType data_;
	SizeType offset_;
	SizeType m_; // partition
}; // class VectorWithOffset

template<typename FieldType>
const FieldType VectorWithOffset<FieldType>::zero_=0;

template<typename FieldType>
FieldType operator*(const Dmrg::VectorWithOffset<FieldType>& v1,
                    const Dmrg::VectorWithOffset<FieldType>& v2)
{
	if (v1.m_ != v2.m_) return 0.0;
	return (v1.data_ * v2.data_);
}

template<typename FieldType,typename FieldType2>
VectorWithOffset<FieldType2> operator*(const FieldType& value,
                                       const VectorWithOffset<FieldType2>& v)
{
	VectorWithOffset<FieldType2> w = v;
	w.data_ *= value;
	return w;
}

template<typename FieldType>
FieldType norm(const Dmrg::VectorWithOffset<FieldType>& v)
{
	return PsimagLite::norm(v.data_);
}

template<typename FieldType>
FieldType norm(const Dmrg::VectorWithOffset<std::complex<FieldType> >& v)
{
	return PsimagLite::norm(v.data_);
}

}
/*@}*/
#endif

