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

/*! \file VectorWithOffsets.h
 *
 *  A class to represent a vector like this 000000 XXXXXXXX 0000XXXXXXX000000000....
 *  offsets_ is where the first X (non-zero element) is in each block of Xs.
 *  data_ contains a vector of vectors for each nonzero part.
 *  size_ is the size of the whole vector
 */
#ifndef VECTOR_WITH_OFFSETS_H
#define VECTOR_WITH_OFFSETS_H
#include "Complex.h"
#include "ProgressIndicator.h"
#include <cassert>
#include "ProgramGlobals.h"

// FIXME: a more generic solution is needed instead of tying
// the non-zero structure to basis
namespace Dmrg {
template<typename ComplexOrRealType>
class VectorWithOffsets {

	typedef VectorWithOffsets<ComplexOrRealType> ThisType;
	static ComplexOrRealType const zero_;

public:

	typedef ComplexOrRealType value_type;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef std::pair<SizeType,SizeType> PairSizeType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;

	VectorWithOffsets()
	    : progress_("VectorWithOffsets"),size_(0),index2Sector_(0)
	{ }

	template<typename SomeBasisType>
	VectorWithOffsets(const typename PsimagLite::Vector<SizeType>::Type& weights,
	                  const SomeBasisType& someBasis)
	    : progress_("VectorWithOffsets"),
	      size_(someBasis.size()),
	      index2Sector_(size_),
	      data_(weights.size()),
	      offsets_(weights.size()+1)
	{
		for (SizeType i=0;i<weights.size();i++) {
			data_[i].resize(weights[i]);
			offsets_[i] = someBasis.partition(i);
			if (weights[i]>0) {
				SizeType qn = someBasis.pseudoEffectiveNumber(offsets_[i]);
				nzMsAndQns_.push_back(PairSizeType(i,qn));
				//firstSector_ = i;
			}
		}

		offsets_[weights.size()]=size_;
		setIndex2Sector();
	}

	void resize(SizeType x)
	{
		size_ = x;
		index2Sector_.resize(x);
		data_.clear();
		offsets_.clear();
		nzMsAndQns_.clear();
	}

	template<typename SomeBasisType>
	void set(const typename PsimagLite::Vector<VectorType>::Type& v,
	         const SomeBasisType& someBasis)
	{
		size_ = someBasis.size();
		nzMsAndQns_.clear();
		data_.clear();
		data_.resize(v.size());
		offsets_.resize(v.size()+1);
		for (SizeType i=0;i<v.size();i++) {
			data_[i] = v[i];
			offsets_[i] = someBasis.partition(i);
			if (v[i].size()>0) {
				SizeType qn = someBasis.pseudoEffectiveNumber(offsets_[i]);
				nzMsAndQns_.push_back(PairSizeType(i, qn));
			}
		}

		offsets_[v.size()]=size_;
		setIndex2Sector();
	}

	template<typename SomeBasisType>
	void populateSectors(const SomeBasisType& someBasis)
	{
		SizeType np = someBasis.partition()-1;
		size_ = someBasis.size();
		nzMsAndQns_.clear();
		data_.clear();
		data_.resize(np);
		offsets_.resize(np+1);
		for (SizeType i=0;i<np;i++) {
			offsets_[i] = someBasis.partition(i);
			SizeType total = someBasis.partition(i+1)-offsets_[i];
			VectorType tmpV(total,0);
			data_[i] = tmpV;
			SizeType qn = someBasis.pseudoEffectiveNumber(offsets_[i]);
			nzMsAndQns_.push_back(PairSizeType(i, qn));
		}

		offsets_[np]=size_;
		setIndex2Sector();
		PsimagLite::OstringStream msg;
		msg<<"Populated "<<np<<" sectors";
		progress_.printline(msg,std::cout);
	}

	template<typename SomeBasisType>
	void populateFromQns(const typename PsimagLite::Vector<SizeType>::Type& qns,
	                     const SomeBasisType& someBasis)
	{
		SizeType np = someBasis.partition()-1;
		size_ = someBasis.size();
		nzMsAndQns_.clear();
		data_.clear();
		data_.resize(np);
		offsets_.resize(np+1);
		for (SizeType i=0;i<np;i++) {
			offsets_[i] = someBasis.partition(i);
		}
		offsets_[np]=size_;

		for (SizeType i=0;i<qns.size();i++) {
			SizeType ip = findPartitionWithThisQn(qns[i],someBasis);
			SizeType total = someBasis.partition(ip+1)-offsets_[ip];
			VectorType tmpV(total,0);
			data_[ip] = tmpV;
			nzMsAndQns_.push_back(PairSizeType(ip, qns[i]));
		}

		setIndex2Sector();
		PsimagLite::OstringStream msg;
		msg<<"populateFromQns "<<qns.size()<<" sectors";
		progress_.printline(msg,std::cout);
	}

	void collapseSectors()
	{
		SizeType np = data_.size();

		typename PsimagLite::Vector<PairSizeType>::Type nzMsAndQns;
		for (SizeType i=0;i<np;i++) {
			if (isZero(data_[i]))
				data_[i].resize(0);
			else
				nzMsAndQns.push_back(nzMsAndQns_[i]);
		}

		nzMsAndQns_ = nzMsAndQns;
		setIndex2Sector();
		PsimagLite::OstringStream msg;
		msg<<"Collapsed. Non-zero sectors now are "<<nzMsAndQns_.size();
		progress_.printline(msg,std::cout);
	}

	void setDataInSector(const VectorType& v,SizeType i0)
	{
		data_[i0] = v;
	}

	SizeType sectors() const { return nzMsAndQns_.size(); }

	SizeType sector(SizeType i) const
	{
		assert(i < nzMsAndQns_.size());
		return nzMsAndQns_[i];
	}

	template<typename SomeBasisType>
	void fromFull(const VectorType& v,const SomeBasisType& someBasis)
	{
		size_ = someBasis.size();

		offsets_.resize(someBasis.partition());
		for (SizeType i=0;i<someBasis.partition();i++)
			offsets_[i] = someBasis.partition(i);
		assert(offsets_[offsets_.size()-1]==size_);

		data_.clear();
		data_.resize(someBasis.partition()-1);

		nzMsAndQns_.clear();
		findPartitions(nzMsAndQns_,v,someBasis);
		for (SizeType jj=0;jj<nzMsAndQns_.size();jj++) {
			SizeType j = nzMsAndQns_[jj];
			//firstSector_ = j;
			SizeType offset = offsets_[j];
			SizeType total = offsets_[j+1]-offset;
			data_[j].resize(total);
			for (SizeType i=0;i<total;i++) data_[j][i] = v[i+offset];
		}
		setIndex2Sector();
	}

	void extract(VectorType& v,SizeType i) const
	{
		v=data_[i];
	}

	SizeType size() const { return size_; }

	SizeType effectiveSize(SizeType i) const { return data_[i].size(); }

	SizeType offset(SizeType i) const { return offsets_[i]; }

	const ComplexOrRealType& fastAccess(SizeType i,SizeType j) const
	{
		assert(i < data_.size());
		assert(j < data_[i].size());
		return data_[i][j];
	}

	ComplexOrRealType& fastAccess(SizeType i,SizeType j)
	{
		assert(i < data_.size());
		assert(j < data_[i].size());
		return data_[i][j];
	}

	const ComplexOrRealType& slowAccess(SizeType i) const
	{
		assert(i<index2Sector_.size());
		int j = index2Sector_[i];
		if (j<0) return zero_;
		return data_[j][i-offsets_[j]];
	}

	ComplexOrRealType& slowAccess(SizeType i)
	{
		int j = index2Sector_[i];
		if (j<0) {
			PsimagLite::String msg("VectorWithOffsets");
			std::cerr<<msg<<" can't build itself dynamically yet (sorry!)\n";
			return data_[0][0];
		}

		return data_[j][i-offsets_[j]];
	}

	template<typename SparseVectorType>
	void toSparse(SparseVectorType& sv) const
	{
		sv.resize(size_);
		for (SizeType jj=0;jj<nzMsAndQns_.size();jj++) {
			SizeType j =  nzMsAndQns_[jj].first;
			for (SizeType i=0;i<data_[j].size();i++)
				sv[i+offsets_[j]] = data_[j][i];
		}
	}

	template<typename IoOutputter>
	void save(IoOutputter& io,const PsimagLite::String& label) const
	{
		io.printline(label);
		PsimagLite::String s="#size="+ttos(size_);
		io.printline(s);
		io.printVector(offsets_,"#offsets");
		s = "#nonzero="+ttos(nzMsAndQns_.size());
		io.printline(s);

		for (SizeType jj=0;jj<nzMsAndQns_.size();jj++) {
			SizeType j =  nzMsAndQns_[jj];
			s="#sector="+ttos(j);
			io.printline(s);
			io.printVector(data_[j],s);
		}
	}

	template<typename IoInputter>
	void load(IoInputter& io,const PsimagLite::String& label,SizeType counter=0)
	{
		PsimagLite::String msg("VectorWithOffsets:");
		io.advance(label,counter);
		int x = 0;
		io.readline(x,"#size=");
		if (x<0)
			throw PsimagLite::RuntimeError(msg + ":load(...): size<0\n");
		size_ = x;
		io.read(offsets_,"#offsets");
		data_.clear();
		data_.resize(offsets_.size());
		io.readline(x,"#nonzero=");
		if (x<0)
			throw PsimagLite::RuntimeError(msg + ":load(...): nonzerosectors<0\n");
		nzMsAndQns_.resize(x);
		for (SizeType jj=0;jj<nzMsAndQns_.size();jj++) {
			io.readline(x,"#sector=");
			if (x<0)
				throw PsimagLite::RuntimeError(msg + ":load(...): sector<0\n");
			if (SizeType(x)>=data_.size())
				throw PsimagLite::RuntimeError(msg + ":load(...): sector too big\n");
			int y = 0;
			io.readline(y, "#qn=");
			if (y < 0)
				throw PsimagLite::RuntimeError(msg + ":load(...): qn<0\n");
			nzMsAndQns_[jj] = PairSizeType(x, y);
			io.read(data_[x],"#sector=");
		}

		setIndex2Sector();
	}

	// We don't have a partitioned basis because we don't
	// have the superblock basis at this point
	// Therefore, partitioning is bogus here
	template<typename IoInputter>
	void loadOneSector(IoInputter& io,
	                   const PsimagLite::String& label,
	                   SizeType counter=0)
	{
		PsimagLite::String msg("VectorWithOffsets:");
		io.advance(label,counter);
		int x = 0;
		io.readline(x,"#size=");
		if (x<0)
			throw PsimagLite::RuntimeError(msg + ":loadOneSector(...): size<0\n");
		size_ = x;

		io.read(offsets_,"#offsets");

		data_.clear();
		data_.resize(offsets_.size());

		io.readline(x,"#nonzero=");
		if (x < 0)
			throw PsimagLite::RuntimeError(msg + ":loadOneSector(...): nonzerosectors<0\n");
		nzMsAndQns_.resize(x);

		for (SizeType jj=0;jj<nzMsAndQns_.size();jj++) {
			io.readline(x,"#sector=");
			if (x<0)
				throw PsimagLite::RuntimeError(msg + ":loadOneSector(...): sector<0\n");
			if (SizeType(x)>=data_.size())
				throw PsimagLite::RuntimeError(msg + ":loadOneSector(...): sector too big\n");
			int y = 0;
			io.readline(y, "#qn=");
			if (y < 0)
				throw PsimagLite::RuntimeError(msg + ":load(...): qn<0\n");
			nzMsAndQns_[jj] = PairSizeType(x, y);
			io.read(data_[x],"#sector=");
		}

		setIndex2Sector();
	}

	VectorWithOffsets<ComplexOrRealType> operator+=(const VectorWithOffsets<ComplexOrRealType>& v)
	{
		if (nzMsAndQns_.size()==0) {
			size_ = v.size_;
			data_ = v.data_;
			offsets_ = v.offsets_;
			nzMsAndQns_ = v.nzMsAndQns_;
			setIndex2Sector();
			return *this;
		}

		for (SizeType ii=0;ii<nzMsAndQns_.size();ii++) {
			SizeType i = nzMsAndQns_[ii].first;
			data_[i] += v.data_[i];
		}

		setIndex2Sector();
		return *this;
	}

	int index2Sector(SizeType i) const
	{
		assert(i < index2Sector_.size());
		return index2Sector_[i];
	}

	template<typename ComplexOrRealType2>
	friend ComplexOrRealType2 norm(const Dmrg::VectorWithOffsets<ComplexOrRealType2>& v);

	template<typename ComplexOrRealType2>
	friend ComplexOrRealType2 norm(const Dmrg::VectorWithOffsets<std::complex<ComplexOrRealType2> >& v);

	template<typename ComplexOrRealType2>
	friend void normalize(Dmrg::VectorWithOffsets<std::complex<ComplexOrRealType2> >& v);

	template<typename ComplexOrRealType3,typename ComplexOrRealType2>
	friend VectorWithOffsets<ComplexOrRealType2> operator*(const ComplexOrRealType3&,
	                                                       const VectorWithOffsets<ComplexOrRealType2>&);

	template<typename ComplexOrRealType2>
	friend VectorWithOffsets<ComplexOrRealType2> operator+(const VectorWithOffsets<ComplexOrRealType2>&,
	                                                       const VectorWithOffsets<ComplexOrRealType2>&);

private:

	void setIndex2Sector()
	{
		if (index2Sector_.size()!=size_) index2Sector_.resize(size_);
		for (SizeType i=0;i<size_;i++) {
			index2Sector_[i] = -1;
			for (SizeType jj=0;jj<nzMsAndQns_.size();jj++) {
				SizeType j = nzMsAndQns_[jj].first;
				if (i<offsets_[j] || i>=offsets_[j+1]) continue;
				index2Sector_[i] = j;
			}
		}
	}

	template<typename SomeBasisType>
	void findPartitions(typename PsimagLite::Vector<SizeType>::Type& p,
	                    const VectorType& v,
	                    const SomeBasisType& someBasis)
	{
		bool found = false;
		p.clear();
		for (SizeType i=0;i<someBasis.partition()-1;i++) {
			if (nonZeroPartition(v,someBasis,i)) {
				found = true;
				p.push_back(i);
			}
		}
		if (!found) {
			PsimagLite::OstringStream msg;
			msg<<"No partition found";
			progress_.printline(msg,std::cout);
		}

	}

	template<typename SomeBasisType>
	bool nonZeroPartition(const VectorType& v,
	                      const SomeBasisType& someBasis,SizeType i)
	{
		typename VectorType::value_type zero = 0;
		for (SizeType j=someBasis.partition(i);j<someBasis.partition(i+1);j++) {
			if (v[j]!=zero) return true;
		}
		return false;
	}

	bool isZero(const VectorType& v) const
	{
		RealType eps = 1e-5;
		for (SizeType i=0;i<v.size();i++)
			if (fabs(PsimagLite::real(v[i]))>eps || fabs(PsimagLite::imag(v[i]))>eps)
				return false;
		return true;
	}

	template<typename SomeBasisType>
	SizeType findPartitionWithThisQn(SizeType qn,
	                                 const SomeBasisType& someBasis) const
	{
		SizeType np = someBasis.partition()-1;
		for (SizeType i=0;i<np;i++) {
			SizeType state = someBasis.partition(i);
			if (SizeType(someBasis.qn(state))==qn) return i;
		}
		throw PsimagLite::RuntimeError("findPartitionWithThisQn\n");
	}

	PsimagLite::ProgressIndicator progress_;
	SizeType size_;
	typename PsimagLite::Vector<int>::Type index2Sector_;
	typename PsimagLite::Vector<VectorType>::Type data_;
	typename PsimagLite::Vector<SizeType>::Type offsets_;
	typename PsimagLite::Vector<PairSizeType>::Type nzMsAndQns_;
}; // class VectorWithOffset

template<typename ComplexOrRealType>
inline ComplexOrRealType norm(const Dmrg::VectorWithOffsets<ComplexOrRealType>& v)
{
	ComplexOrRealType sum=0;
	for (SizeType ii=0;ii<v.nzMsAndQns_.size();ii++) {
		SizeType i = v.nzMsAndQns_[ii];
		ComplexOrRealType tmp = PsimagLite::norm(v.data_[i]);
		sum += tmp*tmp;
	}
	return sqrt(sum);
}

template<typename ComplexOrRealType>
inline ComplexOrRealType norm(const Dmrg::VectorWithOffsets<std::complex<ComplexOrRealType> >& v)
{
	ComplexOrRealType sum=0;
	for (SizeType ii=0;ii<v.nzMsAndQns_.size();ii++) {
		SizeType i = v.nzMsAndQns_[ii];
		ComplexOrRealType tmp = PsimagLite::norm(v.data_[i]);
		sum += tmp*tmp;
	}
	return sqrt(sum);
}

template<typename ComplexOrRealType>
inline void normalize(Dmrg::VectorWithOffsets<std::complex<ComplexOrRealType> >& v)
{
	ComplexOrRealType norma = PsimagLite::norm(v);
	ComplexOrRealType eps = 1e-5;

	if (fabs(norma-1.0)<eps) return;

	PsimagLite::String s(__FILE__);
	s += " " + ttos(__LINE__);
	std::cerr<<s<<" norm= "<<norma<<"\n";
	assert(fabs(norma)>eps);

	for (SizeType i=0;i<v.data_.size();i++)
		for (SizeType j=0;j<v.data_[i].size();j++)
			v.data_[i][j] /= norma;
}

template<typename ComplexOrRealType>
inline ComplexOrRealType operator*(const Dmrg::VectorWithOffsets<ComplexOrRealType>& v1,
                                   const Dmrg::VectorWithOffsets<ComplexOrRealType>& v2)
{
	ComplexOrRealType sum = 0;
	for (SizeType ii=0;ii<v1.sectors();ii++) {
		SizeType i = v1.sector(ii);
		for (SizeType jj=0;jj<v1.sectors();jj++) {
			SizeType j = v2.sector(jj);
			if (i!=j) continue;
			for (SizeType k=0;k<v1.effectiveSize(i);k++)
				sum+= v1.fastAccess(i,k)*PsimagLite::conj(v2.fastAccess(j,k));
		}
	}
	return sum;
}

template<typename ComplexOrRealType,typename ComplexOrRealType2>
inline VectorWithOffsets<ComplexOrRealType2> operator*(const ComplexOrRealType& value,
                                                       const VectorWithOffsets<ComplexOrRealType2>& v)
{
	VectorWithOffsets<ComplexOrRealType2> w = v;

	for (SizeType ii=0;ii<w.nzMsAndQns_.size();ii++) {
		SizeType i = w.nzMsAndQns_[ii];
		w.data_[i] *= value;
	}
	return w;
}

template<typename ComplexOrRealType>
VectorWithOffsets<ComplexOrRealType> operator+(const VectorWithOffsets<ComplexOrRealType>& v1,
                                               const VectorWithOffsets<ComplexOrRealType>& v2)
{
	PsimagLite::String s = "VectorWithOffsets + VectorWithOffsets failed\n";
	if (v1.nzMsAndQns_!=v2.nzMsAndQns_)
		throw PsimagLite::RuntimeError(s.c_str());
	for (SizeType ii=0;ii<v1.nzMsAndQns_.size();ii++) {
		SizeType i = v1.nzMsAndQns_[ii];
		if (v1.data_[i].size()!=v2.data_[i].size())
			throw PsimagLite::RuntimeError(s.c_str());
	}
	VectorWithOffsets<ComplexOrRealType> w = v1;
	w += v2;
	return w;
}

template<typename ComplexOrRealType>
const ComplexOrRealType VectorWithOffsets<ComplexOrRealType>::zero_ = 0;
}
/*@}*/
#endif

