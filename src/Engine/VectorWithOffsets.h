/*
Copyright (c) 2009-2012, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
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

// FIXME: a more generic solution is needed instead of tying the non-zero structure to basis
namespace Dmrg {
	template<typename FieldType>
	class VectorWithOffsets {
		typedef VectorWithOffsets<FieldType> ThisType;
		static FieldType const zero_;
	public:
		typedef FieldType value_type;
		typedef std::pair<size_t,size_t> PairType;
		typedef std::vector<FieldType> VectorType;
			
		VectorWithOffsets() 
		: progress_("VectorWithOffsets",0),size_(0),index2Sector_(0)
		{ }
			
		template<typename SomeBasisType>
		VectorWithOffsets(const std::vector<size_t>& weights,
		                  const SomeBasisType& someBasis)
		: progress_("VectorWithOffsets",0),
		  size_(someBasis.size()),
		  index2Sector_(size_),
		  data_(weights.size()),
		  offsets_(weights.size()+1)
		{
			for (size_t i=0;i<weights.size();i++) {
				data_[i].resize(weights[i]);
				offsets_[i] = someBasis.partition(i);
				if (weights[i]>0) {
					nonzeroSectors_.push_back(i);
					//firstSector_ = i;
				}
			}
			offsets_[weights.size()]=size_;
			setIndex2Sector();
		}

		void resize(size_t x)
		{
			size_ = x;
			index2Sector_.resize(x);
			data_.clear();
			offsets_.clear();
			nonzeroSectors_.clear();
		}

		template<typename SomeBasisType>
		void set(const std::vector<VectorType>& v,//const std::vector<size_t>& weights,
		         const SomeBasisType& someBasis)
		{
			size_ = someBasis.size();
			nonzeroSectors_.clear();
			data_.clear();
			data_.resize(v.size());
			offsets_.resize(v.size()+1);
			for (size_t i=0;i<v.size();i++) {
				data_[i] = v[i];
				offsets_[i] = someBasis.partition(i);
				if (v[i].size()>0) {
					nonzeroSectors_.push_back(i);
					//firstSector_ = i;
				}
				//if (v[i].size()!=0 && weights[i]==0) throw std::runtime_error("VectorWithOffsets::"
				//			"set(...)\n");
			}
			offsets_[v.size()]=size_;
			setIndex2Sector();
		}

		template<typename SomeBasisType>
		void populateSectors(const SomeBasisType& someBasis)
		{
			size_t np = someBasis.partition()-1;
			size_ = someBasis.size();
			nonzeroSectors_.clear();
			data_.clear();
			data_.resize(np);
			offsets_.resize(np+1);
			for (size_t i=0;i<np;i++) {
				offsets_[i] = someBasis.partition(i);
				size_t total = someBasis.partition(i+1)-offsets_[i];
				VectorType tmpV(total,0);
				data_[i] = tmpV;
				nonzeroSectors_.push_back(i);
			}
			offsets_[np]=size_;
			setIndex2Sector();
			std::ostringstream msg;
			msg<<"Populated "<<np<<" sectors";
			progress_.printline(msg,std::cout);
		}

		void collapseSectors()
		{
			size_t np = data_.size();
			nonzeroSectors_.clear();
			for (size_t i=0;i<np;i++) {
				if (isZero(data_[i])) {
					data_[i].resize(0);
				} else {
					nonzeroSectors_.push_back(i);
				}
			}
			setIndex2Sector();
			std::ostringstream msg;
			msg<<"Collapsed. Non-zero sectors now are "<<nonzeroSectors_.size();
			progress_.printline(msg,std::cout);
		}

		void setDataInSector(const VectorType& v,size_t i0)
		{
			data_[i0] = v;
		}

		size_t sectors() const { return nonzeroSectors_.size(); }

		size_t sector(size_t i) const { return nonzeroSectors_[i]; }

		template<typename SomeBasisType>
		void fromFull(const VectorType& v,const SomeBasisType& someBasis)
		{
			size_ = someBasis.size();
			
			offsets_.resize(someBasis.partition());
			for (size_t i=0;i<someBasis.partition();i++)
				offsets_[i] = someBasis.partition(i);
			assert(offsets_[offsets_.size()-1]==size_);
			
			data_.clear();
			data_.resize(someBasis.partition()-1);
			
			nonzeroSectors_.clear();
			findPartitions(nonzeroSectors_,v,someBasis);
			for (size_t jj=0;jj<nonzeroSectors_.size();jj++) {
				size_t j = nonzeroSectors_[jj];
				//firstSector_ = j;
				size_t offset = offsets_[j];
				size_t total = offsets_[j+1]-offset;
				data_[j].resize(total);
				for (size_t i=0;i<total;i++) data_[j][i] = v[i+offset];
			}
			setIndex2Sector();
		}

		void extract(VectorType& v,size_t i) const
		{
			v=data_[i];
		}
		
		size_t size() const { return size_; }
		
		size_t effectiveSize(size_t i) const { return data_[i].size(); }
		 
		size_t offset(size_t i) const { return offsets_[i]; }
		
		const FieldType& fastAccess(size_t i,size_t j) const 
		{
			return data_[i][j];
		}
				
		const FieldType& operator[](size_t i) const //__attribute__((always_inline))
		{
			assert(i<index2Sector_.size());
			int j = index2Sector_[i];
			if (j<0) return zero_;
			return data_[j][i-offsets_[j]];
			/*if (nonzeroSectors_.size()==1) {
				if (i>=offsets_[firstSector_] && i<offsets_[firstSector_+1])
					return data_[firstSector_][i-offsets_[firstSector_]];
				else return zero_;
			}*/
			
// 			for (size_t jj=0;jj<nonzeroSectors_.size();jj++) {
// 				size_t j = nonzeroSectors_[jj];
// 				if (i<offsets_[j] || i>=offsets_[j+1]) continue;
// 				return data_[j][i-offsets_[j]];
// 			}
// 			return zero_;

//				size_t x = nonzeroSectors_.size()/2;
//				size_t j = 0;
//				while(true) {
//					j = nonzeroSectors_[x];
//					if (i<offsets_[j]) {
//						if (x==0) return zero_;
//						x--;
//						continue;
//					}
//					if (i>=offsets_[j+1]) {
//						x++;
//						if (x==nonzeroSectors_.size()) return zero_;
//						continue;
//					}
//					break;
//				}
//				return data_[j][i-offsets_[j]];
		}
		
		FieldType& operator[](size_t i) //__attribute__((always_inline))
		{
// 			for (size_t jj=0;jj<nonzeroSectors_.size();jj++) {
// 				size_t j = nonzeroSectors_[jj];
// 				if (i>=offsets_[j] && i<offsets_[j+1]) {
// 					return data_[j][i-offsets_[j]];
// 				}
// 			}
			int j = index2Sector_[i];
			if (j<0) {
				std::cerr<<"VectorWithOffsets can't build itself dynamically yet (sorry!)\n";
				return data_[0][0];
			}
			return data_[j][i-offsets_[j]];
			//throw std::runtime_error("VectorWithOffsets can't build itself dynamically yet (sorry!)\n");
		}
		
		/*ThisType& operator= (const ThisType& f)
		{
			size_=f.size_;
			data_=f.data_;
			offsets_=f.offsets_;
			nonzeroSectors_=f.nonzeroSectors_;
			return *this;
		}*/

		template<typename SparseVectorType>
		void toSparse(SparseVectorType& sv) const
		{
			sv.resize(size_);
			for (size_t jj=0;jj<nonzeroSectors_.size();jj++) {
				size_t j =  nonzeroSectors_[jj];
				for (size_t i=0;i<data_[j].size();i++)
					sv[i+offsets_[j]] = data_[j][i];
			}
		}

		template<typename IoOutputter>
		void save(IoOutputter& io,const std::string& label) const
		{
			io.printline(label);
			std::string s="#size="+ttos(size_);
			io.printline(s);

//			io.print(label);
//			io.print("#size=",size_);
			io.printVector(offsets_,"#offsets");
			s = "#nonzero="+ttos(nonzeroSectors_.size());
			io.printline(s);
//			io.print("#nonzero=",nonzeroSectors_.size());

			for (size_t jj=0;jj<nonzeroSectors_.size();jj++) {
				size_t j =  nonzeroSectors_[jj];
				s="#sector="+ttos(j);
				io.printline(s);
//				std::string s="#sector="+ttos(j);
//				io.print("#sector",j);
				io.printVector(data_[j],s);
			}
		}
		
		template<typename IoInputter>
		void load(IoInputter& io,const std::string& label,size_t counter=0)
		{
			io.advance(label,counter);
			int x = 0;
			io.readline(x,"#size=");
			if (x<0) throw std::runtime_error("VectorWithOffsets::load(...): size<0\n");
			size_ = x;
			io.read(offsets_,"#offsets");
			data_.clear();
			data_.resize(offsets_.size());
			io.readline(x,"#nonzero=");
			if (x<0) throw std::runtime_error("VectorWithOffsets::load(...): nonzerosectors<0\n");
			nonzeroSectors_.resize(x);
			for (size_t jj=0;jj<nonzeroSectors_.size();jj++) {
				io.readline(x,"#sector=");
				if (x<0) 
					throw std::runtime_error("VectorWithOffsets::load(...): sector<0\n");
				if (size_t(x)>=data_.size()) 
					throw std::runtime_error("VectorWithOffsets::load(...): sector too big\n");
				nonzeroSectors_[jj] = x;
				io.read(data_[x],"#sector=");
			}
			setIndex2Sector();
		}

		VectorWithOffsets<FieldType> operator+=(const VectorWithOffsets<FieldType>& v)
		{
			if (nonzeroSectors_.size()==0) {
				size_ = v.size_;
				data_ = v.data_;
				offsets_ = v.offsets_;
				nonzeroSectors_ = v.nonzeroSectors_;
				setIndex2Sector();
				return *this;
			}
			for (size_t ii=0;ii<nonzeroSectors_.size();ii++) {
				size_t i = nonzeroSectors_[ii];
				data_[i] += v.data_[i];
			}
			setIndex2Sector();
			return *this;
		}
		
		template<typename FieldType2>
		friend FieldType2 std::norm(const Dmrg::VectorWithOffsets<FieldType2>& v);

		template<typename FieldType2>
		friend FieldType2 std::norm(const Dmrg::VectorWithOffsets<std::complex<FieldType2> >& v);
		
		template<typename FieldType2>
		friend std::complex<FieldType2> multiply(const Dmrg::VectorWithOffsets<std::complex<FieldType2> >& v1,
		                                         const Dmrg::VectorWithOffsets<std::complex<FieldType2> >& v2);

		template<typename FieldType2>
		friend FieldType2 multiply(const Dmrg::VectorWithOffsets<FieldType2>& v1,
		                           const Dmrg::VectorWithOffsets<FieldType2>& v2);

		template<typename FieldType2>
		friend void normalize(Dmrg::VectorWithOffsets<std::complex<FieldType2> >& v);
		
		template<typename FieldType3,typename FieldType2>
		friend VectorWithOffsets<FieldType2> operator*(const FieldType3& value,const VectorWithOffsets<FieldType2>& v);
	
		template<typename FieldType2>
		friend VectorWithOffsets<FieldType2> operator+(const VectorWithOffsets<FieldType2>& v1,
							       const VectorWithOffsets<FieldType2>& v2);

	private:
		
		void setIndex2Sector()
		{
			if (index2Sector_.size()!=size_) index2Sector_.resize(size_);
			//assert(size_>0);
			for (size_t i=0;i<size_;i++) {
				index2Sector_[i] = -1;
				for (size_t jj=0;jj<nonzeroSectors_.size();jj++) {
					size_t j = nonzeroSectors_[jj];
					if (i<offsets_[j] || i>=offsets_[j+1]) continue;
					index2Sector_[i] = j;
				}
			}
		}

		template<typename SomeBasisType>
		void findPartitions(std::vector<size_t>& p,const VectorType& v,const SomeBasisType& someBasis)
		{
			bool found = false;
			p.clear();
			for (size_t i=0;i<someBasis.partition()-1;i++) {
				if (nonZeroPartition(v,someBasis,i)) {
					if (found) {
//						std::ostringstream msg;
//						msg<<"More than one partition found";
//						progress_.printline(msg,std::cout);
						//throw std::runtime_error("quiting\n");
					}
					found = true;
					p.push_back(i);
				}
			}
			if (!found) {
				std::ostringstream msg;
				msg<<"No partition found";
				progress_.printline(msg,std::cout);
				//p.push_back(0);
			}

		}
		
		template<typename SomeBasisType>
		bool nonZeroPartition(const VectorType& v,const SomeBasisType& someBasis,size_t i)
		{
			typename VectorType::value_type zero = 0;
			for (size_t j=someBasis.partition(i);j<someBasis.partition(i+1);j++) {
				if (v[j]!=zero) return true;
			}
			return false;
		}
		
		bool isZero(const VectorType& v) const
		{
			double eps = 1e-5;
			for (size_t i=0;i<v.size();i++)
				if (fabs(std::real(v[i]))>eps || fabs(std::imag(v[i]))>eps) return false;
			return true; 
		}

		PsimagLite::ProgressIndicator progress_;
		size_t size_;
		std::vector<int> index2Sector_;
		std::vector<VectorType> data_;
		std::vector<size_t> offsets_;
		std::vector<size_t> nonzeroSectors_;
	}; // class VectorWithOffset

// 	template<typename FieldType>
// 	std::ostream& operator<<(std::ostream& os,const VectorWithOffsets<FieldType>& s)
// 	{
// 		os<<size_<<"\n";
// 		os<<data_;
// 		os<<offsets_;
// 		os<<nonzeroSectors_;
// 		return os;
// 	}
} // namespace Dmrg

namespace std {
	template<typename FieldType>
	inline FieldType norm(const Dmrg::VectorWithOffsets<FieldType>& v)
	{
		FieldType sum=0;
		for (size_t ii=0;ii<v.nonzeroSectors_.size();ii++) {
			size_t i = v.nonzeroSectors_[ii];
			sum += PsimagLite::norm(v.data_[i]);
		}
		return sum;
	}
	
	template<typename FieldType>
	inline FieldType norm(const Dmrg::VectorWithOffsets<std::complex<FieldType> >& v)
	{
		FieldType sum=0;
		for (size_t ii=0;ii<v.nonzeroSectors_.size();ii++) {
			size_t i = v.nonzeroSectors_[ii];
			sum += PsimagLite::norm(v.data_[i]);
		}
		return sum;
	}

}

namespace Dmrg {
	template<typename FieldType>
	inline void normalize(Dmrg::VectorWithOffsets<std::complex<FieldType> >& v)
	{
		FieldType norma = std::norm(v);
		FieldType eps = 1e-5;

		if (fabs(norma-1.0)<eps) return;

		std::string s(__FILE__);
		s += " " + ttos(__LINE__);
		std::cerr<<s<<" norm= "<<norma<<"\n";
		assert(fabs(norma)>eps);

		for (size_t i=0;i<v.data_.size();i++)
			for (size_t j=0;j<v.data_[i].size();j++) 
				v.data_[i][j] /= norma;
	}

	template<typename FieldType>
	inline std::complex<FieldType> multiply(const Dmrg::VectorWithOffsets<std::complex<FieldType> >& v1,
	                                        const Dmrg::VectorWithOffsets<std::complex<FieldType> >& v2)
	{
		std::complex<FieldType> sum=0;
		for (size_t ii=0;ii<v1.nonzeroSectors_.size();ii++) {
			size_t i = v1.nonzeroSectors_[ii];
			sum += v1.data_[i]*v2.data_[i]; // call to * will conj()
		}
		return sum;
	}

	template<typename FieldType>
	inline FieldType multiply(const Dmrg::VectorWithOffsets<FieldType>& v1,
	                          const Dmrg::VectorWithOffsets<FieldType>& v2)
	{
		FieldType sum=0;
		for (size_t ii=0;ii<v1.nonzeroSectors_.size();ii++) {
			size_t i = v1.nonzeroSectors_[ii];
			sum += v1.data_[i]*v2.data_[i]; // call to * will conj()
		}
		return sum;
	}

	//! Isn't this function equal to the prev.? need to check FIXME
	template<typename FieldType>
	inline std::complex<FieldType> operator*(const Dmrg::VectorWithOffsets<std::complex<FieldType> >& v1,
	                                         const Dmrg::VectorWithOffsets<std::complex<FieldType> >& v2)
	{
		std::complex<FieldType> sum = 0;
		for (size_t ii=0;ii<v1.sectors();ii++) {
			size_t i = v1.sector(ii);
			for (size_t jj=0;jj<v1.sectors();jj++) {
				size_t j = v2.sector(jj);
				if (i!=j) continue; //throw std::runtime_error("Not same sector\n");
				size_t offset = v1.offset(i);
				for (size_t k=0;k<v1.effectiveSize(i);k++) 
					sum+= v1[k+offset] * conj(v2[k+offset]);
			}
		}
		return sum;
	}

	template<typename FieldType>
	inline FieldType operator*(const Dmrg::VectorWithOffsets<FieldType>& v1,
	                           const Dmrg::VectorWithOffsets<FieldType>& v2)
	{
		FieldType sum = 0;
		for (size_t ii=0;ii<v1.sectors();ii++) {
			size_t i = v1.sector(ii);
			for (size_t jj=0;jj<v1.sectors();jj++) {
				size_t j = v2.sector(jj);
				if (i!=j) continue; //throw std::runtime_error("Not same sector\n");
				size_t offset = v1.offset(i);
				for (size_t k=0;k<v1.effectiveSize(i);k++)
					sum+= v1[k+offset] * std::conj(v2[k+offset]);
			}
		}
		return sum;
	}

	template<typename FieldType,typename FieldType2>
	inline VectorWithOffsets<FieldType2> operator*(const FieldType& value,const VectorWithOffsets<FieldType2>& v)
	{
		VectorWithOffsets<FieldType2> w = v;

		for (size_t ii=0;ii<w.nonzeroSectors_.size();ii++) {
			size_t i = w.nonzeroSectors_[ii];
			w.data_[i] *= value;
		}
		return w;
	}

	template<typename FieldType>
	VectorWithOffsets<FieldType> operator+(const VectorWithOffsets<FieldType>& v1,
					       const VectorWithOffsets<FieldType>& v2)
	{
		std::string s = "VectorWithOffsets + VectorWithOffsets failed\n";
		if (v1.nonzeroSectors_!=v2.nonzeroSectors_) throw std::runtime_error(s.c_str());
		for (size_t ii=0;ii<v1.nonzeroSectors_.size();ii++) {
			size_t i = v1.nonzeroSectors_[ii];
			if (v1.data_[i].size()!=v2.data_[i].size())
				throw std::runtime_error(s.c_str());
		}
		VectorWithOffsets<FieldType> w = v1;
		w += v2;
		return w;
	}


	template<typename FieldType>
	const FieldType VectorWithOffsets<FieldType>::zero_ = 0;
}
/*@}*/
#endif

