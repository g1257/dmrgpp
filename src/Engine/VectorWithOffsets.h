// BEGIN LICENSE BLOCK
/*
Copyright © 2009 , UT-Battelle, LLC
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
// END LICENSE BLOCK
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
#include "Utils.h"

namespace Dmrg {
	template<typename FieldType>
	class VectorWithOffsets {
			typedef VectorWithOffsets<FieldType> ThisType;
		public:
			typedef FieldType value_type;
			typedef std::pair<size_t,size_t> PairType;
			typedef std::vector<FieldType> VectorType;
			
			VectorWithOffsets() : progress_("VectorWithOffsets",0),zero_(0) { }
			
			template<typename SomeBasisType>
			VectorWithOffsets(const std::vector<size_t>& weights,
					 const SomeBasisType& someBasis)
				: 
				progress_("VectorWithOffsets",0),
				zero_(0),
				size_(someBasis.size()),
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
			}
			
			template<typename SomeBasisType>
			void set(const std::vector<VectorType>& v,//const std::vector<size_t>& weights,
				 const SomeBasisType& someBasis)
			{
				size_ = someBasis.size();
				nonzeroSectors_.clear();
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
			}
			
			template<typename SomeBasisType>
			void populateSectors(const SomeBasisType& someBasis)
			{
				size_t np = someBasis.partition()-1;
				size_ = someBasis.size();
				nonzeroSectors_.clear();
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
				if (offsets_[offsets_.size()-1]!=size_) throw std::runtime_error
					("TST::fromFull(...): internal error\n");
				
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
			}
			
			void extract(VectorType& v,size_t i) const
			{
				v=data_[i];
			}
			
			void print(std::ostream& os,const std::string& label) const
			{
				os<<label<<"\n";
				os<<size_<<"\n";
				for (size_t jj=0;jj<nonzeroSectors_.size();jj++) {
					size_t j =  nonzeroSectors_[jj];
					size_t offset = offsets_[j];
					os<<data_[j].size()<<"\n";
					for (size_t i=0;i<data_[j].size();i++) 
						os<<(i+offset)<<" "<<data_[j][i]<<"\n";
				}
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
				/*if (nonzeroSectors_.size()==1) {
					if (i>=offsets_[firstSector_] && i<offsets_[firstSector_+1])
						return data_[firstSector_][i-offsets_[firstSector_]];
					else return zero_;
				}*/
				
				for (size_t jj=0;jj<nonzeroSectors_.size();jj++) {
					size_t j = nonzeroSectors_[jj];
					if (i>=offsets_[j] && i<offsets_[j+1]) {
						return data_[j][i-offsets_[j]];
					}
				}
				return zero_;
			}
			
			FieldType& operator[](size_t i) //__attribute__((always_inline))
			{
				for (size_t jj=0;jj<nonzeroSectors_.size();jj++) {
					size_t j = nonzeroSectors_[jj];
					if (i>=offsets_[j] && i<offsets_[j+1]) {
						return data_[j][i-offsets_[j]];
					}
				}
				throw std::runtime_error("VectorWithOffsets\n");
			}
			
			ThisType& operator= (const ThisType& f)
			{
				size_=f.size_;
				data_=f.data_;
				offsets_=f.offsets_;
				nonzeroSectors_=f.nonzeroSectors_;
				//firstSector_ = f.firstSector_;
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
			friend void normalize(Dmrg::VectorWithOffsets<std::complex<FieldType2> >& v);
		
		private:
			template<typename SomeBasisType>
			void findPartitions(std::vector<size_t>& p,const VectorType& v,const SomeBasisType& someBasis)
			{
				bool found = false;
				p.clear();
				for (size_t i=0;i<someBasis.partition()-1;i++) {
					if (nonZeroPartition(v,someBasis,i)) {
						if (found) {
							std::ostringstream msg;
							msg<<"More than one partition found";
							progress_.printline(msg,std::cout);
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
					if (fabs(real(v[i]))>eps || fabs(imag(v[i]))>eps) return false;
				return true; 
			}
			
			ProgressIndicator progress_;
			const FieldType zero_;
			size_t size_;
			std::vector<VectorType> data_;
			std::vector<size_t> offsets_;
			std::vector<size_t> nonzeroSectors_;
			//size_t firstSector_;
	}; // class VectorWithOffset
	
	/*template<typename FieldType>
	std::ostream& operator<<(std::ostream& os,const VectorWithOffsets<FieldType>& s)
	{
		s.print(os,"VectorWithOffsets");
		return os;
	}*/
	
} // namespace Dmrg

namespace std {
	template<typename FieldType>
	inline FieldType norm(const Dmrg::VectorWithOffsets<FieldType>& v)
	{
		FieldType sum=0;
		for (size_t ii=0;ii<v.nonzeroSectors_.size();ii++) {
			size_t i = v.nonzeroSectors_[ii];
			sum += std::norm(v.data_[i]);
		}
		return sum;
	}
	
	template<typename FieldType>
	inline FieldType norm(const Dmrg::VectorWithOffsets<std::complex<FieldType> >& v)
	{
		FieldType sum=0;
		for (size_t ii=0;ii<v.nonzeroSectors_.size();ii++) {
			size_t i = v.nonzeroSectors_[ii];
			sum += std::norm(v.data_[i]);
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
		if (fabs(norma-1.0)<eps) {
			//std::cerr<<"VectorWithOffsets::normalize():";
			//std::cerr<<"norm is already one, nothing to do\n";
			return;
		} 
		//std::cerr<<"norm="<<norma<<"\n";
		if (fabs(norma)<eps) throw std::runtime_error("Too small!\n");
		for (size_t i=0;i<v.data_.size();i++) for (size_t j=0;j<v.data_[i].size();j++) 
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
}
/*@}*/
#endif

