// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009, UT-Battelle, LLC
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

/*! \file HamiltonianSymmetrySu2.h
	*
	*  An object of this class is meant to be contained by a Hilbert Space basis
	*  and then this class help creating the outer product states necessary to preserve the SU(2) symmetry
	*  THe counterpart to this class is HamiltonianSymmetryLocal.h for local symmetries only
	*
	*/
#ifndef HAM_SYMM_SU2_H
#define HAM_SYMM_SU2_H

#include "Sort.h" // in PsimagLite
#include "JmPairs.h"
#include "VerySparseMatrix.h"
#include "JmSubspace.h"
#include "ProgramGlobals.h"
#include "CrsMatrix.h"

namespace Dmrg {

	template<typename RealType,typename SparseMatrixType>
	class	HamiltonianSymmetrySu2 {
		public:
			typedef std::pair<size_t,size_t> PairType;
		private:
			typedef JmPairs<PairType> JmPairsType;
			typedef VerySparseMatrix<RealType> VerySparseMatrixType;
			typedef HamiltonianSymmetrySu2<RealType,SparseMatrixType> ThisType;
			typedef JmSubspace<RealType,VerySparseMatrixType,ThisType> JmSubspaceType;
			typedef typename JmSubspaceType::FlavorType FlavorType;
			typedef  BasisData<PairType> BasisDataType;
			
		public:
			typedef PsimagLite::CrsMatrix<RealType> FactorsType;
			
			static size_t const MAX = ProgramGlobals::MaxNumberOfSites;
			
			PairType jmValue(size_t i) const { return jmValues_[i]; }

			void set(const BasisDataType& basisData)
			{
				jmValues_=basisData.jmValues;
				flavors_=basisData.flavors;
				flavorsMax_= *(std::max_element(
						flavors_.begin(),flavors_.end()));
					//flavors_[utils::vectorMax<size_t,std::greater<size_t> >(flavors_,0)];

				size_t electronsMax1 = *(std::max_element(
						basisData.electronsUp.begin(),basisData.electronsUp.end()));
					// basisData.electronsUp[utils::vectorMax<size_t,std::greater<size_t>
					//	>(basisData.electronsUp,0)];

				size_t electronsMax2 = *(std::max_element(
						basisData.electronsDown.begin(),basisData.electronsDown.end()));
					// basisData.electronsDown[utils::vectorMax<size_t,std::greater<size_t>
					//	>(basisData.electronsDown,0)];
				electronsMax_=electronsMax1+electronsMax2;
				jMax_=0;
				jmValues_.maxFirst<std::greater<size_t> >(jMax_);
				jMax_++;
				calcReducedBasis();
			}

			static void findQuantumNumbers(std::vector<size_t> &q,const BasisDataType& basisData) 
			{
				q.resize(basisData.electronsUp.size());
				for (size_t i=0;i<q.size();i++) {
					size_t ne = basisData.electronsUp[i]+basisData.electronsDown[i];
					PairType jmpair = basisData.jmValues[i];
					q[i]=neJmToIndex(ne,jmpair);
				}
			}

			static size_t neJmToIndex(size_t ne,const PairType& jm) 
			{
				std::vector<size_t> v(3);
				v[2]=0;
				v[2]=jm.first;
				v[0]=jm.second;
				double m = jm.second+0.5*(ne-jm.first);
				if (m<0 || m>65535) throw std::runtime_error(" neJmToIndex\n");
				v[0]=size_t(m);
				if (ne<v[0]) throw std::runtime_error(" neJmToIndex 2\n");
				v[1]=ne-v[0];
				return encodeQuantumNumber(v);
			}

			static size_t encodeQuantumNumber(const std::vector<size_t>& v)
			{
				size_t x= v[0] + v[1]*MAX;
				if (v[0]>=MAX || v[1]>=MAX || v[2]>=MAX) throw std::runtime_error("encodeQuantumNumber\n");
				if (v.size()==3) x += v[2]*MAX*MAX;
				return x;
			}

			static std::vector<size_t> decodeQuantumNumber(size_t q)
			{
				std::vector<size_t> v(3);
				v[2] = size_t(q/(MAX*MAX));
				size_t tmp = q - v[2]*MAX*MAX;
				v[1] = size_t(tmp/MAX);
				v[0] = tmp % MAX;
				return v;
			}

			//! targets[0]=nup, targets[1]=ndown,  targets[2]=2j
			static size_t pseudoQuantumNumber(const std::vector<size_t>& v)
			{
				size_t x= (v[0] + v[1]);
				x += v[2]*2*MAX;
				return x;
			}

			void setToProduct(
					const HamiltonianSymmetrySu2& symm1,
					const HamiltonianSymmetrySu2& symm2,
					int pseudoQn,
					const std::vector<size_t>& electrons1,
     					const std::vector<size_t>& electrons2,
					std::vector<size_t>& electrons,
					std::vector<size_t>& quantumNumbers)
			{
				size_t ns = symm1.jmValues_.size();
				size_t ne = symm2.jmValues_.size();

				JmSubspaceType::setToProduct(&symm1,&symm2,electrons1,electrons2);

				findAllowedJm(symm1,symm2,electrons1,electrons2,pseudoQn);;
				createFactors(ns,ne);
				setFlavors(quantumNumbers);
				if (quantumNumbers.size()!=(ns*ne)) {
					std::cerr<<"ns="<<ns<<" ne="<<ne<<" but quantumNumbers.size="<<quantumNumbers.size()<<"\n";
					throw std::runtime_error("Wrong number of quantum numbers\n");
				}
			
				jMax_=0;
				jmValues_.maxFirst<std::greater<size_t> >(jMax_);
				jMax_++;
				calcReducedBasis();
				normalizeFlavors();
				setElectrons(electrons,quantumNumbers);
				electronsMax_ = *(std::max_element(
						electrons.begin(),electrons.end()));
		//				electrons[utils::vectorMax<size_t,std::greater<size_t> >(electrons,0)];
			}

			size_t pseudoEffectiveNumber(size_t nelectrons,size_t jtilde) const
			{
				std::vector<size_t> v(3);
				v[0]=nelectrons;
				v[1]=0;
				v[2]=jtilde;
				return pseudoQuantumNumber(v);
			}

			PairType getJmValue(size_t alpha) const
			{
				return jmValues_[alpha];
			}

			size_t getFlavor(size_t alpha) const
			{
				return flavors_[alpha];
			}

			size_t flavorsMax() const { return flavorsMax_; }

			size_t electronsMax() const { return electronsMax_; }

			size_t jMax() const { return jMax_; }

			template<typename SolverParametersType>
			void calcRemovedIndices(std::vector<size_t>& removedIndices,std::vector<RealType>& eigs,size_t kept,
					const SolverParametersType& solverParams)
			{
				normalizeFlavors();

				// we sort the eigenvalues
				// note: eigenvalues are not ordered because DensityMatrix is diagonalized in blocks
				std::vector<size_t> perm(eigs.size());
				Sort<std::vector<RealType> > sort;
				sort.sort(eigs,perm);
				
				if (eigs.size()<=kept) return;
				size_t target = eigs.size()-kept;

				removedIndices.clear();

				if (solverParams.options.find("inflate")!=std::string::npos)
					inclusiveRemoval(removedIndices,perm,eigs,target);
				else exclusiveRemoval(removedIndices,perm,eigs,target);
			}

			const FactorsType& getFactors() const 
			{
				return factors_;
			}

			void reorder(const std::vector<size_t>& permutationVector)
			{
				// reorder jmValues
				jmValues_.reorder(permutationVector);

				// reorder flavors
				utils::reorder(flavors_,permutationVector);
				utils::reorder(flavorsOld_,permutationVector);
			}

			void truncate(std::vector<size_t> const &removedIndices,const std::vector<size_t>& electrons)
			{
				electronsMax_= * (std::max_element(
						electrons.begin(),electrons.end()));
				//electrons[utils::vectorMax<size_t,std::greater<size_t> >(electrons,0)];

				utils::truncateVector(flavors_,removedIndices);
				flavorsMax_=* (std::max_element(
						flavors_.begin(),flavors_.end()));
				//flavors_[utils::vectorMax<size_t,std::greater<size_t> >(flavors_,0)];

				jmValues_.truncate(removedIndices);
				jMax_=0;
				jmValues_.maxFirst<std::greater<size_t> >(jMax_);
				jMax_++;
				calcReducedBasis();
			}

			size_t size() const {return jmValues_.size(); }

			template<typename IoInputter>
			void load(IoInputter& io) 
			{
				jmValues_.load(io); 
				io.read(flavors_,"#su2flavors");
				io.readline(flavorsMax_,"#su2FlavorsMax=");
				io.readline(electronsMax_,"#su2ElectronsMax=");
				io.readline(jMax_,"#su2JMax=");
				io.read(statesReduced_,"#su2StatesReduced");
				io.read(jvals_,"#su2Jvals");
			}

			template<typename IoOutputter>
			void save(IoOutputter& io) const
			{
				jmValues_.save(io);
				io.printVector(flavors_,"#su2flavors");
				std::string s="#su2FlavorsMax=" + ttos(flavorsMax_)+"\n";
				io.print(s); 
				s="#su2ElectronsMax="+ttos(electronsMax_)+"\n";
				io.print(s); 
				s="#su2JMax="+ttos(jMax_)+"\n";
				io.print(s);
				io.printVector(statesReduced_,"#su2StatesReduced");
				io.printVector(jvals_,"#su2Jvals");
			}

			size_t flavor2Index(size_t f1,size_t f2,size_t ne1,size_t ne2,size_t j1,size_t j2) const
			{
				return JmSubspaceType::flavor(f1,f2,ne1,ne2,j1,j2);
				
			}

			void flavor2Index(std::map<size_t,size_t>& flavorsOldInverse, const PairType& jm) const
			{
				for (size_t i=0;i<flavorsOld_.size();i++) {
					if (jmValues_[i]!=jm) continue;
					flavorsOldInverse[flavorsOld_[i]]=i;
				}
				
			}

			const std::vector<size_t>& flavorsOld() const
			{
				return flavorsOld_;
				
			}

			// reduced:
			size_t reducedIndex(size_t i) const { return statesReduced_[i]; }

			size_t reducedSize() const { return statesReduced_.size(); }

			size_t jVals(size_t i) const { return jvals_[i]; }

			size_t jVals() const { return jvals_.size(); }

		private:

			template<typename JmSubspaceType>
			size_t  setFlavors(std::vector<size_t>& quantumNumbers,JmSubspaceType& jmSubspace,size_t offset)
			{
				// order is important here, electrons must be set after quantumNumbers
				size_t flavors = jmSubspace.numberOfFlavors();
				if (offset==0) {
					quantumNumbers.clear();
					jmValues_.clear();
					flavors_.clear();
				}
				for (size_t i=0;i<flavors;i++ ) {
					PairType jm = jmSubspace.getJmValue();
					quantumNumbers.push_back(neJmToIndex(jmSubspace.getNe(),jm));
					jmValues_.push(jm,i+offset);
					flavors_.push_back(jmSubspace.getFlavor(i));
				}
				offset += flavors;

				return offset;
			}

			void setElectrons(std::vector<size_t>& electrons,const std::vector<size_t>& qns)
			{
				electrons.resize(qns.size());
				for (size_t i=0;i<qns.size();i++) {
					std::vector<size_t> v = decodeQuantumNumber(qns[i]);
					electrons[i]=v[0]+v[1];
				}
			}

			void normalizeFlavors()
			{
				flavorsOld_=flavors_;
				std::vector<size_t> perm(flavors_.size());
				Sort<std::vector<size_t> > sort;
				sort.sort(flavors_,perm);

				size_t counter=0;
				size_t flavorSaved=flavors_[0];
				std::vector<size_t> flavorsTmp(flavors_.size());

				for (size_t i=0;i<flavors_.size();i++) {
					if (flavorSaved!=flavors_[i]) {
						counter++;
						flavorSaved=flavors_[i];
					}
					flavorsTmp[i]=counter;
				}

				for (size_t i=0;i<flavors_.size();i++) 
					flavors_[perm[i]]=flavorsTmp[i];

				flavorsMax_=counter+1;

			}

			void setFlavors(std::vector<size_t>& quantumNumbers) 
			{
				size_t offset=0;
				for (size_t i=0;i<jmSubspaces_.size();i++) {
					offset = setFlavors(quantumNumbers,jmSubspaces_[i],offset);
					jmSubspaces_[i].clear();
				}
			}

			// note: j is actually 2j and m is actually m+j
			// note: this is so that j and m are both always size_t
			void findAllowedJm(const ThisType& symm1,const ThisType& symm2,const std::vector<size_t>& electrons1,
					   const std::vector<size_t>& electrons2,int pseudoQn) 
			{
				size_t ns = symm1.jmValues_.size();
				size_t ne = symm2.jmValues_.size();

				jmSubspaces_.clear();
				for (size_t i=0;i<ns;i++) {
					PairType jm1 = symm1.getJmValue(i);
					for (size_t j=0;j<ne;j++) {
						PairType jm2 = symm2.getJmValue(j);
						size_t nelectrons = electrons1[i]+electrons2[j];
						addAllowedJms(jm1,jm2,i,j,ns,nelectrons,pseudoQn);
					}
				}
			}
			
			void createFactors(size_t ns,size_t ne)
			{
				VerySparseMatrixType factors(ns*ne);
				size_t offset=0;
				for (size_t i=0;i<jmSubspaces_.size();i++) {
					size_t s=0;
					s= jmSubspaces_[i].createFactors(factors,offset);
					offset += s;
				}
				if (factors.nonZero()==0) {
					for (size_t i=0;i<jmSubspaces_.size();i++) {
						std::cerr<<"subspace number "<<i;
						size_t nelectrons=jmSubspaces_[i].getNe();
						std::cerr<<" nelectrons="<<nelectrons;
						PairType jm = jmSubspaces_[i].getJmValue();
						std::cerr<<" pseudo="<<pseudoEffectiveNumber(nelectrons,jm.first);
						
						std::cerr<<" jm=("<<jm.first<<","<<jm.second<<")\n";
						std::cerr<<" heavy="<<jmSubspaces_[i].heavy()<<"\n";
						std::cerr<<"--------------------------------------------\n";
					}
					throw std::runtime_error("HSSU2.h::createFactors(): factors are empty\n");
				}
				factors.sort();
				factors_ = factors;
			}

			// note: j is actually 2j and m is actually m+j
			// note: this is so that j and m are both always size_t
			void addAllowedJms(const PairType& jm1,const PairType& jm2,size_t alpha,size_t beta,size_t ns,size_t nelectrons,
					int pseudoQn)
			{
				int j1 = jm1.first, j2=jm2.first;
				int jinitial = j1-j2;
				if (jinitial<0) jinitial = -jinitial;
				for (int j=jinitial;j<=j1+j2;j++) {
					// go over all hurdles:
					// first hurdle (j1+j2+j is even)
					//if (alpha==6 && beta==0) std::cerr<<__LINE__<<"\n";
					if ((j1+j2+j) %2 !=0) continue;
					// calculate m
					// m = m1+m2+(j1+j2-j)/2
					// note: (j1+j2-j) is even
					// note: (j1+j2-j)>=0
					int m = jm1.second+ jm2.second + int((-j1-j2+j)/2);
					// second hurdle |2m-j| <= j
					int tmp = 2*m-j;
					if (tmp<0) tmp= -tmp;
					if (tmp>j) continue;
					PairType jm(j,m);
					int heavy=1;
					if (pseudoQn>=0 &&  pseudoEffectiveNumber(nelectrons,jm.first) !=size_t(pseudoQn)) {
						heavy=0;
					}
					addJmPair(alpha+beta*ns,jm1,jm2,jm,nelectrons,heavy);
				}
			}

			void addJmPair(size_t index,const PairType& jm1,const PairType& jm2,const PairType& jm,size_t nelectrons,int heavy)
			{
				std::pair<PairType,size_t> triplet;
				triplet.first=jm;
				triplet.second=nelectrons;
				int x = PsimagLite::isInVector(jmSubspaces_,triplet);
				if (x<0) { // add new jmSubspace
					JmSubspace<RealType,VerySparseMatrixType,ThisType> jmSubspace(jm,index,jm1,jm2,nelectrons,heavy);
					jmSubspaces_.push_back(jmSubspace);
				} else {
					jmSubspaces_[x].push(index,jm1,jm2,nelectrons);
				}
			}

			void inclusiveRemoval(std::vector<size_t>& removedIndices,const std::vector<size_t>& perm,
					      const std::vector<RealType>& eigs,size_t target)
			{
				std::vector<size_t> permInverse(perm.size());
				for (size_t i=0;i<permInverse.size();i++) permInverse[perm[i]]=i;

				for (size_t i=0;i<target;i++) {
					if (PsimagLite::isInVector(removedIndices,perm[i])>=0) continue;
					removedIndices.push_back(perm[i]);
				}

				for (size_t i=0;i<target;i++)	{
					for (size_t j=0;j<eigs.size();j++) {
						
						if (flavors_[j]==flavors_[perm[i]] && jmValues_[j].first==jmValues_[perm[i]].first) {
							int x = PsimagLite::isInVector(removedIndices,j);
							if (x<0) {
								std::vector<size_t>::iterator p1 =
										find(removedIndices.begin(),removedIndices.end(),perm[i]);
								if (p1==removedIndices.end()) continue;
								removedIndices.erase(p1);
								if (fabs(eigs[permInverse[j]]-eigs[i])>1e-6) {
									std::cerr<<"ind="<<perm[i]<<" j="<<permInverse[j];
									std::cerr<<" e[ind]="<<eigs[i];
									std::cerr<<" e[j]="<<eigs[permInverse[j]]<<"\n";
									std::cerr<<"flavor="<<flavors_[j];
									std::cerr<<" jm=("<<jmValues_[j].first<<",";
									std::cerr<<jmValues_[j].second<<") ";
									std::cerr<<" jm[ind]=(";
									std::cerr<<jmValues_[perm[i]].first<<",";
									std::cerr<<jmValues_[perm[i]].second<<")\n";
								}
								break;
							}
						}
					}
				}
			}

			void exclusiveRemoval(std::vector<size_t>& removedIndices,const std::vector<size_t>& perm,
					      const std::vector<RealType>& eigs,size_t target)
			{
				std::vector<size_t> permInverse(perm.size());
				for (size_t i=0;i<permInverse.size();i++) permInverse[perm[i]]=i;

				for (size_t i=0;i<target;i++) {
					if (removedIndices.size()>=target) break;
					if (PsimagLite::isInVector(removedIndices,perm[i])>=0) continue;
						removedIndices.push_back(perm[i]);

					for (size_t j=0;j<eigs.size();j++) {

						if (flavors_[j]==flavors_[perm[i]] && jmValues_[j].first==jmValues_[perm[i]].first) {
							int x = PsimagLite::isInVector(removedIndices,j);
							if (x<0) {
								removedIndices.push_back(j);
								if (fabs(eigs[permInverse[j]]-eigs[i])>1e-6) {
									std::cerr<<"ind="<<perm[i]<<" j="<<permInverse[j];
									std::cerr<<" e[ind]="<<eigs[i];
									std::cerr<<" e[j]="<<eigs[permInverse[j]]<<"\n";
									std::cerr<<"flavor="<<flavors_[j]<<" jm=";
									std::cerr<<jmValues_[j].first<<" ";
									std::cerr<<jmValues_[j].second<<" ";
									std::cerr<<" jm[ind]=";
									std::cerr<<jmValues_[perm[i]].first<<" ";
									std::cerr<<jmValues_[perm[i]].second<<"\n";
								}
							}
						}
					}
				}
			}

			void calcReducedBasis()
			{
				jvals_.clear();
				statesReduced_.clear();
				for (size_t i1=0;i1<jmValues_.size();i1++) {
					PairType jm1 = jmValues_[i1];
					if (jm1.first!=jm1.second) continue;
					statesReduced_.push_back(i1);
					int x = PsimagLite::isInVector(jvals_,jm1.first);
					if (x<0) jvals_.push_back(jm1.first);
				}
			}

			JmPairsType jmValues_;
			std::vector<size_t> flavors_,flavorsOld_;
			size_t flavorsMax_,electronsMax_,jMax_;
			FactorsType factors_;
			std::vector<JmSubspaceType> jmSubspaces_;
			// reduced:
			std::vector<size_t> statesReduced_;
			std::vector<size_t> jvals_;
	}; //class HamiltonianSymmetrySu2
} // namespace Dmrg

/*@}*/
#endif
