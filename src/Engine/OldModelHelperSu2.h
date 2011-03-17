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
		
namespace Dmrg {	
		/*const SparseMatrixType& getTcOperator(int i,size_t type) const
		{
			if (type==System) return basis2tc_[i];
			return basis3tc_[i];
		}*/	
	
		const SparseMatrixType& getReducedOperator(char modifier,size_t i,size_t sigma,size_t type) const
		{
			size_t i = i*dof_+sigma;
			if (modifier=='N') {
				if (type==System) return basis2_.getOperatorByIndex(ii).data;
				return basis3_.getOperatorByIndex(ii).data;
			} else {
				return getTcOperator(ii,type);
			} 
		}
		
	//! //! Does matrixBlock= (AB), A belongs to pSprime and B  belongs to pEprime or viceversa (inter)
		void fastOpProdInter(SparseMatrixTemplate<MatrixElementType> const &A,
				SparseMatrixTemplate<MatrixElementType> const &B,
				int type,
				MatrixElementType  &hop,
				SparseMatrixTemplate<MatrixElementType> &matrixBlock,
				bool operatorsAreFermions=true) const
		{
			int const SystemEnviron=1,EnvironSystem=2;
			
			int fermionSign =  (operatorsAreFermions) ? -1 : 1;
			//! work only on partition m
		
			if (type==EnvironSystem)  {
				MatrixElementType hop2 =hop*fermionSign;
				fastOpProdInter(B,A,SystemEnviron,hop2,matrixBlock,operatorsAreFermions);
				return;
			}		
		
			//! work only on partition m
			int m = m_;
			//std::cerr<<"About to allocated block in outer: bs="<<bs<<"\n";
			int offset = basis1_.partition(m);
			//std::cerr<<"doing partition (outer)"<<m<<" with size="<<bs<<"\n";
			int total = basis1_.partition(m+1) - offset;
			int ns = basis2_.size();
			
			int counter=0;
			matrixBlock.resize(total);
			
			const SparseMatrixType& factors = basis1_.getFactors();
			//const SparseMatrixType& factorsInverse = basis1_.getFactorsInverse();

			
			for (int i=0;i<total;i++) {
				// row i of the ordered product basis
				matrixBlock.setRow(i,counter);
				//size_t ii = basis1_.permutation(i+offset);
				std::vector<int> columns;
				std::vector<MatrixElementType> values;
				//for (int k=factorsInverse_.getRowPtr(ii);k<factorsInverse_.getRowPtr(ii+1);k++) {
				for (size_t k=gamma_[i];k<gamma_[i+1];k++) {
// 					int alpha=alpha_[i][k-factorsInverse_.getRowPtr(ii)];
// 					int beta=beta_[i][k-factorsInverse_.getRowPtr(ii)];
					//utils::getCoordinates(alpha,beta,factorsInverse.getCol(k),ns);
					int alpha=alpha_[k];
					int beta=beta_[k];
					//MatrixElementType tmp2=factorsInverse_.getValue(k)*hop;
					MatrixElementType tmp2=factorsInverseValue_[k]*hop;
					if (operatorsAreFermions) tmp2 *= basis2_.fermionicSign(alpha,fermionSign);
					
					for (int k2=A.getRowPtr(alpha);k2<A.getRowPtr(alpha+1);k2++) {
						int alphaPrime = A.getCol(k2);
						MatrixElementType tmp33=tmp2*A.getValue(k2);
						for (int k3=B.getRowPtr(beta);k3<B.getRowPtr(beta+1);k3++) {
							int betaPrime= B.getCol(k3);
							MatrixElementType tmp = tmp33* B.getValue(k3);
							/* fermion signs note:
							here the environ is applied first and has to "cross"
							the system, hence the sign factor pSprime.fermionicSign(alpha)
							*/
							
							size_t jj = alphaPrime + betaPrime*ns;
							
							for (int k4=factors.getRowPtr(jj);k4<factors.getRowPtr(jj+1);k4++) {
								//if (fabs(tmp*factors.getValue(k4))<1e-8) continue;
								int j = basis1_.permutationInverse(factors.getCol(k4)) -offset;
								if (j<0 || j>=total) continue;
								int x = PsimagLite\:\:isInVector(columns,j);
								if (x<0) {
									columns.push_back(j);
									values.push_back(tmp*factors.getValue(k4));
								} else {
									values[x] += tmp*factors.getValue(k4);
								}
							}
						}
					}
				}
				if (columns.size()==0) continue;
				//std::vector<size_t> perm(columns.size());
				//utils::sort<int,double>(columns,perm);
				for (size_t jj=0;jj<columns.size();jj++) {
					matrixBlock.pushCol(columns[jj]);
					matrixBlock.pushValue(values[jj]);
					counter++;
				}
			}
			matrixBlock.setRow(total,counter);
		}
	
		
		//! Does x+= (AB)y, where A belongs to pSprime and B  belongs to pEprime or viceversa (inter)
		//! Has been changed to accomodate for reflection symmetry
		 void fastOpProdInter(	std::vector<MatrixElementType>  &x,
					std::vector<MatrixElementType>  const &y,
					SparseMatrixTemplate<MatrixElementType> const &A,
					SparseMatrixTemplate<MatrixElementType> const &B,
					int type,
					MatrixElementType  &hop,
					bool operatorsAreFermions,size_t angularMomentum,MatrixElementType angularFactor,size_t category,bool flip=false)  const 
		{
			int const SystemEnviron=1,EnvironSystem=2;
			int fermionSign =  (operatorsAreFermions) ? -1 : 1;
			
			if (type==EnvironSystem)  {
				MatrixElementType hop2 =hop*fermionSign;
				fastOpProdInter(x,y,B,A,SystemEnviron,hop2,operatorsAreFermions,angularMomentum,angularFactor,category);
				return;
			}
			
			//! work only on partition m
			int m = m_;
			//std::cerr<<"About to allocated block in outer: bs="<<bs<<"\n";
			int offset = basis1_.partition(m);
			//std::cerr<<"doing partition (outer)"<<m<<" with size="<<bs<<"\n";
			int total = basis1_.partition(m+1) - offset;
			int ns = basis2_.size();
			
			
			int j;
			MatrixElementType tmp,tmp2,tmp3;
			
			for (int i=0;i<total;i++) {
				
				// row i of the ordered product basis
				//utils::getCoordinates(alpha,beta,modelHelper.basis1().permutation(i+offset),ns);
				//size_t ii = basis1_.permutation(i+offset);
				//for (int k=factorsInverse_.getRowPtr(ii);k<factorsInverse_.getRowPtr(ii+1);k++) {
				for (size_t k=gamma_[i];k<gamma_[i+1];k++) {
					int alpha=alpha_[k];
					int beta=beta_[k];
					
					//alpha=alpha_[i][k-factorsInverse_.getRowPtr(ii)];
					//beta=beta_[i][k-factorsInverse_.getRowPtr(ii)];
					
					//utils::getCoordinates(alpha,beta,factorsInverse_.getCol(k),ns);
					//tmp =factorsInverse_.getValue(k)*hop;
					tmp =factorsInverseValue_[k]*hop;
					if (operatorsAreFermions && basis2_.fermionicSign(alpha,fermionSign)<0) tmp = -tmp;		
					for (int k2=A.getRowPtr(alpha);k2<A.getRowPtr(alpha+1);k2++) {
						int alphaPrime = A.getCol(k2);
						tmp2 =tmp*A.getValue(k2);
						for (int k3=B.getRowPtr(beta);k3<B.getRowPtr(beta+1);k3++) {
							int betaPrime= B.getCol(k3);
							tmp3 =tmp2 * B.getValue(k3);
							
							size_t jj = alphaPrime + betaPrime*ns;
							
							for (int k4=factors_.getRowPtr(jj);k4<factors_.getRowPtr(jj+1);k4++) {
						
								/* fermion signs note:
						   		here the environ is applied first and has to "cross"
						   		the system, hence the sign factor pSprime.fermionicSign(alpha,tmp)
						 		*/
								//int j = basis1_.permutationInverse(factors_.getCol(k4)) -offset;
								j = jCached_[jj][k4-factors_.getRowPtr(jj)];
								if (j<0 || j>=total) continue;
								x[i] += y[j]*tmp3*factors_.getValue(k4);
							}
						}
					}
				}
			}
			
		}
		
		//! Let H_{alpha,beta; alpha',beta'} = basis2.hamiltonian_{alpha,alpha'} \delta_{beta,beta'}
		//! Let H_m be  the m-th block (in the ordering of basis1) of H
		//! Then, this function does x += H_m * y
		//! This is a performance critical function
		//! Has been changed to accomodate for reflection symmetry
		void hamiltonianLeftProduct(std::vector<MatrixElementType> &x,std::vector<MatrixElementType> const &y) const 
		{ 
			int m = m_;
			int offset = basis1_.partition(m);
			int bs = basis1_.partition(m+1)-offset;
			SparseMatrixType hamiltonian = basis2_.hamiltonian();
			int ns = basis2_.size();
				
			const SparseMatrixType& factors = basis1_.getFactors();
			//const SparseMatrixType& factorsInverse = basis1_.getFactorsInverse();

			
			for (int i=0;i<bs;i++) {
				//size_t ii = basis1_.permutation(i+offset);
				//for (int k=factorsInverse_.getRowPtr(ii);k<factorsInverse_.getRowPtr(ii+1);k++) {
				for (size_t k=gamma_[i];k<gamma_[i+1];k++) {
					int r=alpha_[k];
					int beta=beta_[k];
					//int r,beta;
					//utils::getCoordinates(r,beta,factorsInverse_.getCol(k),ns);
					
					for (int kk=hamiltonian.getRowPtr(r);kk<hamiltonian.getRowPtr(r+1);kk++) {
						size_t alphaPrime = hamiltonian.getCol(kk);
						size_t jj = alphaPrime + beta*ns;
						MatrixElementType value =   hamiltonian.getValue(kk)*factorsInverseValue_[k];
						for (int k3=factors.getRowPtr(jj);k3<factors.getRowPtr(jj+1);k3++) {
							int j = basis1_.permutationInverse(factors.getCol(k3))-offset;
							if (j<0 || j>=bs) continue;
							x[i] += y[j]*value*factors.getValue(k3);
						}
					}
				}
			}
				
			
		}		
		
		//! Let  H_{alpha,beta; alpha',beta'} = basis2.hamiltonian_{beta,beta'} \delta_{alpha,alpha'}
		//! Let H_m be  the m-th block (in the ordering of basis1) of H
		//! Then, this function does x += H_m * y
		//! This is a performance critical function
		void hamiltonianRightProduct(std::vector<MatrixElementType> &x,std::vector<MatrixElementType> const &y) const 
		{ 
			int m = m_;
			int offset = basis1_.partition(m);
			int bs = basis1_.partition(m+1)-offset;
			SparseMatrixType hamiltonian = basis3_.hamiltonian();
			int ns = basis2_.size();
				
			const SparseMatrixType& factors = basis1_.getFactors();
			//const SparseMatrixType& factorsInverse = basis1_.getFactorsInverse();

			
			for (int i=0;i<bs;i++) {
				//size_t ii = basis1_.permutation(i+offset);
				//for (int k=factorsInverse_.getRowPtr(ii);k<factorsInverse_.getRowPtr(ii+1);k++) {
				for (size_t k=gamma_[i];k<gamma_[i+1];k++) {
					int r=alpha_[k];
					int beta=beta_[k];
					//int r,beta;
					//utils::getCoordinates(r,beta,factorsInverse_.getCol(k),ns);
					
					for (int kk=hamiltonian.getRowPtr(beta);kk<hamiltonian.getRowPtr(beta+1);kk++) {
						size_t betaPrime = hamiltonian.getCol(kk);
						size_t jj = r + betaPrime*ns;
						MatrixElementType value =   hamiltonian.getValue(kk)*factorsInverseValue_[k];
						for (int k3=factors.getRowPtr(jj);k3<factors.getRowPtr(jj+1);k3++) {
							int j = basis1_.permutationInverse(factors.getCol(k3))-offset;
							if (j<0 || j>=bs) continue;
							x[i] += y[j]*value*factors.getValue(k3);
						}
					}
				}
			}
		}
		
		//! Note: USed only for debugging
		void calcHamiltonianPartLeft(SparseMatrixType &matrixBlock) const
		{
			int m  = m_;
			int offset = basis1_.partition(m);
			
			int bs = basis1_.partition(m+1)-offset;
			const SparseMatrixType& hamiltonian = basis2_.hamiltonian();
			int ns = basis2_.size();
			//if (!isHermitian(hamiltonian)) throw std::runtime_error("Not hermitian hamiltonian (left)\n");
			matrixBlock.resize(bs);
			
			const SparseMatrixType& factors = basis1_.getFactors();
			//const SparseMatrixType& factorsInverse = basis1_.getFactorsInverse();

			
			int counter=0;
			for (int i=0;i<bs;i++) {
				matrixBlock.setRow(i,counter);
				//size_t ii = basis1_.permutation(i+offset);
				std::vector<int> columns;
				std::vector<MatrixElementType> values;
				//for (int k=factorsInverse_.getRowPtr(ii);k<factorsInverse_.getRowPtr(ii+1);k++) {
				for (size_t k=gamma_[i];k<gamma_[i+1];k++) {
					int alpha=alpha_[k];
					int beta=beta_[k];	
					//int alpha,beta;
					//utils::getCoordinates(alpha,beta,factorsInverse_.getCol(k),ns);
					for (int k2=hamiltonian.getRowPtr(alpha);k2<hamiltonian.getRowPtr(alpha+1);k2++) {
						MatrixElementType tmp = factorsInverseValue_[k]*hamiltonian.getValue(k2);
						size_t jj = hamiltonian.getCol(k2)+beta*ns;
						for (int k3=factors.getRowPtr(jj);k3<factors.getRowPtr(jj+1);k3++) {
							int j = basis1_.permutationInverse(factors.getCol(k3))-offset;
							if (j<0 || j>=bs) continue;
							int x = PsimagLite\:\:isInVector(columns,j);
							if (x<0) {
								columns.push_back(j);
								values.push_back(tmp*factors.getValue(k3));
							} else {
								values[x] += 	tmp*factors.getValue(k3);
							}
						}
					}
				}
				
				for (size_t jj=0;jj<columns.size();jj++) {
					matrixBlock.pushCol(columns[jj]);
					matrixBlock.pushValue(values[jj]);
					counter++;
				}
			}
			
			matrixBlock.setRow(bs,counter);
			//if (!isHermitian(matrixBlock)) throw std::runtime_error("Not hermitian matrixBlock (left)\n");
			
		}
		
		//! Note: USed only for debugging
		void calcHamiltonianPartRight(SparseMatrixType &matrixBlock) const
		{
			int m  = m_;
			int offset = basis1_.partition(m);
			
			int bs = basis1_.partition(m+1)-offset;
			const SparseMatrixType& hamiltonian= basis3_.hamiltonian();
			int ns = basis2_.size();
			
			//if (!isHermitian(hamiltonian)) throw std::runtime_error("Not hermitian hamiltonian (right)\n");
			
			//std::cerr<<"Adding matrix block of size: "<<ns<<" for option="<<option<<"\n";
			matrixBlock.resize(bs);

			const SparseMatrixType& factors = basis1_.getFactors();
			//const SparseMatrixType& factorsInverse = basis1_.getFactorsInverse();

			
			int counter=0;
			for (int i=0;i<bs;i++) {
				matrixBlock.setRow(i,counter);
				//size_t ii = basis1_.permutation(i+offset);
				std::vector<int> columns;
				std::vector<MatrixElementType> values;
				//for (int k=factorsInverse_.getRowPtr(ii);k<factorsInverse_.getRowPtr(ii+1);k++) {
				for (size_t k=gamma_[i];k<gamma_[i+1];k++) {
					int alpha=alpha_[k];
					int beta=beta_[k];
					//int alpha,beta;
					//utils::getCoordinates(alpha,beta,factorsInverse_.getCol(k),ns);
					for (int k2=hamiltonian.getRowPtr(beta);k2<hamiltonian.getRowPtr(beta+1);k2++) {
						MatrixElementType tmp = factorsInverseValue_[k]*hamiltonian.getValue(k2);
						size_t jj = alpha + hamiltonian.getCol(k2)*ns;
						for (int k3=factors.getRowPtr(jj);k3<factors.getRowPtr(jj+1);k3++) {
							int j = basis1_.permutationInverse(factors.getCol(k3))-offset;
							if (j<0 || j>=bs) continue;
							int x = PsimagLite\:\:isInVector(columns,j);
							if (x<0) {
								columns.push_back(j);
								values.push_back(tmp*factors.getValue(k3));
							} else {
								values[x] += 	tmp*factors.getValue(k3);
							}
						}
					}
				}
				
				for (size_t jj=0;jj<columns.size();jj++) {
					matrixBlock.pushCol(columns[jj]);
					matrixBlock.pushValue(values[jj]);
					counter++;
				}
			}
			matrixBlock.setRow(bs,counter);
			//if (!isHermitian(matrixBlock)) throw std::runtime_error("Not hermitian matrixBlock (right)\n");
		}
		
				
// 		std::vector<MatrixElementType> factorsInverseValue_;
// 		std::vector<std::vector<int> > jCached_;
// 		std::vector<size_t> alpha_,beta_,gamma_;
// 		
// 		
// 		void createTcOperators(std::vector<SparseMatrixType>& basistc,const DmrgBasisWithOperatorsType& basis)
// 		{
// 			for (size_t i=0;i<basistc.size();i++) basistc[i] = transposeConjugate(basis.getOperatorByIndex(i).data);
// 		}
		
// 		void buildAlphaAndBeta()
// 		{
// 			//if (DmrgBasisType::useSu2Symmetry()) 
// 			buildAlphaAndBetaSu2();
// 			//else buildAlphaAndBetaLocal();
// 		}
// 		
// 		void buildAlphaAndBetaSu2()
// 		{
// 			int m = m_;
// 			//std::cerr<<"About to allocated block in outer: bs="<<bs<<"\n";
// 			int offset = basis1_.partition(m);
// 			//std::cerr<<"doing partition (outer)"<<m<<" with size="<<bs<<"\n";
// 			int total = basis1_.partition(m+1) - offset;
// 			int ns = basis2_.size();
// 			size_t counter=0;
// 			SparseMatrixType factorsInverse;
// 			transposeConjugate(factorsInverse,factors_);
// 			for (int i=0;i<total;i++) {
// // 				std::cerr<<"state i="<<i<<" jm="<<basis1_.jmValue(i+offset)<<" electrons="<<basis1_.electrons(i+offset);
// // 				std::cerr<<" flavor="<<basis1_.getFlavor(i+offset)<<"\n";
// 				std::vector<int> alphaTmp,betaTmp;
// 				int ii = basis1_.permutation(i+offset);
// 				
// 				gamma_.push_back(counter);
// 				
// 				for (int k=factorsInverse.getRowPtr(ii);k<factorsInverse.getRowPtr(ii+1);k++) {
// 					
// 					int alpha,beta;
// 					utils::getCoordinates(alpha,beta,factorsInverse.getCol(k),ns);
// 					alpha_.push_back(alpha);
// 					beta_.push_back(beta);
// 					factorsInverseValue_.push_back(factorsInverse.getValue(k));
// 					counter++;
// 					
// 					
// 				}
// 				
// 				
// 			}
// 			gamma_.push_back(counter);
// 							
// 			for (int ii=0;ii<factors_.rank();ii++) {
// 				std::vector<int> jVector;
// 				for (int k4=factors_.getRowPtr(ii);k4<factors_.getRowPtr(ii+1);k4++) {
// 					int j=basis1_.permutationInverse(factors_.getCol(k4)) -offset;
// 					jVector.push_back(j);
// 				}
// 				jCached_.push_back(jVector);
// 			}
// 		}
	
			
// 	template<typename Block_,typename MatrixElementType_,template<typename> class SparseMatrixTemplate,typename OperatorsType_,
// 		typename ReflectionSymmetryType_,
// 		typename ConcurrencyType_>
// 	ClebschGordanCached<MatrixElementType_> DmrgModelHelperSu2<Block_,MatrixElementType_,SparseMatrixTemplate,OperatorsType_,ReflectionSymmetryType_,ConcurrencyType_>::
// 	cgObject_(20);
		
}
