\documentclass{report}
\usepackage[T1]{fontenc}
\usepackage{bera}

\usepackage[pdftex,usenames,dvipsnames]{color}
\usepackage{listings}
\definecolor{mycode}{rgb}{0.9,0.9,1}
\lstset{language=c++,tabsize=1,basicstyle=\scriptsize,backgroundcolor=\color{mycode}}

\usepackage{hyperref}
\usepackage{fancyhdr}
\pagestyle{fancy}
\usepackage{verbatim}
\begin{document}

%\title{The HubbardLanczos Class}
%\author{G.A.}
%\maketitle

\begin{comment}
@o HubbardLanczos.h -t
@{
/*
@i license.txt
*/
@}
\end{comment}

\section{A Hubbard model to use with Lanczos}

HEre is some boilerplate:

@o HubbardLanczos.h -t
@{
#ifndef HUBBARDLANCZOS_H
#define HUBBARDLANCZOS_H

#include "Utils.h"
#include "CrsMatrix.h"
#include "BasisHubbardLanczos.h"

namespace Dmrg {
	@<theClassHere@>
} // namespace
#endif

@}

And the class is:
@d theClassHere
@{
template<typename RealType_,typename IoInputType,typename ParametersType,
	typename GeometryType>
class HubbardLanczos {
	@<privateTypesAndEnums@>
public:
	@<publicTypesAndEnums@>
	@<constructor@>
	@<publicFunctions@>
private:
	@<privateFunctions@>
	@<privateData@>
}; // class HubbardLanczos @}

@d privateTypesAndEnums
@{
typedef psimag::Matrix<RealType_> MatrixType; @}

@d privateData
@{
const ParametersType& mp_;
const GeometryType& geometry_;
size_t nup_;
size_t ndown_;
@}

@d publicTypesAndEnums
@{
typedef BasisHubbardLanczos BasisType;
typedef typename BasisType::WordType WordType;
typedef RealType_ RealType;
typedef CrsMatrix<RealType> SparseMatrixType;
typedef std::vector<RealType> VectorType;
enum {SPIN_UP,SPIN_DOWN};
enum {DESTRUCTOR,CONSTRUCTOR};
@}

@d constructor
@{
HubbardLanczos(IoInputType& io,const ParametersType& mp,GeometryType& geometry)
	: mp_(mp),geometry_(geometry)
{
	std::vector<RealType> qns;
	io.read(qns,"TargetQuantumNumbers");
	if (qns.size()<2) throw std::runtime_error("HubbardLanczos::ctor(...)\n");
	nup_=geometry.numberOfSites()*qns[0];
	ndown_=geometry.numberOfSites()*qns[1];
}
@}

@d publicFunctions
@{
@<setupHamiltonian@>
@<getOperator@>
@}

@d setupHamiltonian
@{
void setupHamiltonian(SparseMatrixType &matrix) const
{
	size_t nsite = geometry_.numberOfSites();
	BasisType basis1(nsite,nup_);
	BasisType basis2(nsite,ndown_);
	setupHamiltonian(matrix,basis1,basis2);
}
@}

@d getOperator
@{
void getOperator(SparseMatrixType& matrix,size_t what,size_t i,size_t flavor) const
{
	throw std::runtime_error("getOperator(...): Unimplemented\n");
}
@}

@d privateFunctions
@{
@<hoppings@>
@<setupHamiltonianP@>
@<countNonZero@>
@<perfectIndex@>
@<doSign1@>
@<doSign2@>
@}

@d setupHamiltonianP
@{
void setupHamiltonian(SparseMatrixType &matrix,const BasisType &basis1,const BasisType& basis2) const
{
	// Calculate diagonal elements AND count non-zero matrix elements
	size_t hilbert1=basis1.size();
	size_t hilbert2=basis2.size();
	MatrixType diag(hilbert2,hilbert1);
	size_t nzero = countNonZero(diag,basis1,basis2);
	
	size_t nsite = geometry_.numberOfSites();
	
	// Setup CRS matrix
	matrix.resize(hilbert1*hilbert2,nzero);
	
	// Calculate off-diagonal elements AND store matrix
	size_t nCounter=0;
	matrix.setRow(0,0);
	for (size_t ispace1=0;ispace1<hilbert1;ispace1++) {
		WordType ket1 = basis1[ispace1];
		for (size_t ispace2=0;ispace2<hilbert2;ispace2++) {
			WordType ket2 = basis2[ispace2];
			// Save diagonal
			matrix.setCol(nCounter,ispace2+ispace1*hilbert2);
			RealType cTemp=diag(ispace2,ispace1);
			matrix.setValues(nCounter,cTemp);
			nCounter++;
			for (size_t i=0;i<nsite;i++) {
				WordType s1i=(ket1 & BasisType::bitmask(i));
				WordType s2i=(ket2 & BasisType::bitmask(i));
				if (s1i>0) s1i=1;
				if (s2i>0) s2i=1;
				
				// Hopping term 
				for (size_t j=0;j<nsite;j++) {
					if (j<i) continue;
					RealType tmp = hoppings(i,j);
					if (tmp==0) continue;
					WordType s1j= (ket1 & BasisType::bitmask(j));
					WordType s2j= (ket2 & BasisType::bitmask(j));
					if (s1j>0) s1j=1;
					if (s2j>0) s2j=1;
					if (s1i+s1j==1) {
						WordType bra1= ket1 ^(BasisType::bitmask(i)|BasisType::bitmask(j));
						size_t temp = perfectIndex(basis1,basis2,bra1,ket2);
						matrix.setCol(nCounter,temp);
						cTemp=hoppings(i,j)*doSign(ket1,ket2,i,j,SPIN_DOWN); // check SIGN FIXME
						if (cTemp==0.0) {
							std::cerr<<"ctemp=0 and hopping="<<hoppings(i,j)<<" and i="<<i<<" and j="<<j<<"\n";
						}
						matrix.setValues(nCounter,cTemp);
						nCounter++;
					}
					if (s2i+s2j==1) {
						WordType bra2= ket2 ^(BasisType::bitmask(i)|BasisType::bitmask(j));
						size_t temp = perfectIndex(basis1,basis2,ket1,bra2);
						matrix.setCol(nCounter,temp);
						cTemp=hoppings(i,j)*doSign(ket1,ket2,i,j,SPIN_UP); // Check SIGN FIXME
						matrix.setValues(nCounter,cTemp);
						nCounter++;					
					}
				}
			}
			matrix.setRow(ispace2+hilbert2*ispace1+1,nCounter);
		}
	}
	matrix.setRow(hilbert1*hilbert2,nCounter);
}
@}

@d hoppings
@{
RealType hoppings(size_t i,size_t j) const
{
	return geometry_(i,0,j,0,0);
}
@}
		
@d countNonZero
@{
size_t countNonZero(MatrixType& diag,const BasisType &basis1,const BasisType& basis2) const
{
	size_t hilbert1=basis1.size();
	size_t hilbert2=basis2.size();
	size_t nsite = geometry_.numberOfSites();

	// Calculate diagonal elements AND count non-zero matrix elements
	size_t nzero = 0;
	for (size_t ispace1=0;ispace1<hilbert1;ispace1++) {
		WordType ket1 = basis1[ispace1];
		for (size_t ispace2=0;ispace2<hilbert2;ispace2++) {
			WordType ket2 = basis2[ispace2];
			RealType s=0;
			for (size_t i=0;i<nsite;i++) {
				WordType s1i=(ket1 & BasisType::bitmask(i));
				WordType s2i=(ket2 & BasisType::bitmask(i));
				if (s1i>0) s1i=1;
				if (s2i>0) s2i=1;
				
				// Hubbard term
				if (s1i>0 && s2i>0 ) s += mp_.hubbardU[i];
						
				// Potential term
				if (s1i>0) s += mp_.potentialV[i];
				if (s2i>0) s += mp_.potentialV[i];
				
				// Hopping term (only count how many non-zero)
				for (size_t j=0;j<nsite;j++) {
					if (j<i) continue;
					RealType tmp = hoppings(i,j);
					if (tmp==0) continue;
					
					WordType s1j= (ket1 & BasisType::bitmask(j));
					WordType s2j= (ket2 & BasisType::bitmask(j));
					if (s1j>0) s1j=1;
					if (s2j>0) s2j=1;
					if (s1i+s1j==1) nzero++;
					if (s2i+s2j==1) nzero++;
				}
			}
			// cout<<"diag of ("<<ispace1<<","<<ispace2<<"): "<<s<<endl;
			diag(ispace2,ispace1)=s;
			nzero++;
		}
	}

	nzero++;
	return nzero;
}
@}

@d perfectIndex
@{
size_t perfectIndex(const BasisType& basis1,const BasisType& basis2,WordType ket1,WordType ket2) const
{
	size_t hilbert2=basis2.size();
	size_t n1 = basis1.perfectIndex(ket1);
	size_t n2 = basis2.perfectIndex(ket2);

	return n2 + n1*hilbert2;
}
@}

@d doSign1
@{
int doSign(WordType a, WordType b,size_t i,size_t j,size_t sector) const
{
	if (i > j) {
		std::cerr<<"FATAL: At doSign\n";
		std::cerr<<"INFO: i="<<i<<" j="<<j<<std::endl;
		std::cerr<<"AT: "<<__FILE__<<" : "<<__LINE__<<std::endl;
		throw std::runtime_error("HubbardLanczos::doSign(...)\n");
	}

	WordType mask = a ^  b;
	mask &= ((1 << (i+1)) - 1) ^ ((1 << j) - 1);
	int s=(BasisType::bitcnt (mask) & 1) ? -1 : 1; // Parity of single occupied between i and j

	if (sector==SPIN_DOWN) { // Is there a down at j?
		if (BasisType::bitmask(j) & b) s = -s;
	}
	if (sector==SPIN_UP) { // Is there an up at i?
		if (BasisType::bitmask(i) & a) s = -s;
	}

	return s;
}
@}

@d doSign2
@{
int doSign(WordType a, WordType b,size_t i,size_t sector) const
{
	size_t nsite = geometry_.numberOfSites();
	if (i<nsite-1) {
		WordType mask = a ^ b;
		mask &= ((1 << (i+1)) - 1) ^ ((1 << nsite) - 1);
		int s=(BasisType::bitcnt (mask) & 1) ? -1 : 1; // Parity of single occupied between i and nsite-1
		//cout<<"sign1: "<<s<<"a="<<a<<" b="<<b<<"sector="<<sector<<endl;
		if (sector==SPIN_UP) { // Is there an up at i?
			if (BasisType::bitmask(i) & a) s = -s;
			//cout<<"sign2: "<<s<<endl;
		}
		return s;
	}

	int s=1;
	if (sector==SPIN_UP) { // Is there an up at i?
		if (BasisType::bitmask(i) & a) s = -s;
	}

	return s;
}
@}
			
\end{document}

