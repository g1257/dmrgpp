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

namespace Dmrg {
	@<theClassHere@>
} // namespace
#endif

@}

And the class is:
@d theClassHere
@{
template<typename ParametersType,typename GeometryType>
class HubbardLanczos {
	@<privateTypesAndEnums@>
public:
	@<publicTypesAndEnums@>
	@<constructor@>
	@<publicFunctions@>
private:
	@<privateFunctions@>
	@<privateData@>
}; // class HubbardLanczos
@}

@d privateTypesAndEnums
@{
typedef unsigned int long long WordType;
@}

@d privateData
@{
const ParametersType& mp_;
const GeometryType& geometry_;
@}

@d publicTypesAndEnums
@{
enum {SPIN_UP,SPIN_DOWN};
enum {DESTRUCTOR,CONSTRUCTOR};
@}

@d constructor
@{
	HubbardLanczos(const ParametersType& mp,GeometryType& geometry) 
	: mp_(mp),geometry_(geometry)
	{
	}
@}

@d publicFunctions
@{
@<setupHamiltonian@>
@<getOperator@>
@}

@d setupHamiltonian
@{
void setupHamiltonian(SparseMatrixType &matrix,const BasisType &basis1,const BasisType& basis2)
{
	// Calculate diagonal elements AND count non-zero matrix elements
	size_t hilbert1=basis1.size();
	size_t hilbert2=basis2.size();
	MatrixType diag(hilbert2,hilbert1);
	size_t nzero = countNonZero(diag,basis1,basis2);
	
	size_t nsite = geometry_.numberOfSites();
	
	
	// Setup CRS matrix
	matrix.initRow(hilbert1*hilbert2+1);
	matrix.initValues(nzero);
	matrix.initCol(nzero);
	
	// Calculate off-diagonal elements AND store matrix
	size_t nCounter=0;
	matrix.setRow(0,0);
	for (size_t ispace1=0;ispace1<hilbert1;ispace1++) {
		WordType ket1 = basis1[ispace1];
		for (size_t ispace2=0;ispace2<hilbert2;ispace2++) {
			WordType ket2 = basis2[ispace2];
			// Save diagonal
			matrix.setCol(nCounter,ispace2+ispace1*hilbert2);
			ComplexType cTemp=diag[ispace2+ispace1*hilbert2];
			matrix.setValues(nCounter,cTemp);
			nCounter++;
			for (size_t i=0;i<nsite;i++) {
				WordType s1i=(ket1 & ether.bitmask[i]);
				WordType s2i=(ket2 & ether.bitmask[i]);
				if (s1i>0) s1i=1;
				if (s2i>0) s2i=1;
				
				// Hopping term 
				for (size_t j=0;j<nsite;j++) {
					if (j<i) continue;
					RealType tmp = dmrgGeometry_(i,0,j,0,0);
					if (tmp==0) continue;
					WordType s1j= (ket1 & ether.bitmask[j]);
					WordType s2j= (ket2 & ether.bitmask[j]);
					if (s1j>0) s1j=1;
					if (s2j>0) s2j=1;
					if (s1i+s1j==1) {
						WordType bra1= ket1 ^(ether.bitmask[i]|ether.bitmask[j]);
						temp = perfectIndex(basis1,basis2,bra1,ket2,ether);
						matrix.setCol(nCounter,temp);
						cTemp=ether.hoppings[i+nsite*j]*doSign(ket1,ket2,i,j,0,ether); // check SIGN FIXME
						if (cTemp==0.0) {
							std::cerr<<"ctemp=0 and hopping="<<ether.hoppings[i+nsite*j]<<" and i="<<i<<" and j="<<j<<"\n";
						}
						matrix.setValues(nCounter,cTemp);
						nCounter++;
					}
					if (s2i+s2j==1) {
						WordType bra2= ket2 ^(ether.bitmask[i]|ether.bitmask[j]);
						size_t temp = perfectIndex(basis1,basis2,ket1,bra2,ether);
						matrix.setCol(nCounter,temp);
						cTemp=ether.hoppings[i+nsite*j]*doSign(ket1,ket2,i,j,1,ether); // Check SIGN FIXME
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

@d getOperator
@{
double calcNiup(int iind,int jind,vector<word> const &basis1,vector<word> const &basis2,vector<Complex> const &gsVector,
Parameters const &ether)
{
	int hilbert1=basis1.size();
	int hilbert2=basis2.size();
	word bra1,bra2,ket1,ket2,s1j,s2j,s1i,s2i;
	int temp,ispace1,ispace2;
	double sum=0;
	double normaliz=0;
		
	for (ispace1=0;ispace1<hilbert1*hilbert2;ispace1++) {	
		normaliz += real(conj(gsVector[ispace1])*gsVector[ispace1]);
	}
	cout<<"Normalization="<<normaliz<<endl;
	
	for (ispace1=0;ispace1<hilbert1;ispace1++) {
		ket1 = basis1[ispace1];
		// cout<<"here ket1="<<ket1<<endl;
		temp = perfectIndex(basis1, ket1,ether);
		
		for (ispace2=0;ispace2<hilbert2;ispace2++) {
			ket2 = basis2[ispace2];
			s1i=(ket1 & ether.bitmask[iind]);
			s2i=(ket2 & ether.bitmask[iind]);
			if (s1i>0) s1i=1;
			if (s2i>0) s2i=1;
			s1j= (ket1 & ether.bitmask[jind]);
			s2j= (ket2 & ether.bitmask[jind]);
			if (s1j>0) s1j=1;
			if (s2j>0) s2j=1;
			if (s1i+s1j==1) {
				bra1= ket1 ^(ether.bitmask[iind]|ether.bitmask[jind]);
				temp = perfectIndex(basis1,basis2,bra1,ket2,ether);
				sum += real(conj(gsVector[ispace1+hilbert1*ispace2])*gsVector[temp])*doSign(ket1,ket2,iind,jind,0,ether); 
			}
			if (s2i+s2j==1) {
				bra2= ket2 ^(ether.bitmask[iind]|ether.bitmask[jind]);
				temp = perfectIndex(basis1,basis2,ket1,bra2,ether);
				sum += real(conj(gsVector[ispace1+hilbert1*ispace2])*gsVector[temp])*doSign(ket1,ket2,iind,jind,1,ether); 							
			}
			if (s1i==1 && iind==jind) {
				sum += real(conj(gsVector[ispace1+hilbert1*ispace2])*gsVector[ispace1+hilbert1*ispace2]);
			}
			if (s2i==1 && iind==jind) {
				sum += real(conj(gsVector[ispace1+hilbert1*ispace2])*gsVector[ispace1+hilbert1*ispace2]); 
			}
			 
		}
	}
	return sum;
}
@}

@d privateFunctions
@{
@<countNonZero@>
@<perfectIndex1@>
@<perfectIndex2@>
@<doSign1@>
@<doSign2@>
@}

@d countNonZero
@{
size_t countNonZero(MatrixType& diag)
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
				WordType s1i=(ket1 & ether.bitmask[i]);
				WordType s2i=(ket2 & ether.bitmask[i]);
				if (s1i>0) s1i=1;
				if (s2i>0) s2i=1;
				
				// Hubbard term
				if (s1i>0 && s2i>0 ) s += ether.U[i];
				
				// Potential term
				if (s1i>0) s += ether.Vg[i];
				if (s2i>0) s += ether.Vg[i];
				
				// Hopping term (only count how many non-zero)
				for (size_t j=0;j<nsite;k++) {
					if (j<i) continue;
					RealType tmp = dmrgGeometry_(i,0,j,0,0);
					if (tmp==0) continue;
					
					s1j= (ket1 & ether.bitmask[j]);
					s2j= (ket2 & ether.bitmask[j]);
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

@d perfectIndex1
@{
size_t perfectIndex(vector<word> const &basis, word state,Parameters const &ether)
{
	size_t nsite=geometry_.numberOfSites();
	for (size_t b=n=0,c=1;state>0;b++,state>>=1)
		if (state&1) n += ether.comb[b+nsite*(c++)];

	return n;
}
@}

@d perfectIndex2
@{
int perfectIndex(vector<word> const &basis1,vector<word> const &basis2,
		word ket1,word ket2,Parameters const &ether)
{
	int hilbert1=basis1.size(), hilbert2=basis2.size();
	int n1,n2,n;
	n1 = perfectIndex(basis1,ket1,ether);
	n2 = perfectIndex(basis2,ket2,ether);
	
	n = n2 + n1*hilbert2;
	return n;
}
@}

@d doSign1
@{
int doSign(word a, word b,int i,int j,int sector,Parameters const &ether)
{
	word mask;
	int s;
	if (i > j) {
		cerr<<"FATAL: At doSign\n";
		cerr<<"INFO: i="<<i<<" j="<<j<<endl;
		cerr<<"AT: "<<__FILE__<<" : "<<__LINE__<<endl;
		exit(1);
	}
	mask = a ^  b;
	mask &= ((1 << (i+1)) - 1) ^ ((1 << j) - 1);
	s=(bitcnt (mask) & 1) ? -1 : 1; // Parity of single occupied between i and j
	
	if (sector==0) { // Is there a down at j?
		if (ether.bitmask[j] & b) s = -s;
	}
	if (sector==1) { // Is there an up at i?
		if (ether.bitmask[i] & a) s = -s;
	}
	/* cout<<"Sector="<<sector<<" Ether bitmask[i]="<<(ether.bitmask[i] & a);
	cout<<" Ether bitmask[j]="<<(ether.bitmask[j] & b)<<" a="<<a<<" b="<<b;
	cout<<"sign="<<s<<endl; */
	return s;
}
@}

@d doSign2
@{
int doSign(word a,word b,int i,int sector,Parameters const &ether,int nsite)
{
	int s;
	word mask;
		
	if (i<nsite-1) {
		mask = a ^ b;
		mask &= ((1 << (i+1)) - 1) ^ ((1 << nsite) - 1);
		s=(bitcnt (mask) & 1) ? -1 : 1; // Parity of single occupied between i and nsite-1
		//cout<<"sign1: "<<s<<"a="<<a<<" b="<<b<<"sector="<<sector<<endl;
		if (sector==1) { // Is there an up at i?
			if (ether.bitmask[i] & a) s = -s;
			//cout<<"sign2: "<<s<<endl;
		}
	} else {
		s=1;
		if (sector==1) { // Is there an up at i?
			if (ether.bitmask[i] & a) s = -s;
		}
	}
	
	/* cout<<"Sector="<<sector<<" Ether bitmask[i]="<<(ether.bitmask[i] & a);
	cout<<" Ether bitmask[j]="<<(ether.bitmask[j] & b)<<" a="<<a<<" b="<<b;
	cout<<"sign="<<s<<endl; */
	return s;
}
@}
			
\end{document}

