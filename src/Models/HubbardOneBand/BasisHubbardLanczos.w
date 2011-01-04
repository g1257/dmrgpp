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

\begin{comment}
@o BasisHubbardLanczos.h -t
@{
/*
@i license.txt
*/
@}
\end{comment}

\section{Basis a Hubbard Model to use with Lanczos}

HEre is some boilerplate:

@o BasisHubbardLanczos.h -t
@{
#ifndef BASISHUBBARDLANCZOS_H
#define BASISHUBBARDLANCZOS_H

#include "Utils.h"

namespace Dmrg {
	@<theClassHere@>
} // namespace
#endif

@}

And the class is:
@d theClassHere
@{
class BasisHubbardLanczos {
public:
	@<publicTypesAndEnums@>
	@<constructor@>
	@<publicFunctions@>

private:
	@<privateFunctions@>
	@<privateData@>
}; // class BasisHubbardLanczos
@<staticDefinitions@>
@}

@d privateData
@{
size_t size_;
size_t npart_;
std::vector<WordType> data_;
@}

@d publicTypesAndEnums
@{
typedef unsigned int long long WordType;
static size_t nsite_;
static psimag::Matrix<size_t> comb_;
static std::vector<WordType> bitmask_; @}

All right, now the constructor:
@d constructor
@{
BasisHubbardLanczos(size_t nsite, size_t npart) : npart_(npart)
{
	if (nsite_>0 && nsite!=nsite_)
		throw std::runtime_error("BasisHubbardLanczos: All basis must have same number of sites\n");
	nsite_ = nsite;
	doCombinatorial();
	doBitmask();

	/* compute size of basis */
	size_t hilbert=1;
	int n=nsite;
	size_t m=1;
	for (;m<=npart;n--,m++)
		hilbert=hilbert*n/m;

	if (data_.size()!=hilbert) {
		data_.clear();
		data_.resize(hilbert);
	}

	if (npart==0) {
		data_[0]=0;
		return;
	}
	
	/* define basis states */
	WordType ket = (1ul<<npart)-1;
	for (size_t i=0;i<hilbert;i++) {
		data_[i] = ket;
		n=m=0;
		for (;(ket&3)!=1;n++,ket>>=1) {
			m += ket&1;
		}
		ket = ((ket+1)<<n) ^ ((1<<m)-1);
	}
	size_ = hilbert;
} @}

@d publicFunctions
@{
@<size@>
@<operatorBracket@>
@<perfectIndex@>
@<bitmask@>
@<bitctn@> @}

@d size
@{
size_t size() const { return size_; } @}

@d operatorBracket
@{
const WordType& operator[](size_t i) const
{
	return data_[i];
} @}

@d perfectIndex
@{
size_t perfectIndex(WordType state) const
{
	size_t n=0;
	for (size_t b=0,c=1;state>0;b++,state>>=1)
		if (state&1) n += comb_(b,c++);

	return n;
} @}

@d bitmask
@{
static const WordType& bitmask(size_t i)
{
	return bitmask_[i];
} @}

@d bitctn
@{
static int bitcnt (WordType b)
{
#if (ULONG_MAX == 0xfffffffful)
	b = (0x55555555 & b) + (0x55555555 & (b >> 1));
	b = (0x33333333 & b) + (0x33333333 & (b >> 2));
	b = (0x0f0f0f0f & b) + (0x0f0f0f0f & (b >> 4));
	b = (0x00ff00ff & b) + (0x00ff00ff & (b >> 8));
	b = (0x0000ffff & b) + (0x0000ffff & (b >> 16));

	return (int) b;
#else
	b = (0x5555555555555555 & b) + (0x5555555555555555 & (b >> 1));
	b = (0x3333333333333333 & b) + (0x3333333333333333 & (b >> 2));
	b = (0x0f0f0f0f0f0f0f0f & b) + (0x0f0f0f0f0f0f0f0f & (b >> 4));
	b = (0x00ff00ff00ff00ff & b) + (0x00ff00ff00ff00ff & (b >> 8));
	b = (0x0000ffff0000ffff & b) + (0x0000ffff0000ffff & (b >> 16));
	b = (0x00000000ffffffff & b) + (0x00000000ffffffff & (b >> 32));

	return (int) b;
#endif
} @}

@d privateFunctions
@{
@<doCombinatorial@>
@<doBitMask@>
@} 

@d doCombinatorial
@{
void doCombinatorial()
{
	/* look-up table for binomial coefficients */
	comb_.resize(nsite_,nsite_);

	for (size_t n=0;n<nsite_;n++)
		for (size_t i=0;i<nsite_;i++)
			comb_(n,i)=0;

	for (size_t n=0;n<nsite_;n++) {
		size_t m = 0;
		int j = n;
		size_t i = 1;
		size_t cnm  = 1;
		for (;m<=n/2;m++,cnm=cnm*j/i,i++,j--)
			comb_(n,m) = comb_(n,n-m) = cnm;
	}
} @}

@d doBitMask
@{
void doBitmask()
{
	bitmask_.resize(nsite_);
	bitmask_[0]=1ul;
	for (size_t i=1;i<nsite_;i++)
		bitmask_[i] = bitmask_[i-1]<<1;
} @}

@d staticDefinitions
@{
size_t BasisHubbardLanczos::nsite_=0;
psimag::Matrix<size_t> BasisHubbardLanczos::comb_;
std::vector<typename BasisHubbardLanczos::WordType> BasisHubbardLanczos::bitmask_; @}
	
\end{document}

