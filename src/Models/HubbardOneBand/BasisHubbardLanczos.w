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

\section{The Hilbert space for a Hubbard Model to use with Lanczos}

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
template<typename WordType>
class BasisHubbardLanczos {
		@<privateTypesAndEnums@>
	public:
		@<constructor@>
	private:
		@<privateFunctions@>
		@<privateData@>
}; // class BasisHubbardLanczos
@}

@d privateData
@{
std::vector<WordType> data_;
@}

All right, now the constructor:
@d constructor
@{
BasisHubbardLanczos(size_t nsite, size_t npart)
{
	/* compute size of basis */
	int hilbert=1;
	for (int n=nsite,int m=1;m<=npart;hilbert=hilbert*n/m,n--,m++);
	
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
	for (i=0;i<hilbert;i++) {
		basis[i] = ket;
		size_t n=m=0;
		for (;(ket&3)!=1;n++,ket>>=1) {
			m += ket&1;
		}
		ket = ((ket+1)<<n) ^ ((1<<m)-1);
	}
}
@}

@d privateFunctions
@{
@<bitctn@>
@} 

@d bitctn
@{
int bitcnt (WordType b)
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
}					
@}


\end{document}

