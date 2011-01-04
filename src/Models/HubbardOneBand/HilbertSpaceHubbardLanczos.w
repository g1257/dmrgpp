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

%\title{The HilbertSpaceHubbardLanczos Class}
%\author{G.A.}
%\maketitle

\begin{comment}
@o HilbertSpaceHubbardLanczos.h -t
@{
/*
@i license.txt
*/
@}
\end{comment}

\section{The Hilbert space for a Hubbard Model to use with Lanczos}


@{
@d bitctn
int bitcnt (word b)
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


@{
@d constructor
void BasisLanczos(size_t nsite, size_t npart)
{
  WordType ket;
  int i, n, m;
	
  /* compute size of basis */
  for (n=nsite,m=1,hilbert=1;m<=npart;hilbert=hilbert*n/m,n--,m++);
  cerr<<"size of basis: "<<hilbert<<endl;

  if (basis.size()!=hilbert) {
  	basis.clear();
	basis.resize(hilbert);
  } 
  
  if (npart==0) {
	  basis[0]=0;
	  return;
  }
  /* define basis states */
  ket = (1ul<<npart)-1;
  for (i=0;i<hilbert;i++){
    basis[i] = ket;
    for (n=m=0;(ket&3)!=1;n++,ket>>=1) {
		// cout<<"ket = "<<ket<<" n="<<n<<endl;
     	 m += ket&1;
	 }
    ket = ((ket+1)<<n) ^ ((1<<m)-1);
  }	
 
}
@}


@o HilbertSpaceHubbardLanczos.h -t
@{
} // namespace Dmrg

#endif
@}
\end{document}

