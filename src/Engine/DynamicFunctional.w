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

%\title{The DynamicFunctional Class}
%\author{G.A.}
%\maketitle

\begin{comment}
@o DynamicFunctional.h -t
@{
/*
@i license.txt
*/
@}
\end{comment}

\section{Dynamic Functional Needed for Dynamic Targetting}
This class implements the functional Eq.~(14) of reference \cite{re:jeckelmann02}.
See also DynamicTargetting, FIXME CROSS REFERENCE HERE.

@o DynamicFunctional.h -t
@{
#ifndef DYN_FUNCTIONAL_H
#define DYN_FUNCTIONAL_H

#include "ProgressIndicator.h"

namespace Dmrg {
	@<theClassHere@>
} // namespace
#endif // DYN_FUNCTIONAL_H

@}

@d theClassHere
@{
template<
	typename RealType,typename SparseMatrixType,typename VectorWithOffsetType>
class DynamicFunctional  {
public:
	@<publicTypedefs@>
	@<constructor@>
	@<publicFunctions@>
private:
	@<privateData@>
}; // class DynamicFunctional
@}

Now comes the constructor which takes 6 arguments.
The first 3 arguments are the system (left-block), environment (right-block), and superblock (system + environment).
As usual, the first 2 are heavy objects---with operators---, and the superblock is light.
The 4th argument is the model object. The 5th argument is a \verb|TargettingStructureType| object 
which is a \verb|TargettingStructureParms| object (as you can see in~@xtargettingstructure@x).
A structure is just a bunch of data bundled together, and you can see this in the file \verb|TargetStructureParams.h|.
The last argument is a \verb|WaveFunctionTransformation| object. More info about this class is
in \verb|WaveFunctionTransformation.h|.
@d constructor
@{
DynamicFunctional(
		const SparseMatrixType& H,
		const VectorWithOffsetType& aVector,
		RealType omega,
		RealType E0,
		RealType eta)
:	@<stackInitialization@>
{}
@}

Now let us look at the private data of this class:
@d privateData
@{
const SparseMatrixType& H_;
const VectorWithOffsetType& aVector_;
RealType omega_;
RealType E0_;
RealType eta_;
ProgressIndicator progress_;
@}

Now we get to the stack initialization of this object.

@d stackInitialization
@{
H_(H),
aVector_(aVector),
omega_(omega),
E0_(E0),
eta_(eta),
progress_("DynamicFunctional",0)
@}

This class needs to tell the minimizer (\verb=Minimizer.w=) class what type of scalar we'll be using%'
The minimizer expects this in type called \verb=FieldType=.
It turns out that here the vectors are complex, but the minimizer only takes real functions, since it is based
on the gsl which uses \verb=gsl_vector= which is real.
Arrrghhh!! To work around this problem we will map the complex vector \verb=v= of size \verb=n=
into a real vector \verb=vReal= of size \verb=2*n=. This is implemented here:
@d packComplex
@{
void packComplexToReal(std::vector<RealType>& svReal,const std::vector<std::complex<RealType> >& sv)
{
	svReal.resize(sv.size()*2);
	size_t j = 0;
	for (size_t i=0;i<sv.size();i++) {
		svReal[j++] = real(sv[i]);
		svReal[j++] = imag(sv[i]);
	}
}
@}

@d packReal
@{
void packRealToComplex(std::vector<std::complex<RealType> >& sv,const std::vector<RealType>& svReal)
{
	sv.resize(svReal.size()/2);
	size_t j = 0;
	for (size_t i=0;i<sv.size();i++) {
		sv[i] = std::complex<RealType>(svReal[j],svReal[j+1]);
		j += 2;
	}
}
@}

All right, so let us tell the minimizer that we are real:
@d publicTypedefs
@{
typedef RealType FieldType; // see documentation
@}


@d publicFunctions
@{
@<packComplex@>
@<packReal@>
@<operatorParens@>
@<size@>
@}

The function below computes $W_{A,\eta}(\omega,\psi) for a fixed
$A$, $\eta$, and $\omega$. See Eq.~(14) of reference \cite{re:jeckelmann02}.
@d operatorParens
@{
template<typename SomeVectorType>
RealType operator()(const SomeVectorType &v) const
{
	throw std::runtime_error("Neeeds implementation (sorry)\n");
}
@}

@d size
@{
size_t size() const {throw std::runtime_error("Neeeds implementation (sorry)\n");; }
@}

\bibliographystyle{plain}
\bibliography{thesis}

\end{document}
