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
	typename RealType,typename SparseMatrixType,typename VectorType>
class DynamicFunctional  {
public:
	@<publicTypedefs@>
	@<constructor@>
	@<publicFunctions@>
private:
	@<privateFunctions@>
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
		const VectorType& aVector,
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
const VectorType& aVector_;
RealType omega_;
RealType E0_;
RealType eta_;
PsimagLite::ProgressIndicator progress_;
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

This class needs to tell the minimizer (\verb=Minimizer.w=) class what type of scalar we'll be using.%'
The minimizer expects this in a type that has to be named \verb=FieldType=.

All right, so let us tell the minimizer that we are real:
@d publicTypedefs
@{
typedef RealType FieldType; // see documentation
//typedef std::complex<RealType> ComplexType;
//typedef typename PsimagLite::Vector<ComplexType>::Type VectorComplexType;
@}


@d publicFunctions
@{
@<operatorParens@>
@<size@>
@}

The function below computes $W_{A,\eta}(\omega,\psi)$ for a fixed
$A$, $\eta$, and $\omega$. See Eq.~(14) of reference \cite{re:jeckelmann02}.
@d operatorParens
@{
template<typename SomeVectorType>
RealType operator()(const SomeVectorType &v) const
{
	VectorType vC(v.size());
	for (SizeType i=0;i<v.size();i++) vC[i] = v[i];
	VectorType x(vC.size(),0.0);

	H_.matrixVectorProduct(x,vC); // x += H_ vC
	RealType sum = utils::square(E0_+omega_) + utils::square(eta_);
	sum *= std::real(vC*vC);
	sum -= 2*(E0_+omega_)*std::real(x*vC);
	sum += std::real(x*x);
	sum += 2*eta_*std::real(aVector_*vC);
	//checkProducts(vC,x);
	return sum;
}
@}

The size of the vectors for this functional is equal to the rank of the Hamiltonian sector we are
considering; the latter is stored in the private member \verb=H_=.
@d size
@{
SizeType size() const {return H_.rank(); }
@}

@d privateFunctions
@{
@}

\bibliographystyle{plain}
\bibliography{thesis}

\end{document}
