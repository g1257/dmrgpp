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

%\title{The Powell Class for Minimization of Functions}
%\author{G.A.}
%\maketitle

\begin{comment}
@o Minimizer.h -t
@{
/*
@i license.txt
*/
@}
\end{comment}

\section{A class to minimize a function without using Derivatives}
This requires the GNU Scientific Library, libraries and devel. headers.

Explain usage here. FIXME.

This is the structure of this file:

@o Minimizer.h -t
@{
#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <iostream>
#include <vector>
#include <stdexcept>
extern "C" {
#include <gsl/gsl_multimin.h>
}

namespace PsimagLite {
	@<MockVector@>
	@<MyFunction@>
	@<theClassHere@>
} // namespace PsimagLite
#endif // MINIMIZER_H

@}

And this is the class here:

@d theClassHere
@{
template<typename RealType,typename FunctionType>
class Minimizer {
	@<privateTypedefsAndConstants@>
public:
	@<constructor@>
	@<destructor@>
	@<publicFunctions@>
private:
	@<privateFunctions@>
	@<privateData@>
}; // class Minimizer
@}


@d privateTypedefsAndConstants
@{
typedef typename FunctionType::FieldType FieldType;
typedef typename Vector<FieldType>::Type VectorType;
typedef Minimizer<RealType,FunctionType> ThisType;

@}

@d privateData
@{
FunctionType& function_;
size_t maxIter_;
const gsl_multimin_fminimizer_type *gslT_;
gsl_multimin_fminimizer *gslS_;
@}

@d constructor
@{
Minimizer(FunctionType& function,size_t maxIter)
		: function_(function),
		  maxIter_(maxIter),
		  gslT_(gsl_multimin_fminimizer_nmsimplex2),
		  gslS_(gsl_multimin_fminimizer_alloc(gslT_,function_.size()))
{
}
@}

@d publicFunctions
@{
@<simplex@>
@}

@d destructor
@{
~Minimizer()
{
	gsl_multimin_fminimizer_free (gslS_);
}
@}


@d simplex
@{
int simplex(VectorType& minVector,RealType delta=1e-3,RealType tolerance=1e-3)
{
	gsl_vector *x;
	/* Starting point,  */
	x = gsl_vector_alloc (function_.size());
	for (size_t i=0;i<minVector.size();i++)
		gsl_vector_set (x, i, minVector[i]);

	gsl_vector *xs;
	xs = gsl_vector_alloc (function_.size());
	for (size_t i=0;i<minVector.size();i++)
		gsl_vector_set (xs, i, delta);

	gsl_multimin_function func;
	func.f= myFunction<FunctionType>;
	func.n = function_.size();
	func.params = &function_;
	gsl_multimin_fminimizer_set (gslS_, &func, x, xs);

	for (size_t iter=0;iter<maxIter_;iter++) {
		int status = gsl_multimin_fminimizer_iterate (gslS_);

		if (status) throw std::runtime_error("Minimizer::simplex(...): Error encountered\n");

		RealType size = gsl_multimin_fminimizer_size(gslS_);
		status = gsl_multimin_test_size(size, tolerance);

		if (status == GSL_SUCCESS) {
			found(minVector,gslS_->x,iter);
			gsl_vector_free (x);
			gsl_vector_free (xs);
			return iter;
		}
	}
	gsl_vector_free (x);
	gsl_vector_free (xs);
	return -1;
}
@}

@d privateFunctions
@{
@<found@>
@}

@d found
@{
void found(VectorType& minVector,gsl_vector* x,size_t iter)
{
	for (size_t i=0;i<minVector.size();i++)
		minVector[i] = gsl_vector_get(x,i);
}
@}

@d MockVector
@{
template<typename FieldType>
class MockVector {
public:
	MockVector(const gsl_vector *v) : v_(v)
	{
	}
	const FieldType& operator[](size_t i) const
	{
		return v_->data[i];
	}
	size_t size() const { return v_->size; }
private:
	const gsl_vector *v_;
}; // class MockVector
@}

@d MyFunction
@{
template<typename FunctionType>
typename FunctionType::FieldType myFunction(const gsl_vector *v, void *params)
{
	MockVector<typename FunctionType::FieldType> mv(v);
	FunctionType* ft = (FunctionType *)params;
	return ft->operator()(mv);
}
@}

\end{document}

