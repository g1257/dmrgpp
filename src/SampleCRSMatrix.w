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

%\title{The SampleCRSMatrix Class}
%\author{G.A.}
%\maketitle

\begin{comment}
@o SampleCRSMatrix.h -t
@{
/*
@i license.txt
*/
@}
\end{comment}

\section{A class to represent a sparse matrix in Compressed Row Storage (CRS)}

HEre is some boilerplate:

@o SampleCRSMatrix.h -t
@{
#ifndef SAMPLE_CRSMATRIX_HEADER_H
#define SAMPLE_CRSMATRIX_HEADER_H

#include <iostream>
#include <vector>
#include <stdexcept>
#include <cstdlib> // <-- just for drand48 and srand48(seed)
#include "Sort.h"

namespace PsimagLite {
	@<theClassHere@>
	@<companions@>
} // namespace Dmrg
#endif

@}

The CRS format puts the subsequent nonzero elements of the matrix rows
in contiguous memory locations. We create 3 vectors: one for complex numbers containing the values of the
matrix entries
and the other two for integers ($colind$ and $rowptr$).
The vector $values$ stores the values of the non-zero elements of the matrix,
as they are traversed in a row-wise fashion.
The $colind$ vector stores the column indices of the elements of the $values$
vector. That is, if $values[k] = a[i][j]$ then $colind[k] = j$.
The $rowptr$ vector stores the locations in the $values$ vector that start
a row, that is $values[k] = a[i][j]$ if $rowptr[i] \le i < rowptr[i + 1]$.
By convention, we define $rowptr[N_{dim}]$ to be equal to the number of non-zero elements,
$n_z$, in the matrix. The storage savings of this approach are significant since instead of
storing $N_{dim}^2$ elements, we need only $2n_z + N_{dim} + 1$ storage locations.\\
To illustrate how the CRS format works, consider the non-symmetric matrix defined by
\begin{equation}
	A=\left[\begin{tabular}{llllll}

	10 &  0 & 0 & 0  & -2 & 0 \\
	3 &  9 &  0 &  0 &  0 &  3 \\
	0 &  7 &  8 &  7 &  0 &  0 \\
	3 &  0 &  8 &  7  & 5 &  0 \\
	0 &   8 &  0 &  9 &  9 & 13 \\
	0 &  4 &  0 &  0 &  2&  -1 \\
\end{tabular}\right]\end{equation}

The CRS format for this matrix is then specified by the arrays:

\begin{tt}
	values = [10 -2  3  9  3  7  8  7  3 ... 9 13  4  2 -1 ]\\
	colind = [ 0  4  0  1  5  1  2  3  0 ... 4  5  1  4  5 ]\\
	rowptr = [ 0  2  5  8 12 16 19 ]\\
\end{tt}

So, this is the class:

@d theClassHere
@{
template<typename T>
class SampleCRSMatrix {
public:
	typedef T value_type;
	@<constructors@>
	@<publicFunctions@>
private:
	@<privateFunctions@>
	@<privateData@>
}; // class SampleCRSMatrix
@}

From what we discussed above the private data is
@d privateData
@{
size_t rank_;
typename Vector<size_t>::Type rowptr_;
typename Vector<size_t>::Type colind_;
typename Vector<T>::Type values_;
@}

The constructor is
@d constructors
@{
@<emptyConstructor@>
@<randomConstructor@>
@<readFromAFile@>
@}

The empty constructor just creates an empty CRS matrix with the rank given:
@d emptyConstructor
@{
SampleCRSMatrix(size_t rank) : rank_(rank),rowptr_(rank+1)
{
}
@}

The random constructor, takes a seed, and the number of non-zero values and fills the matrix
randomly. Note that $T$ must be a real type here:
@d randomConstructor
@{
SampleCRSMatrix(size_t rank,T seed,size_t nonZeros,T maxValue) : rank_(rank),rowptr_(rank+1)
{
	srand48(seed);
	typename Vector<size_t>::Type rows,cols;
	typename Vector<T>::Type vals;
	for (size_t i=0;i<nonZeros;i++) {
		// pick a row
		size_t row = size_t(drand48()*rank);
		// and a column
		size_t col = size_t(drand48()*rank);
		// and a value
		T val = drand48()*maxValue;
		rows.push_back(row);
		cols.push_back(col);
		vals.push_back(val);
	}
	// fill the matrix with this data:
	fillMatrix(rows,cols,vals);
}
@}

This constructor reads the matrix from a file:
@d readFromAFile
@{
template<typename SomeIoInputType>
SampleCRSMatrix(SomeIoInputType& io)
{
	io>>rank_;
	readVector(io,rowptr_);
	readVector(io,colind_);
	readVector(io,values_);
}
@}
The public functions are
@d publicFunctions
@{
@<setRow@>
@<pushCol@>
@<pushValue@>
@<matrixVectorProduct@>
@<rank@>
@<save@>
@}

@d setRow
@{
void setRow(int n,int v)
{
	rowptr_[n]=v;
}
@}

The function below
performs $\vec{x} = \vec{x} + A * \vec{y}$,
 where $\vec{x}$ and $\vec{y}$ are vectors and A is this sparse matrix in
 compressed row format (CRS):
 @d matrixVectorProduct
 @{
void matrixVectorProduct(typename Vector<T>::Type& x, const typename Vector<T>::Type& y) const
{
	for (size_t i = 0; i < y.size(); i++)
		for (size_t j = rowptr_[i]; j < rowptr_[i + 1]; j++)
			x[i] += values_[j] * y[colind_[j]];
}
@}

@d rank
@{
size_t rank() const { return rank_; }
@}

@d pushCol
@{
void pushCol(int i) { colind_.push_back(i); }
@}

@d pushValue
@{
void pushValue(const T& value) { values_.push_back(value); }
@}

@d save
@{
template<typename SomeIoOutputType>
void save(SomeIoOutputType& io) const
{
	io<<rank_<<"\n";
	saveVector(io,rowptr_);
	saveVector(io,colind_);
	saveVector(io,values_);
}
@}

@d privateFunctions
@{
@<saveVector@>
@<readVector@>
@<fillMatrix@>
@}

@d saveVector
@{
template<typename SomeIoOutputType,typename SomeFieldType>
void saveVector(SomeIoOutputType& io,const typename Vector<SomeFieldType>::Type& v) const
{
	io<<v.size()<<"\n";
	for (size_t i=0;i<v.size();i++) {
		io<<v[i]<<" ";
	}
	io<<"\n";
}
@}

@d readVector
@{
template<typename SomeIoInputType,typename SomeFieldType>
void readVector(SomeIoInputType& io,typename Vector<SomeFieldType>::Type& v) const
{
	int size=0;
	io>>size;
	if (size<0) throw std::runtime_error("readVector: size is zero\n");
	v.resize(size);
	for (size_t i=0;i<v.size();i++) {
		io>>v[i];
	}
}
@}

@d fillMatrix
@{
void fillMatrix(typename Vector<size_t>::Type& rows,typename Vector<size_t>::Type& cols,
		typename Vector<T>::Type& vals)
{
	Sort<typename Vector<size_t>::Type > s;
    typename Vector<size_t>::Type iperm(rows.size());
	s.sort(rows,iperm);
	size_t counter = 0;
	size_t prevRow = rows[0]+1;
	for (size_t i=0;i<rows.size();i++) {
		size_t row = rows[i];
		if (prevRow!=row) {
			// add new row
			rowptr_[row] = counter++;
			prevRow = row;
		}
		colind_.push_back(cols[iperm[i]]);
		values_.push_back(vals[iperm[i]]);
	}
	size_t lastNonZeroRow = rows[rows.size()-1];
	for (size_t i=lastNonZeroRow+1;i<=rank_;i++)
		rowptr_[i] = counter;
}
@}

Finally an \verb=operator<<= to print this class:
@d companions
@{
@<operatorDump@>
@}

@d operatorDump
@{
template<typename T>
std::ostream& operator<<(std::ostream& os,const SampleCRSMatrix<T>& m)
{
	m.save(os);
	return os;
}
@}


\end{document}

