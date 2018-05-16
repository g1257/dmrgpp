
/*
// BEGIN LICENSE BLOCK
Copyright (c) 2009 , UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/
/** \ingroup DMRG */
/*@{*/

/*! \file ContinuedFractionCollection.h
 *
 * We need many cont. fractions for the Green's function
 * on different sites, because, you know, we need
 * c_i + c_j and also c_i - c_j, and so on.
 * This class handles the composition
 */

#ifndef CONTINUED_FRACTION_COLL_H
#define CONTINUED_FRACTION_COLL_H
#include <iostream>
#include "TypeToString.h"
#include "ProgressIndicator.h"
#include "FreqEnum.h"

namespace PsimagLite {

template<typename ContinuedFractionType_>
class ContinuedFractionCollection  {
public:

	typedef ContinuedFractionType_ ContinuedFractionType;
	typedef typename ContinuedFractionType::ComplexType ComplexType;
	typedef typename ContinuedFractionType::TridiagonalMatrixType
	TridiagonalMatrixType;
	typedef typename	TridiagonalMatrixType::value_type RealType;
	typedef typename ContinuedFractionType::MatrixType MatrixType;
	typedef typename ContinuedFractionType::PlotDataType PlotDataType;
	typedef typename ContinuedFractionType::PlotParamsType PlotParamsType;

	ContinuedFractionCollection(FreqEnum freqEnum)
	    : freqEnum_(freqEnum), progress_("ContinuedFractionCollection")
	{
	}

	template<typename IoInputType>
	ContinuedFractionCollection(IoInputType& io,SizeType level = 0)
	    : freqEnum_(FREQ_REAL),progress_("ContinuedFractionCollection")
	{
		int n = 0;
		io.readline(n,"#CONTINUEDFRACTIONCOLLECTION=",level);
		if (n<=0) {
			String s = "ContinuedFractionCollection::ctor(...): ";
			s += "Expected a positive number of items, got " +
			        ttos(n);
			throw RuntimeError(s.c_str());
		}
		for (SizeType i=0;i<SizeType(n);i++) {
			ContinuedFractionType cf(io);
			data_.push_back(cf);
			freqEnum_ = cf.freqType();
		}
	}

	template<typename IoOutputType>
	void write(IoOutputType& io) const
	{
		io.write(data_.size(), "#CONTINUEDFRACTIONCOLLECTION");
		for (SizeType i=0;i<data_.size();i++)
			data_[i].write(io, "");
	}

	void push(const ContinuedFractionType& cf)
	{
		data_.push_back(cf);
	}

	void plot(PlotDataType& result,
	          const PlotParamsType& params) const
	{
		for (SizeType i=0;i<data_.size();i++) {
			PlotDataType result1;
			data_[i].plot(result1,params);
			accumulate(result, result1);
		}
	}

	void plotOne(
	        SizeType i,
	        PlotDataType& result,
	        const PlotParamsType& params) const
	{
		data_[i].plot(result,params);
	}

	SizeType size() const { return data_.size(); }

	FreqEnum freqType() const  { return freqEnum_; }

private:

	void accumulate(PlotDataType& v1,const PlotDataType& v2) const
	{
		bool wasEmpty = false;
		if (v1.size()==0) {
			wasEmpty = true;
			v1.resize(v2.size());
		} else {
			if (v1.size()!=v2.size()) {
				String s = "ContinuedFractionCollection::acc...(...)";
				s += " vectors must be of same length\n";
				throw RuntimeError(s.c_str());
			}
		}
		for (SizeType i=0;i<v1.size();i++) {

			if (wasEmpty) {
				v1[i].first = v2[i].first;
				v1[i].second = v2[i].second;
			} else {
				if (v1[i].first!=v2[i].first)
					throw RuntimeError("CF: x coordinate different\n");
				v1[i].second += v2[i].second;
			}
		}
	}

	FreqEnum freqEnum_;
	ProgressIndicator progress_;
	typename Vector<ContinuedFractionType>::Type data_;
}; // class ContinuedFractionCollection
} // namespace PsimagLite
/*@}*/
#endif  //CONTINUED_FRACTION_COLL_H

