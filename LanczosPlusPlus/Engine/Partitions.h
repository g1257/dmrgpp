
/*
// BEGIN LICENSE BLOCK
Copyright (c) 2009-2012, UT-Battelle, LLC
All rights reserved

[Lanczos++, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/

#ifndef PARTITIONS_H
#define PARTITIONS_H
#include "BitManip.h"
#include "Matrix.h"

namespace LanczosPlusPlus {

class Partitions {

public:

	Partitions(SizeType length, SizeType parts)
	    : length_(length)
	{
		PsimagLite::Vector<SizeType>::Type values(parts, 0);

		while (true) {
			if (sumOf(values) == length_)
				partitions_.push_back(values);
			values[0]++;
			if (sumOf(values) > length_) {
				int x = increaseNextIndices(values);
				if (x < 0)
					break;
			}
		}
	}

	SizeType size() const { return partitions_.size(); }

	const PsimagLite::Vector<SizeType>::Type& operator()(SizeType i) const
	{
		return partitions_[i];
	}

private:

	SizeType sumOf(const PsimagLite::Vector<SizeType>::Type& values) const
	{
		SizeType sum = 0;
		for (SizeType i = 0; i < values.size(); i++)
			sum += values[i];
		return sum;
	}

	int increaseNextIndices(PsimagLite::Vector<SizeType>::Type& values) const
	{
		if (0 == values.size() - 1)
			return -1;
		values[0]  = 0;
		SizeType i = 1;
		while (true) {
			values[i]++;
			if (sumOf(values) <= length_)
				return i;
			if (i == values.size() - 1)
				return -1;
			values[i] = 0;
			i++;
		}
	}

	SizeType                                                     length_;
	PsimagLite::Vector<PsimagLite::Vector<SizeType>::Type>::Type partitions_;

}; // class Partitions

} // namespace LanczosPlusPlus
#endif
