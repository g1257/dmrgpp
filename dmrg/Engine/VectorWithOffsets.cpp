#include "VectorWithOffsets.h"

#include <cmath>
#include <complex>

namespace Dmrg {

template <typename ComplexOrRealType>
VectorWithOffsets<ComplexOrRealType>::VectorWithOffsets()
    : progress_("VectorWithOffsets")
    , size_(0)
    , index2Sector_(0)
{
}

template <typename ComplexOrRealType>
void VectorWithOffsets<ComplexOrRealType>::clear()
{
	size_ = 0;
	index2Sector_.clear();
	data_.clear();
	offsets_.clear();
	nzMsAndQns_.clear();
}

template <typename ComplexOrRealType>
void VectorWithOffsets<ComplexOrRealType>::collapseSectors()
{
	SizeType np = data_.size();
	if (np != nzMsAndQns_.size()) {
		PsimagLite::String str("VectorWithOffsets: collapseSectors cannot be called");
		err(str + " on a partially populated vector\n");
	}

	typename PsimagLite::Vector<PairQnType>::Type nzMsAndQns;
	for (SizeType i = 0; i < np; ++i) {
		if (isZero(data_[i])) {
			data_[i].resize(0);
		} else {
			assert(i < nzMsAndQns_.size());
			nzMsAndQns.push_back(nzMsAndQns_[i]);
		}
	}

	nzMsAndQns_ = nzMsAndQns;
	setIndex2Sector();
	PsimagLite::OstringStream                     msgg(std::cout.precision());
	PsimagLite::OstringStream::OstringStreamType& msg = msgg();
	msg << "Collapsed. Non-zero sectors now are " << nzMsAndQns_.size();
	progress_.printline(msgg, std::cout);
}

template <typename ComplexOrRealType>
void VectorWithOffsets<ComplexOrRealType>::setDataInSector(const VectorType& v,
                                                            SizeType         i0)
{
	if (i0 >= data_.size())
		err("VectorWithOffsets: setDataInSector\n");

	data_[i0] = v;
}

template <typename ComplexOrRealType>
SizeType VectorWithOffsets<ComplexOrRealType>::sector(SizeType i) const
{
	assert(i < nzMsAndQns_.size());
	return nzMsAndQns_[i].first;
}

template <typename ComplexOrRealType>
const typename VectorWithOffsets<ComplexOrRealType>::QnType&
VectorWithOffsets<ComplexOrRealType>::qn(SizeType i) const
{
	assert(i < nzMsAndQns_.size());
	return nzMsAndQns_[i].second;
}

template <typename ComplexOrRealType>
void VectorWithOffsets<ComplexOrRealType>::extract(VectorType& v, SizeType i) const
{
	if (i >= data_.size())
		err("VectorWithOffsets: extract\n");

	v = data_[i];
}

template <typename ComplexOrRealType>
SizeType VectorWithOffsets<ComplexOrRealType>::effectiveSize(SizeType i) const
{
	if (i >= data_.size())
		err("VectorWithOffsets: effectiveSize\n");

	return data_[i].size();
}

template <typename ComplexOrRealType>
SizeType VectorWithOffsets<ComplexOrRealType>::offset(SizeType i) const
{
	if (i >= offsets_.size())
		err("VectorWithOffsets: offset\n");

	return offsets_[i];
}

template <typename ComplexOrRealType>
const ComplexOrRealType& VectorWithOffsets<ComplexOrRealType>::fastAccess(SizeType i,
                                                                          SizeType j) const
{
	assert(i < data_.size());
	assert(j < data_[i].size());
	return data_[i][j];
}

template <typename ComplexOrRealType>
ComplexOrRealType& VectorWithOffsets<ComplexOrRealType>::fastAccess(SizeType i,
                                                                    SizeType j)
{
	assert(i < data_.size());
	assert(j < data_[i].size());
	return data_[i][j];
}

template <typename ComplexOrRealType>
typename VectorWithOffsets<ComplexOrRealType>::PairSizeType
VectorWithOffsets<ComplexOrRealType>::sectorAndOffset() const
{
	if (nzMsAndQns_.size() != 1) {
		err("sectorAndOffset() cannot be called if there are more than one "
		    "non trivial sector\n");
	}

	SizeType sector = nzMsAndQns_[0].first;
	assert(sector < offsets_.size());
	return PairSizeType(sector, offsets_[sector]);
}

template <typename ComplexOrRealType>
VectorWithOffsets<ComplexOrRealType>& VectorWithOffsets<ComplexOrRealType>::operator*=(
    const ComplexOrRealType& value)
{
	for (SizeType ii = 0; ii < nzMsAndQns_.size(); ++ii) {
		SizeType i = nzMsAndQns_[ii].first;
		assert(i < data_.size());
		data_[i] *= value;
	}

	return *this;
}

template <typename ComplexOrRealType>
VectorWithOffsets<ComplexOrRealType> VectorWithOffsets<ComplexOrRealType>::operator+=(
    const VectorWithOffsets& v)
{
	if (nzMsAndQns_.size() == 0) {
		size_       = v.size_;
		data_       = v.data_;
		offsets_    = v.offsets_;
		nzMsAndQns_ = v.nzMsAndQns_;
		setIndex2Sector();
		return *this;
	}

	for (SizeType ii = 0; ii < nzMsAndQns_.size(); ii++) {
		SizeType i = nzMsAndQns_[ii].first;
		assert(i < data_.size());
		data_[i] += v.data_[i];
	}

	setIndex2Sector();
	return *this;
}

template <typename ComplexOrRealType>
int VectorWithOffsets<ComplexOrRealType>::index2Sector(SizeType i) const
{
	assert(i < index2Sector_.size());
	return index2Sector_[i];
}

template <typename ComplexOrRealType>
typename PsimagLite::Real<ComplexOrRealType>::Type norm(
    const VectorWithOffsets<ComplexOrRealType>& v)
{
	using RealType = typename PsimagLite::Real<ComplexOrRealType>::Type;
	RealType sum    = 0;
	for (SizeType ii = 0; ii < v.nzMsAndQns_.size(); ii++) {
		SizeType i = v.nzMsAndQns_[ii].first;
		assert(i < v.data_.size());
		RealType tmp = PsimagLite::norm(v.data_[i]);
		sum += tmp * tmp;
	}

	return std::sqrt(sum);
}

template <typename ComplexOrRealType>
void normalize(VectorWithOffsets<ComplexOrRealType>& v)
{
	using RealType = typename PsimagLite::Real<ComplexOrRealType>::Type;
	RealType norma = norm(v);
	RealType eps   = 1e-5;

	if (std::fabs(norma - 1.0) < eps)
		return;

	assert(std::fabs(norma) > eps);

	for (SizeType i = 0; i < v.data_.size(); i++)
		for (SizeType j = 0; j < v.data_[i].size(); j++)
			v.data_[i][j] /= norma;
}

template <typename ComplexOrRealType>
ComplexOrRealType operator*(const VectorWithOffsets<ComplexOrRealType>& v1,
                            const VectorWithOffsets<ComplexOrRealType>& v2)
{
	ComplexOrRealType sum = 0;
	for (SizeType ii = 0; ii < v1.sectors(); ++ii) {
		SizeType i = v1.sector(ii);
		for (SizeType jj = 0; jj < v1.sectors(); ++jj) {
			SizeType j = v2.sector(jj);
			if (i != j)
				continue;
			for (SizeType k = 0; k < v1.effectiveSize(i); ++k)
				sum += v1.fastAccess(i, k) * PsimagLite::conj(v2.fastAccess(j, k));
		}
	}

	return sum;
}

template <typename ComplexOrRealType>
VectorWithOffsets<ComplexOrRealType> operator*(
    const ComplexOrRealType&                    value,
    const VectorWithOffsets<ComplexOrRealType>& v)
{
	VectorWithOffsets<ComplexOrRealType> w = v;

	for (SizeType ii = 0; ii < w.nzMsAndQns_.size(); ++ii) {
		SizeType i = w.nzMsAndQns_[ii].first;
		assert(i < w.data_.size());
		w.data_[i] *= value;
	}

	return w;
}

template <typename ComplexOrRealType>
VectorWithOffsets<ComplexOrRealType> operator+(
    const VectorWithOffsets<ComplexOrRealType>& v1,
    const VectorWithOffsets<ComplexOrRealType>& v2)
{
	PsimagLite::String s = "VectorWithOffsets + VectorWithOffsets failed\n";
	if (v1.nzMsAndQns_ != v2.nzMsAndQns_)
		err(s.c_str());

	for (SizeType ii = 0; ii < v1.nzMsAndQns_.size(); ii++) {
		SizeType i = v1.nzMsAndQns_[ii].first;
		if (i >= v1.data_.size() || i >= v2.data_.size())
			err(s.c_str());
		if (v1.data_[i].size() != v2.data_[i].size())
			err(s.c_str());
	}

	VectorWithOffsets<ComplexOrRealType> w = v1;
	w += v2;
	return w;
}

template <typename ComplexOrRealType>
const ComplexOrRealType& VectorWithOffsets<ComplexOrRealType>::slowAccess(SizeType i) const
{
	assert(i < index2Sector_.size());
	int j = index2Sector_[i];
	if (j < 0)
		return zero_;
	return data_[j][i - offsets_[j]];
}

template <typename ComplexOrRealType>
void VectorWithOffsets<ComplexOrRealType>::setIndex2Sector()
{
	if (index2Sector_.size() != size_)
		index2Sector_.resize(size_);

	for (SizeType i = 0; i < size_; ++i) {
		index2Sector_[i] = -1;
		for (SizeType jj = 0; jj < nzMsAndQns_.size(); ++jj) {
			SizeType j = nzMsAndQns_[jj].first;
			assert(j + 1 < offsets_.size());
			if (i < offsets_[j] || i >= offsets_[j + 1])
				continue;
			index2Sector_[i] = j;
		}
	}
}

template <typename ComplexOrRealType>
bool VectorWithOffsets<ComplexOrRealType>::isZero(const VectorType& v) const
{
	RealType eps = 1e-5;
	for (SizeType i = 0; i < v.size(); ++i) {
		if (std::fabs(PsimagLite::real(v[i])) > eps
		    || std::fabs(PsimagLite::imag(v[i])) > eps)
			return false;
	}

	return true;
}

template <typename ComplexOrRealType>
const ComplexOrRealType VectorWithOffsets<ComplexOrRealType>::zero_ = 0;

template class VectorWithOffsets<double>;
template class VectorWithOffsets<std::complex<double>>;

template double norm(const VectorWithOffsets<double>& v);
template double norm(const VectorWithOffsets<std::complex<double>>& v);
template void normalize(VectorWithOffsets<double>& v);
template void normalize(VectorWithOffsets<std::complex<double>>& v);
template double operator*(const VectorWithOffsets<double>& v1,
                          const VectorWithOffsets<double>& v2);
template std::complex<double> operator*(
    const VectorWithOffsets<std::complex<double>>& v1,
    const VectorWithOffsets<std::complex<double>>& v2);
template VectorWithOffsets<double> operator*(const double&                     value,
                                             const VectorWithOffsets<double>& v);
template VectorWithOffsets<std::complex<double>> operator*(
    const std::complex<double>&                    value,
    const VectorWithOffsets<std::complex<double>>& v);
template VectorWithOffsets<double> operator+(const VectorWithOffsets<double>& v1,
                                             const VectorWithOffsets<double>& v2);
template VectorWithOffsets<std::complex<double>> operator+(
    const VectorWithOffsets<std::complex<double>>& v1,
    const VectorWithOffsets<std::complex<double>>& v2);

} // namespace Dmrg
