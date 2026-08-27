/*
Copyright (c) 2009-2015, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."

*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************

*/
/** \ingroup DMRG */
/*@{*/

/*! \file VectorWithOffsets.h
 * \brief Defines a vector stored as multiple symmetry-sector blocks.
 */
#ifndef VECTOR_WITH_OFFSETS_H
#define VECTOR_WITH_OFFSETS_H
#include "ProgramGlobals.h"
#include "Qn.h"
#include <PsimagLite/Complex.h>
#include <PsimagLite/ProgressIndicator.h>
#include <cassert>
#include <complex>
#include <typeinfo>

// FIXME: a more generic solution is needed instead of tying
// the non-zero structure to basis
namespace Dmrg {

/*!
 * \brief Stores a vector as independently addressable symmetry-sector blocks.
 *
 * The full vector is represented by one data block per basis sector. Only populated
 * sectors are listed in nzMsAndQns_; offsets_ maps each sector to its position in the
 * full vector, and index2Sector_ provides the inverse mapping for element access.
 *
 * Quantum numbers for populated sectors are represented by Dmrg::Qn.
 *
 * \tparam ComplexOrRealType Scalar type stored by the vector.
 */

// forward declaration for friendship
template <typename T> class OffsetVectorAny;
template <typename ComplexOrRealType> class VectorWithOffsets;

template <typename ComplexOrRealType>
typename PsimagLite::Real<ComplexOrRealType>::Type
norm(const VectorWithOffsets<ComplexOrRealType>& v);

template <typename ComplexOrRealType> void normalize(VectorWithOffsets<ComplexOrRealType>& v);

template <typename ComplexOrRealType>
ComplexOrRealType operator*(const VectorWithOffsets<ComplexOrRealType>& v1,
                            const VectorWithOffsets<ComplexOrRealType>& v2);

template <typename ComplexOrRealType>
VectorWithOffsets<ComplexOrRealType>
operator*(const typename VectorWithOffsets<ComplexOrRealType>::value_type& value,
          const VectorWithOffsets<ComplexOrRealType>&                      v);

template <typename ComplexOrRealType>
VectorWithOffsets<ComplexOrRealType> operator+(const VectorWithOffsets<ComplexOrRealType>& v1,
                                               const VectorWithOffsets<ComplexOrRealType>& v2);

template <typename ComplexOrRealType> class VectorWithOffsets {

	using ThisType       = VectorWithOffsets<ComplexOrRealType>;
	using VectorSizeType = typename Qn::VectorSizeType;
	using PairSizeType   = typename Qn::PairSizeType;

	static ComplexOrRealType const zero_;

public:

	using QnType           = Qn;
	using value_type       = ComplexOrRealType;
	using RealType         = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using PairQnType       = std::pair<SizeType, QnType>;
	using VectorType       = typename PsimagLite::Vector<ComplexOrRealType>::Type;
	using VectorVectorType = typename PsimagLite::Vector<VectorType>::Type;

	/*!
	 * \brief CONSTRUCTOR
	 *
	 * Constructs an empty vector with no populated sectors.
	 */
	VectorWithOffsets();

	/*!
	 * \brief Constructs sector storage from one weight per basis sector.
	 *
	 * \param[in] weights Number of stored elements in each sector.
	 * \param[in] someBasis Basis supplying partitions and quantum numbers.
	 */
	template <typename SomeBasisType>
	VectorWithOffsets(const typename PsimagLite::Vector<SizeType>::Type& weights,
	                  const SomeBasisType&                               someBasis)
	    : progress_("VectorWithOffsets")
	    , size_(someBasis.size())
	    , index2Sector_(size_)
	    , data_(weights.size())
	    , offsets_(weights.size() + 1)
	{
		for (SizeType i = 0; i < weights.size(); i++) {
			data_[i].resize(weights[i]);
			offsets_[i] = someBasis.partition(i);
			if (weights[i] > 0) {
				QnType qn = someBasis.pseudoQn(i);
				nzMsAndQns_.push_back(PairQnType(i, qn));
				// firstSector_ = i;
			}
		}

		offsets_[weights.size()] = size_;
		setIndex2Sector();
	}

	/*!
	 * \brief Rejects the single-sector construction path.
	 *
	 * This overload exists for interface compatibility but always reports a fatal error.
	 */
	template <typename SomeBasisType>
	VectorWithOffsets(SizeType weight, SizeType sector, const SomeBasisType& someBasis)
	    : progress_("VectorWithOffsets")
	    , size_(someBasis.size())
	    , index2Sector_(size_)
	    , data_(someBasis.partition() - 1)
	    , offsets_(someBasis.partition())
	{
		for (SizeType i = 0; i < offsets_.size(); ++i)
			offsets_[i] = someBasis.partition(i);

		assert(data_.size() < offsets_.size());
		if (offsets_[data_.size()] != size_)
			err("VectorWithOffsets::ctor(): FATAL: internal error\n");

		assert(sector < data_.size());
		data_[sector].resize(weight);
		QnType qn = someBasis.pseudoQn(sector);
		nzMsAndQns_.push_back(PairQnType(sector, qn));

		setIndex2Sector();
	}

	/*!
	 * \brief Constructs storage for an explicit list of populated sectors.
	 *
	 * \param[in] compactedWeights Number of stored elements for each listed sector.
	 * \param[in] sectors Basis-sector index for each compacted weight.
	 * \param[in] someBasis Basis supplying partitions and quantum numbers.
	 */
	template <typename SomeBasisType>
	VectorWithOffsets(const VectorSizeType& compactedWeights,
	                  const VectorSizeType& sectors,
	                  const SomeBasisType&  someBasis)
	    : progress_("VectorWithOffsets")
	    , size_(someBasis.size())
	    , index2Sector_(size_)
	    , data_(someBasis.partition() - 1)
	    , offsets_(someBasis.partition())
	{
		for (SizeType i = 0; i < offsets_.size(); ++i)
			offsets_[i] = someBasis.partition(i);

		assert(data_.size() < offsets_.size());
		if (offsets_[data_.size()] != size_)
			err("VectorWithOffsets::ctor(): FATAL: internal error\n");

		for (SizeType sectorIndex = 0; sectorIndex < compactedWeights.size();
		     ++sectorIndex) {
			const SizeType sector = sectors[sectorIndex];
			data_[sector].resize(compactedWeights[sectorIndex]);
			QnType qn = someBasis.pseudoQn(sector);
			nzMsAndQns_.push_back(PairQnType(sector, qn));
		}

		setIndex2Sector();
	}

	/*! \brief Resets the object to an empty vector with no sectors. */
	void clear();

	/*!
	 * \brief Replaces the vector with data from one sector.
	 *
	 * The data are swapped into this object, leaving \p v empty.
	 *
	 * \param[in,out] v Sector data transferred into this object.
	 * \param[in] sector Basis-sector index for the data.
	 * \param[in] someBasis Basis supplying partitions and quantum numbers.
	 */
	template <typename SomeBasisType>
	void set(VectorType& v, SizeType sector, const SomeBasisType& someBasis)
	{
		size_ = someBasis.size();
		nzMsAndQns_.clear();
		data_.clear();
		assert(someBasis.partition() > 0);
		SizeType n = someBasis.partition() - 1;
		data_.resize(n);
		offsets_.resize(n + 1);
		for (SizeType i = 0; i < n; ++i)
			offsets_[i] = someBasis.partition(i);

		data_[sector].swap(v);
		QnType qn = someBasis.pseudoQn(sector);
		nzMsAndQns_.push_back(PairQnType(sector, qn));

		offsets_[n] = size_;
		setIndex2Sector();
	}

	/*!
	 * \brief Allocates zero-filled storage for every basis sector.
	 * \param[in] someBasis Basis defining the sectors.
	 */
	template <typename SomeBasisType> void populateSectors(const SomeBasisType& someBasis)
	{
		SizeType np = someBasis.partition() - 1;
		size_       = someBasis.size();
		nzMsAndQns_.clear();
		data_.clear();
		data_.resize(np);
		offsets_.resize(np + 1);
		for (SizeType i = 0; i < np; i++) {
			offsets_[i]      = someBasis.partition(i);
			SizeType   total = someBasis.partition(i + 1) - offsets_[i];
			VectorType tmpV(total, 0);
			data_[i]         = tmpV;
			const QnType& qn = someBasis.pseudoQn(i);
			nzMsAndQns_.push_back(PairQnType(i, qn));
		}

		offsets_[np] = size_;
		setIndex2Sector();
		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "Populated " << np << " sectors";
		progress_.printline(msgg, std::cout);
	}

	/*!
	 * \brief Allocates zero-filled sectors matching another vector's quantum numbers.
	 *
	 * \param[in] v Vector whose populated quantum numbers are replicated.
	 * \param[in] someBasis Destination basis used to locate matching sectors.
	 */
	template <typename SomeBasisType>
	void populateFromQns(const VectorWithOffsets& v, const SomeBasisType& someBasis)
	{
		SizeType np = someBasis.partition() - 1;
		size_       = someBasis.size();
		nzMsAndQns_.clear();
		data_.clear();
		data_.resize(np);
		offsets_.resize(np + 1);
		for (SizeType i = 0; i < np; i++) {
			offsets_[i] = someBasis.partition(i);
		}

		offsets_[np] = size_;

		for (SizeType i = 0; i < v.sectors(); i++) {
			SizeType ip = findPartitionWithThisQn(v.qn(i), someBasis);
			if (ip >= np)
				err("VectorWithOffsets: populateFromQns\n");
			SizeType   total = someBasis.partition(ip + 1) - offsets_[ip];
			VectorType tmpV(total, 0);
			data_[ip] = tmpV;
			nzMsAndQns_.push_back(PairQnType(ip, v.qn(i)));
		}

		setIndex2Sector();
		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "populateFromQns " << v.sectors() << " sectors";
		progress_.printline(msgg, std::cout);
	}

	/*!
	 * \brief Removes sectors whose stored elements are numerically zero.
	 *
	 * All basis sectors must be represented before this operation is called.
	 */
	void collapseSectors();

	/*!
	 * \brief Replaces the data at a basis-sector index.
	 * \param[in] v Replacement data.
	 * \param[in] i0 Basis-sector index in the sector array.
	 */
	void setDataInSector(const VectorType& v, SizeType i0);

	/*! \brief Returns the number of populated sectors. */
	SizeType sectors() const { return nzMsAndQns_.size(); }

	/*!
	 * \brief Returns the basis-sector index of a populated sector.
	 * \param[in] i Index in the compact list of populated sectors.
	 */
	SizeType sector(SizeType i) const;

	/*!
	 * \brief Returns the quantum number of a populated sector.
	 * \param[in] i Index in the compact list of populated sectors.
	 */
	const QnType& qn(SizeType i) const;

	/*!
	 * \brief Builds the sector representation from a full vector.
	 *
	 * Only sectors containing nonzero elements are stored as populated.
	 *
	 * \param[in] v Full vector in basis order.
	 * \param[in] someBasis Basis defining sector boundaries and quantum numbers.
	 */
	template <typename SomeBasisType>
	void fromFull(const VectorType& v, const SomeBasisType& someBasis)
	{
		size_ = someBasis.size();

		offsets_.resize(someBasis.partition());
		for (SizeType i = 0; i < someBasis.partition(); i++)
			offsets_[i] = someBasis.partition(i);
		assert(offsets_[offsets_.size() - 1] == size_);

		data_.clear();
		data_.resize(someBasis.partition() - 1);

		nzMsAndQns_.clear();
		findPartitions(nzMsAndQns_, v, someBasis);
		for (SizeType jj = 0; jj < nzMsAndQns_.size(); jj++) {
			SizeType j = nzMsAndQns_[jj].first;
			assert(j < data_.size());
			SizeType offset = offsets_[j];
			SizeType total  = offsets_[j + 1] - offset;
			data_[j].resize(total);
			for (SizeType i = 0; i < total; i++)
				data_[j][i] = v[i + offset];
		}

		setIndex2Sector();
	}

	/*!
	 * \brief Copies the data stored at a basis-sector index.
	 * \param[out] v Receives the sector data.
	 * \param[in] i Basis-sector index.
	 */
	void extract(VectorType& v, SizeType i) const;

	/*! \brief Returns the size of the corresponding full vector. */
	SizeType size() const { return size_; }

	/*!
	 * \brief Returns the number of stored elements at a basis-sector index.
	 * \param[in] i Basis-sector index.
	 */
	SizeType effectiveSize(SizeType i) const;

	/*!
	 * \brief Returns a sector's starting position in the full vector.
	 * \param[in] i Basis-sector index.
	 */
	SizeType offset(SizeType i) const;

	/*!
	 * \brief Returns a const reference to an element using sector-local indices.
	 * \param[in] i Basis-sector index.
	 * \param[in] j Element index within the sector.
	 */
	const ComplexOrRealType& fastAccess(SizeType i, SizeType j) const;

	/*!
	 * \brief Returns a mutable reference to an element using sector-local indices.
	 * \param[in] i Basis-sector index.
	 * \param[in] j Element index within the sector.
	 */
	ComplexOrRealType& fastAccess(SizeType i, SizeType j);

	PairSizeType sectorAndOffset() const;

	/*!
	 * \brief Copies populated sector data into a full sparse-vector representation.
	 * \param[out] sv Destination sparse vector, resized to size().
	 */
	template <typename SparseVectorType> void toSparse(SparseVectorType& sv) const
	{
		sv.resize(size_);
		for (SizeType jj = 0; jj < nzMsAndQns_.size(); jj++) {
			SizeType j = nzMsAndQns_[jj].first;
			assert(j < data_.size());
			for (SizeType i = 0; i < data_[j].size(); i++)
				sv[i + offsets_[j]] = data_[j][i];
		}
	}

	/*!
	 * \brief Reads the vector from hierarchical storage.
	 * \param[in,out] io Input object.
	 * \param[in] label Storage path containing the vector.
	 */
	template <typename SomeInputType> void read(SomeInputType& io, PsimagLite::String label)
	{
		io.read(size_, label + "/size_");
		if (size_ == 0)
			return;
		io.read(index2Sector_, label + "/index2Sector_");
		SizeType x = 0;
		io.read(x, label + "/data_/Size");
		data_.resize(x);
		bool flag = false;
		for (SizeType i = 0; i < x; ++i) {
			try {
				io.read(data_[i], label + "/data_/" + ttos(i));
				flag = true;
				std::cerr << "VectorWithOffsets: non-zero sector index " << i
				          << " read \n";
			} catch (...) { }
		}

		if (!flag)
			err("VectorWithOffsets: all sectors in data_ are empty (FATAL)\n");

		io.read(offsets_, label + "/offsets_");
		SizeType aSize = 0;
		io.read(aSize, label + "/nzMsAndQns_/Size");
		nzMsAndQns_.resize(aSize, PairQnType(0, QnType::zero()));
		for (SizeType i = 0; i < aSize; ++i) {
			io.read(nzMsAndQns_[i].first, label + "/nzMsAndQns_/" + ttos(i) + "/0");
			nzMsAndQns_[i].second.read(label + "/nzMsAndQns_/" + ttos(i) + "/1", io);
		}
	}

	/*!
	 * \brief Writes the vector to hierarchical storage.
	 * \param[in,out] io Output object.
	 * \param[in] label Storage path for the vector.
	 */
	template <typename SomeIoOutputType>
	void write(SomeIoOutputType& io, const PsimagLite::String& label) const
	{
		io.createGroup(label);
		io.write(size_, label + "/size_");
		io.write(index2Sector_, label + "/index2Sector_");
		io.write(data_, label + "/data_");
		io.write(offsets_, label + "/offsets_");
		io.write(nzMsAndQns_, label + "/nzMsAndQns_");
	}

	// We don't have a partitioned basis because we don't
	// have the superblock basis at this point
	// Therefore, partitioning is bogus here
	/*!
	 * \brief Loads a legacy text representation containing populated sectors.
	 * \param[in,out] io Input object.
	 * \param[in] label Label advanced to before reading.
	 * \param[in] counter Occurrence of \p label to load.
	 */
	template <typename IoInputter>
	void loadOneSector(IoInputter& io, const PsimagLite::String& label, SizeType counter = 0)
	{
		PsimagLite::String msg("VectorWithOffsets:");
		io.advance(label, counter);
		int x = 0;
		io.readline(x, "size=");
		if (x < 0)
			err(msg + ":loadOneSector(...): size<0\n");
		size_ = x;

		io.read(offsets_, "offsets");

		data_.clear();
		data_.resize(offsets_.size());

		io.readline(x, "nonzero=");
		if (x < 0)
			err(msg + ":loadOneSector(...): nonzerosectors<0\n");
		const SizeType nonzero = x;
		nzMsAndQns_.clear();
		nzMsAndQns_.reserve(nonzero);

		for (SizeType jj = 0; jj < nonzero; jj++) {
			io.readline(x, "sector=");
			if (x < 0)
				err(msg + ":loadOneSector(...): sector<0\n");

			QnType y = QnType::zero();
			y.read("qn", io);
			nzMsAndQns_.push_back(PairQnType(x, y));

			if (static_cast<SizeType>(x) >= data_.size())
				err(msg + ":loadOneSector(...): sector too big\n");

			PsimagLite::String s = "data" + ttos(jj);
			io.read(data_[x], s);
		}

		setIndex2Sector();
	}

	/*!
	 * \brief Scales every populated sector in place.
	 * \param[in] value Scalar multiplier.
	 * \return This vector after scaling.
	 */
	VectorWithOffsets& operator*=(const ComplexOrRealType& value);

	/*!
	 * \brief Adds another vector's populated sector data in place.
	 * \param[in] v Vector to add.
	 * \return This vector after addition.
	 */
	VectorWithOffsets operator+=(const VectorWithOffsets& v);

	/*!
	 * \brief Maps a full-vector index to its basis-sector index.
	 * \param[in] i Full-vector index.
	 * \return Basis-sector index, or -1 for an unpopulated sector.
	 */
	int index2Sector(SizeType i) const;

	/*! \brief Returns the serialization name for this vector type. */
	static PsimagLite::String name() { return "vectorwithoffsets"; }

	/*!
	 * \brief Computes the Euclidean norm across all populated sectors.
	 * \param[in] v Vector to measure.
	 */
	template <typename T>
	friend typename PsimagLite::Real<T>::Type norm(const VectorWithOffsets<T>& v);

	/*!
	 * \brief Normalizes a vector in place.
	 * \param[in,out] v Vector to normalize.
	 */
	template <typename T> friend void normalize(VectorWithOffsets<T>& v);

	/*!
	 * \brief Computes the inner product over sectors populated in both vectors.
	 * \param[in] v1 Left vector.
	 * \param[in] v2 Right vector.
	 */
	template <typename T>
	friend T operator*(const VectorWithOffsets<T>& v1, const VectorWithOffsets<T>& v2);

	/*!
	 * \brief Returns a scaled copy of a vector.
	 * \param[in] value Scalar multiplier.
	 * \param[in] v Vector to scale.
	 */
	template <typename T>
	friend VectorWithOffsets<T>
	operator*(const typename VectorWithOffsets<T>::value_type& value,
	          const VectorWithOffsets<T>&                      v);

	/*!
	 * \brief Returns the sum of vectors with matching populated sectors.
	 * \param[in] v1 Left vector.
	 * \param[in] v2 Right vector.
	 */
	template <typename T>
	friend VectorWithOffsets<T> operator+(const VectorWithOffsets<T>& v1,
	                                      const VectorWithOffsets<T>& v2);

	friend class OffsetVectorAny<ComplexOrRealType>;

private:

	// Don't use directly; Use via OffsetVector which is a friend
	const ComplexOrRealType& slowAccess(SizeType i) const;

	void setIndex2Sector();

	template <typename SomeBasisType>
	void findPartitions(typename PsimagLite::Vector<PairQnType>::Type& p,
	                    const VectorType&                              v,
	                    const SomeBasisType&                           someBasis)
	{
		bool found = false;
		p.clear();
		for (SizeType i = 0; i < someBasis.partition() - 1; i++) {
			if (nonZeroPartition(v, someBasis, i)) {
				found            = true;
				const QnType& qn = someBasis.pseudoQn(i);
				p.push_back(PairQnType(i, qn));
			}
		}

		if (!found) {
			PsimagLite::OstringStream                     msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg << "No partition found";
			progress_.printline(msgg, std::cout);
		}
	}

	template <typename SomeBasisType>
	bool nonZeroPartition(const VectorType& v, const SomeBasisType& someBasis, SizeType i)
	{
		typename VectorType::value_type zero = 0;
		for (SizeType j = someBasis.partition(i); j < someBasis.partition(i + 1); ++j) {
			assert(j < v.size());
			if (v[j] != zero)
				return true;
		}

		return false;
	}

	bool isZero(const VectorType& v) const;

	template <typename SomeBasisType>
	SizeType findPartitionWithThisQn(const QnType& qn, const SomeBasisType& someBasis) const
	{
		SizeType np = someBasis.partition() - 1;
		for (SizeType i = 0; i < np; ++i)
			if (someBasis.qnEx(i) == qn)
				return i;

		throw PsimagLite::RuntimeError("findPartitionWithThisQn\n");
	}

	PsimagLite::ProgressIndicator                 progress_;
	SizeType                                      size_;
	typename PsimagLite::Vector<int>::Type        index2Sector_;
	VectorVectorType                              data_;
	typename PsimagLite::Vector<SizeType>::Type   offsets_;
	typename PsimagLite::Vector<PairQnType>::Type nzMsAndQns_;
}; // class VectorWithOffset

extern template class VectorWithOffsets<double>;
extern template class VectorWithOffsets<std::complex<double>>;

extern template double               norm(const VectorWithOffsets<double>& v);
extern template double               norm(const VectorWithOffsets<std::complex<double>>& v);
extern template void                 normalize(VectorWithOffsets<double>& v);
extern template void                 normalize(VectorWithOffsets<std::complex<double>>& v);
extern template double               operator*(const VectorWithOffsets<double>& v1,
                                 const VectorWithOffsets<double>& v2);
extern template std::complex<double> operator*(const VectorWithOffsets<std::complex<double>>& v1,
                                               const VectorWithOffsets<std::complex<double>>& v2);
extern template VectorWithOffsets<double> operator*(const double&                    value,
                                                    const VectorWithOffsets<double>& v);
extern template VectorWithOffsets<std::complex<double>>
operator*(const std::complex<double>& value, const VectorWithOffsets<std::complex<double>>& v);
extern template VectorWithOffsets<double> operator+(const VectorWithOffsets<double>& v1,
                                                    const VectorWithOffsets<double>& v2);
extern template VectorWithOffsets<std::complex<double>>
operator+(const VectorWithOffsets<std::complex<double>>& v1,
          const VectorWithOffsets<std::complex<double>>& v2);
}
/*@}*/
#endif
