#ifndef RestartStruct_H
#define RestartStruct_H

#include "Vector.h"
#include <iostream>
#include "Io/IoSerializerStub.h"
#include "InputNg.h"

namespace Dmrg {

// no longer a struct
struct RestartStruct {

	typedef PsimagLite::Vector<int>::Type VectorIntType;

	RestartStruct()
	    : filename_(""), labelForEnergy_("Energies"), mapStages_(true), sourceTvForPsi_(-1)
	{}

	/* PSIDOC MiscRestartOptions

	 RestartSourceTvForPsi=integer
	 Optional. If present and non-negative, set the g.s. of this targeting to the
	 previous run target vector indicated by this number.

	 RestartMapStages=1 or 0
	 Optional. If present and 0, do not set stages for this targeting. If absent, or
	 present and 1, set stages for this targeting to stages of previous run.

	 RestartMappingTvs (vector of integers)
	 Optional. If present the number of numbers to follow must be provided first,
	 as for all vectors in InputNg legacy, and with the appropriate syntax for Ainur.
	 Sets the i-th target vector for this targeting to the previous run target
	 vector number MappingTvs[i], if MappingTvs[i] is non-negative; skips it if negative.
	 */
	template<typename SomeInputType>
	void read(SomeInputType& io)
	{
		try {
			io.readline(labelForEnergy_, "RestartLabelForEnergy=");
		} catch (std::exception&) {}

		try {
			io.readline(sourceTvForPsi_, "RestartSourceTvForPsi=");
		} catch (std::exception&) {
			sourceTvForPsi_ = -1;
		}

		try {
			int x = 1;
			io.readline(x, "RestartMapStages=");
			mapStages_ = (x > 0);
		} catch (std::exception&) {}

		try {
			io.read(mappingTvs_, "RestartMappingTvs");
		} catch (std::exception&) {}
	}

	void setFilename(PsimagLite::String f) { filename_ = f; }

	PsimagLite::String filename() const { return filename_; }

	PsimagLite::String labelForEnergy() const { return labelForEnergy_; }

	bool mapStages() const { return mapStages_; }

	int mappingTvs(SizeType ind) const
	{
		if (mappingTvs_.size() == 0)
			return ind;

		if (ind >= mappingTvs_.size())
			err("RestartStruct::mappingTvs not provided for ind= " + ttos(ind) + "\n");

		return mappingTvs_[ind];
	}

	int sourceTvForPsi() const { return sourceTvForPsi_; }

	void write(PsimagLite::String label,
	           PsimagLite::IoSerializer& ioSerializer) const
	{
		PsimagLite::String root = label;
		ioSerializer.createGroup(root);
		ioSerializer.write(root + "/filename", filename_);
		ioSerializer.write(root + "/labelForEnergy", labelForEnergy_);
		ioSerializer.write(root + "/mapStages", mapStages_);
		ioSerializer.write(root + "/sourceTvForPsi", sourceTvForPsi_);
		if (mappingTvs_.size() > 0)
			ioSerializer.write(root + "/mappingTvs", mappingTvs_);
	}

	friend std::ostream& operator<<(std::ostream& os, const RestartStruct& c)
	{
	    if (c.filename_ == "") return os;

	    os<<"RestartStruct.filename="<<c.filename_<<"\n";
	    os<<"RestartStruct.labelForEnergy="<<c.labelForEnergy_<<"\n";
		os<<"RestartStruct.mapStages="<<c.mapStages_<<"\n";
		os<<"RestartStruct.sourceTvForPsi="<<c.sourceTvForPsi_<<"\n";
		if (c.mappingTvs_.size() > 0)
			os<<"RestartStruct.mappingTvs="<<c.mappingTvs_<<"\n";
	    return os;
	}

private:

	PsimagLite::String filename_;
	PsimagLite::String labelForEnergy_;
	bool mapStages_;
	int sourceTvForPsi_;
	VectorIntType mappingTvs_;
};

} // namespace Dmrg

#endif // RestartStruct_H

