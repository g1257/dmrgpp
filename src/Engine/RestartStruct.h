#ifndef RestartStruct_H
#define RestartStruct_H

#include "Vector.h"
#include <iostream>
#include "Io/IoSerializerStub.h"
#include "InputNg.h"

namespace Dmrg {

// no longer a struct
struct RestartStruct {

	RestartStruct()
	    : filename_(""),labelForEnergy_("Energy")
	{}

	template<typename SomeInputType>
	void read(SomeInputType& io)
	{
		try {
			io.readline(labelForEnergy_, "RestartLabelForEnergy=");
		} catch (std::exception&) {}
	}

	void setFilename(PsimagLite::String f) { filename_ = f; }

	PsimagLite::String filename() const { return filename_; }

	PsimagLite::String labelForEnergy() const { return labelForEnergy_; }

	SizeType mappingStages(SizeType ind) const
	{
		return ind;
	}

	int mappingTvs(SizeType ind) const
	{
		//if (mappingTvs_.size() == 0)
			return ind;

		//if (ind >= mappingTvs_.size())
		//	err("RestartStruct::mappingTvs not provided for ind= " + ttos(ind) + "\n");

		//return mappingTvs_[ind];
	}

	int sourceTvForPsi() const { return -1; }

	void write(PsimagLite::String label,
	           PsimagLite::IoSerializer& ioSerializer) const
	{
		PsimagLite::String root = label;
		ioSerializer.createGroup(root);
		ioSerializer.write(root + "/filename", filename_);
		ioSerializer.write(root + "/labelForEnergy", labelForEnergy_);
	}

	friend std::ostream& operator<<(std::ostream& os, const RestartStruct& c)
	{
	    if (c.filename_ == "") return os;

	    os<<"RestartStruct.filename="<<c.filename_<<"\n";
	    os<<"RestartStruct.labelForEnergy="<<c.labelForEnergy_<<"\n";
	    return os;
	}

	friend std::istream& operator>>(std::istream& is,RestartStruct& c)
	{
	    is>>c.filename_;
	    is>>c.labelForEnergy_;
	    return is;
	}

private:

	PsimagLite::String filename_;
	PsimagLite::String labelForEnergy_;
};

} // namespace Dmrg

#endif // RestartStruct_H

