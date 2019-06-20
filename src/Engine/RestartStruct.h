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
	    : filename_(""),into_("GroundState"),labelForPsi_("PSI"),labelForEnergy_("Energy")
	{}

	template<typename SomeInputType>
	void read(SomeInputType& io)
	{
		try {
			io.readline(into_, "RestartInto=");
		} catch (std::exception&) {}

		if (into_ != "All" && into_ != "GroundState") {
			PsimagLite::String str = "FATAL: RestartInto=All | GroundState\n";
			throw PsimagLite::RuntimeError(str);
		}

		try {
			io.readline(labelForPsi_, "RestartLabelForPsi=");
		} catch (std::exception&) {}

		try {
			io.readline(labelForEnergy_, "RestartLabelForEnergy=");
		} catch (std::exception&) {}
	}

	void setFilename(PsimagLite::String f) { filename_ = f; }

	PsimagLite::String filename() const { return filename_; }

	PsimagLite::String into() const { return into_; }

	PsimagLite::String labelForPsi() const { return labelForPsi_; }

	PsimagLite::String labelForEnergy() const { return labelForEnergy_; }

	SizeType mappingStages(SizeType ind) const
	{
		return ind;
	}

	SizeType mappingTvs(SizeType ind) const
	{
		return ind;
	}

	void write(PsimagLite::String label,
	           PsimagLite::IoSerializer& ioSerializer) const
	{
		PsimagLite::String root = label;
		ioSerializer.createGroup(root);
		ioSerializer.write(root + "/filename", filename_);
		ioSerializer.write(root + "/into", into_);
		ioSerializer.write(root + "/labelForPsi", labelForPsi_);
		ioSerializer.write(root + "/labelForEnergy", labelForEnergy_);
	}

	friend std::ostream& operator<<(std::ostream& os, const RestartStruct& c)
	{
	    if (c.filename_ == "") return os;

	    os<<"RestartStruct.filename="<<c.filename_<<"\n";
	    os<<"RestartStruct.into="<<c.into_<<"\n";
	    os<<"RestartStruct.labelForPsi="<<c.labelForPsi_<<"\n";
	    os<<"RestartStruct.labelForEnergy="<<c.labelForEnergy_<<"\n";
	    return os;
	}

	friend std::istream& operator>>(std::istream& is,RestartStruct& c)
	{
	    is>>c.filename_;
	    is>>c.into_;
	    is>>c.labelForPsi_;
	    is>>c.labelForEnergy_;
	    return is;
	}

private:

	PsimagLite::String filename_;
	PsimagLite::String into_;
	PsimagLite::String labelForPsi_;
	PsimagLite::String labelForEnergy_;
};

} // namespace Dmrg

#endif // RestartStruct_H

