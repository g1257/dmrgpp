#ifndef INPUTFROMDATAORNOT_H
#define INPUTFROMDATAORNOT_H
#include "Vector.h"
#include "InputNg.h"
#include "Io/IoNg.h"

namespace Dmrg {

template<typename InputCheckType>
class InputFromDataOrNot {

public:

	typedef PsimagLite::InputNg<InputCheckType> InputNgType;
	typedef PsimagLite::IoNg::In IoNgInType;

	InputFromDataOrNot(PsimagLite::String filename, const InputCheckType& inputCheck)
	    : ioWriteable_(0), isData_(false)
	{

		internal(filename);

		ioWriteable_ = (isData_) ? new typename InputNgType::Writeable(inputCheck, data_)
		                         : new typename InputNgType::Writeable(filename, inputCheck);
		//data_ = "";
	}

	~InputFromDataOrNot()
	{
		delete ioWriteable_;
		ioWriteable_ = 0;
	}

	const typename InputNgType::Writeable& ioWriteable() const
	{
		if (ioWriteable_) return *ioWriteable_;
		throw PsimagLite::RuntimeError("InputFromDataOrNot: INTERNAL ERROR (FATAL)\n");
	}

private:

	void internal(PsimagLite::String filename)
	{
		IoNgInType* io = 0;
		try {
			io = new IoNgInType(filename);
			isData_ = true;
		} catch (...) {
			return;
		}

		PsimagLite::String buffer;
		io->read(buffer, "InputBase64Encoded");
		delete io;
		io = 0;

		PsimagLite::PsiBase64::Decode base64decode(buffer);
		data_ = base64decode();
	}

	typename InputNgType::Writeable* ioWriteable_;
	bool isData_;
	PsimagLite::String data_;
};
}
#endif // INPUTFROMDATAORNOT_H
