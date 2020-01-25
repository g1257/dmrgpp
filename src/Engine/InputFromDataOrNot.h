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

	InputFromDataOrNot(PsimagLite::String filename,
	                   const InputCheckType& inputCheck,
	                   bool filenameIsCout)
	    : ioWriteable_(0), isData_(false)
	{

		if (filenameIsCout)
			readFromCout(filename);
		else
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

	void readFromCout(PsimagLite::String filename)
	{
		std::ifstream fin(filename.c_str());
		if (!fin || !fin.good() || fin.bad()) {
			PsimagLite::String s(__FILE__);
			err(s + " Cannot open file " + filename + "\n");
		}

		static const PsimagLite::String search = "PsiApp::echoBase64: ";
		static const SizeType lsearch = search.length();
		PsimagLite::String str;
		bool found = false;
		while (std::getline(fin, str)) {
			if (str.substr(0, lsearch) != search) continue;

			if (std::getline(fin, str)) found = true;
			break;
		}

		if (!found)
			err("Could not find " + search + " in " + filename + "\n");

		PsimagLite::PsiBase64::Decode base64decode(str);
		data_ = base64decode();
		isData_ = true;
	}

	typename InputNgType::Writeable* ioWriteable_;
	bool isData_;
	PsimagLite::String data_;
};
}
#endif // INPUTFROMDATAORNOT_H
