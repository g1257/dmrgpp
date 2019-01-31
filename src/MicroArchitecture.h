#ifndef MICROARCHITECTURE_H
#define MICROARCHITECTURE_H

#include "AllocatorCpu.h"
#include <fstream>

namespace PsimagLite {

class MicroArchitecture {

public:

	MicroArchitecture()
	{
		vendorId_ = readLabel("/proc/cpuinfo", "vendor_id");
		if (vendorId_.find("Intel") != String::npos) {
			vendorId_ = "Intel";
		} else {
			if (vendorId_.find("AMD") != String::npos)
				vendorId_ = "AMD";
			else
				vendorId_ = extractOther(vendorId_);
		}
	}

	const String& vendorId() const { return vendorId_; }

private:

	static String readLabel(String file, String label)
	{
		std::ifstream fin(file.c_str());
		if (!fin || fin.bad() || !fin.good())
			return "";

		String line;
		bool flag = false;
		for (; std::getline(fin, line);) {
			if (line.find(label) != String::npos) {
				flag = true;
				break;
			}
		}

		fin.close();

		if (!flag) return "";
		return line;
	}

	static String extractOther(String vendorId)
	{
		std::string literal = ":";
		size_t it = vendorId.find(literal);
		if (it == String::npos) return "";

		SizeType index = it + 1;
		SizeType ind = index;
		SizeType total = vendorId.length();
		for (; ind < total; ++ind)
			if (vendorId[ind] != ' ' && vendorId[ind] != '\t') break;


		String tmp1 = vendorId.substr(ind, total - ind);
		return tmp1;
	}

	String vendorId_;
};
}
#endif // MICROARCHITECTURE_H
