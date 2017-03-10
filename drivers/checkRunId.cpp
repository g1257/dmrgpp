#include "BitManip.h"
#include <iostream>
#include <cstdlib>
#include "ApplicationInfo.h"

void usage(const char *progName)
{
	std::cerr<<"Usage: "<<progName<<" -i number | -g times\n";
}

bool checkRunId(PsimagLite::String id, PsimagLite::String progName)
{
	PsimagLite::String str(id);
	int l = str.length();
	if (l < 3) {
		std::cerr<<progName<<": too short runId="<<id<<"\n";
		return false;
	}

	PsimagLite::String check("");
	PsimagLite::String number("");
	check += str[l-2];
	check += str[l-1];

	for (int i = 0; i < l - 2; ++i)
		number += str[i];

	unsigned long int x = atol(number.c_str());
	int bits = PsimagLite::BitManip::countKernighan(x);
	std::cout<<number<<" "<<check<<" "<<bits<<"\n";
	//std::cout<<bits<<"\n";
	return (bits == atoi(check.c_str()));
}

int main(int argc, char** argv)
{
	int opt = 0;
	int g = 0;
	PsimagLite::String id;
	while ((opt = getopt(argc, argv,"g:i:")) != -1) {
		switch (opt) {
		case 'g':
			g = atoi(optarg);
			break;
		case 'i':
			id = optarg;
			break;
		default:
			usage(argv[0]);
			return 1;
		}
	}

	if (id == "" && g == 0) {
		usage(argv[0]);
		return 1;
	}

	if (id != "") {
		bool b = checkRunId(id, argv[0]);
		return (b) ? 0 : 3;
	}

	for (int i = 0; i < g; ++i) {
		PsimagLite::ApplicationInfo appInfo("test");
		bool b = checkRunId(appInfo.runId(), argv[0]);
		if (!b) {
			std::cerr<<"Found invalid number\n";
			return 4;
		}
	}

	return 0;
}

