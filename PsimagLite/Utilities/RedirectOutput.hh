#ifndef REDIRECTOUTPUT_HH
#define REDIRECTOUTPUT_HH
#include "PsimagLite.h"
#include <fstream>
#include <iostream>

namespace PsimagLite {

class RedirectOutput {

public:

	RedirectOutput()
	    : cout_buffer_(nullptr)
	    , redirected_(false)
	{ }

	void doIt(const std::string& label, std::ios_base::openmode open_mode, bool unbuffered)
	{
		if (cout_buffer_ || redirected_) {
			err("RedirectOutput::doIt(): ERROR: Already redirected\n");
		}

		cout_stream_.open(label.c_str(), open_mode);
		if (!cout_stream_ || cout_stream_.bad() || !cout_stream_.good()) {
			err("Could not redirect std::cout to " + label + "\n");
		}

		std::cerr << "Standard output sent to " << label << "\n";
		std::cerr.flush();

		assert(!cout_buffer_);
		cout_buffer_ = std::cout.rdbuf(); // save old buf
		std::cout.rdbuf(cout_stream_.rdbuf()); // redirect std::cout to file

		if (unbuffered) {
			std::cout.setf(std::ios::unitbuf);
			cout_stream_.setf(std::ios::unitbuf);
		}

		redirected_ = true;
	}

	bool hasBeenRedirected() const { return redirected_; }

	~RedirectOutput()
	{
		if (!cout_buffer_) {
			return;
		}

		cout_stream_.close();
		std::cout.rdbuf(cout_buffer_);
		cout_buffer_ = nullptr;
	}

	RedirectOutput(const RedirectOutput&) = delete;

	RedirectOutput& operator=(const RedirectOutput&) = delete;

	std::streambuf* cout_buffer_;
	std::ofstream   cout_stream_;
	bool            redirected_;
};
}
#endif // REDIRECTOUTPUT_HH
