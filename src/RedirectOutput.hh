#ifndef REDIRECTOUTPUT_HH
#define REDIRECTOUTPUT_HH
#include <iostream>
#include <fstream>
#include "PsimagLite.h"

namespace PsimagLite {

class RedirectOutput {

public:

	static void setAppName(const std::string& app_name)
	{
		setAppName(app_name, app_name);
	}

	static void setAppName(const std::string& app_name,
	                       const std::string& app_name2)
	{
		app_name_ = app_name;
		app_name2_ = app_name2;
	}

	//  std::ofstream::app or std::ofstream::out
	static void doIt(const std::string& label,
	                 std::ios_base::openmode open_mode,
	                 bool unbuffered)
	{
		cout_stream_.open(label.c_str(), open_mode);
		if (!cout_stream_ || cout_stream_.bad()
		    || !cout_stream_.good()) {
			PsimagLite::String str(app_name_);
			str += ": Could not redirect std::cout to " + label + "\n";
			err(str);
		}

		std::cerr << app_name2_ << ": Standard output sent to ";
		std::cerr << label << "\n";
		std::cerr.flush();
		cout_buffer_ = std::cout.rdbuf(); // save old buf
		std::cout.rdbuf(cout_stream_.rdbuf()); // redirect std::cout to file
		if (unbuffered) {
			std::cout.setf(std::ios::unitbuf);
			cout_stream_.setf(std::ios::unitbuf);
		}

		atexit(restoreCoutBuffer);
		redirected_ = true;
	}

	bool hasBeenRedirected() const { return redirected_; }

private:

	static void restoreCoutBuffer()
	{
		if (cout_buffer_ == 0) {
			return;
		}

		cout_stream_.close();
		std::cout.rdbuf(cout_buffer_);
	}

	static std::string app_name_;
	static std::string app_name2_;
	static std::streambuf* cout_buffer_;
	static std::ofstream cout_stream_;
	static bool redirected_;
};

std::string RedirectOutput::app_name_;
std::string RedirectOutput::app_name2_;
std::streambuf* RedirectOutput::cout_buffer_ = nullptr;
std::ofstream RedirectOutput::cout_stream_;
bool RedirectOutput::redirected_ = false;
}
#endif // REDIRECTOUTPUT_HH
