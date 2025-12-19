#include "RedirectOutput.hh"

namespace PsimagLite {

std::string RedirectOutput::app_name_;
std::string RedirectOutput::app_name2_;
std::streambuf* RedirectOutput::cout_buffer_ = nullptr;
std::ofstream RedirectOutput::cout_stream_;
bool RedirectOutput::redirected_ = false;

}
