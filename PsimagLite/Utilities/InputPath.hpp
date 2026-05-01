#ifndef INPUT_PATH_HPP
#define INPUT_PATH_HPP

#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

namespace PsimagLite {

class InputPath {
public:

	InputPath()
	    : paths_ { std::filesystem::current_path() }
	{ }

	void push(std::filesystem::path path) { paths_.push_back(std::move(path)); }

	std::string findFirst(const std::string& file) const
	{
		auto it = std::find_if(paths_.rbegin(),
		                       paths_.rend(),
		                       [&file](std::filesystem::path const& dir)
		                       { return std::filesystem::exists(dir / file); });

		if (it == paths_.rend())
			throw std::runtime_error("File " + file
			                         + " could not be found within available paths\n");

		return *it / file;
	}

private:

	std::vector<std::filesystem::path> paths_;
};

}

#endif // INPUT_PATH_HPP
