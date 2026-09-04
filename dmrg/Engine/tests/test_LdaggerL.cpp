#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>

#ifndef DMRG_EXECUTABLE_PATH
#error "DMRG_EXECUTABLE_PATH must name the dmrg executable"
#endif

#ifndef DMRG_LIOUVILLIAN_INPUT_PATH
#error "DMRG_LIOUVILLIAN_INPUT_PATH must name a Liouvillian input"
#endif

namespace {

class TemporaryDirectory {
public:

	TemporaryDirectory()
	{
		const auto unique = std::to_string(::getpid()) + "_"
		    + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
		path_ = std::filesystem::temp_directory_path() / ("dmrgpp_ldaggerl_" + unique);
		std::filesystem::create_directories(path_);
	}

	~TemporaryDirectory() { std::filesystem::remove_all(path_); }

	const std::filesystem::path& path() const { return path_; }

private:

	std::filesystem::path path_;
};

struct RunResult {
	int         status;
	std::string output;
};

std::string shellQuote(const std::string& value)
{
	std::string quoted = "'";
	for (const char c : value)
		quoted += (c == '\'') ? "'\\''" : std::string(1, c);
	return quoted + "'";
}

RunResult runDmrg(const std::filesystem::path& input, const std::filesystem::path& workDir)
{
	const std::filesystem::path output  = workDir / "process-output.txt";
	const std::string           command = "cd " + shellQuote(workDir.string()) + " && "
	    + shellQuote(DMRG_EXECUTABLE_PATH) + " -f " + shellQuote(input.string()) + " > "
	    + shellQuote(output.string()) + " 2>&1";
	const int status = std::system(command.c_str());

	std::ifstream      stream(output);
	std::ostringstream buffer;
	buffer << stream.rdbuf();
	return { status, buffer.str() };
}

bool exitedSuccessfully(const int status)
{
	return status != -1 && WIFEXITED(status) && WEXITSTATUS(status) == 0;
}

std::filesystem::path writeHermitianInput(const std::filesystem::path& directory,
                                          const std::string&           solverOptions)
{
	const std::filesystem::path input = directory / "hermitian.ain";
	std::ofstream               stream(input);
	stream << "##Ainur1.0\n"
	       << "TotalNumberOfSites=8;\n"
	       << "NumberOfTerms=1;\n"
	       << "DegreesOfFreedom=1;\n"
	       << "GeometryKind=\"chain\";\n"
	       << "GeometryOptions=\"ConstantValues\";\n"
	       << "dir0:Connectors=[1.0];\n"
	       << "hubbardU=[0.0, ...];\n"
	       << "potentialV=[0.0,...];\n"
	       << "Model=\"HubbardOneBand\";\n"
	       << "SolverOptions=\"" << solverOptions << "\";\n"
	       << "Version=\"test\";\n"
	       << "OutputFile=\"test-hermitian\";\n"
	       << "InfiniteLoopKeptStates=50;\n"
	       << "FiniteLoops=[[@auto, 100, 0]];\n"
	       << "TargetElectronsUp=4;\n"
	       << "TargetElectronsDown=4;\n";
	return input;
}

} // namespace

TEST_CASE("LdaggerL rejects Hermitian models", "[LdaggerL]")
{
	TemporaryDirectory temporaryDirectory;
	const auto         input = writeHermitianInput(temporaryDirectory.path(), "LdAgGeRl");

	const RunResult result = runDmrg(input, temporaryDirectory.path());

	INFO(result.output);
	REQUIRE_FALSE(exitedSuccessfully(result.status));
	REQUIRE(result.output.find("LdaggerL cannot be used with a Hermitian model")
	        != std::string::npos);
}

TEST_CASE("LdaggerL accepts non-Hermitian models", "[LdaggerL]")
{
	TemporaryDirectory temporaryDirectory;

	const RunResult result = runDmrg(DMRG_LIOUVILLIAN_INPUT_PATH, temporaryDirectory.path());

	INFO(result.output);
	REQUIRE(exitedSuccessfully(result.status));
}

TEST_CASE("Hermitian models run without LdaggerL", "[LdaggerL]")
{
	TemporaryDirectory temporaryDirectory;
	const auto         input = writeHermitianInput(temporaryDirectory.path(), "none");

	const RunResult result = runDmrg(input, temporaryDirectory.path());

	INFO(result.output);
	REQUIRE(exitedSuccessfully(result.status));
}
