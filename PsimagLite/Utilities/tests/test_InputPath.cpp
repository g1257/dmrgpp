#include <catch2/catch_test_macros.hpp>

#include <InputPath.hpp>

#include <filesystem>
#include <fstream>

TEST_CASE("InputPath operations", "[InputPath]")
{
	PsimagLite::InputPath inputPath;
	std::string           testFile = "test_input_path.txt";
	std::filesystem::path tempDir  = std::filesystem::temp_directory_path() / "dmrgpp_test_dir";
	std::filesystem::create_directories(tempDir);

	std::filesystem::path filePathInTemp = tempDir / testFile;

	SECTION("Initial path is empty and finds current directory file")
	{
		std::ofstream ofs(testFile);
		ofs << "test";
		ofs.close();

		REQUIRE(std::filesystem::exists(testFile));
		auto found = inputPath.findFirst(testFile);
		CHECK(found == std::filesystem::current_path() / testFile);

		std::filesystem::remove(testFile);
	}

	SECTION("Push new path and find file there")
	{
		std::ofstream ofs(filePathInTemp);
		ofs << "test";
		ofs.close();

		inputPath.push(tempDir);
		auto found = inputPath.findFirst(testFile);
		CHECK(found == filePathInTemp.string());

		std::filesystem::remove(filePathInTemp);
	}

	SECTION("Throws if file not found")
	{
		CHECK_THROWS_AS(inputPath.findFirst("non_existent_file.txt"), std::runtime_error);
	}

	SECTION("LIFO behavior (last pushed path is searched first)")
	{
		std::filesystem::path tempDir2
		    = std::filesystem::temp_directory_path() / "dmrgpp_test_dir2";
		std::filesystem::create_directories(tempDir2);
		std::filesystem::path filePathInTemp2 = tempDir2 / testFile;

		std::ofstream ofs1(filePathInTemp);
		ofs1 << "test1";
		ofs1.close();

		std::ofstream ofs2(filePathInTemp2);
		ofs2 << "test2";
		ofs2.close();

		inputPath.push(tempDir);
		inputPath.push(tempDir2);

		auto found = inputPath.findFirst(testFile);
		CHECK(found == filePathInTemp2.string());

		std::filesystem::remove(filePathInTemp);
		std::filesystem::remove(filePathInTemp2);
		std::filesystem::remove(tempDir2);
	}

	std::filesystem::remove_all(tempDir);
}
