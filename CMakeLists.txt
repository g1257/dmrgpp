cmake_minimum_required(VERSION 3.31)

project(
  dmrgpp
  VERSION 1.1
  LANGUAGES CXX)

option(BUILD_TESTING "Enable building tests" ON)
if(BUILD_TESTING)
  enable_testing()

  include(FetchContent)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.11.0
    SYSTEM FIND_PACKAGE_ARGS NAMES Catch2)
  FetchContent_MakeAvailable(Catch2)
  include(Catch)
endif()

find_program(AWK NAMES awk mawk gawk)

if(AWK MATCHES ".+-NOTFOUND")
  message(
    WARNING
      "awk (and mawk and gawk) could not be found (${AWK}). As a result numeric checks of application level tests are disabled"
  )
elseif()
  message(STATUS "awk found numerical checks for application level tests will be enabled.")
endif()

include(GNUInstallDirs)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if(APPLE)
  set(CMAKE_INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
else()
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()

add_subdirectory(PsimagLite)
add_subdirectory(src)

add_subdirectory(cincuenta)
# add_subdirectory(doc)
