# Add -Wextra in the future
add_library(psimaglite_io TypeToH5.cpp)

target_include_directories(psimaglite_io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_features(psimaglite_io PUBLIC cxx_std_17)
target_link_libraries(psimaglite_io PUBLIC psimaglite::psimaglite)

add_library(psimaglite::io ALIAS psimaglite_io)

# find_package(MPI COMPONENTS CXX)
if(MPI_FOUND)
  target_compile_definitions(psimaglite_io PUBLIC USE_MPI)
  target_link_libraries(psimaglite_io PUBLIC MPI::MPI_CXX)

  # Prefer HDF5 parallel when MPI is enabled
  set(HDF5_PREFER_PARALLEL ON)
endif()

find_package(
  HDF5 REQUIRED
  COMPONENTS CXX HL
  MODULE)
if(MPI_FOUND AND NOT HDF5_IS_PARALLEL)
  message(WARNING "HDF5 that was found \(${HDF5_DIR}\) is not parallel")
endif()
target_link_libraries(psimaglite_io PUBLIC HDF5::HDF5)

find_package(Threads)
if(Threads_FOUND)
  target_compile_definitions(psimaglite_io PUBLIC USE_PTHREADS)
  target_link_libraries(psimaglite_io PUBLIC Threads::Threads)
endif()

if(BUILD_SHARED_LIBS)
  install(TARGETS psimaglite_io LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

# tests below
add_executable(testIoNg testIoNg.cpp)
target_link_libraries(testIoNg psimaglite::io psimaglite::psimaglite)
add_test(NAME testIoNg COMMAND testIoNg)

add_executable(testIoNgBoolean testIoNgBoolean.cpp)
target_link_libraries(testIoNgBoolean psimaglite::io psimaglite::psimaglite)
add_test(NAME testIoNgBoolean0 COMMAND testIoNgBoolean 0)
add_test(NAME testIoNgBoolean1 COMMAND testIoNgBoolean 1)
