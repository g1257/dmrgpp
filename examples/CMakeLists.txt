find_package(BLAS)
find_package(LAPACK)

# testLanczos executable
add_executable(testLanczos
               testLanczos.cpp)

target_link_libraries(
  testLanczos
  psimaglite
  BLAS::BLAS
  LAPACK::LAPACK
)
