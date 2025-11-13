# testLanczos executable
add_executable(testLanczos
               testLanczos.cpp)

target_link_libraries(
  testLanczos
  psimaglite
)

add_test(NAME testLanczos COMMAND testLanczos)
