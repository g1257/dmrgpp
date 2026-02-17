add_library(kronutil KronUtil.cpp util.cpp utilComplex.cpp csc_nnz.cpp)
target_link_libraries(kronutil PUBLIC psimaglite::psimaglite)

if(BUILD_TESTING)
  add_executable(KronUtil_Test1 test1.cpp)
  target_link_libraries(KronUtil_Test1 kronutil)
  add_test(NAME KronUtil_Test1 COMMAND KronUtil_Test1)

  add_executable(KronUtil_Test2 test2.cpp)
  target_link_libraries(KronUtil_Test2 kronutil)
  add_test(NAME KronUtil_Test2 COMMAND KronUtil_Test2)
endif()
