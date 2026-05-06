# Add -Wextra in the future
add_library(ainur AinurSpirit.cpp AinurConvert.cpp)
target_include_directories(ainur PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_features(ainur PUBLIC cxx_std_17)
target_link_libraries(ainur PUBLIC loki)

add_library(psimaglite::ainur ALIAS ainur)

find_package(Boost)
if(Boost_FOUND)
  target_compile_definitions(ainur PUBLIC USE_BOOST)
  target_link_libraries(ainur PUBLIC Boost::boost)
endif()

if(BUILD_SHARED_LIBS)
  install(TARGETS ainur LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

if(BUILD_TESTING)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/function.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
  add_executable(Ainur_Test1 test.cpp)
  target_link_libraries(Ainur_Test1 psimaglite::ainur psimaglite::psimaglite)
  add_test(NAME Ainur_Test1 COMMAND Ainur_Test1 ${CMAKE_CURRENT_SOURCE_DIR}/prelude.ain
                                    ${CMAKE_CURRENT_SOURCE_DIR}/input0.ain)

  add_executable(Ainur_Test2 testReal.cpp)
  target_link_libraries(Ainur_Test2 psimaglite::ainur psimaglite::psimaglite)
  add_test(NAME Ainur_Test2 COMMAND Ainur_Test2 ${CMAKE_CURRENT_SOURCE_DIR}/input1.ain)
endif()
