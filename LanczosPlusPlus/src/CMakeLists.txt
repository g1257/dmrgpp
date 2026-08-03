add_library(lanczos_globals LanczosGlobals.cpp)
target_include_directories(lanczos_globals PUBLIC Engine ${CMAKE_SOURCE_DIR})
target_link_libraries(lanczos_globals PUBLIC psimaglite::psimaglite)

add_executable(
  lanczos
  lanczos.cpp
  LanczosDriver0.cpp
  LanczosDriver1.cpp
  LanczosDriver2.cpp
  LanczosDriver3.cpp
  LanczosDriver4.cpp
  LanczosDriver5.cpp)
target_include_directories(lanczos PRIVATE Engine ${CMAKE_SOURCE_DIR})
target_link_libraries(lanczos PUBLIC lanczos_globals psimaglite::psimaglite)

add_executable(thermal thermal.cpp)
target_link_libraries(thermal psimaglite::psimaglite)

add_executable(lorentzian lorentzian.cpp)

target_link_libraries(lorentzian psimaglite::psimaglite)

install(TARGETS lanczos thermal lorentzian RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
