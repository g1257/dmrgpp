
find_package(GSL)
if(NOT GSL_FOUND)
  MESSAGE(WARNING "libGSL required for cincuenta at runtime.")
endif()

add_subdirectory(src)
