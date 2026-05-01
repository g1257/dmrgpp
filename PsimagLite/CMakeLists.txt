# Directories ordered alphabetically
add_subdirectory(AST)
add_subdirectory(Ainur)
add_subdirectory(Geometry)
add_subdirectory(loki)
add_subdirectory(Utilities)

if(BUILD_TESTING)
  add_subdirectory(examples)
endif()
