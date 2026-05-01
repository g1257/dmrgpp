# Directories ordered alphabetically
# add_subdirectory(AST)
add_subdirectory(Ainur)
#add_subdirectory(Geometry)
add_subdirectory(Io)
add_subdirectory(Loki)
add_subdirectory(Utilities)

if(BUILD_TESTING)
  add_subdirectory(examples)
endif()
