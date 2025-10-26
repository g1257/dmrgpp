add_library(loki INTERFACE)
# This is the point of this, now we don't have to run around managing the
# include directories the loki library just becomes a target that is "linked"
# and actual dependency to it can be reasoned about at the cmake level.
target_include_directories(
  loki INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
                 $<INSTALL_INTERFACE:include>)

install(
  TARGETS loki
  EXPORT PsimagLiteTargets
  INCLUDES
  DESTINATION include)

export(
  TARGETS loki
  FILE "${CMAKE_CURRENT_BINARY_DIR}/PsimagLiteTargets.cmake"
  NAMESPACE psimaglite::
  APPEND # Important: append to existing export set
)
