# This is repeated here but should possibly just use the custom
# command and info from dmrgpp/src

add_library(provenance OBJECT Provenance.cpp)
target_link_libraries(provenance PRIVATE dmrgpp psimaglite)
set(dmftpp_src cincuenta.cpp)
add_executable(dmftpp ${dmftpp_src} $<TARGET_OBJECTS:provenance>)
target_link_libraries(dmftpp PUBLIC dmrgpp kronutil)


add_executable(freeGimp freeGimp.cpp $<TARGET_OBJECTS:provenance>)
target_link_libraries(freeGimp PUBLIC dmrgpp kronutil)
