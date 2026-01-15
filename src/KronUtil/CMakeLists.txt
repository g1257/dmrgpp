add_library(kronutil KronUtil.cpp util.cpp utilComplex.cpp csc_nnz.cpp)
target_link_libraries(kronutil PRIVATE psimaglite::psimaglite)
