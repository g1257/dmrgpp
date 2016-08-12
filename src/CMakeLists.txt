add_library(PsimagLite STATIC MersenneTwister.cpp Matrix.cpp Mpi.cpp ApplicationInfo.cpp Concurrency.cpp ProgressIndicator.cpp Tokenizer.cpp MemResolv.cpp PsimagLite.cpp PsiBase64.cpp SpecialFunctions.cpp)

target_link_libraries( PsimagLite ${CMAKE_THREAD_LIBS_INIT} )
