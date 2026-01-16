# PsimagLite tests
# This will be converted into unit tests using Catch2

set(INPUTS_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# affiniteTest needs -DUSE_PTHREADS
add_executable(affinityTest affinityTest.cpp)
target_link_libraries(affinityTest psimaglite)
# add_test(NAME affinityTest COMMAND affinityTest 4 14)

add_executable(akimaSpline akimaSpline.cpp)
target_link_libraries(akimaSpline psimaglite)
add_test(NAME akimaSpline COMMAND akimaSpline ${INPUTS_DIR}/dataForAkimaSpline.txt 0.5 1.2 10)

add_executable(base64test base64test.cpp)
target_link_libraries(base64test psimaglite)
add_test(NAME base64test COMMAND base64test)

add_executable(calculator calculator.cpp)
target_link_libraries(calculator psimaglite)
add_test(NAME calculator COMMAND calculator "-:4.2:*:%t:3.3")

add_executable(checkRunId checkRunId.cpp)
target_link_libraries(checkRunId psimaglite)
add_test(NAME checkRunId COMMAND checkRunId -g 4)

add_executable(closuresTest closuresTest.cpp)
target_link_libraries(closuresTest psimaglite)
add_test(NAME closuresTest COMMAND closuresTest)

# These continued fraction three tests need data; may test in the future
add_executable(combineContinuedFraction combineContinuedFraction.cpp)
target_link_libraries(combineContinuedFraction psimaglite)
# add_test(NAME combineContinuedFraction COMMAND combineContinuedFraction)

add_executable(continuedFraction continuedFraction.cpp)
target_link_libraries(continuedFraction psimaglite)
# add_test(NAME continuedFraction COMMAND continuedFraction)

add_executable(continuedFractionCollection continuedFractionCollection.cpp)
target_link_libraries(continuedFractionCollection psimaglite)
# add_test(NAME continuedFractionCollection COMMAND continuedFractionCollection)

# The fit test needs data
add_executable(fit fit.cpp)
target_link_libraries(fit psimaglite)
# add_test(NAME fit COMMAND fit)

# Integration need the GSL as a TPL
add_executable(integrator integrator.cpp)
target_link_libraries(integrator psimaglite)
# add_test(NAME integrator COMMAND integrator)

add_executable(internode internode.cpp)
target_link_libraries(internode psimaglite)
add_test(NAME internode COMMAND internode 4)

# Needs -DUSE_PTHREADS
add_executable(isBlasThreaded isBlasThreaded.cpp)
target_link_libraries(isBlasThreaded psimaglite)
# add_test(NAME isBlasThreaded COMMAND isBlasThreaded)

# needs input data
add_executable(kernelPolynomial kernelPolynomial.cpp)
target_link_libraries(kernelPolynomial psimaglite)
# add_test(NAME kernelPolynomial COMMAND kernelPolynomial)

add_executable(linearPrediction linearPrediction.cpp)
target_link_libraries(linearPrediction psimaglite)
add_test(NAME linearPrediction COMMAND linearPrediction -f ${INPUTS_DIR}/series.txt -l Series -p 3 -q 10)

# loadImabalance needs -DUSE_PTHREADS
add_executable(loadImbalance loadImbalance.cpp)
target_link_libraries(loadImbalance psimaglite)
# add_test(NAME loadImbalance COMMAND loadImbalance)

# minimizer needs the GSL
add_executable(minimizer minimizer.cpp)
target_link_libraries(minimizer psimaglite)
# add_test(NAME minimizer COMMAND minimizer)

add_executable(myApp01 myApp01.cpp)
target_link_libraries(myApp01 psimaglite)
add_test(NAME myApp01 COMMAND myApp01)

# nested parallelization needs -DUSE_PTHREADS
add_executable(nested nested.cpp)
target_link_libraries(nested psimaglite)
# add_test(NAME nested COMMAND nested)

add_executable(options options.cpp)
target_link_libraries(options psimaglite)
add_test(NAME options COMMAND options "verbose,fast")

add_executable(randomTest randomTest.cpp)
target_link_libraries(randomTest psimaglite)
add_test(NAME randomTest COMMAND randomTest 1234)

# range makes sense to test when using -DUSE_PTHREADS
add_executable(range range.cpp)
target_link_libraries(range psimaglite)
# add_test(NAME range COMMAND range)

add_executable(sfinaeHasMember sfinaeHasMember.cpp)
target_link_libraries(sfinaeHasMember psimaglite)
add_test(NAME sfinaeHasMember COMMAND sfinaeHasMember)

add_executable(sparseSolverTest sparseSolverTest.cpp)
target_link_libraries(sparseSolverTest psimaglite)
add_test(NAME sparseSolverTest COMMAND sparseSolverTest -n 20)

add_executable(sumDecomposition sumDecomposition.cpp)
target_link_libraries(sumDecomposition psimaglite)
add_test(NAME sumDecomposition COMMAND sumDecomposition 4 3)

add_executable(svd svd.cpp)
target_link_libraries(svd psimaglite)
add_test(NAME svd COMMAND svd)

add_executable(testCRSMatrix testCRSMatrix.cpp)
target_link_libraries(testCRSMatrix psimaglite)
add_test(NAME testCRSMatrix COMMAND testCRSMatrix 10 0.3)

add_executable(testExcitedLanczos testExcitedLanczos.cpp)
target_link_libraries(testExcitedLanczos psimaglite)
add_test(NAME testExcitedLanczos COMMAND testExcitedLanczos 2)

# testGemmR[2] both need -DUSE_PTHREADS to make sense
add_executable(testGemmR testGemmR.cpp)
target_link_libraries(testGemmR psimaglite)
# add_test(NAME testGemmR COMMAND testGemmR)

add_executable(testGemmR2 testGemmR2.cpp)
target_link_libraries(testGemmR2 psimaglite)
# add_test(NAME testGemmR2 COMMAND testGemmR2)

add_executable(testInputNg testInputNg.cpp)
target_link_libraries(testInputNg psimaglite)
add_test(NAME testInputNg COMMAND testInputNg ${INPUTS_DIR}/testInput.ain)

add_executable(testIoNg testIoNg.cpp)
target_link_libraries(testIoNg psimaglite)
add_test(NAME testIoNg COMMAND testIoNg)

add_executable(testIoNgBoolean testIoNgBoolean.cpp)
target_link_libraries(testIoNgBoolean psimaglite)
add_test(NAME testIoNgBoolean0 COMMAND testIoNgBoolean 0)
add_test(NAME testIoNgBoolean1 COMMAND testIoNgBoolean 1)

add_executable(testIsClass testIsClass.cpp)
target_link_libraries(testIsClass psimaglite)
add_test(NAME testIsClass COMMAND testIsClass)

add_executable(testLanczos testLanczos.cpp)
target_link_libraries(testLanczos psimaglite)
add_test(NAME testLanczos COMMAND testLanczos)

add_executable(testLanczosMatrixInFile testLanczosMatrixInFile.cpp)
target_link_libraries(testLanczosMatrixInFile psimaglite)
add_test(NAME testLanczosMatrixInFile COMMAND testLanczosMatrixInFile ${INPUTS_DIR}/matrix.txt)

add_executable(testLapack testLapack.cpp)
target_link_libraries(testLapack psimaglite)
add_test(NAME testLapack COMMAND testLapack)

add_executable(testArnoldiSaI testArnoldiSaI.cpp)
target_link_libraries(testArnoldiSaI psimaglite)
add_test(NAME testArnoldiSaI COMMAND testArnoldiSaI)

# I leave this test on for now but it really needs -DUSE_PTHREADS
# the first number (4) is the matrix rank and the next number (1) is the
# number of pthreads to use
# The same is true for testParallelizer2
add_executable(testParallelSvd testParallelSvd.cpp)
target_link_libraries(testParallelSvd psimaglite)
add_test(NAME testParallelSvd COMMAND testParallelSvd 4 1)

add_executable(testParallelizer2 testParallelizer2.cpp)
target_link_libraries(testParallelizer2 psimaglite)
add_test(NAME testParallelizer2 COMMAND testParallelizer2 10 1)

add_executable(testPredicate testPredicate.cpp)
target_link_libraries(testPredicate psimaglite)
add_test(NAME testPredicateTrue COMMAND testPredicate "+:c:l==8" 3)
add_test(NAME testPredicateFalse COMMAND testPredicate "+:c:l==8" 7)

add_executable(test_redirect_output test_redirect_output.cpp)
target_link_libraries(test_redirect_output psimaglite)
add_test(NAME test_redirect_output COMMAND test_redirect_output)

# See testParallelSvd for use of -DUSE_PTHREADS
add_executable(threads threads.cpp)
target_link_libraries(threads psimaglite)
add_test(NAME threads COMMAND threads 1 4)
