# CMakeLists.txt
#
# Add executables
add_executable(affinityTest affinityTest.cpp)
target_link_libraries(affinityTest psimaglite)

add_executable(akimaSpline akimaSpline.cpp)
target_link_libraries(akimaSpline psimaglite)

add_executable(base64test base64test.cpp)
target_link_libraries(base64test psimaglite)

add_executable(calculator calculator.cpp)
target_link_libraries(calculator psimaglite)

add_executable(checkRunId checkRunId.cpp)
target_link_libraries(checkRunId psimaglite)

add_executable(closuresTest closuresTest.cpp)
target_link_libraries(closuresTest psimaglite)

add_executable(combineContinuedFraction combineContinuedFraction.cpp)
target_link_libraries(combineContinuedFraction psimaglite)

add_executable(continuedFraction continuedFraction.cpp)
target_link_libraries(continuedFraction psimaglite)

add_executable(continuedFractionCollection continuedFractionCollection.cpp)
target_link_libraries(continuedFractionCollection psimaglite)

add_executable(fit fit.cpp)
target_link_libraries(fit psimaglite)

add_executable(integrator integrator.cpp)
target_link_libraries(integrator psimaglite)

add_executable(internode internode.cpp)
target_link_libraries(internode psimaglite)

add_executable(isBlasThreaded isBlasThreaded.cpp)
target_link_libraries(isBlasThreaded psimaglite)

add_executable(kernelPolynomial kernelPolynomial.cpp)
target_link_libraries(kernelPolynomial psimaglite)

add_executable(linearPrediction linearPrediction.cpp)
target_link_libraries(linearPrediction psimaglite)

add_executable(loadImbalance loadImbalance.cpp)
target_link_libraries(loadImbalance psimaglite)

add_executable(minimizer minimizer.cpp)
target_link_libraries(minimizer psimaglite)

add_executable(myApp01 myApp01.cpp)
target_link_libraries(myApp01 psimaglite)

add_executable(nested nested.cpp)
target_link_libraries(nested psimaglite)

add_executable(options options.cpp)
target_link_libraries(options psimaglite)

add_executable(randomTest randomTest.cpp)
target_link_libraries(randomTest psimaglite)

add_executable(range range.cpp)
target_link_libraries(range psimaglite)

add_executable(rungeKuttaTest rungeKuttaTest.cpp)
target_link_libraries(rungeKuttaTest psimaglite)

add_executable(sfinaeHasMember sfinaeHasMember.cpp)
target_link_libraries(sfinaeHasMember psimaglite)

add_executable(sparseSolverTest sparseSolverTest.cpp)
target_link_libraries(sparseSolverTest psimaglite)

add_executable(sumDecomposition sumDecomposition.cpp)
target_link_libraries(sumDecomposition psimaglite)

add_executable(svd svd.cpp)
target_link_libraries(svd psimaglite)

add_executable(testCRSMatrix testCRSMatrix.cpp)
target_link_libraries(testCRSMatrix psimaglite)

add_executable(testExcitedLanczos testExcitedLanczos.cpp)
target_link_libraries(testExcitedLanczos psimaglite)

add_executable(testGemmR testGemmR.cpp)
target_link_libraries(testGemmR psimaglite)

add_executable(testGemmR2 testGemmR2.cpp)
target_link_libraries(testGemmR2 psimaglite)

add_executable(testInputNg testInputNg.cpp)
target_link_libraries(testInputNg psimaglite)

add_executable(testIoNg testIoNg.cpp)
target_link_libraries(testIoNg psimaglite)

add_executable(testIoNgBoolean testIoNgBoolean.cpp)
target_link_libraries(testIoNgBoolean psimaglite)

add_executable(testIsClass testIsClass.cpp)
target_link_libraries(testIsClass psimaglite)

add_executable(testLanczos testLanczos.cpp)
target_link_libraries(testLanczos psimaglite)

add_executable(testLanczosMatrixInFile testLanczosMatrixInFile.cpp)
target_link_libraries(testLanczosMatrixInFile psimaglite)

add_executable(testLapack testLapack.cpp)
target_link_libraries(testLapack psimaglite)

add_executable(testParallelSvd testParallelSvd.cpp)
target_link_libraries(testParallelSvd psimaglite)

add_executable(testParallelizer2 testParallelizer2.cpp)
target_link_libraries(testParallelizer2 psimaglite)

add_executable(testPredicate testPredicate.cpp)
target_link_libraries(testPredicate psimaglite)

add_executable(test_redirect_output test_redirect_output.cpp)
target_link_libraries(test_redirect_output psimaglite)

add_executable(threads threads.cpp)
target_link_libraries(threads psimaglite)

# End of executables

# PsimagLite tests
# This will be converted into unit tests using Catch2
add_test(NAME affinityTest COMMAND affinityTest)
add_test(NAME akimaSpline COMMAND akimaSpline)
add_test(NAME base64test COMMAND base64test)
add_test(NAME binaryIoTest COMMAND binaryIoTest)
add_test(NAME binaryRead COMMAND binaryRead)
add_test(NAME calculator COMMAND calculator)
add_test(NAME checkRunId COMMAND checkRunId)
add_test(NAME closuresTest COMMAND closuresTest)
add_test(NAME combineContinuedFraction COMMAND combineContinuedFraction)
add_test(NAME concurrencyTest COMMAND concurrencyTest)
add_test(NAME continuedFraction COMMAND continuedFraction)
add_test(NAME continuedFractionCollection COMMAND continuedFractionCollection)
add_test(NAME fit COMMAND fit)
add_test(NAME integrator COMMAND integrator)
add_test(NAME internode COMMAND internode)
add_test(NAME isBlasThreaded COMMAND isBlasThreaded)
add_test(NAME kernelPolynomial COMMAND kernelPolynomial)
add_test(NAME linearPrediction COMMAND linearPrediction)
add_test(NAME loadImbalance COMMAND loadImbalance)
add_test(NAME minimizer COMMAND minimizer)
add_test(NAME myApp01 COMMAND myApp01)
add_test(NAME nested COMMAND nested)
add_test(NAME openblasConfig COMMAND openblasConfig)
add_test(NAME options COMMAND options)
add_test(NAME randomTest COMMAND randomTest)
add_test(NAME range COMMAND range)
add_test(NAME rungeKuttaTest COMMAND rungeKuttaTest)
add_test(NAME sfinaeHasMember COMMAND sfinaeHasMember)
add_test(NAME sparseSolverTest COMMAND sparseSolverTest)
add_test(NAME sumDecomposition COMMAND sumDecomposition)
add_test(NAME svd COMMAND svd)
add_test(NAME testCRSMatrix COMMAND testCRSMatrix)
add_test(NAME testExcitedLanczos COMMAND testExcitedLanczos)
add_test(NAME testGemmR COMMAND testGemmR)
add_test(NAME testGemmR2 COMMAND testGemmR2)
add_test(NAME testInputNg COMMAND testInputNg)
add_test(NAME testIoNg COMMAND testIoNg)
add_test(NAME testIoNgBoolean COMMAND testIoNgBoolean)
add_test(NAME testIsClass COMMAND testIsClass)
add_test(NAME testLanczos COMMAND testLanczos)
add_test(NAME testLanczosMatrixInFile COMMAND testLanczosMatrixInFile)
add_test(NAME testLapack COMMAND testLapack)
add_test(NAME testParallelSvd COMMAND testParallelSvd)
add_test(NAME testParallelizer2 COMMAND testParallelizer2)
add_test(NAME testPredicate COMMAND testPredicate)
add_test(NAME test_redirect_output COMMAND test_redirect_output)
add_test(NAME threads COMMAND threads)
