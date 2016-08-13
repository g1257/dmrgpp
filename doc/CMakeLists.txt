add_custom_command( OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/manual.tex
  COMMAND bash -c "cp ${CMAKE_CURRENT_SOURCE_DIR}/*.pdf ${CMAKE_CURRENT_BINARY_DIR}/"
  COMMAND bash -c "cp ${CMAKE_CURRENT_SOURCE_DIR}/*.png ${CMAKE_CURRENT_BINARY_DIR}/"
  COMMAND find ${dmrgpp_SOURCE_DIR} -iname "*.h" -or -iname "*.cpp" | perl ${PsimagLite_SOURCE_DIR}/../scripts/doc.pl ${CMAKE_CURRENT_SOURCE_DIR}/manual.ptex ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} 
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating manual.tex from manual.ptex"
#  VERBATIM
  )

add_custom_target( manual.pdf ALL
  COMMAND pdflatex manual.tex
  COMMAND bibtex manual.aux
  COMMAND pdflatex manual.tex
  COMMAND pdflatex manual.tex
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating manual.pdf"
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/manual.tex
  VERBATIM
  )




