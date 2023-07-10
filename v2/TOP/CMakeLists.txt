cmake_minimum_required(VERSION 3.11.0 FATAL_ERROR)
###########################Package@#####################################
project(PDTOP LANGUAGES C  CXX)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH} )
find_package(LHAPDF REQUIRED)
message(STATUS "LHAPDF_VERSION=${LHAPDF_VERSION} LHAPDF_LIBRARIES=${LHAPDF_LIBRARIES} LHAPDF_INCLUDE_DIRS=${LHAPDF_INCLUDE_DIRS} LHAPDF_DATA_PATH=${LHAPDF_DATA_PATH}")
find_package(GSL REQUIRED)

add_executable(top++  
${CMAKE_CURRENT_SOURCE_DIR}/top++2.0/FixedOrder.cpp
${CMAKE_CURRENT_SOURCE_DIR}/top++2.0/FixedOrder.h
${CMAKE_CURRENT_SOURCE_DIR}/top++2.0/lgamma.cpp
${CMAKE_CURRENT_SOURCE_DIR}/top++2.0/lgamma.h
${CMAKE_CURRENT_SOURCE_DIR}/top++2.0/PartonicFlux.cpp
${CMAKE_CURRENT_SOURCE_DIR}/top++2.0/PartonicFlux.h
${CMAKE_CURRENT_SOURCE_DIR}/top++2.0/psin.cpp
${CMAKE_CURRENT_SOURCE_DIR}/top++2.0/psin.h
${CMAKE_CURRENT_SOURCE_DIR}/top++2.0/Resummation.cpp
${CMAKE_CURRENT_SOURCE_DIR}/top++2.0/Resummation.h
${CMAKE_CURRENT_SOURCE_DIR}/top++2.0/SubtrFlux.cpp
${CMAKE_CURRENT_SOURCE_DIR}/top++2.0/SubtrFlux.h
${CMAKE_CURRENT_SOURCE_DIR}/top++2.0/top++.cpp
${CMAKE_CURRENT_SOURCE_DIR}/top++2.0/Utilities.cpp
${CMAKE_CURRENT_SOURCE_DIR}/top++2.0/Utilities.h)
target_link_libraries(top++ PRIVATE LHAPDF::LHAPDF GSL::gslcblas GSL::gsl)

########################################################################
########################################################################
SET(temp "83300 CABCHSV2023nnlo 1\n")
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/LOCALLHAPDF/pdfsets.index "${temp}") 
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../CABCHSV2023nnlo DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/LOCALLHAPDF/)



ENABLE_TESTING()
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/CABCHSV2023nnlo/pdf.cfg "MSTW2008nlo68cl                    Asymmetric\n") 
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/share/top++best-TEV-pdfvar-MSTW2008nlo68cl.cfg  ${CMAKE_CURRENT_BINARY_DIR}/MSTW2008nlo68cl/top++.cfg COPYONLY)
add_test(NAME top++MSTW2008nlo68cl  COMMAND top++ WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/MSTW2008nlo68cl/)
SET_TESTS_PROPERTIES( top++MSTW2008nlo68cl  PROPERTIES ENVIRONMENT "LHAPDF_DATA_PATH=${LHAPDF_DATA_PATH};")


file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/CABCHSV2023nnlo/pdf.cfg "CABCHSV2023nnlo                    Symmetric\n") 
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/share/top++best-TEV-pdfvar-CABCHSV2023nnlo.cfg  ${CMAKE_CURRENT_BINARY_DIR}/CABCHSV2023nnlo/top++.cfg COPYONLY)
add_test(NAME top++CABCHSV2023nnlo  COMMAND top++  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CABCHSV2023nnlo/)
SET_TESTS_PROPERTIES( top++CABCHSV2023nnlo  PROPERTIES ENVIRONMENT "LHAPATH=${CMAKE_CURRENT_BINARY_DIR}/LOCALLHAPDF/;")



