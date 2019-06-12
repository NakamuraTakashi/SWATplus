cmake_minimum_required(VERSION 2.8)
file(GLOB SWAT_SRCS source_codes/*.f90)
project(swatplus Fortran)
add_executable(swatplus ${SWAT_SRCS})

set_source_files_properties(
  ${SWAT_SRCS}
  PROPERTIES COMPILE_FLAGS "-O3 -ffree-line-length-none")
