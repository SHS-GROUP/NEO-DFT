add_subdirectory(dft)

set(CCHEM_SOURCES
  ${CCHEM_DFT_SOURCES}
  ${CMAKE_CURRENT_SOURCE_DIR}/runtime.cpp
  #"${CMAKE_CURRENT_SOURCE_DIR}/rutime.cpp"
  PARENT_SCOPE)


