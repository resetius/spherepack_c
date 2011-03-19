project(spherepack C)

set(SPHEREPACK_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)

set(source_o0 
  vhagc.c
  vtses.c
  vhsgs.c
  vhses.c
  vtsgs.c
  vhags.c
  vhaes.c
  vhaec.c
  )

set(source_o1
  vtsgc.c
  vtsec.c
  vhsgc.c
  vhsec.c
  hrfft.c
  )

set(source
  vlapgs.c
  vlapec.c
  visgau.c
  visequ.c
  trvsph.c
  sphcom.c
  shses.c
  shpg.c
  shpe.c
  shigc.c
  shaes.c
  sfvpgc.c
  lfp.c
  ivrtec.c
  ivlapes.c
  islapgc.c
  isfvpgc.c
  idvtgs.c
  idvtec.c
  idivec.c
  grades.c
  divgs.c
  alf.c
  vsurf.c
  vshifte.c
  vrtgs.c
  vrtgc.c
  vrtes.c
  vrtec.c
  vlapgc.c
  vlapes.c
  visgeo.c
  trssph.c
  sshifte.c
  slapgs.c
  slapgc.c
  slapes.c
  slapec.c
  shsgs.c
  shsgc.c
  shsec.c
  shigs.c
  shags.c
  shagc.c
  shaec.c
  sfvpgs.c
  sfvpes.c
  sfvpec.c
  pow_di.c
  ivrtgs.c
  ivrtgc.c
  ivrtes.c
  ivlapgs.c
  ivlapgc.c
  ivlapec.c
  islapgs.c
  islapes.c
  islapec.c
  isfvpgs.c
  isfvpes.c
  isfvpec.c
  ihgeod.c
  igradgs.c
  igradgc.c
  igrades.c
  igradec.c
  idvtgc.c
  idvtes.c
  idivgs.c
  idivgc.c
  idives.c
  helmsph.c
  gradgs.c
  gradgc.c
  gradec.c
  geo2math.c
  gaqd.c
  divgc.c
  dives.c
  divec.c
  )

add_library(spherepack ${source_o0} ${source_o1} ${source})

if (MSVC)
  if (MSVC_VERSION GREATER 1400)
    # TODO: add COMPILE_FLAGS_<CONFIGURATION> support in cmake
    set_source_files_properties(${source_o1} PROPERTIES
      COMPILE_FLAGS "/O1")
    set_source_files_properties(${source_o0} PROPERTIES
      COMPILE_FLAGS "/Od")
  endif()
endif ()


set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")

#add_executable(helmsph helmsph.c)
#target_link_libraries(helmsph spherepack m f2c)

