cmake_minimum_required(VERSION 2.8)
project(spherepack_src C)

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
# vis function deps: set_, line_
#  visgau.c
#  visgeo.c
#  vsurf.c
#  visequ.c
  trvsph.c
  sphcom.c
  shses.c
# rand_ function dep
#  shpg.c
#  shpe.c
  shigc.c
  shaes.c
  sfvpgc.c
# fftpack dep
#  lfp.c
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
  vshifte.c
  vrtgs.c
  vrtgc.c
  vrtes.c
  vrtec.c
  vlapgc.c
  vlapes.c
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
  islapec3l.c
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

add_library(spherepack SHARED ${source_o0} ${source_o1} ${source} spherepack.def)

if (MSVC)
  if (MSVC_VERSION GREATER 1400)
    # TODO: add COMPILE_FLAGS_<CONFIGURATION> support in cmake
    set_source_files_properties(${source_o1} PROPERTIES
      COMPILE_FLAGS "/Od")
    set_source_files_properties(${source_o0} PROPERTIES
      COMPILE_FLAGS "/Od")
  endif()
endif ()


set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

if (WIN32)
	target_link_libraries(spherepack f2cembedded)
else ()
	target_link_libraries(spherepack m f2cembedded)
endif ()

#add_executable(helmsph helmsph.c)
#target_link_libraries(helmsph spherepack m f2c)

