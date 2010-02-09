project(spherepack C)

add_library(spherepack
	alf.c                divgc.c     hrfft.c   idvtgs.c   isfvpgc.c  ivlapgc.c  lfp.c     shaes.c  shses.c   sshifte.c  vhses.c   vlapgc.c   vtsec.c
	       divgs.c     idivec.c  igradec.c  isfvpgs.c  ivlapgs.c  lfpt.c    shagc.c  shsgc.c   trssph.c   vhsgc.c   vlapgs.c   vtses.c
	          gaqd.c      idives.c  igrades.c  islapec.c  ivrtec.c     shags.c  shsgs.c   trvsph.c   vhsgs.c   vrtec.c    vtsgc.c
	  geo2math.c  idivgc.c  igradgc.c  islapes.c  ivrtes.c   sfvpec.c  shigc.c  slapec.c  vhaec.c    visequ.c  vrtes.c    vtsgs.c
	      gradec.c    idivgs.c  igradgs.c  islapgc.c  ivrtgc.c   sfvpes.c  shigs.c  slapes.c  vhaes.c    visgau.c  vrtgc.c
	      grades.c    idvtec.c  ihgeod.c   islapgs.c  ivrtgs.c   sfvpgc.c  shpe.c   slapgc.c  vhagc.c    visgeo.c  vrtgs.c
	divec.c              gradgc.c    idvtes.c  isfvpec.c  ivlapec.c  lfim.c     sfvpgs.c  shpg.c   slapgs.c  vhags.c    vlapec.c  vshifte.c
	dives.c              gradgs.c    idvtgc.c  isfvpes.c  ivlapes.c  lfin.c     shaec.c   shsec.c  sphcom.c  vhsec.c    vlapes.c  vsurf.c pow_di.c
       )


add_executable(helmsph helmsph.c)
target_link_libraries(helmsph spherepack m f2c)

