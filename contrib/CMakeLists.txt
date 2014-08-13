project(spherepack_contrib C)

file(GLOB SRC *.c *.h)
add_definitions(-DSYSDEP_H_INCLUDED 
	-DNO_TRUNCATE -DNON_UNIX_STDIO 
	-DNON_POSIX_STDIO -DUSE_CLOCK)

add_library(f2cembedded ${SRC})

