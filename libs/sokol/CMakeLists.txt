project(sokol)

SET(SRCS
	sokol_glue.h
	sokol_app.h
	sokol_gfx.h
	sokol.c
)

add_library(sokol ${SRCS})
target_include_directories(sokol PUBLIC ${PROJECT_SOURCE_DIR})
source_group("" FILES ${SRCS})