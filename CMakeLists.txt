cmake_minimum_required(VERSION 3.14)

project(remez)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

set( CMAKE_CXX_VISIBLITY_PRESET hidden )
set( CMAKE_VISIBILITY_INLINES_HIDDEN 1 )

add_library(
    remez
    remez.c
    remez.h
)

include( GenerateExportHeader )

set( EXPORT_FILE "${CMAKE_CURRENT_SOURCE_DIR}/remez_export.h" )

generate_export_header(
    remez
    BASE_NAME REMEZ
    EXPORT_MACRO_NAME REMEZ_PUBLIC
    NO_EXPORT_MACRO_NAME REMEZ_PRIVATE
    STATIC_DEFINE REMEZ_STATIC_BUILD
    EXPORT_FILE_NAME ${EXPORT_FILE}
)

target_include_directories( remez INTERFACE "$<INSTALL_INTERFACE:include>" )

if( WIN32 )
    install(
        TARGETS remez
        EXPORT remez
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION bin
        ARCHIVE DESTINATION lib
    )

    install( FILES $<TARGET_PDB_FILE:remez> DESTINATION bin OPTIONAL )
else()
    install(
        TARGETS remez
        EXPORT remez
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
    )
endif()

set( PUBLIC_HEADERS remez.h remez_export.h )

install( FILES ${PUBLIC_HEADERS} DESTINATION include/remez )
install( EXPORT remez DESTINATION lib/cmake/remez FILE remez-config.cmake )
