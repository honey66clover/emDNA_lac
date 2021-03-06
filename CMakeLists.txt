# emDNA CMake


# --- CMake gobal directives --- #

cmake_minimum_required (VERSION 2.6)


# --- emDNA project --- #

project (emDNA CXX C)
include (emDNACommon.cmake)
set (EMDNA_MAJOR_VERSION "0")
set (EMDNA_MINOR_VERSION "0")
set (EMDNA_BUGFIX_VERSION "0")
set (EMDNA_VERSION_NUMBER "${EMDNA_MAJOR_VERSION}.${EMDNA_MINOR_VERSION}.${EMDNA_BUGFIX_VERSION}")
set (EMDNA_VERSION "${EMDNA_VERSION_NUMBER}")


# --- build type options --- #
if (NOT CMAKE_BUILD_TYPE)

    # compiler
    # force to discard any user input
    set (a ${CMAKE_CXX_COMPILER})
    mark_as_advanced (CLEAR CMAKE_CXX_COMPILER)
    set (CMAKE_CXX_COMPILER ${a}
        CACHE
        STRING
        "C++ compiler"
        FORCE)

    # set build type as MacOS
    # force to discard any user input
    set (CMAKE_BUILD_TYPE MACOS
        CACHE
        STRING
        "MacOS build settings."
        FORCE)

    # set MacOS c++ compiler flags
    set (MACOS_FLAGS "-O2 -march=native")
    set (MACOS_FLAGS "${MACOS_FLAGS} -std=gnu++11 -stdlib=libc++")
    set (MACOS_FLAGS "${MACOS_FLAGS} -m64 -funroll-loops -msse3 -mssse3")
    set (MACOS_FLAGS "${MACOS_FLAGS} -msse4.1 -msse4.2 -ffast-math")
    set (MACOS_FLAGS "${MACOS_FLAGS} -fstrict-aliasing")
    set (MACOS_FLAGS "${MACOS_FLAGS} -Wno-deprecated -Weffc++")
    set (CMAKE_CXX_FLAGS_MACOS ${MACOS_FLAGS}
        CACHE
        FILEPATH
        "C++ MacOS compiler flags.")

    # set Unix cluster c++ compiler flags
    set (UNIX_FLAGS "-O2 -std=gnu++11")
    set (UNIX_FLAGS "${UNIX_FLAGS} -m64 -march=native -funroll-loops")
    set (UNIX_FLAGS "${UNIX_FLAGS} -msse3 -mssse3 -ffast-math")
    set (UNIX_FLAGS "${UNIX_FLAGS} -fstrict-aliasing")
    set (UNIX_FLAGS "${UNIX_FLAGS} -Wno-deprecated -s")
    set (CMAKE_CXX_FLAGS_UNIX ${UNIX_FLAGS}
        CACHE
        FILEPATH
        "C++ Unix compiler flags.")

	# OS X architectures
	if (APPLE)
		set (CMAKE_OSX_ARCHITECTURES "x86_64"
                CACHE
                STRING
                "Architecture to build against (OS X only)."
                FORCE)
	endif (APPLE)
		
endif (NOT CMAKE_BUILD_TYPE)


# --- source tree directives --- #

# add all source directories
set (LibHeadersDirectories
    src
    src/alglib
    src/DNA
    src/minim
    src/serialization
    src/utils
)
set (LibSourcesDirectories
    src/alglib
	src/DNA
	src/minim
    src/serialization
    src/utils
	)
set (AppHeadersDirectories
    src/app
    )
set (AppSourcesDirectories
    src/app
    )

# headers and cpp files list
append_files (emDNALibHeaders "h" ${LibHeadersDirectories})
append_files (emDNALibCppSources "cpp" ${LibSourcesDirectories})
append_files (emDNAAppHeaders "h" ${AppHeadersDirectories})
append_files (emDNAAppCppSources "cpp" ${AppSourcesDirectories})

# include directories
include_directories (BEFORE src)
include_directories (BEFORE src/alglib)
include_directories (BEFORE src/app)
include_directories (BEFORE src/DNA)
include_directories (BEFORE src/minim)
include_directories (BEFORE src/serialization)
include_directories (BEFORE src/utils)

# DNASim
SET(DNASIM_INCLUDE "/home/juanwei/DNAsim_install/include/DNASim")
SET(DNASIM_LIB_STATIC "/home/juanwei/DNAsim_install/lib/libDNASim.so") 
find_library(DNASIM_LIB_STATIC DNASim_static
	PATHS
	/usr/local/lib/DNASim/lib
	)
if (DNASIM_LIB_STATIC-NOTFOUND)
	message(WARNING "could not find DNASim static library")
endif (DNASIM_LIB_STATIC-NOTFOUND)
find_path(DNASIM_INCLUDE DNASim_Includes.h
	PATHS
	/usr/local/lib/DNASim/include/DNASim
	)
if (DNASIM_INCLUDE-NOTFOUND)
	message(WARNING "could not find DNASim library headers")
endif (DNASIM_INCLUDE-NOTFOUND)
include_directories (BEFORE ${DNASIM_INCLUDE})


# --- target directives --- #

# emDNA library
add_library (emDNALib STATIC
    ${emDNALibHeaders}
    ${emDNALibCppSources}
    )

# library name
set_target_properties(emDNALib
    PROPERTIES
    OUTPUT_NAME emDNA
    COMPILER_FLAGS "-fno-visibility-inlines-hidden -fvisibility=default"
    )

# link
target_link_libraries (emDNALib ${DNASIM_LIB_STATIC})

# emDNA app
add_executable(emDNA ${emDNAAppCppSources})
target_link_libraries(emDNA ${DNASIM_LIB_STATIC} emDNALib)


# --- install directives --- #

# library install
set (emDNALibPermissions OWNER_READ;OWNER_WRITE;OWNER_EXECUTE)
set (emDNALibPermissions ${emDNALibPermissions} GROUP_READ;GROUP_EXECUTE)
set (emDNALibPermissions ${emDNALibPermissions} WORLD_READ;WORLD_EXECUTE)
install (TARGETS emDNALib
    DESTINATION lib
    PERMISSIONS ${emDNALibPermissions}
    )

# headers install
set (emDNAHeadersPermissions OWNER_READ;OWNER_WRITE)
set (emDNAHeadersPermissions ${emDNAHeadersPermissions} GROUP_READ)
set (emDNAHeadersPermissions ${emDNAHeadersPermissions} WORLD_READ)
install (FILES ${emDNALibHeaders}
    DESTINATION include/emDNA
    PERMISSIONS ${emDNAHeadersPermissions}
    )

# app and tools
set (emDNAPermissions OWNER_READ;OWNER_WRITE;OWNER_EXECUTE)
set (emDNAPermissions ${emDNAPermissions} GROUP_READ;GROUP_EXECUTE)
set (emDNAPermissions ${emDNAPermissions} WORLD_READ;WORLD_EXECUTE)
install (TARGETS emDNA
	DESTINATION bin
    PERMISSIONS ${emDNAPermissions}
	)


# --- status messages --- #
message (STATUS "")
message (STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message (STATUS
"C++ compiler (CMAKE_CXX_COMPILER): ${CMAKE_CXX_COMPILER}")
message (STATUS
"C++ compiler flags (CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}): "
"${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")
message (STATUS "Installation path (CMAKE_INSTALL_PREFIX): "
"${CMAKE_INSTALL_PREFIX}")
if (APPLE)
message (STATUS "OS X architecture to build (CMAKE_OSX_ARCHITECTURES): "
"${CMAKE_OSX_ARCHITECTURES}")
message (STATUS "OS X system target (CMAKE_OSX_DEPLOYMENT_TARGET): "
"${CMAKE_OSX_DEPLOYMENT_TARGET}")
endif (APPLE)
