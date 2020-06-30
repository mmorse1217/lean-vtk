# - Try to find the CMakeDemo library
# Once done this will define
#
#  CMakeDemo_FOUND - system has CMakeDemo
#  CMakeDemo_INCLUDE_DIR - CMakeDemo include directory
#  CMakeDemo_LIB - CMakeDemo library directory
#  CMakeDemo_LIBRARIES - CMakeDemo libraries to link

if(CMakeDemo_FOUND)
    return()
endif()

# We prioritize libraries installed in /usr/local with the prefix .../CMakeDemo-*, 
# so we make a list of them here
file(GLOB lib_glob "/usr/local/lib/CMakeDemo-*")
file(GLOB inc_glob "/usr/local/include/CMakeDemo-*")

# Find the library with the name "CMakeDemo" on the system. Store the final path
# in the variable CMakeDemo_LIB
find_library(CMakeDemo_LIB 
    # The library is named "CMakeDemo", but can have various library forms, like
    # libCMakeDemo.a, libCMakeDemo.so, libCMakeDemo.so.1.x, etc. This should
    # search for any of these.
    NAMES CMakeDemo
    # Provide a list of places to look based on prior knowledge about the system.
    # We want the user to override /usr/local with environment variables, so
    # this is included here.
    HINTS
        ${CMakeDemo_DIR}
        ${CMAKEDEMO_DIR}
        $ENV{CMakeDemo_DIR}
        $ENV{CMAKEDEMO_DIR}
        ENV CMAKEDEMO_DIR
    # Provide a list of places to look as defaults. /usr/local shows up because
    # that's the default install location for most libs. The globbed paths also
    # are placed here as well.
    PATHS
        /usr
        /usr/local
        /usr/local/lib
        ${lib_glob}
    # Constrain the end of the full path to the detected library, not including
    # the name of library itself.
    PATH_SUFFIXES 
        lib
)

# Find the path to the file "source_file.hpp" on the system. Store the final
# path in the variables CMakeDemo_INCLUDE_DIR. The HINTS, PATHS, and
# PATH_SUFFIXES, arguments have the same meaning as in find_library().
find_path(CMakeDemo_INCLUDE_DIR source_file.hpp
    HINTS
        ${CMakeDemo_DIR}
        ${CMAKEDEMO_DIR}
        $ENV{CMakeDemo_DIR}
        $ENV{CMAKEDEMO_DIR}
        ENV CMAKEDEMO_DIR
    PATHS
        /usr
        /usr/local
        /usr/local/include
        ${inc_glob}
    PATH_SUFFIXES 
        include
)


# Check that both the paths to the include and library directory were found.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CMakeDemo
    "\nCMakeDemo not found --- You can download it using:\n\tgit clone 
    https://github.com/mmorse1217/cmake-project-template\n and setting the CMAKEDEMO_DIR environment variable accordingly"
    CMakeDemo_LIB CMakeDemo_INCLUDE_DIR)

# These variables don't show up in the GUI version of CMake. Not required but
# people seem to do this...
mark_as_advanced(CMakeDemo_INCLUDE_DIR CMakeDemo_LIB)

# Finish defining the variables specified above. Variables names here follow
# CMake convention.
set(CMakeDemo_INCLUDE_DIRS ${CMakeDemo_INCLUDE_DIR})
set(CMakeDemo_LIBRARIES ${CMakeDemo_LIB})

# If the above CMake code was successful and we found the library, and there is
# no target defined, lets make one.
if(CMakeDemo_FOUND AND NOT TARGET CMakeDemo::CMakeDemo)
    add_library(CMakeDemo::CMakeDemo UNKNOWN IMPORTED)
    # Set location of interface include directory, i.e., the directory
    # containing the header files for the installed library
    set_target_properties(CMakeDemo::CMakeDemo PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${CMakeDemo_INCLUDE_DIRS}"
        )

    # Set location of the installed library
    set_target_properties(CMakeDemo::CMakeDemo PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
        IMPORTED_LOCATION "${CMakeDemo_LIBRARIES}"
        )
endif()
