# It's important to specify the minimum CMake version upfront required by
# CMakeLists.txt. This is so that a user can clearly understand the reason the 
# build will fail before the build actually occurs, instead of searching for the
# CMake function that was used that is causing the failure.
cmake_minimum_required(VERSION 3.10)

# Set the global package-wide C++ standard. This will be inherited by all
# targets specified in the project. One can also specify the C++ standard in a
# target-specific manner, using:
#   set_target_properties(foo PROPERTIES CXX_STANDARD 11)
#   target_compile_features(foo PUBLIC cxx_std_14)
# together.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set the project name and version number. This allows for a user of your
# library or tool to specify a particular version when they include it, as in 
#   find_package(LeanVTK 1.0 REQUIRED)
set(LeanVTK_VERSION_MAJOR 2)
set(LeanVTK_VERSION_MINOR 0)
set(LeanVTK_VERSION_PATCH 0)
set(LeanVTK_VERSION ${LeanVTK_VERSION_MAJOR}.${LeanVTK_VERSION_MINOR})
project(LeanVTK VERSION ${LeanVTK_VERSION})
# Broadcasting the version of the library is always handy
add_compile_definitions("LeanVTK_VERSION_MAJOR=${LeanVTK_VERSION_MAJOR}")
add_compile_definitions("LeanVTK_VERSION_MINOR=${LeanVTK_VERSION_MINOR}")
add_compile_definitions("LeanVTK_VERSION_PATCH=${LeanVTK_VERSION_PATCH}")

# enable unit testing via "make test" once the code has been compiled.
include(CTest)

# The user options
set(BUILD_SHARED OFF CACHE BOOL
    "Build a shared library instead of a static one")

# Function to prepend the subdirectory to source files in subdirectories
FUNCTION(PREPEND var )
   SET(listVar "")
   FOREACH(f ${${var}})
       LIST(APPEND listVar "${CMAKE_CURRENT_SOURCE_DIR}/${f}")
   ENDFOREACH(f)
   SET(${var} "${listVar}" PARENT_SCOPE)
ENDFUNCTION(PREPEND)

# Include source code and headers. This calls the CMakeLists.txt in each
# subdirectory. These can define their own libraries, executables, etc. as targets, 
# but here we define all exportable targets in the root CMakeLists.txt.
add_subdirectory(src)
add_subdirectory(include)

# Add the test directory. It is optional and can be disabled during with
#   cmake -DBUILD_TESTING=OFF ..
# To run unit tests produced here, we only need to run:
#   make test
# or
#   ctest 
# In case your tests are printing to console, you can view their output to
# stdout with:
#   ctest -V
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

# Add the library leanvtk as a target, with the contents of src/ and include/
# as dependencies.
if(${BUILD_SHARED})
    add_library(leanvtk SHARED ${LeanVTK_SRC} ${LeanVTK_INC})
else()
    add_library(leanvtk STATIC ${LeanVTK_SRC} ${LeanVTK_INC})
endif()

# generator expressions are needed for the include directories, since installing headers changes the include path
# Specify that leanvtk requires the files located in the include/ directory at
# compile time. This would normally look like 
#   target_include_directories(leanvtk PUBLIC include/)
# PUBLIC means that other libraries including leanvtk should also include the
# directory include/.
# However, there is a catch. If we are installing the project in
# CMAKE_INSTALL_PREFIX, we can't specify include/ in the build directory: we have 
# copied the contents of include to CMAKE_INSTALL_PREFIX/include and we would
# like  other projects to include this directory instead of include/. The follow
# CMake command handles this. $<BUILD_INTERFACE:...> and
# $<INSTALL_INTERFACE:...> are macros whose values change depending on if we are
# simply building the code or if we are installing it.
target_include_directories(leanvtk PUBLIC
   # headers to include when building from source
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include> 
   $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include> 

   # headers to include when installing  (implicitly prefixes with ${CMAKE_INSTALL_PREFIX}).
   $<INSTALL_INTERFACE:include/leanvtk> 
   )

# Install leanvtk in CMAKE_INSTALL_PREFIX (defaults to /usr/local on linux). 
# To change the install location, run 
#   cmake -DCMAKE_INSTALL_PREFIX=<desired-install-path> ..

# install(...) specifies installation rules for the project. It can specify
# location of installed files on the system, user permissions, build
# configurations, etc. Here, we are only copying files.
# install(TARGETS ...) specifies rules for installing targets. 
# Here, we are taking a target or list of targets (leanvtk) and telling CMake
# the following:
#   - put shared libraries associated with leanvtk in ${LeanVTK_LIB_DEST}
#   - put static libraries associated with leanvtk in ${LeanVTK_LIB_DEST}
#   - put include files associated with leanvtk in ${LeanVTK_INCLUDE_DEST}
# We also need to specify the export that is associated with leanvtk; an export 
# is just a list of targets to be installed.
# So we are associating leanvtk with LeanVTKTargets.
install(
    # targets to install
    TARGETS leanvtk 
    # name of the CMake "export group" containing the targets we want to install
    EXPORT LeanVTKTargets
    # Dynamic, static library and include destination locations after running
    # "make install"
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include/leanvtk
)

# install(FILES ...) simply puts files in a certain place with certain
# properties. We're just copying them to the desired place here.
install(FILES ${LeanVTK_INC} DESTINATION include/leanvtk)

# We now need to install the export LeanVTKTargets that we defined above. This
# is needed in order for another project to import leanvtk using 
#   find_package(LeanVTK)
# find_package(LeanVTK) will look for LeanVTKConfig.cmake to provide
# information about the targets contained in the project LeanVTK. Fortunately,
# this is specified in the export LeanVTKTargets, so we will install this too.
# install(EXPORT ...) will install the information about an export. Here, we
# save it to a file {$LeanVTK_LIB_DEST}/LeanVTKTargets.cmake and prepend 
# everything inside LeanVTKTargets  with the namespace leanvtk::.
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/leanvtk/LeanVTKConfigVersion.cmake"
    VERSION ${LeanVTK_VERSION}
    COMPATIBILITY AnyNewerVersion
    )

export(EXPORT LeanVTKTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/leanvtk/LeanVTKTargets.cmake"
  NAMESPACE leanvtk::
)
set(ConfigPackageLocation lib/cmake/leanvtk)
configure_package_config_file(
  ${CMAKE_SOURCE_DIR}/cmake/LeanVTKConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/leanvtk/LeanVTKConfig.cmake
  INSTALL_DESTINATION ${ConfigPackageLocation})

install(EXPORT LeanVTKTargets
  FILE
    LeanVTKTargets.cmake
  NAMESPACE
    leanvtk::
  DESTINATION
    ${ConfigPackageLocation}
)
install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/leanvtk/LeanVTKConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/leanvtk/LeanVTKConfigVersion.cmake"
  DESTINATION
    ${ConfigPackageLocation}
)
