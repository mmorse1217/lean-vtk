#get_filename_component(SELF_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

include(CMakeFindDependencyMacro)
# Capturing values from configure (optional)
#set(my-config-var @my-config-var@)

# Any extra setup

# Add the targets file
include("${CMAKE_CURRENT_LIST_DIR}/LeanVTKTargets.cmake")
