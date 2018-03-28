# - Try to find NURBS++
# Once done this will define
#  NURBS3_FOUND - System has NURBS++
#  NURBS3_INCLUDE_DIRS - The NURBS++ include directories
#  NURBS3_LIBRARIES - The libraries needed to use NURBS++
#  NURBS++_DEFINITIONS - Compiler switches required for using NURBS++


find_package(PkgConfig)
pkg_check_modules(NURBS++ REQUIRED nurbs++-3.0.12)
set(NURBS3_DEFINITIONS ${NURBS++_CFLAGS_OTHER})

# I have hacked in the nurbs++ here to make it work... Couldn't find the .h files otherwise
find_path(NURBS3_INCLUDE_DIR nurbs.h
          HINTS ${NURBS++_INCLUDEDIR}/nurbs++ ${NURBS++_INCLUDE_DIRS}
         )
## Find the 4 libraries required for NURBS
find_library(NURBS3_LIBRARY1 nurbsf libnurbsf 
             HINTS ${NURBS++_LIBDIR} ${NURBS++_LIBRARY_DIRS})
find_library(NURBS3_LIBRARY2 matrixN libmatrixN 
             HINTS ${NURBS++_LIBDIR} ${NURBS++_LIBRARY_DIRS})
find_library(NURBS3_LIBRARY3 matrix libmatrix 
             HINTS ${NURBS++_LIBDIR} ${NURBS++_LIBRARY_DIRS})
find_library(NURBS3_LIBRARY4 matrixI libmatrixI
             HINTS ${NURBS++_LIBDIR} ${NURBS++_LIBRARY_DIRS})                                       

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set NURBS++_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(Nurbs DEFAULT_MSG
                                  NURBS3_LIBRARY NURBS3_INCLUDE_DIR)

mark_as_advanced(NURBS3_INCLUDE_DIR NURBS3_LIBRARY )

## Add all the libraries to the libraries command - so they are all included by including LIBRARIES
set(NURBS3_LIBRARIES ${NURBS3_LIBRARY1};${NURBS3_LIBRARY2};${NURBS3_LIBRARY3};${NURBS3_LIBRARY4} )
set(NURBS3_INCLUDE_DIRS ${NURBS3_INCLUDE_DIR} )