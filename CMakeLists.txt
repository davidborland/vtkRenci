CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )

PROJECT( vtkMol )

SET( LIBRARY_OUTPUT_PATH "${vtkMol_BINARY_DIR}/lib" )


#######################################
# Include VTK
#######################################

FIND_PACKAGE( VTK REQUIRED )
INCLUDE( ${VTK_USE_FILE} )  


#######################################
# Static or Shared Libraries
#######################################

# Shared libraries are necessary for wrapping Python
OPTION( BUILD_SHARED_LIBS
        "Build with shared libraries."
        ${VTK_BUILD_SHARED_LIBS} )
# Copy the CMake option to a setting with vtkMol_ prefix for use in
# our project.  This name is used in vtkmsConfigure.h.in.
SET( vtkMol_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS} )


#######################################
# Wrap Python
#######################################

IF( VTK_WRAP_PYTHON )

  OPTION( vtkMol_WRAP_PYTHON
          "Wrap classes in Python."
          ON )

  IF( vtkMol_WRAP_PYTHON )
    SET( VTK_WRAP_PYTHON3_INIT_DIR "${vtkMol_SOURCE_DIR}/Wrapping" )
    INCLUDE( ${VTK_CMAKE_DIR}/vtkWrapPython.cmake )
    IF( WIN32 )
      IF( NOT BUILD_SHARED_LIBS )
        MESSAGE( FATAL_ERROR "Python support requires BUILD_SHARED_LIBS to be ON." )
        SET( vtkMol_CAN_BUILD 0 )
      ENDIF( NOT BUILD_SHARED_LIBS )
    ENDIF( WIN32 )
  ENDIF( vtkMol_WRAP_PYTHON )

ELSE( VTK_WRAP_PYTHON )

  IF( vtkMol_WRAP_PYTHON )
    MESSAGE( "Warning. vtkMol_WRAP_PYTHON is ON but the VTK version you have "
             "chosen does not support Python (VTK_WRAP_PYTHON is OFF).  "
             "Please set vtkMol_WRAP_PYTHON to OFF." )
    SET( vtkMol_WRAP_PYTHON OFF )
  ENDIF( vtkMol_WRAP_PYTHON )

ENDIF( VTK_WRAP_PYTHON )


#######################################
# Configure file
#######################################

CONFIGURE_FILE( ${vtkMol_SOURCE_DIR}/vtkMolConfigure.h.in
                ${vtkMol_BINARY_DIR}/vtkMolConfigure.h )
                
INCLUDE_DIRECTORIES( ${vtkMol_BINARY_DIR} )


#######################################
# Include vtkMol code
#######################################

INCLUDE_DIRECTORIES( ${vtkMol_SOURCE_DIR} )

SET( SRC vtkIndependentOpacityScalarsToColors.h vtkIndependentOpacityScalarsToColors.cxx
         vtkPerVertexAmbientOcclusion.h vtkPerVertexAmbientOcclusion.cxx
         vtkSmoothDataArray.h vtkSmoothDataArray.cxx )
ADD_LIBRARY( vtkMol ${SRC} )


#######################################
# Create Python library
#######################################

IF( VTK_WRAP_PYTHON AND vtkMol_WRAP_PYTHON)
  VTK_WRAP_PYTHON3( vtkMolPython vtkMolPython_SRC "${SRC}")
  ADD_LIBRARY( vtkMolPythonD ${vtkMolPython_SRC} )
  ADD_LIBRARY( vtkMolPython MODULE vtkMolPythonInit.cxx )
  TARGET_LINK_LIBRARIES( vtkMolPythonD vtkMol vtkRenderingPythonD )
  TARGET_LINK_LIBRARIES( vtkMolPython vtkMolPythonD )
ENDIF( VTK_WRAP_PYTHON AND vtkMol_WRAP_PYTHON )