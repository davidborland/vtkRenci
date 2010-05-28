/*=========================================================================

  Name:        vtkPerVertexAmbientOcclusion.h

  Author:      David Borland, The Renaissance Computing Institute (RENCI)

  Copyright:   The Renaissance Computing Institute (RENCI)

  License:     Licensed under the RENCI Open Source Software License v. 1.0.
               
               See included License.txt or 
               http://www.renci.org/resources/open-source-software-license
               for details.

=========================================================================*/
// .NAME vtkPerVertexAmbientOcclusion
// .SECTION Description
// vtkPerVertexAmbientOcclusion computes per-vertex ambient occlusion 
// information for a single vtkPolyData input.  The ambient occlusion
// term is stored as point data.

// .SECTION see also
// vtkPolyData

#ifndef __vtkPerVertexAmbientOcclusion
#define __vtkPerVertexAmbientOcclusion

#include "vtkMolConfigure.h"

#include "vtkPolyDataAlgorithm.h"

class vtkMultiThreader;

class VTK_MOL_EXPORT vtkPerVertexAmbientOcclusion : public vtkPolyDataAlgorithm
{
public:
  vtkTypeRevisionMacro(vtkPerVertexAmbientOcclusion,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkPerVertexAmbientOcclusion *New();

  // Description:
  // Set/get the number of rays cast at each polygon to determine occlusion
  vtkSetMacro(NumberOfSamples,int);
  vtkGetMacro(NumberOfSamples,int);

  // Description:
  // Set/get the number of threads to use
  vtkSetMacro(NumberOfThreads,int);
  vtkGetMacro(NumberOfThreads,int);

protected:
  vtkPerVertexAmbientOcclusion();
  ~vtkPerVertexAmbientOcclusion();

  // Usual data generation method
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

  // Number of rays to cast to determine ambient occlusion at each polygon
  int NumberOfSamples;

  // Number of threads to use for computation
  int NumberOfThreads;

  // Use multi-threading to improve performance
  vtkMultiThreader* Threader;
  static VTK_THREAD_RETURN_TYPE ComputeAmbientOcclusion(void* arg);

  // Description:
  // Compute the rotation matrix of a polygon, given its normal
  static void GetPolygonRotation(const double normal[3], double rotation[3][3]);

  // Description:
  // Compute a cosine-weighted sampling of the unit hemisphere
  static void SampleHemisphere(double sample[3]);

private:
  vtkPerVertexAmbientOcclusion(const vtkPerVertexAmbientOcclusion&);  // Not implemented.
  void operator=(const vtkPerVertexAmbientOcclusion&);  // Not implemented.
};

#endif