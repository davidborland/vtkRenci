/*=========================================================================

  Name:        vtkSmoothDataArray.h

  Author:      David Borland, The Renaissance Computing Institute (RENCI)

  Copyright:   The Renaissance Computing Institute (RENCI)

  License:     Licensed under the RENCI Open Source Software License v. 1.0.
               
               See included License.txt or 
               http://www.renci.org/resources/open-source-software-license
               for details.

=========================================================================*/
// .NAME vtkSmoothDataArray
// .SECTION Description
// vtkSmoothDataArray smooths the active point or cell data by solving
// the the diffusion equation on the underlying dataset.  Either 
// isotropic diffusion, equivalent to filtering with a symmetric 
// Gaussian kernel, or anisotropic diffusion, via the original 
// Perona and Malik method, is supported.  Anisotropic diffusion smooths 
// more in homogenous regions, and less across regions of high gradient.

// .SECTION see also
// vtkDataSet

#ifndef __vtkSmoothDataArray_h
#define __vtkSmoothDataArray_h

#include "vtkMolConfigure.h"

#include "vtkDataSetAlgorithm.h"

class vtkCell;

class VTK_MOL_EXPORT vtkSmoothDataArray : public vtkDataSetAlgorithm
{
public:
  vtkTypeRevisionMacro(vtkSmoothDataArray,vtkDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkSmoothDataArray *New();

  // Description:
  // These are convenience methods that call SetInputArrayToProcess
  // to set the array used as the input scalars.  The fieldAssociation comes
  // from the vtkDataObject::FieldAssocations enum.  The fieldAttributeType
  // comes from the vtkDataSetAttributes::AttributeTypes enum.
  virtual void SetInputScalars(int fieldAssociation, const char *name);
  virtual void SetInputScalars(int fieldAssociation, int fieldAttributeType);

  // Description:
  // Set/get the name of the resulting array to create.  If NULL (the
  // default) then the output array will be the input array name appended 
  // with "Smoothed", e.g. "Scalar Data" -> "Scalar Data Smoothed"
  vtkSetStringMacro(ResultArrayName);
  vtkGetStringMacro(ResultArrayName);

  // Description:
  // Set/get the number of iterations to use when smoothing.  More 
  // iterations will result in more smoothing.
  vtkSetMacro(NumberOfIterations,int);
  vtkGetMacro(NumberOfIterations,int);

  // Description:
  // Use anisotropic diffusion or not.
  vtkSetMacro(AnisotropicDiffusion,int);
  vtkGetMacro(AnisotropicDiffusion,int);
  vtkBooleanMacro(AnisotropicDiffusion,int);

  // Description:
  // Parameter controlling how much edge-weighting to use for
  // anisotropic diffusion.
  vtkSetMacro(ConductionConstant,double);
  vtkGetMacro(ConductionConstant,double);

protected:
  vtkSmoothDataArray();
  ~vtkSmoothDataArray();

  // Usual data generation method
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

  char* ResultArrayName;

  int NumberOfIterations;

  int AnisotropicDiffusion;

  double ConductionConstant;

private:
  vtkSmoothDataArray(const vtkSmoothDataArray&);  // Not implemented.
  void operator=(const vtkSmoothDataArray&);  // Not implemented.
};

#endif
