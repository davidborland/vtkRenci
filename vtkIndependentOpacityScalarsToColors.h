/*=========================================================================

  Name:        vtkIndependentOpacityScalarsToColors.h

  Author:      David Borland, The Renaissance Computing Institute (RENCI)

  Copyright:   The Renaissance Computing Institute (RENCI)

  License:     Licensed under the RENCI Open Source Software License v. 1.0.
               
               See included License.txt or 
               http://www.renci.org/resources/open-source-software-license
               for details.

=========================================================================*/
// .NAME vtkIndependentOpacityScalarsToColors
// .SECTION Description
// vtkIndependentOpacityScalarsToColors assumes a two-component data
// array.  The first component is mapped to rgb color via the supplied 
// lookup table (opacity is ignored), and the second component is mapped 
// to opacity via another lookup table (rgb is ignored).
//
// .SECTION see also 
// vtkMergeFields

#ifndef __vtkIndependentOpacityScalarsToColors
#define __vtkIndependentOpacityScalarsToColors

#include "vtkMolConfigure.h"

#include "vtkScalarsToColors.h"

class VTK_MOL_EXPORT vtkIndependentOpacityScalarsToColors : public vtkScalarsToColors
{
public:
  vtkTypeRevisionMacro(vtkIndependentOpacityScalarsToColors,vtkScalarsToColors);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkIndependentOpacityScalarsToColors *New();

  // Description:
  // Set/get the separate lookup tables to use for rgb and opacity
  vtkSetObjectMacro(ColorLookupTable,vtkScalarsToColors);
  vtkGetObjectMacro(ColorLookupTable,vtkScalarsToColors);
  vtkSetObjectMacro(OpacityLookupTable,vtkScalarsToColors);
  vtkGetObjectMacro(OpacityLookupTable,vtkScalarsToColors);

  // Description:
  // Implement parent-class methods
  virtual double* GetRange();
  virtual void SetRange(double min, double max);
  virtual unsigned char* MapValue(double v);
  virtual void GetColor(double v, double rgb[3]);
  virtual double GetOpacity(double v);

  // Description:
  // Build the member lookup tables
  virtual void Build();
  
  // Description:
  // An internal method maps a data array into a 4-component, unsigned char
  // RGBA array, using one component for rgb, and the other for alpha.
  virtual vtkUnsignedCharArray* MapScalars(vtkDataArray* scalars, int colorMode,
                                           int component);

  // Description:
  // An internal method typically not used in applications.
  virtual void MapScalarsThroughTable2(void* input, unsigned char* output,
                                       int inputDataType, int numberOfValues,
                                       int inputIncrement, int outputFormat);
  
protected:
  vtkIndependentOpacityScalarsToColors();
  ~vtkIndependentOpacityScalarsToColors();

  vtkScalarsToColors* ColorLookupTable;
  vtkScalarsToColors* OpacityLookupTable;

private:
  vtkIndependentOpacityScalarsToColors(const vtkIndependentOpacityScalarsToColors&);  // Not implemented.
  void operator=(const vtkIndependentOpacityScalarsToColors&);  // Not implemented.
};

#endif
