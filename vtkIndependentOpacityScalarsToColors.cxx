/*=========================================================================

  Name:        vtkIndependentOpacityScalarsToColors.cxx

  Author:      David Borland, The Renaissance Computing Institute (RENCI)

  Copyright:   The Renaissance Computing Institute (RENCI)

  License:     Licensed under the RENCI Open Source Software License v. 1.0.
               
               See included License.txt or 
               http://www.renci.org/resources/open-source-software-license
               for details.

=========================================================================*/

#include "vtkIndependentOpacityScalarsToColors.h"

#include "vtkLookupTable.h"
#include "vtkObjectFactory.h"
#include "vtkUnsignedCharArray.h"

//-----------------------------------------------------------------------------

vtkCxxRevisionMacro(vtkIndependentOpacityScalarsToColors, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkIndependentOpacityScalarsToColors);

//-----------------------------------------------------------------------------

vtkIndependentOpacityScalarsToColors::vtkIndependentOpacityScalarsToColors()
{
  // Start with default lookup tables so we don't have to do a bunch of 
  // non-NULL checks later
  this->ColorLookupTable = vtkLookupTable::New();
  this->OpacityLookupTable = vtkLookupTable::New();
}

vtkIndependentOpacityScalarsToColors::~vtkIndependentOpacityScalarsToColors()
{
  this->SetColorLookupTable(NULL);
  this->SetOpacityLookupTable(NULL);
}

//-----------------------------------------------------------------------------

double* vtkIndependentOpacityScalarsToColors::GetRange()
{
  // Return the range of the color lookup table
  return this->ColorLookupTable->GetRange();
}

void vtkIndependentOpacityScalarsToColors::SetRange(double min, double max)
{
  // Set range of both tables
  this->ColorLookupTable->SetRange(min, max);
  this->OpacityLookupTable->SetRange(min, max);
}

unsigned char* vtkIndependentOpacityScalarsToColors::MapValue(double v)
{
  // Only mapping rgb here
  return this->ColorLookupTable->MapValue(v);
}

void vtkIndependentOpacityScalarsToColors::GetColor(double v, double rgb[3])
{
  this->ColorLookupTable->GetColor(v, rgb);
}

double vtkIndependentOpacityScalarsToColors::GetOpacity(double v)
{
  return this->OpacityLookupTable->GetOpacity(v);
}

//-----------------------------------------------------------------------------

void vtkIndependentOpacityScalarsToColors::Build()
{
  this->ColorLookupTable->Build();
  this->OpacityLookupTable->Build();
}

vtkUnsignedCharArray* vtkIndependentOpacityScalarsToColors::MapScalars(vtkDataArray* scalars, 
                                                                       int colorMode, 
                                                                       int component)
{
  vtkUnsignedCharArray* newColors = NULL;
 
  // Check that the internal lookup tables are set
  if (this->ColorLookupTable == NULL)
    {
    vtkWarningMacro(<<"ColorLookupTable not set.");
    }

  if (this->OpacityLookupTable == NULL)
    {
    vtkWarningMacro(<<"OpacityLookupTable not set.");
    }

  if (component == -1)
    {
    newColors = vtkUnsignedCharArray::New();
    newColors->SetNumberOfComponents(4);
    newColors->SetNumberOfTuples(scalars->GetNumberOfTuples());

    // Map color and opacity
    if (this->ColorLookupTable)
      {
      vtkUnsignedCharArray* color = this->ColorLookupTable->MapScalars(scalars, colorMode, 0);
      for (int i = 0; i < color->GetNumberOfTuples(); i++) 
        {
        newColors->SetComponent(i, 0, color->GetComponent(i, 0));
        newColors->SetComponent(i, 1, color->GetComponent(i, 1));
        newColors->SetComponent(i, 2, color->GetComponent(i, 2));
        }
      color->Delete();
      }

    if (this->OpacityLookupTable)
      {
      vtkUnsignedCharArray* opacity = this->OpacityLookupTable->MapScalars(scalars, colorMode, 1);
      for (int i = 0; i < opacity->GetNumberOfTuples(); i++) 
        {
        newColors->SetComponent(i, 3, opacity->GetComponent(i, 3));
        }
      opacity->Delete();
      }
    }
  else
    {
    if (component == 0)
      {
      newColors = this->ColorLookupTable->MapScalars(scalars, colorMode, component);
      }
    else if (component == 1) 
      {
      newColors = this->OpacityLookupTable->MapScalars(scalars, colorMode, component);
      }
    else 
      {
      vtkWarningMacro(<<"Component is not 0 or 1.  Returning component 0");
      newColors = this->ColorLookupTable->MapScalars(scalars, colorMode, 0);
      }
    }
    
  return newColors;
}

void vtkIndependentOpacityScalarsToColors::MapScalarsThroughTable2(void *input, unsigned char *output,
                                                                   int inputDataType, int numberOfValues,
                                                                   int inputIncrement, int outputFormat)
{
  // Not necessary, as the actual mapping is performed by the 
  // ColorLookupTable and OpacityLookupTable membesr
}

//-----------------------------------------------------------------------------

void vtkIndependentOpacityScalarsToColors::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "ColorLookupTable: ";
  this->ColorLookupTable->PrintSelf(os,indent.GetNextIndent());

  os << indent << "OpacityLookupTable: ";
  this->OpacityLookupTable->PrintSelf(os,indent.GetNextIndent());
}