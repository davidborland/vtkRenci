/*=========================================================================

  Name:        vtkSmoothDataArray.cxx

  Author:      David Borland, The Renaissance Computing Institute (RENCI)

  Copyright:   The Renaissance Computing Institute (RENCI)

  License:     Licensed under the RENCI Open Source Software License v. 1.0.
               
               See included License.txt or 
               http://www.renci.org/resources/open-source-software-license
               for details.

=========================================================================*/

#include "vtkSmoothDataArray.h"

#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkDataSet.h"
#include "vtkFieldData.h"
#include "vtkGenericCell.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"

//-----------------------------------------------------------------------------

vtkCxxRevisionMacro(vtkSmoothDataArray, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkSmoothDataArray);

template<class data_type>
void vtkSmoothDataArrayDoComputePoints(vtkDataSet *structure,
                                       data_type *smoothed,
                                       data_type *temp,                                       
                                       int anisotropic,
                                       double coefficient);

template<class data_type>
void vtkSmoothDataArrayDoComputeCells(vtkDataSet *structure,
                                      data_type *smoothed,
                                      data_type *temp,
                                      int anisotropic,
                                      double coefficient);

//-----------------------------------------------------------------------------

vtkSmoothDataArray::vtkSmoothDataArray()
{
  this->ResultArrayName = NULL;
  this->NumberOfIterations = 100;
  this->AnisotropicDiffusion = 0;
  this->ConductionConstant = 0.05;

  this->SetInputScalars(vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS,
                        vtkDataSetAttributes::SCALARS);
}

vtkSmoothDataArray::~vtkSmoothDataArray()
{
  this->SetResultArrayName(NULL);
}

//-----------------------------------------------------------------------------

void vtkSmoothDataArray::SetInputScalars(int fieldAssociation, const char *name)
{
  if (   (fieldAssociation != vtkDataObject::FIELD_ASSOCIATION_POINTS)
      && (fieldAssociation != vtkDataObject::FIELD_ASSOCIATION_CELLS)
      && (fieldAssociation != vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS) )
    {
    vtkErrorMacro(<<"Input scalars must be associated with points or cells.");
    return;
    }

  this->SetInputArrayToProcess(0, 0, 0, fieldAssociation, name);
}

void vtkSmoothDataArray::SetInputScalars(int fieldAssociation, int fieldAttributeType)
{
  if (   (fieldAssociation != vtkDataObject::FIELD_ASSOCIATION_POINTS)
      && (fieldAssociation != vtkDataObject::FIELD_ASSOCIATION_CELLS)
      && (fieldAssociation != vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS) )
    {
    vtkErrorMacro(<<"Input scalars must be associated with points or cells.");
    return;
    }

  this->SetInputArrayToProcess(0, 0, 0, fieldAssociation, fieldAttributeType);
}

//-----------------------------------------------------------------------------

static int vtkSmoothDataArrayHasArray(vtkFieldData *fieldData,
                                      vtkDataArray *array)
{
  int numarrays = fieldData->GetNumberOfArrays();
  for (int i = 0; i < numarrays; i++)
    {
    if (array == fieldData->GetArray(i))
      {
      return 1;
      }
    }
  return 0;
}

static void vtkSmoothDataArrayComputeCentroid(vtkCell* cell, double centroid[3]) 
{
  vtkPoints* points = cell->GetPoints();
  int numPoints = points->GetNumberOfPoints();

  centroid[0] = centroid[1] = centroid[2] = 0.0;
  
  for (int j = 0; j < numPoints; j++)
    {
    double* p = points->GetPoint(j);
    centroid[0] += p[0];
    centroid[1] += p[1];
    centroid[2] += p[2];
    }
  centroid[0] /= numPoints;
  centroid[1] /= numPoints;
  centroid[2] /= numPoints;
}

//-----------------------------------------------------------------------------

int vtkSmoothDataArray::RequestData(
  vtkInformation *vtkNotUsed(request),                                        
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // Get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // Get the input and ouptut
  vtkDataSet *input
    = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkDataSet *output
    = vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkDataArray *scalars = this->GetInputArrayToProcess(0, inputVector);

  vtkDebugMacro(<<"Smoothing data array");

  // Check scalar data
  if (scalars == NULL)
    {
    vtkErrorMacro(<<"No input scalars.");
    return 0;
    }
  if (scalars->GetNumberOfComponents() != 1)
    {
    vtkErrorMacro(<<"Input scalars must have one component.");
    return 0;
    }

  // Set which type of scalar data
  int fieldAssociation;
  if (vtkSmoothDataArrayHasArray(input->GetPointData(), scalars))
    {
    fieldAssociation = vtkDataObject::FIELD_ASSOCIATION_POINTS;
    }
  else if (vtkSmoothDataArrayHasArray(input->GetCellData(), scalars))
    {
    fieldAssociation = vtkDataObject::FIELD_ASSOCIATION_CELLS;
    }
  else
    {
    vtkErrorMacro("Input scalars do not seem to be either point or cell scalars.");
    return 0;
    }
  
  // Copy input data to output
  output->CopyStructure(input);
  output->GetPointData()->PassData(input->GetPointData());
  output->GetCellData()->PassData(input->GetCellData());

  // Create the output data array
  vtkDataArray* smoothed = vtkDataArray::CreateDataArray(scalars->GetDataType());
  smoothed->SetNumberOfComponents(1);
  smoothed->SetNumberOfTuples(scalars->GetNumberOfTuples());
  smoothed->DeepCopy(scalars);
  char* name;
  if (this->ResultArrayName)
    {
    int nameLength = strlen(this->ResultArrayName);
    name = new char[nameLength];
    sprintf_s(name, nameLength, "%s", this->ResultArrayName);
    }
  else 
    {
    int nameLength = strlen(scalars->GetName()) + strlen(" Smoothed") + 1;
    name = new char[nameLength];
    sprintf_s(name, nameLength, "%s%s", scalars->GetName(), " Smoothed");
    }
  smoothed->SetName(name);

  // Create a temporary data array for calculations
  vtkDataArray* temp = vtkDataArray::CreateDataArray(scalars->GetDataType());
  temp->SetNumberOfComponents(1);
  temp->SetNumberOfTuples(scalars->GetNumberOfTuples());

  if (fieldAssociation == vtkDataObject::FIELD_ASSOCIATION_POINTS)
    {
    for (int i = 0; i < this->NumberOfIterations; i++)
      {
      switch(scalars->GetDataType())
        {
        vtkTemplateMacro(vtkSmoothDataArrayDoComputePoints(
                         input,
                         static_cast<VTK_TT *>(smoothed->GetVoidPointer(0)),
                         static_cast<VTK_TT *>(temp->GetVoidPointer(0)),
                         this->AnisotropicDiffusion,
                         this->ConductionConstant));
        }
      }

    output->GetPointData()->AddArray(smoothed);
    output->GetPointData()->SetActiveScalars(name);
    }
  else 
    {
    for (int i = 0; i < this->NumberOfIterations; i++)
      {
      switch(scalars->GetDataType())
        {
        vtkTemplateMacro(vtkSmoothDataArrayDoComputeCells(
                         input,
                         static_cast<VTK_TT *>(smoothed->GetVoidPointer(0)),
                         static_cast<VTK_TT *>(temp->GetVoidPointer(0)),
                         this->AnisotropicDiffusion,
                         this->ConductionConstant));
        }
      }

    output->GetCellData()->AddArray(smoothed);
    output->GetCellData()->SetActiveScalars(name);
    }

  delete [] name;    
  temp->Delete();

  return 1;
}

//-----------------------------------------------------------------------------

template<class data_type>
void vtkSmoothDataArrayDoComputePoints(vtkDataSet *structure,
                                       data_type *smoothed,
                                       data_type *temp,
                                       int anisotropic,
                                       double coefficient)
{

  vtkIdType numberOfCells = structure->GetNumberOfCells();
  vtkIdType numberOfPoints = structure->GetNumberOfPoints();

  double* weight = new double[numberOfPoints];

  vtkGenericCell* cell = vtkGenericCell::New();

  // Initialize temp array
  for (int i = 0; i < numberOfPoints; i++)
    {
    temp[i] = 0.0;
    weight[i] = 0.0;
    }

  // Compute the gradients along each edge
  for (int i = 0; i < numberOfCells; i++) 
    {
    // Get the edges for this cell
    structure->GetCell(i, cell);
    int numEdges = cell->GetNumberOfEdges();

    for (int j = 0; j < numEdges; j++) 
      {
      // Get the edge
      vtkCell* edge = cell->GetEdge(j);
      
      // Get the two points defining this edge
      int id1 = edge->GetPointId(0);
      int id2 = edge->GetPointId(1);

      // Compute the distance between the two points
      double p1[3];
      double p2[3];
      structure->GetPoint(id1, p1);
      structure->GetPoint(id2, p2);
      double d = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));

      // Perform smoothing based on gradient
      double K = coefficient;
      double g = smoothed[id2] - smoothed[id1];
      double c = anisotropic ? 1.0 / (1.0 + (g / K) * (g / K)) : 1.0;

      double s = (g * c) / d;
      temp[id1] += s;
      temp[id2] -= s;

      weight[id1] += 1.0 / d;
      weight[id2] += 1.0 / d;
      }
    }

  // Normalize and add to current 
  for (int i = 0; i < numberOfPoints; i++) 
    {
    temp[i] = smoothed[i] + temp[i];// / weight[i];
    smoothed[i] = temp[i];
    }

  // Clean up
  delete [] weight;
  cell->Delete();
}

template<class data_type>
void vtkSmoothDataArrayDoComputeCells(vtkDataSet *structure,
                                      data_type *smoothed,
                                      data_type *temp,
                                      int anisotropic,
                                      double coefficient)
{
  vtkIdType numberOfCells = structure->GetNumberOfCells();

  double* weight = new double[numberOfCells];

  vtkGenericCell* cell = vtkGenericCell::New();
  vtkGenericCell* neighbor = vtkGenericCell::New();

  vtkIdList* cellIds = vtkIdList::New();

  // Initialize temp array
  for (int i = 0; i < numberOfCells; i++)
    {
    temp[i] = 0.0;
    weight[i] = 0;
    }

  // Compute the gradients along each edge
  for (int i = 0; i < numberOfCells; i++) 
    {
    structure->GetCell(i, cell);

    // Compute the centroid of this cell
    double p1[3];
    vtkSmoothDataArrayComputeCentroid(cell, p1);

    // Loop over the edges for this cell
    int numEdges = cell->GetNumberOfEdges();
    for (int j = 0; j < numEdges; j++) 
      {
      vtkCell* edge = cell->GetEdge(j);

      // Get the cells that share this edge
      structure->GetCellNeighbors(i, edge->GetPointIds(), cellIds);
      int numCellNeighbors = cellIds->GetNumberOfIds();

      for (int k = 0; k < numCellNeighbors; k++) 
        {
        vtkIdType id = cellIds->GetId(k);
        structure->GetCell(id, neighbor);

        // Compute the centroid of this neighbor cell
        double p2[3];
        vtkSmoothDataArrayComputeCentroid(neighbor, p2);

        // Compute the distance between the two cells
        double d = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));   

        // Perform smoothing based on gradient
        double K = coefficient;
        double g = smoothed[id] - smoothed[i];
        double c = anisotropic ? 1.0 / (1.0 + (g / K) * (g / K)) : 1.0;

        temp[i] += (g * c) / d;
        weight[i] += 1.0 / d;
        }
      }
    }

  // Normalize and add to current 
  for (int i = 0; i < numberOfCells; i++) 
    {
    temp[i] = smoothed[i] + temp[i] / weight[i];
    smoothed[i] = temp[i];
    }

  // Clean up
  delete [] weight;
  cell->Delete();
  neighbor->Delete();
  cellIds->Delete();
}

//-----------------------------------------------------------------------------

void vtkSmoothDataArray::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "ResultArrayName: " << this->ResultArrayName << "\n";
  os << indent << "NumberOfIterations: " << this->NumberOfIterations << "\n";
  os << indent << "AnisotropicDiffusion: " << this->AnisotropicDiffusion << "\n";
  os << indent << "ConductionConstant: " << this->ConductionConstant << "\n";
}