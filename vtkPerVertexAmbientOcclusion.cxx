/*=========================================================================

  Name:        vtkPerVertexAmbientOcclusion.cxx

  Author:      David Borland, The Renaissance Computing Institute (RENCI)

  Copyright:   The Renaissance Computing Institute (RENCI)

  License:     Licensed under the RENCI Open Source Software License v. 1.0.
               
               See included License.txt or 
               http://www.renci.org/resources/open-source-software-license
               for details.

=========================================================================*/

#include "vtkPerVertexAmbientOcclusion.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkGenericCell.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkMultiThreader.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolygon.h"
#include "vtkTriangleStrip.h"

//-----------------------------------------------------------------------------

vtkCxxRevisionMacro(vtkPerVertexAmbientOcclusion, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkPerVertexAmbientOcclusion);

//-----------------------------------------------------------------------------

struct vtkPerVertexAmbientOcclusionThreadStruct
{
  vtkPerVertexAmbientOcclusion* Filter;
  vtkPolyData* PolyData;
  double* AOArray;
  double RayLength;
};

//-----------------------------------------------------------------------------

vtkPerVertexAmbientOcclusion::vtkPerVertexAmbientOcclusion()
{
  this->NumberOfSamples = 64;
  this->Threader = vtkMultiThreader::New();
  this->NumberOfThreads = this->Threader->GetNumberOfThreads();
}

vtkPerVertexAmbientOcclusion::~vtkPerVertexAmbientOcclusion()
{
  this->Threader->Delete();
}

//-----------------------------------------------------------------------------

int vtkPerVertexAmbientOcclusion::RequestData(
  vtkInformation* vtkNotUsed(request),                                        
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // Get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // Get the input and ouptut
  vtkPolyData *input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkDebugMacro(<<"Computing per-vertex ambient occlusion.");

  // Null input check
  if (!input)
    {
    vtkErrorMacro(<<"No input data.");
    return 0;
    }

  // Points check  
  if (input->GetNumberOfPoints() <= 0)
    {
    vtkErrorMacro(<<"No data to generate ambient occlusion for.");
    return 0;
    }

  // If there is nothing to do, pass the data through
  if (input->GetNumberOfPolys() <= 0 && input->GetNumberOfStrips() <= 0)
    {
    output->CopyStructure(input);
    output->GetPointData()->PassData(input->GetPointData());
    output->GetCellData()->PassData(input->GetCellData());
    output->GetFieldData()->PassData(input->GetFieldData());
    return 1;
    }

  // Create temporary data structure, decomposing triangle strips if necessary
  vtkPolyData* tempData = vtkPolyData::New();
  tempData->SetPoints(input->GetPoints());

  vtkCellArray* polys;
  if (input->GetNumberOfStrips() > 0) 
    {
    polys = vtkCellArray::New();

    if (input->GetNumberOfPolys() > 0) 
      {
      polys = vtkCellArray::New();
      polys->DeepCopy(input->GetPolys());
      }
    else
      {
      polys = vtkCellArray::New();
      polys->Allocate(polys->EstimateSize(input->GetNumberOfStrips(), 5));
      }

    vtkIdType npts;
    vtkIdType* pts;
    for (input->GetStrips()->InitTraversal(); input->GetStrips()->GetNextCell(npts, pts); )
      {
      vtkTriangleStrip::DecomposeStrip(npts, pts, polys);
      }

    tempData->SetPolys(polys);
    polys->Delete();
    }
  else 
    {
    tempData->SetPolys(input->GetPolys());
    }

  // Copy data to output
  output->SetVerts(input->GetVerts());
  output->SetLines(input->GetLines());
  output->SetPolys(tempData->GetPolys());
  output->SetPoints(tempData->GetPoints());
  output->GetPointData()->PassData(input->GetPointData());
  output->GetCellData()->PassData(input->GetCellData());
  output->SetFieldData(input->GetFieldData());

  // Build links to enable gather operation for points after computing
  // ambient occlusion for each cell
  tempData->BuildLinks(); 

  // Ambient occlusion array to hold values per polygon
  vtkDoubleArray* aoCell = vtkDoubleArray::New();
  aoCell->SetName("Ambient Occlusion Cell");
  aoCell->SetNumberOfComponents(1);
  aoCell->SetNumberOfTuples(tempData->GetNumberOfCells());
  double* aoCellData = aoCell->GetPointer(0);

  // Set up the struct to pass info to the threads
  vtkPerVertexAmbientOcclusionThreadStruct aoStruct;
  aoStruct.Filter = this;
  aoStruct.PolyData = tempData;
  aoStruct.AOArray = aoCellData;
  aoStruct.RayLength = tempData->GetLength();

  // Some data access methods must be called once from a single thread before they
  // can safely be used. Call those now
  vtkGenericCell* cell = vtkGenericCell::New();
  tempData->GetCell(0, cell);
  cell->Delete();

  // Run multiple threads to do the computation
  this->Threader->SetNumberOfThreads(this->NumberOfThreads);
  this->Threader->SetSingleMethod(vtkPerVertexAmbientOcclusion::ComputeAmbientOcclusion, &aoStruct);
  this->Threader->SingleMethodExecute();
  
  // Create array for ambient occlusion point data
  vtkDoubleArray* ao = vtkDoubleArray::New();
  ao->SetName("Ambient Occlusion");
  ao->SetNumberOfComponents(1);
  ao->SetNumberOfTuples(tempData->GetNumberOfPoints());
  double* aoData = ao->GetPointer(0);

  // Distribute the ambient occlusion data to points from cells
  // XXX:  Should we use vtkCellDataToPointData instead?
  for (int i = 0; i < tempData->GetNumberOfPoints(); i++)
    {
    unsigned short numCells;
    vtkIdType* cells;

    tempData->GetPointCells(i, numCells, cells);

    double sum = 0;
    for (int j = 0; j < numCells; j++) 
      {
      sum += aoCellData[cells[j]];
      }

    aoData[i] = (double)sum / numCells;
    }

  output->GetPointData()->AddArray(ao);
  output->GetPointData()->SetActiveScalars("Ambient Occlusion");

  // Clean up
  tempData->Delete();
  aoCell->Delete();
  ao->Delete();

  return 1;
}

//-----------------------------------------------------------------------------

VTK_THREAD_RETURN_TYPE vtkPerVertexAmbientOcclusion::ComputeAmbientOcclusion(void* arg)
{
  // Get the info out of the input structure
  int threadID = ((vtkMultiThreader::ThreadInfo *)(arg))->ThreadID;
  int numberOfThreads = ((vtkMultiThreader::ThreadInfo *)(arg))->NumberOfThreads;

  vtkPerVertexAmbientOcclusionThreadStruct* threadStruct =
      (vtkPerVertexAmbientOcclusionThreadStruct*)((vtkMultiThreader::ThreadInfo*)arg)->UserData;

  vtkPerVertexAmbientOcclusion* self = threadStruct->Filter;
  vtkPolyData* polyData = threadStruct->PolyData;
  double* aoArray = threadStruct->AOArray;
  double rayLength = threadStruct->RayLength;

  int numberOfSamples = self->GetNumberOfSamples();
  int numberOfCells = polyData->GetNumberOfCells();

  vtkGenericCell* poly = vtkGenericCell::New();

  // Cast rays to compute the ambient occlusion
  vtkIdTypeArray* ptArray = vtkIdTypeArray::New();
  int progressCount = 1;
  const int progressRate = 20;
  for (int i = threadID; i < numberOfCells; i += numberOfThreads)
    {
    // Update progress
    if (threadID == 0 && i > progressCount * numberOfCells / progressRate) 
      {
      double progress = (double)progressCount / progressRate;
      self->UpdateProgress((double)progressCount / progressRate);
      progressCount++;

      if (self->GetAbortExecute())
        {
        break;
        }
      }

    // Check cell type
    polyData->GetCell(i, poly);

    // Compute center and normal
    double p1[3];
    double n[3];
    ptArray->SetArray(poly->GetPointIds()->GetPointer(0), poly->GetNumberOfPoints(), 1);
    vtkPolygon::ComputeCentroid(ptArray, polyData->GetPoints(), p1);
    vtkPolygon::ComputeNormal(ptArray, polyData->GetPoints(), n);

    // Get rotation of this polygon from the normal
    double rotMat[3][3];
    vtkPerVertexAmbientOcclusion::GetPolygonRotation(n, rotMat);
    
    // Cast a number of sample rays
    int numIntersections = 0;
    for (int r = 0; r < numberOfSamples; r++)
      {
      // Obtain a cosine-weighted sampling of the hemisphere
      double p2[3];
      vtkPerVertexAmbientOcclusion::SampleHemisphere(p2);
      
      // Rotate
      vtkMath::Multiply3x3(rotMat, p2, p2);

      // Translate and scale
      p2[0] = p1[0] + p2[0] * rayLength;
      p2[1] = p1[1] + p2[1] * rayLength;
      p2[2] = p1[2] + p2[2] * rayLength;

      // Loop over all polygons checking for ray intersection
      // XXX: Should use an acceleration structure, such as a BSP tree...
      for (int j = 0; j < numberOfCells; j++)
        {
        if (j == i) continue;

        polyData->GetCell(j, poly);

        // Check for intersection
        double t;
        double x[3];
        double pcoords[3];
        int subId;
        if (poly->IntersectWithLine(p1, p2, 0.0, t, x, pcoords, subId))
          {
          numIntersections++;
          break;
          }
        }
      }

    aoArray[i] = (double)numIntersections / numberOfSamples;
    }

  poly->Delete();

  return VTK_THREAD_RETURN_VALUE;
}

//-----------------------------------------------------------------------------

void vtkPerVertexAmbientOcclusion::GetPolygonRotation(const double normal[3], double rotation[3][3])
{
  // Get the quaternion rotating the up vector to this normal
  double up[3] = { 0, 0, 1 };
  double v[3];
  v[0] = (up[0] + normal[0]);
  v[1] = (up[1] + normal[1]);
  v[2] = (up[2] + normal[2]);
  vtkMath::Normalize(v);

  double axis[3];
  vtkMath::Cross(v, normal, axis);
  double dot = vtkMath::Dot(v, normal);

  double quat[4];
  quat[0] = dot;
  quat[1] = axis[0];
  quat[2] = axis[1];
  quat[3] = axis[2];

  // Convert the quaternion to a rotation matrix
  vtkMath::QuaternionToMatrix3x3(quat, rotation);
}

void vtkPerVertexAmbientOcclusion::SampleHemisphere(double sample[3])
{
  double r1 = (double)rand() / RAND_MAX;
  double r2 = (double)rand() / RAND_MAX;

  double v = sqrt(r1);
  double theta = 2 * vtkMath::DoublePi() * r2;

  double x = v * cos(theta);
  double y = v * sin(theta);
  double z = 1 - x * x - y * y;
  z = z > 0.0 ? z : 0.0;

  sample[0] = x;
  sample[1] = y;
  sample[2] = z;
}

//-----------------------------------------------------------------------------

void vtkPerVertexAmbientOcclusion::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "NumberOfSamples: " << this->NumberOfSamples << "\n";
  os << indent << "NumberOfThreads: " << this->NumberOfThreads << "\n";

  os << indent << "Threader: ";
  this->Threader->PrintSelf(os,indent.GetNextIndent());
}