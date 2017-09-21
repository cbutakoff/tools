/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/



//
// GenerateCubesFromLabels
//   Usage: GenerateCubesFromLabels InputVolume Startlabel Endlabel
//          where
//          InputVolume is a meta file containing a 3 volume of
//            discrete labels.
//          StartLabel is the first label to be processed
//          EndLabel is the last label to be processed
//          NOTE: There can be gaps in the labeling. If a label does
//          not exist in the volume, it will be skipped.
//
//
#include <vtkMetaImageReader.h>
#include <vtkImageAccumulate.h>
#include <vtkImageWrapPad.h>
#include <vtkMaskFields.h>
#include <vtkThreshold.h>
#include <vtkTransformFilter.h>
#include <vtkGeometryFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>

#include <vtkTransform.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>
#include <sstream>

#include <vtkDataSetReader.h>
#include <vtkDataSetWriter.h>
#include <vtkDataSetTriangleFilter.h>
#include <vtkImageData.h>
#include <vtkCallbackCommand.h>

 
#include "vtkPolyDataAlgorithm.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"
 


class vtkTestProgressReportFilter : public vtkPolyDataAlgorithm
{
public:
  static vtkTestProgressReportFilter *New();
  vtkTypeMacro(vtkTestProgressReportFilter,vtkAlgorithm);
 
protected:
  vtkTestProgressReportFilter(){}
  ~vtkTestProgressReportFilter() VTK_OVERRIDE {}
 
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) VTK_OVERRIDE;
 
private:
  vtkTestProgressReportFilter(const vtkTestProgressReportFilter&) VTK_DELETE_FUNCTION;
  void operator=(const vtkTestProgressReportFilter&) VTK_DELETE_FUNCTION;
 
};


 
vtkStandardNewMacro(vtkTestProgressReportFilter);
 
int vtkTestProgressReportFilter::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
 
  // Get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
 
 
  // Get the input and ouptut
  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
 
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
 
  for(vtkIdType i = 0; i < input->GetNumberOfPoints(); i++)
  {
    this->UpdateProgress(static_cast<double>(i)/input->GetNumberOfPoints());
  }
 
  output->ShallowCopy(input);
 
  return 1;
}






void beholder(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData) {
    vtkTestProgressReportFilter* testFilter = static_cast<vtkTestProgressReportFilter*> (caller);
    std::cout << caller->GetClassName()<<"> Progress: " << testFilter->GetProgress() << std::endl;
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        cout << "Usage: " << argv[0] << " InputVolume StartLabel EndLabel" << endl;
        return EXIT_FAILURE;
    }


    vtkSmartPointer<vtkCallbackCommand> progressCallback =
            vtkSmartPointer<vtkCallbackCommand>::New();
    progressCallback->SetCallback(beholder);

    // Create all of the classes we will need
    vtkSmartPointer<vtkDataSetReader> reader =
            vtkSmartPointer<vtkDataSetReader>::New();
    vtkSmartPointer<vtkImageAccumulate> histogram =
            vtkSmartPointer<vtkImageAccumulate>::New();
    vtkSmartPointer<vtkMaskFields> scalarsOff =
            vtkSmartPointer<vtkMaskFields>::New();
    vtkSmartPointer<vtkThreshold> selector =
            vtkSmartPointer<vtkThreshold>::New();
    vtkSmartPointer<vtkDataSetWriter> writer =
            vtkSmartPointer<vtkDataSetWriter>::New();
    vtkSmartPointer<vtkDataSetTriangleFilter> triangulator =
            vtkSmartPointer<vtkDataSetTriangleFilter>::New();


    reader->AddObserver(vtkCommand::ProgressEvent, progressCallback);
    selector->AddObserver(vtkCommand::ProgressEvent, progressCallback);
    triangulator->AddObserver(vtkCommand::ProgressEvent, progressCallback);


    // Define all of the variables
    unsigned int startLabel = atoi(argv[2]);
    if (startLabel > VTK_SHORT_MAX) {
        std::cout << "ERROR: startLabel is larger than " << VTK_SHORT_MAX << std::endl;
        return EXIT_FAILURE;
    }
    unsigned int endLabel = atoi(argv[3]);
    if (endLabel > VTK_SHORT_MAX) {
        std::cout << "ERROR: endLabel is larger than " << VTK_SHORT_MAX << std::endl;
        return EXIT_FAILURE;
    }
    std::string filePrefix = "Cubes";

    // Generate cubes from labels
    // 1) Read the meta file
    // 2) Generate a histogram of the labels
    // 3) Convert point data to cell data
    // 4) Output each cube model into a separate file

    reader->SetFileName(argv[1]);

    histogram->SetInputConnection(reader->GetOutputPort());
    histogram->SetComponentExtent(0, endLabel, 0, 0, 0, 0);
    histogram->SetComponentOrigin(0, 0, 0);
    histogram->SetComponentSpacing(1, 1, 1);
    histogram->Update();


    // Copy the scalar point data of the volume into the scalar cell data

    selector->SetInputConnection(reader->GetOutputPort());
    selector->SetInputArrayToProcess(0, 0, 0,
            vtkDataObject::FIELD_ASSOCIATION_POINTS,
            vtkDataSetAttributes::SCALARS);



    // Strip the scalars from the output
    scalarsOff->SetInputConnection(selector->GetOutputPort());
    scalarsOff->CopyAttributeOff(vtkMaskFields::POINT_DATA,
            vtkDataSetAttributes::SCALARS);
    scalarsOff->CopyAttributeOff(vtkMaskFields::CELL_DATA,
            vtkDataSetAttributes::SCALARS);


    triangulator->SetInputConnection(scalarsOff->GetOutputPort());

    writer->SetInputConnection(triangulator->GetOutputPort());
    writer->SetFileTypeToBinary();

    for (unsigned int i = startLabel; i <= endLabel; i++) {
        // see if the label exists, if not skip it
        double frequency =
                histogram->GetOutput()->GetPointData()->GetScalars()->GetTuple1(i);
        if (frequency == 0.0) {
            continue;
        }

        // select the cells for a given label
        selector->ThresholdBetween(i, i);

        // output the polydata
        std::stringstream ss;
        ss << filePrefix << i << ".vtk";
        cout << argv[0] << " writing " << ss.str() << endl;

        writer->SetFileName(ss.str().c_str());
        writer->Write();

    }
    return EXIT_SUCCESS;
} 