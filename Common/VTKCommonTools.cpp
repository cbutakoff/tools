/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#include <vtkCleanPolyData.h>
#include <vtkThreshold.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkTriangleFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkDelaunay2D.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataSetWriter.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageShrink3D.h>
#include <vtkImageData.h>
#include <vtkShortArray.h>
#include <vtkPolyData.h>
#include <vtkDataSet.h>
#include <vtkPoints.h>
#include <vtkStripper.h>
#include <vtkCutter.h>
#include <vtkPointLocator.h>
#include <vtkPlane.h>
#include <vtkStringArray.h>
#include <vtkSmartPointer.h>
#include <vtkDataArray.h>
#include <vtkDataSetReader.h>
#include <vtkFeatureEdges.h>
#include <vtkTransformFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPLYReader.h>
#include <vtkSTLReader.h>
#include <vtkPolyDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPLYWriter.h>
#include <vtkSTLWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkIVWriter.h>
#include <vtkCellLocator.h>
#include <vtkUnstructuredGridWriter.h>

#include "VTKCommonTools.h"
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <stdexcept>
#include <string>
#include <algorithm>


#include <VTKCommandLineProgress.h>

//---------------------------------------------------
// OPERATIONS
//---------------------------------------------------



//------------------------------------------------------------------

vtkPolyData* CommonTools::LoadShapeFromFile(const char *shapeFileName)
//------------------------------------------------------------------
{
    vtkPolyData* shapePt = NULL;
    CommonTools::VTKSurfaceMeshFormats vtkSurfaceFormat;
    vtkSurfaceFormat = GetTypeOfVTKData(shapeFileName);

    switch (vtkSurfaceFormat) {
        case CommonTools::PLYType:
        {
            vtkSmartPointer<vtkPLYReader> reader = vtkSmartPointer<vtkPLYReader>::New();
            reader -> SetFileName(shapeFileName);
            reader -> Update();
            if (reader->GetOutput()) {
                shapePt = vtkPolyData::New();
                shapePt -> ShallowCopy(reader->GetOutput());
            }
        }
            break;
        case CommonTools::STLType:
        {
            vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
            reader -> SetFileName(shapeFileName);
            reader -> Update();
            if (reader->GetOutput()) {
                shapePt = vtkPolyData::New();
                shapePt -> ShallowCopy(reader->GetOutput());
            }
        }
            break;
        case CommonTools::VTKPolyDataType:
        {
            vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
            reader -> SetFileName(shapeFileName);
            reader -> Update();
            if (reader->GetOutput()) {
                shapePt = vtkPolyData::New();
                shapePt -> ShallowCopy(reader->GetOutput());
            }
        }
            break;
        case CommonTools::VTKXMLPolyDataType:
        {
            vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
            reader -> SetFileName(shapeFileName);
            reader -> Update();
            if (reader->GetOutput()) {
                shapePt = vtkPolyData::New();
                shapePt -> ShallowCopy(reader->GetOutput());
            }
        }
            break;

        default:
            break;
    }

    return shapePt;
}
//------------------------------------------------------------------

vtkUnstructuredGrid* CommonTools::LoadVolumeFromFile(const char *volumeFileName)
//------------------------------------------------------------------
{
    vtkUnstructuredGrid* VolPt = NULL;
    CommonTools::VTKVolumeMeshFormats vtkVolumeFormat;
    vtkVolumeFormat = GetTypeOfVTKVolumeData(volumeFileName);

    switch (vtkVolumeFormat) {
        case CommonTools::VTKUnstructuredGridType:
        {
            vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
            reader -> SetFileName(volumeFileName);
            reader->ReadAllScalarsOn();
            reader->ReadAllVectorsOn();
            reader -> Update();
            if (reader->GetOutput()->GetNumberOfPoints()) {
                VolPt = vtkUnstructuredGrid::New();
                VolPt -> ShallowCopy(reader->GetOutput());
            }
        }
            break;
        default:
            break;
    }

    return VolPt;
}

//------------------------------------------------------------------

void CommonTools::SaveVolumeToFile(
        vtkUnstructuredGrid *volumePt,
        const char *volumeFileName,
        const char *header)
//------------------------------------------------------------------
{
    std::string str_volumeFileName = volumeFileName;

    std::string ext = str_volumeFileName.substr(str_volumeFileName.find_last_of("."));
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    if (ext == ".vtk") {
        vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
        writer -> SetFileTypeToBinary();
        writer -> SetFileName(volumeFileName);
        writer -> SetInputData(volumePt);
        if (header != NULL) {
            writer->SetHeader(header);
        }
        writer -> Write();
    }
}
//


//------------------------------------------------------------------

void CommonTools::SaveShapeToFile(
        vtkPolyData *shapePt,
        const char *shapeFileName,
        const char *header)
//------------------------------------------------------------------
{
    std::string str_volumeFileName = shapeFileName;

    std::string ext = str_volumeFileName.substr(str_volumeFileName.find_last_of("."));
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);


    if (ext == ".ply") {
        vtkSmartPointer<vtkPLYWriter> writer = vtkSmartPointer<vtkPLYWriter>::New();
        writer -> SetFileName(shapeFileName);
        writer -> SetFileTypeToBinary();
        writer -> SetInputData(shapePt);
        writer -> Write();
    }
    if (ext == ".stl") {
        vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
        writer -> SetFileName(shapeFileName);
        writer -> SetFileTypeToBinary();
        writer -> SetInputData(shapePt);
        writer -> Write();
    }
    if (ext == ".vtk") {
        vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        writer -> SetFileName(shapeFileName);
        writer -> SetFileTypeToBinary();
        writer -> SetInputData(shapePt);
        writer -> Write();
    }
    if (ext == ".vtp") {
        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer -> SetFileName(shapeFileName);
        writer -> SetInputData(shapePt);
        writer -> SetDataModeToAscii();
        writer -> Write();
    }
    if (ext == ".iv") {
        vtkSmartPointer<vtkIVWriter> writer = vtkSmartPointer<vtkIVWriter>::New();
        writer -> SetFileName(shapeFileName);
        writer -> SetInputData(shapePt);
        writer -> Write();
    }

}


//-------------------------------------------------------------------------------

CommonTools::VTKVolumeMeshFormats CommonTools::GetTypeOfVTKVolumeData(const char *volumeFileName)
//-------------------------------------------------------------------------------
{
    CommonTools::VTKVolumeMeshFormats type = CommonTools::UnknownVolumeType;

    std::string str_volumeFileName = volumeFileName;

    std::string ext = str_volumeFileName.substr(str_volumeFileName.find_last_of("."));
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

    if (ext == ".vtk") {
        vtkSmartPointer<vtkDataReader> reader = vtkSmartPointer<vtkDataReader>::New();
        reader->SetFileName(volumeFileName);
        if (reader->IsFileUnstructuredGrid())
            type = CommonTools::VTKUnstructuredGridType;
    }

    return type;
}

//------------------------------------------------------------------

bool CommonTools::CheckSaveFileExtension(const char *shapeFileName)
//------------------------------------------------------------------
{
    bool bRes = false;

    if (strlen(shapeFileName) < 4) {
        return false;
    }

    bRes |= strcmp(shapeFileName + strlen(shapeFileName) - 4, ".stl") == 0;
    bRes |= strcmp(shapeFileName + strlen(shapeFileName) - 4, ".vtk") == 0;
    bRes |= strcmp(shapeFileName + strlen(shapeFileName) - 4, ".vtp") == 0;
    bRes |= strcmp(shapeFileName + strlen(shapeFileName) - 3, ".iv") == 0;

    return bRes;
}


//-------------------------------------------------------------------------------

CommonTools::VTKSurfaceMeshFormats CommonTools::GetTypeOfVTKData(const char *shapeFileName)
//-------------------------------------------------------------------------------
{
    CommonTools::VTKSurfaceMeshFormats type = CommonTools::UnknownType;


    std::string str_volumeFileName = shapeFileName;

    std::string ext = str_volumeFileName.substr(str_volumeFileName.find_last_of("."));
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);



    if (ext == ".stl") {
        vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
        reader->SetFileName(shapeFileName);
        reader->Update();
        if (reader->GetOutput() != NULL)
            type = CommonTools::STLType;
    }
    if (ext == ".vtk") {
        vtkSmartPointer<vtkDataReader> reader = vtkSmartPointer<vtkDataReader>::New();
        reader->SetFileName(shapeFileName);
        if (reader->IsFilePolyData())
            type = CommonTools::VTKPolyDataType;
        //if(reader->IsFileUnstructuredGrid())
        //  type = CommonTools::VTKUnstructuredGridType;
    }
    if (ext == ".ply") {
        //vtkSmartPointer<vtkDataReader> reader = vtkSmartPointer<vtkDataReader>::New();
        //reader->SetFileName (shapeFileName);
        //if(reader->IsFilePolyData())
        type = CommonTools::PLYType;
        //if(reader->IsFileUnstructuredGrid())
        //  type = CommonTools::VTKUnstructuredGridType;
    }
    if (ext == ".vtp") {
        vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        reader->SetFileName(shapeFileName);
        if (reader->GetOutput())
            type = CommonTools::VTKXMLPolyDataType;
    }

    return type;
}

void CommonTools::CommandLineProgressIndicator(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData) {
    vtkTestProgressReportFilter* testFilter = static_cast<vtkTestProgressReportFilter*> (caller);
    std::cout << caller->GetClassName() << "> Progress: " << testFilter->GetProgress() << std::endl;
}

vtkSmartPointer<vtkCallbackCommand>
CommonTools::AssociateProgressFunction(vtkAlgorithm* algorithm, vtkCallbackCommand* existing_function) {
    vtkSmartPointer<vtkCallbackCommand> progressCallback;

    if (existing_function == NULL) {
        progressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
        progressCallback->SetCallback(CommonTools::CommandLineProgressIndicator);
    } else
        progressCallback = existing_function;


    algorithm->AddObserver(vtkCommand::ProgressEvent, progressCallback);
    return progressCallback;
}