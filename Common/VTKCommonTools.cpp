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



void CommonTools::SaveVolMeshBSC(const char* infile, const char* outfile_prefix, float scale, bool correct_orientation)
{
    vtkSmartPointer<vtkDataSetReader> reader = vtkSmartPointer<vtkDataSetReader>::New();
    reader->SetFileName(infile);
    reader->Update();
    
    vtkDataSet* volmesh = reader->GetOutput(); 

    CommonTools::SaveVolMeshBSC(volmesh, outfile_prefix, scale);
}





void CommonTools::SaveVolMeshBSC(vtkDataSet* volmesh, const char* outfile_prefix, float scale, bool correct_orientation)
{
    std::string outfile_nodes(outfile_prefix);
    std::string outfile_elements(outfile_prefix);
    std::string outfile_gradient(outfile_prefix);
    std::string outfile_elementtype(outfile_prefix);
    outfile_nodes += ".node";
    outfile_elements += ".ele";
    outfile_gradient += ".grad";
    outfile_elementtype += ".elem_type";    

    std::cout<<outfile_nodes<<std::endl;
    std::cout<<outfile_elements<<std::endl;
    std::cout<<"DON'T FORGET TO CHECK ELEMENT ORIENTATION!"<<std::endl;
    
    //vtkSmartPointer<vtkDataSetReader> reader = vtkSmartPointer<vtkDataSetReader>::New();
    //reader->SetFileName(infile);
    //reader->Update();
    
    //vtkDataSet* volmesh = reader->GetOutput(); 
    
    std::ofstream node_file;
    node_file.open(outfile_nodes.c_str(),ios::trunc);
    node_file<< "COORDINATES"<<std::endl;
    
    for(int i=0; i<volmesh->GetNumberOfPoints(); i++)
    {
        double *pt = volmesh->GetPoint(i);
        node_file<<i+1<<" "<< pt[0]*scale<<" "<< pt[1]*scale<< " "<<pt[2]*scale<<LINEBREAK;
    }

    node_file<< "END_COORDINATES"<<std::endl;
    node_file.close();
    
    
    std::ofstream ele_file;
    ele_file.open(outfile_elements.c_str(),ios::trunc);
    ele_file<< "ELEMENTS"<<std::endl;

    /*https://cgns.github.io/CGNS_docs_current/sids/conv.html */
    double hex_order[8] = {0,1,3,2,4,5,6,7};
    double wedge_order[6] = {2-1,3-1,1-1,6-1,4-1,5-1};
    double pyramid_order[5] = {1-1,4-1,2-1,3-1,5-1};

    for(int i=0; i<volmesh->GetNumberOfCells(); i++)
    {
        vtkCell *cell = volmesh->GetCell(i);

        if(cell->GetNumberOfPoints()==4)
        {
            if( correct_orientation )
                ele_file<<i+1<<" "<< cell->GetPointId(3)+1 <<" "<< cell->GetPointId(1)+1 << " "
                    <<cell->GetPointId(2)+1<<" "<<cell->GetPointId(0)+1<<LINEBREAK; //flipping elements
            else
                ele_file<<i+1<<" "<< cell->GetPointId(0)+1 <<" "<< cell->GetPointId(1)+1 << " "
                    <<cell->GetPointId(2)+1<<" "<<cell->GetPointId(3)+1<<LINEBREAK; //flipping elements
        }

        else if(cell->GetNumberOfPoints()==8)
        {
            ele_file<<i+1;
            if( correct_orientation )
                for(int k=0; k<4; k++)
                    ele_file <<" "<< cell->GetPointId(hex_order[k])+1;
            else
                for(int k=0; k<4; k++)
                    ele_file <<" "<< cell->GetPointId(k)+1;

            ele_file<<LINEBREAK; 
        }

        else if(cell->GetNumberOfPoints()==6)
        {

            ele_file<<i+1;
            if(correct_orientation)
                for(int k=0; k<6; k++)
                    ele_file <<" "<< cell->GetPointId(wedge_order[k])+1;
            else
                for(int k=0; k<6; k++)
                    ele_file <<" "<< cell->GetPointId(k)+1;

            ele_file<<LINEBREAK; 
        }
        else if(cell->GetNumberOfPoints()==5)
        {
            ele_file<<i+1;
            if(correct_orientation)
                for(int k=0; k<cell->GetNumberOfPoints(); k++)
                    ele_file <<" "<< cell->GetPointId(pyramid_order[k])+1;
            else
                for(int k=0; k<cell->GetNumberOfPoints(); k++)
                    ele_file <<" "<< cell->GetPointId(k)+1;

            ele_file<<LINEBREAK; 
        }
        else
            cout<<"element "<<i<<" has "<<cell->GetNumberOfPoints()<<" points and is not supported."<<endl;

    }

    ele_file<< "END_ELEMENTS"<<std::endl;
    ele_file.close();
    
    vtkFloatArray* grads = (vtkFloatArray*)volmesh->GetPointData()->GetArray("Fibers");
    if(grads!=NULL)
    {
        std::cout<<outfile_gradient<<std::endl;
        
        std::ofstream grad_file;
        grad_file.open(outfile_gradient.c_str(),ios::trunc);

        for(int i=0; i<grads->GetNumberOfTuples(); i++)
        {
            double *tuple = grads->GetTuple(i);
            double length =  sqrt( tuple[0]*tuple[0]+ tuple[1]*tuple[1]+ tuple[2]*tuple[2]);
            grad_file<< i+1<<" "<< tuple[0]/length <<" "<< tuple[1]/length << " " << tuple[2]/length<<LINEBREAK;
        }

        grad_file.close();
    }

    {
        std::ofstream type_file;
        std::cout<<"Saving element types (number of vertices per element)"<<std::endl;
        type_file.open(outfile_elementtype.c_str(),ios::trunc);


        for(int i=0; i<volmesh->GetNumberOfCells(); i++)
        {
            vtkCell *cell = volmesh->GetCell(i);

            type_file<<i+1<<" "<< cell->GetNumberOfPoints()<<LINEBREAK; 

        }

        type_file.close();

    }

}
