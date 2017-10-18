/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
// Given a mesh (unstructured grid) and 3 images with x,y,z fibers, 
// create fibers at mesh vertices by sampling the images
// images are loaded sequentially, to reduce memory usage

#include <vtkSmartPointer.h>
#include <vtkDataSetReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkProbeFilter.h>
#include <vtkImageData.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>

#include <VTKCommonTools.h>
#include <vtkCallbackCommand.h>

#include <iostream>

void AddCoordinate(vtkUnstructuredGrid* mesh, vtkImageData* image, vtkFloatArray* array, int coord_id);

int main(int argc, char** argv)
{
    if(argc<5)
    {
        std::cout<<"Usage: AddFibers2Mesh mesh.vtk image_x.vtk image_y.vtk image_z.vtk out.vtk"<<std::endl;
        exit(-1);
    }

    const char * arrayname = "Fibers";
    
    int c=1;
    const char* meshfilename = argv[c++];
    const char* xfilename = argv[c++];
    const char* yfilename = argv[c++];
    const char* zfilename = argv[c++];
    const char* outfilename = argv[c++];
    const char * filenames[3] = {xfilename, yfilename, zfilename};

    
    std::cout<<"Mesh: "<<meshfilename<<std::endl;
    std::cout<<"Output: "<<outfilename<<std::endl;
    std::cout<<"Image x: "<<filenames[0]<<std::endl;
    std::cout<<"Image y: "<<filenames[1]<<std::endl;
    std::cout<<"Image z: "<<filenames[2]<<std::endl;
    
    
    
    vtkSmartPointer<vtkDataSetReader> mesh_reader = vtkSmartPointer<vtkDataSetReader>::New();
    CommonTools::AssociateProgressFunction(mesh_reader);

    std::cout<<"Reading mesh: "<<meshfilename<<std::endl;
    mesh_reader->SetFileName(meshfilename);
    mesh_reader->Update();


    
    vtkUnstructuredGrid* mesh = vtkUnstructuredGrid::SafeDownCast(mesh_reader->GetOutput());
    
    
    vtkSmartPointer<vtkFloatArray> fibers = vtkSmartPointer<vtkFloatArray>::New();
    fibers->SetName("Fibers");
    fibers->SetNumberOfComponents(3);
    fibers->SetNumberOfTuples(mesh->GetNumberOfPoints());
    
    for(int i=0; i<3; i++)
    {
        vtkSmartPointer<vtkDataSetReader> image_reader = vtkSmartPointer<vtkDataSetReader>::New();
        CommonTools::AssociateProgressFunction(mesh_reader);

        std::cout<<"Reading image"<<filenames[i]<<std::endl;
        image_reader->SetFileName(filenames[i]);
        image_reader->Update();

        
        AddCoordinate(mesh, vtkImageData::SafeDownCast(image_reader->GetOutput()), fibers, i);
    }
    
    mesh->GetPointData()->AddArray(fibers);
    
    vtkSmartPointer<vtkUnstructuredGridWriter> wr = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    CommonTools::AssociateProgressFunction(wr);
    std::cout<<"Output filename: "<<outfilename<<std::endl;
    wr->SetFileName(outfilename);
    wr->SetFileTypeToBinary();
    wr->SetInputData(mesh);
    wr->Write();
}




void AddCoordinate(vtkUnstructuredGrid* mesh, vtkImageData* image, vtkFloatArray* array, int coord_id)
{
    vtkSmartPointer<vtkProbeFilter> prober = vtkSmartPointer<vtkProbeFilter>::New();
    CommonTools::AssociateProgressFunction(prober);

    std::cout<<"Sampling the image"<<std::endl;
    prober->SetInputData(mesh);
    prober->SetSourceData(image);
    prober->Update();
    
    vtkDataArray *samples = prober->GetOutput()->GetPointData()->GetScalars();
    
    std::cout<<"Copying the samples"<<std::endl;
    for(vtkIdType i=0; i<array->GetNumberOfTuples(); i++)
    {
        double tuple[3];
        array->GetTuple(i, tuple);
        tuple[coord_id] = samples->GetTuple1(i);
        array->SetTuple(i, tuple);
    }
}
