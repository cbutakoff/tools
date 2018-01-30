/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
// Given a mesh (unstructured grid) and and images , 
// create a pointdata array by sampling the image with image values

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

void AddScalarArray(vtkUnstructuredGrid* mesh, vtkImageData* image, vtkFloatArray* array);

int main(int argc, char** argv)
{
    if(argc<4)
    {
        std::cout<<"Usage: SampleImage2UG image.vtk mesh.vtk out.vtk array_name"<<std::endl;
        std::cout<<"Image is stored in a float array for now"<<std::endl;
        std::cout<<"The image is sampled at vertex positions"<<std::endl;
        
        exit(-1);
    }

    
    int c=1;
    const char* imagefilename = argv[c++];
    const char* meshfilename = argv[c++];
    const char* outfilename = argv[c++];
    const char * arrayname = argv[c++];

    
    std::cout<<"Mesh: "<<meshfilename<<std::endl;
    std::cout<<"Output: "<<outfilename<<std::endl;
    std::cout<<"Image: "<<imagefilename<<std::endl;
    std::cout<<"Array to create: "<<arrayname<<std::endl;
    
    
    
    
    vtkSmartPointer<vtkDataSetReader> mesh_reader = vtkSmartPointer<vtkDataSetReader>::New();
    CommonTools::AssociateProgressFunction(mesh_reader);

    std::cout<<"Reading mesh: "<<meshfilename<<std::endl;
    mesh_reader->SetFileName(meshfilename);
    mesh_reader->Update();


    
    vtkUnstructuredGrid* mesh = vtkUnstructuredGrid::SafeDownCast(mesh_reader->GetOutput());
    
    
    vtkSmartPointer<vtkFloatArray> scalars = vtkSmartPointer<vtkFloatArray>::New();
    scalars->SetName(arrayname);
    scalars->SetNumberOfComponents(1);
    scalars->SetNumberOfTuples(mesh->GetNumberOfPoints());
    
    vtkSmartPointer<vtkDataSetReader> image_reader = vtkSmartPointer<vtkDataSetReader>::New();
    CommonTools::AssociateProgressFunction(mesh_reader);

    std::cout<<"Reading image"<<std::endl;
    image_reader->SetFileName(imagefilename);
    image_reader->Update();


    AddScalarArray(mesh, vtkImageData::SafeDownCast(image_reader->GetOutput()), scalars);
    
    mesh->GetPointData()->AddArray(scalars);
    
    vtkSmartPointer<vtkUnstructuredGridWriter> wr = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    CommonTools::AssociateProgressFunction(wr);
    std::cout<<"Output filename: "<<outfilename<<std::endl;
    wr->SetFileName(outfilename);
    wr->SetFileTypeToBinary();
    wr->SetInputData(mesh);
    wr->Write();
}




void AddScalarArray(vtkUnstructuredGrid* mesh, vtkImageData* image, vtkFloatArray* array)
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
        array->SetTuple1(i, samples->GetTuple1(i));
    }
}
