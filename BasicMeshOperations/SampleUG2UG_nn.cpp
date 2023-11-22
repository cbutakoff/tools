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
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkProbeFilter.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellCenters.h>
#include <vtkResampleWithDataSet.h>
#include "vtkOutputWindow.h"
#include <vtkStaticCellLocator.h>
#include <vtkShortArray.h>

#include <iostream>

#include <VTKCommonTools.h>
#include <vtkCallbackCommand.h>


int main(int argc, char** argv)
{
    if(argc<4)
    {
        std::cout<<"Usage: SampleUG2UG meshsrc.vtu meshtgt.vtu arrayname out.vtu"<<std::endl;        
        std::cout<<"Try to use as recent VTK as possible, this uses vtkStaticCellLocator and it used to have bugs. Since VTK 9.3 seems to be ok, and vtkCellLocator is excruciatingly slow since vtk 9.0."<<std::endl;        
        exit(-1);
    }

    vtkOutputWindow::GetInstance()->SetDisplayModeToAlwaysStdErr();
    
    int c=1;
    const char* mesh_src_filename = argv[c++];
    const char* mesh_dst_filename = argv[c++];
    const char* arrayname         = argv[c++];
    const char* output_filename   = argv[c++];
    
    std::cout<<"Src:   "<<mesh_src_filename<<std::endl;
    std::cout<<"Dst:   "<<mesh_dst_filename<<std::endl;
    std::cout<<"Array: "<<arrayname<<std::endl;
    std::cout<<"Out:   "<<output_filename<<std::endl;
    
    
    
    vtkSmartPointer<vtkXMLUnstructuredGridReader> src_reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    CommonTools::AssociateProgressFunction(src_reader);
    std::cout<<"Reading mesh: "<<mesh_src_filename<<std::endl;
    src_reader->SetFileName(mesh_src_filename);
    src_reader->Update();

    vtkSmartPointer<vtkXMLUnstructuredGridReader> dst_reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    CommonTools::AssociateProgressFunction(dst_reader);
    std::cout<<"Reading mesh: "<<mesh_dst_filename<<std::endl;
    dst_reader->SetFileName(mesh_dst_filename);
    dst_reader->Update();

    
    vtkSmartPointer<vtkCellCenters> cc =     vtkSmartPointer<vtkCellCenters>::New();
    CommonTools::AssociateProgressFunction(cc);
    std::cout<<"Running vtkCellCenters"<<std::endl;
    cc->SetInputConnection( dst_reader->GetOutputPort() );
    cc->CopyArraysOff ();
    cc->Update();

    
    printf("Bulding locator\n");
    vtkSmartPointer<vtkStaticCellLocator> loc = vtkSmartPointer<vtkStaticCellLocator>::New();
    loc->SetDataSet(src_reader->GetOutput());
    loc->BuildLocator();

    auto aha_ahamesh = src_reader->GetOutput()->GetCellData()->GetArray(arrayname);

    vtkSmartPointer<vtkShortArray> aha = vtkSmartPointer<vtkShortArray>::New();
    aha->SetName(arrayname);
    aha->SetNumberOfComponents(1);
    aha->SetNumberOfTuples(dst_reader->GetOutput()->GetNumberOfCells());

    printf("Main loop\n");
    for( vtkIdType i=0; i<cc->GetOutput()->GetNumberOfPoints(); i++ ){
        if( i%1000 == 0 ) printf("Processed %lld/%lld cells\r",i, cc->GetOutput()->GetNumberOfPoints());
        double cp[3];
        vtkIdType cellId;
        int subId;
        double dist2;
        auto pt = cc->GetOutput()->GetPoint(i);
        loc->FindClosestPoint(pt, cp, cellId, subId, dist2);
        if(cellId<0){
            printf ("Cell not found for point (%f,%f,%f), cellid=%lld", pt[0], pt[1], pt[2],cellId);
        }

        aha->SetTuple1(i, aha_ahamesh->GetTuple1(cellId));
    }
    printf("\n");


    dst_reader->GetOutput()->GetCellData()->AddArray(aha);



    vtkSmartPointer<vtkXMLUnstructuredGridWriter> wr = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    CommonTools::AssociateProgressFunction(wr);
    std::cout<<"Output filename: "<<output_filename<<std::endl;
    wr->SetFileName(output_filename);
    wr->EncodeAppendedDataOff();
    wr->SetInputData(dst_reader->GetOutput());
    wr->Write();
}



