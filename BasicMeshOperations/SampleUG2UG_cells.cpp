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

#include <iostream>

#include <VTKCommonTools.h>
#include <vtkCallbackCommand.h>


int main(int argc, char** argv)
{
    if(argc<4)
    {
        std::cout<<"Usage: SampleUG2UG meshsrc.vtu meshtgt.vtu arrayname out.vtu"<<std::endl;        
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

    auto ug = vtkSmartPointer<vtkUnstructuredGrid>::New();
    ug->SetPoints( cc->GetOutput()->GetPoints() );

    vtkSmartPointer<vtkResampleWithDataSet> res = vtkSmartPointer<vtkResampleWithDataSet>::New();
    CommonTools::AssociateProgressFunction(res);
    std::cout<<"Running vtkResampleWithDataSet"<<std::endl;
    res->SetInputData( ug );
    res->SetSourceConnection( src_reader->GetOutputPort() );
    res->Update();

    auto dst_mesh = vtkDataSet::SafeDownCast(dst_reader->GetOutput());
    auto res_mesh = vtkDataSet::SafeDownCast(res->GetOutput());
    dst_mesh->GetCellData()->AddArray( res_mesh->GetPointData()->GetArray(arrayname) );

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> wr = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    CommonTools::AssociateProgressFunction(wr);
    std::cout<<"Output filename: "<<output_filename<<std::endl;
    wr->SetFileName(output_filename);
    wr->EncodeAppendedDataOff();
    wr->SetInputData(dst_mesh);
    wr->Write();
}



