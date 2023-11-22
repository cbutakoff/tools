/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
// Passes cell arrays from one unstructured grid to another , 
// given that the two are the same up to element order.

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
#include <vector>
#include <string>
#include <numeric>

#include <VTKCommonTools.h>
#include <vtkCallbackCommand.h>


std::vector<vtkIdType> argsort_cells(vtkDataSet* mesh);
std::vector<vtkIdType> argsort(std::vector<vtkIdType> v);
void pass_array(std::vector<vtkIdType> order_src, std::vector<vtkIdType> order_dst, 
              std::string array_name, vtkDataSet* mesh_src, vtkDataSet* mesh_dst );

template<class T>
void print_vector(std::vector<T> v)
{
   for( const auto& x: v) std::cout<<x<<",";
}


int main(int argc, char** argv)
{
    if(argc<4)
    {
        std::cout<<"Usage: SampleUG2UG meshsrc.vtu meshtgt.vtu out.vtu arrayname1 arrayname2 ..."<<std::endl;        
        std::cout<<"Try to use as recent VTK as possible, this uses vtkStaticCellLocator and it used to have bugs. Since VTK 9.3 seems to be ok, and vtkCellLocator is excruciatingly slow since vtk 9.0."<<std::endl;        
        exit(-1);
    }

    vtkOutputWindow::GetInstance()->SetDisplayModeToAlwaysStdErr();
    
    int c=1;
    const char* mesh_src_filename = argv[c++];
    const char* mesh_dst_filename = argv[c++];
    const char* output_filename   = argv[c++];
    int narrays = argc-c;
    if (narrays<0) {
        std::cout<<"No arrays specified"<<std::endl;        
        exit(-1);
    }

    std::vector<std::string> arraynames;
    while(c<argc) arraynames.push_back(argv[c++]);



    std::cout<<"Src:   "<<mesh_src_filename<<std::endl;
    std::cout<<"Dst:   "<<mesh_dst_filename<<std::endl;
    std::cout<<"Out:   "<<output_filename<<std::endl;
    std::cout<<"Arrays: "<<std::endl;
    for (const auto& name: arraynames) std::cout<<" - "<<name<<std::endl;
    
    
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


    //Identify max number of nodes in a cell
    const int maxcellsize_src = src_reader->GetOutput()->GetMaxCellSize();
    const int maxcellsize_dst = dst_reader->GetOutput()->GetMaxCellSize();
    if(maxcellsize_src!=maxcellsize_dst) {
        std::cout<<"Meshes have different element types"<<std::endl;        
        exit(-1);
    }
    if(src_reader->GetOutput()->GetNumberOfCells() != dst_reader->GetOutput()->GetNumberOfCells()){
        std::cout<<"Meshes have different number of elements"<<std::endl;        
        exit(-1);
    }
    if(src_reader->GetOutput()->GetNumberOfPoints() != dst_reader->GetOutput()->GetNumberOfPoints()){
        std::cout<<"Meshes have different number of nodes"<<std::endl;        
        exit(-1);
    }

    //Sort cells lexicographically
    std::cout<<"Identifying connectivity correspondence"<<std::endl;
    std::cout<<"Source mesh"<<std::endl;
    std::vector<vtkIdType> order_src = argsort_cells(src_reader->GetOutput());
    std::cout<<"Destination mesh"<<std::endl;
    std::vector<vtkIdType> order_dst = argsort_cells(dst_reader->GetOutput());

    for( const auto& arrayname : arraynames) {
      std::cout<<"Passing array "<<arrayname<<std::endl;
      pass_array( order_src, order_dst, arrayname, src_reader->GetOutput(), dst_reader->GetOutput() );
    }


    vtkSmartPointer<vtkXMLUnstructuredGridWriter> wr = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    CommonTools::AssociateProgressFunction(wr);
    std::cout<<"Output filename: "<<output_filename<<std::endl;
    wr->SetFileName(output_filename);
    wr->EncodeAppendedDataOff();
    wr->SetInputData(dst_reader->GetOutput());
    wr->Write();
}


void pass_array(std::vector<vtkIdType> order_src, 
              std::vector<vtkIdType> order_dst, 
              std::string array_name,
              vtkDataSet* mesh_src,
              vtkDataSet* mesh_dst
              ) {
   auto array_src = mesh_src->GetCellData()->GetArray(array_name.c_str());
   mesh_dst->GetCellData()->AddArray(array_src);

   auto array_dst = mesh_dst->GetCellData()->GetArray(array_name.c_str());
   auto order_dst_indices = argsort(order_dst);
   
   std::cout<<"Source indices"<<std::endl;
   print_vector(order_src);
   std::cout<<std::endl;
   std::cout<<"DST indices"<<std::endl;
   print_vector(order_dst);
   std::cout<<std::endl;
   std::cout<<"order_dst_indices"<<std::endl;
   print_vector(order_dst_indices);
   std::cout<<std::endl;


   for (const auto& idx_src: order_src){
      auto idx_dst = order_dst[order_dst_indices[idx_src]];
      array_dst->SetTuple1(idx_dst, array_src->GetTuple1(idx_src));
      std::cout<<idx_src<<" -> "<<idx_dst<<std::endl;
   }  
}

std::vector<vtkIdType> 
argsort(std::vector<vtkIdType> v) {
   // assumes v is an array of inides, i.e. it';s size is N and it has values from 0 to N-1
   // sort is in ascending order
   std::vector<vtkIdType> indices(v.size());
   std::iota(indices.begin(),indices.end(),0); //Fill with values from 0 to N-1
   std::sort(indices.begin(),indices.end(), [&](vtkIdType i, vtkIdType j){return v[i]>v[j];} );

   return indices;
}


std::vector<vtkIdType> 
argsort_cells(vtkDataSet* mesh) {
   using row=std::vector<vtkIdType>;
   std::vector<row> connectivity;

   vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
   for( vtkIdType cellId=0; cellId<mesh->GetNumberOfCells(); cellId++ ){
      mesh->GetCellPoints(cellId, ptIds);
      row r;
      for(vtkIdType ptid=0; ptid<ptIds->GetNumberOfIds(); ptid++){
         r.push_back(ptIds->GetId(ptid));
      }
      connectivity.push_back(r);

      std::cout<<cellId<<" : ";
      print_vector(r);
      std::cout<<std::endl;

   }

   std::vector<vtkIdType> indices(mesh->GetNumberOfCells());
   std::iota(indices.begin(),indices.end(),0); //Fill with values from 0 to N-1
   std::sort(indices.begin(),indices.end(), [&](vtkIdType i, vtkIdType j){return connectivity[i]<connectivity[j];} );

   for( const auto& v:indices ){
      std::cout<<v<<" : ";
      for( const auto& x: connectivity[v]) std::cout<<x<<",";
      std::cout<<std::endl;
   }


   return indices;
}


