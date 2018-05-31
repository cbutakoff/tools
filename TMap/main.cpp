#include "TMap.h"

#include <stdlib.h>
#include <stdio.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

int main(int argc, char* argv[])
{
    const char* outputmesh_filename = "out.vtk";
    const char* inputmesh_filename = "sample.vtk";

    vtkSmartPointer<vtkPolyDataReader> pdr = vtkSmartPointer<vtkPolyDataReader>::New();
    pdr->SetFileName(inputmesh_filename);
    pdr->Update();

    
    
    
//    TMapUtils tt;
//    tt.SetInputMesh(pdr->GetOutput());
//    tt.AddLandmarkConstraint(1,1,0);
//    tt.AddLandmarkConstraint(2,0,1);
//    tt.ZeroBC();
//    tt.MeshFromBC();
//
//
//    vtkSmartPointer<vtkPolyDataWriter> wr = vtkSmartPointer<vtkPolyDataWriter>::New();
//    wr->SetFileName("ttt.vtk");
//    wr->SetInputData(tt.GetOutput());
//    wr->Write();    
//
//    exit(0);
    
    TMap T;
   
    try
    {
        T.SetInput( pdr->GetOutput() );
    
        T.AddFixedBoundaryLandmarkConstraints();
        T.AddLandmarkConstraint(293,0.5,0.5); //for sample.vtk
        T.AddLandmarkConstraint(550,-0.5,-0.5); //for sample.vtk
        T.AddLandmarkConstraint(250,-0.5,0.5); //for sample.vtk
        T.AddLandmarkConstraint(487,0.5,-0.5); //for sample.vtk


        T.SetMaxNumberOfIterations(500);
        T.Update();
        
        vtkSmartPointer<vtkPolyDataWriter> wr = vtkSmartPointer<vtkPolyDataWriter>::New();
        wr->SetFileName(outputmesh_filename);
        wr->SetInputData(T.GetOutput());
        wr->Write();        
    }
    catch( const std::logic_error& err )
    {
        std::cout<<err.what()<<std::endl;
    }
        
 
}
