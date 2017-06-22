#include "TMap.h"

#include <stdlib.h>
#include <stdio.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

int main(int argc, char* argv[])
{
    //const char* inputmesh_filename = "/home/costa/Dropbox/Ideas/T-Map/mediummesh.vtk";
    const char* outputmesh_filename = "/home/costa/Dropbox/Ideas/T-Map/out.vtk";
    //const char* inputmesh_filename = "/home/costa/Dropbox/Ideas/atrial_trunks/flat_1.vtk";
    const char* inputmesh_filename = "/home/costa/Dropbox/Ideas/T-Map/sample.vtk";
    //const char* inputmesh_filename = "/home/costa/Dropbox/Ideas/T-Map/smallmesh1.vtk";
    //const char* inputmesh_filename = "/home/costa/Data/code/CPP/bin/tools-private/bin/tri.vtk";

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


        //T.AddLandmarkConstraint(6,0.2,0.2); //for mediummesh.vtk
        //T.AddLandmarkConstraint(751,0.5,0.5); //for flat_1.vtk
        //T.AddLandmarkConstraint(7385,-0.5,-0.5); //for flat_1.vtk
        //T.AddLandmarkConstraint(5250,-0.5,0.5); //for flat_1.vtk
        //T.AddLandmarkConstraint(4070,0.5,-0.5); //for flat_1.vtk
        //T.AddLandmarkConstraint(761,0.05,-0.15); //for atrium.vtk
        //T.AddLandmarkConstraint(11825,-0.5,0.0); //for atrium.vtk
        //T.AddLandmarkConstraint(12345,-0.5,0.4); //for atrium.vtk
        //T.AddLandmarkConstraint(2608,0.5,0.25); //for atrium.vtk
        //T.AddLandmarkConstraint(3767,0.5,-0.25); //for atrium.vtk
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