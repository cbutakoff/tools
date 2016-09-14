/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>

int main(int argc, char* argv[])
{
    if (argc < 4) {
        std::cout<<"Boolean surface mesh operations. v1.0"<<std::endl;
        std::cout<<"Usage: SurfaceMeshBoolean mesh1.vtk mesh2.vtk operation output.vtk"<<std::endl;
        std::cout<<"Operation: union, difference, intersect"<<std::endl;
        return -1;
    }
    
    int c=1;
    char* inshape1 = argv[c++];
    char* inshape2 = argv[c++];
    char* operation = argv[c++];
    char* outshape = argv[c++];
    
    std::cout << "The following parameters will be used:" << std::endl;
    std::cout << "Mesh 1:" << inshape1 << std::endl;
    std::cout << "Mesh 2:" << inshape2 << std::endl;
    std::cout << "Output mesh:" << outshape << std::endl;    
    std::cout << "Operation:" << operation << std::endl;    
    
    
    vtkSmartPointer<vtkPolyDataReader> reader1 =
      vtkSmartPointer<vtkPolyDataReader>::New();
    reader1->SetFileName(inshape1);
    reader1->Update();
    vtkPolyData* input1 = reader1->GetOutput();
 
    vtkSmartPointer<vtkPolyDataReader> reader2 =
      vtkSmartPointer<vtkPolyDataReader>::New();
    reader2->SetFileName(inshape2);
    reader2->Update();
    vtkPolyData* input2 = reader2->GetOutput();
    
    vtkSmartPointer<vtkBooleanOperationPolyDataFilter> booleanOperation =
        vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
    
    
    booleanOperation->SetInputData( 0, input1 );
    booleanOperation->SetInputData( 1, input2 );
    
    if(strcmp(operation,"union")==0)
        booleanOperation->SetOperationToUnion();    
    else if(strcmp(operation,"difference")==0)
        booleanOperation->SetOperationToDifference();
    else if(strcmp(operation,"intersect")==0)
        booleanOperation->SetOperationToIntersection();
    else 
    {
        std::cout<<"Unknown operation"<<std::endl;
        return -1;
    }
    
    

    booleanOperation->Update();
    
    vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
    cleaner->SetInputData(booleanOperation->GetOutput());
    cleaner->Update();
    
    
    vtkSmartPointer<vtkPolyDataWriter> writer =
      vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(outshape);
    writer->SetInputData(cleaner->GetOutput());
    writer->Write();
    
    
    return 0;
}
