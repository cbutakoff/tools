/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
/* Reads VTK polydata, and applies union boolean operation to all disconnected 
 meshes in the same file
 */

#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <CommonTools.h>

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cout<<"Reads VTK polydata, and applies union boolean operation to all disconnected meshes in the same file. v1.0"<<std::endl;
        std::cout<<"Usage: CreateSurfaceMeshEnvelope mesh.vtk output.vtk"<<std::endl;
        std::cout<<"Note: also reads stl and ply"<<std::endl;
        return -1;
    }
    
    int c=1;
    char* inshape = argv[c++];
    char* outshape = argv[c++];
    
    std::cout << "The following parameters will be used:" << std::endl;
    std::cout << "Mesh:" << inshape << std::endl;
    std::cout << "Output mesh:" << outshape << std::endl;    
    std::cout << "Operation: union"<< std::endl;    
    
    vtkSmartPointer<vtkPolyData> whole_mesh =  vtkSmartPointer<vtkPolyData>::Take(
            CommonTools::LoadShapeFromFile(inshape)          );
    
    //separate mesh into subparts
    vtkSmartPointer<vtkPolyDataConnectivityFilter> conn = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    conn->SetInputData(whole_mesh);
    conn->Update();
    
    
    const int n_regions = conn->GetNumberOfExtractedRegions();
    
    std::cout<<"Number of subparts: "<<n_regions<<std::endl;
    
    if(n_regions<2)
    {
        std::cout<<"Too few subparts, nothing to do. Quitting"<<std::endl;
        return -1;
    }
    

    vtkSmartPointer<vtkPolyData> union_mesh = vtkSmartPointer<vtkPolyData>::New(); //the result
    
    //put the first region into the union_mesh
    std::cout<<"Running connectivity"<<std::endl;
    conn->AddSpecifiedRegion(0);
    conn->SetExtractionModeToSpecifiedRegions();
    conn->Update();
    

    union_mesh->DeepCopy(conn->GetOutput());
    conn->DeleteSpecifiedRegion(0);
    
    vtkSmartPointer<vtkBooleanOperationPolyDataFilter> booleanOperation =
        vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
    booleanOperation->SetOperationToUnion();        
    
    for(int i=1; i<n_regions; i++)
    {
        std::cout<<"Adding subpart "<<i<<std::endl;
        conn->AddSpecifiedRegion(i);
        conn->SetExtractionModeToSpecifiedRegions();
        conn->Update();
        
        

        booleanOperation->SetInputData( 0, conn->GetOutput() );
        booleanOperation->SetInputData( 1, union_mesh );
        booleanOperation->Update();
        
        union_mesh->DeepCopy(booleanOperation->GetOutput());        
        conn->DeleteSpecifiedRegion(i);
    }
    
    //clean duplicate vertices to join the mesh
    std::cout<<"Cleaning"<<std::endl;
    vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
    cleaner->SetInputData(union_mesh);
    cleaner->Update();
    
    
    std::cout<<"Saving"<<std::endl;
    CommonTools::SaveShapeToFile(cleaner->GetOutput(), outshape);

   
    
    return 0;
}
