/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#include <vtkSmartPointer.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vector>
#include <vtkCellArray.h>
#include <vtkType.h>
#include <vtkCellData.h>

int main(int argc, char **argv)
{
    std::cout<<"Removes elements that are not tetrahedra"<<std::endl;
    std::cout<<"Usage: CheckTetraMesh inputmesh.vtk outputmessh.vtk"<<std::endl;
    
    if(argc<=1) exit(-1);
    
    vtkSmartPointer<vtkDataSetReader> rdr = vtkSmartPointer<vtkDataSetReader>::New();
    rdr->SetFileName(argv[1]);
    rdr->Update();
    
    vtkUnstructuredGrid *mesh = (vtkUnstructuredGrid *)rdr->GetOutput();
    
    vtkSmartPointer<vtkCellArray> array =vtkSmartPointer<vtkCellArray>::New();
    
    //verify that all the cells are tetrahedra
    //and create an indicator array of point ids belonging to cells
    //std::vector<int> indic(mesh->GetNumberOfPoints,0);
    
    for(int i=0; i<mesh->GetNumberOfCells(); i++)
    {
        if(mesh->GetCell(i)->GetNumberOfPoints()!=4)
        {
            std::cout<<"Cell "<<i<<" is not tetrahedron. Number of vertices: "<<mesh->GetCell(i)->GetNumberOfPoints()<<std::endl;
        }
        else
            array->InsertNextCell(mesh->GetCell(i));
    }

    
    mesh->SetCells(VTK_TETRA,array);
    mesh->GetCellData()->Initialize();
    
    std::cout<<"Checking tetra orientation. For (ABCD) tetra it is (AB x AC).AD:"<<std::endl;
    //test tetra orientation
    //formula for (ABCD) tetra is (AB x AC).AD
    unsigned long dp_positive=0; 
    unsigned long dp_negative=0; 
    
    for(int i=0; i<mesh->GetNumberOfCells(); i++){
        double A[3], B[3], C[3], D[3];
        mesh->GetCell(i)->GetPoints()->GetPoint(0,A);
        mesh->GetCell(i)->GetPoints()->GetPoint(1,B);
        mesh->GetCell(i)->GetPoints()->GetPoint(2,C);
        mesh->GetCell(i)->GetPoints()->GetPoint(3,D);
        
        
        double AB[3], AC[3], AD[3];
        for(int j=0; j<3; j++) 
        {
            AB[j]=B[j]-A[j];
            AC[j]=C[j]-A[j];
            AD[j]=D[j]-A[j];
        }
        
        //cross product
        double ABxAC[3]; 
        ABxAC[0] = AB[1]*AC[2] - AB[2]*AC[1];
        ABxAC[1] = AB[2]*AC[0] - AB[0]*AC[2];
        ABxAC[2] = AB[0]*AC[1] - AB[1]*AC[0];

        //dot product
        const double dp = ABxAC[0]*AD[0]+ABxAC[1]*AD[1]+ABxAC[2]*AD[2];
        if(dp>=0) dp_positive++;
        else dp_negative++;
    }
    
    std::cout<<"Number positive dot products: "<<dp_positive<<std::endl;
    std::cout<<"Number negative dot products: "<<dp_negative<<std::endl;
    
    
    
    vtkSmartPointer<vtkDataSetWriter> wr = vtkSmartPointer<vtkDataSetWriter>::New();
    wr->SetFileName(argv[2]);
    wr->SetInputData(mesh);
    wr->Write();

    
    
    return 0;
}
