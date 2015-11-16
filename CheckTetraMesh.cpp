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
    
    vtkSmartPointer<vtkDataSetWriter> wr = vtkSmartPointer<vtkDataSetWriter>::New();
    wr->SetFileName(argv[2]);
    wr->SetInputData(mesh);
    wr->Write();

    
    
    return 0;
}
