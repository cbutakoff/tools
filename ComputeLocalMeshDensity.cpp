/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#include "CommonTools.h"
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>
#include <vtkFloatArray.h>
#include <vtkTriangle.h>
#include <vtkCellData.h>
#include <vtkType.h>

#include <set>

void usage(char *exe) {
    std::cout << "Compute local density of the mesh based on areas." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "-i <mesh.vtk> \t - input mesh" << std::endl;
    std::cout << "-o <mesh.vtk> \t - output mesh" << std::endl;

    exit(0);
}

int main(int argc, char **argv) {
    std::cout << "Version 1.0" << std::endl;
    if (argc < 2) usage(argv[0]);

    char* inmesh_filename = NULL;
    char* outmesh_filename = NULL;

    for (int c = 1; c < argc; c++) {
        if (strcmp(argv[c], "-i") == 0) {
            inmesh_filename = argv[++c];
        } else if (strcmp(argv[c], "-o") == 0) {
            outmesh_filename = argv[++c];
        }
    }

    vtkSmartPointer<vtkPolyData> inmesh = vtkSmartPointer<vtkPolyData>::Take(
            CommonTools::LoadShapeFromFile(inmesh_filename));

    vtkSmartPointer<vtkFloatArray> areas = vtkSmartPointer<vtkFloatArray>::New();
    areas->SetName("Area");
    areas->SetNumberOfComponents(1);
    areas->SetNumberOfTuples(inmesh->GetNumberOfCells());

    for (int i = 0; i < inmesh->GetNumberOfCells(); i++) {
        vtkCell* cell = inmesh->GetCell(i);

        vtkTriangle* triangle = dynamic_cast<vtkTriangle*> (cell);
        double p0[3];
        double p1[3];
        double p2[3];
        triangle->GetPoints()->GetPoint(0, p0);
        triangle->GetPoints()->GetPoint(1, p1);
        triangle->GetPoints()->GetPoint(2, p2);

        double area = vtkTriangle::TriangleArea(p0, p1, p2);
        areas->SetValue(i, area);
    }

    vtkSmartPointer<vtkPolyData> outmesh = vtkSmartPointer<vtkPolyData>::New();
    outmesh -> DeepCopy(inmesh);
    outmesh -> GetCellData() -> AddArray(areas);

    //compute densities
    vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkFloatArray> densities = vtkSmartPointer<vtkFloatArray>::New();
    densities->SetName("Density");
    densities->SetNumberOfComponents(1);
    densities->SetNumberOfTuples(inmesh->GetNumberOfCells());

    
    for (int cellid = 0; cellid < outmesh->GetNumberOfCells(); cellid++) {
        outmesh->GetCellPoints(cellid, cellPointIds);
        std::set<vtkIdType> neighbors;


        for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
            vtkSmartPointer<vtkIdList> idList =
                    vtkSmartPointer<vtkIdList>::New();
            idList->InsertNextId(cellPointIds->GetId(i));

            //get the neighbors of the cell
            vtkSmartPointer<vtkIdList> neighborCellIds =
                    vtkSmartPointer<vtkIdList>::New();

            outmesh->GetCellNeighbors(cellid, idList, neighborCellIds);

            for (vtkIdType j = 0; j < neighborCellIds->GetNumberOfIds(); j++) {
                neighbors.insert(neighborCellIds->GetId(j));
            }

        }

        double maxarea = 0;
        std::set<vtkIdType>::iterator it = neighbors.begin();
        for(; it != neighbors.end(); it++)
        {
            double area = areas->GetValue(*it);
            if( maxarea<area ) maxarea = area;
        }
        densities->SetValue(cellid, areas->GetValue(cellid)/maxarea);
        
    }

    outmesh->GetCellData()->AddArray(densities);
    
    std::cout<<"new code"<<std::endl;
    
    CommonTools::SaveShapeToFile(outmesh, outmesh_filename);



}

