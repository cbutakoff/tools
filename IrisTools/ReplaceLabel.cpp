//reclassifies tiny regions and region 0
//uses majority voting among the neighboring cells
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkUnstructuredGrid.h"
#include "vtkShortArray.h"
#include "vtkDataArray.h"
#include "vtkCellData.h"
#include "vtkDataSet.h"
#include "vtkIdList.h"
#include "vtkUnstructuredGridWriter.h"
#include <vtkType.h>

//std
#include <set>
#include <iostream>

typedef std::set<vtkIdType> CellNeighborsType;
void GetCellNeighbors(vtkDataSet* mesh, vtkIdType cellId, CellNeighborsType &neighbors);

int main(int argc, char** argv) {
    std::cout << "Usage: ReplaceLabel inUGmesh.vtk outUGmesh.vtk <label2find> <label2substitute>" << std::cout;
    if (argc < 5) return -1;

    int c = 1;
    const char* in_meshfilename = argv[c++];
    const char* out_meshfilename = argv[c++];
    const int material_id2replace = atoi(argv[c++]);
    const int material_id2new = atoi(argv[c++]);

    std::cout << "Input mesh: " << in_meshfilename << std::endl;
    std::cout << "Output mesh: " << out_meshfilename << std::endl;
    std::cout << "Material id to replace: " << material_id2replace << std::endl;
    std::cout << "Material new label: " << material_id2new << std::endl;

    std::cout << "Reading mesh" << std::endl;
    vtkSmartPointer<vtkUnstructuredGridReader> rd = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    rd->SetFileName(in_meshfilename);
    rd->Update();

    std::cout << "Npoints :" << rd->GetOutput()->GetNumberOfPoints() << std::endl;
    std::cout << "Ncells :" << rd->GetOutput()->GetNumberOfCells() << std::endl;



    vtkUnstructuredGrid* mesh = rd->GetOutput();
    mesh->GetCellData()->SetActiveScalars("Material");

    vtkShortArray* material_ids = dynamic_cast<vtkShortArray*> (mesh->GetCellData()->GetArray("Material"));
    double material_range[2];
    material_ids->GetRange(material_range);
    std::cout << "Range of materials: " << material_range[0] << "-" << material_range[1] << std::endl;

    if (material_id2replace < material_range[0] || material_id2replace > material_range[1]) {
        std::cout << "Requested material out of range." << std::endl;
        return -1;
    }

    ///////////////////////////////////////////////////////////////////
    // 
    // main processing
    for(vtkIdType i=0; i<mesh->GetNumberOfCells(); i++)
    {
        if(material_ids->GetTuple1(i)==material_id2replace) //find cell matching material id
        {
            material_ids->SetTuple1(i, material_id2new);
        }
    }
    
    //saving
    vtkSmartPointer<vtkUnstructuredGridWriter> wr = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    wr->SetFileName(out_meshfilename);
    wr->SetInputData(mesh);
    wr->SetFileTypeToBinary();
    wr->Write();
        
    return 0;
}

void GetCellNeighbors(vtkDataSet* mesh, vtkIdType cellId, std::set<vtkIdType> &neighbors) {

    vtkSmartPointer<vtkIdList> cellPointIds =
            vtkSmartPointer<vtkIdList>::New();

    mesh->GetCellPoints(cellId, cellPointIds);
    
    vtkSmartPointer<vtkIdList> idList =
            vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> neighborCellIds =
            vtkSmartPointer<vtkIdList>::New();
    
    for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
        idList->Reset();
        neighborCellIds->Reset();
        idList->InsertNextId(cellPointIds->GetId(i));

        //get the neighbors of the cell
        mesh->GetCellNeighbors(cellId, idList, neighborCellIds);

        for (vtkIdType j = 0; j < neighborCellIds->GetNumberOfIds(); j++) {
            neighbors.insert(neighborCellIds->GetId(j));
        }

    }
}