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

//std
#include <set>
#include <iostream>
#include <vtkType.h>

typedef std::set<vtkIdType> CellNeighborsType;
void GetCellNeighbors(vtkDataSet* mesh, vtkIdType cellId, CellNeighborsType &neighbors);

int main(int argc, char** argv) {
    std::cout << "Usage: Reclassify inUGmesh.vtk outUGmesh.vtk <material_id>" << std::cout;
    if (argc < 4) return -1;

    int c = 1;
    const char* in_meshfilename = argv[c++];
    const char* out_meshfilename = argv[c++];
    const int material_id2replace = atoi(argv[c++]);

    std::cout << "Input mesh: " << in_meshfilename << std::endl;
    std::cout << "Output mesh: " << out_meshfilename << std::endl;
    std::cout << "Material id to reclassify: " << material_id2replace << std::endl;

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
    CellNeighborsType neighbors;
    std::vector<int> materials( material_range[1]+1, 0 ); //increment elements for every material found in the neighborhood
    
    for(vtkIdType i=0; i<mesh->GetNumberOfCells(); i++)
    {
        neighbors.clear();
        if(material_ids->GetTuple1(i)==material_id2replace) //find cell matching material id
        {
            //reset material array
            for(int i=0; i<materials.size(); i++) materials[i]=0;
            
            GetCellNeighbors(mesh, i, neighbors);
            for(CellNeighborsType::iterator it=neighbors.begin(); it!=neighbors.end(); it++)
            {
                const int material = material_ids->GetTuple1(*it);
                materials[material] = materials[material]+1;
            }
            
            materials[material_id2replace] = 0; //to avoid assigning the same material
            materials[1] = 0; //to avoid assigning the unclassified region
            
            if(neighbors.size()>0)
            {
                //find the most popular material
                int voted_material = 0; 
                int n_votes = materials[0];
                for(int material=1; material<materials.size(); material++)
                {
                    //std::cout<<materials[material]<<" ";
                    if(materials[material]>n_votes) 
                    {
                        voted_material = material;
                        n_votes = materials[material];
                    }
                }
                
                material_ids->SetTuple1(i, voted_material);
                std::cout<<"Correcting element "<<i<<" to "<<voted_material<<std::endl;
            }
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