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
#include <vtkPolyData.h>
#include <vtkShortArray.h>
#include <vtkCell.h>
#include <vtkPointLocator.h>
#include <vtkType.h>
#include <vtkPointData.h>
#include <vtkCellDataToPointData.h>
#include <vector>
#include <vtkCellLocator.h>
#include <vtkPolyDataReader.h>
#include <vtkCellData.h>
#include <vtkCell.h>
#include <fstream>

    
typedef struct __bsc_entry {
    int id;
    int pt1;
    int pt2;
    int pt3;
    int tet_id;
} BscEntry;


    
int main(int argc, char** argv)
{
    std::cout<<"exec volmesh.vtk surfmesh.vtk out_labels.txt arrayname volmesh_out.vtk"<<std::endl;
    std::cout<<"Labels must be celldata"<<std::endl;
    std::cout<<"Skip extension in volmesh_out to get mesh in Elmer format"<<std::endl;
    
    if(argc<4) return -1;
    
    const char* volmesh_file = argv[1];
    const char* surfmesh_file = argv[2];
    const char* outfile = argv[3];
    const char* array_name = argv[4];
    const char* volmeshout = argv[5];
    
    bool vtkoutput_mesh = false;
    
    if(volmeshout!=NULL)
    {
        const char *ext = volmeshout+strlen(volmeshout)-4;
        if(strcmp(ext,".vtk")==0)  //if we are requesting .node .element volumetric mesh    
            vtkoutput_mesh = true;
    }

    vtkSmartPointer<vtkDataSetReader> vol_rdr = vtkSmartPointer<vtkDataSetReader>::New();
    vol_rdr->SetFileName(volmesh_file);
    vol_rdr->Update();
    vtkUnstructuredGrid* volmesh = (vtkUnstructuredGrid*)vol_rdr->GetOutput();
    
    vtkSmartPointer<vtkPolyDataReader> poly_rdr = vtkSmartPointer<vtkPolyDataReader>::New();
    poly_rdr->SetFileName(surfmesh_file);
    poly_rdr->Update();
    vtkPolyData* surfmesh = poly_rdr->GetOutput();
    
    //read cell scalars from the surface mesh (labels)
    vtkShortArray* scalars = (vtkShortArray*)surfmesh->GetCellData()->GetArray(array_name);
    
    vtkSmartPointer<vtkPointLocator> ptloc = vtkSmartPointer<vtkPointLocator>::New();
    ptloc->SetDataSet(volmesh);
    ptloc->BuildLocator();
    
    vtkSmartPointer<vtkCellLocator> cellloc = vtkSmartPointer<vtkCellLocator>::New();
    cellloc->SetDataSet(volmesh);
    cellloc->BuildLocator();
    
    if(volmeshout!=NULL)
    {
        //create PointData array in the volumetric mesh with labels
        
        //pass the cell array with labels to point array
        vtkSmartPointer<vtkCellDataToPointData> c2pdata = vtkSmartPointer<vtkCellDataToPointData>::New();
        c2pdata->SetInputData(surfmesh);
        c2pdata->PassCellDataOn();
        c2pdata->Update();

        vtkShortArray* pscalars = (vtkShortArray*)c2pdata->GetOutput()->GetPointData()->GetArray(array_name);
        
        //to store ids with the tetra mesh points
        vtkSmartPointer<vtkShortArray> volmesh_regions = vtkSmartPointer<vtkShortArray>::New();
        volmesh_regions->SetName(array_name);
        volmesh_regions->SetNumberOfComponents(1);
        volmesh_regions->SetNumberOfValues(volmesh->GetNumberOfPoints());

        //fill the array of scalars to store with volumetric mesh
        for(int i=0; i<volmesh->GetNumberOfPoints();i++)
        {
            volmesh_regions->SetValue(i,0);
        }

        //for every cell of the surface mesh
        for(int i=0; i<surfmesh->GetNumberOfPoints(); i++)
        {
            vtkIdType ptid = ptloc->FindClosestPoint(surfmesh->GetPoint(i));
            volmesh_regions->SetValue(ptid, pscalars->GetValue(i));
        }
        
        volmesh->GetPointData()->AddArray(volmesh_regions);

        //:~ create PointData array in the volumetric mesh with labels
    }
    
    std::vector<BscEntry> labeldata;
    
    //create std::vector with information about the cell and point ids for every label
    //for every cell of the surface mesh
    for(int i=0; i<surfmesh->GetNumberOfCells(); i++)
    {
        vtkCell* cell = surfmesh->GetCell(i);
        BscEntry entry;
        
        if( cell->GetNumberOfPoints()!=3 )
        {
            std::cout<<"Face does not have 3 vertices. Id= "<<i<<std::endl;
            continue;
        }

        entry.id = scalars->GetValue(i);
        entry.pt1 = ptloc->FindClosestPoint(cell->GetPoints()->GetPoint(0));
        entry.pt2 = ptloc->FindClosestPoint(cell->GetPoints()->GetPoint(1));
        entry.pt3 = ptloc->FindClosestPoint(cell->GetPoints()->GetPoint(2));
        
        double pt1[3];
        volmesh->GetPoint(entry.pt1, pt1);
        double pt2[3];
        volmesh->GetPoint(entry.pt2, pt2);
        double pt3[3];
        volmesh->GetPoint(entry.pt3, pt3);
        double c[3];
        c[0] = (pt1[0]+pt2[0]+pt3[0])/3;
        c[1] = (pt1[1]+pt2[1]+pt3[1])/3;
        c[2] = (pt1[2]+pt2[2]+pt3[2])/3;
        
        vtkIdType cellid = cellloc->FindCell(c);
        
        if(cellid==-1)
        {
            std::cout<<"Cell not found. Look for an error. Closest point:"<<c[0]<<" "<<c[1]<<" "<<c[2]<<" "<<std::endl;
        }
        entry.tet_id = cellid;
        
        labeldata.push_back(entry);
    }

    //to store labels for cells
    vtkSmartPointer<vtkShortArray> volmesh_regions_cells = vtkSmartPointer<vtkShortArray>::New();
    volmesh_regions_cells->SetName(array_name);
    volmesh_regions_cells->SetNumberOfComponents(1);
    volmesh_regions_cells->SetNumberOfValues(volmesh->GetNumberOfCells());
    
    for(int i1=0; i1<volmesh->GetNumberOfCells(); i1++)
    {
        volmesh_regions_cells->SetValue(i1,0); //reset
    }
    
    //store the labels
    std::ofstream file(outfile);
    for(int i=0; i<labeldata.size(); i++)
    {
        BscEntry entry = labeldata[i];
        file<<entry.id<<" "<<
                entry.pt1+1<<" "<<
                entry.pt2+1<<" "<<
                entry.pt3+1<<" "<<
                entry.tet_id+1<<"\x0D\x0A";
        
        if(volmeshout!=NULL)
        {
            volmesh_regions_cells->SetValue(entry.tet_id, entry.id);
        }
    
    }
    
    if(volmeshout!=NULL)
    {
        volmesh->GetCellData()->AddArray(volmesh_regions_cells);
         
        if(vtkoutput_mesh)
        {
            vtkSmartPointer<vtkDataSetWriter> wrwr = vtkSmartPointer<vtkDataSetWriter>::New();
            wrwr->SetFileTypeToBinary();
            wrwr->SetFileName(volmeshout);
            wrwr->SetInputData(volmesh);
            wrwr->Update();
        }
        else //Write elmer mesh
        {
            std::string header_filename = "mesh.header";
            std::ofstream header_file(header_filename.c_str());
            header_file<<volmesh->GetNumberOfPoints()<<" "<<volmesh->GetNumberOfCells()<<" "<<labeldata.size()<<std::endl;
            header_file<<"2"<<std::endl;
            header_file<<"303 "<<labeldata.size()<<std::endl;
            header_file<<"504 "<<volmesh->GetNumberOfCells()<<std::endl;

            std::string node_filename = "mesh.nodes";
            std::ofstream node_file(node_filename.c_str());
            for(int i=0; i<volmesh->GetNumberOfPoints(); i++)
            {
                double* pt = volmesh->GetPoint(i);
                node_file<<i+1<<" -1 "<<pt[0]<<" "<<pt[1]<<" "<<pt[2]<<std::endl;
            }

            std::string ele_filename = "mesh.elements";
            std::ofstream ele_file(ele_filename.c_str());
            for(int i=0; i<volmesh->GetNumberOfCells(); i++)
            {
                vtkCell* cell = volmesh->GetCell(i);
                ele_file<<i+1<<" 1 504 "<<cell->GetPointId(0)+1<<" "<<
                        cell->GetPointId(1)+1<<" "<<
                        cell->GetPointId(2)+1<<" "<<
                        cell->GetPointId(3)+1<<std::endl;
            }
            
            std::string bound_filename =  "mesh.boundary";
            std::ofstream bound_file(bound_filename.c_str());
            for(int i=0; i<labeldata.size(); i++)
            {
                BscEntry entry = labeldata[i];
                bound_file<<i+1<<" "<<
                        entry.id<<" "<<
                        entry.tet_id+1<<" 0 303 "<<
                        entry.pt1+1<<" "<<
                        entry.pt2+1<<" "<<
                        entry.pt3+1<<std::endl;
            }
        }
    }

    return 0;
}
