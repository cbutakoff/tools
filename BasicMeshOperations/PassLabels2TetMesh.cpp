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

#include <vtkCallbackCommand.h>
#include <VTKCommonTools.h>
#include <vtkType.h>

    
typedef struct __bsc_entry {
    int nvertices;
    vtkIdType id;
    vtkIdType pt1;
    vtkIdType pt2;
    vtkIdType pt3;
    vtkIdType pt4;
    vtkIdType tet_id;
} BscEntry;


//#define LINEBREAK "\x0D\x0A"
//#define LINEBREAK "\x0A"

    
int main(int argc, char** argv)
{
    std::cout<<"exec volmesh.vtk surfmesh.vtk out_labels.txt arrayname volmesh_out.vtk scale correct_orientation(0|1)"<<std::endl;
    std::cout<<"Scale - rescale mesh by this factor, !!! for now works only for ALYA, for everything else put 1"<<std::endl;
    std::cout<<"correct_orientation - 1 or 0 whether to correct cell orientation or not"<<std::endl;
    std::cout<<"Labels must be celldata"<<std::endl;
    std::cout<<"Add extension .bsc to volmesh_out to get mesh in Alya format"<<std::endl;
    std::cout<<"Skip extension in volmesh_out to get mesh in Elmer format"<<std::endl;
    std::cout<<"out_labels.txt -- faces with labels, format: label point1 point2 point3 tetra_id"<<std::endl;
     
    
    
    if(argc<4) return -1;
    
    const char* volmesh_file = argv[1];
    const char* surfmesh_file = argv[2];
    const char* outfile = argv[3];
    const char* array_name = argv[4];
    const char* volmeshout = argv[5];
    
    float scale = atof(argv[6]);
    bool correct_orientation = atoi(argv[7])==1;
    
    
    std::cout<<"Volumetric mesh: "<<volmesh_file<<std::endl;
    std::cout<<"Surface mesh: "<<surfmesh_file<<std::endl;
    std::cout<<"Labels: "<<outfile<<std::endl;
    std::cout<<"Label array: "<<array_name<<std::endl;
    std::cout<<"Output mesh: "<<volmeshout<<std::endl;
    std::cout<<"Scale: "<<scale<<std::endl;
    if(correct_orientation)
        std::cout<<"Correcting orientation: ON"<<std::endl;
    else
        std::cout<<"Correcting orientation: OFF"<<std::endl;
    
    
    bool vtkoutput_mesh = false;
    bool alyaoutput_mesh = false;
        
    if(volmeshout!=NULL)
    {
        const char *ext = volmeshout+strlen(volmeshout)-4;
        if(strcmp(ext,".vtk")==0)  //if we are requesting .node .element volumetric mesh    
            vtkoutput_mesh = true;
        else if(strcmp(ext,".bsc")==0)
            alyaoutput_mesh = true;
    }

    vtkSmartPointer<vtkDataSetReader> vol_rdr = vtkSmartPointer<vtkDataSetReader>::New();
    CommonTools::AssociateProgressFunction(vol_rdr);

    vol_rdr->SetFileName(volmesh_file);
    vol_rdr->Update();
    vtkUnstructuredGrid* volmesh = (vtkUnstructuredGrid*)vol_rdr->GetOutput();
    
    
    if(volmesh->GetCell(0)->GetNumberOfPoints()<4)
    {
        std::cout<<"Supplied volumetric mesh has cells with : "<<volmesh->GetCell(0)->GetNumberOfPoints()<<" vertices"<<std::endl;
        std::cout<<"Make sure your mesh is actually volumetric"<<std::endl;
        exit(-1);
    }

    
    
    vtkSmartPointer<vtkPolyDataReader> poly_rdr = vtkSmartPointer<vtkPolyDataReader>::New();
    CommonTools::AssociateProgressFunction(poly_rdr);

    poly_rdr->SetFileName(surfmesh_file);
    poly_rdr->Update();
    vtkPolyData* surfmesh = poly_rdr->GetOutput();
    
    //read cell scalars from the surface mesh (labels)
    auto scalars = surfmesh->GetCellData()->GetArray(array_name);
    
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

        auto pscalars = c2pdata->GetOutput()->GetPointData()->GetArray(array_name);
        
        //to store ids with the tetra mesh points
        vtkSmartPointer<vtkShortArray> volmesh_regions = vtkSmartPointer<vtkShortArray>::New();
        volmesh_regions->SetName(array_name);
        volmesh_regions->SetNumberOfComponents(1);
        volmesh_regions->SetNumberOfValues(volmesh->GetNumberOfPoints());

        //fill the array of scalars to store with volumetric mesh
        for(vtkIdType i=0; i<volmesh->GetNumberOfPoints();i++)
        {
            volmesh_regions->SetTuple1(i,0);
        }

        //for every cell of the surface mesh
        for(vtkIdType i=0; i<surfmesh->GetNumberOfPoints(); i++)
        {
            vtkIdType ptid = ptloc->FindClosestPoint(surfmesh->GetPoint(i));
            volmesh_regions->SetTuple1(ptid, pscalars->GetTuple1(i));
        }
        
        volmesh->GetPointData()->AddArray(volmesh_regions);

        //:~ create PointData array in the volumetric mesh with labels
    }
    
    std::vector<BscEntry> labeldata;
    
    //create std::vector with information about the cell and point ids for every label
    //for every cell of the surface mesh
    std::cout<<"Looking for boundaries"<<std::endl;
    for(vtkIdType i=0; i<surfmesh->GetNumberOfCells(); i++)
    {
        if( i%10000 == 0 )
            std::cout<<"Cell "<<i<<"/"<<surfmesh->GetNumberOfCells()<<"\r"<<std::flush;
        
        vtkCell* cell = surfmesh->GetCell(i);
        BscEntry entry;
        
        //if( cell->GetNumberOfPoints()!=3 )
        //{
        //    std::cout<<"Face does not have 3 vertices. Id= "<<i<<std::endl;
        //    continue;
        //}

        entry.id = scalars->GetTuple1(i);
        entry.pt1 = ptloc->FindClosestPoint(cell->GetPoints()->GetPoint(0));
        entry.pt2 = ptloc->FindClosestPoint(cell->GetPoints()->GetPoint(1));
        entry.pt3 = ptloc->FindClosestPoint(cell->GetPoints()->GetPoint(2));
        if( cell->GetNumberOfPoints()==4 )
        {
            entry.pt4 = ptloc->FindClosestPoint(cell->GetPoints()->GetPoint(3));
            entry.nvertices=4;
        }   
        else if(cell->GetNumberOfPoints()==3)
        {
            entry.nvertices=3;
        }
        else 
        {
            cout<<"cell "<<i<<" has "<<cell->GetNumberOfPoints()<<" vertices, it is unsupported"<<endl;
            exit(-1);
        }
        


        double pt1[3];
        double pt2[3];
        double pt3[3];
        double pt4[3];
        volmesh->GetPoint(entry.pt1, pt1);
        volmesh->GetPoint(entry.pt2, pt2);
        volmesh->GetPoint(entry.pt3, pt3);
    
        if( entry.nvertices==4 )
        {
            volmesh->GetPoint(entry.pt4, pt4);            
        }
        else
        {
            pt4[0]=0;
            pt4[1]=0;
            pt4[2]=0;
        }

        double c[3];
        c[0] = (pt1[0]+pt2[0]+pt3[0]+pt4[0])/entry.nvertices;
        c[1] = (pt1[1]+pt2[1]+pt3[1]+pt4[1])/entry.nvertices;
        c[2] = (pt1[2]+pt2[2]+pt3[2]+pt4[2])/entry.nvertices;
        
        vtkIdType cellid = cellloc->FindCell(c);
        
        if(cellid==-1)
        {
            std::cout<<"Cell not found. Look for an error. Closest point:"<<c[0]<<" "<<c[1]<<" "<<c[2]<<" "<<std::endl;
        }
        entry.tet_id = cellid;
        
        labeldata.push_back(entry);
    }
    std::cout<<std::endl;
    
    
    //to store labels for cells
    vtkSmartPointer<vtkShortArray> volmesh_regions_cells = vtkSmartPointer<vtkShortArray>::New();
    volmesh_regions_cells->SetName(array_name);
    volmesh_regions_cells->SetNumberOfComponents(1);
    volmesh_regions_cells->SetNumberOfValues(volmesh->GetNumberOfCells());
    
    for(vtkIdType i1=0; i1<volmesh->GetNumberOfCells(); i1++)
    {
        volmesh_regions_cells->SetTuple1(i1,0); //reset
    }
    
    //store the labels
    std::ofstream file(outfile);
    for(vtkIdType i=0; i<labeldata.size(); i++)
    {
        BscEntry entry = labeldata[i];
        file<<entry.id<<" "<<
                entry.pt1+1<<" "<<
                entry.pt2+1<<" "<<
                entry.pt3+1<<" ";

        if( entry.nvertices==4 )
            file<<entry.pt4+1<<" ";

        file<<entry.tet_id+1<<LINEBREAK;
        
        if(volmeshout!=NULL)
        {
            volmesh_regions_cells->SetTuple1(entry.tet_id, entry.id);
        }
    
    }

    //store the labels 1
    std::ofstream file1((std::string(outfile)+".1").c_str());
    for(vtkIdType i=0; i<labeldata.size(); i++)
    {
        BscEntry entry = labeldata[i];
        file1<<i+1<<" "<<entry.id<<endl;        
    }


    
    if(volmeshout!=NULL)
    {
        volmesh->GetCellData()->AddArray(volmesh_regions_cells);
         
        
        if(vtkoutput_mesh)
        {   
            vtkSmartPointer<vtkDataSetWriter> wrwr = vtkSmartPointer<vtkDataSetWriter>::New();
            CommonTools::AssociateProgressFunction(wrwr);
            wrwr->SetFileTypeToBinary();
            wrwr->SetFileName(volmeshout);
            wrwr->SetInputData(volmesh);
            wrwr->Update();
        }
        else if(alyaoutput_mesh)
        {
            std::cout<<"Writing alya format"<<std::endl;
            CommonTools::SaveVolMeshBSC(volmesh, volmeshout, scale, correct_orientation);      
            
            //add the boundary
            char boundary_file[256];
            char boundary_elemtype_file[256];
            sprintf(boundary_file, "%s.bound", volmeshout);
            sprintf(boundary_elemtype_file, "%s.bound_type", volmeshout);
            std::ofstream file(boundary_file);
            std::ofstream file_elemtype(boundary_elemtype_file);

            for(vtkIdType i=0; i<labeldata.size(); i++)
            {
                BscEntry entry = labeldata[i];
                file<<i+1<<" "<<
                        entry.pt1+1<<" "<<
                        entry.pt2+1<<" "<<
                        entry.pt3+1<<" ";

                if( entry.nvertices==4 )
                    file<<entry.pt4+1;
        

                file<<LINEBREAK;

                file_elemtype<<i+1<<" "<<entry.nvertices<<LINEBREAK;

            }            
            


            

        }
        else //Write elmer mesh
        {
            std::cout<<"Writing elmer file"<<std::endl;
            std::string header_filename = "mesh.header";
            std::ofstream header_file(header_filename.c_str());
            header_file<<volmesh->GetNumberOfPoints()<<" "<<volmesh->GetNumberOfCells()<<" "<<labeldata.size()<<std::endl;
            header_file<<"2"<<std::endl;
            header_file<<"303 "<<labeldata.size()<<std::endl;
            header_file<<"504 "<<volmesh->GetNumberOfCells()<<std::endl;

            std::string node_filename = "mesh.nodes";
            std::ofstream node_file(node_filename.c_str());
            for(vtkIdType i=0; i<volmesh->GetNumberOfPoints(); i++)
            {
                double* pt = volmesh->GetPoint(i);
                node_file<<i+1<<" -1 "<<pt[0]<<" "<<pt[1]<<" "<<pt[2]<<std::endl;
            }

            std::string ele_filename = "mesh.elements";
            std::ofstream ele_file(ele_filename.c_str());
            for(vtkIdType i=0; i<volmesh->GetNumberOfCells(); i++)
            {
                vtkCell* cell = volmesh->GetCell(i);
                ele_file<<i+1<<" 1 504";
                for(int k=0; k<cell->GetNumberOfPoints(); k++)                
                    ele_file<<" "<<cell->GetPointId(k)+1;
                ele_file<<std::endl;
            }
            
            std::string bound_filename =  "mesh.boundary";
            std::ofstream bound_file(bound_filename.c_str());
            for(vtkIdType i=0; i<labeldata.size(); i++)
            {
                BscEntry entry = labeldata[i];
                bound_file<<i+1<<" "<<
                        entry.id<<" "<<
                        entry.tet_id+1<<" 0 303 "<<
                        entry.pt1+1<<" "<<
                        entry.pt2+1<<" "<<
                        entry.pt3+1;
                if( entry.nvertices==4 )
                    bound_file<<" "<<entry.pt4+1;

                bound_file<<std::endl;
            }
        }
    }

    return 0;
}
