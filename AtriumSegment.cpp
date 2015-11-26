/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file
    \brief Mesh thresholding based on some scalar field to extract the required number of regions. Searches for a threshold.
 */

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <string>
#include <iostream>
#include <set>
#include <limits>
#include <string>
#include <stack>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkShortArray.h>
#include <vtkType.h>
#include <vtkIdList.h>
#include <vtkPolyDataReader.h>
#include <vtkWriter.h>
#include <vtkPolyDataWriter.h>
#include <vtkThreshold.h>
#include <vtkConnectivityFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellLocator.h>
#include <vtkCell.h>
#include <vtkLongArray.h>
#include <vtkDoubleArray.h>

#include "MinHeap/MinHeap.h"

#define MARKER_SET "MarkerSet"



typedef struct __CellDescriptor
{
    unsigned long label;
    unsigned long cell_id;
} CellDescriptor;




//return all faces that are local minima. Arrayname is the array containing the height function
//the markers are stored in cell array in the mesh "MarkerSet"
void GenerateMarkerSetAsLocalMin(vtkPolyData* mesh, const char* arrayname);
int GenerateMarkerSetByIterativeThresholding(vtkPolyData* mesh, const char* arrayname, const int nregions);
int GenerateMarkerSetByIterativeThresholdingSteps(vtkPolyData* mesh, const char* arrayname, const int nregions);
void GetCellNeighbors(vtkPolyData* mesh, vtkIdType cellId, std::set<vtkIdType>& neighbors);

void FastMarchingGrowingFromBiggestTriangle(vtkPolyData* mesh, const char* arrayname);
void FastMarchingCumulativeFunction(vtkPolyData* mesh, const char* arrayname);

void FastMarchingWatershed(vtkPolyData* mesh, const char* arrayname);

void GrowHeap(vtkPolyData* mesh, unsigned long cellid, MinHeap<CellDescriptor, double>& heap);




int main(int argc, char *argv[]) {

    if (argc < 2) {
        std::cout << "Atrial mesh thresholding v1.0" << std::endl;
        std::cout << "Parameters: " << std::endl;
        std::cout << "-i input_mesh.vtk" << std::endl;
        std::cout << "-o output_mesh.vtk" << std::endl;
        std::cout << "-r number_of_regions \t -how many regions are to be extracted (default 5)" << std::endl;
        std::cout << "-a array_name \t\t -name of the cell array that is to be used for thresholding" << std::endl;
        return -1;
    }

    char* mesh_filename=NULL;
    char* output_filename=NULL;
    std::string array_name = "area_change"; //default
    int nregions = 5; //number of regions to expect (normally 5- 4 PVs + Appendage)
    

    for(int c=1; c<argc; c++)
    {
            if( strcmp(argv[c],"-i")==0 )
            {
                    mesh_filename = argv[++c];
            }
            else if( strcmp(argv[c],"-o")==0 )
            {
                    output_filename = argv[++c];
            }
            else if( strcmp(argv[c],"-a")==0 )
            {
                    array_name = argv[++c];
            }
            else if( strcmp(argv[c],"-r")==0 )
            {
                    nregions = atoi(argv[++c]);
            }
    }    


    if(mesh_filename==NULL || output_filename==NULL){
        std::cout<<"Input and output meshes are required"<<std::endl;
        return -1;
    }
    
    std::cout<<"Input mesh: "<<mesh_filename<<std::endl;
    std::cout<<"Output mesh: "<<output_filename<<std::endl;
    std::cout<<"Array for thresholding: "<<array_name<<std::endl;
    std::cout<<"Required number of regions: "<<nregions<<std::endl;
    
    vtkSmartPointer<vtkPolyDataReader> rd = vtkSmartPointer<vtkPolyDataReader>::New();
    rd->SetFileName(mesh_filename);
    rd->Update();
    vtkPolyData* mesh = rd->GetOutput();
    
    //GenerateMarkerSetByIterativeThresholding(mesh, array_name.c_str(), nregions);
    FastMarchingCumulativeFunction(mesh, array_name.c_str());
    //GenerateMarkerSetByIterativeThresholdingSteps(mesh, "GrowingOrder", nregions);
    //GenerateMarkerSetAsLocalMin(mesh, array_name.c_str());
    //FastMarchingWatershed(mesh, array_name.c_str());
    
    vtkSmartPointer<vtkPolyDataWriter> wr = vtkSmartPointer<vtkPolyDataWriter>::New();
    wr->SetFileName(output_filename);
    wr->SetInputData(mesh);
    wr->Write();
    
    
    return 0;
}

void GenerateMarkerSetAsLocalMin(vtkPolyData* mesh, const char* arrayname) {
    vtkDataArray* f = mesh->GetCellData()->GetArray(arrayname);
    
    if(f==NULL)
    {
        std::cout<<"Cell array "<<arrayname<<" was not found"<<std::endl;
        throw;
    }

    vtkSmartPointer<vtkShortArray> labels = vtkSmartPointer<vtkShortArray>::New();
    labels->SetName(MARKER_SET);
    labels->SetNumberOfComponents(1);
    labels->SetNumberOfTuples(f->GetNumberOfTuples());

    //for each cell, compare f with the neighborhood, and mark the cell if the value is 
    //strictly smaller
    
    const unsigned long ncells = mesh->GetNumberOfCells();
    
    std::cout<<"Building marker set: "<<std::endl;
    unsigned long label_id = 1;
    for (int i = 0; i < ncells; i++) {
        std::set<vtkIdType> neighbors;
        GetCellNeighbors(mesh, i, neighbors);

        if(i%100)
            std::cout<<i<<"/"<<ncells<<"\r";
        
        //std::cout << "Cell neighbor ids are: " << std::endl;
        
        double f_mean = 0;
        for(std::set<vtkIdType>::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
        {
            f_mean += f->GetTuple1(*it1);
            //std::cout << " " << *it1;
        }
        f_mean /= neighbors.size();
        //std::cout << std::endl;

        if(f->GetTuple1(i)<f_mean) labels->SetTuple1(i,label_id++);
        else labels->SetTuple1(i,0);       
    }
    std::cout<<std::endl;
    
    mesh->GetCellData()->AddArray(labels);
}



//from: http://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/CellPointNeighbors
//does not return the cell cellId
void GetCellNeighbors(vtkPolyData* mesh, vtkIdType cellId, std::set<vtkIdType>& neighbors) {
    // Find all cells connected to point 0
    mesh->BuildLinks();
    vtkSmartPointer<vtkIdList> cellPointIds =
            vtkSmartPointer<vtkIdList>::New();
    mesh->GetCellPoints(cellId, cellPointIds);

    // neighbor cells may be listed multiple times
    // use std::set instead of std::list if you want a unique list of neighbors
    //std::list<vtkIdType> neighbors;

    /*For each vertex of the cell, we calculate which cells uses that point.
     So if we make this, for each vertex, we have all the neighbors.
     In the case we use ''cellPointIds'' as a parameter of ''GeteCellNeighbors'',
     we will obtain an empty set. Because the only cell that is using that set of points
     is the current one. That is why we have to make each vertex at time.*/

    for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
        vtkSmartPointer<vtkIdList> idList =
                vtkSmartPointer<vtkIdList>::New();
        idList->InsertNextId(cellPointIds->GetId(i));

        //get the neighbors of the cell
        vtkSmartPointer<vtkIdList> neighborCellIds =
                vtkSmartPointer<vtkIdList>::New();

        mesh->GetCellNeighbors(cellId, idList, neighborCellIds);

        for (vtkIdType j = 0; j < neighborCellIds->GetNumberOfIds(); j++) {
            neighbors.insert(neighborCellIds->GetId(j));
        }
    }
}


//return 0 if ok, return -1 if could not find nregions
int GenerateMarkerSetByIterativeThresholding(vtkPolyData* mesh, const char* arrayname, const int nregions) {
    vtkDataArray* f = mesh->GetCellData()->GetArray(arrayname);
    
    if(f==NULL)
    {
        std::cout<<"Cell array "<<arrayname<<" was not found"<<std::endl;
        throw;
    }

    vtkSmartPointer<vtkShortArray> labels = vtkSmartPointer<vtkShortArray>::New();
    labels->SetName(MARKER_SET);
    labels->SetNumberOfComponents(1);
    labels->SetNumberOfTuples(f->GetNumberOfTuples());

    //for each cell, compare f with the neighborhood, and mark the cell if the value is 
    //strictly smaller
    
    const unsigned long ncells = mesh->GetNumberOfCells();
    
    //find min and max of the f
    double* range = f->GetRange(0);
    const double f_min = range[0];
    const double f_max = range[1];
    
    std::cout<<"The values in the function are in the range: ["<<f_min<<";"<<f_max<<"]"<<std::endl;


    mesh->GetCellData()->SetActiveScalars(arrayname);
    
    vtkSmartPointer<vtkThreshold> th = vtkSmartPointer<vtkThreshold>::New();
    th->SetInputData(mesh);
    //th->SetAttributeModeToUseCellData();
    th->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, arrayname);
    
    vtkSmartPointer<vtkConnectivityFilter> conn = vtkSmartPointer<vtkConnectivityFilter>::New();
    conn->SetExtractionModeToAllRegions();
    conn->ColorRegionsOn();    
    
    double T_old = f_max;
    double T = f_min+(f_max-f_min)/2; //thresholds
    
    bool regions_found = false;
    
    while( !regions_found && fabs(T-T_old)>1e-5 )
    {
        th->ThresholdByLower(T);
        th->Update();
        
        conn->SetInputData(th->GetOutput());
        conn->Update();    
        
        const unsigned long n_extr_regions = conn->GetNumberOfExtractedRegions();
        std::cout<<"T="<<T<<"; components="<<n_extr_regions<<std::endl;
        
        if(n_extr_regions == nregions) 
        {   
            regions_found = true;
        }
        else if(n_extr_regions > nregions)
        {
            T_old = T;
            T = f_min+(T-f_min)/2;
        }
        else 
        {
            double t_min = T>T_old?T_old:T;
            double t_max = T>T_old?T:T_old;
            T_old = T;
            T = t_min+(t_max-t_min)/2;
        }
    }

    if(regions_found)
    {
        //mark the cells based on the threshold filter (using locator)
        for (int i = 0; i < ncells; i++) labels->SetTuple1(i,0);

        vtkSmartPointer<vtkCellLocator> in_target_locator = vtkSmartPointer<vtkCellLocator>::New();
        in_target_locator->SetDataSet(mesh);
        in_target_locator->BuildLocator();

        for( int region=0; region<conn->GetNumberOfExtractedRegions(); region++)
        {
            const int label = region+1;


            vtkSmartPointer<vtkThreshold> th1 = vtkSmartPointer<vtkThreshold>::New();
            th1->SetInputData(conn->GetOutput());
            th1->ThresholdBetween(region,region);
            th1->Update();
            vtkUnstructuredGrid* submesh = th1->GetOutput();

            for( int i=0; i<submesh->GetNumberOfCells(); i++)
            {
                double cell_center_p[3];
                double pt[3];
                double weights[3];
                int not_used;
                submesh->GetCell(i)->GetParametricCenter(cell_center_p);
                submesh->GetCell(i)->EvaluateLocation(not_used,cell_center_p,pt,weights);

                double closestpoint[3];
                vtkIdType cellid;
                int subid;
                double dist2;
                in_target_locator->FindClosestPoint(pt,closestpoint,cellid,subid, dist2);

                labels->SetValue(cellid, label);
            }    
        }
    
    
        mesh->GetCellData()->AddArray(labels);
    }
    
    return regions_found?0:-1;
}



void FastMarchingGrowingFromBiggestTriangle(vtkPolyData* mesh, const char* arrayname)
{
    vtkDataArray* f = mesh->GetCellData()->GetScalars(arrayname);
    
    //find the triangle with max f
    unsigned long max_f_ind = 0;
    double f_max = f->GetTuple1(0);
    for(int i=1; i<f->GetNumberOfTuples(); i++)
    {
        if(f->GetTuple1(i)>f_max)
        {
            f_max = f->GetTuple1(i);
            max_f_ind = i;
        }           
    }
    
    
    //start growing the region from the found triangle
    MinHeap<unsigned long, double> heap;
    heap.Insert(0, max_f_ind); //put initial thing


    
    vtkSmartPointer<vtkLongArray> labels = vtkSmartPointer<vtkLongArray>::New();
    labels->SetNumberOfComponents(1);
    labels->SetNumberOfTuples(mesh->GetNumberOfCells());
    labels->SetName("GrowingOrder");
    
    //initialize
    for(int i=0; i<labels->GetNumberOfTuples(); i++) labels->SetTuple1(i,-1);
    
    unsigned long c=0;
    while (heap.Size()>0)
    {
        //get the cell
        unsigned long cell_id = heap.GetMin();
        heap.DeleteMin();
        
        labels->SetTuple1(cell_id, c++); 
        
        //get cell neighbors
        std::set<vtkIdType> neighbors;
        GetCellNeighbors(mesh, cell_id, neighbors);
    
        //add cells without label to the heap
        for(std::set<vtkIdType>::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
        {
            //std::cout<<"Neighbor "<<*it1<<std::endl;
            if(labels->GetTuple1(*it1)<0)
            {
                unsigned long neighbor_id = *it1;
                const double fdiff = fabs(f->GetTuple1(neighbor_id)-f->GetTuple1(cell_id));
                heap.Insert(fdiff, neighbor_id);
                labels->SetTuple1(*it1, 0);
            }
        }
        
        if(c%100==0)
        {
            std::cout<<"Cell: "<<c<<"/"<<mesh->GetNumberOfCells()<<". Heap: "<<heap.Size()<<"\n";
        }      
    }
    
    mesh->GetCellData()->AddArray(labels);
    
    std::cout<<std::endl;
}





void FastMarchingCumulativeFunction(vtkPolyData* mesh, const char* arrayname)
{
    vtkDataArray* f = mesh->GetCellData()->GetScalars(arrayname);
    
    //find the triangle with max f
    unsigned long max_f_ind = 0;
    double f_max = f->GetTuple1(0);
    for(int i=1; i<f->GetNumberOfTuples(); i++)
    {
        if(f->GetTuple1(i)>f_max)
        {
            f_max = f->GetTuple1(i);
            max_f_ind = i;
        }           
    }
    
    
    //start growing the region from the found triangle
    std::stack<unsigned long> cell_stack;
    cell_stack.push(max_f_ind); //put the initial cell


    
    vtkSmartPointer<vtkDoubleArray> labels = vtkSmartPointer<vtkDoubleArray>::New();
    labels->SetNumberOfComponents(1);
    labels->SetNumberOfTuples(mesh->GetNumberOfCells());
    labels->SetName("GrowingOrder");

    labels->SetTuple1(max_f_ind, f_max); 

    
    //initialize
    for(int i=0; i<labels->GetNumberOfTuples(); i++) labels->SetTuple1(i,1);
    
    unsigned long c=0;
    while (!cell_stack.empty())
    {
        //get the cell
        unsigned long cell_id = cell_stack.top();
        cell_stack.pop();
        
        
        //get cell neighbors
        std::set<vtkIdType> neighbors;
        GetCellNeighbors(mesh, cell_id, neighbors);
    
        //add cells without label to the heap
        for(std::set<vtkIdType>::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
        {
            //std::cout<<"Neighbor "<<*it1<<std::endl;
            if(labels->GetTuple1(*it1)>0)
            {
                unsigned long neighbor_id = *it1;
                const double fval = f->GetTuple1(neighbor_id)+f->GetTuple1(cell_id);
                cell_stack.push(neighbor_id);
                labels->SetTuple1(*it1, fval);
            }
        }
        
        c++;
        if(c%1000==0)
        {
            std::cout<<"Cell: "<<c<<"/"<<mesh->GetNumberOfCells()<<"\n";
        }      
    }
    
    mesh->GetCellData()->AddArray(labels);
    
    std::cout<<std::endl;
}






int GenerateMarkerSetByIterativeThresholdingSteps(vtkPolyData* mesh, const char* arrayname, const int nregions)
{
    vtkDataArray* f = mesh->GetCellData()->GetArray(arrayname);
    
    if(f==NULL)
    {
        std::cout<<"Cell array "<<arrayname<<" was not found"<<std::endl;
        throw;
    }

    vtkSmartPointer<vtkShortArray> labels = vtkSmartPointer<vtkShortArray>::New();
    labels->SetName(MARKER_SET);
    labels->SetNumberOfComponents(1);
    labels->SetNumberOfTuples(f->GetNumberOfTuples());

    //for each cell, compare f with the neighborhood, and mark the cell if the value is 
    //strictly smaller
    
    const unsigned long ncells = mesh->GetNumberOfCells();
    
    //find min and max of the f
    double* range = f->GetRange(0);
    const double f_min = range[0];
    const double f_max = range[1];
    
    std::cout<<"The values in the function are in the range: ["<<f_min<<";"<<f_max<<"]"<<std::endl;


    mesh->GetCellData()->SetActiveScalars(arrayname);
    
    vtkSmartPointer<vtkThreshold> th = vtkSmartPointer<vtkThreshold>::New();
    th->SetInputData(mesh);
    //th->SetAttributeModeToUseCellData();
    th->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, arrayname);
    
    vtkSmartPointer<vtkConnectivityFilter> conn = vtkSmartPointer<vtkConnectivityFilter>::New();
    conn->SetExtractionModeToAllRegions();
    conn->ColorRegionsOn();    
    
    double step = -0.1;
    double T = f_max+step; //thresholds
    
    bool regions_found = false;
    
    while( !regions_found && fabs(T-f_min)>1e-5 )
    {
        th->ThresholdByLower(T);
        th->Update();
        
        conn->SetInputData(th->GetOutput());
        conn->Update();    
        
        const unsigned long n_extr_regions = conn->GetNumberOfExtractedRegions();
        std::cout<<"T="<<T<<"; components="<<n_extr_regions<<std::endl;
        
        if(n_extr_regions >= nregions) 
        {   
            regions_found = true;
        }
        else 
        {
            T += step;
        }
    }

    if(regions_found)
    {
        //mark the cells based on the threshold filter (using locator)
        for (int i = 0; i < ncells; i++) labels->SetTuple1(i,0);

        vtkSmartPointer<vtkCellLocator> in_target_locator = vtkSmartPointer<vtkCellLocator>::New();
        in_target_locator->SetDataSet(mesh);
        in_target_locator->BuildLocator();

        for( int region=0; region<conn->GetNumberOfExtractedRegions(); region++)
        {
            const int label = region+1;


            vtkSmartPointer<vtkThreshold> th1 = vtkSmartPointer<vtkThreshold>::New();
            th1->SetInputData(conn->GetOutput());
            th1->ThresholdBetween(region,region);
            th1->Update();
            vtkUnstructuredGrid* submesh = th1->GetOutput();

            for( int i=0; i<submesh->GetNumberOfCells(); i++)
            {
                double cell_center_p[3];
                double pt[3];
                double weights[3];
                int not_used;
                submesh->GetCell(i)->GetParametricCenter(cell_center_p);
                submesh->GetCell(i)->EvaluateLocation(not_used,cell_center_p,pt,weights);

                double closestpoint[3];
                vtkIdType cellid;
                int subid;
                double dist2;
                in_target_locator->FindClosestPoint(pt,closestpoint,cellid,subid, dist2);

                labels->SetValue(cellid, label);
            }    
        }
    
    
        mesh->GetCellData()->AddArray(labels);
    }
    
    return regions_found?0:-1;
}




void FastMarchingWatershed(vtkPolyData* mesh, const char* arrayname)
{
    vtkDataArray* label = mesh->GetCellData()->GetArray(MARKER_SET);
    vtkDataArray* f = mesh->GetCellData()->GetArray("area_change");
      
    MinHeap<CellDescriptor, double> heap;
    
    for(int i=0; i<mesh->GetNumberOfCells(); i++)
    {
        if(label->GetTuple1(i)==1)
            GrowHeap(mesh, i, heap);
    }
    
    unsigned long cc = 0;
    while (heap.Size()>0)
    {
        CellDescriptor data = heap.GetMin();
        heap.DeleteMin();
        
        const unsigned long vi = data.cell_id;
        if(label->GetTuple1(vi)==0)
        {
            label->SetTuple1(vi, data.label);
            GrowHeap(mesh, vi, heap);
        }
     
        if(cc++%100)
        {
            std::cout<<"Heap: "<<heap.Size()<<"\n";
        }
    }
    std::cout<<std::endl;
}


void GrowHeap(vtkPolyData* mesh, unsigned long cellid, MinHeap<CellDescriptor, double> & heap)
{
    CellDescriptor data;
    vtkDataArray* label = mesh->GetCellData()->GetArray(MARKER_SET);
    vtkDataArray* f = mesh->GetCellData()->GetArray("area_change");
    
    //get cell neighbors
    std::set<vtkIdType> neighbors;
    GetCellNeighbors(mesh, cellid, neighbors);

    //add cells without label to the heap
    for(std::set<vtkIdType>::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
    {
        if(label->GetTuple1(*it1)==0)
        {
            data.cell_id = *it1;
            data.label = label->GetTuple1(cellid);
            const double fdiff = fabs(f->GetTuple1(*it1)-f->GetTuple1(cellid));
            heap.Insert(fdiff, data);
        }
    }    
}
