/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkShortArray.h>
#include <vtkDataArray.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>


#include <armadillo>
#include <stack>
#include <set>
#include <unordered_map>

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <limits>
#include <queue>

typedef arma::Col<char> LabelArrayType;
typedef std::set<vtkIdType> NeighborSetType;
typedef std::unordered_map<vtkIdType, double> PropagationFrontListType;

void CalculateDTPerBlobID(vtkPolyData* mesh, vtkIdType blob_id, const char* blob_array, arma::rowvec& dt);

void GetPointNeighbors(vtkPolyData* mesh, vtkIdType id, NeighborSetType& neighbors, LabelArrayType& labels);

//if the vertex exists - update, if not - insert
inline vtkIdType PopSmallestFromFront(PropagationFrontListType& front);

int main(int argc, char** argv)
{
    std::cout<<"Calculate distance transform on a mesh from blobs to other vertices using Fast Marching v1.0"<<std::endl;
    
    if(argc<4)
    {
        std::cout<<std::endl<<"Usage: DistanceTransformMesh <input_mesh.vtk> <blob_point_shortarray> <output_mesh.vtk>"<<std::endl<<std::endl;
        std::cout<<"blob_point_shortarray - must be vtkShortArray stored with points."
                <<"Values are continuous indices from 0 to N"<<std::endl;
        std::cout<<"Blobs = regions of the mesh with the same id. "
                <<"This function calculates the distance from these blobs to every vertex"<<std::endl;
        std::cout<<"Vertices that do not belong to any blob should be marked "
                <<"with -1 in the blob_point_shortarray."<<std::endl;
        std::cout<<"The distances are stored in point array DT as float vectors "
                <<"(vtkFloatArray), i-th component corresponds to i-th blob"<<std::endl;
        return -1;
    }
    
    //const char* filename = "/home/costa/Dropbox/code/DistTransform/build/vein_rpv_sup.vtk";
    //const char* blob_array = "blob";
    const char* input_filename = argv[1];
    const char* blob_array = argv[2];
    const char* output_filename = argv[3];

    std::cout<<"Processing mesh: "<<input_filename<<std::endl;
    std::cout<<"Array with blob ids: "<<blob_array<<std::endl;
    std::cout<<"Output will be stored in: "<<output_filename<<std::endl;
    
    vtkSmartPointer<vtkPolyDataReader> rd = vtkSmartPointer<vtkPolyDataReader>::New();
    rd->SetFileName(input_filename);
    rd->Update();
    vtkPolyData* mesh = rd->GetOutput();
    
    vtkShortArray* blob_ids = (vtkShortArray*)mesh->GetPointData()->GetArray(blob_array);
    
    if(blob_ids==NULL)
    {
        std::cout<<"Point array "<<blob_array<<" not found. Aborting"<<std::endl;
        return -1;
    }
    
    short* range = blob_ids->GetValueRange();
    const vtkIdType n_blobs = range[1]+1;
    const vtkIdType n_points = mesh->GetNumberOfPoints();

    arma::mat dts(n_blobs, n_points);
    
    for(vtkIdType i=0; i<n_blobs; i++)
    {
        std::cout<<"Processing blob "<<i<<"/"<<n_blobs<<std::endl;
        
        arma::rowvec dt(n_points);
        CalculateDTPerBlobID(mesh, i, blob_array, dt);
        dts.row(i) = dt;
    }

    std::cout<<"Saving "<<output_filename<<std::endl;
    
    vtkSmartPointer<vtkFloatArray> dts_vtk = vtkSmartPointer<vtkFloatArray>::New();
    dts_vtk->SetNumberOfComponents(n_blobs);
    dts_vtk->SetNumberOfTuples(n_points);
    dts_vtk->SetName("DT");
    
    for(vtkIdType i=0; i<n_points; i++)
    {
        double *tuple = dts.colptr(i); 
        dts_vtk->SetTuple(i, tuple);
    }
    
    mesh->GetPointData()->AddArray(dts_vtk);
    
    vtkSmartPointer<vtkPolyDataWriter> wr = vtkSmartPointer<vtkPolyDataWriter>::New();
    wr->SetFileName(output_filename);
    wr->SetFileTypeToBinary();
    wr->SetInputData(mesh);
    wr->Write();
    
    
    return 0;
}


inline void UpdateFront(vtkPolyData* mesh, arma::rowvec& dt, LabelArrayType& labels, PropagationFrontListType& front, vtkIdType i){
    double Pi[3];
    mesh->GetPoint(i, Pi);
    
    //distance stored at i
    const double Di = dt[i];
    
    //  Get neighbors k of i that are not-accepted
    NeighborSetType neighbors;
    GetPointNeighbors(mesh, i, neighbors, labels);
    
    //  For every k 
    for(NeighborSetType::iterator it = neighbors.begin(); it != neighbors.end(); it++)
    {        
        //      calculate distance Dik to k 
        const vtkIdType k = *it;
        double Pk[3];
        mesh->GetPoint(k, Pk);
        const double Dik = sqrt((Pk[0]-Pi[0])*(Pk[0]-Pi[0]) + 
        (Pk[1]-Pi[1])*(Pk[1]-Pi[1]) + 
        (Pk[2]-Pi[2])*(Pk[2]-Pi[2]) );
        
        //      calculate the overall distance of vertex k: Dk = Di + Dik
        const double Dk_new = Di + Dik;
        
        //      If new Dk is smaller than old Dk - replace 
        if(Dk_new < dt[k]) dt[k] = Dk_new;
        
        //      If k was labeled as far, update the label to considered.
        labels[k] = 1;
        
        //      Add k to the list
        front[k] = dt[k];
    }
}

void CalculateDTPerBlobID(vtkPolyData* mesh, vtkIdType blob_id, const char* blob_array, arma::rowvec& dt)
{
    dt.fill(std::numeric_limits<double>::max());
    LabelArrayType labels(dt.size()); // 2-far, 1-front, 0-accepted 
    labels.fill(2);
    
    vtkShortArray* blob_ids = (vtkShortArray*)mesh->GetPointData()->GetArray(blob_array);
    
    //find the vertices that belong to the blob in question and add them to the queue
    std::queue<vtkIdType> stack; //will contain accepted vertices
    
    const long int n_points = mesh->GetNumberOfPoints();
    for(vtkIdType i=0; i<n_points; i++)
    {
        const vtkIdType local_blob_id = blob_ids->GetTuple1(i);
        if(local_blob_id==blob_id) 
        {
            stack.push(i);
            labels[i] = 0; //distance 
            dt[i] = 0; //accepted
        }
    }

    //create the front (first layer of "considered vertices") 
    //Initialize the state of all vertices to 2 (unexplored), except for the boundary point S to 1 (front).
    //interior of the blob is 0 - accepted
    PropagationFrontListType front;
    
    std::cout<<"Creating propagation front"<<std::endl;
    while(!stack.empty())
    {
        //get the topmost accepted vertex
        const vtkIdType i = stack.front();
        stack.pop();
       
        //get the point i
        UpdateFront(mesh, dt, labels, front, i);
    }
    
    
    //until front is empty
    std::cout<<"Running fast marching"<<std::endl;
    while(!front.empty())
    {
        std::cout<<"Remaining elements: "<<front.size()<<"\r"<<std::flush;
                
        //  Pop i with smallest distance
        vtkIdType i = PopSmallestFromFront(front);
        
        //  Label i as accepted.
        labels[i] = 0;

        //get the point i
        UpdateFront(mesh, dt, labels, front, i);
    }
    
    std::cout<<std::endl;
   
}



void GetPointNeighbors(vtkPolyData* mesh, vtkIdType id, NeighborSetType& neighbors, LabelArrayType& labels)
{
    vtkSmartPointer<vtkIdList> cell_ids = vtkSmartPointer<vtkIdList>::New();
    mesh->GetPointCells(id, cell_ids);
    
    neighbors.clear();
    
    vtkSmartPointer<vtkIdList> pt_ids = vtkSmartPointer<vtkIdList>::New();
    for(vtkIdType i=0; i<cell_ids->GetNumberOfIds(); i++)
    {
        pt_ids->Reset();
        mesh->GetCellPoints(cell_ids->GetId(i), pt_ids);
        
        for(vtkIdType j=0; j<pt_ids->GetNumberOfIds(); j++)
        {
            const vtkIdType nb_id = pt_ids->GetId(j);
            if( nb_id!=id && labels[nb_id]!=0 ) //not current point and not accepted
                neighbors.insert(nb_id);
        }
        
    }
}




inline vtkIdType PopSmallestFromFront(PropagationFrontListType& front)
{
    double min_dist = std::numeric_limits<double>::max();
    vtkIdType min_id = 0;
    
    for(auto& it : front)
    {
        if(it.second < min_dist)
        {
            min_dist = it.second;
            min_id = it.first;
        }
    }
    
    front.erase(min_id);
    return min_id;
}
