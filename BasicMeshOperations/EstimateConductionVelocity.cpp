/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
// Given a mesh (unstructured grid) and 3 images with x,y,z fibers, 
// create fibers at mesh vertices by sampling the images
// images are loaded sequentially, to reduce memory usage

#include <vtkSmartPointer.h>
#include <vtkDataSetReader.h>
#include <vtkXMLMultiBlockDataReader.h>
#include <vtkDataSetWriter.h>
#include <vtkProbeFilter.h>
#include <vtkCell.h>
#include <vtkPoints.h>
#include <vtkGradientFilter.h>
#include <vtkPointData.h>
#include <vtkProbeFilter.h>
#include <vtkPolyData.h>
#include <vtkFloatArray.h>
#include <vtkCompositeDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataArray.h>
#include <vtkCellData.h>

#include <VTKCommonTools.h>
#include <vtkCallbackCommand.h>

#include <iostream>

#include <Eigen/Dense>
#include <vtkMultiBlockDataSet.h>

//WOrks on point data
//Adds float array "IsoDerivative"
//velocity by averaging within each cell
void EstimateConductionVelocity(vtkDataSet* mesh, const char* isoch_array);

int main(int argc, char** argv)
{
    if(argc<3)
    {
        std::cout<<"Usage: EstimateConductionVelocity mesh.vtk isochrone_array_name outmesh.vtk"<<std::endl;
        exit(-1);
    }

    const char * arrayname = "Fibers";
    
    int c=1;
    const char* meshfilename = argv[c++];
    const char* isochrone_array_name = argv[c++];
    const char* outfilename = argv[c++];
    
    std::cout<<"Mesh: "<<meshfilename<<std::endl;
    std::cout<<"Output: "<<outfilename<<std::endl;
    std::cout<<"Isochrone array: "<<isochrone_array_name<<std::endl;
    
    vtkDataSet* mesh;
        
    vtkSmartPointer<vtkXMLMultiBlockDataReader> mesh_reader_xml = vtkSmartPointer<vtkXMLMultiBlockDataReader>::New();
    vtkSmartPointer<vtkDataSetReader> mesh_reader = vtkSmartPointer<vtkDataSetReader>::New();

    if( mesh_reader_xml->CanReadFile(meshfilename) )
    {
        CommonTools::AssociateProgressFunction(mesh_reader_xml);
        std::cout<<"Reading mesh: "<<meshfilename<<std::endl;
        mesh_reader_xml->SetFileName(meshfilename);
        mesh_reader_xml->Update();

        
        vtkMultiBlockDataSet* block = vtkMultiBlockDataSet::SafeDownCast(mesh_reader_xml->GetOutput());        
    
        mesh = vtkDataSet::SafeDownCast( block->GetBlock(0) );
        
    }
    else
    {
        CommonTools::AssociateProgressFunction(mesh_reader);
        std::cout<<"Reading mesh: "<<meshfilename<<std::endl;
        mesh_reader->SetFileName(meshfilename);
        mesh_reader->Update();        
        mesh  = mesh_reader->GetOutput();
    }
    
    EstimateConductionVelocity(mesh, isochrone_array_name);
    
    
    vtkSmartPointer<vtkDataSetWriter> wr = vtkSmartPointer<vtkDataSetWriter>::New();
    CommonTools::AssociateProgressFunction(wr);
    std::cout<<"Output filename: "<<outfilename<<std::endl;
    wr->SetFileName(outfilename);
    wr->SetFileTypeToBinary();
    wr->SetInputData(mesh);
    wr->Write();
}




void EstimateConductionVelocity(vtkDataSet* mesh, const char* isoch_array)
{
    vtkDataArray* isoc = mesh->GetPointData()->GetArray(isoch_array);
    vtkSmartPointer<vtkFloatArray> velocity = vtkSmartPointer<vtkFloatArray>::New();
    velocity->SetName("Velocity");
    velocity->SetNumberOfComponents(1);
    velocity->SetNumberOfTuples(mesh->GetNumberOfCells());
        
    
    
    for(vtkIdType cellid=0; cellid<mesh->GetNumberOfCells(); cellid++)
    {
        if( cellid % 100000 == 0 )
        {
            std::cout<<"Processing point "<<cellid<<"/"<<mesh->GetNumberOfCells()<<"\r"<<std::flush;
        }        
        

        vtkCell *cell = mesh->GetCell(cellid);
        
        Eigen::Vector3d direction;

        //calculate the gradient within the cell
        double pcenter[3];
        int subId = cell->GetParametricCenter(pcenter);
        

        //save the isochrone values in the same order as the vertices
        Eigen::VectorXd iso_values(cell->GetNumberOfPoints());
        for(int i=0; i<iso_values.size(); i++) iso_values[i] = isoc->GetTuple1(cell->GetPointId(i));
                
        cell->Derivatives(0, pcenter, iso_values.data(), 1, direction.data());
        
        direction.normalize();
        //std::cout<<"Cell: "<<cellid<<std::endl;
        //std::cout<<"Gradient: "<<direction<<std::endl;
        
 
        //project all the vertices onto the gradient direction
        //get cell center
        Eigen::Vector3d center;
        Eigen::VectorXd weights(cell->GetNumberOfPoints());
        cell->EvaluateLocation(subId, pcenter, center.data(), weights.data());
        
        //prepare storage for the projection coordinates.
        //center of the cell will be the origin
        Eigen::VectorXd vertex_projections( cell->GetNumberOfPoints() );
        
        //find the longest projection of the direction on the edges
        for(int ptid = 0; ptid < cell->GetNumberOfPoints(); ptid++)
        {
            Eigen::Vector3d vertex; 
            mesh->GetPoint(cell->GetPointId(ptid), vertex.data());

            Eigen::Vector3d vertex_vector = vertex-center;
            vertex_projections[ptid] = vertex_vector.dot(direction);
        }

        //find the projections with the smallest and largest coordinate
        int max_vertex_id = 0;
        int min_vertex_id = 1;
        double max_time = iso_values[max_vertex_id];
        double min_time = iso_values[min_vertex_id];      
        
        //std::cout<<"Projection Coordinates: "<<vertex_projections<<std::endl;
        //std::cout<<"Time values: "<<iso_values<<std::endl;
        
        for(int i=0; i<iso_values.size(); i++)
        {
            if( max_time<iso_values[i] )
            {
                max_time = iso_values[i];
                max_vertex_id = i;
            }

            if( min_time>iso_values[i] )
            {
                min_time = iso_values[i];
                min_vertex_id = i;
            }
                
        }
        
        //std::cout<<"Min id: "<<min_vertex_id<<", max id: "<<max_vertex_id<<std::endl; 
        //std::cout<<"Min time:"<<min_time<<", max time: "<<max_time<<std::endl;
        
        double v=0;
        double max_time_diff = max_time-min_time;
        double coordinate_diff = fabs(vertex_projections[max_vertex_id] - vertex_projections[min_vertex_id]);
        if(max_time_diff>0)
        {
            v = coordinate_diff/max_time_diff;
            //std::cout<<"Velocity:"<<v<<std::endl;
        }
        
        velocity->SetTuple1(cellid, v);
    }
    std::cout<<std::endl;
    
    mesh->GetCellData()->AddArray(velocity);
}


