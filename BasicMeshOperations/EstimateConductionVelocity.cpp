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
//velocity by sampling -- very slow
void EstimateConductionVelocityOld(vtkDataSet* mesh, const char* isoch_array);
double GetDifferentiationStep(vtkDataSet* mesh);

//by default use order 5 (no reason)
//http://web.media.mit.edu/~crtaylor/calculator.html
void CreateSamplingPositions( Eigen::Vector3d& center, Eigen::Vector3d& direction, double h, vtkPoints* pts );
double CalculateDerivative( vtkDataArray* samples, double h );

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
        cell->GetParametricCenter(pcenter);
        
        
        Eigen::VectorXd iso_values(cell->GetNumberOfPoints());
        for(int i=0; i<iso_values.size(); i++) iso_values[i] = isoc->GetTuple1(cell->GetPointId(i));
                
        cell->Derivatives(0, pcenter, iso_values.data(), 1, direction.data());
        
        direction.normalize();
        //std::cout<<"Cell: "<<cellid<<std::endl;
        //std::cout<<"Gradient: "<<direction<<std::endl;
        
        
        //find the longest projection of the direction on the edges
        double proj_length = 0;
        Eigen::Vector3d pt0l, pt1l; //values for the longest projection
        double t0l, t1l; //values for the longest projection

        for(int edgeid = 0; edgeid < cell->GetNumberOfEdges(); edgeid++)
        {
            vtkCell* edge = cell->GetEdge(edgeid);
            Eigen::Vector3d pt0, pt1; 
            mesh->GetPoint(edge->GetPointId(0), pt0.data());
            mesh->GetPoint(edge->GetPointId(1), pt1.data());

            Eigen::Vector3d edge_vector = (pt0-pt1).normalized();
            const double dp = fabs(edge_vector.dot(direction));
            
            if( dp>proj_length )
            {
                proj_length = dp;
                //std::cout<<"Pt id0: "<<edge->GetPointId(0)<<std::endl;
                //std::cout<<"Pt id1: "<<edge->GetPointId(1)<<std::endl;
                t0l = isoc->GetTuple1(edge->GetPointId(0));
                t1l = isoc->GetTuple1(edge->GetPointId(1));
                pt0l = pt0;
                pt1l = pt1;
            }
        }

        //std::cout<<"Max pt0: "<<pt0l<<std::endl;
        //std::cout<<"Max pt1: "<<pt1l<<std::endl;
        //std::cout<<"Max t0: "<<t0l<<std::endl;
        //std::cout<<"Max t1: "<<t1l<<std::endl;
        
        //use the longest projection to estimate velocity along the chosen edge
        const double edge_length = (pt0l-pt1l).norm();
        const double time_edge = fabs(t0l-t1l); //ignore the sign, velocity is taken along the gradient, so the time is always increasing and velocity always positive
        
        double v=0;
        if(time_edge>0)
        {
            const double v_edge = edge_length/time_edge;
            //std::cout<<"Velocoity edge: "<<v_edge<<std::endl;
            //std::cout<<"Dot product:"<<proj_length<<std::endl;
            v = v_edge/proj_length; //from right triangle: v - velocity, v_e - velocity along the edge, proj_length - cos of the angle between
            //std::cout<<"Velocity:"<<v<<std::endl;
        }
        
        velocity->SetTuple1(cellid, v);
    }
    std::cout<<std::endl;
    
    mesh->GetCellData()->AddArray(velocity);
}





void EstimateConductionVelocityOld(vtkDataSet* mesh, const char* isoch_array)
{
    std::cout<<"Estimating differentiation step"<<std::endl;
    const double h = GetDifferentiationStep(mesh);
    
    std::cout<<"Edge length:"<<h<<std::endl;
    
    //calculate the gradient of the isochrones
    vtkSmartPointer<vtkGradientFilter> grad_filt = vtkSmartPointer<vtkGradientFilter>::New();
    CommonTools::AssociateProgressFunction(grad_filt);
    grad_filt->SetInputData(mesh);
    grad_filt->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, isoch_array);
    grad_filt->ComputeDivergenceOff();
    grad_filt->ComputeGradientOn();
    grad_filt->ComputeQCriterionOff();
    grad_filt->ComputeVorticityOff();
    grad_filt->Update();
    
    std::cout<<"Getting gradients"<<std::endl;
    vtkDataArray *grad = grad_filt->GetOutput()->GetPointData()->GetArray("Gradients");
    //std::cout<<grad->GetNumberOfTuples()<<std::endl;
    
    
    std::cout<<"Setting up prober"<<std::endl;
    vtkSmartPointer<vtkProbeFilter> prober = vtkSmartPointer<vtkProbeFilter>::New();
    prober->SetSourceData(mesh);
    prober->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS, isoch_array);
    
    
    std::cout<<"Creating array"<<std::endl;
    vtkSmartPointer<vtkFloatArray> derivs = vtkSmartPointer<vtkFloatArray>::New();
    derivs->SetName("Velocity");
    derivs->SetNumberOfComponents(1);
    derivs->SetNumberOfTuples(mesh->GetNumberOfPoints());
    
    //for every mesh vertex 
    for(vtkIdType ptid = 0; ptid<mesh->GetNumberOfPoints(); ptid++ )
    {
        if( ptid % 1000 == 0 )
        {
            std::cout<<"Processing point "<<ptid<<"/"<<mesh->GetNumberOfPoints()<<"\r"<<std::flush;
        }
        

        Eigen::Vector3d center;
        mesh->GetPoint(ptid, center.data());
        
        Eigen::Vector3d direction;
        grad->GetTuple(ptid, direction.data());
        
        vtkSmartPointer<vtkPoints> samples = vtkSmartPointer<vtkPoints>::New();
        CreateSamplingPositions(center, direction, h, samples);
//        for(int i=0; i<samples->GetNumberOfPoints(); i++)
//            std::cout<<samples->GetPoint(i)<<std::endl;
        
        vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
        pd->SetPoints(samples);
        prober->SetInputData(pd);
        prober->Update();
        
        const double d = CalculateDerivative(prober->GetOutput()->GetPointData()->GetScalars(), h);
        derivs->SetTuple1(ptid, d);
    }
    std::cout<<std::endl;
    
    mesh->GetPointData()->AddArray(derivs);
}




double GetDifferentiationStep(vtkDataSet* mesh)
{
    //For now I assume that all the mesh edges have equal length
    Eigen::Vector3d pt0, pt1;

    mesh->GetPoint( mesh->GetCell(0)->GetPointIds()->GetId(0), pt0.data() );
    mesh->GetPoint( mesh->GetCell(0)->GetPointIds()->GetId(1), pt1.data() );

    return (pt0-pt1).norm();
}




void CreateSamplingPositions( Eigen::Vector3d& center, Eigen::Vector3d& direction, double h, vtkPoints* pts )
{
    double displacement[] = {-2,2,-1,1};
    pts->SetNumberOfPoints(2);

    direction.normalize();
    
    for(int i=0; i<2; i++)
    {
        Eigen::Vector3d pt = center + direction*h*displacement[i];
        pts->SetPoint(i, pt.data());
    }
}




double CalculateDerivative( vtkDataArray* samples, double h )
{
    double result = 0;
    
    double diff = (samples->GetTuple1(0)-samples->GetTuple1(1));
    result += 4*h/diff;

    //std::cout<<diff<<std::endl;
    
    diff = (samples->GetTuple1(2)-samples->GetTuple1(3));
    result += 2*h/diff;
    //std::cout<<diff<<std::endl;
    
    return result/2;
}
