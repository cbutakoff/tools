/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#include "UmbrellaWeightedOrder2Smoother.h"
#include <vtkSmartPointer.h>
#include <vtkPolyDataWriter.h>



double UmbrellaWeightedOrder2Smoother::TriangleWeightScaleDependent(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2) const
{
    return 1/( (v1-v2).norm() );
}




void UmbrellaWeightedOrder2Smoother::GetVertexNeighbors( vtkIdType vertexId, VertexIDArrayType& neighbors)  {
    
    std::set<vtkIdType> connectedVertices;
    
    //get all cells that vertex 'seed' is a part of 
    vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
    m_originalMesh->GetPointCells(vertexId, cellIdList);

    //loop through all the cells that use the seed point 
    for (vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++) {

        vtkCell* cell = m_originalMesh->GetCell(cellIdList->GetId(i));

        for (vtkIdType e = 0; e < cell->GetNumberOfEdges(); e++) {
            vtkCell* edge = cell->GetEdge(e);

            vtkIdList* pointIdList = edge->GetPointIds();

            if (pointIdList->GetId(0) == vertexId || pointIdList->GetId(1) == vertexId) {
                if (pointIdList->GetId(0) == vertexId) {
                    connectedVertices.insert(pointIdList->GetId(1));
                }
                else {
                    connectedVertices.insert(pointIdList->GetId(0));
                }
            }
        }
    }
    
    neighbors.clear();
//    std::cout<<"Neighbors of "<<vertexId<<" "<<std::endl;
    for( std::set<vtkIdType>::const_iterator it=connectedVertices.begin(); it!=connectedVertices.end(); it++ )
    {
        neighbors.push_back( (*it) );
//        std::cout<<(*it)<<" "<<std::endl;
    }
}

static int idx = 0;
void UmbrellaWeightedOrder2Smoother::Update()
{
    vtkSmartPointer<vtkPolyDataWriter> wr = vtkSmartPointer<vtkPolyDataWriter>::New();
    wr->SetFileName("input_mesh.vtk");
    wr->SetInputData(m_originalMesh);
    wr->Write();

    CalculateConnectivity();
    
    //start iterative updater
    bool converged = false;
    
    Eigen::VectorXd U(m_C.size());
    Eigen::VectorXd U2(m_C.size());
    U.fill(-1); //initialize. Some values will not be initialized, so checking will be needed
    U2.fill(-1);
    
    int iter = 0;
    
    while (!converged)
    {
        //std::cout<<"Smoothing iter "<<iter<<std::endl;
        
        //calculate U for inner and boundary vertices
        MatrixUType U, U2;
        Eigen::SparseMatrix<double> weights;
        CalculateU(U, weights);
        
        //calculate U^2 only for interior vertices
        CalculateU2(U2, U, weights);
        
        //update only the interior vertices
        const double diff = UpdateMeshPoints(U2);
        

        char filename[100];
        sprintf(filename,"mesh%03d.vtk",idx++);
        m_originalMesh->BuildCells();
        vtkSmartPointer<vtkPolyDataWriter> wr = vtkSmartPointer<vtkPolyDataWriter>::New();
        wr->SetFileName(filename);
        wr->SetInputData(m_originalMesh);
        wr->Write();

        iter++;
        
        if( diff<m_tolerance || iter>=m_maxIter )
        {
            converged = true;
        }
        
    }
}


void UmbrellaWeightedOrder2Smoother::ClassifyVertex(VertexConnectivityType& v)
{
    VertexIDArrayType::iterator cover_it = std::find(m_coverVertexIDs->begin(), m_coverVertexIDs->end(), v.originalID);
    
    if( cover_it != m_coverVertexIDs->end() ) //the point is on the cover. See if at least one neighbor is from the exterior
    {
        v.vertexClass = vcInterior; //interior if not boundary
        
        
        //verify if it could be on the boundary
        for(VertexIDArrayType::const_iterator it = v.connectedVertices.begin(); 
                it!=v.connectedVertices.end(); it++)
        {
            const vtkIdType id = (*it);
            VertexIDArrayType::iterator neighb_it = std::find(m_coverVertexIDs->begin(), m_coverVertexIDs->end(), id);
            
            if( neighb_it == m_coverVertexIDs->end() )
            {
                v.vertexClass = vcBoundary;
                break;
            }
        }
    }
    else
    {
        v.vertexClass = vcExterior;
    }
}



void UmbrellaWeightedOrder2Smoother::CalculateConnectivity()
{
    //go over all cover vertices
    m_C.clear();
    
    for( VertexIDArrayType::const_iterator it = m_coverVertexIDs->begin(); it!=m_coverVertexIDs->end(); it++ )
    {
        //get neighbors for this vertex
        VertexConnectivityType vc; 
        vc.vertexClass = vcExterior; //default
        vc.originalID = (*it); 
        GetVertexNeighbors( (*it), vc.connectedVertices );
        ClassifyVertex( vc );
        
        m_C.push_back(vc);
    }

}



void UmbrellaWeightedOrder2Smoother::CalculateU(MatrixUType& U, Eigen::SparseMatrix<double>& weights)
{
    weights.resize(m_originalMesh->GetNumberOfPoints(), m_originalMesh->GetNumberOfPoints());
    U.resize(Eigen::NoChange, m_C.size());
    
    vtkIdType vertex_index=0;    
    for( std::vector<VertexConnectivityType>::const_iterator it=m_C.begin(); it!=m_C.end(); it++, vertex_index++ )
    {
        Eigen::VectorXd w(it->connectedVertices.size()); //weights

        Eigen::Vector3d v;
        m_originalMesh->GetPoint( it->originalID, v.data() );
        
        //iterate over all the neighbors and calculate weights
        vtkIdType neighb_index=0;
        Eigen::Vector3d vi;

        Eigen::Vector3d Un; Un.fill(0);
        
        for( VertexIDArrayType::const_iterator neighb_it = it->connectedVertices.begin();
                neighb_it != it->connectedVertices.end(); neighb_it++, neighb_index++)
        {
            const vtkIdType neighb_id = (*neighb_it);
            m_originalMesh->GetPoint( neighb_id, vi.data() );

            w(neighb_index) = TriangleWeightScaleDependent(v,vi);
            
            //save only upper triangular matrix
//            std::cout<<"Adding connectivity : "<<it->originalID<<", "<<neighb_id<<std::endl;
            weights.coeffRef(std::min(it->originalID, neighb_id), std::max(it->originalID, neighb_id)) = w(neighb_index);
            
            Un += vi*w(neighb_index);
        }
        
        const double Wt = w.sum();
        
        U.col(vertex_index) = -v + Un/Wt;
    }
}


void UmbrellaWeightedOrder2Smoother::CalculateU2(MatrixUType& U2, const MatrixUType& U, const Eigen::SparseMatrix<double>& weights)
{
    U2.resize(Eigen::NoChange, m_C.size());

        
    vtkIdType vertex_index=0;
    for( std::vector<VertexConnectivityType>::const_iterator it=m_C.begin(); it!=m_C.end(); it++, vertex_index++ )
    {
        Eigen::VectorXd w(it->connectedVertices.size()); //weights

        Eigen::Vector3d Uv = U.col(vertex_index);
        
        
        //iterate over all the neighbors and calculate weights
        vtkIdType neighb_index=0;
        Eigen::Vector3d Uvi;

        Eigen::Vector3d Un; Un.fill(0);
        
        for( VertexIDArrayType::const_iterator neighb_it = it->connectedVertices.begin();
                neighb_it != it->connectedVertices.end(); neighb_it++, neighb_index++)
        {
            const vtkIdType neighb_id = (*neighb_it);

            Uvi = U.col(neighb_index);
            w(neighb_index) = weights.coeff( std::min(it->originalID, neighb_id), std::max(it->originalID, neighb_id) );
            
            Un += Uvi*w(neighb_index);
        }
        
        const double Wt = w.sum();
        
        U2.col(vertex_index) = -Uv + Un/Wt;
    }    
}



double UmbrellaWeightedOrder2Smoother::UpdateMeshPoints( const MatrixUType& U2 )
{
    double vertex_difference = 0; //how much the update displaced the vertices

    vtkIdType vertex_index=0;
    
    for( std::vector<VertexConnectivityType>::const_iterator it=m_C.begin(); it!=m_C.end(); it++, vertex_index++ )
    {
        if(it->vertexClass==vcInterior)
        {
            Eigen::Vector3d Pi;
            m_originalMesh->GetPoint( it->originalID, Pi.data() );
            
            //iterate over all the neighbors and calculate valences
            vtkIdType neighb_index=0;
            double Vsum = 0; //valence factor
            for( VertexIDArrayType::const_iterator neighb_it = it->connectedVertices.begin();
                    neighb_it != it->connectedVertices.end(); neighb_it++, neighb_index++)
            {
                const vtkIdType neighb_id = (*neighb_it);
                
                //find neighb_id in m_C and get its valence
                vtkIdType neighb_C_id = FindVertexConnectivityInfo(neighb_id);
                Vsum += 1/m_C.at(neighb_C_id).connectedVertices.size();
                
            }
            
            const double V = 1 + Vsum/it->connectedVertices.size();
            
            Eigen::Vector3d Pi_new = Pi - U2.col(vertex_index)/V;
            m_originalMesh->GetPoints()->SetPoint( it->originalID, Pi_new.data() );
            
            vertex_difference += (Pi-Pi_new).norm();
        }
    }
    
    return vertex_difference;
}




vtkIdType UmbrellaWeightedOrder2Smoother::FindVertexConnectivityInfo( vtkIdType id ) //uses ids within mesh, compares to originalID
{
    vtkIdType vertex_index=0;
    for( std::vector<VertexConnectivityType>::const_iterator it=m_C.begin(); it!=m_C.end(); it++, vertex_index++ )
    {
        if(it->originalID==id)
            break;
    }    
    
    return vertex_index;
}
