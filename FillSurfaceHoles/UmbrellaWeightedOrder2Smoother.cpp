/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#include "UmbrellaWeightedOrder2Smoother.h"
#include <vtkSmartPointer.h>



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
    for( std::set<vtkIdType>::const_iterator it=connectedVertices.begin(); it!=connectedVertices.end(); it++ )
        neighbors.push_back( (*it) );
}

void UmbrellaWeightedOrder2Smoother::Update()
{
    CalculateConnectivity();
    
    vtkPoints* pts = m_originalMesh->GetPoints();
    //start iterative updater
    bool converged = false;
    
    Eigen::VectorXd U(m_C.size());
    Eigen::VectorXd U2(m_C.size());
    U.fill(-1); //initialize. Some values will not be initialized, so checking will be needed
    U2.fill(-1);
    while (!converged)
    {
        //calculate U for inner and boundary vertices
        MatrixUType U, U2;
        Eigen::SparseMatrix<double> weights;
        CalculateU(U, weights);
        
        //calculate U^2 only for interior vertices
        CalculateU2(U2, U, weights);
        
        //update only the interior vertices
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
                neighb_it != it->connectedVertices.end(); it++, neighb_index++)
        {
            const vtkIdType neighb_id = (*neighb_it);
            m_originalMesh->GetPoint( neighb_id, vi.data() );

            w(neighb_index) = TriangleWeightScaleDependent(v,vi);
            
            //save only upper triangular matrix
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
                neighb_it != it->connectedVertices.end(); it++, neighb_index++)
        {
            const vtkIdType neighb_id = (*neighb_it);

            Uvi = U.col(neighb_id);
            w(neighb_index) = weights.coeff( std::min(it->originalID, neighb_id), std::max(it->originalID, neighb_id) );
            
            Un += Uvi*w(neighb_index);
        }
        
        const double Wt = w.sum();
        
        U2.col(vertex_index) = -Uv + Un/Wt;
    }    
}



void UmbrellaWeightedOrder2Smoother::UpdateMeshPoints( const MatrixUType& U2 )
{
    
}
