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



double UmbrellaWeightedOrder2Smoother::CalculateEdgeWeight(vtkIdType v_third, vtkIdType v1_edge, vtkIdType v2_edge) const
{
    Eigen::VectorXd v_third_p, v1_edge_p, v2_edge_p;
    m_originalMesh->GetPoint( v_third, v_third_p.data() );
    m_originalMesh->GetPoint( v1_edge, v1_edge_p.data() );
    m_originalMesh->GetPoint( v2_edge, v2_edge_p.data() );
    
    return CalculateEdgeWeight( v_third_p, v1_edge_p, v2_edge_p );
}


double UmbrellaWeightedOrder2Smoother::CalculateEdgeWeight(const Eigen::VectorXd& v_third, const Eigen::VectorXd& v1_edge, const Eigen::VectorXd& v2_edge) const
{
    double value = -1;
    switch( m_weightingType )
    {
        case vwCotangent:
            value = EdgeWeightCotangent(v_third, v1_edge, v2_edge);
            break;
        case vwInvEdgeLength:
            value = EdgeWeightInvEdgeLength(v_third, v1_edge, v2_edge);
            break;
        default:
            std::cout<<"Undefined vertex weighting type";
            throw "Undefined vertex weighting type";
    }
}

double UmbrellaWeightedOrder2Smoother::EdgeWeightCotangent(const Eigen::VectorXd& v_third, const Eigen::VectorXd& v1_edge, const Eigen::VectorXd& v2_edge) const
{
    const VectorType v13 = (v2_edge-v_third).normalized();
    const VectorType v12 = (v1_edge-v_third).normalized();
    const double angle = std::acos(v13.dot(v12));
    const double w23 = 1/std::tan(angle);
    return w23;
}




double UmbrellaWeightedOrder2Smoother::EdgeWeightInvEdgeLength(const Eigen::VectorXd& v_third, const Eigen::VectorXd& v1_edge, const Eigen::VectorXd& v2_edge) const
{
    const double w23 = 1/( (v1_edge-v2_edge).norm() );
}







void UmbrellaWeightedOrder2Smoother::GetVertexNeighbors( vtkIdType vertexId, VertexIDArrayType& neighbors) const {
    
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
    
    SparseMatrixDoubleType W;
    
    CalculateEdgeWeightMatrix(W);

    std::cout<<W<<std::endl;
    
    char filename[100];
    sprintf(filename,"mesh%03d.vtk",idx++);
    m_originalMesh->BuildCells();

    wr->SetFileName(filename);
    wr->SetInputData(m_originalMesh);
    wr->Write();

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


//
//void UmbrellaWeightedOrder2Smoother::CalculateU(MatrixUType& U, Eigen::SparseMatrix<double>& weights)
//{
//    weights.resize(m_originalMesh->GetNumberOfPoints(), m_originalMesh->GetNumberOfPoints());
//    U.resize(Eigen::NoChange, m_C.size());
//    
//    vtkIdType vertex_index=0;    
//    for( VertexConnectivityArrayType::const_iterator it=m_C.begin(); it!=m_C.end(); it++, vertex_index++ )
//    {
//        Eigen::VectorXd w(it->connectedVertices.size()); //weights
//
//        Eigen::Vector3d v;
//        m_originalMesh->GetPoint( it->originalID, v.data() );
//        
//        //iterate over all the neighbors and calculate weights
//        vtkIdType neighb_index=0;
//        Eigen::Vector3d vi;
//
//        Eigen::Vector3d Un; Un.fill(0);
//        
//        for( VertexIDArrayType::const_iterator neighb_it = it->connectedVertices.begin();
//                neighb_it != it->connectedVertices.end(); neighb_it++, neighb_index++)
//        {
//            const vtkIdType neighb_id = (*neighb_it);
//            m_originalMesh->GetPoint( neighb_id, vi.data() );
//
//            w(neighb_index) = TriangleWeightScaleDependent(v,vi);
//            
//            //save only upper triangular matrix
////            std::cout<<"Adding connectivity : "<<it->originalID<<", "<<neighb_id<<std::endl;
//            weights.coeffRef(std::min(it->originalID, neighb_id), std::max(it->originalID, neighb_id)) = w(neighb_index);
//            
//            Un += vi*w(neighb_index);
//        }
//        
//        const double Wt = w.sum();
//        
//        U.col(vertex_index) = -v + Un/Wt;
//    }
//}

//
//void UmbrellaWeightedOrder2Smoother::CalculateU2(MatrixUType& U2, const MatrixUType& U, const Eigen::SparseMatrix<double>& weights)
//{
//    U2.resize(Eigen::NoChange, m_C.size());
//
//        
//    vtkIdType vertex_index=0;
//    for( VertexConnectivityArrayType::const_iterator it=m_C.begin(); it!=m_C.end(); it++, vertex_index++ )
//    {
//        Eigen::VectorXd w(it->connectedVertices.size()); //weights
//
//        Eigen::Vector3d Uv = U.col(vertex_index);
//        
//        
//        //iterate over all the neighbors and calculate weights
//        vtkIdType neighb_index=0;
//        Eigen::Vector3d Uvi;
//
//        Eigen::Vector3d Un; Un.fill(0);
//        
//        for( VertexIDArrayType::const_iterator neighb_it = it->connectedVertices.begin();
//                neighb_it != it->connectedVertices.end(); neighb_it++, neighb_index++)
//        {
//            const vtkIdType neighb_id = (*neighb_it);
//
//            Uvi = U.col(neighb_index);
//            w(neighb_index) = weights.coeff( std::min(it->originalID, neighb_id), std::max(it->originalID, neighb_id) );
//            
//            Un += Uvi*w(neighb_index);
//        }
//        
//        const double Wt = w.sum();
//        
//        U2.col(vertex_index) = -Uv + Un/Wt;
//    }    
//}
//

//
//double UmbrellaWeightedOrder2Smoother::UpdateMeshPoints( const MatrixUType& U2 )
//{
//    double vertex_difference = 0; //how much the update displaced the vertices
//
//    vtkIdType vertex_index=0;
//    
//    for( VertexConnectivityArrayType::const_iterator it=m_C.begin(); it!=m_C.end(); it++, vertex_index++ )
//    {
//        if(it->vertexClass==vcInterior)
//        {
//            Eigen::Vector3d Pi;
//            m_originalMesh->GetPoint( it->originalID, Pi.data() );
//            
//            //iterate over all the neighbors and calculate valences
//            vtkIdType neighb_index=0;
//            double Vsum = 0; //valence factor
//            for( VertexIDArrayType::const_iterator neighb_it = it->connectedVertices.begin();
//                    neighb_it != it->connectedVertices.end(); neighb_it++, neighb_index++)
//            {
//                const vtkIdType neighb_id = (*neighb_it);
//                
//                //find neighb_id in m_C and get its valence
//                vtkIdType neighb_C_id = FindVertexConnectivityInfo(neighb_id);
//                Vsum += 1/m_C.at(neighb_C_id).connectedVertices.size();
//                
//            }
//            
//            const double V = 1 + Vsum/it->connectedVertices.size();
//            
//            Eigen::Vector3d Pi_new = Pi - U2.col(vertex_index)/V;
//            m_originalMesh->GetPoints()->SetPoint( it->originalID, Pi_new.data() );
//            
//            vertex_difference += (Pi-Pi_new).norm();
//        }
//    }
//    
//    return vertex_difference;
//}
//



vtkIdType UmbrellaWeightedOrder2Smoother::FindVertexConnectivityLocalID( vtkIdType id ) const//uses ids within mesh, compares to originalID
{
    vtkIdType vertex_index=0;
    for( VertexConnectivityArrayType::const_iterator it=m_C.begin(); it!=m_C.end(); it++, vertex_index++ )
    {
        if(it->originalID==id)
            break;
    }    
    
    return vertex_index;
}



void UmbrellaWeightedOrder2Smoother::CalculateEdgeWeightMatrix( SparseMatrixDoubleType& W ) const
{
    const vtkIdType n_pts = m_originalMesh->GetNumberOfPoints();
    
    //Create cotan weight n x n matrix 
    W.resize(n_pts, n_pts);
    
    //for every vertex
    for( VertexIDArrayType::const_iterator vert_it = m_coverVertexIDs->begin();
            vert_it!=m_coverVertexIDs->end(); vert_it++)
    {
        const vtkIdType v1id = (*vert_it);
        
        //get connectivity
        const vtkIdType v1id_in_C = FindVertexConnectivityLocalID( v1id );
        const VertexIDArrayType &v1_neighbors = m_C.at(v1id_in_C).connectedVertices;

        //for every edge find the 2 triangles
        for( VertexIDArrayType::const_iterator neighb_it = v1_neighbors.begin(); neighb_it!=v1_neighbors.end(); neighb_it++ )
        {
            const vtkIdType v2id = (*neighb_it);
            //edge formed by v1id and v2id
            //find triangles by finding the vertices shared by both
            
            //get neighbors of v2id
            const vtkIdType v2id_in_C = FindVertexConnectivityLocalID(v2id);
            const VertexIDArrayType &v2_neighbors = m_C.at(v2id_in_C).connectedVertices;

            std::cout<<"W("<<v1id<<", "<<v2id<<")="<<W.coeff(v1id, v2id)<<std::endl<<std::flush;
            if( W.coeff(v1id, v2id)==0 ) //to avoid recalculating the weight
            {
                //get intersection of  v1_neighbors and v2_neighbors
                VertexIDArrayType common_vertex_ids;
                if( IntersectVectors(v1_neighbors, v2_neighbors, common_vertex_ids) )
                {
                    for( VertexIDArrayType::const_iterator common_v_it = common_vertex_ids.begin();
                            common_v_it != common_vertex_ids.end(); common_v_it++)
                    {
                        const vtkIdType V_common_id = (*common_v_it);

                        const double w = CalculateEdgeWeight(V_common_id, v1id, v2id);
                        W.coeffRef( v1id, v2id ) += w;
                        W.coeffRef( v2id, v1id ) += w;
                    }
                }
            }
                   
        }                
    }    
}






bool UmbrellaWeightedOrder2Smoother::IntersectVectors( const VertexIDArrayType& a, const VertexIDArrayType& b, VertexIDArrayType& c ) const
{
    std::set<VertexIDArrayType::value_type> t1;

    t1.insert<VertexIDArrayType::const_iterator> (a.begin(), a.end());
    
    c.clear();
    c.reserve( a.size() );
    for( VertexIDArrayType::const_iterator b_it = b.begin(); b_it<b.end(); b_it++ )
    {
        std::set<VertexIDArrayType::value_type>::iterator found_id_it = t1.find( (*b_it) );
        if( found_id_it != t1.end() )
            c.push_back( (*found_id_it) );
    }
    
    return c.size()>0;
}
