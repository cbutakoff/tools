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



vtkIdType UmbrellaWeightedOrder2Smoother::FindVertexConnectivityLocalID( vtkIdType id ) //uses ids within mesh, compares to originalID
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
        const vtkIdType vertex_ID_in_C = FindVertexConnectivityLocalID(v1id);
        VertexIDArrayType &neighbors = m_C.at(vertex_ID_in_C).connectedVertices;

        //for every edge find the 2 triangles
        
        //for each triangle calculate weight
        
    }
    
    
    //for every triangle calculate cot of every angle and fill the matrix W
//    
//    { //Create a scope for triplet list (to release it afterwards)
//        typedef Eigen::Triplet<double> T;
//        std::vector<T> tripletList;
//
//        VectorType v1, v2, v3;
//        for( HoleCoverType::const_iterator it = m_coverFaces->begin(); it!=m_coverFaces->end(); it++)
//        {
//            const TriangleCellType& tri = *it;
//            for(int j=0; j<3; j++)
//            {
//                const int id1 = tri.id[j];
//                const int id2 = tri.id[(j+1)%3];
//                const int id3 = tri.id[(j+2)%3];
//                m_coverVertices->GetPoint(id1, v1.data());
//                m_coverVertices->GetPoint(id2, v2.data());
//                m_coverVertices->GetPoint(id3, v3.data());
//
//                //angle between v[3],v[1] and v[1],v[2]
////                const VectorType v13 = (v3-v1).normalized();
////                const VectorType v12 = (v2-v1).normalized();
////                const double angle = std::acos(v13.dot(v12));
////                const double w23 = 1/std::tan(angle);
//    //            W.coeffRef(id2,id3) += w23;
//    //            W.coeffRef(id3,id2) += w23;
//                
//#ifdef USE_COTANGENT_WEIGHTS                
//                const double w23 = TriangleWeightCotangent(v1, v2, v3);
//#else
//                const double w23 = TriangleWeightScaleDependent(v1, v2, v3);
//#endif
//
//                tripletList.push_back(T(id2, id3, 23));
//                tripletList.push_back(T(id3, id2, w23));
//            }
//        }
//
//        //fill W from triplets
//        W.setFromTriplets(tripletList.begin(), tripletList.end());
//    }m_C
}