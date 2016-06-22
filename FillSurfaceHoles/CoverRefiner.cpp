/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#include "CoverRefiner.h"
#include "HoleFillerDefines.h"


void CoverRefiner::Update()
{
    //SparseIDMatrixType conn(coverVertices->GetNumberOfPoints(), coverVertices->GetNumberOfPoints());
    SparseIDMatrixType conn(10000, 10000);
    
    //std::cout<<"Cover size: "<<localCover.size()<<std::endl;
    for (HoleCoverType::const_iterator it = localCover.begin(); it != localCover.end(); ++it) {
        const vtkIdType id0 = (*it).id[0];
        const vtkIdType id1 = (*it).id[1];
        const vtkIdType id2 = (*it).id[2];
        
        conn.coeffRef(std::min(id0, id1), std::max(id0, id1)) += 1;
        conn.coeffRef(std::min(id1, id2), std::max(id1, id2)) += 1;
        conn.coeffRef(std::min(id0, id2), std::max(id0, id2)) += 1;
        
        //        std::cout<<"Adding ("<<id0<<", "<<id1<<")"<<std::endl;
        //        std::cout<<"Adding ("<<id1<<", "<<id2<<")"<<std::endl;
        //        std::cout<<"Adding ("<<id0<<", "<<id2<<")"<<std::endl;
        
        
        //std::cout<<"Adding "<<std::min(id0, id1)<<", "<< std::max(id0, id1)<<std::endl;
    }
    //create storage for weights
    //this will be synchronized with coverVertices
    std::vector<double> sigmas(coverVertices->GetNumberOfPoints()); 
    
   
    //------------------------------------------
    //
    //  Step 1:    
    //   for each vertex on hole boundary calculate s(vi) = average adjacent edge lengths, !!do not consider cover edges!!   
    //
    //
    //--------------------------------------------------------------
    VectorType pt0, pt1;
    for( vtkIdType vertexId=0; vertexId<coverVertices->GetNumberOfPoints(); vertexId++ )
    {
        //get the point id in the original mesh
        const vtkIdType originalVertexID = boundaryVertexIDs.at(vertexId);
        
        //find vertex neighbors
        std::set<vtkIdType> vertexNeighbors;
        GetVertexNeighbors(mesh, originalVertexID, vertexNeighbors);
        
        //        std::cout<<"Vertex "<<originalVertexID<<" connected to ";
        //        for(std::set<vtkIdType>::const_iterator it=vertexNeighbors.begin(); it!=vertexNeighbors.end(); it++)
        //            std::cout<<(*it)<<" ";
        //        std::cout<<std::endl;                             
        double sigma = 0;
        
        mesh->GetPoint(originalVertexID, pt0.data());
        
        for(std::set<vtkIdType>::const_iterator it=vertexNeighbors.begin(); it!=vertexNeighbors.end(); it++)
        {
            mesh->GetPoint((*it), pt1.data());
            sigma += (pt1-pt0).norm();
        }
        
        sigmas.at(vertexId) = sigma/vertexNeighbors.size();
    }
    
    //    for(int i=0;i<sigmas.size();i++)
    //        std::cout<<sigmas[i]<<" "<<std::endl;
    
  
    
    //------------------------------------------
    //
    //  Step 1.5:    
    //   relax all edges of the cover
    //
    //
    //--------------------------------------------------
    while( RelaxAllCoverEdges(localCover, coverVertices, conn) ) {};

    
        
    //------------------------------------------
    //
    //  Step 2:    
    //   
    //    
    //2: for each (vi,vj,vk) of the cover 
    //   get centroid vc
    //   s(vc) = (s(vi)+s(vj)+s(vk))/3
    //   (originally) check for all m=i,j,k sqrt(2)*||vc-vm||>max(s(vc),s(vm))
    //   check for all m=i,j,k 2*||vc-vm||^2 > max(s(vc)^2,s(vm)^2)
    //          replace (vi,vj,vk) with (vc,vj,vk), (vi,vc,vk), (vi,vj,vc)
    //          relax edges (vi,vj), (vi,vk), (vj,vk)
    //          split_created = true
    //   if !split_created end
    //
    //--------------------------------------------------
    
    while ( true ) 
    {
        bool TriangleSplitted = false; //algorithm stops when no more splits are required
        
        for(HoleCoverType::iterator coverIt = localCover.begin(); coverIt!=localCover.end(); ++coverIt)    
        {
            
            //calculate centroid
            const vtkIdType idVi = (*coverIt).id[0];
            const vtkIdType idVj = (*coverIt).id[1];
            const vtkIdType idVk = (*coverIt).id[2];
            
            VectorType Vc;
            double Svc;
            
            //            std::cout<<"Splitting verifying: "<<idVi<<" "<<idVj<<" "<<idVk<<std::endl;
            if( IsTriangleSplitRequired(coverVertices, sigmas, idVi, idVj, idVk, Vc, Svc) )
            {  //create new triangles
                //erase old triangle
                
                //                std::cout<<"Splitting confirmed: "<<idVi<<" "<<idVj<<" "<<idVk<<std::endl;
                
                localCover.erase(coverIt);
                //add new vertex
                const vtkIdType idVc = coverVertices->InsertNextPoint(Vc.data());
                
                
                //replace (vi,vj,vk) with (vc,vj,vk), (vi,vc,vk), (vi,vj,vc)
                TriangleCellType tri1; tri1.id[0]=idVc; tri1.id[1]=idVj; tri1.id[2]=idVk;
                TriangleCellType tri2; tri2.id[0]=idVi; tri2.id[1]=idVc; tri2.id[2]=idVk;
                TriangleCellType tri3; tri3.id[0]=idVi; tri3.id[1]=idVj; tri3.id[2]=idVc;
                //relax edges (vi,vj), (vi,vk), (vj,vk)
                //              TriangleSplitted = true;
                localCover.push_back(tri1);
                if( CheckForDuplicateTriangles(localCover) )
                {
                    for( HoleCoverType::const_iterator it1=localCover.begin(); it1!=localCover.end(); it1++  )
                    {                
                        std::cout<<(*it1).id[0]<<" "<<(*it1).id[1]<<" "<<(*it1).id[2]<<" "<<std::endl;
                    }
                }
                localCover.push_back(tri2);
                if( CheckForDuplicateTriangles(localCover) )
                {
                    for( HoleCoverType::const_iterator it1=localCover.begin(); it1!=localCover.end(); it1++  )
                    {                
                        std::cout<<(*it1).id[0]<<" "<<(*it1).id[1]<<" "<<(*it1).id[2]<<" "<<std::endl;
                    }
                }
                localCover.push_back(tri3);
                if( CheckForDuplicateTriangles(localCover) )
                {
                    for( HoleCoverType::const_iterator it1=localCover.begin(); it1!=localCover.end(); it1++  )
                    {                
                        std::cout<<(*it1).id[0]<<" "<<(*it1).id[1]<<" "<<(*it1).id[2]<<" "<<std::endl;
                    }
                }
                
                
                //Update matrices
                sigmas.push_back(Svc);
                
                //new point added - need resize
                //conn.conservativeResize(coverVertices->GetNumberOfPoints(), coverVertices->GetNumberOfPoints());
                std::cout<<"conn matrix size: "<<conn.rows()<<", "<<conn.cols()<<std::endl;
                std::cout<<"max id to store: "<<std::max(tri1.id[0], std::max(tri1.id[1], std::max(tri1.id[3], std::max(
                        tri2.id[0], std::max( tri2.id[1], std::max( tri2.id[2], std::max(
                        tri3.id[0], std::max( tri3.id[1], tri3.id[2])         ))) 
                        ) ) ) )                         <<std::endl;
                for(int i=0; i<3; i++)
                {
                    conn.coeffRef(std::min(tri1.id[i], tri1.id[(i+1)%3]), std::max(tri1.id[i], tri1.id[(i+1)%3])) = 2;
                    conn.coeffRef(std::min(tri2.id[i], tri2.id[(i+1)%3]), std::max(tri2.id[i], tri2.id[(i+1)%3])) = 2;
                    conn.coeffRef(std::min(tri3.id[i], tri3.id[(i+1)%3]), std::max(tri3.id[i], tri3.id[(i+1)%3])) = 2;
                }
                std::cout<<"Stored"<<std::endl;
                //                std::cout<<conn<<std::endl;
                
                //relax edges (vi,vj), (vi,vk), (vj,vk)
                EdgeType edge; edge.v0 = idVi; edge.v1 = idVj; 
                EdgeType candidateEdge;
                
                std::cout<<"relaxing after split"<<std::endl;
                if( FindConnectedVertices(coverVertices, localCover, edge, candidateEdge) )
                    RelaxEdgeIfPossible(edge, candidateEdge, coverVertices, localCover, conn);
                std::cout<<"relaxed after split"<<std::endl;
                
                TriangleSplitted = true;
                
                //                SaveIsolatedCover(localCover, coverVertices, "refined.vtk");
            }
            
        }
        
        
        std::cout<<"Cover size: "<<localCover.size()<<std::endl;
        std::cout<<"Sigmas size: "<<sigmas.size()<<std::endl;
        std::cout<<"Connectivity: "<<conn.rows()<<", "<<conn.cols()<<", "<<conn.size()<<std::endl;
        
        //------------------------------------------
        //
        //  Step 3: If no splits were performed - finish
        //   
        if(!TriangleSplitted) break;
        
        
        
        
        //------------------------------------------
        //
        //  Step 4:    
        //   
        //  Relax all cover edges
        std::cout<<"relaxing started"<<std::endl;
        while( RelaxAllCoverEdges(localCover, coverVertices, conn) ) {};
        std::cout<<"relaxing complete"<<std::endl;
        
    }    
}