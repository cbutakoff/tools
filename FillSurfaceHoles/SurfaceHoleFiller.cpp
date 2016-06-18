/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
#include "SurfaceHoleFiller.h"


#include <vtkCell.h>
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>
#include <vtkDataSetAttributes.h>
#include <vtkCellData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkDataWriter.h>
#include <vtkWriter.h>
#include <vtkAlgorithm.h>
#include <vtkAdjacentVertexIterator.h>


#include <iostream>
#include <limits>
#include <cmath>
#include <float.h>



#include <stack>

//To create initial cover, by default the triangle areas augmented by dihedral angles are used.
//Uncomment this to return to just triangle areas. Not recommended, becasue sometimes the cover can become crap
//see the paper for explanation
//#define USE_TRIANGLEAREA_WEIGHT



void SurfaceHoleFiller::Update()
{
    vtkSmartPointer<vtkPolyDataNormals> normal_gen = vtkSmartPointer<vtkPolyDataNormals>::New();
    normal_gen->SetInputData(this->m_inputMesh);
    normal_gen->SplittingOff();
    normal_gen->ConsistencyOn();
    normal_gen->Update();

    m_outputMesh = normal_gen->GetOutput();


    HoleBoundaryType unordered_edges;
    FindHoles(m_outputMesh, unordered_edges);

    ArrayOfBoundariesType hole_boundaries;
    SplitHoles(unordered_edges, hole_boundaries);


    std::cout << "Found " << hole_boundaries.size() << " holes" << std::endl;

    
    //fill the holes
    ArrayOfCoversType covers(hole_boundaries.size());

    for (int i = 0; i < hole_boundaries.size(); i++) {
        std::cout << "Filling hole " << i << "/" << hole_boundaries.size() << std::endl;

        //        std::cout<<"Boundary "<<i<<std::endl;
        //        for(int j=0; j<hole_boundaries[i].size(); j++)
        //        {
        //            std::cout<<hole_boundaries[i][j].v0<<" "<<hole_boundaries[i][j].v1<<std::endl;
        //        }

        HoleCoverType cover;

        vtkSmartPointer<vtkPolyData> refinedCover = vtkSmartPointer<vtkPolyData>::New();
        InitialCoverTriangulation(m_outputMesh, hole_boundaries[i], cover);
        
        //updates the mesh inside
        RefineCover(m_outputMesh, hole_boundaries[i], cover);

//            vtkSmartPointer<vtkPolyDataWriter> wr = vtkSmartPointer<vtkPolyDataWriter>::New();
//            wr->SetFileName("closed.vtk");
//            wr->SetInputData(m_outputMesh);
//            wr->Write();

    }

}



//======================================================================
//
//
//      Other functions
//
//
//
//======================================================

void SurfaceHoleFiller::FindHoles(vtkPolyData *mesh, HoleBoundaryType& unordered_edges) const {

    const vtkIdType nPts = mesh->GetNumberOfPoints();
    SparseShortMatrixType adj(nPts, nPts); //adjacency matrix, counts number of times an edge appears

    unordered_edges.clear();

    for (vtkIdType cellid = 0; cellid < mesh->GetNumberOfCells(); cellid++) {
        vtkCell* cell = mesh->GetCell(cellid);

        //std::cout<<cellid<<": "<<cell->GetPointId(0)<<" "<<cell->GetPointId(1)<<" "<<cell->GetPointId(2)<<std::endl;

        for (vtkIdType edgeid = 0; edgeid < cell->GetNumberOfEdges(); edgeid++) {
            vtkCell* edge = cell->GetEdge(edgeid);
            const vtkIdType pt0 = edge->GetPointId(0);
            const vtkIdType pt1 = edge->GetPointId(1);
            if (pt0 < pt1) //to conserve memory due to symmetry
                adj.coeffRef(pt0, pt1) += 1;
            else
                adj.coeffRef(pt1, pt0) -= 1; //to keep track of orientation
        }
    }

    //std::cout << adj << std::endl;
    //iterate over every edge that adj(i,j)==1, these are  the boundaries

    for (vtkIdType k=0; k<adj.outerSize(); ++k) {        
        for (SparseShortMatrixType::InnerIterator it(adj,k); it; ++it)
            if ( it.value() != 0) {
                EdgeType edge;

                if ( it.value() > 0) {
                    edge.v0 = it.row();
                    edge.v1 = it.col();
                } else //orientation was reversed
                {
                    edge.v0 = it.col();
                    edge.v1 = it.row();
                }

                //                double n[3];
                //                mesh->GetPointData()->GetNormals()->GetTuple(edge.v0,n);
                //                edge.n0[0] = n[0];
                //                edge.n0[1] = n[1];
                //                edge.n0[2] = n[2];
                //                mesh->GetPointData()->GetNormals()->GetTuple(edge.v1,n);
                //                edge.n1[0] = n[0];
                //                edge.n1[1] = n[1];
                //                edge.n1[2] = n[2];
                unordered_edges.push_back(edge);
                //                cout << "(" << i2.index1() << "," << i2.index2()
                //                     << ":" << *i2 << ")  "<< endl;
            }

    }
}

void SurfaceHoleFiller::SplitHoles(HoleBoundaryType& unordered_edges, ArrayOfBoundariesType& hole_boundaries) const {
    std::vector<bool> visited_edges(unordered_edges.size(), false);

    for (vtkIdType edgeid = 0; edgeid < unordered_edges.size(); edgeid++) {
        if (!visited_edges[edgeid]) {
            const EdgeType &edgec = unordered_edges[edgeid];

            vtkIdType front = edgec.v1; //front vertex in the queue

            visited_edges[edgeid] = true;

            HoleBoundaryType bdry;
            bdry.push_back(edgec);

            //            for(int k=0; k<visited_edges.size(); k++)
            //                std::cout<<" | "<<k<<" : "<<visited_edges[k]?1:0;
            //            std::cout<<std::endl;

            vtkIdType j = edgeid;
            while (j < unordered_edges.size()) {
                //std::cout<<"checking edge "<<j<<" - "<<visited_edges[j]<<std::endl;
                if (!visited_edges[j]) {
                    const EdgeType &edge = unordered_edges[j];
                    if (edge.v0 == front) //if this edge back vertex matches current front vertex, add the edge
                    {
                        bdry.push_back(edge);
                        visited_edges[j] = true;
                        front = edge.v1;

                        j = edgeid; //reset counter
                    }
                }
                j++;
            }

            hole_boundaries.push_back(bdry);
        }
    }
}

void SurfaceHoleFiller::EdgesToVertexArray(const HoleBoundaryType& ordered_boundary, VertexIDArrayType& vertices) const {
    vertices.clear();

    if (ordered_boundary.empty()) {
        return; //nothing to do
    }

    //    vertices.push_back(ordered_boundary[0].v0);
    //    vertices.push_back(ordered_boundary[0].v1);
    //
    //    for (vtkIdType i = 1; i < ordered_boundary.size()-1; i++) {
    //        vertices.push_back(ordered_boundary[i].v1);
    //    }

    //save in reversed order, for correct triangle orientation
    vertices.push_back(ordered_boundary[ordered_boundary.size() - 2].v1);
    vertices.push_back(ordered_boundary[ordered_boundary.size() - 2].v0);

    for (vtkIdType i = ordered_boundary.size() - 3; i >= 0; i--) {
        vertices.push_back(ordered_boundary[i].v0);
    }



    //vertices.push_back(ordered_boundary[0].v0);
}

double SurfaceHoleFiller::TriangleWeightFunctionArea(const VectorType& u, const VectorType& v, const VectorType& w) const {
    //triangle area   
    return 0.5 * ((v - u).cross(w - u)).norm();
}

void SurfaceHoleFiller::InitialCoverTriangulation(vtkPolyData* mesh, HoleBoundaryType& ordered_boundary, HoleCoverType& cover) const {
    //create an array of vertices, this will create an array ending with vertex 0
    VertexIDArrayType vertices;
    EdgesToVertexArray(ordered_boundary, vertices);

    //print vertices
    //    std::cout<<"Vertices: ";
    //    for(int i=0; i<vertices.size(); i++)
    //    {   
    //        std::cout<<vertices[i]<<" ";
    //    }
    //    std::cout<<std::endl;

    //initialize some variables, consistently with th paper
    //n - number of vertices+1 (vn = v0)
    //number of unique vertices = number of edges
    vtkIdType n = vertices.size();

    //1. Create upper triangular weight matrix W
    using namespace boost::numeric::ublas;

#ifdef USE_TRIANGLEAREA_WEIGHT    
    triangular_matrix<double, upper> W(n, n);
#else //otherwise use angle/area pair. Use complex just for convenience, 
    //but the operations will have to be redefined
    triangular_matrix<AreaAngleMeasureType, upper> W(n, n);
#endif


    //initialize weight matrix
    for (int i = 0; i < n - 1; i++) W(i, i + 1) = 0;

    VectorType vi, vi1, vi2, vm, vk;

    //create storage for the minimum ms
    TriangularIDMatrixType O(n, n);

    //fill O with -1 for debugging
    for (int i = 1; i < n; i++)
        for (int j = i + 1; j < n; j++)
            O(i, j) = -1;


    for (int i = 0; i < n - 2; i++) {
        mesh->GetPoint(vertices[i], vi.data());
        mesh->GetPoint(vertices[i + 1], vi1.data());
        mesh->GetPoint(vertices[i + 2], vi2.data());

#ifdef USE_TRIANGLEAREA_WEIGHT    
        W(i, i + 2) = TriangleWeightFunctionArea(vi, vi1, vi2);
#else
        //find the adjacent triangles on the original mesh
        const vtkIdType va1id = FindThirdVertexId(mesh, vertices[i], vertices[i + 1]);
        const vtkIdType va2id = FindThirdVertexId(mesh, vertices[i + 1], vertices[i + 2]);

        VectorType va1;
        VectorType va2;
        mesh->GetPoint(va1id, va1.data());
        mesh->GetPoint(va2id, va2.data());


        //initialize matrix 0 with the neighbors of the boundary elements
        O(i, i + 1) = va1id;
        O(i + 1, i + 2) = va2id;

        //triangle1: vi, vi1, vi2 -> correct orientation -> vi2, vi1, vi (the boundary orientation corresponds to triangles outside of the hole)
        //triangle2: vi, vi1, va1 - this is in the original mesh, so it has correct orientation
        //triangle3: vi1, vi2, va2 - this is in the original mesh, so it has correct orientation
        //
        //transform the dihedral angle cosines so that it is larger for larger angles (by multiplying by -1)
        //
        //dihedral angle with the triangle adjacent to edge (i, i+1)
        const double dih_angle1 = -CalculateDihedralAngleCos(vi, vi1, vi2, va1, vi1, vi);
        //dihedral angle with the triangle adjacent to edge (i+1, i+2)
        const double dih_angle2 = -CalculateDihedralAngleCos(vi, vi1, vi2, va2, vi2, vi1);

        //        std::cout<<"Angle 1 :"<<vertices[i]<<" "<<vertices[i+1]<<" "<<vertices[i+2]<<" vs "<<
        //                va1id<<" "<<vertices[i+1]<<" "<<vertices[i]<<std::endl;
        //        std::cout<<"Angle 2 :"<<vertices[i]<<" "<<vertices[i+1]<<" "<<vertices[i+2]<<" vs "<<
        //                va2id<<" "<<vertices[i+2]<<" "<<vertices[i+1]<<std::endl;
        //        std::cout<<"Angles: "<<dih_angle1<<" "<<dih_angle2<<std::endl;


        W(i, i + 2) = AreaAngleMeasureType(std::max(dih_angle1, dih_angle2),
                TriangleWeightFunctionArea(vi, vi1, vi2));

        O(i, i + 2) = i + 1;

#endif
    }

    //one more special case edge between first and last point
    O(0, n - 1) = FindThirdVertexId(mesh, vertices[0], vertices[n - 1]);

    //==========================================
    //
    // !!!warning!!!!
    // Notice that for the O(i,i+1) for all i, including O(0, n-1) have the ids 
    // of the triangles in the original mesh, not the ids in the vertices array.
    // Every other element of O is actually an id in the vertices array.
    // I could not come up with a more elegant solution. Therefore there will be 
    // some "if" statements ahead


    //2. Look for minimal triangulation
    vtkIdType j = 2;


    while (j < n - 1) {
        j = j + 1;
        for (vtkIdType i = 0; i < n - j; i++) {
            vtkIdType k = i + j;

            //find smallest W(i,m)+W(m,k)+F(vi,vm,vk) for i<m<k
            vtkIdType m_min = 0;

#ifdef USE_TRIANGLEAREA_WEIGHT    
            double W_min = DBL_MAX;
#else
            AreaAngleMeasureType W_min = AreaAngleMeasureType(DBL_MAX, 0);
#endif

            for (vtkIdType m = i + 1; m < k; m++) {
                mesh->GetPoint(vertices[i], vi.data());
                mesh->GetPoint(vertices[m], vm.data());
                mesh->GetPoint(vertices[k], vk.data());

#ifdef USE_TRIANGLEAREA_WEIGHT    
                double Wik = W(i, m) + W(m, k) + TriangleWeightFunctionArea(vi, vm, vk);
#else

                //calculate dihedral angles of triangles 
                //(imk) vs (i,Oim,m)
                //(imk) vs (m,Omk,k)
                VectorType vOim;
                VectorType vOmk;

                //std::cout<<"Triangle "<<vertices[i]<<" "<<vertices[m]<<" "<<vertices[k]<<std::endl;

                //This is a bit messy
                if (abs(i - m) != 1) //nonconsecutive boundary vertices - O has id of the vertices array
                {
                    //std::cout<<"Oim:"<<vertices[i]<<" "<<vertices[O(i,m)]<<" "<<vertices[m]<<std::endl;
                    mesh->GetPoint(vertices[O(i, m)], vOim.data());
                } else //consecutive boundary vertices - O has the id of the mesh cell array
                {
                    //std::cout<<"Oim:"<<vertices[i]<<" "<<O(i,m)<<" "<<vertices[m]<<std::endl;
                    mesh->GetPoint(O(i, m), vOim.data()); //this vertex is not on the boundary, can't find elegant solution
                }

                if (abs(m - k) != 1) {
                    mesh->GetPoint(vertices[O(m, k)], vOmk.data());
                    //std::cout<<"Omk:"<<vertices[m]<<" "<<vertices[O(m,k)]<<" "<<vertices[k]<<std::endl<<std::flush;
                } else {
                    mesh->GetPoint(O(m, k), vOmk.data()); //this vertex is not on the boundary, can't find elegant solution
                    //std::cout<<"Omk:"<<vertices[m]<<" "<<O(m,k)<<" "<<vertices[k]<<std::endl<<std::flush;
                }

                //calculate dihedral angles from vertex coordinates
                const double dih_angle1 = -CalculateDihedralAngleCos(vi, vm, vk, vi, vOim, vm);
                double dih_angle2 = -CalculateDihedralAngleCos(vi, vm, vk, vm, vOmk, vk);
                //std::cout<<"Angles: "<<dih_angle1<<" "<<dih_angle2<<std::endl;

                if (i == 0 && k == n - 1) //these are consecutive vertices - O - has id within mesh cell array
                {
                    VectorType vOik;
                    //std::cout<<"Oik:"<<O(i,k)<<std::endl<<std::flush;

                    mesh->GetPoint(O(i, k), vOik.data());
                    const double dih_angle3 = -CalculateDihedralAngleCos(vi, vm, vk, vi, vOik, vk);

                    //std::cout<<"Last triangle:"<<vertices[i]<<" "<<O(i,k)<<" "<<vertices[k]<<std::endl<<std::flush;

                    dih_angle2 = std::max(dih_angle2, dih_angle3);
                }

                AreaAngleMeasureType Aimk(std::max(dih_angle1, dih_angle2), TriangleWeightFunctionArea(vi, vm, vk));

                //Calculate the combined weight as a sum of 
                //W(i, m) + W(m, k) + TriangleWeightFunctionArea(vi, vm, vk);  
                AreaAngleMeasureType Wik = SumAreaTriangleMeasures(
                        SumAreaTriangleMeasures(W(i, m), W(m, k)), Aimk);
#endif

#ifdef USE_TRIANGLEAREA_WEIGHT    
                if (Wik < W_min)
#else
                if (AreaTriangleMeasureLess(Wik, W_min))
#endif
                {
                    W_min = Wik;
                    m_min = m;
                }

            }

            W(i, k) = W_min;
            O(i, k) = m_min;
            //std::cout<<"("<<i<<","<<k<<") Best: "<<vertices[i]<<" "<<vertices[m_min]<<" "<<vertices[k]<<std::endl;
        }
    }

    //4. Create the covering mesh
    //vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    //copy the original cells
    //for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); i++)
    //    cells->InsertNextCell(mesh->GetCell(i));

    //create new cells
    PopulateCover(cover, 0, n - 1, O, vertices);

    //for(vtkIdType j=0; j<cover.size(); j++)
    //{
    //    std::cout<<"Verifying: "<<cover[j].id[0]<<" "<<cover[j].id[1]<<" "<<cover[j].id[2]<<std::endl;
    //}


    //Add the cover to the mesh
    //mesh->SetPolys(cells);
    //mesh->GetCellData()->Initialize();
    //mesh->BuildCells();
}

void SurfaceHoleFiller::PopulateCover(HoleCoverType& cover, vtkIdType i, vtkIdType k, 
        const TriangularIDMatrixType& O, const VertexIDArrayType& vertices) const {
    if (i + 2 == k) {
        //cover->InsertNextCell(3);
        //cover->InsertCellPoint(vertices[i]); //readjusted the orientation wrt the paper
        //cover->InsertCellPoint(vertices[i + 1]);
        //cover->InsertCellPoint(vertices[k]);
        TriangleCellType triangle;
        triangle.id[0] = vertices[i];
        triangle.id[1] = vertices[i + 1];
        triangle.id[2] = vertices[k];
        cover.push_back(triangle);
        //std::cout<<"Adding: "<<vertices[i]<<" "<<vertices[i+1]<<" "<<vertices[k]<<std::endl;
    } else {
        const vtkIdType o = O(i, k);
        if (o != i + 1)
            PopulateCover(cover, i, o, O, vertices);

        //cover->InsertNextCell(3);
        //cover->InsertCellPoint(vertices[i]); //readjusted the orientation wrt the paper
        //cover->InsertCellPoint(vertices[o]);
        //cover->InsertCellPoint(vertices[k]);
        TriangleCellType triangle;
        triangle.id[0] = vertices[i];
        triangle.id[1] = vertices[o];
        triangle.id[2] = vertices[k];
        cover.push_back(triangle);
        //std::cout<<"Adding: "<<vertices[i]<<" "<<vertices[o]<<" "<<vertices[k]<<std::endl;

        if (o != k - 1)
            PopulateCover(cover, o, k, O, vertices);

    }
}



//make sure the vertices have correct order\
//https://en.wikipedia.org/wiki/Dihedral_angle

double SurfaceHoleFiller::CalculateDihedralAngleCos(const VectorType& v1, const VectorType& v2, const VectorType& v3,
        const VectorType& u1, const VectorType& u2, const VectorType& u3) const {
    const VectorType n1 = ((v3 - v2).cross(v1 - v2)).normalized();
    const VectorType n2 = ((u3 - u2).cross(u1 - u2)).normalized();

    //    std::cout<<"V1: "<<v1(0)<<" "<<v1(1)<<" "<<v1(2)<<" "<<std::endl;
    //    std::cout<<"V2: "<<v2(0)<<" "<<v2(1)<<" "<<v2(2)<<" "<<std::endl;
    //    std::cout<<"V3: "<<v3(0)<<" "<<v3(1)<<" "<<v3(2)<<" "<<std::endl;
    //    std::cout<<"U1: "<<u1(0)<<" "<<u1(1)<<" "<<u1(2)<<" "<<std::endl;
    //    std::cout<<"U2: "<<u2(0)<<" "<<u2(1)<<" "<<u2(2)<<" "<<std::endl;
    //    std::cout<<"U3: "<<u3(0)<<" "<<u3(1)<<" "<<u3(2)<<" "<<std::endl;
    //    std::cout<<"n1: "<<n1(0)<<" "<<n1(1)<<" "<<n1(2)<<" "<<std::endl;
    //    std::cout<<"n2: "<<n2(0)<<" "<<n2(1)<<" "<<n2(2)<<" "<<std::endl;

    return n1.dot(n2);
}


//this will be exected only for boundaries, so there is only 1 triangle with the edge (p1,p2)

vtkIdType SurfaceHoleFiller::FindThirdVertexId(vtkPolyData* mesh, vtkIdType p1, vtkIdType p2) const {
    vtkSmartPointer<vtkIdList> n = vtkSmartPointer<vtkIdList>::New();
    mesh->GetPointCells(p1, n);

    vtkIdType cell_id;
    vtkIdList* ptids;

    for (vtkIdType i = 0; i < n->GetNumberOfIds(); i++) {
        cell_id = n->GetId(i);
        ptids = mesh->GetCell(cell_id)->GetPointIds();

        for (vtkIdType j = 0; j < ptids->GetNumberOfIds(); j++) {
            if (p2 == ptids->GetId(j)) //found the cell sharing 2 points
            {
                j = ptids->GetNumberOfIds();
                i = n->GetNumberOfIds();
            }
        }

    }

    //cell_id has the identified cell id;
    //ptids - has ids of all the cell points
    vtkIdType thirdid = VTK_ID_MAX;
    for (vtkIdType i = 0; i < ptids->GetNumberOfIds(); i++) {
        thirdid = ptids->GetId(i);
        if (thirdid != p1 && thirdid != p2)
            break;
    }

    return thirdid;
}





//does m1+m2, as described in the paper
//(a,b)+(c,d) = (max(a,c),b+d)

AreaAngleMeasureType SurfaceHoleFiller::SumAreaTriangleMeasures(const AreaAngleMeasureType& m1,
        const AreaAngleMeasureType& m2) const {
    return AreaAngleMeasureType(std::max(std::real(m1), std::real(m2)),
            std::imag(m1) + std::imag(m2));
}

//checks if m1 < m2, lexicographically
//(a,b)<(c,d) iff (a<c or (a==c and b<d))

bool SurfaceHoleFiller::AreaTriangleMeasureLess(const AreaAngleMeasureType& m1,
        const AreaAngleMeasureType& m2) const {
    bool less = false;
    if (std::real(m1) < std::real(m2))
        less = true;
    else if (std::imag(m1) < std::imag(m2))
        less = true;

    return less;
}




void SurfaceHoleFiller::SaveIsolatedCover(const HoleCoverType& localCover, vtkPoints * coverVertices, const char* filename) const {
    vtkSmartPointer<vtkCellArray> coverCells = vtkSmartPointer<vtkCellArray>::New();
    for( HoleCoverType::const_iterator it = localCover.begin(); it!=localCover.end(); ++it )
    {
        coverCells->InsertNextCell(3);
        coverCells->InsertCellPoint((*it).id[0]);
        coverCells->InsertCellPoint((*it).id[1]);
        coverCells->InsertCellPoint((*it).id[2]);
    }
    
    vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
    pd->SetPoints(coverVertices);
    pd->SetPolys(coverCells);
    
    vtkSmartPointer<vtkPolyDataWriter> wr = vtkSmartPointer<vtkPolyDataWriter>::New();
    wr->SetFileName(filename);
    wr->SetInputData(pd);
    wr->Write();
}

//static int cover_id = 0; //for debugging output
bool SurfaceHoleFiller::RelaxEdgeIfPossible(const EdgeType& edge, const EdgeType& candidateEdge,
        vtkPoints* coverVertices, HoleCoverType& localCover, SparseShortMatrixType& conn) const {
    VectorType edgeV0, edgeV1, candV0, candV1; 
    
    coverVertices->GetPoint(edge.v0, edgeV0.data());
    coverVertices->GetPoint(edge.v1, edgeV1.data());
    coverVertices->GetPoint(candidateEdge.v0, candV0.data());
    coverVertices->GetPoint(candidateEdge.v1, candV1.data());
    
    const double edge_length = (edgeV0-edgeV1).squaredNorm();
    const double cand_edge_length = (candV0-candV1).squaredNorm();
    
    bool swap_performed = false;
    
    if(cand_edge_length<edge_length) //check the circumference criterion
    {
        if( IsPointInCircle(edgeV0, edgeV1, candV0, candV1) ) //do the swap
        {
            //do the swap here. Erase old triangles
            HoleCoverType::iterator it1;
            HoleCoverType::iterator it2;
            it1 = FindTriangleByPointIds(localCover, edge.v0, edge.v1, candidateEdge.v0);
            it2 = FindTriangleByPointIds(localCover, edge.v0, edge.v1, candidateEdge.v1);
            
            TriangleCellType newtriangle1;
            TriangleCellType newtriangle2;
            
            //identify vertex order and create new triangles
            //forward order 
            bool triangles_created = false;
            for(int i=0; i<3; i++)
            {
//                std::cout<<"Triangle: "<<edge.v0<<" "<<edge.v1<<" "<<candidateEdge.v0<<std::endl;
//                std::cout<<"Checking: "<<(*it1).id[i]<<" "<<(*it1).id[(i+1)%3]<<" "<<(*it1).id[(i+2)%3]<<std::endl;
                if( ( (*it1).id[i]==edge.v0 && (*it1).id[(i+1)%3]==edge.v1 && (*it1).id[(i+2)%3]==candidateEdge.v0) )
                {
                    newtriangle1.id[0] = candidateEdge.v0;
                    newtriangle1.id[1] = edge.v0;
                    newtriangle1.id[2] = candidateEdge.v1;
                    newtriangle2.id[0] = candidateEdge.v1;
                    newtriangle2.id[1] = edge.v1;
                    newtriangle2.id[2] = candidateEdge.v0;
                    triangles_created = true;
                    break;
                }
            }
            
            if(!triangles_created)
            {
                for(int i=0; i<3; i++)
                {
//                    std::cout<<"Triangle: "<<edge.v0<<" "<<edge.v1<<" "<<candidateEdge.v0<<std::endl;
//                    std::cout<<"Checking: "<<(*it1).id[(i+2)%3]<<" "<<(*it1).id[(i+1)%3]<<" "<<(*it1).id[i]<<std::endl;
                    if( ( (*it1).id[(i+2)%3]==edge.v0 && (*it1).id[(i+1)%3]==edge.v1 && (*it1).id[i]==candidateEdge.v0) )
                    {
                        newtriangle1.id[0] = candidateEdge.v1;
                        newtriangle1.id[1] = edge.v0;
                        newtriangle1.id[2] = candidateEdge.v0;
                        newtriangle2.id[0] = candidateEdge.v0;
                        newtriangle2.id[1] = edge.v1;
                        newtriangle2.id[2] = candidateEdge.v1;
                        triangles_created = true;
                        break;
                    }
                }
            }
            
            if(triangles_created)
            {
                //create triangles
                localCover.erase(it1); //erase the second trianle
                localCover.erase(it2);
                
                //add the new triangles
                localCover.push_back(newtriangle1);
                localCover.push_back(newtriangle2);
                
//                std::cout<<"Added triangles"<<std::endl;
//                std::cout<<"T1: "<<newtriangle1.id[0]<<" "<<newtriangle1.id[1]<<" "<<newtriangle1.id[2]<<std::endl;
//                std::cout<<"T2: "<<newtriangle2.id[0]<<" "<<newtriangle2.id[1]<<" "<<newtriangle2.id[2]<<std::endl;
//                std::cout<<"edge: "<<edge.v0<<" "<<edge.v1<<std::endl;
//                std::cout<<"cand: "<<candidateEdge.v0<<" "<<candidateEdge.v1<<std::endl;
                //update connectivity matrix
                conn.coeffRef( std::min(edge.v0, edge.v1), std::max(edge.v0, edge.v1) ) = 0;
                conn.coeffRef( std::min(candidateEdge.v0, candidateEdge.v1), std::max(candidateEdge.v0, candidateEdge.v1) ) = 2;
                
                swap_performed = true;
            }
        }
    }
    
    return swap_performed;
}

bool SurfaceHoleFiller::RelaxAllCoverEdges(HoleCoverType& localCover, 
        vtkPoints * coverVertices, SparseShortMatrixType& conn) const {

    
    std::stack<EdgeType> EdgeStack;
    
    for (int k=0; k<conn.outerSize(); ++k)
        for (SparseShortMatrixType::InnerIterator it(conn,k); it; ++it)  {
            if ( it.value() == 2) { //only interior edges (i.e. 2 cover triangles share it)
                EdgeType e;
                e.v0 = it.row();
                e.v1 = it.col();
                EdgeStack.push(e);
            }
            //std::cout<<(*i2)<<" ";
        }
    
    
    //      
    //   for each element in the queue.
    //      take the intersecting edge vertices
    //      fit a circle to the original 2 vertices plus one of the intersecting ones
    //      check if the 2nd intersecting vertex is inside. If it is swap the edge only if it becomes shorter
    //          compare the length of the old edge and the proposed edge ->swap ->update conn matrix
    //          
    
    
    bool swap_performed = false;
    
    while(!EdgeStack.empty())
    {
        EdgeType edge = EdgeStack.top();
        EdgeStack.pop();
        
        //Using the upper triangular connectivity matrix find the 2 vertices that are connected to the edge
        EdgeType candidateEdge;
        if( FindConnectedVertices(coverVertices, localCover, edge, candidateEdge) )
        {
            
            //verify if swap is needed, first check the edge length criterion, then circumference
            //std::cout<<"edge "<<edge.v0<<" "<<edge.v1<<std::endl;
            //std::cout<<"cand "<<candidateEdge.v0<<" "<<candidateEdge.v1<<std::endl;
            swap_performed = RelaxEdgeIfPossible(edge, candidateEdge, coverVertices, localCover, conn);
            
        }
    }
    
    return swap_performed;
}

bool SurfaceHoleFiller::IsTriangleSplitRequired(vtkPoints* coverVertices, const std::vector<double>& sigmas, 
        const vtkIdType idVi, const vtkIdType idVj, const vtkIdType idVk, VectorType& Vc, double& Svc) const {
    VectorType Vi, Vj, Vk;
    coverVertices->GetPoint( idVi, Vi.data() );
    coverVertices->GetPoint( idVj, Vj.data() );
    coverVertices->GetPoint( idVk, Vk.data() );
    
    Vc = (Vi+Vj+Vk)/3;
    Svc = ( sigmas[idVi] + sigmas[idVj] + sigmas[idVk] )/3;
    
    //check for all m=i,j,k 2*||vc-vm||^2 > max(s(vc)^2,s(vm)^2)
    const double Svc2 = Svc*Svc;
    const double Svi2 = sigmas[idVi]*sigmas[idVi];
    const double Svj2 = sigmas[idVj]*sigmas[idVj];
    const double Svk2 = sigmas[idVk]*sigmas[idVk];
    const double dist2_Vc_Vi = 2*(Vc-Vi).squaredNorm();
    const double dist2_Vc_Vj = 2*(Vc-Vj).squaredNorm();
    const double dist2_Vc_Vk = 2*(Vc-Vk).squaredNorm();
    
//    std::cout<<"Splitting triangle ("<<idVi<<" "<<idVj<<" "<<idVk<<std::endl;
//    std::cout<<"S(vi,vj,vk)="<<sigmas[idVi]<<" "<<sigmas[idVj]<<" "<<sigmas[idVk]<<std::endl;
//    std::cout<<"S(vc)="<<Svc<<std::endl;
//    std::cout<<"||Vc-Vi||="<<std::sqrt(dist2_Vc_Vi)<<std::endl;
//    std::cout<<"||Vc-Vj||="<<std::sqrt(dist2_Vc_Vj)<<std::endl;
//    std::cout<<"||Vc-Vk||="<<std::sqrt(dist2_Vc_Vk)<<std::endl;
    
    const bool split_required = ( dist2_Vc_Vi > std::max(Svc2, Svi2) &&
            dist2_Vc_Vj > std::max(Svc2, Svj2) && 
            dist2_Vc_Vk > std::max(Svc2, Svk2) );
    
//    std::cout<<"Split ?: "<<split_required<<std::endl;
    return split_required;
}

void SurfaceHoleFiller::RefineCover(vtkPolyData* mesh, const HoleBoundaryType& ordered_boundary, 
        const HoleCoverType& cover) const {

    //create vertex storage
    VertexIDArrayType boundaryVertexIDs;
    EdgesToVertexArray(ordered_boundary, boundaryVertexIDs);

    //Relabel the ids of the cover to ids from 0 to n_vertices
    //the correspondence with the original is given by
    //i -> boundaryVertexIDs[i]
    HoleCoverType localCover;
    for( HoleCoverType::const_iterator it = cover.begin(); it!=cover.end(); ++it )
    {
        TriangleCellType cell; 
        
        //for each id, find its position within boundaryVertexIDs
        for(int j=0; j<3; j++)
        {
            vtkIdType i;
            for(i=0; i<boundaryVertexIDs.size(); i++)
                if( (*it).id[j] == boundaryVertexIDs[i] ) break;
            
            cell.id[j] = i;
        }
        
        localCover.push_back(cell);
    }
    
    
    //copy boundary vertices to point storage
    vtkSmartPointer<vtkPoints> coverVertices = vtkSmartPointer<vtkPoints>::New();
    for(vtkIdType i=0; i<boundaryVertexIDs.size(); i++)
    {
        coverVertices->InsertNextPoint(mesh->GetPoint( boundaryVertexIDs[i] ));
    }

    //Save the cover for debugging
//    char name[100];
//    sprintf(name,"initial_cover%d.vtk",cover_id++);
//    SaveIsolatedCover(localCover, coverVertices, name);
            
    
    //build upper triangular! vertex connectivity matrix for the cover
    SparseShortMatrixType conn(coverVertices->GetNumberOfPoints(), coverVertices->GetNumberOfPoints());

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
        const vtkIdType originalVertexID = boundaryVertexIDs[vertexId];
    
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

        sigmas[vertexId] = sigma/vertexNeighbors.size();
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
                localCover.push_back(tri2);
                localCover.push_back(tri3);
                

                //Update matrices
                sigmas.push_back(Svc);
                
                //new point added - need resize
                conn.conservativeResize(coverVertices->GetNumberOfPoints(), coverVertices->GetNumberOfPoints());
                for(int i=0; i<3; i++)
                {
                    conn.coeffRef(std::min(tri1.id[i], tri1.id[(i+1)%3]), std::max(tri1.id[i], tri1.id[(i+1)%3])) = 2;
                    conn.coeffRef(std::min(tri2.id[i], tri2.id[(i+1)%3]), std::max(tri2.id[i], tri2.id[(i+1)%3])) = 2;
                    conn.coeffRef(std::min(tri3.id[i], tri3.id[(i+1)%3]), std::max(tri3.id[i], tri3.id[(i+1)%3])) = 2;
                }

//                std::cout<<conn<<std::endl;

                //relax edges (vi,vj), (vi,vk), (vj,vk)
                EdgeType edge; edge.v0 = idVi; edge.v1 = idVj; 
                EdgeType candidateEdge;
                if( FindConnectedVertices(coverVertices, localCover, edge, candidateEdge) )
                    RelaxEdgeIfPossible(edge, candidateEdge, coverVertices, localCover, conn);
                
                TriangleSplitted = true;
                
//                SaveIsolatedCover(localCover, coverVertices, "refined.vtk");
            }

        }
        
        
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
        while( RelaxAllCoverEdges(localCover, coverVertices, conn) ) {};
        
    }

    

    
    
    //    //Create edge storage. Put only interior edges (no boundary)
//    for (SparseIDMatrixType::iterator1 it1 = conn.begin1(); it1 != conn.end1(); it1++) {
//        for (SparseIDMatrixType::iterator2 it2 = it1.begin(); it2 != it1.end(); it2++) {
//            if (*it2 == 2) {
//                EdgeType edge;
//                edge.v0 = it2.index1();
//                edge.v1 = it2.index2();
//            }
//        }
//    }    

    //get the cover with original ids
    //the correspondence with the original is given by
    //i -> boundaryVertexIDs[i]
    
//    SaveIsolatedCover(localCover, coverVertices, "refined.vtk");
    
    
//    std::cout<<"Copying final refined cover"<<std::endl;

    vtkCellArray* mesh_cells = m_outputMesh->GetPolys();
    vtkPoints* mesh_points = m_outputMesh->GetPoints();
    
    const vtkIdType nOriginalVertices = m_outputMesh->GetNumberOfPoints();
   
    //insert new points, skip old ones, they are in the beginning
    for( vtkIdType ptid = boundaryVertexIDs.size(); ptid<coverVertices->GetNumberOfPoints(); ptid++)
    {
        const vtkIdType id = mesh_points->InsertNextPoint( coverVertices->GetPoint(ptid) );
        boundaryVertexIDs.push_back(id);
    }
    
    for( HoleCoverType::const_iterator it = localCover.begin(); it!=localCover.end(); ++it )
    {
        //TriangleCellType cell; 
        
        //for each id, find its position within boundaryVertexIDs
        mesh_cells->InsertNextCell(3);
        for(int j=0; j<3; j++)
        {
            const vtkIdType local_id = (*it).id[j];
            
            const vtkIdType targetId = boundaryVertexIDs[ local_id ];

            mesh_cells->InsertCellPoint( targetId );
            
//            std::cout<<"Id "<<local_id<<" converted to "<<targetId<<std::endl;
        }
        
        //refinedCover.push_back(cell);
    }    
    
    m_outputMesh->GetCellData()->Initialize();
    m_outputMesh->GetPointData()->Initialize();
    m_outputMesh->BuildCells();
    
//    vtkSmartPointer<vtkPolyDataWriter> wr = vtkSmartPointer<vtkPolyDataWriter>::New();
//    wr->SetFileName("integrated.vtk");
//    wr->SetInputData(m_outputMesh);
//    wr->Write();
}


void SurfaceHoleFiller::GetVertexNeighbors(vtkPolyData *mesh, vtkIdType vertexId, 
        std::set<vtkIdType>& connectedVertices) const {

    connectedVertices.clear();
    
    //get all cells that vertex 'seed' is a part of 
    vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
    mesh->GetPointCells(vertexId, cellIdList);

    //loop through all the cells that use the seed point 
    for (vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++) {

        vtkCell* cell = mesh->GetCell(cellIdList->GetId(i));

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
}


//conn - upper triangular connectivity matrix
bool SurfaceHoleFiller::FindConnectedVertices(vtkPoints* vertices, 
        const HoleCoverType& localCover, const EdgeType& edge, 
        EdgeType& intersectingEdge) const {

    std::vector<vtkIdType> ids;
    bool found_edge = false;
    
//    std::cout<<"Edge orthogonal to ("<<edge.v0<<" "<<edge.v1<<")"<<std::endl;
    
    //iterate through the triangles and find the 2 sharing the edge
    for( HoleCoverType::const_iterator it = localCover.begin(); it!=localCover.end(); it++ )
    {
        int mask[] = {0,0,0}; //0 - vertex not used, 1 - vertex used
        for(int i=0; i<3; i++)
            if( (*it).id[i]==edge.v0 )
            {
                mask[i] = 1;
                for(int j=0; j<3; j++)
                    if( (*it).id[j]==edge.v1 )
                    {
                        mask[j]=1;
                        
                        //triangle found
                        for(int k=0; k<3; k++)
                            if(mask[k]==0)
                            {
                                ids.push_back( (*it).id[k] );
                                k=4; //continue to the next triangle
                                j=4;
                                i=4;
                                found_edge = true; 
                            }
                    }
            }
    }

    if (found_edge)
    {
        intersectingEdge.v0 = ids[0];
        intersectingEdge.v1 = ids[1];
//        std::cout<<"("<<ids.size()<<") Triangles adjacent to ("<<intersectingEdge.v0<<" "<<intersectingEdge.v1<<")"<<std::endl;
    }
    else
    {
        intersectingEdge.v0 = -1;
        intersectingEdge.v1 = -1;
//        std::cout<<"("<<ids.size()<<") No alternative edge"<<std::endl;
    }
    
    return found_edge;
}




bool SurfaceHoleFiller::IsPointInCircle(const VectorType& pt0, 
        const VectorType& pt1, const VectorType& pt2, 
        const VectorType& ptcheck) const {
    //https://en.wikipedia.org/wiki/Circumscribed_circle

    typedef Eigen::Vector2d Vector2DType;

    //map points to 2 dimensions    
    Eigen::Matrix<double, 2, 3> T3Dto2D;
    VectorType v10 = (pt1-pt0).normalized();
    VectorType v  = pt2-pt0;
    VectorType n = v.cross( v10 );
    VectorType v20 = v10.cross(n).normalized();
    
    
    T3Dto2D.row(0) = v10;
    T3Dto2D.row(1) = v20;
    
    Vector2DType A;
    A.fill(0); // = T3Dto2D * (pt0-pt0);
    
    Eigen::Matrix<double, 2, 1> B = T3Dto2D * (pt1-pt0);
    Eigen::Matrix<double, 2, 1> C = T3Dto2D * (pt2-pt0);
    Eigen::Matrix<double, 2, 1> v1 = T3Dto2D * (ptcheck-pt0);
    
    Eigen::Matrix4d M;
//    if det(M)<0 - inside
//    if det(M)=0 - on the circle
//    if det(M)>0 - outside
    M(0,0) = v1.squaredNorm();
    M(1,0) = A.squaredNorm();
    M(2,0) = B.col(0).squaredNorm();
    M(3,0) = C.col(0).squaredNorm();
    
//    std::cout<<"B: "<<B<<std::endl;
//    std::cout<<"T3Dto2D * (pt1-pt0): "<<T3Dto2D * (pt1-pt0)<<std::endl;
//    std::cout<<"C: "<<C<<std::endl;
//    std::cout<<"v1: "<<v1<<std::endl;
    
    M(0,1) = v1[0]; 
    M(1,1) = A[0];
    M(2,1) = B[0];
    M(3,1) = C[0];
    
    M(0,2) = v1[1];
    M(1,2) = A[1];
    M(2,2) = B[1];
    M(3,2) = C[1];

    M(0,3) = 1;
    M(1,3) = 1;
    M(2,3) = 1;
    M(3,3) = 1;
    
//    std::cout<<"3DPoints 1: "<<pt0<<" | "<<pt1<<" | "<<pt2<<" | "<<ptcheck<<std::endl;
//    std::cout<<"M: "<<M<<std::endl;
    
    
//    M = [v*v' p(1) p(2) 1; ...
//     A*A' A(1) A(2) 1; ...
//     B*B' B(1) B(2) 1; ...
//     C*C' C(1) C(2) 1];
    
    const double det = M.determinant();

//    std::cout<<"Det: "<<det<<std::endl;
//    exit(0);

    
    if(det > 0)
        return false;
    else
        return true;    
}





HoleCoverType::iterator SurfaceHoleFiller::FindTriangleByPointIds(HoleCoverType& localCover, 
        vtkIdType id0, vtkIdType id1, vtkIdType id2) const {   
    HoleCoverType::iterator retval;
    
    std::vector<vtkIdType> argumentIDs;
    argumentIDs.push_back(id0);
    argumentIDs.push_back(id1);
    argumentIDs.push_back(id2);
    std::sort(argumentIDs.begin(), argumentIDs.end());

    std::vector<vtkIdType> IDs;

//    std::cout<<"To find: "<<argumentIDs[0]<<" "<<argumentIDs[1]<<" "<<argumentIDs[2]<<std::endl;
    
    for(retval = localCover.begin(); retval!=localCover.end(); retval++)
    {
        IDs.clear();
        IDs.push_back((*retval).id[0]);
        IDs.push_back((*retval).id[1]);
        IDs.push_back((*retval).id[2]);
        std::sort(IDs.begin(), IDs.end());

//        std::cout<<"Comparing to: "<<IDs[0]<<" "<<IDs[1]<<" "<<IDs[2]<<std::endl;

        
        std::vector<vtkIdType> v_intersection;
        std::set_intersection(argumentIDs.begin(), argumentIDs.end(),
                              IDs.begin(), IDs.end(),
                              std::back_inserter(v_intersection));         
        
        if(v_intersection.size()==3) //bingo
            break;
    }
    
    return retval;
}
