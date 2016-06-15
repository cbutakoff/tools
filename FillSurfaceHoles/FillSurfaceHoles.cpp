//Remesh-like hole filling
//based on Eurographics Symposium on Geometry Processing(2003), Filling Holes in Meshes, Peter Liepa
//and Filling Gaps in the Boundary of a Polyhedron, Gill Barequet, Micha Sharir
#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkCell.h>
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>
#include <vtkDataSetAttributes.h>
#include <vtkCellData.h>
#include <vtkPolyDataWriter.h>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <Eigen/Dense>

#include <vector>
#include <iostream>
#include <limits>
#include <cmath>
#include <float.h>
#include <vtkDataWriter.h>
#include <vtkWriter.h>
#include <vtkAlgorithm.h>


//#define USE_TRIANGLEAREA_WEIGHT
typedef std::complex<double> AreaAngleMeasureType;

typedef Eigen::Vector3d VectorType;
typedef boost::numeric::ublas::triangular_matrix<vtkIdType, boost::numeric::ublas::upper> TriangularIDMatrixType;

typedef struct __edge {
    vtkIdType v0; //vertices
    vtkIdType v1;
} EdgeType;

typedef std::vector<EdgeType> HoleBoundaryType;
typedef std::vector<HoleBoundaryType> ArrayOfBoundariesType;

typedef std::vector<vtkIdType> VertexArrayType;

//returns the unordered boundary of all the holes
void FindHoles(vtkPolyData *mesh, HoleBoundaryType& unordered_edges);
void SplitHoles(HoleBoundaryType& unordered_edges, ArrayOfBoundariesType& hole_boundaries);

//the last vertex is equal to the first vertex to simulate closed loop.
//so the number of returned vertices is actually the number of edges +1
void EdgesToVertexArray(const HoleBoundaryType& ordered_boundary, VertexArrayType& vertices);


void InitialTriangulation(vtkPolyData* mesh, HoleBoundaryType& ordered_boundary);

//for explanation of O, see the paper. This is recursive function Trace
void PopulateCover(vtkCellArray* cover, vtkIdType i, vtkIdType k, const TriangularIDMatrixType& O, const VertexArrayType& vertices);

inline double TriangleWeightFunctionArea(const VectorType& u, const VectorType& v, const VectorType& w);


//calculates cosine of the dihedral angle. Vertices must have correct orientation
//v1, v2, v3 - triangle 1
//u1, u2, u3 - triangle 2
inline double CalculateDihedralAngleCos(const VectorType& v1, const VectorType& v2, const VectorType& v3,
                                        const VectorType& u1, const VectorType& u2, const VectorType& u3);



inline vtkIdType FindThirdVertexId(vtkPolyData* mesh, vtkIdType p1, vtkIdType p2);


//does m1+m2, as described in the paper
//(a,b)+(c,d) = (max(a,c),b+d)
inline AreaAngleMeasureType SumAreaTriangleMeasures(const AreaAngleMeasureType& m1, 
        const AreaAngleMeasureType& m2);

//checks if m1 < m2, lexicographically
//(a,b)<(c,d) iff (a<c or (a==c and b<d))
inline bool AreaTriangleMeasureLess(const AreaAngleMeasureType& m1, 
        const AreaAngleMeasureType& m2);

//======================================================================
//
//
//      Main function
//
//
//
//======================================================

int main(int argc, char **argv) {
    std::cout << "Hole filling based on Eurographics Symposium on Geometry Processing(2003), Filling Holes in Meshes, Peter Liepa"<<std::endl;
    std::cout << "usage: FillSurfaceHoles -i mesh.vtk -o mesh.vtk -method [Barequet|Liepa_flat|Liepa_smooth]" << std::endl;
    std::cout << "version 1.0" << std::endl;
    std::cout << "-------------------------------------------------------------------" << std::endl;

    if (argc < 3) {
        return -1;
    }

    char *in_filename=NULL;
    char *out_filename=NULL;

    enum AlgorithmType {Liepa_flat, Liepa_smooth};
    AlgorithmType filling_algorithm = Liepa_flat;
    
    for (int c = 1; c < argc; c++) {
        if (strcmp(argv[c], "-i") == 0) {
            in_filename = argv[++c];
        } else if (strcmp(argv[c], "-o") == 0) {
            out_filename = argv[++c];
        } else if (strcmp(argv[c], "-method") == 0) {
            const char* methodtype = argv[++c];
            if(strcmp(methodtype, "Liepa_flat")==0) 
                filling_algorithm = Liepa_flat;
            if(strcmp(methodtype, "Liepa_smooth")==0) 
                filling_algorithm = Liepa_smooth;
            else
                std::cout<<"Unknown filling algorithm, using default."<<std::endl;
        }
    }

    std::cout<<"Input file: "<<in_filename<<std::endl;
    std::cout<<"Output file: "<<out_filename<<std::endl;
    std::cout << "-------------------------------------------------------------------" << std::endl;

    vtkSmartPointer<vtkPolyDataReader> rdr = vtkSmartPointer<vtkPolyDataReader>::New();
    rdr->SetFileName(in_filename);
    
    if (!rdr->IsFilePolyData())
    {
        std::cout<<"Input file is not plydata. Verify your file. Aborting"<<std::endl;
        return -1;
    }

    rdr->Update();

    vtkSmartPointer<vtkPolyDataNormals> normal_gen = vtkSmartPointer<vtkPolyDataNormals>::New();
    normal_gen->SetInputData(rdr->GetOutput());
    normal_gen->SplittingOff();
    normal_gen->ConsistencyOn();
    normal_gen->Update();

    vtkPolyData *mesh = normal_gen->GetOutput();


    HoleBoundaryType unordered_edges;
    FindHoles(mesh, unordered_edges);

    ArrayOfBoundariesType hole_boundaries;
    SplitHoles(unordered_edges, hole_boundaries);


    std::cout << "Found " << hole_boundaries.size() << " holes" << std::endl;


    //fill the holes
    for (int i = 0; i < hole_boundaries.size(); i++) {
        std::cout << "Filling hole " << i << "/" << hole_boundaries.size() << std::endl;

        std::cout<<"Boundary "<<i<<std::endl;
        for(int j=0; j<hole_boundaries[i].size(); j++)
        {
            std::cout<<hole_boundaries[i][j].v0<<" "<<hole_boundaries[i][j].v1<<std::endl;
        }

        InitialTriangulation(mesh, hole_boundaries[i]);
        
        if( filling_algorithm == Liepa_smooth )
        {
            //add smoothing remeshing here
        }

    }


    //save the result
    vtkSmartPointer<vtkPolyDataWriter> wr = vtkSmartPointer<vtkPolyDataWriter>::New();
    wr->SetFileName(out_filename);
    wr->SetInputData(mesh);
    wr->Write();

    return 0;
}



//======================================================================
//
//
//      Other functions
//
//
//
//======================================================

void FindHoles(vtkPolyData *mesh, HoleBoundaryType& unordered_edges) {
    using namespace boost::numeric::ublas;

    const vtkIdType nPts = mesh->GetNumberOfPoints();
    coordinate_matrix<short> adj(nPts, nPts); //adjacency matrix, counts number of times an edge appears

    unordered_edges.clear();

    for (vtkIdType cellid = 0; cellid < mesh->GetNumberOfCells(); cellid++) {
        vtkCell* cell = mesh->GetCell(cellid);

        //std::cout<<cellid<<": "<<cell->GetPointId(0)<<" "<<cell->GetPointId(1)<<" "<<cell->GetPointId(2)<<std::endl;

        for (vtkIdType edgeid = 0; edgeid < cell->GetNumberOfEdges(); edgeid++) {
            vtkCell* edge = cell->GetEdge(edgeid);
            const vtkIdType pt0 = edge->GetPointId(0);
            const vtkIdType pt1 = edge->GetPointId(1);
            if (pt0 < pt1) //to conserve memory due to symmetry
                adj(pt0, pt1) += 1;
            else
                adj(pt1, pt0) -= 1; //to keep track of orientation
        }
    }

    //std::cout << adj << std::endl;
    //iterate over every edge that adj(i,j)==1, these are  the boundaries
    typedef coordinate_matrix<short>::iterator1 i1_t;
    typedef coordinate_matrix<short>::iterator2 i2_t;

    for (i1_t i1 = adj.begin1(); i1 != adj.end1(); ++i1) {
        for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
            if (*i2 != 0) {
                EdgeType edge;

                if (*i2 > 0) {
                    edge.v0 = i2.index1();
                    edge.v1 = i2.index2();
                } else //orientation was reversed
                {
                    edge.v0 = i2.index2();
                    edge.v1 = i2.index1();
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

void SplitHoles(HoleBoundaryType& unordered_edges, ArrayOfBoundariesType& hole_boundaries) {
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

void EdgesToVertexArray(const HoleBoundaryType& ordered_boundary, VertexArrayType& vertices) {
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
    vertices.push_back(ordered_boundary[ordered_boundary.size()-2].v1);
    vertices.push_back(ordered_boundary[ordered_boundary.size()-2].v0);

    for (vtkIdType i = ordered_boundary.size()-3; i>=0; i--) {
        vertices.push_back(ordered_boundary[i].v0);
    }


    
    //vertices.push_back(ordered_boundary[0].v0);
}

inline double TriangleWeightFunctionArea(const VectorType& u, const VectorType& v, const VectorType& w) {
    //triangle area   
    return 0.5 * ((v - u).cross(w - u)).norm();
}

void InitialTriangulation(vtkPolyData* mesh, HoleBoundaryType& ordered_boundary) {
    //create an array of vertices, this will create an array ending with vertex 0
    VertexArrayType vertices;
    EdgesToVertexArray(ordered_boundary, vertices);

    //print vertices
    std::cout<<"Vertices: ";
    for(int i=0; i<vertices.size(); i++)
    {   
        std::cout<<vertices[i]<<" ";
    }
    std::cout<<std::endl;
    
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

    VectorType vi(3);
    VectorType vi1(3);
    VectorType vi2(3);
    VectorType vm(3);
    VectorType vk(3);

    //create storage for the minimum ms
    TriangularIDMatrixType O(n, n);
    
    //fill O with -1 for debugging
    for(int i=1; i<n; i++)
        for(int j=i+1; j<n; j++)
            O(i,j) = -1;
    
    
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

        VectorType va1(3);
        VectorType va2(3);
        mesh->GetPoint(va1id, va1.data());
        mesh->GetPoint(va2id, va2.data());

        
        //initialize matrix 0 with the neighbors of the boundary elements
        O(i, i+1) = va1id;
        O(i+1, i+2) = va2id;

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
        std::cout<<"Angle 1 :"<<vertices[i]<<" "<<vertices[i+1]<<" "<<vertices[i+2]<<" vs "<<
                va1id<<" "<<vertices[i+1]<<" "<<vertices[i]<<std::endl;
        std::cout<<"Angle 2 :"<<vertices[i]<<" "<<vertices[i+1]<<" "<<vertices[i+2]<<" vs "<<
                va2id<<" "<<vertices[i+2]<<" "<<vertices[i+1]<<std::endl;
        std::cout<<"Angles: "<<dih_angle1<<" "<<dih_angle2<<std::endl;
        
        
        W(i, i + 2) = AreaAngleMeasureType( std::max(dih_angle1, dih_angle2),
                TriangleWeightFunctionArea(vi, vi1, vi2) );
        
        O(i, i+2) = i+1;

#endif
    }
    
    //one more special case edge between first and last point
    O(0, n-1) = FindThirdVertexId(mesh, vertices[0], vertices[n-1]);

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
            AreaAngleMeasureType W_min = DBL_MAX;
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
                VectorType vOim(3);
                VectorType vOmk(3);
                std::cout<<"Triangle "<<vertices[i]<<" "<<vertices[m]<<" "<<vertices[k]<<std::endl;

                
                if(abs(i-m)!=1)
                {
                    std::cout<<"Oim:"<<vertices[i]<<" "<<vertices[O(i,m)]<<" "<<vertices[m]<<std::endl;
                    mesh->GetPoint(vertices[O(i,m)], vOim.data());
                }
                else
                {
                    std::cout<<"Oim:"<<vertices[i]<<" "<<O(i,m)<<" "<<vertices[m]<<std::endl;
                    mesh->GetPoint(O(i,m), vOim.data()); //this vertex is not on the boundary, can't find elegant solution
                }

                if(abs(m-k)!=1)
                {
                    mesh->GetPoint(vertices[O(m,k)], vOmk.data());
                    std::cout<<"Omk:"<<vertices[m]<<" "<<vertices[O(m,k)]<<" "<<vertices[k]<<std::endl<<std::flush;
                }
                else
                {
                    mesh->GetPoint(O(m,k), vOmk.data()); //this vertex is not on the boundary, can't find elegant solution
                    std::cout<<"Omk:"<<vertices[m]<<" "<<O(m,k)<<" "<<vertices[k]<<std::endl<<std::flush;
                }
                
                const double dih_angle1 = -CalculateDihedralAngleCos(vi, vm, vk, vi, vOim, vm);
                double dih_angle2 = -CalculateDihedralAngleCos(vi, vm, vk, vm, vOmk, vk);
                std::cout<<"Angles: "<<dih_angle1<<" "<<dih_angle2<<std::endl;

                if(i==0 && k==n-1)
                {
                    VectorType vOik(3);
                    std::cout<<"Oik:"<<O(i,k)<<std::endl<<std::flush;
                    mesh->GetPoint(O(i,k), vOik.data());
                    const double dih_angle3 = -CalculateDihedralAngleCos(vi, vm, vk, vi, vOik, vk);
                    std::cout<<"Last triangle:"<<vertices[i]<<" "<<O(i,k)<<" "<<vertices[k]<<std::endl<<std::flush;
                    dih_angle2 = std::max( dih_angle2, dih_angle3 );
                }
                
                AreaAngleMeasureType Aimk(std::max(dih_angle1, dih_angle2), TriangleWeightFunctionArea(vi, vm, vk));
                
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
            std::cout<<"("<<i<<","<<k<<") Best: "<<vertices[i]<<" "<<vertices[m_min]<<" "<<vertices[k]<<std::endl;
        }
    }

    //4. Create the covering mesh
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    //copy the original cells
    for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); i++)
        cells->InsertNextCell(mesh->GetCell(i));

    //create new cells
    PopulateCover(cells, 0, n - 1, O, vertices);

    //Add the cover to the mesh
    mesh->SetPolys(cells);
    mesh->GetCellData()->Initialize();
    mesh->BuildCells();
}

void PopulateCover(vtkCellArray* cover, vtkIdType i, vtkIdType k, const TriangularIDMatrixType& O, const VertexArrayType& vertices) {
    if (i + 2 == k) {
        cover->InsertNextCell(3);
        cover->InsertCellPoint(vertices[i]); //readjusted the orientation wrt the paper
        cover->InsertCellPoint(vertices[i + 1]);
        cover->InsertCellPoint(vertices[k]);
    } else {
        const vtkIdType o = O(i, k);
        if (o != i + 1)
            PopulateCover(cover, i, o, O, vertices);

        cover->InsertNextCell(3);
        cover->InsertCellPoint(vertices[i]); //readjusted the orientation wrt the paper
        cover->InsertCellPoint(vertices[o]);
        cover->InsertCellPoint(vertices[k]);

        if (o != k - 1)
            PopulateCover(cover, o, k, O, vertices);

    }
}



//make sure the vertices have correct order\
//https://en.wikipedia.org/wiki/Dihedral_angle
inline double CalculateDihedralAngleCos(const VectorType& v1, const VectorType& v2, const VectorType& v3,
                                        const VectorType& u1, const VectorType& u2, const VectorType& u3)
{
    const VectorType n1 = ((v3-v2).cross(v1-v2)).normalized();
    const VectorType n2 = ((u3-u2).cross(u1-u2)).normalized();
    
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
inline vtkIdType FindThirdVertexId(vtkPolyData* mesh, vtkIdType p1, vtkIdType p2)
{
    vtkSmartPointer<vtkIdList> n = vtkSmartPointer<vtkIdList>::New();
    mesh->GetPointCells(p1, n);
    
    vtkIdType cell_id;
    vtkIdList* ptids;
    
    for(vtkIdType i=0; i<n->GetNumberOfIds(); i++)
    {
        cell_id = n->GetId(i);
        ptids = mesh->GetCell( cell_id )->GetPointIds();
        
        for( vtkIdType j=0; j<ptids->GetNumberOfIds(); j++)
        {
            if(p2 == ptids->GetId(j)) //found the cell sharing 2 points
            {
                j=ptids->GetNumberOfIds();
                i=n->GetNumberOfIds();
            }
        }
        
    }
    
    //cell_id has the identified cell id;
    //ptids - has ids of all the cell points
    vtkIdType thirdid = VTK_ID_MAX;
    for(vtkIdType i=0; i<ptids->GetNumberOfIds(); i++)
    {
        thirdid = ptids->GetId(i);
        if(thirdid!=p1 && thirdid!=p2)
            break;
    }
    
    return thirdid;
}





//does m1+m2, as described in the paper
//(a,b)+(c,d) = (max(a,c),b+d)
inline AreaAngleMeasureType SumAreaTriangleMeasures(const AreaAngleMeasureType& m1, 
        const AreaAngleMeasureType& m2)
{
    return AreaAngleMeasureType( std::max( std::real(m1), std::real(m2) ), 
            std::imag(m1) + std::imag(m2) );
}

//checks if m1 < m2, lexicographically
//(a,b)<(c,d) iff (a<c or (a==c and b<d))
inline bool AreaTriangleMeasureLess(const AreaAngleMeasureType& m1, 
        const AreaAngleMeasureType& m2)
{
    bool less = false;
    if( std::real(m1) < std::real(m2) )
        less = true;
    else if( std::imag(m1) < std::imag(m2) )
        less = true;
    
    return less;
}
