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


#define USE_TRIANGLEAREA_WEIGHT


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

inline double TriangleWeightFunction(const VectorType& u, const VectorType& v, const VectorType& w);


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
        //        std::cout<<"Boundary "<<i<<std::endl;
        //        for(int j=0; j<hole_boundaries[i].size(); j++)
        //        {
        //            std::cout<<hole_boundaries[i][j].v0<<" "<<hole_boundaries[i][j].v1<<std::endl;
        //        }

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

    vertices.push_back(ordered_boundary[0].v0);
    vertices.push_back(ordered_boundary[0].v1);

    for (vtkIdType i = 1; i < ordered_boundary.size(); i++) {
        vertices.push_back(ordered_boundary[i].v1);
    }

    vertices.push_back(ordered_boundary[0].v0);
}

inline double TriangleWeightFunction(const VectorType& u, const VectorType& v, const VectorType& w) {
#ifdef USE_TRIANGLEAREA_WEIGHT
    //triangle area   
    return 0.5 * ((v - u).cross(w - u)).norm();
#endif
}

void InitialTriangulation(vtkPolyData* mesh, HoleBoundaryType& ordered_boundary) {
    //create an arra of vertices
    VertexArrayType vertices;
    EdgesToVertexArray(ordered_boundary, vertices);

    //initialize some variables, consistently with th paper
    //n - number of vertices+1 (vn = v0)
    //number of unique vertices = number of edges
    vtkIdType n = vertices.size();

    //1. Create upper triangular weight matrix W
    using namespace boost::numeric::ublas;
    triangular_matrix<double, upper> W(n, n);



    //initialize weight matrix
    for (int i = 0; i < n - 1; i++) W(i, i + 1) = 0;

    VectorType vi(3);
    VectorType vi1(3);
    VectorType vi2(3);
    VectorType vm(3);
    VectorType vk(3);

    for (int i = 0; i < n - 2; i++) {
        mesh->GetPoint(vertices[i], vi.data());
        mesh->GetPoint(vertices[i + 1], vi1.data());
        mesh->GetPoint(vertices[i + 2], vi2.data());
        W(i, i + 3) = TriangleWeightFunction(vi, vi1, vi2);
    }

    //2. Look for minimal triangulation
    vtkIdType j = 2;


    //create storage for the minimum ms
    TriangularIDMatrixType O(n, n);

    while (j < n - 1) {
        j = j + 1;
        for (vtkIdType i = 0; i < n - j; i++) {
            vtkIdType k = i + j;

            //find smallest W(i,m)+W(m,k)+F(vi,vm,vk) for i<m<k
            vtkIdType m_min = 0;
            double W_min = DBL_MAX;
            for (vtkIdType m = i + 1; m < k; m++) {
                mesh->GetPoint(vertices[i], vi.data());
                mesh->GetPoint(vertices[m], vm.data());
                mesh->GetPoint(vertices[k], vk.data());
                double Wik = W(i, m) + W(m, k) + TriangleWeightFunction(vi, vm, vk);

                if (Wik < W_min) {
                    W_min = Wik;
                    m_min = m;
                }
            }
            W(i, k) = W_min;
            O(i, k) = m_min;
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
        cover->InsertCellPoint(vertices[k]); //readjusted the orientation wrt the paper
        cover->InsertCellPoint(vertices[i + 1]);
        cover->InsertCellPoint(vertices[i]);
    } else {
        const vtkIdType o = O(i, k);
        if (o != i + 1)
            PopulateCover(cover, i, o, O, vertices);

        cover->InsertNextCell(3);
        cover->InsertCellPoint(vertices[k]); //readjusted the orientation wrt the paper
        cover->InsertCellPoint(vertices[o]);
        cover->InsertCellPoint(vertices[i]);

        if (o != k - 1)
            PopulateCover(cover, o, k, O, vertices);

    }
}



