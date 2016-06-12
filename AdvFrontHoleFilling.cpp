//The advancing‐front mesh generation method revisited - ReadCube
#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkCell.h>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <vector>


typedef double PointType[3];

typedef struct __edge{
    vtkIdType v0; //vertices
    vtkIdType v1;
    PointType n0; //normals at the vertices
    PointType n1;
} EdgeType;

typedef std::vector<EdgeType> HoleBoundaryType; 
    
//returns the unordered boundary of all the holes
void FindHoles(vtkPolyData *mesh, HoleBoundaryType& boundary);

int main(int argc, char **argv)
{
    const char *filename = argv[1];
    
    vtkSmartPointer<vtkPolyDataReader> rdr = vtkSmartPointer<vtkPolyDataReader>::New();
    rdr->SetFileName(filename);
    rdr->Update();
    
    vtkPolyData *mesh = rdr->GetOutput();
    
    HoleBoundaryType boundary; 
    FindHoles(mesh, boundary);
    
    return 0;
}

void FindHoles(vtkPolyData *mesh, HoleBoundaryType& boundary)
{
    using namespace boost::numeric::ublas;    
    
    const vtkIdType nPts = mesh->GetNumberOfPoints();
    coordinate_matrix<short> adj(nPts, nPts); //adjacency matrix, counts number of times an edge appears
    
    boundary.clear();
    
    for(vtkIdType cellid = 0; cellid<mesh->GetNumberOfCells(); cellid++)
    {
        vtkCell* cell = mesh->GetCell(cellid);
        
        for(vtkIdType edgeid = 0; edgeid<cell->GetNumberOfEdges(); edgeid++)
        {
            vtkCell* edge = cell->GetEdge(edgeid);
            const vtkIdType pt0 = edge->GetPointId(0);
            const vtkIdType pt1 = edge->GetPointId(1);
            if(pt0<pt1) //to conserve memory due to symmetry
                adj(pt0, pt1) += 1;
            else
                adj(pt1, pt0) += 1;
        }
    }
    
    //std::cout << adj << std::endl;
    //iterate over every edge that adj(i,j)==1, these are  the boundaries
    typedef coordinate_matrix<short>::iterator1 i1_t;
    typedef coordinate_matrix<short>::iterator2 i2_t;

    for (i1_t i1 = adj.begin1(); i1 != adj.end1(); ++i1) {
        for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
            if(*i2==1)
            {
                EdgeType edge;
                edge.v0 = i2.index1();
                edge.v1 = i2.index2();
                edge.n0[0] = 0;
                edge.n0[1] = 0;
                edge.n0[2] = 0;
                edge.n1[0] = 0;
                edge.n1[1] = 0;
                edge.n1[2] = 0;
                boundary.push_back(edge);
                //cout << "(" << i2.index1() << "," << i2.index2()
                //     << ":" << *i2 << ")  "<< endl;
            }
        
    }    
}


