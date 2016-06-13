//The advancing‐front mesh generation method revisited - ReadCube
#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkCell.h>
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <vector>
#include <vtkDataSetAttributes.h>
#include <iostream>


typedef double PointType[3];

typedef struct __edge{
    vtkIdType v0; //vertices
    vtkIdType v1;
    PointType n0; //normals at the vertices
    PointType n1;
} EdgeType;

typedef std::vector<EdgeType> HoleBoundaryType; 
typedef std::vector<HoleBoundaryType> ArrayOfBoundariesType; 
    

//returns the unordered boundary of all the holes
void FindHoles(vtkPolyData *mesh, HoleBoundaryType& unordered_edges);
void SplitHoles(HoleBoundaryType& unordered_edges, ArrayOfBoundariesType& hole_boundaries);

int main(int argc, char **argv)
{
    const char *filename = argv[1];
    
    vtkSmartPointer<vtkPolyDataReader> rdr = vtkSmartPointer<vtkPolyDataReader>::New();
    rdr->SetFileName(filename);
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
    
    for(int i=0; i<hole_boundaries.size(); i++)
    {
        std::cout<<"Boundary "<<i<<std::endl;
        for(int j=0; j<hole_boundaries[i].size(); j++)
        {
            std::cout<<hole_boundaries[i][j].v0<<" "<<hole_boundaries[i][j].v1<<std::endl;
        }
    }
    
    return 0;
}

void FindHoles(vtkPolyData *mesh, HoleBoundaryType& unordered_edges)
{
    using namespace boost::numeric::ublas;    
    
    const vtkIdType nPts = mesh->GetNumberOfPoints();
    coordinate_matrix<short> adj(nPts, nPts); //adjacency matrix, counts number of times an edge appears
    
    unordered_edges.clear();
    
    for(vtkIdType cellid = 0; cellid<mesh->GetNumberOfCells(); cellid++)
    {
        vtkCell* cell = mesh->GetCell(cellid);
        
        //std::cout<<cellid<<": "<<cell->GetPointId(0)<<" "<<cell->GetPointId(1)<<" "<<cell->GetPointId(2)<<std::endl;
        
        for(vtkIdType edgeid = 0; edgeid<cell->GetNumberOfEdges(); edgeid++)
        {
            vtkCell* edge = cell->GetEdge(edgeid);
            const vtkIdType pt0 = edge->GetPointId(0);
            const vtkIdType pt1 = edge->GetPointId(1);
            if(pt0<pt1) //to conserve memory due to symmetry
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
            if(*i2!=0)
            {
                EdgeType edge;
        
                if( *i2 > 0 )
                {
                    edge.v0 = i2.index1();
                    edge.v1 = i2.index2();
                }
                else //orientation was reversed
                {
                    edge.v0 = i2.index2();
                    edge.v1 = i2.index1();
                }
                    
                double n[3];
                mesh->GetPointData()->GetNormals()->GetTuple(edge.v0,n);
                edge.n0[0] = n[0];
                edge.n0[1] = n[1];
                edge.n0[2] = n[2];
                mesh->GetPointData()->GetNormals()->GetTuple(edge.v1,n);
                edge.n1[0] = n[0];
                edge.n1[1] = n[1];
                edge.n1[2] = n[2];
                unordered_edges.push_back(edge);
//                cout << "(" << i2.index1() << "," << i2.index2()
//                     << ":" << *i2 << ")  "<< endl;
            }
        
    }    
}




void SplitHoles(HoleBoundaryType& unordered_edges, ArrayOfBoundariesType& hole_boundaries)
{
    std::vector<bool> visited_edges(unordered_edges.size(), false);
    
    for(vtkIdType edgeid=0; edgeid<unordered_edges.size(); edgeid++)
    {
        if(!visited_edges[edgeid])
        {
            const EdgeType &edgec = unordered_edges[edgeid];

            vtkIdType front = edgec.v1; //front vertex in the queue

            visited_edges[edgeid] = true;

            HoleBoundaryType bdry;
            bdry.push_back(edgec);

//            for(int k=0; k<visited_edges.size(); k++)
//                std::cout<<" | "<<k<<" : "<<visited_edges[k]?1:0;
//            std::cout<<std::endl;

            vtkIdType j=edgeid;
            while (j<unordered_edges.size())
            {
                //std::cout<<"checking edge "<<j<<" - "<<visited_edges[j]<<std::endl;
                if(!visited_edges[j])
                {
                    const EdgeType &edge = unordered_edges[j];
                    if(edge.v0==front) //if this edge back vertex matches current front vertex, add the edge
                    {
                        bdry.push_back(edge);
                        visited_edges[j]=true;
                        front = edge.v1;

                        j=edgeid; //reset counter
                    }
                }
                j++;
            }

            hole_boundaries.push_back(bdry);
        }
    }
}

