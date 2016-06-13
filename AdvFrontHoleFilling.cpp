//The advancing‐front mesh generation method revisited - ReadCube
#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkCell.h>
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>
#include <vtkDataSetAttributes.h>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>

#include <vector>
#include <iostream>
#include <limits>
#include <cmath>

#include "MinHeap/MinHeap.h"


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

//boundary edges have to be ordered. Clockwise when looking from the top, normal pointing outwards
//this coincides with normal orientation of the triangles (the edges are follow corkscrew rule with respect to the normal)
void FillHole(vtkPolyData* mesh, HoleBoundaryType& ordered_boundary);


inline void CrossProduct(const PointType& u, const PointType& v, PointType& cp);
inline double DotProduct(const PointType& u, const PointType& v);

//calculates (u x v).w
inline double MixedProduct(const PointType& u, const PointType& v, const PointType& w);

//create vector from the points: v = p2-p1
inline void MakeVector(const PointType& p1, const PointType& p2, PointType& v);
//same as above but from point ids
inline void MakeVector(vtkPolyData* mesh, vtkIdType p1id, vtkIdType p2id, PointType& v);

//angle between two vectors (arccos)
inline double CalculateAngle(const PointType& v1, const PointType& v2);

//norm of the vector
inline double CalculateVectorNorm(const PointType& v);

//angles[i] is the angle between edge[i] and edge[i+1]
void CalculateHoleAngles(vtkPolyData* mesh, HoleBoundaryType& ordered_boundary, std::vector<double> angles);

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
    
    //print the boundaries (point ids of edges)
//    for(int i=0; i<hole_boundaries.size(); i++)
//    {
//        std::cout<<"Boundary "<<i<<std::endl;
//        for(int j=0; j<hole_boundaries[i].size(); j++)
//        {
//            std::cout<<hole_boundaries[i][j].v0<<" "<<hole_boundaries[i][j].v1<<std::endl;
//        }
//    }
    
    std::cout<<"Found "<<hole_boundaries.size()<<" holes"<<std::endl;

                
    //fill the holes
    for(int i=0; i<hole_boundaries.size(); i++)
    {
        std::cout<<"Filling hole "<<i<<"/"<<hole_boundaries.size()<<std::endl;
        FillHole(mesh, hole_boundaries[i]);
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




//angles[i] is the angle between edge[i] and edge[i+1]
//assumes correct edge orientation
void CalculateHoleAngles(vtkPolyData* mesh, HoleBoundaryType& ordered_boundary, std::vector<double> angles)
{
    for(vtkIdType edgeid = 0; edgeid<ordered_boundary.size(); edgeid++)
    {
        //check if edge[edgeid].v1 is concave or not
        //assuming clockwise orientation looking from the top, corresponds to 
        //triangles having counterclockwise orientation looking from top.
        //
        //if e1, e2 are two edges, n-normal, then
        // (e1 x e2).n < 0 - angle is concave
        // (e1 x e2).n > 0 - angle is convex
        // (e1 x e2) = 0 - angle is pi
        //
        const EdgeType &e1 = ordered_boundary[edgeid];
        const EdgeType &e2 = ordered_boundary[edgeid+1>ordered_boundary.size()?edgeid+1:0];

        PointType v1;
        PointType v2;
        MakeVector(mesh, e1.v0, e1.v1, v1);
        MakeVector(mesh, e2.v0, e2.v1, v2);
        
        //get average normal for the vertex e1.v1;
        PointType n;
        n[0] = (e1.n0[0] + e1.n1[0] + e2.n1[0])/3;
        n[1] = (e1.n0[1] + e1.n1[1] + e2.n1[1])/3;
        n[2] = (e1.n0[2] + e1.n1[2] + e2.n1[2])/3;
        
        //calculate mixed product
        const double mp = MixedProduct(v1,v2,n);
        double angle = CalculateAngle(v1,v2);

        const double pi = boost::math::constants::pi<double>();
        
        if(fabs(mp)<std::numeric_limits<double>::epsilon())
            angle = pi;
        else if( mp < 0 )
            angle = 2*pi - angle;
        
        angles[edgeid] = angle;
    }
}








void CrossProduct(const PointType& u, const PointType& v, PointType& cp)
{
    cp[0] = u[1]*v[2]-u[2]*v[1];
    cp[1] = u[2]*v[0]-u[0]*v[2];
    cp[2] = u[0]*v[1]-u[1]*v[0];
}


double DotProduct(const PointType& u, const PointType& v)
{
    return u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
}

//calculates (u x v).w
double MixedProduct(const PointType& u, const PointType& v, const PointType& w)
{
    PointType cp;
    
    CrossProduct(u,v,cp);
    return DotProduct(cp, w);
}


//create vector from the points: v = p2-p1
void MakeVector(const PointType& p1, const PointType& p2, PointType& v)
{
    v[0] = p2[0]-p1[0];
    v[1] = p2[1]-p1[1];
    v[2] = p2[2]-p1[2];
}


inline void MakeVector(vtkPolyData* mesh, vtkIdType p1id, vtkIdType p2id, PointType& v)
{
    PointType p1;
    PointType p2;
    mesh->GetPoint(p1id,p1);
    mesh->GetPoint(p2id,p2);
    
    MakeVector(p1,p2,v);
}


//angle between two vectors (arccos)
double CalculateAngle(const PointType& v1, const PointType& v2)
{
    const double dp = DotProduct(v1,v2);
    const double v1norm = CalculateVectorNorm(v1);
    const double v2norm = CalculateVectorNorm(v2);
    return std::acos(dp/(v1norm*v2norm));
}


//norm of the vector
double CalculateVectorNorm(const PointType& v)
{
    return std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}




void FillHole(vtkPolyData* mesh, HoleBoundaryType& ordered_boundary)
{
    bool finished=false;
    
    //create front boundary
    HoleBoundaryType front(ordered_boundary);
    
    
    while( !finished )
    {
        std::vector<double> angles;
        CalculateHoleAngles(mesh, front, angles);
    }
}
