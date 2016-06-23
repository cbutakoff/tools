/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
//patch smoothing by solving 
//U^2 (v) = 0
//where U is the second order weighted umbrella operator with cotangent weights
//v - each of the patch vertices
//The patch must be a topological disk embedded in 3D
//For details see paper of Liepa, Filling Holes in Meshes
//
//Solution using iterative process from :
//Interactive Multi-Resolution Modeling on Arbitrary Meshes, Leif Kobbelt

#ifndef __UmbrellaWeightedOrder2Smoother_h
#define __UmbrellaWeightedOrder2Smoother_h

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include "HoleFillerDefines.h"



//uncomment to use cotangent weights during smoothing instead of scale-dependent
//#define USE_COTANGENT_WEIGHTS


//updates the input mesh
class UmbrellaWeightedOrder2Smoother
{
public:
    typedef std::vector<vtkIdType> VertexIDArrayType;
    
    void SetInputMesh(vtkPolyData *mesh) {m_originalMesh = mesh;};
    void SetCoverVertexIds(VertexIDArrayType* ids) {m_coverVertexIDs = ids;};
    
    void SetTolerance( double eps ) {m_tolerance = eps; };
    double GetTolerance() {return m_tolerance; };
    
    void SetMaxIter( int i ) {m_maxIter = i; };
    int GetMaxIter() {return m_maxIter; };
    
    void Update();
    
    
    UmbrellaWeightedOrder2Smoother():m_originalMesh(NULL), m_coverVertexIDs(NULL), m_maxIter(50), m_tolerance(1e-6) {};
    
    
protected:
    typedef enum {vcInterior, vcBoundary, vcExterior} VertexClassType;
 
    typedef struct __connectivity
    {
        vtkIdType originalID;
        VertexIDArrayType connectedVertices;
        VertexClassType vertexClass;
    } VertexConnectivityType;    
    
    
    typedef Eigen::Matrix<double, 3, Eigen::Dynamic> MatrixUType;
    
    double TriangleWeightScaleDependent(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2) const;    
    void GetVertexNeighbors(vtkIdType vertexId, VertexIDArrayType& neighbors);
    
    void CalculateConnectivity(); //updates m_C
    void ClassifyVertex(VertexConnectivityType& v);
    
    //Fills U and weights
    void CalculateU( MatrixUType& U, Eigen::SparseMatrix<double>& weights ); 
    //uses U and weights and fills U2
    void CalculateU2( MatrixUType& U2, const MatrixUType& U, const Eigen::SparseMatrix<double>& weights );
    
    double UpdateMeshPoints( const MatrixUType& U2 );
    
    vtkIdType FindVertexConnectivityInfo( vtkIdType id ); //uses ids within mesh, compares to originalID
    
private:
    vtkPolyData *m_originalMesh;
    VertexIDArrayType* m_coverVertexIDs;
    
    double m_tolerance;
    int m_maxIter;

    
    std::vector<VertexConnectivityType> m_C;
};


#endif
