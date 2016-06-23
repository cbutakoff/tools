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


#ifndef __UmbrellaWeightedOrder2Smoother_h
#define __UmbrellaWeightedOrder2Smoother_h

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include "HoleFillerDefines.h"


//uncomment to use cotangent weights during smoothing instead of scale-dependent
//#define USE_COTANGENT_WEIGHTS


class UmbrellaWeightedOrder2Smoother
{
public:
    void SetInputFaces(const HoleCoverType *coverFaces) {m_coverFaces = coverFaces;};
    void SetInputVertices(vtkPoints* coverVertices) {m_coverVertices = coverVertices; };
    void SetInputBoundaryIds( const VertexIDArrayType *boundaryIds ) {m_boundaryIds=boundaryIds;};
    
    void Update();
    
    vtkSmartPointer<vtkPoints> GetOutputVertices() {return m_smoothedCoverPoints; }; //the output will be deleted upon class destruction, copy it
    
    UmbrellaWeightedOrder2Smoother():m_coverFaces(NULL), m_boundaryIds(NULL) {};
protected:
    
    void CalculateWeightMatrix();
    void SumSparseMatrixCols( const SparseDoubleMatrixType& m, Eigen::VectorXd& s);
    void FormRightHandSide();
    void AddBoundaryToWeigtMatrix();
    void CreateOutput();
    void TestBoundaryConservation(); //check that the smoothing solution did not modify the boundary
    
    //returns cotangent of the angle between v1v2 and v1v3
    double TriangleWeightCotangent(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2, const Eigen::VectorXd& v3) const;
    double TriangleWeightScaleDependent(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2, const Eigen::VectorXd& v3) const;
private:
    const HoleCoverType *m_coverFaces;
    const VertexIDArrayType *m_boundaryIds;
    vtkPoints* m_coverVertices;
    vtkSmartPointer<vtkPoints> m_smoothedCoverPoints;
    
    
    SparseDoubleMatrixType m_C; //sparse matrix of the system of equations
    Eigen::MatrixXd m_b; //solve Cx = b, where x and b have 3 coordinates (3 columns)
    Eigen::MatrixXd m_x;
};


#endif
