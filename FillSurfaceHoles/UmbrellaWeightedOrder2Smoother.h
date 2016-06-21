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


class UmbrellaWeightedOrder2Smoother
{
public:
    void SetInputFaces(const HoleCoverType *coverFaces) {m_coverFaces = coverFaces;};
    void SetInputVertices(vtkPoints* coverVertices) {m_coverVertices = coverVertices; };
    void SetInputBOundaryIds( const VertexIDArrayType *boundaryIds ) {m_boundaryIds=boundaryIds;};
    
    void Update();
    
    vtkSmartPointer<vtkPoints> GetOutputVertices() {return m_smoothedCoverPoints; }; //the output will be deleted upon class destruction, copy it
    
    UmbrellaWeightedOrder2Smoother():m_coverFaces(NULL), m_boundaryIds(NULL) {};
protected:
    
    void CalculateWeightMatrix();
    
private:
    const HoleCoverType *m_coverFaces;
    const VertexIDArrayType *m_boundaryIds;
    vtkPoints* m_coverVertices;
    vtkSmartPointer<vtkPoints> m_smoothedCoverPoints;
    
    SparseDoubleMatrixType m_C; //sparse matrix of the system of equations
};


#endif
