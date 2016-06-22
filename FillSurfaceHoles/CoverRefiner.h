/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
//patch refiner


#ifndef __CoverRefiner_h
#define __CoverRefiner_h

#include <vtkPoints.h>
#include <vtkSmartPointer.h>

#include "HoleFillerDefines.h"


class CoverRefiner //inplace refinement
{
public:
    void SetInputFaces(HoleCoverType *coverFaces) {m_coverFaces = coverFaces;};
    void SetInputVertices(vtkPoints* coverVertices) {m_coverVertices = coverVertices; };
    void SetInputBoundaryIds( const VertexIDArrayType *boundaryIds ) {m_boundaryIds=boundaryIds;};
    
    //this has to be called from exterior
    void InitializeVertexWeights( vtkPolyData* mesh, const VertexIDArrayType *originalBoundaryIds );
    
    void Update();
        
    CoverRefiner():m_coverFaces(NULL), m_boundaryIds(NULL) {};
protected:
    
    void GetVertexNeighbors(vtkPolyData *mesh, vtkIdType vertexId, std::set<vtkIdType>& connectedVertices);
    bool RelaxAllCoverEdges(Eigen::MatrixXi& conn);    
    bool FindConnectedVertices(const EdgeType& edge, EdgeType& intersectingEdge);
    bool RelaxEdgeIfPossible(const EdgeType& edge, const EdgeType& candidateEdge, Eigen::MatrixXi& conn);
    bool CheckForDuplicateTriangles();
    bool IsTriangleSplitRequired(const vtkIdType idVi, const vtkIdType idVj, const vtkIdType idVk, VectorType& Vc, double& Svc);
    
    //check if ptcheck is inside a circle defined by the 3 points
    bool IsPointInCircle(const VectorType& pt0, const VectorType& pt1, const VectorType& pt2, const VectorType& ptcheck) const;
    HoleCoverType::iterator FindTriangleByPointIds( vtkIdType id0, vtkIdType id1, vtkIdType id2);     
private:
    HoleCoverType *m_coverFaces;
    const VertexIDArrayType *m_boundaryIds;
    vtkPoints* m_coverVertices;    

    std::vector<double> m_sigmas;
};


#endif
