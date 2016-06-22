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
    
    void InitializeVertexWeights( const vtkPolyData* mesh, const VertexIDArrayType *originalBoundaryIds );
    
    void Update();
        
    CoverRefiner():m_coverFaces(NULL), m_boundaryIds(NULL) {};
protected:
    
 
private:
    HoleCoverType *m_coverFaces;
    const VertexIDArrayType *m_boundaryIds;
    vtkPoints* m_coverVertices;    

    std::vector<double> m_sigmas;
};


#endif
