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
    void SetInput(const HoleCoverType *cover, const VertexIDArrayType *boundaryIds);
    
    void Update();
    
    std::auto_ptr<HoleCoverType> GetOutput() {return m_smoothCover; }; //the output will be deleted upon class destruction, copy it
    
    UmbrellaWeightedOrder2Smoother():m_originalCover(NULL), m_boundaryIds(NULL) {};
protected:
private:
    const HoleCoverType *m_originalCover;
    const VertexIDArrayType *m_boundaryIds;
    std::auto_ptr<HoleCoverType> m_smoothCover;
};


#endif
