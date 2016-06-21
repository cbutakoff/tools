/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#include "UmbrellaWeightedOrder2Smoother.h"
#include <Eigen/Sparse>

void UmbrellaWeightedOrder2Smoother::SetInput(const HoleCoverType* cover, const VertexIDArrayType* boundaryIds)
{
    this->m_originalCover = cover;
    this->m_boundaryIds = boundaryIds;
}


void UmbrellaWeightedOrder2Smoother::Update()
{
    
}

