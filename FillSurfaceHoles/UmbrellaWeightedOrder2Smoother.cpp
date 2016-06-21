/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#include "UmbrellaWeightedOrder2Smoother.h"
#include <Eigen/Sparse>



void UmbrellaWeightedOrder2Smoother::CalculateWeightMatrix()
{
    const vtkIdType n_pts = m_coverVertices->GetNumberOfPoints();
    const vtkIdType n_faces = m_coverFaces->size();
    
    //Create cotan weight n x n matrix 
    SparseDoubleMatrixType W(n_pts, n_pts);
    
 
}


void UmbrellaWeightedOrder2Smoother::Update()
{
    if (this->m_coverFaces == NULL )
    {
        std::cout<<"UmbrellaWeightedOrder2Smoother:Error - cover not set. Aborting."<<std::endl;
        return;
    }

    if (this->m_boundaryIds == NULL )
    {
        std::cout<<"UmbrellaWeightedOrder2Smoother:Error - boundary array not set. Aborting."<<std::endl;
        return;
    }

    
    CalculateWeightMatrix();
    
//    C = ComputeWeightMatrix(v.T, f.T);
//            
//    C[b,:] = 0;
//    for i in range(b.size):
//        C[b[i],b[i]] = 1;
//    
//    n = v.shape[0]
//    Rx = np.zeros(n); 
//    Rx[b] = v[b,0];
//    Ry = np.zeros(n); 
//    Ry[b] = v[b,1];
//    Rz = np.zeros(n); 
//    Rz[b] = v[b,2];
//    
//    result = np.zeros( (Rx.size, 3) )
//    result[:,0] = linalg_sp.spsolve(C, Rx); #x
//    result[:,1] = linalg_sp.spsolve(C, Ry); #y
//    result[:,2] = linalg_sp.spsolve(C, Rz); #y
//
//    return result
}

