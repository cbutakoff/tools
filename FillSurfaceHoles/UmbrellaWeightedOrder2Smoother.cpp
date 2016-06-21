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
    
    //for every triangle calculate cot of every angle and fill the matrix W
    VectorType v1, v2, v3;
    for( HoleCoverType::const_iterator it = m_coverFaces->begin(); it!=m_coverFaces->end(); it++)
    {
        const TriangleCellType& tri = *it;
        for(int j=0; j<3; j++)
        {
            const int id1 = j;
            const int id2 = (j+1)%3;
            const int id3 = (j+2)%3;
            m_coverVertices->GetPoint(id1, v1.data());
            m_coverVertices->GetPoint(id2, v2.data());
            m_coverVertices->GetPoint(id3, v3.data());
            
            //angle between v[3],v[1] and v[1],v[2]
            const VectorType v31 = (v1-v3).normalized();
            const VectorType v12 = (v2-v1).normalized();
            const double angle = std::acos(v32.dot(v12));
            const double w23 = std::atan(angle);
        }
    }
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

