/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#include "UmbrellaWeightedOrder2Smoother.h"
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <vtkSmartPointer.h>


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
            const int id1 = tri.id[j];
            const int id2 = tri.id[(j+1)%3];
            const int id3 = tri.id[(j+2)%3];
            m_coverVertices->GetPoint(id1, v1.data());
            m_coverVertices->GetPoint(id2, v2.data());
            m_coverVertices->GetPoint(id3, v3.data());
            
            //angle between v[3],v[1] and v[1],v[2]
            const VectorType v13 = (v3-v1).normalized();
            const VectorType v12 = (v2-v1).normalized();
            const double angle = std::acos(v13.dot(v12));
            const double w23 = 1/std::tan(angle);
            W.coeffRef(id2,id3) += w23;
            W.coeffRef(id3,id2) += w23;
        }
    }
    
    
//    std::cout<<"Weight matrix: "<<std::endl;
//    std::cout<<W<<std::endl;
    
    //compute the system matrix
//    WW = W.sum(axis=0);
//    WWDiag = sparse.dia_matrix((WW, 0), shape=(n,n) );
//    WWDiagInv = sparse.dia_matrix((1/WW, 0), shape=(n,n) );
//    
//    C = W*WWDiagInv*W.transpose()-2*W+WWDiag
    Eigen::VectorXd WW;
    SumSparseMatrixCols(W, WW);
    Eigen::VectorXd WWInv = WW.cwiseInverse();
    
    //std::cout<<"Sums: "<<WW<<std::endl;
        
    //Wt = SparseDoubleMatrixType(W.transpose()); //Eigen does not support direct W*W'
    m_C = W * WWInv.asDiagonal() * SparseDoubleMatrixType(W.transpose()) - 2*W; // + WW.asDiagonal();
    
    //since adding diagonal matrix to the sparse is not supported by Eigen, add diagonal stuff manually
    for( Eigen::Index i=0; i<WW.innerSize(); i++)
        m_C.coeffRef(i,i) += WW(i);
    
//    std::cout<<"C: "<<m_C<<std::endl;
    
}




void UmbrellaWeightedOrder2Smoother::SumSparseMatrixCols( const SparseDoubleMatrixType& m, Eigen::VectorXd& s )
{
    s.resize(m.outerSize());
    
    Eigen::Index k;
    
    for (k = 0; k < m.outerSize(); ++k){
        s[k] = 0;
//        std::cout<<"k= "<<k<<std::endl;
        for (SparseDoubleMatrixType::InnerIterator it(m, k); it; ++it) {
            s[k] += it.value();
//            std::cout<<"V= "<<it.value()<<std::endl;
//            std::cout<<"s[k]= "<<s[k]<<std::endl;

        }
    }
    
//    std::cout<<"s: "<<s<<std::endl;
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

    
//    C = ComputeWeightMatrix(v.T, f.T);
    CalculateWeightMatrix();

//            
//    C[b,:] = 0;
//    for i in range(b.size):
//        C[b[i],b[i]] = 1;
    AddBoundaryToWeigtMatrix(); //reset rows of C to 1 for the boundaries

//    
//    n = v.shape[0]
//    Rx = np.zeros(n); 
//    Rx[b] = v[b,0];
//    Ry = np.zeros(n); 
//    Ry[b] = v[b,1];
//    Rz = np.zeros(n); 
//    Rz[b] = v[b,2];
//the boundary vertices should stay at the same place
    FormRightHandSide();
    
    
//  Solve the system  
//    result = np.zeros( (Rx.size, 3) )
//    result[:,0] = linalg_sp.spsolve(C, Rx); #x
//    result[:,1] = linalg_sp.spsolve(C, Ry); #y
//    result[:,2] = linalg_sp.spsolve(C, Rz); #y
//
//    return result

//    std::cout<<"Matrix C: "<<std::endl<<m_C<<std::endl;
//    std::cout<<"Matrix b: "<<std::endl<<m_b<<std::endl;
    
    Eigen::BiCGSTAB < Eigen::SparseMatrix<double> > cg;
    cg.compute(m_C);
    m_x = cg.solve(m_b);
    std::cout << "Solving sparse system using BiCGSTAB" << std::endl;
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error() << std::endl;     
    
//    std::cout<<"Matrix x: "<<std::endl<<m_x<<std::endl;

    TestBoundaryConservation();
    
    CreateOutput();
}


void UmbrellaWeightedOrder2Smoother::TestBoundaryConservation()
{
    //compare calculated m_x to the m_coverVertices[m_boundaryIds] to make sure they are the same
    //throw an error otherwise
    Eigen::RowVector3d v;
    double error = 0;
    for(Eigen::Index i=0; i<m_boundaryIds->size(); i++)
    {
        m_coverVertices->GetPoint(i, v.data());
        const Eigen::Index k=(*m_boundaryIds)[i];

        error += (m_x.row(k) - v).norm();
    }        

    const double eps = m_boundaryIds->size() * std::numeric_limits<double>::epsilon();
    
    if(error>eps )
    {
        std::cout<<"Linear system solver did converge to an accurate result"<<std::endl;
        std::cout<<"Matrix C: "<<std::endl<<m_C<<std::endl;
        std::cout<<"Matrix b: "<<std::endl<<m_b<<std::endl;
        std::cout<<"Matrix x: "<<std::endl<<m_x<<std::endl;    
    }
    
    assert( error <= eps );
    
}





void UmbrellaWeightedOrder2Smoother::FormRightHandSide()
{
    //intialize to 0
    m_b.resize(m_C.rows(), 3);
    m_b.fill(0);
    
    Eigen::RowVector3d v;
    for(Eigen::Index i=0; i<m_boundaryIds->size(); i++)
    {
        m_coverVertices->GetPoint(i, v.data());
        const Eigen::Index k=(*m_boundaryIds)[i];

        m_b.row(k) = v;
    }
}




void UmbrellaWeightedOrder2Smoother::AddBoundaryToWeigtMatrix()
{
    //reset rows of C to 1 for the boundaries
    
    
    for(Eigen::Index i=0; i<m_boundaryIds->size(); i++)
    {
        const Eigen::Index k=(*m_boundaryIds)[i];

        //reset the k-th row
        for (Eigen::Index j = 0; j < m_C.outerSize(); ++j)
            for (SparseDoubleMatrixType::InnerIterator it(m_C, j); it; ++it) {
                if(it.row()==k)
                    m_C.coeffRef(it.row(), it.col()) = 0;
            }
        
        
        m_C.coeffRef(k,k) = 1;
    }    
    
    
    m_C.prune(0,0); //remove zeros
}


void UmbrellaWeightedOrder2Smoother::CreateOutput()
{
    m_smoothedCoverPoints = vtkSmartPointer<vtkPoints>::New();
    
    //generate vtkPoints
    Eigen::RowVector3d v;
    for( Eigen::Index i=0; i<m_x.rows(); i++)
    {
        v = m_x.row(i);
        m_smoothedCoverPoints->InsertNextPoint( v.data() );
    }
}