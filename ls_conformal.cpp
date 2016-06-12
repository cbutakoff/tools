/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file 
    \brief Least squares conformal map. Based on "Least squares conformal maps for automatic texture atlas generation" by Bruno Levy
 */
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkCell.h>
#include <vtkPointData.h>
#include <vtkShortArray.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkFeatureEdges.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPointLocator.h>

#include "CommonTools.h"

#include <math.h>
#include <armadillo>



#include <stdio.h>
#include <stdlib.h>
#include <vtkCellArray.h>

#include <fstream>

#include <Eigen/SparseCore>
#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>


#define USE_REGULARIZATION 
//#define USE_APEX
#define USE_EIGENSOLVER



//#define RUN_TESTS
void ComputeW(vtkPolyData* mesh, arma::mat& W_real, arma::mat& W_imag, arma::colvec& dTs);
vtkPolyData* ComputeGradient(vtkPolyData* mesh, arma::mat& W_real, arma::mat& W_imag, arma::colvec& dTs, const char *ptarray);
void ComputeM(vtkPolyData* mesh, arma::mat& W_real, arma::mat& W_imag, arma::colvec& dTs, arma::sp_mat& M_real, arma::sp_mat& M_imag);
void ExtractBoundaryIds(vtkPolyData* mesh, arma::uvec& boundary_ids);
void SplitM(arma::uvec& fixed_ids, arma::sp_mat& M_real, arma::sp_mat& M_imag, arma::sp_mat& Mf_real, arma::sp_mat& Mf_imag, arma::sp_mat& Mp_real, arma::sp_mat& Mp_imag);
void GeneratePinnedPointsR2(long npoints, arma::colvec& Up1, arma::colvec& Up2);

void ComputeB(arma::sp_mat& Mp_real, arma::sp_mat& Mp_imag, arma::colvec& Up1, arma::colvec& Up2, arma::colvec& b);

void SortClosedPolyline(vtkPolyData* polyline, vtkPoints* sorted_points);

void SolveSparseSystemEigen(Eigen::SparseMatrix<double>& A, arma::colvec& b, arma::colvec& x);

void AssembleMatrix(arma::sp_mat& A, arma::sp_mat& B, Eigen::SparseMatrix<double>& C);
void AssembleMatrix(arma::sp_mat& A, arma::sp_mat& B, arma::sp_mat& C);
void SPMAT_ElemIndex(arma::sp_mat& A, unsigned int linear_ind, unsigned int& row, unsigned int& col);

void FlattenTriangle(arma::colvec3& p0, arma::colvec3& p1, arma::colvec3& p2, arma::colvec2& q0, arma::colvec2& q1, arma::colvec2& q2);
double TriangleGradientMatrix(arma::colvec2& q0, arma::colvec2& q1, arma::colvec2& q2, arma::mat& T,  bool normalize=true);
double TriangleGradientMatrix(vtkPolyData* mesh, vtkIdType pt1id, vtkIdType pt2id, vtkIdType pt3id, arma::mat& T,  bool normalize=true);


void CreateRegularizationTerm(vtkPolyData* mesh, arma::uvec& fixed_vert_ids, 
        arma::colvec& Up1, arma::colvec& Up2, 
        Eigen::SparseMatrix<double>& A, arma::colvec& B);

vtkIdType FindTriangle3rdPoint(vtkPolyData* mesh, vtkIdType cellId, vtkIdType p1, vtkIdType p2);



int main(int argc, char **argv) {
    
    const char *filename = argv[1];


    //---------------------------------------
    //
    // Prepare the mesh
    //

    vtkSmartPointer<vtkPolyData> mesh = vtkSmartPointer<vtkPolyData>::Take(CommonTools::LoadShapeFromFile(filename));

    std::cout << "Generating normals" << std::endl;

    vtkSmartPointer<vtkPolyDataNormals> normal_gen = vtkSmartPointer<vtkPolyDataNormals>::New();
    normal_gen->SetInputData(mesh);
    normal_gen->SplittingOff();
    normal_gen->ComputePointNormalsOff();
    normal_gen->ComputeCellNormalsOn();
    normal_gen->Update();

    //-----------------------------------------------------------------------
    //
    // calculate the weights W both real and complex part (for every triangle 3 values)
    // associated to every vertex
    arma::mat W_real;
    arma::mat W_imag;
    arma::colvec dTs;

    std::cout << "Computing W" << std::endl;
    ComputeW(normal_gen->GetOutput(), W_real, W_imag, dTs);


    //test gradient
#ifdef RUN_TESTS
    vtkSmartPointer<vtkPolyData> grad_mesh = vtkSmartPointer<vtkPolyData>::Take(
            ComputeGradient(normal_gen->GetOutput(), W_real, W_imag, dTs, "Ids"));
    CommonTools::SavePolydata(grad_mesh, "/home/costa/Copy/Ideas/LSConformal/test_gradient.vtk");
#endif

    arma::sp_mat M_real;
    arma::sp_mat M_imag;

    std::cout << "Computing M" << std::endl;
    ComputeM(mesh, W_real, W_imag, dTs, M_real, M_imag);


    //free up space
    W_real.reset();
    W_imag.reset();
    //:~

    //---------------------------
    // 
    //  Prepare the boundary
    //

    //extract boundary point ids
    arma::uvec boundary_ids;

    std::cout << "Extracting boundary" << std::endl;
    ExtractBoundaryIds(mesh, boundary_ids);

    //generate target points for the pinned points
    arma::colvec Up1;
    arma::colvec Up2;
    const long nboundary_points = boundary_ids.size(); //i expect in the future to have labels in "boundary"

    std::cout << "Generating boundary" << std::endl;
    GeneratePinnedPointsR2(nboundary_points, Up1, Up2);


#ifdef USE_APEX
    //insert apex
    boundary_ids.insert_rows(boundary_ids.size(),1,true);
    boundary_ids(boundary_ids.size()-1)=38;
    Up1.insert_rows(Up1.size(),1,true);
    Up2.insert_rows(Up2.size(),1,true);
#endif

    
    
    //reorder Up1 and Up2 according to the pointID 
    arma::uvec order = arma::sort_index(boundary_ids);
    Up1 = Up1(order);
    Up2 = Up2(order);


    //std::cout<<boundary_ids<<std::endl;
    
    //------------------------------------------
    //
    //  Calculate the remainder of matrices
    //

    //split matrix M
    arma::sp_mat Mp_real;
    arma::sp_mat Mp_imag;
    arma::sp_mat Mf_real;
    arma::sp_mat Mf_imag;

    std::cout << "Splitting M" << std::endl;
    SplitM(boundary_ids, M_real, M_imag, Mf_real, Mf_imag, Mp_real, Mp_imag);

    //free up space
    M_real.reset();
    M_imag.reset();
    //:~



    //create matrices A and b
    Eigen::SparseMatrix<double> A;
    arma::colvec b;

    std::cout << "Computing A" << std::endl;
    AssembleMatrix(Mf_real, Mf_imag, A);

//    std::ofstream ff;
//    arma::mat AAA(A);
//    ff.open("/home/costa/Copy/Ideas/LSConformal/A.txt");
//    AAA.print(ff);
//    ff.close();
    //    ff.open("/home/costa/Copy/Ideas/LSConformal/Aold.txt");
    //    Aold.print(ff);
    //    ff.close();

    Mf_real.reset();
    Mf_imag.reset();

    std::cout << "Computing B" << std::endl;
    ComputeB(Mp_real, Mp_imag, Up1, Up2, b);


    //    std::ofstream ff1;
    //    ff1.open("/home/costa/Copy/Ideas/LSConformal/b.txt");
    //    b.print(ff1);
    //    ff1.close();
    //    ff1.open("/home/costa/Copy/Ideas/LSConformal/b_old.txt");
    //    b_old.print(ff1);
    //    ff1.close();

    Mp_real.reset();
    Mp_imag.reset();

    
    //include regularization
#ifdef USE_REGULARIZATION    
    CreateRegularizationTerm(mesh, boundary_ids, Up1, Up2, A, b);
#endif
    
    
    //normalize the rows of the system
    std::cout << "Normalizing rows of the system" << std::endl;
    for(int i=0; i<A.rows(); i++)
    {
        double norm = 0;
        for (int j=0; j<A.cols(); j++)
            norm += A.coeff(i,j)*A.coeff(i,j);
        
        
        
        norm = sqrt( norm + b(i)*b(i) );
        assert(norm>0);

        //std::cout<<"b("<<i<<")="<<b(i)<<" n "<<norm<<std::endl;
        b(i) = b(i)/norm;
        for (int j=0; j<A.cols(); j++)
            if( A.coeff(i,j)!=0 ) A.coeffRef(i,j) = A.coeff(i,j)/norm;

    }
    
    //std::cout<<"b "<<b<<std::endl;

    //solve the system of equations

    //arma::mat A1(A);


    std::cout << "Solving system " << A.rows() << "x" << A.cols() << std::endl;
    arma::colvec Uf;
    
#ifdef USE_EIGENSOLVER    
    SolveSparseSystemEigen(A, b, Uf);
#else        
    //try iterative algorithm to make it better
    Eigen::VectorXd b1(b.size());
    Eigen::VectorXd x(A.cols());

    for(int i=0; i<b.size(); i++) b1(i)=b(i);
    for(int i=0; i<x.size(); i++) x(i)=0;


    Eigen::SparseMatrix<double> G = A.transpose()*A; 
    Eigen::VectorXd c = -A.transpose()*b1;

    const double norm_b = b1.norm();
    double threshold = 1e-5 * norm_b;
        
    Eigen::VectorXd g = -(G*x+c);
    Eigen::VectorXd r = g;
    
    //std::cout<<"Error: "<<g.norm()<<" / "<<threshold<<std::endl;
    while( g.norm()>threshold ){
        Eigen::VectorXd p = G*r;
        const double rho_sq = p.norm();
        const double rho = rho_sq*rho_sq;
        const double sigma = r.dot(p);
        const double tau = g.dot(r);
        const double t = tau/sigma;
        x += t*r;
        g -= t*p;
        const double gamma = (t*t*rho-tau)/tau;
        r += gamma*r+g;
        //std::cout<<"Error: "<<g.norm()<<std::endl;
    }

    //copy the reuslt back
    Uf.set_size(x.size());
    for(int i=0; i<Uf.size(); i++) Uf(i)=x(i);
#endif    
    
 
    
    
    
    
    
    //////////////////////////////////////////////////////////
    //
    //  Saving result
    //
    
    std::cout << "Uf : " << Uf.size() << std::endl;

    arma::colvec Uf1 = Uf(arma::span(0, Uf.size() / 2 - 1));
    arma::colvec Uf2 = Uf(arma::span(Uf.size() / 2, Uf.size() - 1));


    //----------------------------------------------
    //
    //  create the new polydata
    //
    //---------------------------------------------
    std::cout << "Saving the result" << std::endl;

    vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
    pd->DeepCopy(mesh);

    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    pts->SetNumberOfPoints(mesh->GetNumberOfPoints());

    arma::Col<unsigned char> mask(mesh->GetNumberOfPoints());
    mask.fill(0);
    mask(boundary_ids).fill(1);

    //create the points
    long free_index = 0;
    long pinned_index = 0;
    for (int ptid = 0; ptid < mesh->GetNumberOfPoints(); ptid++) {
        if (mask(ptid) > 0) {
            pts->SetPoint(ptid, Up1(pinned_index), Up2(pinned_index), 0);
            pinned_index++;
        } else {
            pts->SetPoint(ptid, Uf1(free_index), Uf2(free_index), 0);
            free_index++;
        }
    }

    pd->SetPoints(pts);

    CommonTools::SavePolydata(pd, argv[2]);


    return 0;
}



/*! \brief Create the rows of the system matrix corresponding to regularization 
 *  This has to be appended manually to the system matrix.
 *  The matrix has one row per interior edge (not boundary) and expresses the variation of the gradient across 
 *  that edge. The gradient is (approximately) projected onto the tangent plane to both triangles 
 *  The paper probably has a minor mistake. It treats the difference between the gradients as if the 
 *  both triangles were in the same plane. If we consider gradient projections to the tangent plane these 
 *  formulas would be correct. There might be a scaling factor involved (in case of tangent plane, 
 *  because the Y axis would form an angle with the plane and both triangles would form the same angle 
 *  with the plane) but then it probably would disappear once the equation "scale*variation = 0" is 
 *  formed and since X axis is not involved.
 * 
 *  see: Constrained Texture Mapping for Polygonal Meshes, Bruno Levy
 * 
 *  @param mesh - the mesh
 *  @param A - the matrix. Make sure it already has correct number of columns
 *  @param b - right hand side vector
 * 
 *  @return updated A and b 
 */
void CreateRegularizationTerm(vtkPolyData* mesh, arma::uvec& fixed_vert_ids, 
        arma::colvec& Up1, arma::colvec& Up2,
        Eigen::SparseMatrix<double>& A, arma::colvec& b)
{
    const int npoints = mesh->GetNumberOfPoints();
    arma::SpMat<char> edge_mask(npoints, npoints);
    
    mesh->BuildLinks();
    
    //to store triplets (row,col,value)
    std::vector<unsigned long> el_rows;  
    std::vector<unsigned long> el_cols;
    std::vector<double> el_values;
    
    unsigned long row = 0;
    
    
    //now we will solve Ax-b =0 the unknowns are only the nonfixed vertices (fixed ones are excluded)
    //they are ordered. We need to create the mapping between the real point ids and 
    //the indices within the array of unknowns which is smaller.
    //
    // E.g. vertices 0, 1, 2, 3. vertex 1 is fixed. Then unknowns will be 0,2,3
    // The correspondence will be 0 -> 0; 2 -> 1; 3 -> 2
    arma::Col<long int> vert_ids( mesh->GetNumberOfPoints() );
    vert_ids.fill(-1);
    
    arma::Col<unsigned char> mask( mesh->GetNumberOfPoints(), arma::fill::zeros );
    mask(fixed_vert_ids).fill(1);    

    long c_free=0;
    long c_pinned=0;
    for (int i = 0; i < mask.size(); i++) {
        if(mask(i)==0)
            vert_ids(i) = c_free++;
        else
            vert_ids(i) = c_pinned++;
    }

    
    
    //std::cout<<"free_vert_ids "<<vert_ids<<std::endl;
    //std::cout<<"pinned vertices "<<fixed_vert_ids<<std::endl;
    
    std::vector<double> regBu1_v;  //for coords Up1 
    std::vector<double> regBu2_v;  //for coords Up2  

    
    //create the triplets to fill the system
    for (int cellId = 0; cellId < mesh->GetNumberOfCells(); cellId++) {
        vtkCell* triangle = mesh->GetCell(cellId);
        
        for (int edgeId = 0; edgeId < triangle->GetNumberOfEdges(); edgeId++) {
            vtkCell* edge = triangle->GetEdge(edgeId);
            const vtkIdType p1 = edge->GetPointId(0);
            const vtkIdType p2 = edge->GetPointId(1);
            
            //the edge has not been processed yet
            if( edge_mask(p1, p2)==0 )
            {
                //get adjacent triangle
                vtkSmartPointer<vtkIdList> neighborCellIds = vtkSmartPointer<vtkIdList>::New();
                mesh->GetCellEdgeNeighbors( cellId, p1, p2, neighborCellIds);
                
                //make sure it's not an edge
                if( neighborCellIds->GetNumberOfIds()>0 )
                {
                    const vtkIdType adjCellId = neighborCellIds->GetId(0);
                    //now cellId has the cell and adjCellId is the adjacent cell 
                    //sharing the edge (p1, p2)
                    
                    //now make sure we don't mess up the point order. We need 2 triangles
                    //(p3, p1, p2) of cellId and (q3, p1 ,p2) of adjCellId
                    
                    //start with cellId
                    const vtkIdType p3 = FindTriangle3rdPoint( mesh, cellId, p1, p2 );
                    const vtkIdType q3 = FindTriangle3rdPoint( mesh, adjCellId, p1, p2);
                                        
                    arma::mat T; //grad of the cellId
                    arma::mat adjT; //grad of the adjCellId
                    (void) TriangleGradientMatrix(mesh, p3, p1, p2, T, true);
                    (void) TriangleGradientMatrix(mesh, q3, p1, p2, adjT, true);
                    
                    regBu1_v.push_back(0);
                    regBu2_v.push_back(0);
                    
                    if( mask(p1)==0 )
                    {
                        el_rows.push_back(row);
                        el_cols.push_back( vert_ids(p1) );
                        el_values.push_back(T(1,1)+adjT(1,1));                    
                    }
                    else
                    {
                        //std::cout<<"p1: "<<vert_ids(p1)<<" / "<<Up1.size()<<std::endl;
                        regBu1_v[row] -= (T(1,1)+adjT(1,1))*Up1( vert_ids(p1) );
                        regBu2_v[row] -= (T(1,1)+adjT(1,1))*Up2( vert_ids(p1) );
                    }

                    if( mask(p2)==0 )
                    {
                        el_rows.push_back(row);
                        el_cols.push_back( vert_ids(p2) );
                        el_values.push_back(T(1,2)+adjT(1,2));                    
                    }
                    else
                    {
                        //std::cout<<"p2: "<<vert_ids(p2)<<" / "<<Up1.size()<<std::endl;
                        regBu1_v[row] -= (T(1,2)+adjT(1,2))*Up1( vert_ids(p2) );
                        regBu2_v[row] -= (T(1,2)+adjT(1,2))*Up2( vert_ids(p2) );
                    }

                    if( mask(p3)==0 )
                    {
                        el_rows.push_back(row);
                        el_cols.push_back( vert_ids(p3) );
                        el_values.push_back(T(1,0));                    
                    }
                    else
                    {
                        //std::cout<<"p3: "<<vert_ids(p3)<<" / "<<Up1.size()<<std::endl;
                        regBu1_v[row] -= T(1,0)*Up1( vert_ids(p3) );
                        regBu2_v[row] -= T(1,0)*Up2( vert_ids(p3) );
                    }

                    if( mask(q3)==0 )
                    {
                        el_rows.push_back(row);
                        el_cols.push_back( vert_ids(q3) );
                        el_values.push_back(adjT(1,0));                    
                    }
                    else
                    {
                        //std::cout<<"q3: "<<q3<<"->"<<vert_ids(q3)<<" / "<<Up1.size()<<std::endl;
                        regBu1_v[row] -= adjT(1,0)*Up1( vert_ids(q3) );
                        regBu2_v[row] -= adjT(1,0)*Up2( vert_ids(q3) );
                    }

                    row++;

                    edge_mask(p1, p2)==1;
                    edge_mask(p2, p1)==1;
                }
            }
        }    
    }

    //edge_mask.reset();
    const long nrows_reg = el_rows[el_rows.size()-1]+1; 
    const long nrows = A.rows();
    const long ncols = 2*(mesh->GetNumberOfPoints()-fixed_vert_ids.size()); //A.cols();
    const long half_cols = ncols/2; //this is integer because the number of columns = 2*number of vertices
    
    std::cout<<"cols check: "<<ncols<<" vs "<<A.cols()<<std::endl;
    
    
    //nrows_reg*2 -- each coordinate has separate rows
//    A.conservativeResize( nrows + nrows_reg*2, ncols );
//
//    arma::colvec regBu1(nrows_reg, arma::fill::zeros);  //for coords Up1 
//    arma::colvec regBu2(nrows_reg, arma::fill::zeros);  //for coords Up2  
//
//    
//    for (int i = 0; i < el_rows.size(); i++) {
//        A.coeffRef(el_rows[i]+nrows, el_cols[i]) = el_values[i];
//        A.coeffRef(el_rows[i]+nrows+nrows_reg, el_cols[i]+half_cols) = el_values[i];
//    }
//
//    for(int i = 0; i<regBu1_v.size(); i++ ) regBu1(i) = regBu1_v[i];
//    for(int i = 0; i<regBu2_v.size(); i++ ) regBu2(i) = regBu2_v[i];
//    
//    b = arma::join_vert(b, arma::join_vert(regBu1,regBu2));

    
    A.conservativeResize( nrows + nrows_reg, ncols );

    arma::colvec regBu1(nrows_reg, arma::fill::zeros);  //for coords Up1 
    arma::colvec regBu2(nrows_reg, arma::fill::zeros);  //for coords Up2  
        
    for (int i = 0; i < el_rows.size(); i++) {
        A.coeffRef(el_rows[i]+nrows, el_cols[i]) = el_values[i];
        A.coeffRef(el_rows[i]+nrows, el_cols[i]+half_cols) = el_values[i];
    }
    

    for(int i = 0; i<regBu1_v.size(); i++ )
    {
        regBu1(i) = regBu1_v[i];
        regBu2(i) = regBu2_v[i];
    }
    

    b = arma::join_vert(b, regBu1+regBu2);
}














/*! \brief Calculate the gradient matrix for a 2D triangle (q0,q1,q2)
 * T = y2-y3  y3-y1  y1-y2 == TX
 *     x2-X3  x3-x1  x1-x2 == TY
 * See "Least Squares Conformal Maps for Automatic Texture Atlas Generation" by Levy
 * 
 *  @param normalize - (default true) if false, the matrix is not normalized by the area (if you plan to do it later)
 * 
 *  @return W 2x3 matrix, the return value is twice the triangle area. 
 */
double TriangleGradientMatrix(arma::colvec2& q0, arma::colvec2& q1, arma::colvec2& q2, arma::mat& T, bool normalize)
{
    T.set_size(2,3);
    T(0,0) = q1(1) - q2(1);
    T(0,1) = q2(1) - q0(1);
    T(0,2) = q0(1) - q1(1);

    T(1,0) = q1(0) - q2(0);
    T(1,1) = q2(0) - q0(0);
    T(1,2) = q0(0) - q1(0);

    
    //double the area of the triangle
    double area2 =  (q0(0)*q1(1) - q0(1)*q1(0)) + 
                    (q1(0)*q2(1) - q1(1)*q2(0)) + 
                    (q2(0)*q0(1) - q2(1)*q0(0));    
    
    if(normalize)
        T/=area2;
    
    return area2;
}




/*! \brief Flatten the 3D triangle (p0,p1,p2)
 * Axes: v0 = p1-p0, v0 x (p2-p1) x v0
 * Origin: p0
 * 
 *  @return q0, q1, q2
 */
void FlattenTriangle(arma::colvec3& p0, arma::colvec3& p1, arma::colvec3& p2, arma::colvec2& q0, arma::colvec2& q1, arma::colvec2& q2)
{
    //project the points to the local coordinate system
    //define the axes
    arma::colvec3 v0 = p2 - p1;
    double pt12_dist = arma::norm(v0, 2);
    v0 = arma::normalise(v0, 2);
    arma::colvec3 v1 = arma::cross( arma::cross(v0, p2-p0), v0 );
    v1 = arma::normalise(v1, 2);

    //in local coordinates v0, v1
    //pt1 will have coordinates (0,0)
    //pt2 will have coordinates ( ||pt2-pt1||, 0  ), because it lies on v0
    //pt3 will have coordinates ( v0.(pt3-pt1), v1.(pt3-pt1) )
    q0.fill(0);    
    q1.fill(0);
    
    q1(0) = pt12_dist;

    q2(0) = arma::dot(v0, p2 - p0);
    q2(1) = arma::dot(v1, p2 - p0);    
}



/*! \brief Create a matrix C from A and B
 * C = | A -B |
 *     | B  A |
 *
 *  @return x
 */
void AssembleMatrix(arma::sp_mat& A, arma::sp_mat& B, Eigen::SparseMatrix<double>& C) {
    const unsigned int n_colA = A.n_cols;
    const unsigned int n_rowA = A.n_rows;
    const unsigned int n_colB = B.n_cols;
    const unsigned int n_rowB = B.n_rows;

    const unsigned int nzA = A.col_ptrs[n_colA];
    const unsigned int nzB = B.col_ptrs[n_colB];

    C.setZero();
    C.resize(n_rowA + n_rowB, n_colA + n_colB);

    //fill in matrix A
    unsigned int row;
    unsigned int col;
    for (int i = 0; i < nzA; i++) {
        SPMAT_ElemIndex(A, i, row, col);
        C.coeffRef(row, col) = A.values[i];
        C.coeffRef(row + n_rowA, col + n_colA) = A.values[i];
    }

    //fill in matrix B
    for (int i = 0; i < nzB; i++) {
        SPMAT_ElemIndex(B, i, row, col);
        C.coeffRef(row, col + n_colA) = -B.values[i];
        C.coeffRef(row + n_rowA, col) = B.values[i];
    }
}




/*! \brief Create a matrix C from A and B
 * C = | A -B |
 *     | B  A |
 *
 *  @return x
 */
void AssembleMatrix(arma::sp_mat& A, arma::sp_mat& B, arma::sp_mat& C) {
    const unsigned int n_colA = A.n_cols;
    const unsigned int n_rowA = A.n_rows;
    const unsigned int n_colB = B.n_cols;
    const unsigned int n_rowB = B.n_rows;

    const unsigned int nzA = A.col_ptrs[n_colA];
    const unsigned int nzB = B.col_ptrs[n_colB];

    C.reset();
    C.set_size(n_rowA + n_rowB, n_colA + n_colB);

    //fill in matrix A
    unsigned int row;
    unsigned int col;
    for (int i = 0; i < nzA; i++) {
        SPMAT_ElemIndex(A, i, row, col);
        C(row, col) = A.values[i];
        C(row + n_rowA, col + n_colA) = A.values[i];
    }

    //fill in matrix B
    for (int i = 0; i < nzB; i++) {
        SPMAT_ElemIndex(B, i, row, col);
        C(row, col + n_colA) = -B.values[i];
        C(row + n_rowA, col) = B.values[i];
    }
}





void SPMAT_ElemIndex(arma::sp_mat& A, unsigned int linear_ind, unsigned int& row, unsigned int& col) {

    row = A.row_indices[linear_ind];
    col = 0;
    while (!(A.col_ptrs[col] <= linear_ind && A.col_ptrs[col + 1] > linear_ind)) ++col;
}






/*! \brief Solve Ax=b using UMFPACK
 *
 *  @return x
 */
void SolveSparseSystemEigen(Eigen::SparseMatrix<double>& A, arma::colvec& b, arma::colvec& x) {

//    Eigen::SparseMatrix<double> A1(A.n_rows, A.n_cols);
//
//    const unsigned int nzA = A.col_ptrs[A.n_cols]; //number of nonzero elements
//    unsigned int row;
//    unsigned int col;
//    for (int i = 0; i < nzA; i++) {
//        SPMAT_ElemIndex(A, i, row, col);
//        A1.coeffRef(row, col) = A.values[i];
//    }
//    
    A.makeCompressed();

    Eigen::VectorXd b1(b.size());
    for(int i=0; i<b.size(); i++) b1(i) = b(i);

//    std::cout << A1 << std::endl;
//    std::cout << b1 << std::endl;

    Eigen::SparseQR< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
    //Eigen::ConjugateGradient< Eigen::SparseMatrix<double> > solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        // decomposition failed
        std::cout << "Eigen solver, decomposition failed";
    }
    
    Eigen::VectorXd x1 = solver.solve(b1);
    if (solver.info() != Eigen::Success) {
        // solving failed
        std::cout << "Eigen solver, solving failed";
    }

//    std::cout << "Solution: " << x1 << std::endl;
    x.set_size(x1.size());
    for(int i=0; i<x.size(); i++) x(i)=x1(i);

}

/*! \brief Generate vector b as in (4)
 *
 *  @return b
 */
void ComputeB(arma::sp_mat& Mp_real, arma::sp_mat& Mp_imag, arma::colvec& Up1, arma::colvec& Up2, arma::colvec& b) {

    arma::sp_mat M;

    AssembleMatrix(Mp_real, Mp_imag, M);

    b = -M * arma::join_vert(Up1, Up2);
}

/*! \brief Generate Up1, Up2 in eq (4).  
 *
 *  @param npoints - number of pointsarma::colvec y
 *  @return Up1, Up2 - x,y coordinates of the pinned points (boundary in the flattening)
 */
void GeneratePinnedPointsR2(long npoints, arma::colvec& Up1, arma::colvec& Up2) {
    //put the points on the circumference of radius 1
    arma::colvec angle = arma::linspace<arma::colvec>(0, 2 * M_PI, npoints + 1).subvec(0, npoints - 1);
    const double r = 1.0;


    Up1 = r * arma::cos(angle);
    Up2 = r * arma::sin(angle);
}





/*! \brief Split matrix M into free Mf and pinned Mp based on point ids
 *
 *  @return Mf and Mp
 */
void SplitM(arma::uvec& fixed_ids, arma::sp_mat& M_real, arma::sp_mat& M_imag, 
        arma::sp_mat& Mf_real, arma::sp_mat& Mf_imag, 
        arma::sp_mat& Mp_real, arma::sp_mat& Mp_imag) {
    //count number of pinned points
    const long n_pinned = fixed_ids.size();


    Mp_real.set_size(M_real.n_rows, n_pinned);
    Mp_imag.set_size(M_real.n_rows, n_pinned);
    Mf_real.set_size(M_real.n_rows, M_real.n_cols - n_pinned);
    Mf_imag.set_size(M_real.n_rows, M_real.n_cols - n_pinned);


    arma::Col<unsigned char> mask(M_real.n_cols);
    mask.fill(0);
    mask(fixed_ids).fill(1);



    //split the matrix
    long free_index = 0;
    long pinned_index = 0;
    for (int ptid = 0; ptid < mask.size(); ptid++) {
        if (mask(ptid) > 0) {
            Mp_real.col(pinned_index) = M_real.col(ptid);
            Mp_imag.col(pinned_index) = M_imag.col(ptid);
            pinned_index++;
        } else {
            Mf_real.col(free_index) = M_real.col(ptid);
            Mf_imag.col(free_index) = M_imag.col(ptid);
            free_index++;
        }

    }

}

/*! \brief Sort the vertices of the closed polyline
 *
 *  @return indicator set, 1-s are on the boundary (later store boundary Id for more boundaries)
 */
void SortClosedPolyline(vtkPolyData* polyline, vtkPoints* sorted_points) {



    const long npoints = polyline->GetNumberOfPoints();

    arma::Mat<unsigned char> A(npoints, npoints, arma::fill::zeros); //adjacency matrix

    for (int lineid = 0; lineid < polyline->GetNumberOfCells(); lineid++) {
        const vtkIdType pt1id = polyline->GetCell(lineid)->GetPointIds()->GetId(0);
        const vtkIdType pt2id = polyline->GetCell(lineid)->GetPointIds()->GetId(1);

        A(pt1id, pt2id) = 1;
        A(pt2id, pt1id) = 1;
    }

    //extract vertices in a sorted way
    vtkIdType prevId = 0; //this is point that we added before the last added point

    //insert first 2 points
    sorted_points->InsertNextPoint(polyline->GetPoint(prevId));


    int currId = 0;
    while (A(prevId, currId) == 0)
        ++currId;

    sorted_points->InsertNextPoint(polyline->GetPoint(currId));



    while (currId != 0) {
        int i = 0;
        while (A(currId, i) == 0 || i == prevId)
            i++;

        prevId = currId;
        currId = i;


        if (i != 0)
            sorted_points->InsertNextPoint(polyline->GetPoint(i));

    }

    //CommonTools::SavePoints(sorted_points,"/home/costa/Copy/Ideas/LSConformal/pts.vtk");
}

/*! \brief Find the boundary and return point ids (only 1 contour for now)
 *
 *  @return indicator set, 1-s are on the boundary (later store boundary Id for more boundaries)
 */
void ExtractBoundaryIds(vtkPolyData* mesh, arma::uvec& boundary_ids) {
    vtkSmartPointer<vtkFeatureEdges> fe = vtkSmartPointer<vtkFeatureEdges>::New();
    fe->SetInputData(mesh);
    fe->BoundaryEdgesOn();
    fe->NonManifoldEdgesOff();
    fe->FeatureEdgesOff();
    fe->ManifoldEdgesOff();
    fe->Update();

    vtkSmartPointer<vtkPolyDataConnectivityFilter> connect =
            vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    connect->SetInputData(fe->GetOutput());
    connect->Update();

    const int ncontours = connect->GetNumberOfExtractedRegions();
    const int contour_id = 1; //work with only 1 contour for now



    connect->AddSpecifiedRegion(contour_id - 1);
    connect->SetExtractionModeToSpecifiedRegions();
    connect->Update();
    vtkPolyData *edges = connect->GetOutput();


    boundary_ids.set_size(edges->GetNumberOfPoints());


    vtkSmartPointer<vtkPointLocator> loc = vtkSmartPointer<vtkPointLocator>::New();
    loc->SetDataSet(mesh);
    loc->BuildLocator();

    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    SortClosedPolyline(edges, pts);



    for (int ptid = 0; ptid < pts->GetNumberOfPoints(); ptid++) {
        boundary_ids(ptid) = loc->FindClosestPoint(pts->GetPoint(ptid));
        ;
    }

}

/*! \brief Construct matrix M from W and dT
 *
 *  @return M_real and M_imag
 */
void ComputeM(vtkPolyData* mesh, arma::mat& W_real, arma::mat& W_imag, arma::colvec& dTs, arma::sp_mat& M_real, arma::sp_mat& M_imag) {
    const long int nfaces = mesh->GetNumberOfCells();

    M_real.set_size(nfaces, mesh->GetNumberOfPoints());
    M_imag.set_size(nfaces, mesh->GetNumberOfPoints());

    for (vtkIdType faceid = 0; faceid < nfaces; faceid++) {
        vtkCell *cell = mesh->GetCell(faceid);
        const vtkIdType pt1id = cell->GetPointId(0);
        const vtkIdType pt2id = cell->GetPointId(1);
        const vtkIdType pt3id = cell->GetPointId(2);

        const double dts_sqrt = sqrt(dTs(faceid));

        M_real(faceid, pt1id) = W_real(faceid, 0) / dts_sqrt;
        M_real(faceid, pt2id) = W_real(faceid, 1) / dts_sqrt;
        M_real(faceid, pt3id) = W_real(faceid, 2) / dts_sqrt;

        M_imag(faceid, pt1id) = W_imag(faceid, 0) / dts_sqrt;
        M_imag(faceid, pt2id) = W_imag(faceid, 1) / dts_sqrt;
        M_imag(faceid, pt3id) = W_imag(faceid, 2) / dts_sqrt;
    }
}

/*! \brief calculate the gradient using the matrix W input not modified
 *  du/dx = 1/dT(i)  ( -W_imag(i,:)*u(i) )  --- i-triangle index, u-scalars at vertices
 *  du/dy = 1/dT(i)  ( W_real(i,:)*u(i) )
 *
 *  @param mesh vtkPolyData with generated FACE normals
 *  @param W_real armadillo matrix (of size nfaces x 3)
 *  @param W_imag armadillo matrix (of size nfaces x 3)
 *  @param dTs areas of triangles
 *  @param ptarray point array that contains scalars (vtkShortArray)
 *  @return gradient saved inside the mesh. The saved polydata contains cell centroids with vectors. Use glyphs to visualize
 */
vtkPolyData* ComputeGradient(vtkPolyData* mesh, arma::mat& W_real, arma::mat& W_imag, arma::colvec& dTs, const char *ptarray) {
    vtkDoubleArray* face_normals = static_cast<vtkDoubleArray*> (
            mesh->GetCellData()->GetNormals());


    vtkShortArray *scalars = dynamic_cast<vtkShortArray*> (mesh->GetPointData()->GetArray(ptarray));

    const long int nfaces = mesh->GetNumberOfCells();

    vtkSmartPointer<vtkFloatArray> gradient = vtkSmartPointer<vtkFloatArray>::New();
    gradient->SetNumberOfComponents(3);
    gradient->SetNumberOfTuples(nfaces);
    gradient->SetName("Gradient");

    vtkSmartPointer<vtkPoints> centroids = vtkSmartPointer<vtkPoints>::New();
    centroids->SetNumberOfPoints(nfaces);


    for (vtkIdType faceid = 0; faceid < nfaces; faceid++) {
        vtkCell *cell = mesh->GetCell(faceid);
        const vtkIdType pt1id = cell->GetPointId(0);
        const vtkIdType pt2id = cell->GetPointId(1);
        const vtkIdType pt3id = cell->GetPointId(2);

        arma::colvec3 pt1_R3;
        arma::colvec3 pt2_R3;
        arma::colvec3 pt3_R3;

        mesh->GetPoint(pt1id, pt1_R3.colptr(0));
        mesh->GetPoint(pt2id, pt2_R3.colptr(0));
        mesh->GetPoint(pt3id, pt3_R3.colptr(0));

        //create the local coordinate system
        //define the axes
        arma::colvec3 v0 = pt2_R3 - pt1_R3;
        double pt12_dist = arma::norm(v0, 2);
        v0 = arma::normalise(v0, 2);
        arma::colvec3 n;
        face_normals->GetTuple(faceid, n.colptr(0));
        arma::colvec3 v1 = arma::cross(n, v0);


        //du/dx = 1/dT(i)  ( -W_imag(i,:)*u(i) )  --- i-triangle index, u-scalars at vertices
        //du/dy = 1/dT(i)  ( W_real(i,:)*u(i) )
        arma::colvec3 u;
        u(0) = scalars->GetValue(pt1id);
        u(1) = scalars->GetValue(pt2id);
        u(2) = scalars->GetValue(pt3id);


        const double dudx = arma::dot(-W_imag.row(faceid), u) / dTs(faceid);
        const double dudy = arma::dot(W_real.row(faceid), u) / dTs(faceid);

        arma::colvec3 grad = dudx * v0 + dudy*v1;

        arma::colvec3 cntr = (pt1_R3 + pt2_R3 + pt3_R3) / 3;
        centroids->SetPoint(faceid, cntr.memptr());
        gradient->SetTuple(faceid, grad.memptr());
    }

    //create the new polydata
    vtkPolyData* result = vtkPolyData::New();
    result->SetPoints(centroids);
    result->GetPointData()->AddArray(gradient);

    return result;
}

/*! \brief calculate the weights W both real and complex part 
 *         (for every triangle 3 values)
 *
 *  @param mesh vtkPolyData with generated FACE normals
 *  @param W_real armadillo matrix (of size nfaces x 3)
 *  @param W_imag armadillo matrix (of size nfaces x 3)
 *  @param dTs areas of triangles
 *  @return W_real, W_imag, dTs
 */
void ComputeW(vtkPolyData* mesh, arma::mat& W_real, arma::mat& W_imag, arma::colvec& dTs) {
    vtkDoubleArray* face_normals = static_cast<vtkDoubleArray*> (
            mesh->GetCellData()->GetNormals());

    const long int nfaces = mesh->GetNumberOfCells();

    W_real.set_size(nfaces, 3);
    W_imag.set_size(nfaces, 3);
    dTs.set_size(nfaces);

    for (vtkIdType faceid = 0; faceid < nfaces; faceid++) {
        vtkCell *cell = mesh->GetCell(faceid);
        const vtkIdType pt1id = cell->GetPointId(0);
        const vtkIdType pt2id = cell->GetPointId(1);
        const vtkIdType pt3id = cell->GetPointId(2);

        arma::mat T;
        dTs(faceid) = TriangleGradientMatrix(mesh, pt1id, pt2id, pt3id, T, false);

        W_real(faceid, arma::span::all) = -T.row(1);
        W_imag(faceid, arma::span::all) = -T.row(0);
    }


}













/*! \brief Return the Id of the 3rd point of the triangle
 *
 *  @param mesh - triangular mesh
 *  @param cellId - id of the triangular cell in which to search
 *  @param p1 - id of one point
 *  @param p2 - id of the other point
 *  @return id of the third point
 */
vtkIdType FindTriangle3rdPoint(vtkPolyData* mesh, vtkIdType cellId, vtkIdType p1, vtkIdType p2)
{
    vtkIdType npts;
    vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

    mesh->GetCellPoints (cellId, ptIds);
    
    vtkIdType p3 = 0;
    for (int i = 0; i < ptIds->GetNumberOfIds(); i++) {
        const vtkIdType id = ptIds->GetId(i);
        
        if( id!=p1 && id!=p2 ){
            p3 = id;
            break;
        }        
    }
        
    return p3;
}




/*! \brief Calculate the gradient matrix for a 3D triangle (q0,q1,q2)
 * Projects the triangle to 2D and calls the corresponding function
 * 
 *  @param normalize - (default true) if false, the matrix is not normalized by the area (if you plan to do it later)
 * 
 *  @return W 2x3 matrix, the return value is twice the triangle area. 
 */
double TriangleGradientMatrix(vtkPolyData* mesh, vtkIdType pt1id, vtkIdType pt2id, vtkIdType pt3id, arma::mat& T,  bool normalize)
{
    arma::colvec3 pt1_R3;
    arma::colvec3 pt2_R3;
    arma::colvec3 pt3_R3;

    mesh->GetPoint(pt1id, pt1_R3.colptr(0));
    mesh->GetPoint(pt2id, pt2_R3.colptr(0));
    mesh->GetPoint(pt3id, pt3_R3.colptr(0));

    //project the points to the local coordinate system
    //define the axes
    arma::colvec2 pt1_R2, pt2_R2, pt3_R2;
    FlattenTriangle(pt1_R3, pt2_R3, pt3_R3, pt1_R2, pt2_R2, pt3_R2);

    return TriangleGradientMatrix(pt1_R2, pt2_R2, pt3_R2, T, normalize);   
}





