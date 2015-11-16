/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file
    \brief Generalized Ricci Flow. Using notation from the paper in TPAMI
 * there were complilation errors marked with //!!!
 */
#include "CommonTools.h"

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkType.h>


#include <cfloat>
#include <math.h>
#include <armadillo>
#include <vtkIdList.h>

//------------------------------------------------------------------



void CalculateEdgeLengths(vtkPolyData* mesh, arma::sp_mat& edge_lengths);
void CalculateRadii(vtkPolyData* mesh, const arma::sp_mat& edge_lengths, arma::vec& radii);
void CalculateInversiveDistanceE2(const arma::sp_mat& l, const arma::vec& gamma, arma::sp_mat& I);
void CalculateGeneralizedEdgeLengthsE2(const arma::vec& gamma, const arma::sp_mat& I, arma::sp_mat& gen_edge_lengths);
void CalculateFaceAnglesE2(vtkPolyData* mesh, const arma::vec& edge_lengths, arma::vec& angles);
void CalculateAngleDerivativesE2(vtkPolyData* mesh, const arma::sp_mat& edge_lengths, const arma::vec& radii);


//reorders indices in1, ind2, ind3 (one of them is i) so that ind1 is i
void SortIndices(int i, int& ind1, int& ind2, int& ind3);

int main(int argc, char *argv[]) {
    std::cout << "GeneralizedRicciFlow for surfaces v0.1" << std::endl;

    if (argc < 1) {
        std::cout << "Usage: GeneralizedRicciFlow mesh.vtk" << std::endl;
        return -1;
    }

    char* inshape = argv[1];

    vtkSmartPointer<vtkPolyData> mesh = vtkSmartPointer<vtkPolyData>::Take(
            CommonTools::LoadShapeFromFile(inshape));


    //////////////////////////////////////////////////////////
    //
    //
    //      Calculate initial circle packing metric
    //
    //
    /////////////////////////////////////////////////////////////
    arma::sp_mat l;
    CalculateEdgeLengths(mesh, l);

    arma::vec gamma;
    CalculateRadii(mesh, l, gamma);

    arma::sp_mat I; //inversive distance
    CalculateInversiveDistanceE2(l, gamma, I);


    //////////////////////////////////////////////////////////
    //
    //
    //      Discrete Ricci Flow
    //
    //
    /////////////////////////////////////////////////////////////
    bool converged = false;
    while (!converged) {
        arma::sp_mat gen_edge_lengths;
        CalculateGeneralizedEdgeLengthsE2(gamma, I, gen_edge_lengths);

        //compute corner angles using cosine laws
        arma::vec angles;
//!!!        CalculateFaceAnglesE2(mesh, gen_edge_lengths, angles);
    }


    return 0;
}



/////////////////////////////////////////////////////////
//
//   Eq. (4)
//
//
//
//
//

void CalculateFaceAnglesE2(vtkPolyData* mesh, const arma::vec& edge_lengths, arma::vec& angles) {
    const int ncells = mesh->GetNumberOfCells();
    angles.set_size(ncells);

    vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();

    for (vtkIdType i = 0; i < ncells; i++) {
        pointIdList->Reset();
        mesh->GetCellPoints(i, pointIdList);

        int ind_i = pointIdList->GetId(0);
        int ind_j = pointIdList->GetId(1);
        int ind_k = pointIdList->GetId(2);

        SortIndices(i, ind_i, ind_j, ind_k);

        const int l_ij = edge_lengths(ind_i, ind_j);
        const int l_ki = edge_lengths(ind_k, ind_i);
        const int l_kj = edge_lengths(ind_k, ind_j);

        angles(i) = acos((l_ij * l_ij + l_ki * l_ki - l_kj * l_kj) / (2 * l_ki * l_ij));
    }
}




/////////////////////////////////////////////////////////
//
//   Eq. (8)
//
//
//
//
//

void CalculateGeneralizedEdgeLengthsE2(const arma::vec& gamma, const arma::sp_mat& I, arma::sp_mat& gen_edge_lengths) {
    const int npts = gamma.n_elem;
    gen_edge_lengths.set_size(npts, npts);

    //the measure is symmetric, process only upper triangle
    for (int i = 0; i < npts; i++) {
        for (int j = i + 1; j < npts; j++) {
            if (I(i, j) > 0) { //only for existing edges
                const double l_ij = sqrt(gamma(i) * gamma(i) + gamma(j) * gamma(j) + 2 * gamma(i) * gamma(j) * I(i, j));
                gen_edge_lengths(i, j) = l_ij;
                gen_edge_lengths(j, i) = l_ij;
            }
        }
    }

}




/////////////////////////////////////////////////////////
//
//   For every pair of vertices calculate the inversive distance
//
//
//
//
//

void CalculateInversiveDistanceE2(const arma::sp_mat& l, const arma::vec& gamma, arma::sp_mat& I) {
    const int npts = gamma.n_elem;
    I.set_size(npts, npts);

    //the measure is symmetric, process only upper triangle
    for (int i = 0; i < npts; i++) {
        for (int j = i + 1; j < npts; j++) {
            if (l(i, j) > 0) { //only for existing edges
                const double I_ij = (l(i, j) * l(i, j) - gamma(i) * gamma(i) - gamma(j) * gamma(j)) / (2 * gamma(i) * gamma(j));
                I(i, j) = I_ij;
                I(j, i) = I_ij;
            }
        }
    }

}




/////////////////////////////////////////////////////////
//
//   For every vertex calculate all the radii and take minimum
//
//
//
//
//

void CalculateRadii(vtkPolyData* mesh, const arma::sp_mat& edge_lengths, arma::vec& radii) {
    const int npts = mesh->GetNumberOfPoints();
    const int ncells = mesh->GetNumberOfCells();
    radii.set_size(npts);

    vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();

    //for every face calculate gamma_ijk
    arma::vec gamma_percell(ncells);
    for (vtkIdType i = 0; i < ncells; i++) {
        pointIdList->Reset();
        mesh->GetCellPoints(cellIdList->GetId(i), pointIdList);

//!!!        assert(pointIdList->GetNumberOfIds() == 3);

        int ind_i = pointIdList->GetId(0);
        int ind_j = pointIdList->GetId(1);
        int ind_k = pointIdList->GetId(2);

        SortIndices(i, ind_i, ind_j, ind_k);

        const double l_ij = edge_lengths(ind_i, ind_j);
        const double l_ki = edge_lengths(ind_k, ind_i);
        const double l_jk = edge_lengths(ind_j, ind_k);

        gamma_percell(i) = (l_ij + l_ki - l_jk) / 2;
    }

    //for every vertex get minimum gamma    
    for (vtkIdType i = 0; i < npts; i++) {
        cellIdList->Reset();
        mesh->GetPointCells(i, cellIdList);

        radii(i) = DBL_MAX;
        for (vtkIdType j = 0; j < cellIdList->GetNumberOfIds(); j++) {
            if (gamma_percell(j) < radii(i)) radii(i) = gamma_percell(j); //get minimum of gammas from surrounding cells
        }
    }

}



/////////////////////////////////////////////////////////
//
//
//
//
//
//
//

void CalculateEdgeLengths(vtkPolyData* mesh, arma::sp_mat& edge_lengths) {
    vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();

    const int npts = mesh->GetNumberOfPoints();
    edge_lengths.set_size(npts, npts);

    for (vtkIdType i = 0; i < npts; i++) {
        cellIdList->Reset();
        mesh->GetPointCells(i, cellIdList);
        double *pt1 = mesh->GetPoint(i);


        for (vtkIdType j = 0; j < cellIdList->GetNumberOfIds(); j++) {
            pointIdList->Reset();
            mesh->GetCellPoints(cellIdList->GetId(i), pointIdList);
            for (vtkIdType k = 0; k < pointIdList->GetNumberOfIds(); k++) {
                vtkIdType pt2_id = pointIdList->GetId(k);
                if (pt2_id != i) { //not the same point
                    if (edge_lengths(i, pt2_id) != 0) { //is not calculated yet
                        double *pt2 = mesh->GetPoint(pt2_id);

                        const double dx = pt1[0] - pt2[0];
                        const double dy = pt1[1] - pt2[1];
                        const double dz = pt1[2] - pt2[2];
                        const double d = sqrt(dx * dx + dy * dy + dz * dz);

                        edge_lengths(i, pt2_id) = d;
                        edge_lengths(pt2_id, i) = d;
                    }
                }
            }
        }
    }
}



/////////////////////////////////////////////////////////
//
//
//
//  Reorders indices in1, ind2, ind3 (one of them is i) so that ind1 is i. 
//  Keeps the orientation of triangle
//
//
//
//

void SortIndices(int i, int& ind1, int& ind2, int& ind3) {
    if (ind1 == i)
        return;
    else if (ind2 == i) {
        ind3 = ind1;
        ind2 = ind3;
        ind1 = i;
    } else if (ind3 == i) {
        ind3 = ind2;
        ind2 = ind1;
        ind1 = i;
    }

}





/////////////////////////////////////////////////////////
//
//
//
//  Calculate dTheta/du, eq 10
//
//
//

void CalculateAngleDerivativesE2(vtkPolyData* mesh, const arma::sp_mat& edge_lengths, const arma::vec& radii) {
    //calculate the power center for each triangle as follows
    //for 2 circles centered at c1 and c2 solve for p:
    //|p-c1|^2 - r1^2 = |p-c2|^2 - r2^2    -- this is power line
    //add one more condition for the third circle to get one point
    //  Funny part: the paper says that the circle of radius sqrt(|c1-p|^2 - r1^2)
    //  is perpendicular to all three circles. But this works only when the 3 circles do not intersect
    const int npts = mesh->GetNumberOfPoints();
    const int ncells = mesh->GetNumberOfCells();

    vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();

    //for every face calculate gamma_ijk
    arma::vec gamma_percell(ncells);
    for (vtkIdType i = 0; i < ncells; i++) {
        pointIdList->Reset();
        mesh->GetCellPoints(cellIdList->GetId(i), pointIdList);


        const int ind_i = pointIdList->GetId(0);
        const int ind_j = pointIdList->GetId(1);
        const int ind_k = pointIdList->GetId(2);

        const double ri = radii[ind_i];
        const double rj = radii[ind_j];
        const double rk = radii[ind_k];

        arma::vec ci(3), cj(3), ck(3);
        mesh->GetPoint(ind_i, ci.memptr());
        mesh->GetPoint(ind_j, cj.memptr());
        mesh->GetPoint(ind_k, ck.memptr());

        //create local 2D coordinates, and map vertices
        arma::vec vji = cj - ci;
        arma::vec v   = ck - ci;

        arma::vec n = arma::cross(vji, v); //normal
        arma::vec vki = arma::cross(n, vji); //creates local y' axis

        const double lenX = arma::norm(vji,2);
        const double lenY = arma::norm(vki,2);
        
        vji /= lenX; //normalize axes
        vki /= lenY;
        
        //2D coordinates
        arma::vec vi(3,0); //(0,0)
        arma::vec vj(3);
        arma::vec vk(3);

        vj(0) = lenX;
        vj(1) = 0.0;
        vk(0) = arma::dot(v, vji);
        vk(1) = arma::dot(v, vki);
        
        //form matrices for system of linear equations and solve
        //A = [2*(c2-c1)'; 2*(c2-c3)'];
        //B = [ c2'*c2 - c1'*c1 + r1^2 - r2^2; c2'*c2 - c3'*c3 + r3^2 - r2^2];
        //P = A\B;
        arma::mat A(2,2);
        arma::vec B(2);
        
        A.insert_rows( 0, 2.0*(cj-ci) );
        A.insert_rows( 1, 2.0*(cj-ck) );
        B(0) = arma::dot(cj, cj) - arma::dot(ci, ci) + ri*ri - rj*rj;
        B(1) = arma::dot(cj, cj) - arma::dot(ck, ck) + rk*rk - rj*rj;
        arma::vec P = arma::solve( A, B ); 
        
        //calculate distance from P to all the edges: hi, hj, hk
        //http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
        const double hi = arma::norm( arma::cross(cj-ck, ck-P), 2 ) / arma::norm(cj-ck, 2);
        const double hj = arma::norm( arma::cross(ci-ck, ck-P), 2 ) / arma::norm(ci-ck, 2);
        const double hk = arma::norm( arma::cross(cj-ci, ci-P), 2 ) / arma::norm(cj-ci, 2);
        
        const double li = edge_lengths(ind_j, ind_k);
        const double lj = edge_lengths(ind_i, ind_k);
        const double lk = edge_lengths(ind_i, ind_j);
        
        const double dTi_duj = hk/lk;
        const double dTj_duk = hi/li;
        const double dTk_dui = hj/lj;

        const double dTi_dui = - dTi_duj - dTk_dui;
        const double dTj_duj = - dTj_duk - dTi_duj;
        const double dTk_duk = - dTj_duk - dTk_dui;


    }
    
    
    
}

