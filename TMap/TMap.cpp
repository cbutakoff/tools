#include "TMap.h"

#include <vtkIdList.h>
#include <vtkTriangle.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>

#include <vtkType.h>
#include <complex>
#include <set>
#include <tuple>

#include<Eigen/SparseLU>

//#define DEBUG_MESSAGES


void TMapUtils::ValidateInputs() const
{
    if(m_InputMesh.GetPointer()==NULL)
        throw std::logic_error("Input mesh is not set");
    
}



/* calculates the BC with f given via 2 meshes. f maps the mesh from SetInputMesh to
 * tgtMesh passed as a prameter
 */
void TMapUtils::CalculateBC(vtkPolyData* tgtMesh)
{
    ValidateInputs();

    if(tgtMesh==NULL)
        throw std::logic_error("Target mesh is not set");

    
    const vtkIdType nTriangles = m_InputMesh->GetNumberOfCells();
    
    m_BC.resize(nTriangles);
    
    DenseMatrix2x2Type A; //left side matrix in (5.2)
    DenseMatrix2x2Type B; //right side matrix in (5.2)
    DenseMatrix2x2Type nabla_f;
    double v0[3], v1[3], v2[3];
    double fv0[3], fv1[3], fv2[3];

    
    for(vtkIdType i=0; i<nTriangles; i++)
    {
        if( m_InputMesh->GetCell(i)->GetNumberOfPoints()!= 3)
        {
            char error[100];
            sprintf(error,"Cell %lld is not a triangle",i);
            throw std::logic_error(error);
        }
            
        m_InputMesh->GetPoint( m_InputMesh->GetCell(i)->GetPointIds()->GetId(0), v0 );
        m_InputMesh->GetPoint( m_InputMesh->GetCell(i)->GetPointIds()->GetId(1), v1 );
        m_InputMesh->GetPoint( m_InputMesh->GetCell(i)->GetPointIds()->GetId(2), v2 );

        tgtMesh->GetPoint( tgtMesh->GetCell(i)->GetPointIds()->GetId(0), fv0 );
        tgtMesh->GetPoint( tgtMesh->GetCell(i)->GetPointIds()->GetId(1), fv1 );
        tgtMesh->GetPoint( tgtMesh->GetCell(i)->GetPointIds()->GetId(2), fv2 );

        /* remember we use only X, Y cordinates*/
        A(0,0) = v1[0] - v0[0];
        A(0,1) = v1[1] - v0[1];
        A(1,0) = v2[0] - v0[0];
        A(1,1) = v2[1] - v0[1];
        
        B(0,0) = fv1[0] - fv0[0]; //note that there should be 
        B(0,1) = fv1[1] - fv0[1]; //no normalization in the equation
        B(1,0) = fv2[0] - fv0[0]; //you can easily check that the 
        B(1,1) = fv2[1] - fv0[1]; //eq does not satisfy identity transform
        
        nabla_f = A.colPivHouseholderQr().solve(B);
        
        /* nabla_f = [[at, ct], [bt,dt]] */
        const double at = nabla_f(0,0);
        const double bt = nabla_f(1,0);
        const double ct = nabla_f(0,1);
        const double dt = nabla_f(1,1);
        
        m_BC(i) = ComplexType(at-dt, ct+bt) / ComplexType(at+dt, ct-bt);
    }
}




/* Adjacency matrix.  */
void TMapUtils::Precalculate()
{
    ValidateInputs();

    if( HasInputMeshChanged() )
    {
        ResetModificationFlag();

        this->m_CellsSharedByVertex.clear();
        this->m_CellsAdjacent2Cell.clear();
        this->m_AdjacencyMatrix.setZero();
        
        const vtkIdType n_points = m_InputMesh->GetNumberOfPoints();
        const vtkIdType n_cells = m_InputMesh->GetNumberOfCells();

        this->m_AdjacencyMatrix.resize(n_points, n_points);
        this->m_Areas.resize(n_cells,1);
        //this->m_Cotangents.resize(n_cells,3);

        this->m_CellsSharedByVertex.resize(n_points);
        this->m_CellsAdjacent2Cell.resize(n_cells);
        //for( auto vert_vector : this->m_CellsSharedByVertex ) vert_vector.clear(); /* clean up the array*/

        for( vtkIdType i=0; i<n_cells; i++)
        {
            /* Fill in the adjacency matrix */
            for( vtkIdType j=0; j< this->m_InputMesh->GetCell(i)->GetNumberOfEdges(); j++)
            {
                StopOnError(this->m_InputMesh->GetCell(i)->GetEdge(j)->GetNumberOfPoints()!=2, "Edge has wrong number of points");
                const vtkIdType pt0id = this->m_InputMesh->GetCell(i)->GetEdge(j)->GetPointId(0);
                const vtkIdType pt1id = this->m_InputMesh->GetCell(i)->GetEdge(j)->GetPointId(1);

                if( pt0id<pt1id )
                    this->m_AdjacencyMatrix.coeffRef( pt0id, pt1id ) += 1;
                else
                    this->m_AdjacencyMatrix.coeffRef( pt1id, pt0id ) += 1;
            }


            StopOnError(this->m_InputMesh->GetCell(i)->GetNumberOfPoints()!=3, "Triangle has wrong number of points");

            double vi[3], vj[3], vk[3];
            this->m_InputMesh->GetPoint(this->m_InputMesh->GetCell(i)->GetPointId(0), vi);
            this->m_InputMesh->GetPoint(this->m_InputMesh->GetCell(i)->GetPointId(1), vj);
            this->m_InputMesh->GetPoint(this->m_InputMesh->GetCell(i)->GetPointId(2), vk);

            m_Areas[i] = vtkTriangle::TriangleArea(vi, vj, vk);


            /* Fill in the matrix that stores which vertex is contained in which triangles 
             * push this triangle into the vectors of each of the cell vertices
             */
            this->m_CellsSharedByVertex[this->m_InputMesh->GetCell(i)->GetPointId(0)].push_back(i);
            this->m_CellsSharedByVertex[this->m_InputMesh->GetCell(i)->GetPointId(1)].push_back(i);
            this->m_CellsSharedByVertex[this->m_InputMesh->GetCell(i)->GetPointId(2)].push_back(i);
            
            /* Calculate the cotangents */
            //this->m_Cotangents(i,0) = 1.0/tan(this->GetAngle(vi,vj,vk));
            //this->m_Cotangents(i,1) = 1.0/tan(this->GetAngle(vj,vk,vi));
            //this->m_Cotangents(i,2) = 1.0/tan(this->GetAngle(vk,vi,vj));
        }

        //calculate the neighbors for every cell
        for( vtkIdType i=0; i<n_cells; i++)
        {
            vtkIdList* ids = this->m_InputMesh->GetCell(i)->GetPointIds();
            std::unordered_set<vtkIdType> neighbors;

            //for every vertex add its cells
            for( vtkIdType j=0; j<ids->GetNumberOfIds(); j++ )
                for( auto cell_id : this->m_CellsSharedByVertex[ ids->GetId(j) ] )
                    neighbors.insert(cell_id);

            for( auto cell_id : neighbors )
                m_CellsAdjacent2Cell[i].push_back(cell_id);        
        }

    }
    

}


/*Transforms global vertex IDs to local vertex ids inside the triangle = 0,1,2*/
/* lazy implementation */
std::tuple<vtkIdType, vtkIdType, vtkIdType> 
TMapUtils::VertexIdGlobal2Local(vtkIdType i, vtkIdType j, vtkIdType k, vtkIdType cellId)
{
    auto ptids = this->m_InputMesh->GetCell(cellId)->GetPointIds();
    std::vector<vtkIdType> ids( ptids->GetNumberOfIds() );
    for( vtkIdType i=0; i<ptids->GetNumberOfIds(); i++) ids[i] = ptids->GetId(i);

    //i
    auto it = std::find(ids.begin(),ids.end(),i);
    
    if(it==ids.end()) 
    {
        char message[100];
        sprintf(message, "VertexIdGlobal2Local: Cannot find vertex %ld in triangle %ld", (long int)i, (long int)cellId);
        throw std::logic_error( message );
    }

    const vtkIdType iloc = std::distance(ids.begin(), it);
    
    //j
    it = std::find(ids.begin(),ids.end(),j);
    
    if(it==ids.end()) 
    {
        char message[100];
        sprintf(message, "VertexIdGlobal2Local: Cannot find vertex %ld in triangle %ld", (long int)j, (long int)cellId);
        throw std::logic_error( message );
    }

    const vtkIdType jloc = std::distance(ids.begin(), it);
    
    //k
    it = std::find(ids.begin(),ids.end(),k);
    
    if(it==ids.end()) 
    {
        char message[100];
        sprintf(message, "VertexIdGlobal2Local: Cannot find vertex %ld in triangle %ld", (long int)k, (long int)cellId);
        throw std::logic_error( message );
    }

    const vtkIdType kloc = std::distance(ids.begin(), it);
    
    return std::tuple<vtkIdType, vtkIdType, vtkIdType> (iloc, jloc, kloc);
}





/* Return the angle formed by vectors v0v1 and v0v2*/
//double TMapUtils::GetAngle(const double* v0, const double* v1, const double* v2) const
//{
//    /* vectors */
//    Eigen::Vector3d v0v1({v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]});
//    Eigen::Vector3d v0v2({v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]});
//    v0v1.normalize();
//    v0v2.normalize();
//    const double dp = v0v1.dot(v0v2);
//    return acos(dp);
//}



void TMapUtils::MeshEdgePointIds(std::unordered_set<vtkIdType>& edgePtIds) 
{
    Precalculate();
    
    edgePtIds.clear();
    
    for (vtkIdType k=0; k<this->m_AdjacencyMatrix.outerSize(); ++k)
        for (AdjacencyMatrixType::InnerIterator it(this->m_AdjacencyMatrix,k); it; ++it)
        {
            if(it.value()==1) //one adjacent triangle = edge
            {
                edgePtIds.insert(it.row());   // one point
                edgePtIds.insert(it.col());   // another point
            }
        }    
}








 void TMapUtils::MeshFromBC()
 {
    Precalculate();
     
    if( m_LandmarkConstraints.size()<2 )
    {
        throw std::logic_error("Too few landmark constraints");
    }
    
#ifdef DEBUG_MESSAGES
    std::cout<<"Landmarks: ";
    for(auto& lmk : m_LandmarkConstraints)
        std::cout<<lmk.Id<<" "<<lmk.x<<" "<<lmk.y<<std::endl;
#endif
    
    
    if( m_OutputMesh.GetPointer()==NULL)
        m_OutputMesh = vtkSmartPointer<vtkPolyData>::New();
    
    m_OutputMesh->DeepCopy(m_InputMesh);
    
    Eigen::VectorXd x, y;
    TMapUtils::CalculateVerticesFromBC(x,y);
    
#ifdef DEBUG_MESSAGES
    std::cout<<"X: "<<x<<std::endl;
    std::cout<<"y: "<<y<<std::endl;
#endif

    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    pts->SetNumberOfPoints(m_InputMesh->GetNumberOfPoints());
    for( vtkIdType i=0; i<m_OutputMesh->GetNumberOfPoints(); i++)
    {
       pts->SetPoint(i, x(i), y(i), 0);
    }
    
    m_OutputMesh->SetPoints(pts);        
 }
 
 
 
 
void TMapUtils::CalculateAlphasFromBC(TMapUtils::DenseMatrixType& alpha)
{
    alpha.resize(this->m_BC.size(), 3);
    
    for(vtkIdType i=0; i<this->m_BC.size(); i++)
    {
        const double rho = m_BC(i).real();
        const double tau = m_BC(i).imag();
        const double denom = 1 - rho*rho - tau*tau;
        
        alpha(i,0) = ( (rho-1)*(rho-1) + tau*tau )/denom; //alpha1
        alpha(i,1) = -2*tau / denom; //alpha2
        alpha(i,2) = (1 + 2*rho + rho*rho + tau*tau )/denom; //alpha3
    }
}



void TMapUtils::CalculateVerticesFromBC(Eigen::VectorXd& x, Eigen::VectorXd& y)
{
    /* For each vertex create a system of equations, eq (5.12) and (5.13) */
    const vtkIdType npoints = m_InputMesh->GetNumberOfPoints();

    SysEqMatrixType L(npoints, npoints); /*Left side*/
    DenseMatrixType R(npoints, 2);  /* Right side - x,y */

    
    L.reserve(Eigen::VectorXi::Constant(npoints,3)); /*reserve 3 entries per row*/
    
    Eigen::Vector3d A,B;

    DenseMatrixType alphas;
    CalculateAlphasFromBC(alphas);
    
#ifdef DEBUG_MESSAGES
    std::cout<<"Alphas: "<<alphas<<std::endl;
#endif
    
    /* create a sorted array of ids of the constrainted landmarks for fast search */
    std::set<vtkIdType> landmark_ids;
    for(const auto& lmk : m_LandmarkConstraints )
        landmark_ids.insert(lmk.Id);
    /* set created */
    
    
    /* For each vertex */
    for(vtkIdType vertId = 0; vertId<npoints; vertId++)
    {
        //if(vertId%1000==0)
        //    std::cout<<"Creating system, vertex:"<<vertId<<"\r";
        
        /* Check if the vertex is part of the landmark constraints */
        auto search4lmk_id = landmark_ids.find(vertId);
        
        if( search4lmk_id == landmark_ids.end() ) /* not landmark, create a regular equation */
        {
	    double Ai_sum = 0;
	    double Bi_sum = 0;

            /* For each triangle containing that vertex */
            for( const vtkIdType& triId : m_CellsSharedByVertex[vertId] )
            {
                const vtkIdType i = vertId;
                vtkIdType j;
                vtkIdType k;


                this->Get_jk( m_InputMesh->GetCell(triId)->GetPointIds(), i, j, k );


                CalculateAB(i,j,k,triId,A,B);

                //double Ki, Li;
                //std::tie(Ki, Li) = CalculateKL(i,j,k,triId);
                
                const auto Ai = A(0); 
                const auto Aj = A(1); 
                const auto Ak = A(2); 

                const auto Bi = B(0); 
                const auto Bj = B(1); 
                const auto Bk = B(2); 

                const auto alpha1 = alphas(triId,0);
                const auto alpha2 = alphas(triId,1);
                const auto alpha3 = alphas(triId,2);
	
		//the sum of As and Bs has to be zero, verify
		if(Ai+Aj+Ak>0.001)
 		     throw std::logic_error("As in the triangle do not sum to 0");
 		if(Bi+Bj+Bk>0.001)
                {
                    std::cout<<"(i,j,k)="<<i<<","<<j<<","<<k<<std::endl;
                    std::cout<<"B(i,j,k)="<<Bi<<","<<Bj<<","<<Bk<<std::endl;
 		    throw std::logic_error("Bs in the triangle do not sum to 0");
                }
		Ai_sum += Ai;
		Bi_sum += Bi;

                /* Assemble one line of the system */
                L.coeffRef(vertId, i) += (Ai * (Ai*alpha1 + alpha2*Bi) + Bi * (alpha2*Ai + alpha3*Bi))*(2*m_Areas[triId]);
                L.coeffRef(vertId, j) += (Ai * (Aj*alpha1 + alpha2*Bj) + Bi * (alpha2*Aj + alpha3*Bj))*(2*m_Areas[triId]);
                L.coeffRef(vertId, k) += (Ai * (Ak*alpha1 + alpha2*Bk) + Bi * (alpha2*Ak + alpha3*Bk))*(2*m_Areas[triId]);            

#ifdef DEBUG_MESSAGES
    std::cout<<"Alphas: "<<alpha1<<", "<<alpha2<<", "<<alpha3<<std::endl;
#endif
                
                //L.coeffRef(vertId, i) += Ki * (Ai*alpha1 + alpha2*Bi) + Li * (alpha2*Ai + alpha3*Bi);
                //L.coeffRef(vertId, j) += Ki * (Aj*alpha1 + alpha2*Bj) + Li * (alpha2*Aj + alpha3*Bj);
                //L.coeffRef(vertId, k) += Ki * (Ak*alpha1 + alpha2*Bk) + Li * (alpha2*Ak + alpha3*Bk);            

                R(vertId, 0) = 0;
                R(vertId, 1) = 0;

            }
	
 	    //if(Ai_sum>0.1)
            //	throw std::logic_error("As around the vertex do not sum to 0");
 	    //if(Bi_sum>0.1)
            //	throw std::logic_error("Bs around the vertex do not sum to 0");

        }
 
    }
    
    std::cout<<std::endl;
    
    /* Add landmark constraints */
    for( vtkIdType i=0; i<m_LandmarkConstraints.size(); i++ )
    {
        const auto& lmk = m_LandmarkConstraints[i];

        //for (int k=0; k<npoints; k++)
        //    L.coeffRef(lmk.Id, k) = 0;
        
        L.coeffRef(lmk.Id, lmk.Id) = 1;
        R(lmk.Id, 0) = lmk.x;
        R(lmk.Id, 1) = lmk.y;
    }
    
    
    L.makeCompressed();
    
#ifdef DEBUG_MESSAGES
    std::cout<<"L: "<<L.toDense()<<std::endl;
    std::cout<<"R: "<<R<<std::endl;
#endif
    
    /* Solve the system */
    /* L has to be compressed and in Eigen::ColMajor ordering*/
    Eigen::SparseLU< Eigen::SparseMatrix<SysEqMatrixType::Scalar, Eigen::ColMajor, SysEqMatrixType::Index>, Eigen::COLAMDOrdering<SysEqMatrixType::Index> >   solver;

    std::cout<<"Solving the system: analyzePattern"<<std::endl;
    
    // Compute the ordering permutation vector from the structural pattern of A
    solver.analyzePattern(L);
    // Compute the numerical factorization 
    std::cout<<"Solving the system: factorize"<<std::endl;
    solver.factorize(L); 
    if(solver.info()!=Eigen::Success) {
        throw std::logic_error("Solver factorization failed");
    }    
    
    //Use the factors to solve the linear system 
    std::cout<<"Solving the system: solve"<<std::endl;
    DenseMatrixType result = solver.solve(R);     
    
    if(solver.info()!=Eigen::Success) {
        throw std::logic_error("Solving failed");
    }    
    
    x = result.col(0);
    y = result.col(1);
}


/* Ki = (hj-hi)*ctg(theta_k) + (hk-hi)*ctg(theta_j) */
/* Li = (gj-gi)*ctg(theta_k) + (gk-hi)*ctg(theta_j) */
//std::tuple<double, double> TMapUtils::CalculateKL(vtkIdType i, vtkIdType j, vtkIdType k, vtkIdType cellId)
//{
//    double vi[3], vj[3], vk[3];
//    this->m_InputMesh->GetPoint(i, vi);
//    this->m_InputMesh->GetPoint(j, vj);
//    this->m_InputMesh->GetPoint(k, vk);
//
//    const double hi = vi[1];
//    const double hj = vj[1];
//    const double hk = vk[1];
//    const double gi = vi[0];
//    const double gj = vj[0];
//    const double gk = vk[0];
//    
//    vtkIdType ilocal, jlocal, klocal;
//    std::tie(ilocal, jlocal, klocal) = VertexIdGlobal2Local(i,j,k,cellId);
//    
//    const double theta_k = this->m_Cotangents(cellId, klocal); /* angle cotangent on vertex k */
//    const double theta_j = this->m_Cotangents(cellId, jlocal); /* angle cotangent on vertex k */
//
//    //normalize the hj-hi by the length of the vector
//    const double vji_length = sqrt( (hj-hi)*(hj-hi)+(gj-gi)*(gj-gi) );
//    const double vki_length = sqrt( (hk-hi)*(hk-hi)+(gk-gi)*(gk-gi) );
//    
//    const double Ki = (hj-hi)*theta_k/vji_length + (hk-hi)*theta_j/vki_length; 
//    const double Li = (gj-gi)*theta_k/vji_length + (gk-hi)*theta_j/vki_length;
//    
//    return std::tuple<double, double>(Ki, Li);
//}




void TMapUtils::CalculateAB(vtkIdType i, vtkIdType j, vtkIdType k, vtkIdType cellId, Eigen::Vector3d& A, Eigen::Vector3d& B)
{
    double vi[3], vj[3], vk[3];
    this->m_InputMesh->GetPoint(i, vi);
    this->m_InputMesh->GetPoint(j, vj);
    this->m_InputMesh->GetPoint(k, vk);

    const double hi = vi[1];
    const double hj = vj[1];
    const double hk = vk[1];
    const double gi = vi[0];
    const double gj = vj[0];
    const double gk = vk[0];

    A(0) = (hj-hk)/(2*m_Areas[cellId]); //Ai
    A(1) = (hk-hi)/(2*m_Areas[cellId]); //Aj
    A(2) = (hi-hj)/(2*m_Areas[cellId]); //Ak

    B(0) = -(gj-gk)/(2*m_Areas[cellId]); //Bi
    B(1) = -(gk-gi)/(2*m_Areas[cellId]); //Bj
    B(2) = -(gi-gj)/(2*m_Areas[cellId]); //Bk
    
#ifdef DEBUG_MESSAGES
    std::cout<<"Area: "<<m_Areas[cellId]<<std::endl;
    std::cout<<"A: "<<A<<std::endl;
    std::cout<<"B: "<<B*(2*m_Areas[cellId])<<std::endl;
    std::cout<<"(gi, gj, gk)="<<gi<<","<<gj<<","<<gk<<std::endl;
#endif
    
}




void TMapUtils::Get_jk(vtkIdList* ijk, vtkIdType i, vtkIdType& j, vtkIdType& k)
{
    if( ijk->GetId(0) == i )
    {
        j = ijk->GetId(1);
        k = ijk->GetId(2);
    }
    else if( ijk->GetId(1) == i )
    {
        j = ijk->GetId(2);
        k = ijk->GetId(0);
    } 
    else
    {
        j = ijk->GetId(0);
        k = ijk->GetId(1);
    }
}



void TMapUtils::SetFixedBoundaryLandmarkConstraints()
{
    std::unordered_set<vtkIdType> edgePtIds;
    MeshEdgePointIds(edgePtIds);
    
    for( auto edgePtId : edgePtIds )
    {
        double *pt = m_InputMesh->GetPoint(edgePtId);
        m_LandmarkConstraints.push_back({edgePtId,pt[0],pt[1]});
    }
            

}



void TMapUtils::ZeroBC()
{
    const vtkIdType nTriangles = m_InputMesh->GetNumberOfCells();

    m_BC.resize(nTriangles);
    for( vtkIdType i=0; i<nTriangles; i++ )
        m_BC(i) = ComplexType(0,0);
}


bool TMapUtils::HasInputMeshChanged()
{
    return m_InputMeshModified || (m_InputMesh->GetMTime() != m_InputMeshMTime);
}

void TMapUtils::ResetModificationFlag()
{
    m_InputMeshModified=false;
    m_InputMeshMTime = m_InputMesh->GetMTime();
}



void TMapUtils::SmoothBC()
{
    //m_CellsSharedByCell
    
    /* For each triangle, get it's neighbors and average independently 
     the modulo and argument. The cell itself is not included */
    
    Eigen::VectorXd BC_mod(m_BC.size()); 
    Eigen::VectorXd BC_arg(m_BC.size()); 
           
    for( vtkIdType i=0; i<m_BC.size(); i++) //precalculate
    {
        BC_mod[i] = std::abs( m_BC[i] );
        BC_arg[i] = std::arg( m_BC[i] );
    }
    
    for( vtkIdType cell_id=0; cell_id<m_CellsAdjacent2Cell.size(); cell_id++)        
    {
        double mod = 0;
        double arg = 0;
        
        ComplexType t(0,0);
        
        const auto n_neighbors = m_CellsAdjacent2Cell[cell_id].size()-1; //exclude the cell cell_id
        
        for( const auto& neigb_cell_id : m_CellsAdjacent2Cell[cell_id])
            if(neigb_cell_id!=cell_id)
            {
                mod += BC_mod[neigb_cell_id]/n_neighbors;
                arg += BC_arg[neigb_cell_id]/n_neighbors;
            }
        
        m_BC[cell_id] = std::polar(mod, arg);
        
    }
    
}





void TMap::Update()
{
    //The complete t-map algorithm
    //Requires input: mesh, landmark constraints set beforehand
    double error=10;
    
    //create a working copy of the mesh
    vtkSmartPointer<vtkPolyData> mesh = vtkSmartPointer<vtkPolyData>::New();
    mesh->DeepCopy(m_InputMesh);
    
    
    //1. Set BC mu to 0    
    m_Utils.ZeroBC();    
    TMapUtils::BCVectorType BC_old = m_Utils.GetBC(); //create a copy
    
    int iter = 0;
    while( iter++ < m_MaxIter )
    {
        std::cout<<"Iteration "<<iter<<"; Error "<<error<<std::endl;
        
        //2. Calculate mesh f with mu 
        std::cout<<"Calculating mesh from BC"<<std::endl;
        m_Utils.MeshFromBC();

        char fname[100];
        sprintf(fname,"iter%04d.vtk",iter);
        vtkSmartPointer<vtkPolyDataWriter> wr = vtkSmartPointer<vtkPolyDataWriter>::New();
        wr->SetFileName(fname);
        wr->SetInputData(m_Utils.GetOutput());
        wr->Write();
        
        //3. Recalculate mu
        std::cout<<"Calculating BC"<<std::endl;
        m_Utils.CalculateBC( m_Utils.GetOutput() );

        TMapUtils::BCVectorType& BC = m_Utils.GetBC(); 
        error = (BC-BC_old).norm();
        BC_old = BC; //copy

        if( error<=m_Tolerance ) break; //probably not the most elegant solution, will change later on
        
        
        //4. Smooth mu
        std::cout<<"Smoothing BC"<<std::endl;        
        m_Utils.SmoothBC();
        
        //5. Reset mu's to constant, eq (5.26)
        m_Utils.SetConstantBC();
        
        //6. go to 2 until change in mu is not small
        
#ifdef DEBUG_MESSAGES
        std::cout<<"BC: "<<BC<<std::endl;
        std::cout<<"BCold: "<<BC_old<<std::endl;
#endif
        
    }
    
    
    m_OutputMesh = m_Utils.GetOutput();
        
}





void TMapUtils::SetConstantBC()
{   //eq (5.26)
    //calculate the mean abs BC 
    double mean_mod = 0;
    const auto ncells = m_BC.size();
    
    
    for(vtkIdType i=0; i<m_BC.size(); i++)
    {
        const double mod = std::abs(m_BC[i]);
        mean_mod += mod/m_BC.size();
        
        //normalize the BCs meanwhile
        m_BC[i] /= mod; 
    }

    //set BC to a constant magnitude
    for(vtkIdType i=0; i<m_BC.size(); i++)
    {
        m_BC[i] *= mean_mod; 
    }
    
}
