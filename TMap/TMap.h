/* Implmenetation of Teichmuller map from
 * LM Lui, Teichmuller mapping (T-Map) and its applications to landmark matching registration
 * SIAM J IMG SCI 7(1) 2014
 */

#ifndef __TMap_h
#define __TMap_h

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <complex>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <stdexcept>

#include <unordered_set>
#include <vector>


/*
 * Steps:
 * 1. SetInput
 * 2. Precalculate
 * 3. SetLandmarks
 * 4. LandmarkConstrainedTMap
 * If anything happens throws std::logic_error
 */
class TMapUtils
{
public:
    typedef struct __landmark {
        vtkIdType Id; /* Id of the vertex for which the coordinates are specified*/
        double x; /* Target x,y coordinates */
        double y;
    } LandmarkType;    
    
    typedef std::complex<double> ComplexType;
    typedef Eigen::VectorXcd BCVectorType;
    typedef Eigen::MatrixXd DenseMatrixType;
    typedef Eigen::Matrix2d DenseMatrix2x2Type;
    typedef Eigen::SparseMatrix<unsigned char, Eigen::ColMajor, vtkIdType> AdjacencyMatrixType;
    typedef std::vector<std::vector<vtkIdType> > NeighborIDMatrixType; 
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor, vtkIdType> SysEqMatrixType;     
    
    /* Input mesh has to be flat, only in XY plane, Z coordinate will be ignored */
    void SetInputMesh(vtkPolyData* mesh) { m_InputMesh = mesh; m_InputMeshModified = true;};

    /* Calculate Beltrami coefficients for each triangle 
     * calculates the BC with f given via 2 meshes. f maps the mesh from SetInputMesh to
     * tgtMesh passed as a prameter
     */
    void CalculateBC(vtkPolyData* tgtMesh); 
    BCVectorType& GetBC() {return m_BC; };
    void ZeroBC();
    void SmoothBC();
    void SetConstantBC();
    
    /* Reconstruct the mesh from the stored BCs
     * the mesh will have the same topology as the input mesh
     */
    void MeshFromBC(); 
    

    /* Build adjacency matrix, triangle areas. Make sure it's called if the 
     * mesh has been modified. You can call SetInputMesh() to force execution 
     * of this method since it checks m_Modified flag */
    void Precalculate();


    
    void SetFixedBoundaryLandmarkConstraints();
    void ClearLandmarkConstraints() {m_LandmarkConstraints.clear();};
    void AddLandmarkConstraint(vtkIdType ptId, double x, double y) {m_LandmarkConstraints.push_back({ptId,x,y});};

    vtkPolyData* GetOutput() {return m_OutputMesh;};
    
    
    TMapUtils():m_InputMeshModified(false) {};
    ~TMapUtils() {};
protected:

    bool HasInputMeshChanged();
    void ResetModificationFlag();
    
    /* Find edge of the mesh */
    void MeshEdgePointIds(std::unordered_set<vtkIdType>& edgePtIds);
    
    /* Throw exception if inputs are not specified. Stores number of edges >=0
     * For each triangle containing an edge add 1 
     * Only upper triangle is used
     */
    void ValidateInputs() const;
    
    void StopOnError(bool condition, const char* message) {if(condition) throw message;};
    
    void CalculateAlphasFromBC(TMapUtils::DenseMatrixType& alpha);
    
    void CalculateVerticesFromBC(Eigen::VectorXd& x, Eigen::VectorXd& y); /* Calculate vertex positions x,y from m_BC including the landmark constraints */

    /* Return the angle formed by vectors v0v1 and v0v2*/
    //double GetAngle(const double* v0, const double* v1, const double* v2) const;
    
    /* For vertex indices i,j,k, returns Ai, Aj, Ak*/
    void CalculateAB(vtkIdType i, vtkIdType j, vtkIdType k, vtkIdType cellId, Eigen::Vector3d& A, Eigen::Vector3d& B);
    
    //std::tuple<double, double> CalculateKL(vtkIdType i, vtkIdType j, vtkIdType k, vtkIdType cellId);    
    
    /*Transforms global vertex IDs to local vertext ids inside the triangle = 0,1,2*/
    std::tuple<vtkIdType, vtkIdType, vtkIdType> VertexIdGlobal2Local(vtkIdType i, vtkIdType j, vtkIdType k, vtkIdType cellId);
    
    /* from the array of 3 triangle vertex indices, given vertex id i, return indices j and k conserving order 
     * e.g. ijk = [1 2 3], i=2 return j = 3, k= 1
     */
    void Get_jk(vtkIdList* ijk, vtkIdType i, vtkIdType& j, vtkIdType& k);
    
    
private:

    
    vtkSmartPointer<vtkPolyData> m_OutputMesh;
    vtkSmartPointer<vtkPolyData> m_InputMesh;
    vtkMTimeType m_InputMeshMTime;
    
    
    BCVectorType m_BC; /* Beltrami coefficients */
    
    AdjacencyMatrixType m_AdjacencyMatrix; /*only upper triangle is used*/
    Eigen::VectorXd m_Areas; /* Triangle areas */
    //Eigen::MatrixXd m_Cotangents; /* Nx3 matrix, of angle cotangents, in the same order as vertices*/

    
    std::vector<LandmarkType> m_LandmarkConstraints; /*n x 2 matrix. */
    
    
    NeighborIDMatrixType m_CellsSharedByVertex; /* i-th element is a list of triangles containing the vertex i*/
    NeighborIDMatrixType m_CellsAdjacent2Cell; /* i-th element is a list of triangles sharing an edge with the cell i*/
    
    bool m_InputMeshModified;
};



class TMap
{
public:
    void SetInput(vtkPolyData* mesh) { m_InputMesh = mesh; m_Utils.SetInputMesh(m_InputMesh);};
    vtkPolyData* GetOutput() {return m_OutputMesh;};

    void Update();
    
    
    void AddFixedBoundaryLandmarkConstraints() {m_Utils.SetFixedBoundaryLandmarkConstraints(); };
    void ClearLandmarkConstraints() {m_Utils.ClearLandmarkConstraints(); };
    void AddLandmarkConstraint(vtkIdType ptId, double x, double y) {m_Utils.AddLandmarkConstraint(ptId,x,y);};
    
    void SetMaxNumberOfIterations(int niter) {m_MaxIter = niter;};
    
    TMap():m_MaxIter(100), m_Tolerance(1e-6) {};
    ~TMap() {};
protected:
private:
    vtkSmartPointer<vtkPolyData> m_OutputMesh;
    vtkSmartPointer<vtkPolyData> m_InputMesh;
    
    TMapUtils m_Utils;
    int m_MaxIter;
    double m_Tolerance;

};





#endif
