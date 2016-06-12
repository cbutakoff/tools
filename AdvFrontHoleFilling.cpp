//The advancing‐front mesh generation method revisited - ReadCube
#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>

#include <vector>


typedef double PointType[3];

typedef struct __edge{
    vtkIdType v0; //vertices
    vtkIdType v1;
    PointType n0; //normals at the vertices
    PointType n1;
} EdgeType;

typedef std::vector<PointType> HoleBoundaryType; 
    
//returns the unordered boundary of all the holes
void FindHoles(vtkPolyData *mesh, HoleBoundaryType& boundary);

int main(int argc, char **argv)
{
    const char *filename = argv[1];
    
    vtkSmartPointer<vtkPolyDataReader> rdr = vtkSmartPointer<vtkPolyDataReader>::New();
    rdr->SetFileName(filename);
    rdr->Update();
    
    vtkPolyData *mesh = rdr->GetOutput();
    
    HoleBoundaryType boundary; 
    FindHoles(mesh, boundary);
    
    return 0;
}

void FindHoles(vtkPolyData *mesh, HoleBoundaryType& boundary)
{
    boundary.clear();
    
    //iterate over every triangle and edge
}
