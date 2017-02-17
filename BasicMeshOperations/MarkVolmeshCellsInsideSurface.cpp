// VTK headers
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkCell.h>
#include <vtkIdList.h>
#include <vtkCellData.h>
#include <vtkIntArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkCellCenters.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkPointData.h>

#include "VTKCommonTools.h"

void ReadUnstructuredGrid(vtkUnstructuredGrid* outImage, const char *fileName);
void WriteUnstructuredData(vtkUnstructuredGrid* in_shape, const char *fileName);
void ReadPolydata(vtkPolyData* out_shape, const char *fileName);



int main(int argc, char** argv) {
    
    if (argc != 5 ) {
        std::cout << " Usage: " << std::endl;
        std::cout << argv[0] << " vol_mesh.vtk surf_mesh.vtk output_vol_mesh.vtk elements.txt" << std::endl;
        std::cout << "surf_mesh.vtk must be a closed surface" << std::endl;
        return 0;
    }
   
    const char* input1 = argv[1];
    const char* input2 = argv[2];
    const char* output = argv[3];
    const char* elements = argv[4];
    
    // -------------------------------------------------------------------------
    // Reading input poly data
    vtkSmartPointer<vtkUnstructuredGrid> volume = vtkSmartPointer<vtkUnstructuredGrid>::New();
    ReadUnstructuredGrid(volume, input1);
    
    vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::Take(
        CommonTools::LoadShapeFromFile( input2 ) );
    
    // -------------------------------------------------------------------------
    
    vtkSmartPointer<vtkIntArray> isIn = vtkSmartPointer<vtkIntArray>::New();
    isIn -> SetName("Is_within_distance");
    
    unsigned int ncells = volume->GetNumberOfCells();
    isIn->SetNumberOfValues(ncells);
    
    for (unsigned int i = 0; i < ncells; i++)
    {
        isIn -> SetValue(i, 0);
    }
    volume -> GetCellData() -> AddArray(isIn);

    // -------------------------------------------------------------------------
    
    std::ofstream file;
    file.open(elements);


    vtkSmartPointer<vtkCellCenters> centers = vtkSmartPointer<vtkCellCenters>::New();
    centers->SetInputData(volume);
    centers->Update();
    
    
    vtkSmartPointer<vtkSelectEnclosedPoints> inside_slector = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
    inside_slector->SetInputData(centers->GetOutput());
    inside_slector->SetSurfaceData(surface);
    inside_slector->CheckSurfaceOn();    
    inside_slector->Update();

    
    
    for (int i = 0; i < volume->GetNumberOfCells(); i++)
    {
        isIn->SetValue(i, inside_slector->IsInside(i));
    }
    
    centers->GetOutput()->GetPointData()->AddArray(isIn);
    CommonTools::SaveShapeToFile(centers->GetOutput(), "centers.vtk");

    file.close();
    
    WriteUnstructuredData(volume, output);
    
}



// -----------------------------------------------------------------------------
// Read poly data
void ReadPolydata(vtkPolyData* out_shape, const char *fileName)
{
    std::cout << "Reading vtk polydata file: " << fileName << " ... " << std::endl;
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader -> SetFileName (fileName);
    reader -> Update();
    out_shape->ShallowCopy(reader->GetOutput());
    std::cout << "Vtk poly data file has been read. " << std::endl;
    std::cout << "Number of vertices in mesh: " << out_shape->GetNumberOfPoints() << std::endl;
}
// -----------------------------------------------------------------------------
// Unstructured grid reader
void ReadUnstructuredGrid(vtkUnstructuredGrid* outImage, const char *fileName)
{
    std::cout << "Reading vtk unstructured grid data from file: " << fileName << " ... " << std::endl;
    vtkUnstructuredGridReader* reader = vtkUnstructuredGridReader::New();
    reader->SetFileName(fileName);
    reader->Update();
    outImage->ShallowCopy(reader->GetOutput());
    std::cout << "Vtk unstructured grid data file has been read." << std::endl; 
}
// -----------------------------------------------------------------------------
// Unstructured data writer
void WriteUnstructuredData(vtkUnstructuredGrid* in_shape, const char *fileName)
{
    std::cout << "Writing unstructured grid data to file: " << fileName << " ... " << std::endl;
    vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(in_shape);
#else
    writer->SetInputData(in_shape);
#endif
    writer -> SetFileName (fileName);
    //writer->SetFileTypeToBinary();
    writer->Write();
    std::cout << "Unstructured data written to file: " << fileName<< " ." << std::endl;  
}
