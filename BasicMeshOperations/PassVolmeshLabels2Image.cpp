#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkType.h>
#include <vtkCellLocator.h>
#include <vtkDataArray.h>
#include <vtkFieldData.h>
#include <vtkCellData.h>
#include <vtkCell.h>

//return false if not found

int main(int argc, char **argv)
{
    if (argc < 3) {
        std::cerr << "Usage: PassVolmeshLabels2Image inputImageFile inputMeshFile outputImageFile" << std::endl;
        return EXIT_FAILURE;
    }

    const int outside_label = 0; //to label outside pixels
    const char* label_array = "Material";
    
    //parsing parameters
    int c = 1;
    const char* imagefilename = argv[c++];
    const char* meshfilename = argv[c++];
    const char* outputfilename = argv[c++];

    
    std::cout<<"Loading mesh "<<meshfilename<<std::endl;
    vtkSmartPointer<vtkXMLUnstructuredGridReader> rdr = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    rdr->SetFileName(meshfilename);
    rdr->Update();
    
    vtkUnstructuredGrid *mesh = (vtkUnstructuredGrid *)rdr->GetOutput();    
    
    
}



