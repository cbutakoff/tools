
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>

#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkType.h>

#include <cgnslib.h>

#include <cmath>

using namespace std;

void write_cgns(vtkUnstructuredGrid *mesh, const char* filename) {
    cout<<"Storing tetrahedral mesh in CGNS format.."<<endl;
    
    int64_t nr_nodes = mesh->GetNumberOfPoints();
    int64_t nr_elem = mesh->GetNumberOfCells();

    cgsize_t size[3] = {nr_nodes, nr_elem, 0};

    double* xc = new double[nr_nodes];
    double* yc = new double[nr_nodes];
    double* zc = new double[nr_nodes];
    cgsize_t* tets = new cgsize_t[nr_elem * 4];

    cout<<"Extracting nodes: "<<nr_nodes<<flush<<endl;
    for (uint64_t i = 0; i < nr_nodes; ++i) {
        if( i%100000 == 0 )
            cout<<"Progress "<<i*100/nr_nodes<<"%\r"<<flush;
        double pt[3];
        mesh->GetPoint(i, pt);
        xc[i] = pt[0];
        yc[i] = pt[1];
        zc[i] = pt[2];
    }
    cout<<endl;

    cout<<"Extracting elements "<<nr_elem<<endl<<flush;
    for (size_t i = 0, j = 0; i < nr_elem; ++i, j += 4) {
        if( i%100000 == 0 )
            cout<<"Progress "<<i*100/nr_elem<<"%\r"<<flush;

        vtkCell* cell = mesh->GetCell(i);
        if(cell->GetNumberOfPoints()!=4) 
        {
            cout<<"non tetra id:"<<i<<endl;
            throw;
        }
        tets[j + 0] = cell->GetPointId(0) + 1;
        tets[j + 1] = cell->GetPointId(1) + 1;
        tets[j + 2] = cell->GetPointId(2) + 1;
        tets[j + 3] = cell->GetPointId(3) + 1;
    }
    cout<<endl;

     
    int fn, bn, zn, cn, sn;
    if (cg_open(filename, CG_MODE_WRITE, &fn)) {
        cg_error_exit();
    }
    if (cg_base_write(fn, "Base", 3, 3, &bn)) {
        cg_error_exit();
    }
    if (cg_zone_write(fn, bn, "Zone", size, Unstructured, &zn)) {
        cg_error_exit();
    }
    if (cg_coord_write(fn, bn, zn, RealDouble, "CoordinateX", xc, &cn)) {
        cg_error_exit();
    }
    if (cg_coord_write(fn, bn, zn, RealDouble, "CoordinateY", yc, &cn)) {
        cg_error_exit();
    }
    if (cg_coord_write(fn, bn, zn, RealDouble, "CoordinateZ", zc, &cn)) {
        cg_error_exit();
    }

    if (cg_section_write(fn, bn, zn, "Tetra", CGNS_ENUMV(TETRA_4), 1,  nr_elem, 0,  tets, &sn)) {
        cg_error_exit();
    }

    if (cg_close(fn)) {
        cg_error_exit();
    }


    delete xc;
    delete yc;
    delete zc;
    delete tets;

    printf(" done\n");
    
}


int main(int argc, char *argv[])
{
    const char* vtufile = argv[1];
    const char* cgnsfile = argv[2];

    cout<<"Reading "<<vtufile<<endl;
    vtkSmartPointer<vtkXMLUnstructuredGridReader> rd =   vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    rd->SetFileName(vtufile);
    rd->Update();

    cout<<"calling writer"<<flush<<endl;
    write_cgns(rd->GetOutput(), cgnsfile);

    return 0;
}