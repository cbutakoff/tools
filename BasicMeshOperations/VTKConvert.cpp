/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file
    \brief Convert vtk polydata between formats
 */
#include "VTKCommonTools.h"

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <string.h>
#include <stdio.h>
#include <iosfwd>
#include <string>

#include <iostream>
#include <fstream>

#include <vtkDataSetReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkShortArray.h>
#include <vtkUnstructuredGridWriter.h>

//#define LINEBREAK "\x0D\x0A"
#define LINEBREAK "\x0A"

//------------------------------------------------------------------
void SaveVolMeshVTK(const char* infile, const char* outfile_prefix, float scale);
void SaveVolMeshTetgen(const char* infile, const char* outfile_prefix, float scale);
void SaveVolMeshBSC(const char* infile, const char* outfile_prefix, float scale);
void SaveSurfMeshOFF(const char* infile, const char* outfile_prefix, float scale);
void SaveVolMesh2Mesh(const char* infile, const char* outfile_prefix, float scale, const char* array_name);
void SaveSurfMeshTetgen(const char* infile, const char* outfile_prefix, float scale);

int main(int argc, char *argv[]) {

    if (argc < 2) {
        std::cout << "VTK PolyData Converter v1.0 " << std::endl;
        std::cout << "Usage: VTKConvert inshape outshape scale [array name]" << std::endl;
        std::cout << "types determined from extension (.node - tetgen, .bsc - mesh for BSC)" << std::endl;
        std::cout << ".off - for surface meshes" << std::endl;
        std::cout << ".poly - for surface meshes, tetgen surface mesh" << std::endl;
        std::cout << ".mesh - for vol meshes. One more parameter is needed with array name for groups (must exist for both points and cells)" << std::endl;
        std::cout << ".vtkbin - for vol meshes. Saves as binary (to convert vtk asciii to vtk binary)" << std::endl;
        std::cout << "scale - scale factor to apply to points (e.g. 1)" << std::endl;
        return -1;
    }

    int c=1;
    char* inshape = argv[c++];
    char* outshape = argv[c++];
    float scale = atof(argv[c++]);


    char *ext;
    ext = outshape+strlen(outshape)-4;
    if(strcmp(ext,"node")==0)  //if we are requesting .node .element volumetric mesh
    {
        std::cout<<"Processing volumetric mesh for tetgen"<<std::endl;
        SaveVolMeshTetgen(inshape, outshape, scale);
    }
    else if(strcmp(ext,"poly")==0)  //if we are requesting .node .element volumetric mesh
    {
        std::cout<<"Generating surface tetgen mesh"<<std::endl;
        SaveSurfMeshTetgen(inshape, outshape, scale);
    }
    else if(strcmp(ext,".bsc")==0)
    {
        std::cout<<"Processing volumetric mesh for bsc"<<std::endl;
        SaveVolMeshBSC(inshape, outshape, scale);
    }
    else if(strcmp(ext,"kbin")==0)
    {
        std::cout<<"Converting UG vtk ascii to vtk bin"<<std::endl;
        SaveVolMeshVTK(inshape, outshape, scale);
    }
    else if(strcmp(ext,".off")==0)
    {
        std::cout<<"Processing volumetric mesh for off format"<<std::endl;
        SaveSurfMeshOFF(inshape, outshape, scale);
    }
    else if(strcmp(ext,".off")==0)
    {
        std::cout<<"Processing volumetric mesh for off format"<<std::endl;
        SaveSurfMeshOFF(inshape, outshape, scale);
    }
    else if(strcmp(ext,"mesh")==0)
    {
        std::cout<<"Processing volumetric mesh for .mesh format"<<std::endl;
        const char* array = argv[c++];
        std::cout<<"Array: "<<array<<std::endl;
        SaveVolMesh2Mesh(inshape, outshape, scale, array);
    }
    else //standard polygonal meshes
    {
        vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::Take(
            CommonTools::LoadShapeFromFile(inshape));

        CommonTools::SaveShapeToFile(pd, outshape);
    }

    return 0;
}



void SaveVolMeshVTK(const char* infile, const char* outfile_prefix, float scale)
{
    std::string outfile_nodes(outfile_prefix);
    std::string outfile; //output filename
    outfile = outfile_nodes.substr(0,outfile_nodes.length()-3);
    
    std::cout<<outfile<<std::endl;
    
    vtkSmartPointer<vtkDataSetReader> reader = vtkSmartPointer<vtkDataSetReader>::New();
    reader->SetFileName(infile);
    reader->Update();
    
    vtkSmartPointer<vtkUnstructuredGridWriter> wr = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    wr->SetInputData(reader->GetOutput());
    wr->SetFileTypeToBinary();
    wr->SetFileName(outfile.c_str());
    wr->Write();
}






void SaveVolMeshTetgen(const char* infile, const char* outfile_prefix, float scale)
{
    std::string outfile_nodes(outfile_prefix);
    std::string outfile_elements;
    outfile_elements = outfile_nodes.substr(0,outfile_nodes.length()-5) + ".ele";
    
    std::cout<<outfile_nodes<<std::endl;
    std::cout<<outfile_elements<<std::endl;
    
    vtkSmartPointer<vtkDataSetReader> reader = vtkSmartPointer<vtkDataSetReader>::New();
    reader->SetFileName(infile);
    reader->Update();
    
    vtkDataSet* volmesh = reader->GetOutput(); 
    
    std::ofstream node_file;
    node_file.open(outfile_nodes.c_str(),ios::trunc);
    node_file<< volmesh->GetNumberOfPoints() <<" 3 0 0"<<std::endl;
    
    for(int i=0; i<volmesh->GetNumberOfPoints(); i++)
    {
        double *pt = volmesh->GetPoint(i);
        node_file<<i+1<<" "<< pt[0]*scale<<" "<< pt[1]*scale<< " "<<pt[2]*scale<<std::endl;
    }
    node_file.close();
    
    
    std::ofstream ele_file;
    ele_file.open(outfile_elements.c_str(),ios::trunc);
    ele_file<< volmesh->GetNumberOfCells() <<" 4 0"<<std::endl;

    for(int i=0; i<volmesh->GetNumberOfCells(); i++)
    {
        vtkCell *cell = volmesh->GetCell(i);
        ele_file<<i+1<<" "<< cell->GetPointId(0)+1 <<" "<< cell->GetPointId(1)+1 << " "
                <<cell->GetPointId(2)+1<<" "<<cell->GetPointId(3)+1<<std::endl;
    }
    ele_file.close();
}


void SaveVolMeshBSC(const char* infile, const char* outfile_prefix, float scale)
{
    std::string outfile_nodes(outfile_prefix);
    std::string outfile_elements(outfile_prefix);
    std::string outfile_gradient(outfile_prefix);
    outfile_nodes += ".node";
    outfile_elements += ".ele";
    outfile_gradient += ".grad";
    
    std::cout<<outfile_nodes<<std::endl;
    std::cout<<outfile_elements<<std::endl;
    
    vtkSmartPointer<vtkDataSetReader> reader = vtkSmartPointer<vtkDataSetReader>::New();
    reader->SetFileName(infile);
    reader->Update();
    
    vtkDataSet* volmesh = reader->GetOutput(); 
    
    std::ofstream node_file;
    node_file.open(outfile_nodes.c_str(),ios::trunc);
    node_file<< "COORDINATES"<<std::endl;
    
    for(int i=0; i<volmesh->GetNumberOfPoints(); i++)
    {
        double *pt = volmesh->GetPoint(i);
        node_file<<i+1<<" "<< pt[0]*scale<<" "<< pt[1]*scale<< " "<<pt[2]*scale<<LINEBREAK;
    }

    node_file<< "END_COORDINATES"<<std::endl;
    node_file.close();
    
    
    std::ofstream ele_file;
    ele_file.open(outfile_elements.c_str(),ios::trunc);
    ele_file<< "ELEMENTS"<<std::endl;

    for(int i=0; i<volmesh->GetNumberOfCells(); i++)
    {
        vtkCell *cell = volmesh->GetCell(i);
        ele_file<<i+1<<" "<< cell->GetPointId(0)+1 <<" "<< cell->GetPointId(1)+1 << " "
                <<cell->GetPointId(3)+1<<" "<<cell->GetPointId(2)+1<<LINEBREAK; //flip the element for alya
    }

    ele_file<< "END_ELEMENTS"<<std::endl;
    ele_file.close();
    
    vtkFloatArray* grads = (vtkFloatArray*)volmesh->GetPointData()->GetArray("Fibers");
    if(grads!=NULL)
    {
        std::cout<<outfile_gradient<<std::endl;
        
        std::ofstream grad_file;
        grad_file.open(outfile_gradient.c_str(),ios::trunc);

        for(int i=0; i<grads->GetNumberOfTuples(); i++)
        {
            double *tuple = grads->GetTuple(i);
            double length =  sqrt( tuple[0]*tuple[0]+ tuple[1]*tuple[1]+ tuple[2]*tuple[2]);
            grad_file<< i+1<<" "<< tuple[0]/length <<" "<< tuple[1]/length << " " << tuple[2]/length<<LINEBREAK;
        }

        grad_file.close();
    }
}


void SaveSurfMeshTetgen(const char* infile, const char* outfile_prefix, float scale)
{
    std::cout<<outfile_prefix<<std::endl;
    
    vtkSmartPointer<vtkDataSetReader> reader = vtkSmartPointer<vtkDataSetReader>::New();
    reader->SetFileName(infile);
    reader->Update();
    
    vtkDataSet* mesh = reader->GetOutput(); 
    
    std::ofstream node_file;
    node_file.open(outfile_prefix,ios::trunc);
    node_file<< "# Part 1 - node list"<<std::endl;
    node_file<< "# node count, 3 dim, no attribute, no boundary marker"<<std::endl;
    
    node_file<< mesh->GetNumberOfPoints() <<" 3 0 0"<<std::endl;
    
    for(int i=0; i<mesh->GetNumberOfPoints(); i++)
    {
        double *pt = mesh->GetPoint(i);
        node_file<< i+1<<" "<<pt[0]*scale<<" "<< pt[1]*scale<< " "<<pt[2]*scale<<std::endl;
    }

    node_file<< "# Part 2 - facet list"<<std::endl;
    node_file<< "# facet count, no boundary marker"<<std::endl;
    node_file<< mesh->GetNumberOfCells() <<" 0"<<std::endl;

    for(int i=0; i<mesh->GetNumberOfCells(); i++)
    {
        vtkCell *cell = mesh->GetCell(i);
        if(cell->GetNumberOfPoints()>2)
        {
            node_file<< " 1            # 1 polygon, no hole, no boundary marker"<<std::endl;

            node_file<<cell->GetNumberOfPoints();
            for(int j=0; j<cell->GetNumberOfPoints(); j++)
            {
                    node_file << " "<< cell->GetPointId(j)+1;
            }

            node_file<<std::endl;
        }
    }
    
    node_file<<"# Part 3 - hole list"<<std::endl;
    node_file<<"0            # no hole"<<std::endl;
    node_file<<"# Part 4 - region list"<<std::endl;
    node_file<<"0            # no region"<<std::endl;
    
    node_file.close();
    
    
}




void SaveSurfMeshOFF(const char* infile, const char* outfile_prefix, float scale)
{
    std::cout<<outfile_prefix<<std::endl;
    
    vtkSmartPointer<vtkDataSetReader> reader = vtkSmartPointer<vtkDataSetReader>::New();
    reader->SetFileName(infile);
    reader->Update();
    
    vtkDataSet* volmesh = reader->GetOutput(); 
    
    std::ofstream node_file;
    node_file.open(outfile_prefix,ios::trunc);
    node_file<< "OFF"<<std::endl;
    node_file<< volmesh->GetNumberOfPoints() <<" "<<
            volmesh->GetNumberOfCells() <<" 0"<<std::endl;
    
    for(int i=0; i<volmesh->GetNumberOfPoints(); i++)
    {
        double *pt = volmesh->GetPoint(i);
        node_file<< pt[0]*scale<<" "<< pt[1]*scale<< " "<<pt[2]*scale<<std::endl;
    }

    for(int i=0; i<volmesh->GetNumberOfCells(); i++)
    {
        vtkCell *cell = volmesh->GetCell(i);
        if(cell->GetNumberOfPoints()>2)
        {
            node_file<<cell->GetNumberOfPoints();
            for(int j=0; j<cell->GetNumberOfPoints(); j++)
                    node_file << " "<< cell->GetPointId(j);

            node_file<<std::endl;
        }
    }
    
    node_file.close();
    
    
}

void SaveVolMesh2Mesh(const char* infile, const char* outfile_prefix, float scale, const char* array_name)
{
    std::string outfile_nodes(outfile_prefix);
    
    std::cout<<outfile_nodes<<std::endl;
    
    vtkSmartPointer<vtkDataSetReader> reader = vtkSmartPointer<vtkDataSetReader>::New();
    reader->SetFileName(infile);
    reader->Update();
    
    vtkDataSet* volmesh = reader->GetOutput(); 
    
    std::ofstream node_file;
    node_file.open(outfile_nodes.c_str(),ios::trunc);
    node_file<< "MeshVersionFormatted 2"<<std::endl;
    node_file<< "Vertices"<<std::endl;
    node_file<< volmesh->GetNumberOfPoints() <<std::endl;
    
    vtkShortArray *pscal = (vtkShortArray*)volmesh->GetPointData()->GetArray(array_name);
    vtkShortArray *cscal = (vtkShortArray*)volmesh->GetCellData()->GetArray(array_name);
    
    for(int i=0; i<volmesh->GetNumberOfPoints(); i++)
    {
        double *pt = volmesh->GetPoint(i);
        node_file<< pt[0]*scale<<" "<< pt[1]*scale<< " "<<pt[2]*scale<<" "<<pscal->GetValue(i)<<std::endl;
    }
    
    
    node_file<< "Tetrahedra"<<std::endl;
    node_file<< volmesh->GetNumberOfCells() <<std::endl;

    for(int i=0; i<volmesh->GetNumberOfCells(); i++)
    {
        vtkCell *cell = volmesh->GetCell(i);
        node_file<< cell->GetPointId(0)+1 <<" "<< cell->GetPointId(1)+1 << " "
                <<cell->GetPointId(2)+1<<" "<<cell->GetPointId(3)+1<<" "<<cscal->GetValue(i)<<std::endl;
    }

    node_file<< "End"<<std::endl;

    node_file.close();
}
