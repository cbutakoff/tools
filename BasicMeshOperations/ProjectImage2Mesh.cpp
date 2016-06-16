/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file 
    \brief Given an image and a surface mesh in the same coordinate system, sample image along the normals
 */
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataNormals.h>

#include <vtkDataSetReader.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkProbeFilter.h>


#include <vtkFloatArray.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCell.h>
#include <vtkCellCenters.h>

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

vtkFloatArray* SamplePoints(double start, double end, double step, vtkPoints* points, vtkDataArray* normals, vtkImageData* image,
        const char* spositions_filename = NULL);


//samples_per_pt - number of samples per point, samples should be arranged as tuples
void StoreSamplesCSV(vtkDataArray* samples, const char* filename);

int main(int argc, char **argv) {
    std::cout << "Sample image around the points or cells of the mesh, along normal. Version 0.9" << std::endl;

    if (argc < 2) {
        std::cout << "Params: " << std::endl;
        std::cout << "-image <image.vtk> \t\t- image" << std::endl;
        std::cout << "-imesh <mesh.vtk> \t\t- input surface mesh" << std::endl;
        std::cout << "-omesh <mesh.vtk> \t\t- output surface mesh" << std::endl;
        std::cout << "-sp <mesh.vtk> \t\t- [optional] outputs sampling points" << std::endl;
        std::cout << "-csv <samples.csv> \t\t- [optional] store the samples in a matrix"<<std::endl; 
        std::cout << "-check_normals \t\t- try to impose that normals point inwards" << std::endl;
        std::cout << "-points|cells \t\t- sample for points or cells" << std::endl;
        std::cout << "-range <float> <float> <float>\t\t- initial and final distance for sampling (can be negative), step" << std::endl;
        std::cout << "Example parameters:-image image.vtk -imesh mesh.vtk -cells -range -2 2 0.5 -omesh out.vtk -- samples the "
                "image image.vtk along the normals to cell centers, to both sides of the surface. The result is stored in the out.vtk in float vector array 'Samples'" << std::endl;
        return -1;
    }

    char *image_filename = NULL;
    char *imesh_filename = NULL;
    char *omesh_filename = NULL;
    char *spositions_filename = NULL; //for sampling points
    char *csv_filename = NULL; 

    bool check_normals = false;

    typedef enum __sampling_location {
        slPoints, slCells
    } TypeSamplingLocation;
    TypeSamplingLocation sampling = slPoints;


    //sampling start, end and step
    float start = 1.0;
    float end = 1.0;
    float step = 1.0;


    //basic parameter parsing
    for (int c = 1; c < argc; c++) {
        if (strcmp(argv[c], "-image") == 0) {
            image_filename = argv[++c];
        } else if (strcmp(argv[c], "-imesh") == 0) {
            imesh_filename = argv[++c];
        } else if (strcmp(argv[c], "-omesh") == 0) {
            omesh_filename = argv[++c];
        } else if (strcmp(argv[c], "-sp") == 0) {
            spositions_filename = argv[++c];
        } else if (strcmp(argv[c], "-csv") == 0) {
            csv_filename = argv[++c];           
        } else if (strcmp(argv[c], "-check_normals") == 0) {
            check_normals = true;
        } else if (strcmp(argv[c], "-points") == 0) {
            sampling = slPoints;
        } else if (strcmp(argv[c], "-cells") == 0) {
            sampling = slCells;
        } else if (strcmp(argv[c], "-range") == 0) {
            start = atof(argv[++c]);
            end = atof(argv[++c]);
            step = atof(argv[++c]);
        }
    }

    //Output the parameters
    std::cout << "The following parameters will be used:" << std::endl;
    std::cout << "Image:" << image_filename << std::endl;
    std::cout << "Input mesh:" << imesh_filename << std::endl;
    std::cout << "Output mesh:" << omesh_filename << std::endl;
    std::cout << "Normal checking:" << (check_normals ? "enabled" : "disabled") << std::endl;
    std::cout << "Sampling for:" << (sampling == slPoints ? "points" : sampling == slCells ? "cells" : "unrecognized") << std::endl;
    std::cout << "Sampling from " << start << " to " << end << " with step " << step << std::endl;
    
    if( spositions_filename != NULL)
        std::cout << "Sampled points saved in:" << spositions_filename << std::endl;
        
    if(csv_filename != NULL )  
        std::cout << "Samples saved in ASCII in:" << csv_filename << std::endl;

    
    if (image_filename == NULL || imesh_filename == NULL || omesh_filename == NULL) {
        cout << "Filenames are not specified. Aborting." << std::endl;
        return -1;
    }

    if( step==0 ){
        cout << "Step cannot be 0." << std::endl;
        return -1;        
    }
    
    //read the image
    std::cout<<"Reading image"<<std::endl;
    vtkSmartPointer<vtkDataSetReader> imrdr = vtkSmartPointer<vtkDataSetReader>::New();
    imrdr->SetFileName(image_filename);
    imrdr->Update();
    vtkImageData* image = dynamic_cast<vtkImageData*> (imrdr->GetOutput());

    //read the input mesh
    std::cout<<"Reading mesh"<<std::endl;
    vtkSmartPointer<vtkPolyDataReader> meshrdr = vtkSmartPointer<vtkPolyDataReader>::New();
    meshrdr->SetFileName(imesh_filename);
    meshrdr->Update();
    vtkPolyData* imesh = meshrdr->GetOutput();

    //generate normals
    std::cout<<"Generating normals"<<std::endl;
    vtkSmartPointer<vtkPolyDataNormals> norm_gen = vtkSmartPointer<vtkPolyDataNormals>::New();
    norm_gen->SetInputData(imesh);
    norm_gen->SplittingOff();
 //   norm_gen->ComputeCellNormalsOn();
    norm_gen->ComputePointNormalsOn();
    norm_gen->Update();

    vtkDataArray* pnormals = norm_gen->GetOutput()->GetPointData()->GetNormals();
//    vtkDataArray* cnormals = norm_gen->GetOutput()->GetCellData()->GetNormals();

    //check normals of requested
    if (check_normals) {
        std::cout<<"Checking normals"<<std::endl;
        
        //find centroid
        double c[] = {0, 0, 0};
        for (int i = 0; i < imesh->GetNumberOfPoints(); i++) {
            double pt[3];
            imesh->GetPoint(i, pt);
            c[0] += pt[0];
            c[1] += pt[1];
            c[2] += pt[2];
        }
        c[0] /= imesh->GetNumberOfPoints();
        c[1] /= imesh->GetNumberOfPoints();
        c[2] /= imesh->GetNumberOfPoints();

        //take one normal and check if it points inwards
        double n[3];
        pnormals->GetTuple(0, n);
        double pt[3];
        imesh->GetPoint(0, pt);
        //dot product of c-pt and n
        double dp = (c[0] - pt[0]) * n[0] + (c[1] - pt[1]) * n[1] + (c[2] - pt[2]) * n[2];
        
//        std::cout<<"c: "<<c[0]<<" "<<c[1]<<" "<<c[2]<<std::endl;
//        std::cout<<"pt: "<<pt[0]<<" "<<pt[1]<<" "<<pt[2]<<std::endl;
//        std::cout<<"n: "<<n[0]<<" "<<n[1]<<" "<<n[2]<<std::endl;
//        std::cout<<"dp: "<<dp<<std::endl;
        
        if (dp < 0) //normals point outside
        {
            std::cout<<"Flipping normals required"<<std::endl;
            norm_gen->FlipNormalsOn();
        }
    }


    //create sampling positions
    if (sampling == slPoints) {
        std::cout<<"Sampling for points"<<std::endl;
        norm_gen->ComputePointNormalsOn();
        norm_gen->Update();
        pnormals = norm_gen->GetOutput()->GetPointData()->GetNormals();


        vtkSmartPointer<vtkFloatArray> samples =
                vtkSmartPointer<vtkFloatArray>::Take(
                SamplePoints(start, end, step, imesh->GetPoints(), pnormals,
                image, spositions_filename));

        samples->SetName("Samples");
        imesh->GetPointData()->AddArray(samples);
        //save samples if required
        if(csv_filename!=NULL) StoreSamplesCSV(samples, csv_filename);
        
    } else if (sampling == slCells) {
        //sampling positions
        norm_gen->ComputeCellNormalsOn();
        norm_gen->Update();
        vtkDataArray* cnormals = norm_gen->GetOutput()->GetCellData()->GetNormals();

        std::cout<<"Sampling for cells"<<std::endl;
        vtkSmartPointer<vtkCellCenters> centers = vtkSmartPointer<vtkCellCenters>::New();
        centers->SetInputData(imesh);
        centers->Update();

        vtkSmartPointer<vtkFloatArray> samples =
                vtkSmartPointer<vtkFloatArray>::Take(
                SamplePoints(start, end, step, centers->GetOutput()->GetPoints(), cnormals,
                image, spositions_filename));


        samples->SetName("Samples");
        imesh->GetCellData()->AddArray(samples); //add data to the cells

        //save samples if required
        if(csv_filename!=NULL) StoreSamplesCSV(samples, csv_filename);
    }

    
    
    //save the result
    std::cout<<"Writing results"<<std::endl;
    vtkSmartPointer<vtkPolyDataWriter> pdwr = vtkSmartPointer<vtkPolyDataWriter>::New();
    pdwr->SetFileName(omesh_filename);
    pdwr->SetFileTypeToBinary();
    pdwr->SetInputData(imesh);
    pdwr->Write();


    return 0;
}


//given starting, ending positions, step; sampling centers and normals
vtkFloatArray* SamplePoints(double start, double end, double step, vtkPoints* points, vtkDataArray* normals, vtkImageData* image,
        const char* spositions_filename) {
    vtkSmartPointer<vtkPoints> pos = vtkSmartPointer<vtkPoints>::New();

    //calculate number of elements in a tuple. Can be done using one formula, 
    //i'm doing with a loop just to be sure.
    int n_elements = 0;
    for (double j = start; j <= end; j += step) n_elements++;
    
    //pos->SetNumberOfPoints(n_elements*points->GetNumberOfPoints());
    
    std::cout<<"Creating sampling positions"<<std::endl;
    for (int i = 0; i < points->GetNumberOfPoints(); i++) {
        double pt[] = {0, 0, 0};
        points->GetPoint(i, pt);

        double n[3];
        normals->GetTuple(i, n); //also take cell normals

        for (double j = start; j <= end; j += step) {
            double pt1[3];
            pt1[0] = pt[0] + j * n[0];
            pt1[1] = pt[1] + j * n[1];
            pt1[2] = pt[2] + j * n[2];
            pos->InsertNextPoint(pt1);
            
//            std::cout<<"pt:"<<pt[0]<<" "<<pt[1]<<" "<<pt[2]<<std::endl;
//            std::cout<<"n:"<<n[0]<<" "<<n[1]<<" "<<n[2]<<std::endl;
//            std::cout<<"pt1:"<<pt1[0]<<" "<<pt1[1]<<" "<<pt1[2]<<std::endl;
//            char aa[20];
//            std::cin>>aa; 
       }
    }

    //create a fake polydata
    vtkSmartPointer<vtkPolyData> pos_pd = vtkSmartPointer<vtkPolyData>::New();
    pos_pd->SetPoints(pos);



    //sample the image
    std::cout<<"Calling probe"<<std::endl;
    vtkSmartPointer<vtkProbeFilter> probe = vtkSmartPointer<vtkProbeFilter>::New();
    probe->SetInputData(pos_pd);
    probe->SetSourceData(image);
    probe->Update();
    vtkDataArray* samples = probe->GetOutput()->GetPointData()->GetScalars();

    //fill the original polydata with the vectors containing the samples
    vtkFloatArray* vaSamples = vtkFloatArray::New();


    vaSamples->SetNumberOfComponents(n_elements);
    vaSamples->SetName("Samples");

    
    std::cout<<"Storing samples"<<std::endl;
    unsigned long linear_index = 0;
    for (int i = 0; i < points->GetNumberOfPoints(); i++) //iterate over cells
    {
        double tuple[n_elements];
        int k = 0;
        for (double j = start; j <= end; j += step, k++) {
            tuple[k] = samples->GetTuple1(linear_index++);
        }

        vaSamples->InsertNextTuple(tuple);
    }

    
    if (spositions_filename != NULL) {
        pos_pd->GetPointData()->AddArray(samples);
        
        vtkSmartPointer<vtkPolyDataWriter> pdwr = vtkSmartPointer<vtkPolyDataWriter>::New();
        pdwr->SetFileName(spositions_filename);
        pdwr->SetInputData(pos_pd);
        pdwr->SetFileTypeToBinary();
        pdwr->Write();
    }
    

    return vaSamples;
}


void StoreSamplesCSV(vtkDataArray* samples, const char* filename)
{
    std::ofstream file(filename);
    for (int i = 0; i < samples->GetNumberOfTuples(); i++) //iterate over cells
    {
        int j;
        for (j = 0; j < samples->GetNumberOfComponents()-1; j++) {
            file<<samples->GetTuple(i)[j]<<", ";
        }
        file<<samples->GetTuple(i)[j]<<std::endl;
    }
}
