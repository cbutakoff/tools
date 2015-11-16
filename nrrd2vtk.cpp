/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file 
    \brief Convert nrrd unsigned int images to vtk. Works both ways
    */
#include <vtkImageData.h>

#include <vtkSmartPointer.h>

#include "CommonTools.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>


#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>


int main(int argc, char **argv)
{
	std::cout << "nrrd2vtk. Version 1.3. Works bboth ways."<< std::endl;

	if(argc<4)
	{
		std::cout << "nrrd2vtk input.nrrd output.vtk uint|short|uchar"<< std::endl;		
		return 0;
	}

	const char* inp = argv[1];
	const char* outp = argv[2];
	const char* type = argv[3];
	
	if(strcmp(type, "uint")==0)
	{
		typedef itk::Image< unsigned int, 3 >         ImageType;	
		typedef ImageType::Pointer	    ImageTypePointer;


		typedef itk::ImageFileReader<ImageType> FileReaderType;
		FileReaderType::Pointer reader = FileReaderType::New();
		reader->SetFileName(inp);
		reader->Update();
	

		typedef itk::ImageFileWriter<ImageType> FileWriterType;
		FileWriterType::Pointer writer = FileWriterType::New();
		writer->SetFileName(outp);
		writer->SetInput( reader->GetOutput() );
		writer->Update();
	}
	else if(strcmp(type, "short")==0)
	{
		typedef itk::Image< short, 3 >         ImageType;	
		typedef ImageType::Pointer	    ImageTypePointer;


		typedef itk::ImageFileReader<ImageType> FileReaderType;
		FileReaderType::Pointer reader = FileReaderType::New();
		reader->SetFileName(inp);
		reader->Update();
	

		typedef itk::ImageFileWriter<ImageType> FileWriterType;
		FileWriterType::Pointer writer = FileWriterType::New();
		writer->SetFileName(outp);
		writer->SetInput( reader->GetOutput() );
		writer->Update();
	}
	else if(strcmp(type, "uchar")==0)
	{
		typedef itk::Image< unsigned char, 3 >         ImageType;	
		typedef ImageType::Pointer	    ImageTypePointer;


		typedef itk::ImageFileReader<ImageType> FileReaderType;
		FileReaderType::Pointer reader = FileReaderType::New();
		reader->SetFileName(inp);
		reader->Update();
	

		typedef itk::ImageFileWriter<ImageType> FileWriterType;
		FileWriterType::Pointer writer = FileWriterType::New();
		writer->SetFileName(outp);
		writer->SetInput( reader->GetOutput() );
		writer->Update();
	}


/*	ConnectorType::Pointer connector = ConnectorType::New();
	connector->SetInput(reader->GetOutput());
	connector->Update();


	
	CommonTools::SaveImage( connector->GetOutput(), outp );
*/
	return 0;
}

