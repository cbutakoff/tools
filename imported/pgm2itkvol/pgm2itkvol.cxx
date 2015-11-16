/*
 author  Rubén Cárdenes, Nov. 2012 
 adapted from Michel Couprie's pink software
/* Michel Couprie - janvier 2000 */

#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <stdlib.h>
#include <mcimage.h>
#include <mccodimage.h>

#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkImageFileWriter.h>

#define VERBOSE

#define UINTDATA unsigned int *
/* =============================================================== */
int main(int argc, char **argv)
/* =============================================================== */
{
  int32_t rs, cs, ds, N, ret;
  struct xvimage * image;

  if (argc != 3)
  {
    fprintf(stderr, "usage: %s in.pgm out.raw \n", argv[0]);
    exit(1);
  }

  char* output_file = argv[2];
	
  image = readimage(argv[1]);
  if (image == NULL)
  {
    fprintf(stderr, "%s: readimage failed\n", argv[0]);
    exit(1);
  }

  if ( (datatype(image) != VFF_TYP_1_BYTE) && (datatype(image) != VFF_TYP_4_BYTE) && (datatype(image) != VFF_TYP_FLOAT))
  {
    fprintf(stderr, "%s: only byte|long|float images supported\n", argv[0]);
    exit(1);
  }

  rs = rowsize(image);
  cs = colsize(image);
  ds = depth(image);
  N = rs * cs * ds;
	
#ifdef VERBOSE
  printf("rs = %d ; cs = %d ; ds = %d ; N = rs * cs * ds = %d\n", rs, cs, ds, N);
#endif

		
	
	
  if(datatype(image) == VFF_TYP_1_BYTE)
  {
	  printf("unsigned char\n");

	  typedef unsigned char               PixelType;
	  typedef itk::Image< PixelType,  3 >  ImageType3D;
	  typedef itk::ImageFileWriter< ImageType3D >  WriterType3D;
	  typedef itk::ImageRegionIterator<ImageType3D> IteratorType;
	  
	  ImageType3D::Pointer input_image = ImageType3D::New();
	  ImageType3D::RegionType region;
	  ImageType3D::SizeType image_size;
	  image_size[0] = rs;
	  image_size[1] = cs;
	  image_size[2] = ds;
	  region.SetSize(image_size);
	  input_image->SetRegions(region);
	  input_image->Allocate();
	  
	  ImageType3D::SpacingType spacing;
	  spacing[0] = image->xdim;
	  spacing[1] = image->ydim;
	  spacing[2] = image->zdim;
	  
	  ImageType3D::PointType origin;
	  origin[0] = image->origin_x;
	  origin[1] = image->origin_y;
	  origin[2] = image->origin_z;
	  
	  input_image->SetSpacing(spacing);
	  input_image->SetOrigin(origin);
	  
	  IteratorType It( input_image, input_image->GetLargestPossibleRegion() );
	  int i=0;
	  unsigned char* I = UCHARDATA(image);
	  for( It.GoToBegin(); !It.IsAtEnd(); ++It  ) {
		  It.Set(I[i]) ;
		  i++;
	  }	  
	  
	  WriterType3D::Pointer writer = WriterType3D::New();
	  writer->SetInput(input_image);
	  writer->SetFileName(output_file);
	  writer->Update();
	  
  }
  else if(datatype(image) == VFF_TYP_2_BYTE)
  {
  	  printf("unsigned short\n");

  }
  else if(datatype(image) == VFF_TYP_4_BYTE)
  {
  	  printf("int\n");

  }
  else if(datatype(image) == VFF_TYP_FLOAT)
  {
	  printf("float\n");

	  typedef float                        PixelType;
	  typedef itk::Image< PixelType,  3 >  ImageType3D;
	  typedef itk::ImageFileWriter< ImageType3D >  WriterType3D;
	  typedef itk::ImageRegionIterator<ImageType3D> IteratorType;
	  
	  ImageType3D::Pointer input_image = ImageType3D::New();
	  ImageType3D::RegionType region;
	  ImageType3D::SizeType image_size;
	  image_size[0] = rs;
	  image_size[1] = cs;
	  image_size[2] = ds;
	  region.SetSize(image_size);
	  input_image->SetRegions(region);
	  input_image->Allocate();
	  
	  ImageType3D::SpacingType spacing;
	  spacing[0] = image->xdim;
	  spacing[1] = image->ydim;
	  spacing[2] = image->zdim;
	  
	  ImageType3D::PointType origin;
	  origin[0] = image->origin_x;
	  origin[1] = image->origin_y;
	  origin[2] = image->origin_z;
	  
	  input_image->SetSpacing(spacing);
	  input_image->SetOrigin(origin);
	  
	  IteratorType It( input_image, input_image->GetLargestPossibleRegion() );
	  int i=0;
	  float* I = FLOATDATA(image);
	  for( It.GoToBegin(); !It.IsAtEnd(); ++It  ) {
		  It.Set(I[i]) ;
		  i++;
	  }	  
	  
	  WriterType3D::Pointer writer = WriterType3D::New();
	  writer->SetInput(input_image);
	  writer->SetFileName(output_file);
	  writer->Update();
  
  }
  

  return 0;
} /* main */


