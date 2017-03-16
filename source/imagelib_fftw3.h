#ifndef IMAGELIB_H
#define IMAGELIB_H

#include "stemtypes_fftw3.h"
#include <vector>
#include <string>
#include <iostream>
#include <memory>
//#include "boost/shared_ptr.hpp"

/**************************************************************
 * Here is how to use the new image writing routines
 *
 * ImageIOPtr imageio = ImageIOPtr(new CImageIO(nx,ny))
 *
 * imageIO->WriteRealImage((void**)real_image,filename);
 * imageIO->WriteComplexImage((void**)complex_image,filename);
 *
 **************************************************************
 * Reading an image works like this:
 * 
 * ImageIOPtr imageio = ImageIOPtr(new CImageIO(nx,ny))
 * imageio->ReadImage((void **)pix,nx,ny,fileName);
 *
 * Note that header parameters are persistent on any object.  You
 *   should use the various Set* functions to set parameters as
 *   necessary.  You should not need to read values from this class - 
 *   only set them.  They will be recorded to any file saved from this
 *   this object.
 **************************************************************/

class CImageIO {
  int m_headerSize;  // first byte of image will be size of image header (in bytes)
                   // This is the size without the data, parameters, and comment!!!
  int m_paramSize;   // number of additional parameters
  int m_commentSize; // length of comment string
  int m_nx,m_ny;
  int m_complexFlag;
  int m_dataSize;    // size of one data element in bytes (e.g. complex double: 16)
  int m_version;     // The version flag will later help to find out how to 
                   // distinguish between images produced by different versions of stem
  double m_t;        // thickness
  double m_dx,m_dy;    // size of one pixel
  std::vector<double> m_params;  // array for additional parameters
  std::string m_comment;   // comment of prev. specified length
  char m_buf[200];  // General purpose temporary text buffer
public:
  CImageIO(int nx, int ny);
  CImageIO(int nx, int ny, double t, double dx, double dy,
           std::vector<double> params=std::vector<double>(), 
		   std::string comment="");

  void WriteRealImage(void **pix, const char *fileName);
  void WriteComplexImage(void **pix, const char *fileName);
  void ReadImage(void **pix, int nx, int ny, const char *fileName);
  
  //void WriteImage( std::string fileName);
        
  void SetComment(std::string comment);
  void SetThickness(double thickness);
  void SetParams(std::vector<double> params);
  void SetParameter(int index, double value);
  void SetResolution(double resX, double resY);
private:
  void WriteData(void **pix, const char *fileName);
  // reads in the header; returns the byte offset at which we should start reading image data.
  void ReadHeader(const char *fileName);
};

//typedef boost::shared_ptr<CImageIO> ImageIOPtr;
typedef std::shared_ptr<CImageIO> ImageIOPtr;

#endif
