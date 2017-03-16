/*
QSTEM - image simulation for TEM/STEM/CBED
    Copyright (C) 2000-2010  Christoph Koch
	Copyright (C) 2010-2013  Christoph Koch, Michael Sarahan

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "stdio.h"
#include <string.h>
#include "data_containers.h"

WAVEFUNC::WAVEFUNC(int x, int y, float_tt resX, float_tt resY) :
detPosX(0),
detPosY(0),
iPosX(0),
iPosY(0),
thickness(0.0),
nx(x),
ny(y),
resolutionX(resX),
resolutionY(resY)
{
	char waveFile[256];
	const char *waveFileBase = "mulswav";
#if FLOAT_PRECISION == 1
	diffpat = float2D(nx,ny,"diffpat");
	avgArray = float2D(nx,ny,"avgArray");
#else
	diffpat = double2D(nx,ny,"diffpat");
	avgArray = double2D(nx,ny,"avgArray");
#endif

	m_imageIO=ImageIOPtr(new CImageIO(nx, ny, thickness, resolutionX, resolutionY));


#if FLOAT_PRECISION == 1
	wave = complex2Df(nx, ny, "wave");
	fftPlanWaveForw = fftwf_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_FORWARD, FFTW_ESTIMATE);
	fftPlanWaveInv = fftwf_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_BACKWARD, FFTW_ESTIMATE);
#else
	wave = complex2D(nx, ny, "wave");
	fftPlanWaveForw = fftw_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_FORWARD,
		fftMeasureFlag);
	fftPlanWaveInv = fftw_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_BACKWARD,
		fftMeasureFlag);
#endif

	sprintf(waveFile,"%s.img",waveFileBase);
	strcpy(fileout,waveFile);
	sprintf(fileStart,"mulswav.img");
}

void WAVEFUNC::WriteWave(const char *fileName, const char *comment,
	std::vector<double>params)
{
	m_imageIO->SetComment(comment);
	m_imageIO->SetResolution(resolutionX, resolutionY);
	m_imageIO->SetParams(params);
	m_imageIO->SetThickness(thickness);
	m_imageIO->WriteComplexImage((void **)wave, fileName);
}

void WAVEFUNC::WriteDiffPat(const char *fileName, const char *comment,
	std::vector<double>params)
{
	m_imageIO->SetComment(comment);
	m_imageIO->SetResolution(1.0/(nx*resolutionX), 1.0/(ny*resolutionY));
	m_imageIO->SetParams(params);
	m_imageIO->SetThickness(thickness);
	m_imageIO->WriteRealImage((void**)diffpat, fileName);
}

void WAVEFUNC::WriteAvgArray(const char *fileName, const char *comment,
	std::vector<double>params)
{
	m_imageIO->SetComment(comment);
	m_imageIO->SetResolution(1.0/(nx*resolutionX), 1.0/(ny*resolutionY));
	m_imageIO->SetParams(params);
	m_imageIO->SetThickness(thickness);
	m_imageIO->WriteRealImage((void **)avgArray, fileName);
}

void WAVEFUNC::ReadWave(const char *fileName)
{
	// printf("Debug Wavefunc::ReadWave\n");
	m_imageIO->ReadImage((void **)wave, nx, ny, fileName);
}

void WAVEFUNC::ReadDiffPat(const char *fileName)
{
	m_imageIO->ReadImage((void **)diffpat, nx, ny, fileName);
}

void WAVEFUNC::ReadAvgArray(const char *fileName)
{
	m_imageIO->ReadImage((void **)avgArray, nx, ny, fileName);
}



Detector::Detector(int nx, int ny, float_tt resX, float_tt resY) :
  error(0),
  shiftX(0),
  shiftY(0),
  Navg(0),
  thickness(0)
{
#if FLOAT_PRECISION == 1
	image = float2D(nx,ny,"ADFimag");
	image2 = float2D(nx,ny,"ADFimag");
#else
	image = double2D(nx,ny,"ADFimag");
	image2 = double2D(nx,ny,"ADFimag");
#endif
	m_imageIO=ImageIOPtr(new CImageIO(nx, ny, thickness, resX, resY, std::vector<double>(2+nx*ny), "STEM image"));
}

void Detector::WriteImage(const char *fileName)
{
	m_imageIO->SetThickness(thickness);
	m_imageIO->WriteRealImage((void **)image, fileName);
}

void Detector::SetThickness(float_tt t)
{
	thickness=t;
}

void Detector::SetParameter(int index, double value)
{
	m_imageIO->SetParameter(index, value);
}

void Detector::SetParams(std::vector<double> params)
{
	m_imageIO->SetParams(params);
}

void Detector::SetComment(const char *comment)
{
	m_imageIO->SetComment(comment);
}
