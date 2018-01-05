#include "QSTEM.h"
#include "data_containers.h"
#include "stemlib.h"
#include "stemutil.h"
#include <string.h>
#include <assert.h>
#include <memory>
#include <cmath>

#ifndef _WIN32
#ifdef __cplusplus
#include <cmath>
#else
#include <math.h>
#endif
#else
#include <math.h>
#endif

#define RAD2DEG 57.2958
int fftMeasureFlag = FFTW_ESTIMATE;
using namespace shapes;

QSTEM::QSTEM(std::string mode)
{

  if (mode=="STEM") muls.mode = STEM;
  else if (mode=="TEM") muls.mode = TEM;
  else if (mode=="CBED") muls.mode = CBED;
  else if (mode=="NBED") muls.mode = NBED;
  else if (mode=="TOMO") muls.mode = TOMO;
  else if (mode=="REFINE") muls.mode = REFINE;

  //std::vector<DetectorPtr> detectors;
  //muls.detectors.push_back(detectors);

  muls.printLevel = 0;
  muls.saveLevel = 0;

  muls.nx = 0; // resolution x
  muls.ny = 0; // resolution y
  muls.slices = 0;

  muls.complete_pixels=0;

  muls.atoms=NULL;
  muls.trans = NULL;
  muls.diffpat = NULL;

  muls.ax = 0; // box side length (set by set_box)
  muls.by = 0; // box side length (set by set_box)
  muls.c = 0;  // box height (set by set_box)
  muls.cellDiv = 1; // cell z division (set by set_box)
  muls.equalDivs = 1; // equal cell divisions (fixed)
  muls.nonPeriodZ = 1; // periodic in z (set by set_box)
  muls.nonPeriod = 0; // periodic in x and y (set by set_box)

  muls.tds = 0; // turn on tds
  muls.avgCount = 0; // current run # for tds
  muls.avgRuns = 1; // max run # for tds
  muls.Einstein = 1; // use Einstein model (fixed)
  muls.atomRadius = 5.0; // atomic radius (fixed)

  muls.fileBase[0] = '\0'; // file base (unused)
  muls.fileWaveIn[0] = '\0'; // input wave file (unused)
  muls.atomPosFile[0] = '\0'; // cfg file (unused)
  muls.folder[0] = '\0'; // folder (unused)
  muls.phononFile[0] = '\0'; // phonons file (unused)
  muls.cfgFile[0] = '\0'; // cfg file (unused)

  muls.nCellX = 1; // cell repetition (fixed)
  muls.nCellY = 1; // cell repetition (fixed)
  muls.nCellZ = 1; // cell repetition (fixed)
  muls.ctiltx = 0; // crystal tilt (fixed)
  muls.ctilty = 0;	// crystal tilt (fixed)
  muls.ctiltz = 0; // crystal tilt (fixed)
  muls.adjustCubeSize = 0; // adjust box size (fixed)
  muls.cubex = 0; // adjusted box side length (fixed)
  muls.cubey = 0; // adjusted box side length (fixed)
  muls.cubez = 0; // adjusted box side length (fixed)
  muls.xOffset = 0; // atom offset (fixed)
  muls.yOffset = 0; // atom offset (fixed)
  muls.czOffset = 0; // atom offset (fixed)
  muls.cz = NULL;
  muls.centerSlices = 0;

  muls.bandlimittrans = 0; //  (fixed)
  muls.readPotential = 0; // read potential from file (fixed)
  muls.savePotential = 0; // save potential to file (fixed)
  muls.saveTotalPotential = 0; // save projected potential to file (fixed)
  muls.plotPotential = 0; // plot potential (fixed)
  muls.fftpotential = 1; // fftpotential (fixed)
  muls.potential3D = 1; // potential 3d (fixed)

  muls.totalSliceCount = 0;
  muls.storeSeries = 0; // ??? (fixed)
  muls.scatFactor = WEICK_KOHL;

  muls.webUpdate = 0; // ??? (never used)
  muls.imageGamma = 1.0;
  muls.showProbe = 0;

  muls.ismoth = 1;
  muls.u2 = NULL;
  muls.u2avg = NULL;
  muls.detectorNum = 0;

  muls.beamCurrent = 1;  // pico Ampere
  muls.dwellTime = 1;    // msec
  muls.electronScale = muls.beamCurrent*muls.dwellTime*MILLISEC_PICOAMP;


}

QSTEM::~QSTEM()
{
	
	if (muls.atoms!=NULL){
		free(muls.atoms);
	}
	
	if (muls.trans != NULL){
		fftwf_free(muls.trans[0][0]);
		for(int i = 0; i < muls.slices; i++)
			fftwf_free(muls.trans[i]);
		fftwf_free(muls.trans);
	}
		
}

void QSTEM::set_atoms(int natom, int atomKinds, const std::vector< std::vector<double> > & pos,
                        const std::vector<double> & occ, const std::vector<double> & q,
                        const std::vector<int> & Znum)
{
  if (muls.atoms!=NULL){
      free(muls.atoms);
  }
  
  muls.atoms = (atom *)malloc(natom*sizeof(atom));

  for (int i=0; i<natom; i++){
    muls.atoms[i].x = pos[i][0];
    muls.atoms[i].y = pos[i][1];
    muls.atoms[i].z = pos[i][2];
    muls.atoms[i].dw = 0.; //dw[i];
    muls.atoms[i].occ = occ[i];
    muls.atoms[i].q = q[i];
    muls.atoms[i].Znum = Znum[i];
  }
  muls.natom = natom;
  muls.atomKinds = atomKinds;
}

void QSTEM::set_positions(int natom,const std::vector< std::vector<double> > & pos)
{
  for (int i=0; i<natom; i++){
    muls.atoms[i].x = pos[i][0];
    muls.atoms[i].y = pos[i][1];
    muls.atoms[i].z = pos[i][2];
  }
}

void QSTEM::set_box(const std::vector<double> & box, int nonPeriod, int nonPeriodZ, float cellDiv)
{
  muls.ax = box[0];
  muls.by = box[1];
  muls.c = box[2];
  muls.nonPeriod = nonPeriod;
  muls.nonPeriodZ = nonPeriodZ;
  muls.cellDiv = cellDiv;
}

void QSTEM::get_resolution(float* resolutionX, float* resolutionY)
{
  (*resolutionX)=muls.resolutionX;
  (*resolutionY)=muls.resolutionY;
}

void QSTEM::get_potential_samples(int* potNx, int* potNy, int* slices)
{
  (*potNx)=muls.potNx;
  (*potNy)=muls.potNy;
  (*slices)=muls.slices;
}

void QSTEM::get_potential_extent(float* potOffsetX,float* potOffsetY,float* potSizeX,float* potSizeY)
{
  (*potOffsetX)=muls.potOffsetX;
  (*potOffsetY)=muls.potOffsetY;
  (*potSizeX)=muls.potSizeX;
  (*potSizeY)=muls.potSizeY;
}

void QSTEM::get_probe_samples(int* nx, int* ny)
{
  (*nx)=muls.nx;
  (*ny)=muls.ny;
}

void QSTEM::get_scan_range(float* scanXStart, float* scanXStop, int* scanXN, float* scanYStart, float* scanYStop, int* scanYN)
{
  (*scanXStart)=muls.scanXStart;
  (*scanXStop)=muls.scanXStop;
  (*scanXN)=muls.scanXN;
  (*scanYStart)=muls.scanYStart;
  (*scanYStop)=muls.scanYStop;
  (*scanYN)=muls.scanYN;
}

void QSTEM::get_energy(float* v0)
{
  (*v0)=muls.v0;
}

void QSTEM::allocate_potential()
{
  if (muls.trans != NULL){
    fftwf_free(muls.trans[0][0]);
    fftwf_free(muls.trans[0]);
    fftwf_free(muls.trans);
    fftwf_destroy_plan(muls.fftPlanPotForw);
    fftwf_destroy_plan(muls.fftPlanPotInv);
    muls.trans = NULL;
  }

  muls.sliceThickness = muls.c / (double)muls.slices;

  int potDimensions[2];
 	potDimensions[0] = muls.potNx;
  potDimensions[1] = muls.potNy;
	muls.trans = complex3Df(muls.slices,muls.potNx,muls.potNy,"trans");

  muls.fftPlanPotForw = fftwf_plan_many_dft(2,potDimensions, muls.slices,muls.trans[0][0], NULL,
	                                          1, muls.potNx*muls.potNy,muls.trans[0][0], NULL,
                                            1, muls.potNx*muls.potNy, FFTW_FORWARD, fftMeasureFlag);
  muls.fftPlanPotInv = fftwf_plan_many_dft(2,potDimensions, muls.slices,muls.trans[0][0], NULL,
	                                          1, muls.potNx*muls.potNy,muls.trans[0][0], NULL,
                                            1, muls.potNx*muls.potNy, FFTW_BACKWARD, fftMeasureFlag);
}

void QSTEM::build_potential(int slices)
{

  if ((muls.mode == STEM) || (muls.mode == CBED)) {
		muls.potOffsetX = muls.scanXStart - 0.5*muls.nx*muls.resolutionX;
		muls.potOffsetY = muls.scanYStart - 0.5*muls.ny*muls.resolutionY;
		muls.potNx = (int)((muls.scanXStop-muls.scanXStart)/muls.resolutionX);
		muls.potNy = (int)((muls.scanYStop-muls.scanYStart)/muls.resolutionY);
		muls.scanXStop = muls.scanXStart+muls.resolutionX*muls.potNx;
		muls.scanYStop = muls.scanYStart+muls.resolutionY*muls.potNy;
		muls.potNx += muls.nx;
		muls.potNy += muls.ny;
		muls.potSizeX = muls.potNx*muls.resolutionX;
		muls.potSizeY = muls.potNy*muls.resolutionY;
	}
	else {
    muls.scanXStart = muls.ax/2.0;
    muls.scanYStart = muls.by/2.0;
    muls.scanXN = 1;
    muls.scanYN = 1;
    muls.scanXStop = muls.scanXStart;
    muls.scanYStop = muls.scanYStart;
		muls.potNx = muls.nx;
		muls.potNy = muls.ny;
		muls.potSizeX = muls.potNx*muls.resolutionX;
		muls.potSizeY = muls.potNy*muls.resolutionY;
		muls.potOffsetX = muls.scanXStart - 0.5*muls.potSizeX;
		muls.potOffsetY = muls.scanYStart - 0.5*muls.potSizeY;
	}

  muls.slices=slices;
 	muls.nlayer=slices;
  muls.outputInterval = muls.slices;

  allocate_potential();


  make3DSlices(&muls,muls.slices,muls.atomPosFile,NULL);

}

void QSTEM::set_potential(const std::vector< std::vector< std::vector< std::vector<double> > > > & potential,
                  const std::vector<int> & size, const std::vector<double> & extent)
{

  muls.potNx=size[0];
  muls.potNy=size[1];
  muls.slices=size[2];
  muls.nlayer=muls.slices;
  muls.outputInterval = muls.slices;
  muls.resolutionX=(extent[1]-extent[0])/((double)size[0]);
  muls.resolutionY=(extent[3]-extent[2])/((double)size[1]);
  muls.sliceThickness=(extent[5]-extent[4])/((double)size[2]);

  if (muls.cz == NULL) {
    muls.cz = float1D(muls.nlayer,"cz");
  }
  muls.cz[0]=muls.sliceThickness;

  muls.nlayer=muls.slices;

  muls.potOffsetX = 0;
  muls.potOffsetY = 0;
  //muls.scanXStop = muls.scanXStart+muls.resolutionX*muls.potNx;
  //muls.scanYStop = muls.scanYStart+muls.resolutionY*muls.potNy;

  muls.potSizeX = muls.potNx*muls.resolutionX;
  muls.potSizeY = muls.potNy*muls.resolutionY;

  if (muls.mode == TEM){
    muls.scanXStart = muls.ax/2.0;
    muls.scanYStart = muls.by/2.0;
    muls.scanXN = 1;
    muls.scanYN = 1;
    muls.scanXStop = muls.scanXStart;
    muls.scanYStop = muls.scanYStart;
  }

  allocate_potential();


  for (int iz=0;iz<muls.slices;iz++){
    for (int ix=0;ix<muls.potNx;ix++){
      for (int iy=0;iy<muls.potNy;iy++){
        for (int ic=0;ic<2;ic++){
          muls.trans[iz][ix][iy][ic] = potential.at(ix).at(iy).at(iz).at(ic);
        }
      }
    }
  }
}

std::vector <std::vector <std::vector <std::vector <double> > > > QSTEM::get_potential_or_transfunc(float* resolutionX,float* resolutionY,float* sliceThickness)
{

    using namespace std;
    vector<vector<vector<vector<double>>>> trans_vec(muls.potNx,vector<vector<vector<double>>>(muls.potNy,vector<vector<double>>(muls.slices,vector<double>(2))));

    for (int iz=0;iz<muls.slices;iz++){
      for (int ix=0;ix<muls.potNx;ix++){
        for (int iy=0;iy<muls.potNy;iy++){
          for (int ic=0;ic<2;ic++){
            trans_vec.at(ix).at(iy).at(iz).at(ic) = muls.trans[iz][ix][iy][ic];
          }
        }
      }
    }

    (*resolutionX)=muls.resolutionX;
    (*resolutionY)=muls.resolutionY;
    (*sliceThickness)=muls.sliceThickness;

    return trans_vec;
}

void QSTEM::calculate_transfunc()
{
  muls.mulsRepeat1 = 1;
	muls.mulsRepeat2 = 1;

  initSTEMSlices(&muls,muls.slices);
}

void QSTEM::set_scan_range(float scanXStart,float scanXStop,int scanXN,float scanYStart,float scanYStop,int scanYN)
{
  muls.scanXStart = scanXStart;
  muls.scanXStop = scanXStop;
  muls.scanXN = scanXN;
  muls.scanYStart = scanYStart;
  muls.scanYStop = scanYStop;
  muls.scanYN = scanYN;
}

void QSTEM::free_wave(){
  if (wave != NULL){
    fftwf_free(wave->wave[0]);
    fftwf_free(wave->wave);
    fftwf_free(wave->diffpat[0]);
    fftwf_free(wave->diffpat);
    fftwf_free(wave->avgArray[0]);
    fftwf_free(wave->avgArray);
    fftwf_destroy_plan(wave->fftPlanWaveForw);
    fftwf_destroy_plan(wave->fftPlanWaveInv);
  }
}

void QSTEM::set_wave(std::vector <std::vector <std::vector <double> > >wave_vec,float v0,int nx,int ny,float resolutionX,float resolutionY)
{
  free_wave();

  muls.nx = nx;
  muls.ny = ny;

  if (resolutionX <= 0){
    muls.resolutionX = muls.ax / (double)muls.nx;
  }
  else {
    muls.resolutionX = resolutionX;
  }
  if (resolutionX < 0){
    muls.resolutionY = muls.by / (double)muls.ny;
  }
  else {
    muls.resolutionY = resolutionY;
  }

  muls.v0 = v0;
  wave = WavePtr(new WAVEFUNC(muls.nx,muls.ny,muls.resolutionX,muls.resolutionY));
    for (int ix=0;ix<muls.nx;ix++){
      for (int iy=0;iy<muls.ny;iy++){
        for (int ic=0;ic<2;ic++){
          wave->wave[ix][iy][ic] = wave_vec.at(ix).at(iy).at(ic);
      }
    }
  }
}

void QSTEM::build_wave(int type,float v0,int nx,int ny,float resolutionX,float resolutionY)
{
  free_wave();
  
  muls.nx = nx;
  muls.ny = ny;
  if (resolutionX <= 0){
    muls.resolutionX = muls.ax / (double)muls.nx;
  }
  else {
    muls.resolutionX = resolutionX;
  }
  if (resolutionX < 0){
    muls.resolutionY = muls.by / (double)muls.ny;
  }
  else {
    muls.resolutionY = resolutionY;
  }

  muls.v0 = v0;

  wave = WavePtr(new WAVEFUNC(muls.nx,muls.ny,muls.resolutionX,muls.resolutionY));
  //wave = new WAVEFUNC(muls.nx,muls.ny,muls.resolutionX,muls.resolutionY);

  if (type==0){
    for (int ix=0;ix<muls.nx;ix++){
      for (int iy=0;iy<muls.ny;iy++){
        wave->wave[ix][iy][0]=1.;
        wave->wave[ix][iy][1]=0;
      }
    }
  }
}

void QSTEM::build_probe(float v0,float alpha,int nx,int ny,float resolutionX,float resolutionY,std::unordered_map<std::string,float> abberations)
{
  free_wave();

  muls.nx = nx;
  muls.ny = ny;
  muls.resolutionX = resolutionX;
  muls.resolutionY = resolutionY;

  muls.df0 = abberations["a20"];
  muls.Cs = abberations["a40"];
  muls.C5 = abberations["a60"];
  muls.astigMag = abberations["a22"];

  muls.a33 = abberations["a33"];
  muls.a31 = abberations["a31"];
  muls.a44 = abberations["a44"];
  muls.a42 = abberations["a42"];
  muls.a55 = abberations["a55"];
  muls.a53 = abberations["a53"];
  muls.a51 = abberations["a51"];
  muls.a66 = abberations["a66"];
  muls.a64 = abberations["a64"];
  muls.a62 = abberations["a62"];
  muls.astigAngle = abberations["phi22"] / (float)RAD2DEG;
  muls.phi33 = abberations["phi33"] / (float)RAD2DEG;
  muls.phi31 = abberations["phi31"] / (float)RAD2DEG;
  muls.phi44 = abberations["phi44"] / (float)RAD2DEG;
  muls.phi42 = abberations["phi42"] / (float)RAD2DEG;
  muls.phi55 = abberations["phi55"] / (float)RAD2DEG;
  muls.phi53 = abberations["phi53"] / (float)RAD2DEG;
  muls.phi51 = abberations["phi51"] / (float)RAD2DEG;
  muls.phi66 = abberations["phi66"] / (float)RAD2DEG;
  muls.phi64 = abberations["phi64"] / (float)RAD2DEG;
  muls.phi62 = abberations["phi62"] / (float)RAD2DEG;

  muls.v0 = v0;
  muls.alpha = alpha;

  muls.gaussScale = 0.05f;
	muls.gaussFlag = 0;
  muls.aAIS = 0; // AIS aperture
  muls.dE_E = 0.0;
	muls.dI_I = 0.0;
	muls.dV_V = 0.0;
	muls.Cc = 0.0;
  float dE_E0;
  dE_E0 = sqrt(muls.dE_E*muls.dE_E+muls.dI_I*muls.dI_I+muls.dV_V*muls.dV_V);
  muls.dE_EArray = (double *)malloc((muls.avgRuns+1)*sizeof(double));
  muls.dE_EArray[0] = 0.0;

  wave = WavePtr(new WAVEFUNC(muls.nx,muls.ny,muls.resolutionX,muls.resolutionY));
  probe(&muls, wave, muls.nx/2*muls.resolutionX, muls.ny/2*muls.resolutionY);

}

std::vector <std::vector <std::vector <double> > > QSTEM::get_wave(float* resolutionX,float* resolutionY,float* v0)
{

  using namespace std;
  vector<vector<vector<double>>> wave_vec(muls.nx,vector<vector<double>>(muls.ny,vector<double>(2)));

    for (int ix=0;ix<muls.nx;ix++){
      for (int iy=0;iy<muls.ny;iy++){
        for (int ic=0;ic<2;ic++){
          wave_vec.at(ix).at(iy).at(ic) = wave->wave[ix][iy][ic];
        }
      }
    }

    (*resolutionX)=muls.resolutionX;
    (*resolutionY)=muls.resolutionY;
    (*v0)=muls.v0;

    return wave_vec;
}

void QSTEM::create_detectors(std::vector<std::string> names,std::vector<std::vector<double>> properties, int detectorNum)
{

  muls.detectors.clear();
  int tCount = (int)(ceil((double)((muls.slices * muls.cellDiv) / muls.outputInterval)));

  for (int islice=0; islice<=tCount; islice++)
  {

    std::vector<DetectorPtr> detectors;
    for (int i=0; i<detectorNum; i++){

      DetectorPtr det = DetectorPtr(new Detector(muls.scanXN, muls.scanYN,
        (muls.scanXStop-muls.scanXStart)/(float)muls.scanXN,
        (muls.scanYStop-muls.scanYStart)/(float)muls.scanYN));

      det->rInside = properties[i][0];
      det->rOutside = properties[i][1];
      det->name = names[i].c_str();
      det->shiftX = properties[i][2];
      det->shiftY = properties[i][3];

      det->k2Inside = (float)(sin(det->rInside*0.001)/(wavelength(muls.v0)));
      det->k2Outside = (float)(sin(det->rOutside*0.001)/(wavelength(muls.v0)));

      det->k2Inside *= det->k2Inside;
      det->k2Outside *= det->k2Outside;
      detectors.push_back(det);

    }
    muls.detectors.push_back(detectors);

  }
  muls.detectorNum = detectorNum;

}

void QSTEM::run(int displayProgInterval)
{
  muls.displayProgInterval = displayProgInterval;

  double timer, total_time=0;

  //timer=cputim();
  if (muls.mode == STEM){
    muls.complete_pixels = 0;

    int ix, iy;
    for (int i=0; i < (muls.scanXN * muls.scanYN); i++){
      ix = i / muls.scanYN;
      iy = i % muls.scanYN;

      WavePtr wave_tmp;
      wave_tmp = WavePtr(new WAVEFUNC(muls.nx,muls.ny,muls.resolutionX,muls.resolutionY));
      for (int wix=0;wix<muls.nx;wix++){
        for (int wiy=0;wiy<muls.ny;wiy++){
          for (int wic=0;wic<2;wic++){
            wave_tmp->wave[wix][wiy][wic] = wave->wave[wix][wiy][wic];
          }
        }
      }

      //wave->iPosX =(int)(ix*(muls.scanXStop-muls.scanXStart)/
      //          ((float)muls.scanXN*muls.resolutionX));
      //wave->iPosY = (int)((iy*(muls.scanYStop-muls.scanYStart))/
      //           ((float)muls.scanYN*muls.resolutionY));
      wave_tmp->iPosX =(int)(muls.scanXStart/muls.resolutionX+(ix*(muls.scanXStop-muls.scanXStart)/
                ((float)muls.scanXN*muls.resolutionX))-((float)muls.nx/2.)-muls.potOffsetX/muls.resolutionX);
      wave_tmp->iPosY =(int)(muls.scanYStart/muls.resolutionY+(iy*(muls.scanYStop-muls.scanYStart)/
                          ((float)muls.scanYN*muls.resolutionY))-((float)muls.ny/2.)-muls.potOffsetY/muls.resolutionY);

      //if (wave->iPosX > muls.potNx-muls.nx){
      //  wave->iPosX = muls.potNx-muls.nx;
      //}

      //if (wave->iPosY > muls.potNy-muls.ny){
      //  wave->iPosY = muls.potNy-muls.ny;
      //}

      wave_tmp->detPosX=ix;
      wave_tmp->detPosY=iy;

      runMulsSTEM(&muls,wave_tmp);

      muls.complete_pixels++;
      /*
      if (muls.displayProgInterval > 0) if ((muls.complete_pixels) % muls.displayProgInterval == 0)
      {
        total_time += cputim()-timer;
        printf("Pixels complete: (%d/%d), int.=%.3f, avg time per pixel: %.2fsec\n",
          muls.complete_pixels, muls.scanXN*muls.scanYN, wave->intIntensity,
          (total_time)/muls.complete_pixels);
        timer=cputim();
      }
      */
    }
  }
  else if (muls.mode == CBED){
    wave->iPosX =(int)(muls.scanXStart/muls.resolutionX-((float)muls.nx/2.)-muls.potOffsetX/muls.resolutionX);
    wave->iPosY =(int)(muls.scanYStart/muls.resolutionY-((float)muls.ny/2.)-muls.potOffsetY/muls.resolutionY);
    runMulsSTEM(&muls,wave);
  }
  else{
    runMulsSTEM(&muls,wave);
  }
}

std::vector <std::vector <double> > QSTEM::read_detector(int detNum)
{

  std::vector<DetectorPtr> detectors;
  int tCount = (int)(ceil((double)((muls.slices * muls.cellDiv) / muls.outputInterval)));

  detectors = muls.detectors[tCount];
  std::vector<std::vector<double>> img_vec(muls.scanXN,std::vector<double>(muls.scanYN));

  for (int ix=0;ix<muls.scanXN;ix++){
    for (int iy=0;iy<muls.scanYN;iy++){
      img_vec.at(ix).at(iy) = detectors[detNum]->image[ix][iy];
    }
  }

  return img_vec;
}
