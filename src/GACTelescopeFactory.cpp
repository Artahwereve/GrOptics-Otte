/*
VERSION4.0
30May2016
*/
/*!  GSegSCTelescopeFactory.cpp

     Charlie Duke
     Grinnell College
     May 2011

 */
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <list>
#include <iterator>
#include <algorithm>
#include <bitset>
#include <iomanip>
#include <assert.h>

#include "TROOT.h"
#include "TGraph.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "Math/Vector3D.h"

using namespace std;

#include "GDefinition.h"

#include "GPilot.h"
#include "GUtilityFuncts.h"

#include "AOpticsManager.h" 
#include "GTelescope.h"
#include "GSegSCTelescope.h"
#include "GACTelescope.h"
#include "GSCTelescope.h"
#include "GTelescopeFactory.h"
#include "GSCTelescopeFactory.h"
#include "GSegSCTelescopeFactory.h"
#include "GACTelescopeFactory.h"

#include "GReadSegSCStd.h"
#include "GReadACStd.h"

#define DEBUG(x) *oLog << #x << " = " << x << endl
#define DEBUGS(x) *oLog << "       "<< #x << " = " << x << endl

// define useful units
static const double cm = AOpticsManager::cm();
static const double mm = AOpticsManager::mm();
static const double um = AOpticsManager::um();
static const double nm = AOpticsManager::nm();
static const double  m = AOpticsManager::m();

/*!  /brief SegSCStdOptics structure stores details of a standard 
     Davis-Cotton telescope
 */
ACStdOptics::ACStdOptics() {

  bool debug = false;
  if (debug) {
    *oLog << "  -- ACStdOptics::ACStdOptics " << endl;
  }

  
    //geoStruct = 0;

    //igeoStd = 0;
    stdType         = AC;
    defoc           = 0.0;
    foclength       = 0.0;
    camRadius       = 0.0;
    plateScaleFactor = 0.0;
    radius          = 0.0;
    rotation_offset = 0.0;
    xoff            = 0.0;
    yoff            = 0.0;
    nbr_mir         = 0;
    nbinsx          = 0;
    nbinsy          = 0;
    gridOption      = 0;
    gridFile        = "";
    grid            = 0;
    //rayTracerEnum   = RTDC;
    fAvgTransitTime = 0.0;
    
    cameraRadius = 0.0;

    oStr = oLog;
    iPrtMode = 0;
};
/************** end of SegSCStdOptics ***********************/

ACStdOptics::ACStdOptics(const ACStdOptics &sco) {

  (void) sco; // unused
  bool debug = false;
  if (debug) {
    *oLog << "  -- ACStdOptics::ACStdOptics " << endl;
  }
  

};

/************** end of SegSCStdOptics ***********************/
ACStdOptics::~ACStdOptics() {
  bool debug = false;
  if (debug) {
    *oLog << "  -- ACStdOptics::~ACStdOptics" << endl;
  }
    SafeDelete(grid);
      for (unsigned i = 0;i<facet.size();i++) {
        SafeDelete(facet[i]);
      }
};
/********************* end of ~SegSCStdOptics *****************/

void ACStdOptics::printACStdOptics() {

    bool debug = false;
    if (debug) {
      *oStr << "  -- DCStdOptics::printDCStdOptics" << endl;
      *oStr << "       iPrtMode " << iPrtMode << endl;
    }
    if (iPrtMode == 0) return;
    *oStr << " -- GDCTelescopeFactory  DCStdOptics::printDCStdOptics" << endl;
    *oStr << "       stdType:        " << getTelType(stdType) << endl;
    *oStr << "       defoc           " << defoc << endl;
    *oStr << "       foclength       " << foclength << endl;
    *oStr << "       camRadius       " << camRadius << endl;
    *oStr << " plateScaleFactor      " << plateScaleFactor << endl;
    *oStr << "       radius          " << radius << endl;
    *oStr << "       rotation_offset " << rotation_offset << endl;
    *oStr << "       xoff            " << xoff << endl;
    *oStr << "       yoff            " << yoff << endl;
    *oStr << "       nbr_mir         " << nbr_mir << endl;
    *oStr << "       nbinsx          " << nbinsx  << endl;
    *oStr << "       nbinsy          " << nbinsy  << endl;
    *oStr << "       gridOption      " << gridOption << endl;
    *oStr << "       gridFile        " << gridFile << endl;
    *oStr << "       rayTracer       " << getRayTracerType(rayTracerEnum)
      << endl;
    *oStr << "       facet.size()    " << facet.size() << endl;
    *oStr << "       avgTransitTime  " << fAvgTransitTime <<endl;
    if (grid == 0) {
      *oStr << "      grid not instantiated " << endl;
    }

    if ( iPrtMode > 2 ) {
      for (unsigned i=0;i<facet.size(); i++) {
        facet[i]->printDCStdFacet(*oStr);
        *oStr << endl;
      }
    }

    if (iPrtMode > 3 ) {

    }
    *oStr << endl;
};
/************** end of printSegSCStdOptics ***********************/

GACTelescopeFactory::
GACTelescopeFactory(GReadACStd &acReader,
		       const string &editPilotFile) {
  
  printParameters = false;
  bool debug = false;
  sPilotEdit = editPilotFile;
  if (debug) {
    *oLog << "  -- GACTelescopeFactory Constructor:  " << sPilotEdit << endl;
  }
  iNumACTelMade = 0;
  readAC = 0;
  pi = 0;
  sPilotEdit = "";
  opt = 0;
  mVReflWaveLgts = 0;
  mVCoeffs = 0;
  ACTel = 0;

  readAC = &(acReader);
  readAC->setACTelescopeFactory(this);

  sPilotEdit = editPilotFile;

  pi = new GPilot(sPilotEdit);
 
  if (debug) {
    *oLog << "    -- end of GACTelescopeFactory constructor" << endl;
  }
  
};

/************** end of GACTelescopeFactory ***********************/

GACTelescopeFactory::~GACTelescopeFactory() {
 
  bool debug = false;
  if (debug) {
    *oLog << "  -- GACTelescopeFactory::~GACTelescopeFactory" << endl;
  }

  for (itmStdOp=mStdOptics.begin();
     itmStdOp!=mStdOptics.end(); itmStdOp++) {
    // WHY CAN'T I DELETE THIS CLASS? thinks I have already deleted the class?
    // was created in reader: SCFac->mStdOptics[i] = new ACStdOptics();
    //ACStdOptics *tmpSeg;
    //tmpSeg = itmStdOp->second;
    //*oLog << "tmpSeg->stdNum " <<tmpSeg->stdNum << endl;
    //delete tmpSeg;
    SafeDelete( itmStdOp->second );
    //*oLog << " (itmStdOp->second)->stdNum " <<  (itmStdOp->second)->stdNum << endl;
    //*oLog << "itmStdOp->second " << itmStdOp->second << endl;
  }

  SafeDelete(pi);
  SafeDelete(readAC);
// reflection coeffs owned by the factory so delete here
  SafeDelete(mVReflWaveLgts);
  SafeDelete(mVCoeffs);
    
  map<int,GGeometryBase*>::iterator iter;
    for (iter = mDCGeo.begin(); iter != mDCGeo.end();iter++) {
      SafeDelete(iter->second);
    }
 
};
/************** end of ~GACTelescopeFactory ***********************/

GACTelescope* GACTelescopeFactory::makeTelescope(const int &id,
                                                     const int &std) {
  
  int debug = true;
  if (debug) {
    *oLog << " -- GACTelescopeFactory::makeTelescope" << endl;
    *oLog << "      telID  = " << id << endl;
    *oLog << "      telStd = " << std << endl;
  }

  Int_t idTel = id;
  Int_t iStdID = std;
  iNumACTelMade++; // increment number of AC telescopes made by factory
  
  // get parameters for this telescope
  itmStdOp = mStdOptics.find(iStdID);
  assert(itmStdOp != mStdOptics.end());
  // make pointer to working stdoptics structure, easier typing
  opt = mStdOptics[iStdID];

  // make the telescope
  ACTel = new GACTelescope;

  ACTel->iPrtMode = opt->iPrtMode;

  
  ACTel->setTelID(idTel);
  ACTel->setStdID(iStdID);
  
  ACTel->eTelType = opt->stdType;
  ACTel->fAvgTransitTime = opt->fAvgTransitTime;
  ACTel->fRotationOffset = opt->fRotationOffset;
  ACTel->fPlateScaleFactor = opt->fPlateScaleFactor;

  // general telescope parameters
  ACTel->fF = (opt->fF);
  ACTel->fAlpha = opt->fAlpha;
  ACTel->fQ = opt->fQ;

  // primary parameters
  ACTel->fRpMax = (opt->fRpMax);
  ACTel->fRpMin = (opt->fRpMin);
  ACTel->fZp = (opt->fZp);
  ACTel->fNp = opt->iNParP;
  ACTel->fzp = opt->fzp;
  // make primary poly coefficients
  ACTel->fP = new Double_t[ACTel->fNp];
  Double_t fF = (ACTel->fF)*m;
  
  (ACTel->fP)[0] = TMath::Power(fF,  1)* ((opt->fzp[0]));

  
  if (printParameters) {
    *oLog << endl;
    *oLog << "  calculation of primary parameters " << endl;
    *oLog << "           fzp[i]: input parameters from configuration file " << endl;
    *oLog << "           fP[i]  =  TMath::Power(fF,  iPowerF)* fzp[i] " << endl;  
    *oLog << "  i     iPowerF     fzp[i]           fP[i] " << endl;
    *oLog <<"  0" << "        " << "1" << "       " << opt->fzp[0] << "       " 
          << (ACTel->fP)[0] << endl;
  }
  int iPowerF = -1;
  for (int i = 1;i < ACTel->fNp; i++) {
    (ACTel->fP)[i] = TMath::Power(fF,  iPowerF)* ((opt->fzp[i]));
    if (printParameters) {
      *oLog << "  " <<  i << "       " << iPowerF
            << "       " << opt->fzp[i] << "       " << (ACTel->fP)[i] << endl;
    }
    iPowerF = iPowerF -2;
  }
  if(printParameters) *oLog << endl;
  
  // secondary parameters
  ACTel->fRsMax = (opt->fRsMax);
  ACTel->fRsMin = (opt->fRsMin);
  ACTel->fZs = (opt->fZs);
  ACTel->fNs = opt->fNs;
  ACTel->fS = opt->fS;
  ACTel->iNParS = opt->iNParS;
  ACTel->fNs = opt->iNParS;
  ACTel->fzs = opt->fzs;
  // make secondary poly coefficients
  ACTel->fS = new Double_t[ACTel->iNParS];
  fF = (ACTel->fF)*m;

  (ACTel->fS)[0] = TMath::Power(fF,  1)* ((opt->fzs[0]));

  if (printParameters) {
    *oLog << "  calculation of secondary parameters " << endl;
    *oLog << "           fzs[i]: input parameters from configuration file " << endl;
    *oLog << "           fS[i]  =  TMath::Power(fF,  iPowerF)* fzs[i] " << endl;  
    *oLog << "  i     iPowerF     fzs[i]           fS[i] " << endl;
    *oLog <<"  0" << "        " << "1" << "       " << opt->fzs[0] << "       " 
          << (ACTel->fS)[0] << endl;
  }

  iPowerF = -1;
  for (int i = 1;i < ACTel->fNs; i++) {
    (ACTel->fS)[i] = TMath::Power(fF,  iPowerF)* ((opt->fzs[i]));
    if (printParameters) {
      *oLog << "  " <<  i << "       " << iPowerF
            << "       " << opt->fzs[i] << "       " << (ACTel->fS)[i] << endl;
    }
    iPowerF = iPowerF -2;
  }
  if (printParameters) *oLog << endl;

  // primary segment details
  ACTel->iNumP1Mirrors = opt->iNumP1Mirrors;
  ACTel->iNumP2Mirrors = opt->iNumP2Mirrors;
  ACTel->vSegP1 = opt->vSegP1;
  ACTel->vSegP2 = opt->vSegP2;

 // secondary segment S1details
  ACTel->iNumS1Mirrors = opt->iNumS1Mirrors;
  ACTel->iNumS2Mirrors = opt->iNumS2Mirrors;
  ACTel->vSegS1 = opt->vSegS1;
  ACTel->vSegS2 = opt->vSegS2;

  // focal plane
  ACTel->fKappa1 = opt->fKappa1;
  ACTel->fKappa2 = opt->fKappa2;
  ACTel->fRf     = opt->fRf;
  ACTel->fZf     = opt->fZf;

  // primary baffle
  ACTel->bpBaffleFlag = opt->bpBaffleFlag;
  ACTel->fpBRadOffset = opt-> fpBRadOffset;
  ACTel->fpBLen = opt-> fpBLen;
  ACTel->fpBZOffset = opt-> fpBZOffset;
  ACTel->fpBTilt = opt-> fpBTilt;

  // secondary baffle
  ACTel->bsBaffleFlag = opt->bsBaffleFlag;
  ACTel->fsBRadOffset = opt-> fsBRadOffset;
  ACTel->fsBLen = opt-> fsBLen;
  ACTel->fsBZOffset = opt-> fsBZOffset;
  ACTel->fsBTilt = opt-> fsBTilt;
 
  // camera
  ACTel->bCameraFlag = opt->bCameraFlag;
  ACTel->fPixelSize   = opt->fPixelSize;
  ACTel->fSubCells   = opt->fSubCells;
  ACTel->fMAPMTWidth  = opt->fMAPMTWidth;
  ACTel->fMAPMTLength = opt->fMAPMTLength;
  ACTel->fInputWindowThickness   = opt->fInputWindowThickness;
  ACTel->fMAPMTOffset = opt->fMAPMTOffset;
  ACTel->fMAPMTGap    = opt->fMAPMTGap;
  ACTel->fMAPMTRefIndex   = opt->fMAPMTRefIndex;
  ACTel->bSingleMAPMTmodule = opt->bSingleMAPMTmodule;

  // entrance window

  ACTel->bEntranceWindowFlag      = opt->bEntranceWindowFlag;
  ACTel->bEntranceWindowAbsFlag   = opt->bEntranceWindowAbsFlag;
  ACTel->fEntranceWindowThickness = opt->fEntranceWindowThickness;
  ACTel->fEntranceWindowN         = opt->fEntranceWindowN;
  ACTel->fEntranceWindowAbsLength = opt->fEntranceWindowAbsLength;
  ACTel->fEntranceWindowOffset    = opt->fEntranceWindowOffset;

  //

  ACTel->fFocalSurfaceXOffset     = opt->fFocalSurfaceXOffset;
  ACTel->fFocalSurfaceYOffset     = opt->fFocalSurfaceYOffset;
  ACTel->fFocalSurfaceZOffset     = opt->fFocalSurfaceZOffset;
  ACTel->fFocalSurfacePhiOffset   = opt->fFocalSurfacePhiOffset;
  ACTel->fFocalSurfaceThetaOffset = opt->fFocalSurfaceThetaOffset;
  ACTel->fFocalSurfacePsiOffset   = opt->fFocalSurfacePsiOffset;
 
  ACTel->buildTelescope();

  return ACTel;
};
/************** end of makeTelescope ***********************/

void GACTelescopeFactory::editWorkingTelescope(GACTelescope *ACTel1) {
  //ACTel1->printTelescope();  // unused
  bool debug = false;
  if (debug) {
    *oLog << " -- GACTelescopeFactory::editWorkingTelescope" << endl;
  }

 
};
/************** end of editWorkingTelescope ***********************/

void GACTelescopeFactory::printStdTelescope(const int &iStd,
                                               const int &mode,ostream &oStr) {
  (void) iStd;  // unused
  (void) mode; // unused
  oStr << "unused oStr in parameter list" << endl;
  // DO NOT USE.
  bool debug = false;
  if (debug) {
    *oLog << " -- GACTelescopeFactory::printStdTelescope" << endl;
  }
  
};
/************** end of :printStdTelescope ***********************/

void GACTelescopeFactory::setPrintMode(ostream &oStr,
                                          const int prtMode) {
  oStr << "unused ostream in parameter list" << endl;
  (void) prtMode; // unused
  bool debug = false;
  if (debug) {
    *oLog << " -- GACTelescopeFactory::setPrintMode" << endl;
  }
 
}; 
/************** end of setPrintMode  ***********************/



