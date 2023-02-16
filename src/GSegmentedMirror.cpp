// New Schwartschild-Couder optical system with segmented mirrors and the
// telescope structure.
// Nov/21/2012 Akira Okumura <oxon@mac.com>

#include "TGeoSphere.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TObjArray.h"
#include "TMath.h"
#include "TVector3.h"

#include "TGeoBBox.h"
#include "TGeoCompositeShape.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoMatrix.h"
#include "TGeoTrd1.h"
#include "TGeoTube.h"
#include "TGeoXtru.h"
#include "TPolyLine3D.h"
#include "TRandom.h"
#include "TROOT.h"

#include "ABorderSurfaceCondition.h"
#include "AGeoAsphericDisk.h"
#include "AGlassCatalog.h"
#include "ALens.h"
#include "AMirror.h"
#include "AObscuration.h"
#include "AOpticsManager.h"
#include "ARay.h"
#include "ARayArray.h"
#include "ARayShooter.h"
#include "GSegmentedMirror.h"

#include <iostream>
#include <string>

// define useful units
static const Double_t cm = AOpticsManager::cm();
static const Double_t mm = AOpticsManager::mm();
static const Double_t um = AOpticsManager::um();
static const Double_t nm = AOpticsManager::nm();
static const Double_t  m = AOpticsManager::m();

//______________________________________________________________________________
SegmentedMirror::SegmentedMirror(Double_t rmin, Double_t rmax,
                                 Double_t phimin, Double_t phimax)
  : fRmin(rmin), fRmax(rmax), fPhimin(phimin), fPhimax(phimax),
    fMargin(0), fPositionErrorX(0), fPositionErrorY(0), fPositionErrorZ(0),
    fRotationErrorXY(0), fRotationErrorSagittal(0), fRotationErrorTangential(0),
    fRoughness(0)
{
}

//______________________________________________________________________________
TGeoCombiTrans* SegmentedMirror::BuildMirrorCombiTrans(const char* name,AGeoAsphericDisk* disk, Bool_t isPrimary)
{

  // Double_t cphi = (fPhimax + fPhimin)/2.;
  // Double_t cr = (fRmax + fRmin)/2.;
  // Double_t cz = isPrimary ? disk->CalcF2(cr) : disk->CalcF1(cr);

  // // TGeoRotation rot1("", -90. + cphi + fRotationErrorXY, fRotationErrorSagittal, 0);
  // // TGeoRotation rot2("", -cphi, 0., 0.);
  // // TGeoRotation rot3("", 0., fRotationErrorTangential, 0.);
  // // TGeoRotation rot4("", cphi, 0., 0.);

  // // TGeoTranslation tr4(cr*TMath::Cos(cphi*TMath::DegToRad()) + fPositionErrorX, cr*TMath::Sin(cphi*TMath::DegToRad()) + fPositionErrorY, cz + fPositionErrorZ);
 
  // TGeoRotation rot1("", 0, 0, 0);

  // TGeoTranslation tr4(0, 0, 0);
  // TGeoCombiTrans* combi = new TGeoCombiTrans(tr4, rot1);
  // // combi->RegisterYourself();
  // return combi;
  // TGeoHMatrix* hmat = new TGeoHMatrix((*combi) * (*transZ));

  // return combi;
  Double_t kF = 148.5 * cm;        // focal length
    Double_t kMirrorR = kF * 2;      // the radius of curvature
    Double_t kMirrorD = 15 * cm;     // facet diameter, circular mirror
    Double_t kMirrorT = 1.2954 * cm; // mirror thickness
    Double_t theta = TMath::ASin(kMirrorD / 2. / kMirrorR) * TMath::RadToDeg();

    const int kNMirror = 84;
    
    Double_t xy[kNMirror][3] = {{30.2514, -8.001, 3.3274}, {30.2514, 8.001, 3.3274},
        {44.1198, -16.002, 7.62}, {44.1198, 0, 6.7056}, {44.1198, 16.002, 7.62},
        {57.9628, -24.003, 13.8938}, {57.9628, -8.001, 12.0142}, {57.9628, 8.001, 12.0142}, {57.9628, 24.003, 13.8938},
        {71.8312, -32.004, 22.5298}, {71.8312, -16.002, 19.5072}, { 71.8312, 0, 18.5166}, {71.8312, 16.002, 19.5072}, {71.8312, 32.004, 22.5298},
        {22.054769255679297, 22.197980900044563, 3.3274}, {8.19663074432071, 30.198980900044567, 3.3274}, {35.91803851135859, 30.20786760988867, 7.62}, {22.059900000000003, 38.20886760988867, 6.7056}, {8.201761488641417, 46.209867609888676, 7.62}, {49.76860776703789, 38.19575727447665, 13.8938}, {35.9104692556793, 46.19675727447665, 12.0142}, {22.052330744320713, 54.19775727447666, 12.0142}, {8.194192232962127, 62.19875727447666, 13.8938}, {63.63187702271718, 46.20564398432076, 22.5298}, {49.77373851135859, 54.20664398432076, 19.5072}, {35.9156, 62.20764398432076, 18.5166}, {22.05746148864142, 70.20864398432076, 19.5072}, {8.199322977282833, 78.20964398432076, 22.5298},
        {-8.196630744320698, 30.19898090004457, 3.3274}, {-22.05476925567929, 22.197980900044573, 3.3274}, {-8.201761488641399, 46.20986760988868, 7.62}, {-22.0599, 38.208867609888685, 6.7056}, {-35.918038511358574, 30.207867609888687, 7.62}, {-8.194192232962102, 62.19875727447666, 13.8938}, {-22.05233074432069, 54.197757274476665, 12.0142}, {-35.91046925567928, 46.196757274476674, 12.0142}, {-49.76860776703787, 38.19575727447668, 13.8938}, {-8.199322977282804, 78.20964398432076, 22.5298}, {-22.057461488641394, 70.20864398432077, 19.5072}, {-35.9156, 62.207643984320775, 18.5166}, {-49.77373851135857, 54.20664398432078, 19.5072}, {-63.63187702271716, 46.20564398432079, 22.5298},
        {-30.2514, 8.001, 3.3274}, {-30.2514, -8.001, 3.3274}, {-44.1198, 16.002, 7.62}, {-44.1198, 0, 6.7056}, {-44.1198, -16.002, 7.62}, {-57.9628, 24.003, 13.8938}, {-57.9628, 8.001, 12.0142}, {-57.9628, -8.0001, 12.0142}, {-57.9628, -24.003, 13.8938}, {-71.8312, 32.004, 22.5298}, {-71.8312, 16.002, 19.5072}, {-71.8312, 0, 18.5166}, {-71.8312, -16.002, 19.5072}, {-71.8312, -32.004, 22.5298},
        {-22.054769255679304, -22.197980900044556, 3.3274}, {-8.196630744320723, -30.19898090004456, 3.3274}, {-35.9180385113586, -30.20786760988866, 7.62}, {-22.0599, -38.208867609888664, 6.7056}, {-8.201761488641438, -46.20986760988867, 7.62}, {-49.7686077670379, -38.195757274476634, 13.8938}, {-35.91046925567932, -46.196757274476646, 12.0142}, {-22.052330744320734, -54.19775727447665, 12.0142}, {-8.194192232962152, -62.19875727447666, 13.8938}, {-63.63187702271719, -46.20564398432073, 22.5298}, {-49.77373851135861, -54.20664398432074, 19.5072}, {-35.9156, -62.20764398432075, 18.5166}, {-22.057461488641444, -70.20864398432076, 19.5072}, {-8.199322977282861, -78.20964398432076, 22.5298},
        {8.196630744320707, -30.198980900044564, 3.3274}, {22.054769255679293, -22.197980900044566, 3.3274}, {8.201761488641413, -46.20986760988867, 7.62}, {22.0599, -38.20886760988867, 6.7056}, {35.91803851135859, -30.207867609888673, 7.62}, {8.19419223296212, -62.198757274476655, 13.8938}, {22.052330744320706, -54.19775727447666, 12.0142}, {35.910469255679295, -46.19675727447665, 12.0142}, {49.76860776703788, -38.195757274476655, 13.8938}, {8.199322977282826, -78.20964398432076, 22.5298}, {22.05746148864141, -70.20864398432076, 19.5072}, {35.9156, -62.20764398432076, 18.5166}, {49.77373851135859, -54.20664398432076, 19.5072}, {63.63187702271717, -46.20564398432076, 22.5298}
    };
    // clang-format on

    for (int i = 0; i < kNMirror; i++) {
      char* mirname = Form("%d",i);
      if(strcmp(name,mirname)==0){
        std::cout<< mirname << std::endl;
        Double_t x = xy[i][0] * cm;
        Double_t y = xy[i][1] * cm;
        Double_t z = xy[i][2] * cm;
        Double_t r2d = TMath::RadToDeg();
        Double_t r2 = TMath::Power(x, 2) + TMath::Power(y, 2);
      
        // each mirror center is relocated from the origin (0, 0, 0) to (x, y, z)
        TGeoTranslation* trans =
            new TGeoTranslation(Form("mirTrans%d", i), x, y, z);
        trans->RegisterYourself();

        // and is rotated to compose a DC optics
        Double_t phi = (TMath::ATan2(y, x)) * r2d;
        theta = (TMath::ATan2(TMath::Sqrt(r2), 2 * kF - z))* r2d;
        TGeoRotation* rot = new TGeoRotation("", phi - 90., theta, 0);

        // make a matrix from translation and rotation matrices
        TGeoTranslation* transZ = new TGeoTranslation(0, 0, kMirrorR);
        TGeoCombiTrans* combi = new TGeoCombiTrans(*trans, *rot);
        combi->RegisterYourself();
        return combi;
        TGeoHMatrix* hmat = new TGeoHMatrix((*combi) * (*transZ));
      }
    }
    
}

//______________________________________________________________________________
void SegmentedMirror::SetMargin(Double_t margin)
{
  // Set the margin around the mirror
  fMargin = margin;
}

//______________________________________________________________________________
void SegmentedMirror::SetPositionErrors(Double_t x, Double_t y, Double_t z)
{
  // Set the position errors along X/Y/Z axes
  fPositionErrorX = x;
  fPositionErrorY = y;
  fPositionErrorZ = z;
}

//______________________________________________________________________________
void SegmentedMirror::SetRotationErrors(Double_t xy, Double_t sag, Double_t tan)
{
  // Set the rotation errors in XY/sagittal/tangential planes
  fRotationErrorXY = xy;
  fRotationErrorSagittal = sag;
  fRotationErrorTangential = tan;
}

//______________________________________________________________________________
void SegmentedMirror::SetRougness(Double_t roughness)
{
  // Set the mirror roughness in unit of degree.
  fRoughness = roughness;
}

SectorSegmentedMirror::SectorSegmentedMirror(Double_t rmin, Double_t rmax, Double_t phimin, Double_t phimax)
  : SegmentedMirror(rmin, rmax, phimin, phimax)
{
}

//______________________________________________________________________________
AMirror* SectorSegmentedMirror::BuildMirror(const char* name,
                                            AGeoAsphericDisk* disk,
                                            Bool_t isPrimary)
{
  Double_t zmin = disk->GetOrigin()[2] - disk->GetDZ();
  Double_t zmax = disk->GetOrigin()[2] + disk->GetDZ();
  Double_t dphi = (fPhimax - fPhimin)/2.;

  // (R, Z) = (cr, cz) is used as the origin of the rotation error matrix
  Double_t cr = (fRmax + fRmin)/2.;
  Double_t cz = isPrimary ? disk->CalcF2(cr) : disk->CalcF1(cr);

  // Define the base of the segment shape
  TGeoTubeSeg* seg1 = new TGeoTubeSeg(Form("%s_seg1", name),
                                      fRmin + fMargin, fRmax - fMargin,
                                      (zmax - zmin)/2., -dphi + 90., dphi + 90.);

  // To be used for the margin cut of the side edges
  TGeoTubeSeg* seg2 = new TGeoTubeSeg(Form("%s_seg2", name),
                                      0, fRmax - fMargin,
                                      (zmax - zmin)/2., -dphi + 90., dphi + 90.);

  TGeoTranslation* tr1
    = new TGeoTranslation(Form("%s_tr1", name), 0, -cr, (zmax + zmin)/2. - cz);
  tr1->RegisterYourself();

  Double_t d = fMargin/TMath::Sin(dphi*TMath::DegToRad());
  TGeoTranslation* tr2
    = new TGeoTranslation(Form("%s_tr2", name), 0, d - cr, (zmax + zmin)/2. - cz);
  tr2->RegisterYourself();

  TGeoTranslation* tr3 = new TGeoTranslation(Form("%s_tr3", name), 0, -cr, -cz);
  tr3->RegisterYourself();

  TGeoCompositeShape* comp
    = new TGeoCompositeShape(Form("%s_comp", name),
                             Form("%s:%s*%s:%s*%s:%s",
                                  seg1->GetName(), tr1->GetName(),
                                  seg2->GetName(), tr2->GetName(),
                                  disk->GetName(), tr3->GetName()));

  AMirror* mirror = new AMirror(Form("%s_mirror", name), comp);
  
  return mirror;
}

//______________________________________________________________________________
TetragonSegmentedMirror::TetragonSegmentedMirror(Double_t rmin, Double_t rmax, Double_t phimin, Double_t phimax) :
  SegmentedMirror(rmin, rmax, phimin, phimax)
{
}

//______________________________________________________________________________
AMirror* TetragonSegmentedMirror::BuildMirror(const char* name,
                                              AGeoAsphericDisk* disk,
                                              Bool_t isPrimary)
{
  Double_t zmin = disk->GetOrigin()[2] - disk->GetDZ();
  Double_t zmax = disk->GetOrigin()[2] + disk->GetDZ();
  Double_t dphi = (fPhimax - fPhimin)/2.;

  // (R, Z) = (cr, cz) is used as the origin of the rotation error matrix
  Double_t cr = (fRmax + fRmin)/2.;
  Double_t cz = isPrimary ? disk->CalcF2(cr) : disk->CalcF1(cr);
  
  Double_t kF = 148.5 * cm;        // focal length
  Double_t kMirrorR = kF * 2;      // the radius of curvature
  Double_t kMirrorD = 15 * cm;     // facet diameter, circular mirror
  Double_t kMirrorT = 1.2954 * cm; // mirror thickness
  Double_t theta = TMath::ASin(kMirrorD / 2. / kMirrorR) * TMath::RadToDeg();

  // clang-format on
    TGeoSphere* mirSphere = new TGeoSphere("mirSphere", kMirrorR, kMirrorR + kMirrorT, 180. - theta, 180.);
    // AGeoAsphericDisk * mirSphere = new AGeoAsphericDisk("mirSphere", 0, 0, 0, 0, (7.5 * cm));
    AMirror* mirror = new AMirror("mirror", mirSphere);
    return mirror;
  //TGeoSphere* mirSphere = new TGeoSphere("mirSphere", kMirrorR, kMirrorR + kMirrorT, 180. - theta, 180.);
  
}
//______________________________________________________________________________
PentagonSegmentedMirror::PentagonSegmentedMirror(Double_t rmin, Double_t rmax, Double_t phimin, Double_t phimax) :
  SegmentedMirror(rmin, rmax, phimin, phimax)
{
}

//______________________________________________________________________________
AMirror* PentagonSegmentedMirror::BuildMirror(const char* name,
                                              AGeoAsphericDisk* disk,
                                              Bool_t isPrimary)
{
  Double_t zmin = disk->GetOrigin()[2] - disk->GetDZ();
  Double_t zmax = disk->GetOrigin()[2] + disk->GetDZ();
  Double_t dphi = (fPhimax - fPhimin)/2.;

  // (R, Z) = (cr, cz) is used as the origin of the rotation error matrix
  Double_t cr = (fRmax + fRmin)/2.;
  Double_t cz = isPrimary ? disk->CalcF2(cr) : disk->CalcF1(cr);

  // Define the base of the segment shape
  // ____
  // ___ ----____
  // C  ---- B    Rmax  
  //p2---p1/
  //    / /
  // --p0/
  // ----____
  //     A   ---- Rmin
  //
  Double_t d2r = TMath::DegToRad();
  Double_t ax = fRmin*TMath::Cos((90. - dphi)*d2r);
  Double_t ay = fRmin*TMath::Sin((90. - dphi)*d2r);
  //Double_t bx = fRmax*TMath::Cos((90. - dphi)*d2r);
  //Double_t by = fRmax*TMath::Sin((90. - dphi)*d2r);
  Double_t cy = fRmax;

  Double_t p0x = ax - fMargin*(1./TMath::Sin(dphi*d2r) - 1.)/TMath::Tan((90. - dphi)*d2r);
  Double_t p0y = ay + fMargin;

  Double_t t1 = TMath::Tan((90. - dphi)*d2r);
  Double_t t2 = TMath::Tan((90. - dphi/2.)*d2r);
  Double_t s1 = TMath::Sin(dphi*d2r);
  Double_t s2 = TMath::Sin((90. - dphi/2.)*d2r);
  Double_t p1x = 1/(t1 + 1/t2)*(cy - fMargin/s2 + t1*ax - ay - fMargin/s1);
  Double_t p1y = t1*(p1x - ax) + ay + fMargin/s1;
  Double_t p2x = 0;
  Double_t p2y = cy - fMargin/s2;

  TGeoXtru* xtru1 = new TGeoXtru(2);
  xtru1->SetName(Form("%s_xtru1", name));
  Double_t x[5] = {p0x, p1x, p2x, -p1x, -p0x};
  Double_t y[5] = {p0y, p1y, p2y, p1y, p0y};
  xtru1->DefinePolygon(5, x, y);
  xtru1->DefineSection(0, zmin);
  xtru1->DefineSection(1, zmax);

  TGeoTranslation* tr1 = new TGeoTranslation(Form("%s_tr1", name), 0, -cr, -cz);
  tr1->RegisterYourself();
  
  TGeoCompositeShape* comp
    = new TGeoCompositeShape(Form("%s_comp", name),
                             Form("%s:%s*%s:%s",
                                  disk->GetName(), tr1->GetName(),
                                  xtru1->GetName(), tr1->GetName()));
  
  AMirror* mirror = new AMirror(Form("%s_mirror", name), comp);
  
  return mirror;
}
