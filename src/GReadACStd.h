/*
VERSION4.0
30May2016
*/
/*! \brief  GReadSegSCStd class: concrete class for reading 
      standard telescope configurations for use by telescope factory

 */

#ifndef GREADACSTD
#define GREADACSTD

class GPilot; 
class GACTelescopeFactory;

class GReadACStd {

  // set factory as a friend or the reverse?
 private: 

  // used by all telescope types for various printing modes
  //ostream *oPrtStrm;
  //int iPrtMode;

  GACTelescopeFactory *ACFac;

  // for reading pilot file, can be appended files
  string spilotfile;   //!< pilot file
  GPilot *pi;        //!< pilot file reader
  vector<string> tokens; //!< vector of pilot record string tokens
  string flagline;       //!< record flag for pilot file reading

  //map<int, vector<double> > *mVReflWaveLgts;
  //map<int, vector<double> > *mVCoeffs;

  int iStdNum;  //!< active telescope number used in filling entries

  ACStdOptics *opt;  //!< working SCStdOptics from DCFac.

  void setupACFactory();

  void getReflCoeff();

  void getPolyCoeffs();

  void readBasicRecord(const vector<string> &tokens,
                       const MirSeg &eMirPS,
                       ACStdOptics *opt);
 
void readChangeRecord(const vector<string> &tokens,
                       const MirSeg &eMirPS,
                       ACStdOptics *opt);
 public:

  GReadACStd(const string &pilotfile,GACTelescopeFactory *ACFactory );

  GReadACStd(const string &pilotfile);

  ~GReadACStd();

  void setACTelescopeFactory(GACTelescopeFactory *ACFactory);

  //void setPrintMode(ostream &oStr=cout,const int prtMode=0);

};

#endif
