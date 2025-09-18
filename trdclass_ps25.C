#define trdclass_ps25_cxx
#include "trdclass_ps25.h"
#include "PlotLib.C"
#include "GNN/gnn_model.h"
#include "GNN/gnn_model.cpp"
#include "GNN/toGraph.cpp"

#define NPRT 1000
#define USE_TRK
#define MAX_PRINT 10
//
#define USE_GNN  1
#define USE_FIT  1
#define USE_CLUST 1
#define USE_PULSE 0
//
#define USE_125_RAW
//#define USE_250_PULSE
#define MAX_CLUST 500
#define MAX_NODES 100
#define USE_MAXPOS 1
//
#define SAVE_TRACK_HITS
#define SAVE_PDF
//#define WRITE_CSV
#define DEBUG 0
//
int runSwitch=6080; //Change setup from Quad-GEMTRD to MMG1-TRD
int timeSwitchRun1=6155; //RunNum where timing window was changed to 250
int timeSwitchRun2=6210; //RunNum where timing window was changed back to 200
int runAluminum=6357; //Beginning of runs with al cathodes
double firstTimeWin=200.;
double secondTimeWin=250.;
int noRadList[] = {6314,6315,6316,6317,6318,6319,6320,6389,6390,6391,6392,6393,6394};
int listSize = sizeof(noRadList) / sizeof(noRadList[0]);

//-- For single evt clustering display, uncomment BOTH:
//#define SHOW_EVTbyEVT
//#define SHOW_EVT_DISPLAY

void WriteToCSV(std::ofstream &csvFile, float v1, float v2, float v3, float v4, float v5, float v6, float v7, float v8, float v9, float v10) {
  csvFile<<v1<<","<<v2<<","<<v3<<","<<v4<<","<<v5<<","<<v6<<","<<v7<<","<<v8<<","<<v9<<","<<v10<<std::endl;
}

//===================================
//      TRD DAQ Channel Mapping
//===================================

// -- Triple-GEMTRD mapping --
int Get3GEMChan(int ch, int slot, int runNum) {
  int cardNumber = ch/24;
  int cardChannel = ch-cardNumber*24;
  float dchan = cardChannel+cardNumber*24+(slot-3)*72.;
  if ((slot==6 && ch>23) || slot==7 || slot==8 || (slot==9 && ch<48)) {
    if ( (runNum<runAluminum && ((dchan-240.)==16. || (dchan-240.)==24. || (dchan-240.)==51. || (dchan-240.)==52. || (dchan-240.)==55. || (dchan-240.)==179. || (dchan-240.)==193. || (dchan-240.)==212. || (dchan-240.)==224.)) || (runNum>=runAluminum && ((dchan-240.)==0. || (dchan-240.)==16. || (dchan-240.)==2. || (dchan-240.)==24. || (dchan-240.)==51. || (dchan-240.)==52. || (dchan-240.)==55. || (dchan-240.)==56. || (dchan-240.)==179. || (dchan-240.)==193. || (dchan-240.)==126. || (dchan-240.)==212. || (dchan-240.)==224.))) { return -1; } else { return dchan-240.; }
  }
  return -1;
}

// -- Quad-GEMTRD mapping --
int Get4GEMChan(int ch, int slot, int runNum) {
  int cardNumber = ch/24;
  int cardChannel = ch-cardNumber*24;
  float dchan = cardChannel+cardNumber*24+(slot-3)*72.;
  if (runNum<runSwitch) {
    if (slot==3 || slot==4 || slot==5 || (slot==6 && ch<24)) {
      /*if (dchan==16. || dchan==23. || dchan==35. || dchan==45. || dchan==224. || dchan==226.) { return -1; } else { */return dchan; }
    //}
  }
  return -1;
}

// -- MMG1TRD mapping --
int GetMMG1Chan(int ch, int slot, int runNum) {
  int cardNumber = ch/24;
  int cardChannel = ch-cardNumber*24;
  float dchan = cardChannel+cardNumber*24+(slot-3)*72.;
  if (runNum>runSwitch) {
    if (slot==3 || slot==4 || slot==5 || (slot==6 && ch<24)) {
      if ((runNum<runAluminum && (dchan==16. || dchan==45. || dchan==71. || dchan==136. || dchan==215. || dchan==224. || dchan==226.)) || (runNum>=runAluminum && (dchan==16. || dchan==45. || dchan==71. || dchan==136. || dchan==215. || dchan==224. || dchan==226.))) { return -1; } else { return dchan; }
    }
  }
  return -1;
}

// -- uRWELL+GEM mapping --
int GetRWELLXChan(int ch, int slot, int runNum) {
  if (slot>10) slot=(slot-2);
  int cardNumber = ch/24;
  int cardChannel = ch-cardNumber*24;
  cardChannel = 24-cardChannel; //new style electronics have inverted polarity
  float dchan = cardChannel+cardNumber*24+(slot-3)*72.;
    if (slot==10 || (slot==11 && ch<48)) {
      if ((runNum<runAluminum && (dchan-504.)==99.) || (runNum>runAluminum && (dchan-504.)==102.)) { return -1; } else { return dchan-504.; }
    }
  return -1;
}

int GetRWELLYChan(int ch, int slot, int runNum) {
  if (slot>10) slot=(slot-2);
  int cardNumber = ch/24;
  int cardChannel = ch-cardNumber*24;
  cardChannel = 24-cardChannel; //new style electronics have inverted polarity
  float dchan = cardChannel+cardNumber*24+(slot-3)*72.;
    if ((slot==11 && ch>47)) {
      /*if ((dchan-480.)==42.) { return -1; } else {*/ return dchan-576.;//}
    }
  return -1;
}

//============ END DAQ Channel Mapping ============

//============ ADC Gain Calibrations ============
float Get3GEMCalib(float amp, int ch, int runNum) {
  float calibrated_amp = amp;
  if ((runNum>=6295 && runNum<=6297) || (runNum>=6315 && runNum<=6317)) { //LG Cu
    double gemCoefficients[] = {};
    calibrated_amp = amp*gemCoefficients[ch];
  }
  else if ((runNum>=6301 && runNum<=6304) || (runNum>=6318 && runNum<=6320)) { //HG Cu
    double gemCoefficients[] = {};
    calibrated_amp = amp*gemCoefficients[ch];
  }
  else if ((runNum>=6381 && runNum<=6385) || (runNum>=6392 && runNum<=6394)) { //LG Al
    double gemCoefficients[] = {};
    calibrated_amp = amp*gemCoefficients[ch];
  }
  else if ((runNum>=6389 && runNum<=6391) || (runNum>=6385 && runNum<=6388)) { //HG Al
    double gemCoefficients[] = {};
    calibrated_amp = amp*gemCoefficients[ch];
  }
  return calibrated_amp;
}

float GetMMGCalib(float amp, int ch, int runNum) {
  float calibrated_amp = amp;
  if ((runNum>=6295 && runNum<=6297) || (runNum>=6315 && runNum<=6317)) { //LG Cu
    double mmgCoefficients[] = {};
    calibrated_amp = amp*mmgCoefficients[ch];
  }
  else if ((runNum>=6301 && runNum<=6304) || (runNum>=6318 && runNum<=6320)) { //HG Cu
    double mmgCoefficients[] = {};
    calibrated_amp = amp*mmgCoefficients[ch];
  }
  else if ((runNum>=6381 && runNum<=6385) || (runNum>=6392 && runNum<=6394)) { //LG Al
    double mmgCoefficients[] = {};
    calibrated_amp = amp*mmgCoefficients[ch];
  }
  else if ((runNum>=6389 && runNum<=6391) || (runNum>=6385 && runNum<=6388)) { //HG Al
    double mmgCoefficients[] = {};
    calibrated_amp = amp*mmgCoefficients[ch];
  }
  return calibrated_amp;
}

//TODO: uRWELL Gain Calib

//============ END ADC Gain Calibrations ============

//--UNDO SRS strip-to-mm conversion from JANA2
/*int CalcAPVChannel(int peakPosition) {
  double pitch = 0.4;
  return (peakPosition - (-0.5*(102.4 - pitch))) / pitch;
}
*/
void trdclass_ps25::Loop() {
  
  if (fChain == 0) return;

  //==================================================================================================
  //            B o o k    H i s t o g r a m s
  //==================================================================================================
  
  gErrorIgnoreLevel = kBreak; // Suppress warning messages from empty chi^2 fit data
  TList *HistList = new TList();
  #ifdef WRITE_CSV
    std::ofstream csvFile("EventByEvent.csv");
  #endif
  //-----------------  (canvas 0) Event Display ----------
  #ifdef SHOW_EVT_DISPLAY
    /*char c0Title[256];
    sprintf(c0Title,"Event_Display_Run=%d",RunNum);
    TCanvas *c0 = new TCanvas("DISP",c0Title,1100,200,1500,1300);
    c0->Divide(4,3); c0->cd(1);*/
    //-----------------  canvas 2 FPGA Display ----------
    char c2Title[256];
    sprintf(c2Title,"FPGA_Event_Display_Run=%d",RunNum);
    TCanvas *c2 = new TCanvas("FPGA",c2Title,1000,100,1500,1000);
    c2->Divide(5,2); c2->cd(1);
/*    char c3Title[256];
    sprintf(c3Title,"DQM_Event_Display_Run=%d",RunNum);
    TCanvas *c3 = new TCanvas("PID",c3Title,1000,100,1100,900);
    c3->Divide(5,1); //c3->cd(1);
*/
  #endif
  
  //RunNum-based timing window change
  double histTime = -1;
  if (RunNum<=timeSwitchRun1 || RunNum>=timeSwitchRun2) { histTime = firstTimeWin; } else { histTime = secondTimeWin; }
  // Track fit in time
  TF1 fx1("fx1","pol1",40,150);
  TF1 fx2("fx2","pol1",40,150);
  f125_fit = new TH2F("f125_fit","Triple GEM-TRD Track Fit; Time Response (*8ns) ; X Channel",(int)histTime,0.5,histTime+0.5,240,-0.5,239.5);
  mmg1_f125_fit = new TH2F("mmg1_f125_fit","MMG1-TRD Track Fit; Time Response (*8ns) ; X Channel",(int)histTime,0.5,histTime+0.5,240,-0.5,239.5);
  urw_f125_fit = new TH2F("urw_f125_fit","uRWELL-TRD Track Fit; Time Response (*8ns) ; X Channel",(int)histTime,0.5,histTime+0.5,120,-0.5,119.5);
  //-- TRD - GEMTRKR alignment --------
  double xx1=-37., yy1=-55.,  xx2=53., yy2=44.;
  double aa=(yy2-yy1)/(xx2-xx1);
  double bb=yy1-aa*xx1;
  TF1 ftrk("ftrk","[0]*x+[1]",-55.,55.);
  ftrk.SetParameter(0,aa);
  ftrk.SetParameter(1,bb);
  TF1 ftrkr("ftrkr","(x-[1])/[0]",0.,255.);
  ftrkr.SetParameter(0,aa);
  ftrkr.SetParameter(1,bb);
  //-----------------------
  float z1 = 0.;
  float z2 = 1043.; //1047.75;
  float z3 = 1125.; //1114.;
  float zgem = 593.; //584;
  float zmmg1 = 303.; //298.5;
  float zurw = 880.; //857;
  
  //----------------  GEM TRK fiducial area selection (box cut) [mm] ----------------------
  double xbc1=-50., xbc2=+50., ybc1=-50., ybc2=+50.;
  //----------------------------------------------------------------------------------
  
  //---- Set ADC Thresholds ----
  float TGEM_THRESH=130.; //175
  float MMG1_THRESH=110.; //150
  float URW_THRESH=120.;
  float TRKR_THRESH=700.;
  
  hcount = new TH1D("hcount","Count",3,0,3);  HistList->Add(hcount);
  hcount->SetStats(0);   hcount->SetFillColor(38);   hcount->SetMinimum(1.);
  #if ROOT_VERSION_CODE > ROOT_VERSION(6,0,0)
    hcount->SetCanExtend(TH1::kXaxis);
  #else
    hcount->SetBit(TH1::kCanRebin);
  #endif
  
  htgem_nhits = new TH1F("hgem_nhits","Triple GEM-TRD X Hits (fADC)",100,0,100);  HistList->Add(htgem_nhits);
  hmmg1_nhits = new TH1F("hmmg1_nhits","MMG1-TRD X Hits (fADC)",100,0,100);  HistList->Add(hmmg1_nhits);
  hurw_nxhits = new TH1F("hurw_nxhits","uRWELL-TRD X Hits (fADC)",100,0,100);  HistList->Add(hurw_nxhits);
  hurw_nyhits = new TH1F("hurw_nyhits","uRWELL-TRD Y Hits (fADC)",100,0,100);  HistList->Add(hurw_nyhits);
  hgt1_nhits = new TH1F("hgt1_nhits","GEM-TRKR1 Hits (SRS)",12,0,12);  HistList->Add(hgt1_nhits);
  hgt2_nhits = new TH1F("hgt2_nhits","GEM-TRKR2 Hits (SRS)",12,0,12);  HistList->Add(hgt2_nhits);
  hgt3_nhits = new TH1F("hgt3_nhits","GEM-TRKR3 Hits (SRS)",12,0,12);  HistList->Add(hgt3_nhits);
  
  cout<<"**************************RunNum="<<RunNum<<endl;
  int nx0=100;
  int mfac=25; //60
  int ny0=256; //240;
  double Ymin=0.;    double Ymax=ny0*0.4;
  double Xmin=0.;    double Xmax=30.; //double Xmax=26.;
  mhevt  = new TH2F("mhevt","MMG1-TRD Event Display; z pos [mm]; y pos [mm]",nx0+mfac,Xmin,Xmax,ny0,Ymin,Ymax); mhevt->SetStats(0); mhevt->SetMaximum(10.); mhevt->SetMinimum(0.);
  mhevtc = new TH2F("mhevtc","Clustering; FADC bins; MMG1 strips",nx0+mfac,-0.5,(nx0+mfac)-0.5,ny0,-0.5,ny0-0.5);  mhevtc->SetStats(0);   mhevtc->SetMinimum(0.07); mhevtc->SetMaximum(40.);
  mhevtf = new TH2F("mhevtf","MMG1: Clusters for FPGA; z pos [mm]; y pos [mm]",nx0+mfac,Xmin,Xmax,ny0,Ymin,Ymax);  mhevtf->SetStats(0); mhevtf->SetMaximum(10.);
  
  hevt = new TH2F("hevt","GEM-TRD Event display; z pos [time *8ns]; y pos [chan #]",nx0,Xmin,Xmax,ny0,Ymin,Ymax); hevt->SetStats(0); hevt->SetMaximum(10.); hevt->SetMinimum(0.);
  hevtc = new TH2F("hevtc","Clustering; FADC bins; GEM strips",nx0,-0.5,nx0-0.5,ny0,-0.5,ny0-0.5); hevtc->SetStats(0);   hevtc->SetMinimum(0.07); hevtc->SetMaximum(40.);
  hevtf = new TH2F("hevtf","GEM: Clusters for FPGA; z pos [mm]; y pos [mm]",nx0,Xmin,Xmax,ny0,Ymin,Ymax);  hevtf->SetStats(0); hevtf->SetMaximum(10.);
  #if (USE_PULSE>0)
    hevtk = new TH2F("hevtk","Event display; z pos [mm]; y pos [mm]",nx0,Xmin,Xmax,ny0,Ymin,Ymax); /*hevtk->SetStats(0);*/ hevtk->SetMaximum(10.);
    hevtck = new TH2F("hevtck","Clustering; FADC bins; GEM strips",nx0,-0.5,nx0-0.5,ny0,-0.5,ny0-0.5);
  #endif
  
  //////
  f125_el = new TH1F("f125_el","Triple GEM-TRD f125 Peak Amp; ADC Amplitude ; Counts ",410,0.,4100.);          HistList->Add(f125_el);
  f125_el_max = new TH1F("f125_el_max","GEM-TRD f125 Max Amp; ADC Amplitude ; Counts ",410,0.,4100.);          HistList->Add(f125_el_max);
  f125_el_max_late = new TH1F("f125_el_max_late","GEM-TRD f125 Max Amp for Electrons ; ADC Amplitude ; Counts ",410,0.,4100.);           HistList->Add(f125_el_max_late);
  mmg1_f125_el = new TH1F("mmg1_f125_el","MMG1-TRD f125 Peak Amp; ADC Amplitude ; Counts ",410,0.,4100.);       HistList->Add(mmg1_f125_el);
  urw_f125_el_x = new TH1F("urw_f125_el_x","uRWELL-TRD X f125 Peak Amp; ADC Amplitude ; Counts ",410,0.,4100.);       HistList->Add(urw_f125_el_x);
  urw_f125_el_y = new TH1F("urw_f125_el_y","uRWELL-TRD Y f125 Peak Amp; ADC Amplitude ; Counts ",410,0.,4100.);       HistList->Add(urw_f125_el_y);
  mmg1_f125_el_max = new TH1F("mmg1_f125_el_max","MMG1-TRD f125 Max Amp; ADC Amplitude ; Counts ",410,0.,4100.);           HistList->Add(mmg1_f125_el_max);
  mmg1_f125_el_max_late = new TH1F("mmg1_f125_el_max_late","MMG1-TRD f125 Max Amp for Electrons ; ADC Amplitude ; Counts ",410,0.,4100.);           HistList->Add(mmg1_f125_el_max_late);
  urw_f125_el_xmax = new TH1F("urw_f125_el_xmax","uRWELL-TRD X f125 Max Amp; ADC Amplitude ; Counts ",410,0.,4100.);           HistList->Add(urw_f125_el_xmax);
  urw_f125_el_xmax_late = new TH1F("urw_f125_el_xmax_late","uRWELL-TRD X f125 Max Amp for Electrons ; ADC Amplitude ; Counts ",410,0.,4100.);           HistList->Add(urw_f125_el_xmax_late);
  urw_f125_el_ymax = new TH1F("urw_f125_el_ymax","uRWELL-TRD Y f125 Max Amp; ADC Amplitude ; Counts ",410,0.,4100.);           HistList->Add(urw_f125_el_ymax);

  // --- SRS ---
  hgemtrkr_1_max_xch = new TH1F("hgemtrkr_1_max_xch"," GEM-TRKR1 Max X Position ; X Chan [mm]",256,-0.2,102.2);  HistList->Add(hgemtrkr_1_max_xch);
  hgemtrkr_1_max_xamp = new TH1F("hgemtrkr_1_max_xamp"," GEM-TRKR1 Max X Amp ; ADC Amp ",410,0.,4100.);  HistList->Add(hgemtrkr_1_max_xamp);
  hgemtrkr_2_max_xch = new TH1F("hgemtrkr_2_max_xch"," GEM-TRKR2 Max X Position ; X Chan [mm]",256,-0.2,102.2);  HistList->Add(hgemtrkr_2_max_xch);
  hgemtrkr_2_max_xamp = new TH1F("hgemtrkr_2_max_xamp"," GEM-TRKR2 Max X Amp ; ADC Amp ",410,0.,4100.);  HistList->Add(hgemtrkr_2_max_xamp);
  hgemtrkr_3_max_xch = new TH1F("hgemtrkr_3_max_xch"," GEM-TRKR3 Max X Position ; X Chan [mm]",256,-0.2,102.2);  HistList->Add(hgemtrkr_3_max_xch);
  hgemtrkr_3_max_xamp = new TH1F("hgemtrkr_3_max_xamp"," GEM-TRKR3 Max X Amp ; ADC Amp ",410,0.,4100.);  HistList->Add(hgemtrkr_3_max_xamp);
  //--GEMTracker 1
  hgemtrkr_1_peak_xy = new TH2F("hgemtrkr_1_peak_xy","GEM-TRKR1 Peak X-Y Correlation; Peak X [mm]; Peak Y [mm] ",256,-0.2,102.2,256,-0.2,102.2);    hgemtrkr_1_peak_xy->SetStats(0); HistList->Add(hgemtrkr_1_peak_xy);
  hgemtrkr_1_max_xy = new TH2F("hgemtrkr_1_max_xy","GEM-TRKR1 X-Y Correlation for Max Hits; Peak X [mm]; Peak Y [mm] ",256,-0.2,102.2,256,-0.2,102.2);    hgemtrkr_1_max_xy->SetStats(0); HistList->Add(hgemtrkr_1_max_xy);
  hgemtrkr_1_peak_x = new TH1F("hgemtrkr_1_peak_x"," GEM-TRKR1 Peak X Position; X [mm] ",256,-0.2,102.2);  HistList->Add(hgemtrkr_1_peak_x);
  hgemtrkr_1_peak_y = new TH1F("hgemtrkr_1_peak_y"," GEM-TRKR1 Peak Y Position; Y [mm] ",256,-0.2,102.2);  HistList->Add(hgemtrkr_1_peak_y);
  hgemtrkr_1_peak_x_height = new TH1F("hgemtrkr_1_peak_x_height"," GEM-TRKR1 Peak Amplitudes in X; ADC Value ",410,0.,4100.);  HistList->Add(hgemtrkr_1_peak_x_height);
  hgemtrkr_1_peak_y_height = new TH1F("hgemtrkr_1_peak_y_height"," GEM-TRKR1 Peak Amplitudes in Y ; ADC Value ",410,0.,4100.);  HistList->Add(hgemtrkr_1_peak_y_height);
  //--GEMTracker 2
  hgemtrkr_2_peak_xy = new TH2F("hgemtrkr_2_peak_xy","GEM-TRKR2 Peak X-Y Correlation; Peak X [mm]; Peak Y [mm] ",256,-0.2,102.2,256,-0.2,102.2); hgemtrkr_2_peak_xy->SetStats(0); HistList->Add(hgemtrkr_2_peak_xy);
  hgemtrkr_2_max_xy = new TH2F("hgemtrkr_2_max_xy","GEM-TRKR2 X-Y Correlation for Max Hits; Peak X [mm]; Peak Y [mm] ",256,-0.2,102.2,256,-0.2,102.2); hgemtrkr_2_max_xy->SetStats(0);   HistList->Add(hgemtrkr_2_max_xy);
  hgemtrkr_2_peak_x = new TH1F("hgemtrkr_2_peak_x"," GEM-TRKR2 Peak X Position; X [mm] ",256,-0.2,102.2);  HistList->Add(hgemtrkr_2_peak_x);
  hgemtrkr_2_peak_y = new TH1F("hgemtrkr_2_peak_y"," GEM-TRKR2 Peak Y Position ; Y [mm] ",256,-0.2,102.2);  HistList->Add(hgemtrkr_2_peak_y);
  hgemtrkr_2_peak_x_height = new TH1F("hgemtrkr_2_peak_x_height"," GEM-TRKR2 Peak Amplitudes in X; ADC Value ",410,0.,4100.);  HistList->Add(hgemtrkr_2_peak_x_height);
  hgemtrkr_2_peak_y_height = new TH1F("hgemtrkr_2_peak_y_height"," GEM-TRKR2 Peak Amplitudes in Y ; ADC Value ",410,0.,4100.);  HistList->Add(hgemtrkr_2_peak_y_height);
  //--GEMTracker 3
  hgemtrkr_3_peak_xy = new TH2F("hgemtrkr_3_peak_xy","GEM-TRKR3 Peak X-Y Correlation; Peak X [mm]; Peak Y [mm] ",256,-0.2,102.2,256,-0.2,102.2); hgemtrkr_3_peak_xy->SetStats(0); HistList->Add(hgemtrkr_3_peak_xy);
  hgemtrkr_3_max_xy = new TH2F("hgemtrkr_3_max_xy","GEM-TRKR3 X-Y Correlation for Max Hits; Peak X [mm]; Peak Y [mm] ",256,-0.2,102.2,256,-0.2,102.2); hgemtrkr_3_max_xy->SetStats(0);   HistList->Add(hgemtrkr_3_max_xy);
  hgemtrkr_3_peak_x = new TH1F("hgemtrkr_3_peak_x"," GEM-TRKR3 Peak X Position ; X [mm] ",256,-0.2,102.2);  HistList->Add(hgemtrkr_3_peak_x);
  hgemtrkr_3_peak_y = new TH1F("hgemtrkr_3_peak_y"," GEM-TRKR3 Peak Y Position ; Y [mm] ",256,-0.2,102.2);  HistList->Add(hgemtrkr_3_peak_y);
  hgemtrkr_3_peak_x_height = new TH1F("hgemtrkr_3_peak_x_height"," GEM-TRKR3 Peak Amplitudes in X; ADC Value ",410,0.,4100.);  HistList->Add(hgemtrkr_3_peak_x_height);
  hgemtrkr_3_peak_y_height = new TH1F("hgemtrkr_3_peak_y_height"," GEM-TRKR3 Peak Amplitudes in Y ; ADC Value ",410,0.,4100.);  HistList->Add(hgemtrkr_3_peak_y_height);
  
  mmg1_peak_y = new TH1F("mmg1_peak_y"," MMG1-TRD Peak Y Position (SRS) ; Y [mm] ",256,-0.2,102.2);  HistList->Add(mmg1_peak_y);
  hmmg1_peak_y_height = new TH1F("hmmg1_peak_y_height"," MMG1-TRD Peak Amplitudes in Y ; ADC Value ",410,0.,4100.);  HistList->Add(hmmg1_peak_y_height);
  tgem_peak_y = new TH1F("tgem_peak_y","Triple GEM-TRD Peak Y Position (SRS) ; Y [mm] ",256,-0.2,102.2);  HistList->Add(tgem_peak_y);
  htgem_peak_y_height = new TH1F("htgem_peak_y_height","Triple GEM-TRD Peak Amplitudes in Y ; ADC Value ",410,0.,4100.);  HistList->Add(htgem_peak_y_height);
  hgemtrkr_1_tgem = new TH2F("hgemtrkr_1_tgem","GEM-TRKR1 & Triple GEM-TRD Y Correlation; Triple GEM-TRD Y [mm]; GEM-TRKR1 Y [mm] ",256,-0.2,102.2,256,-0.2,102.2); hgemtrkr_1_tgem->SetStats(0); HistList->Add(hgemtrkr_1_tgem);
  hgemtrkr_1_mmg1 = new TH2F("hgemtrkr_1_mmg1","GEM-TRKR1 & MMG1 Y Correlation; MMG1-TRD Y [mm]; GEM-TRKR1 Y [mm] ",256,-0.2,102.2,256,-0.2,102.2); hgemtrkr_1_mmg1->SetStats(0); HistList->Add(hgemtrkr_1_mmg1);
  
  //--External Tracking
  f125_el_tracker_hits = new TH1F("f125_el_tracker_hits","GEM-TRD Track Extr. Hits; X Chan [mm]",256,-0.2,102.2);   HistList->Add(f125_el_tracker_hits);
  f125_el_tracker_eff = new TH1F("f125_el_tracker_eff","GEM-TRD Track Eff.; X Chan [mm]",256,-0.2,102.2);   HistList->Add(f125_el_tracker_eff);
  tgem_residuals = new TH1F("tgem_residuals","Triple GEM-TRD Residual Hits; X Chan [mm] (actual - expected)",250,-25.,25.);   HistList->Add(tgem_residuals);
  tgem_residualscorr = new TH1F("tgem_residualscorr","Triple GEM-TRD Residual Hits WITH CORR.; X Chan [mm] (actual - expected)",250,-25.,25.);   HistList->Add(tgem_residualscorr);
  tgem_residual_ch = new TH2F("tgem_residual_ch","Triple GEM-TRD Residual Hits vs Chan; X Chan [mm] (actual); X Chan [mm] (actual - expected)",175, 11.5, 81.5,35,-7.5,6.5); tgem_residual_ch->SetStats(0); HistList->Add(tgem_residual_ch);
  tgem_residual_chcorr = new TH2F("tgem_residual_chcorr","Triple GEM-TRD Residual Hits vs Chan WITH CORR.; X Chan [mm] (actual); X Chan [mm] (actual - expected)",175, 11.5, 81.5,35,-7.5,6.5); tgem_residual_chcorr->SetStats(0); HistList->Add(tgem_residual_chcorr);
  mmg1_f125_el_tracker_hits = new TH1F("mmg1_f125_el_tracker_hits","MMG1-TRD Track Extr. Hits; X Chan [mm]",256,-0.2,102.2);   HistList->Add(mmg1_f125_el_tracker_hits);
  mmg1_f125_el_tracker_eff = new TH1F("mmg1_f125_el_tracker_eff","MMG1-TRD Track Eff.; X Chan [mm]",256,-0.2,102.2);   HistList->Add(mmg1_f125_el_tracker_eff);
  mmg1_residuals = new TH1F("mmg1_residuals","MMG1-TRD Residual Hits; X Chan [mm] (actual - expected)",250,-25.,25.);   HistList->Add(mmg1_residuals);
  mmg1_residualscorr = new TH1F("mmg1_residualscorr","MMG1-TRD Residual Hits WITH CORR.; X Chan [mm] (actual - expected)",250,-25.,25.);   HistList->Add(mmg1_residualscorr);
  mmg1_residual_ch = new TH2F("mmg1_residual_ch","MMG1-TRD Residual Hits vs Chan; X Chan [mm] (actual); X Chan [mm] (actual - expected)",125, 10.5, 60.5,35,-7.5,6.5); mmg1_residual_ch->SetStats(0); HistList->Add(mmg1_residual_ch);
  mmg1_residual_chcorr = new TH2F("mmg1_residual_chcorr","MMG1-TRD Residual Hits vs Chan WITH CORR.; X Chan [mm] (actual); X Chan [mm] (actual - expected)",125, 10.5, 60.5,35,-7.5,6.5); mmg1_residual_chcorr->SetStats(0); HistList->Add(mmg1_residual_chcorr);
  urw_f125_x_tracker_hits = new TH1F("urw_f125_x_tracker_hits","uRWELL-TRD X Track Extr. Hits; X Chan [mm]",128,-0.4,102.);   HistList->Add(urw_f125_x_tracker_hits);
  urw_f125_x_tracker_eff = new TH1F("urw_f125_x_tracker_eff","uRWELL-TRD X Track Eff.; X Chan [mm]",128,-0.4,102.);   HistList->Add(urw_f125_x_tracker_eff);
  urw_x_residuals = new TH1F("urw_x_residuals","uRWELL-TRD X Residual Hits; X Chan [mm] (actual - expected)",250,-25.,25.);   HistList->Add(urw_x_residuals);
  urw_x_residualscorr = new TH1F("urw_x_residualscorr","uRWELL-TRD X Residual Hits WITH CORR.; X Chan [mm] (actual - expected)",250,-25.,25.);   HistList->Add(urw_x_residualscorr);
  urw_x_residual_ch = new TH2F("urw_x_residual_ch","uRWELL-TRD Residual Hits vs Chan; X Chan [mm] (actual); X Chan [mm] (actual - expected)",95, 12.5, 88.5,15,-4.5,7.5); urw_x_residual_ch->SetStats(0); HistList->Add(urw_x_residual_ch);
  urw_x_residual_chcorr = new TH2F("urw_x_residual_chcorr","uRWELL-TRD Residual Hits vs Chan WITH CORR.; X Chan [mm] (actual); X Chan [mm] (actual - expected)",95, 12.5, 88.5,15,-4.5,7.5); urw_x_residual_chcorr->SetStats(0); HistList->Add(urw_x_residual_chcorr);
  
  //-- Amplitude Histos --
  f125_el_amp2d = new TH2F("f125_el_amp2d","Triple GEM-TRD ADC Amp in Time (After Selections); Time Response (*8ns) ; X Channel ",(int)histTime,0.5,histTime+.5,240,-0.5,239.5); f125_el_amp2d->SetStats(0); HistList->Add(f125_el_amp2d);
  f125_el_amp2d_max = new TH2F("f125_el_amp2d_max","Triple GEM-TRD Max ADC Amp in Time; Time Response (*8ns) ; X [mm] ",(int)histTime,0.5,histTime+.5,256,-0.2,102.2); f125_el_amp2d_max->SetStats(0); HistList->Add(f125_el_amp2d_max);
  f125_xVSamp = new TH2F("f125_xVSamp","Triple GEM-TRD X Channel vs ADC Amp; X Channel; ADC Value",240,-0.5,239.5,410,0.,4100.); f125_xVSamp->SetStats(0); HistList->Add(f125_xVSamp);
  f125_xVSamp_max = new TH2F("f125_xVSamp_max","Triple GEM-TRD X Channel vs Max ADC Amp; X Channel; ADC Value",256,-0.2,102.2,410,0.,4100.); f125_xVSamp_max->SetStats(0); HistList->Add(f125_xVSamp_max);
  f125_timeVSamp = new TH2F("f125_timeVSamp","Triple GEM-TRD X ADC Amp in Time; Time Response (*8ns); ADC Value",(int)histTime,0.5,histTime+.5,410,0.,4100.); f125_timeVSamp->SetStats(0); HistList->Add(f125_timeVSamp);
  f125_timeVSamp_max = new TH2F("f125_timeVSamp_max","Triple GEM-TRD X Max ADC Amp in Time; Time Response (*8ns); ADC Value",(int)histTime,0.5,histTime+.5,410,0.,4100.); f125_timeVSamp_max->SetStats(0); HistList->Add(f125_timeVSamp_max);
  urw_f125_x_amp2d = new TH2F("urw_f125_x_amp2d","uRWELL-TRD X ADC Amp in Time (After Selections); Time Response (*8ns) ; X Channel ",(int)histTime,0.5,histTime+.5,120,-0.5,119.5); urw_f125_x_amp2d->SetStats(0); HistList->Add(urw_f125_x_amp2d);
  urw_f125_x_amp2d_max = new TH2F("urw_f125_x_amp2d_max","uRWELL-TRD X Max ADC Amp in Time; Time Response (*8ns) ; X [mm] ",(int)histTime,0.5,histTime+.5,128,-0.2,102.2); urw_f125_x_amp2d_max->SetStats(0); HistList->Add(urw_f125_x_amp2d_max);
  urw_f125_xVSamp = new TH2F("urw_f125_xVSamp","uRWELL-TRD X Channel vs ADC Amp; X Channel; ADC Value",120,-0.5,119.5,410,0.,4100.); urw_f125_xVSamp->SetStats(0); HistList->Add(urw_f125_xVSamp);
  urw_f125_xVSamp_max = new TH2F("urw_f125_xVSamp_max","uRWELL-TRD X Channel vs Max ADC Amp; X Channel; ADC Value",128,-0.4,102.4,410,0.,4100.); urw_f125_xVSamp_max->SetStats(0); HistList->Add(urw_f125_xVSamp_max);
  urw_f125_x_timeVSamp = new TH2F("urw_f125_x_timeVSamp","uRWELL-TRD X ADC Amp in Time (After Selections); Time Response (*8ns); ADC Value",(int)histTime,0.5,histTime+.5,410,0.,4100.); urw_f125_x_timeVSamp->SetStats(0); HistList->Add(urw_f125_x_timeVSamp);
  urw_f125_x_timeVSamp_max = new TH2F("urw_f125_x_timeVSamp_max","uRWELL-TRD X Max ADC Amp in Time; Time Response (*8ns); ADC Value",(int)histTime,0.5,histTime+.5,410,0.,4100.); urw_f125_x_timeVSamp_max->SetStats(0); HistList->Add(urw_f125_x_timeVSamp_max);
  mmg1_f125_el_amp2d = new TH2F("mmg1_f125_el_amp2d","MMG1-TRD ADC Amp in Time (After Selections); Time Response (*8ns) ; X Channel ",(int)histTime,0.5,histTime+.5,240,-0.5,239.5); mmg1_f125_el_amp2d->SetStats(0); HistList->Add(mmg1_f125_el_amp2d);
  mmg1_f125_el_amp2d_max = new TH2F("mmg1_f125_el_amp2d_max","MMG1-TRD Max ADC Amp in Time; Time Response (*8ns) ; X [mm] ",(int)histTime,0.5,histTime+.5,256,-0.2,102.2); mmg1_f125_el_amp2d_max->SetStats(0); HistList->Add(mmg1_f125_el_amp2d_max);
  mmg1_f125_xVSamp = new TH2F("mmg1_f125_xVSamp","MMG1-TRD X Channel vs ADC Amp; X Channel; ADC Value",240,-0.5,239.5,410,0.,4100.); mmg1_f125_xVSamp->SetStats(0); HistList->Add(mmg1_f125_xVSamp);
  mmg1_f125_xVSamp_max = new TH2F("mmg1_f125_xVSamp_max","MMG1-TRD X Channel vs Max ADC Amp; X Channel; ADC Value",256,-0.2,102.2,410,0.,4100.); mmg1_f125_xVSamp_max->SetStats(0); HistList->Add(mmg1_f125_xVSamp_max);
  mmg1_f125_timeVSamp = new TH2F("mmg1_f125_timeVSamp","MMG1-TRD X ADC Amp in Time; Time Response (*8ns); ADC Value",(int)histTime,0.5,histTime+.5,410,0.,4100.); mmg1_f125_timeVSamp->SetStats(0); HistList->Add(mmg1_f125_timeVSamp);
  mmg1_f125_timeVSamp_max = new TH2F("mmg1_f125_timeVSamp_max","MMG1-TRD X Max ADC Amp in Time; Time Response (*8ns); ADC Value",(int)histTime,0.5,histTime+.5,410,0.,4100.); mmg1_f125_timeVSamp_max->SetStats(0); HistList->Add(mmg1_f125_timeVSamp_max);
  urw_f125_y_amp2d = new TH2F("urw_f125_y_amp2d","uRWELL-TRD Y ADC Amp in Time; Time Response (*8ns) ; Y Channel ",(int)histTime,0.5,histTime+.5,120,-0.5,119.5); urw_f125_y_amp2d->SetStats(0); HistList->Add(urw_f125_y_amp2d);
  urw_f125_y_amp2d_max = new TH2F("urw_f125_y_amp2d_max","uRWELL-TRD Y Max ADC Amp in Time; Time Response (*8ns) ; Y [mm] ",(int)histTime,0.5,histTime+.5,128,-0.2,102.2); urw_f125_y_amp2d_max->SetStats(0); HistList->Add(urw_f125_y_amp2d_max);
  urw_f125_yVSamp = new TH2F("urw_f125_yVSamp","uRWELL-TRD Y Channel vs ADC Amp; Y Channel; ADC Value",120,-0.5,119.5,410,0.,4100.); urw_f125_yVSamp->SetStats(0); HistList->Add(urw_f125_yVSamp);
  
  f125_el_amp2ds = new TH2F("f125_el_amp2ds","Triple GEM-TRD ADC Amp in Time (Before Selections); Time Response (*8ns) ; X Channel ",(int)histTime,0.5,histTime+.5,240,-0.5,239.5); f125_el_amp2ds->SetStats(0); HistList->Add(f125_el_amp2ds);
  urw_f125_x_amp2ds = new TH2F("urw_f125_x_amp2ds","uRWELL-TRD X ADC Amp in Time (Before Selections); Time Response (*8ns) ; X Channel ",(int)histTime,0.5,histTime+.5,120,-0.5,119.5); urw_f125_x_amp2ds->SetStats(0); HistList->Add(urw_f125_x_amp2ds);
  mmg1_f125_el_amp2ds = new TH2F("mmg1_f125_el_amp2ds","MMG1-TRD ADC Amp in Time (Before Selections); Time Response (*8ns) ; X Channel ",(int)histTime,0.5,histTime+.5,240,-0.5,239.5); mmg1_f125_el_amp2ds->SetStats(0); HistList->Add(mmg1_f125_el_amp2ds);
  urw_f125_y_amp2ds = new TH2F("urw_f125_y_amp2ds","uRWELL-TRD Y ADC Amp in Time (Before Selections); Time Response (*8ns) ; Y Channel ",(int)histTime,0.5,histTime+.5,120,-0.5,119.5); urw_f125_y_amp2ds->SetStats(0); HistList->Add(urw_f125_y_amp2ds);

  //-- 2D Hit Displays --
  htgem_xy = new TH2F("htgem_xy","Triple GEM-TRD X-Y Hit Display; X (fADC) [mm]; Y (SRS) [mm] ",256,-0.2,102.2,256,-0.2,102.2); htgem_xy->SetStats(0);  HistList->Add(htgem_xy);
  htgem_max_xy = new TH2F("htgem_max_xy","Triple GEM-TRD X-Y Max Hit Display; X (fADC) [mm]; Y (SRS) [mm] ",256,-0.2,102.2,256,-0.2,102.2); htgem_max_xy->SetStats(0);  HistList->Add(htgem_max_xy);
  hmmg1_xy = new TH2F("hmmg1_xy","MMG1-TRD X-Y Hit Display; X (fADC) [mm]; Y (SRS) [mm] ",256,-0.2,102.2,256,-0.2,102.2); hmmg1_xy->SetStats(0);  HistList->Add(hmmg1_xy);
  hmmg1_max_xy = new TH2F("hmmg1_max_xy","MMG1-TRD X-Y Max Hit Display; X (fADC) [mm]; Y (SRS) [mm] ",256,-0.2,102.2,256,-0.2,102.2); hmmg1_max_xy->SetStats(0);  HistList->Add(hmmg1_max_xy);
  hurw_xy = new TH2F("hurw_xy","uRWELL-TRD X-Y Hit Display; X (fADC) [mm]; Y (fADC) [mm] ",128,-0.4,102.,128,-0.4,102.); hurw_xy->SetStats(0);  HistList->Add(hurw_xy);
  hurw_max_xy = new TH2F("hurw_max_xy","uRWELL-TRD X-Y Max Hit Display; X (fADC) [mm]; Y (fADC) [mm] ",128,-0.4,102.,128,-0.4,102.); hurw_max_xy->SetStats(0);  HistList->Add(hurw_max_xy);

  //-- X,Y Correlations --
  tgem_mmg1_xcorr = new TH2F("tgem_mmg1_xcorr","Triple GEM-TRD X Correlation With MMG1-TRD; Triple GEM-TRD X [mm]; MMG1-TRD X [mm] ",256,-0.2,102.2,256,-0.2,102.2); HistList->Add(tgem_mmg1_xcorr);
  urw_tgem_xcorr = new TH2F("urw_tgem_xcorr","Triple GEM-TRD X Correlation With uRWELL-TRD; uRWELL-TRD X [mm]; Triple GEM-TRD X [mm]",128,-0.4,102.,256,-0.2,102.2); urw_tgem_xcorr->SetStats(0); HistList->Add(urw_tgem_xcorr);
  urw_mmg1_xcorr = new TH2F("urw_mmg1_xcorr","MMG1-TRD X Correlation With uRWELL-TRD; uRWELL-TRD X [mm];  MMG1-TRD X [mm]",128,-0.4,102.,256,-0.2,102.2); urw_mmg1_xcorr->SetStats(0); HistList->Add(urw_mmg1_xcorr);
  tgem_mmg1_ycorr = new TH2F("tgem_mmg1_ycorr","Triple GEM-TRD Y Correlation With MMG1-TRD; Triple GEM-TRD Y [mm]; MMG1-TRD Y [mm] ",256,-0.2,102.2,256,-0.2,102.2); HistList->Add(tgem_mmg1_ycorr);
  tgem_mmg1_max_xcorr = new TH2F("tgem_mmg1_max_xcorr","Triple GEM-TRD Max X Correlation With MMG1-TRD; Triple GEM-TRD X [mm]; MMG1-TRD X [mm] ",256,-0.2,102.2,256,-0.2,102.2); HistList->Add(tgem_mmg1_max_xcorr);
  tgem_urw_max_xcorr = new TH2F("tgem_urw_max_xcorr","Triple GEM-TRD Max X Correlation With uRWELL-TRD; uRWELL-TRD X [mm]; Triple GEM-TRD X [mm]",128,-0.4,102.,256,-0.2,102.2); HistList->Add(tgem_urw_max_xcorr);
  urw_mmg1_max_xcorr = new TH2F("urw_mmg1_max_xcorr","uRWELL-TRD Max X Correlation With MMG1-TRD; uRWELL-TRD X [mm]; MMG1-TRD X [mm] ",128,-0.4,102.,256,-0.2,102.2); HistList->Add(urw_mmg1_max_xcorr);
  tgem_gt1_xcorr = new TH2F("tgem_gt1_xcorr","Triple GEM-TRD X Correlation With GEMTRKR-1; Triple GEM-TRD X [mm]; GEMTRKR-1 X [mm] ",256,-0.2,102.2,256,-0.2,102.2); HistList->Add(tgem_gt1_xcorr);
  tgem_gt2_xcorr = new TH2F("tgem_gt2_xcorr","Triple GEM-TRD X Correlation With GEMTRKR-2; Triple GEM-TRD X [mm]; GEMTRKR-2 X [mm] ",256,-0.2,102.2,256,-0.2,102.2); HistList->Add(tgem_gt2_xcorr);
  tgem_gt3_xcorr = new TH2F("tgem_gt3_xcorr","Triple GEM-TRD X Correlation With GEMTRKR-3; Triple GEM-TRD X [mm]; GEMTRKR-3 X [mm] ",256,-0.2,102.2,256,-0.2,102.2); HistList->Add(tgem_gt3_xcorr);
  mmg1_gt1_xcorr = new TH2F("mmg1_gt1_xcorr","MMG1-TRD X Correlation With GEMTRKR-1; MMG1-TRD X [mm]; GEMTRKR-1 X [mm] ",256,-0.2,102.2,256,-0.2,102.2); HistList->Add(mmg1_gt1_xcorr);
  mmg1_gt2_xcorr = new TH2F("mmg1_gt2_xcorr","MMG1-TRD X Correlation With GEMTRKR-2; MMG1-TRD X [mm]; GEMTRKR-2 X [mm] ",256,-0.2,102.2,256,-0.2,102.2); HistList->Add(mmg1_gt2_xcorr);
  mmg1_gt3_xcorr = new TH2F("mmg1_gt3_xcorr","MMG1-TRD X Correlation With GEMTRKR-3; MMG1-TRD X [mm]; GEMTRKR-3 X [mm] ",256,-0.2,102.2,256,-0.2,102.2); HistList->Add(mmg1_gt3_xcorr);
  urw_gt1_xcorr = new TH2F("urw_gt1_xcorr","uRWell-TRD X Correlation With GEMTRKR-1; uRWell-TRD X [mm]; GEMTRKR-1 X [mm] ",128,-0.4,102.,256,-0.2,102.2); HistList->Add(urw_gt1_xcorr);
  urw_gt2_xcorr = new TH2F("urw_gt2_xcorr","uRWell-TRD X Correlation With GEMTRKR-2; uRWell-TRD X [mm]; GEMTRKR-2 X [mm] ",128,-0.4,102.,256,-0.2,102.2); HistList->Add(urw_gt2_xcorr);
  urw_gt3_xcorr = new TH2F("urw_gt3_xcorr","uRWell-TRD X Correlation With GEMTRKR-3; uRWell-TRD X [mm]; GEMTRKR-3 X [mm] ",128,-0.4,102.,256,-0.2,102.2); HistList->Add(urw_gt3_xcorr);
  
  hgemClusterDiff_el = new TH1F("hgemClusterDiff_el","GEM Cluster Distance from MMG1 Track; Distance [mm]",200,-20.,20.); HistList->Add(hgemClusterDiff_el);
  hmmg1ClusterDiff_el = new TH1F("hmmg1ClusterDiff_el","MMG1 Cluster Distance from GEM Track; Distance [mm]",200,-20.,20.); HistList->Add(hmmg1ClusterDiff_el);
  hgemPulseDiff_el = new TH1F("hgemPulseDiff_el","GEM Pulse Distance from MMG1 Track; Distance [mm]",150,-30.,30.); HistList->Add(hgemPulseDiff_el);
  hmmg1PulseDiff_el = new TH1F("hmmg1PulseDiff_el","MMG1 Pulse Distance from GEM Track; Distance [mm]",150,-30.,30.); HistList->Add(hmmg1PulseDiff_el);
  hurwPulseDiff_el = new TH1F("hurwPulseDiff_el","uRWELL Pulse Distance from GEM Track; Distance [mm]",75,-30.,30.); HistList->Add(hurwPulseDiff_el);
  hurwPulseDiff_mmg = new TH1F("hurwPulseDiff_mmg","uRWELL Pulse Distance from MMG Track; Distance [mm]",75,-30.,30.); HistList->Add(hurwPulseDiff_mmg);
  hTrackDiff = new TH1F("hTrackDiff","Tracker 1 & 3 X Pos Difference; Distance [mm]",275,-100.,10.); HistList->Add(hTrackDiff);
  
  hClusterMaxdEdx_el = new TH1F("hClusterMaxdEdx_el","GEM Max Cluster Energy; Cluster Energy ; Counts ",100,0.,4960); HistList->Add(hClusterMaxdEdx_el);
  hClusterTotaldEdx_el = new TH1F("hClusterTotaldEdx_el","GEM Total Cluster Energy; Total Cluster Energy ; Counts ",100,0.,4960); HistList->Add(hClusterTotaldEdx_el);
  hmmg1ClusterMaxdEdx_el = new TH1F("hmmg1ClusterMaxdEdx_el","MMG1 Max Cluster Energy; Cluster Energy ; Counts ",100,0.,4960); HistList->Add(hmmg1ClusterMaxdEdx_el);
  hmmg1ClusterTotaldEdx_el = new TH1F("hmmg1ClusterTotaldEdx_el","MMG1 Total Cluster Energy; Total Cluster Energy ; Counts ",100,0.,4960); HistList->Add(hmmg1ClusterTotaldEdx_el);
  
  //=========================================
  
  TFile* fHits;
  #ifdef SAVE_TRACK_HITS
    #if ANALYZE_MERGED
    char hitsFileName[256]; sprintf(hitsFileName, "RootOutput/ps25/merged/trd_singleTrackHits_Run_%06d_%06dEntries.root",RunNum,nEntries);
    #else
    char hitsFileName[256]; sprintf(hitsFileName, "RootOutput/ps25/trd_singleTrackHits_Run_%06d.root", RunNum);
    #endif
    fHits = new TFile(hitsFileName, "RECREATE");
    tgem_zHist = new TH1F("tgem_zHist", "tgem_zHist", 20, 80., histTime);
    mmg1_zHist = new TH1F("mmg1_zHist", "mmg1_zHist", 20, 80., histTime);
    urw_zHist = new TH1F("urw_zHist", "urw_zHist", 20, 80., histTime);
    //-- Triple GEM-TRD
    EVENT_VECT_GEM = new TTree("gem_hits","Triple-GEM TTree with single track hit info");
    EVENT_VECT_GEM->Branch("event_num",&event_num,"event_num/I");
    EVENT_VECT_GEM->Branch("nhit",&tgem_nhit,"tgem_nhit/I");
    EVENT_VECT_GEM->Branch("nclu",&tgem_nclu,"tgem_nclu/I");
    EVENT_VECT_GEM->Branch("ntracks",&tgem_ntracks,"tgem_ntracks/I");
    EVENT_VECT_GEM->Branch("xpos",&tgem_xpos);
    EVENT_VECT_GEM->Branch("zpos",&tgem_zpos);
    EVENT_VECT_GEM->Branch("dedx",&tgem_dedx);
    EVENT_VECT_GEM->Branch("parID",&tgem_parID);
    EVENT_VECT_GEM->Branch("zHist",&tgem_zHist_vect);
    EVENT_VECT_GEM->Branch("xposc",&clu_xpos);
    EVENT_VECT_GEM->Branch("zposc",&clu_zpos);
    EVENT_VECT_GEM->Branch("dedxc",&clu_dedx);
    EVENT_VECT_GEM->Branch("widthc",&clu_width);
    EVENT_VECT_GEM->Branch("xposc_max",&clu_xpos_max,"clu_xpos_max/f");
    EVENT_VECT_GEM->Branch("zposc_max",&clu_zpos_max,"clu_zpos_max/f");
    EVENT_VECT_GEM->Branch("dedxc_max",&clu_dedx_max,"clu_dedx_max/f");
    EVENT_VECT_GEM->Branch("widthc_max",&clu_width_max,"clu_width_max/f");
    EVENT_VECT_GEM->Branch("dedxc_tot",&clu_dedx_tot,"clu_dedx_tot/f");
    EVENT_VECT_GEM->Branch("xch_max",&tgem_xch_max,"tgem_xch_max/f");
    EVENT_VECT_GEM->Branch("amp_max",&tgem_amp_max,"tgem_amp_max/f");
    EVENT_VECT_GEM->Branch("time_max",&tgem_time_max,"tgem_time_max/f");
    EVENT_VECT_GEM->Branch("chi2",&tgem_chi2cc);
    EVENT_VECT_GEM->Branch("Fint",&tgem_integral);
    EVENT_VECT_GEM->Branch("a0",&a0);
    EVENT_VECT_GEM->Branch("a1",&a1);
    
    //-- MMG1-TRD
    EVENT_VECT_MMG1 = new TTree("mmg1_hits","MMG1 TTree with single track hit info");
    EVENT_VECT_MMG1->Branch("event_num",&event_num,"event_num/I");
    EVENT_VECT_MMG1->Branch("nhit",&mmg1_nhit,"mmg1_nhit/I");
    EVENT_VECT_MMG1->Branch("nclu",&mmg1_nclu,"mmg1_nclu/I");
    EVENT_VECT_MMG1->Branch("ntracks",&mmg1_ntracks,"mmg1_ntracks/I");
    EVENT_VECT_MMG1->Branch("xpos",&mmg1_xpos);
    EVENT_VECT_MMG1->Branch("zpos",&mmg1_zpos);
    EVENT_VECT_MMG1->Branch("dedx",&mmg1_dedx);
    EVENT_VECT_MMG1->Branch("parID",&mmg1_parID);
    EVENT_VECT_MMG1->Branch("zHist",&mmg1_zHist_vect);
    EVENT_VECT_MMG1->Branch("xposc",&mmg1_clu_xpos);
    EVENT_VECT_MMG1->Branch("zposc",&mmg1_clu_zpos);
    EVENT_VECT_MMG1->Branch("dedxc",&mmg1_clu_dedx);
    EVENT_VECT_MMG1->Branch("widthc",&mmg1_clu_width);
    EVENT_VECT_MMG1->Branch("xposc_max",&mmg1_clu_xpos_max,"mmg1_clu_xpos_max/f");
    EVENT_VECT_MMG1->Branch("zposc_max",&mmg1_clu_zpos_max,"mmg1_clu_zpos_max/f");
    EVENT_VECT_MMG1->Branch("dedxc_max",&mmg1_clu_dedx_max,"mmg1_clu_dedx_max/f");
    EVENT_VECT_MMG1->Branch("widthc_max",&mmg1_clu_width_max,"mmg1_clu_width_max/f");
    EVENT_VECT_MMG1->Branch("dedxc_tot",&mmg1_clu_dedx_tot,"mmg1_clu_dedx_tot/f");
    EVENT_VECT_MMG1->Branch("xch_max",&mmg1_xch_max,"mmg1_xch_max/f");
    EVENT_VECT_MMG1->Branch("amp_max",&mmg1_amp_max,"mmg1_amp_max/f");
    EVENT_VECT_MMG1->Branch("time_max",&mmg1_time_max,"mmg1_time_max/f");
    EVENT_VECT_MMG1->Branch("chi2",&mmg1_chi2cc);
    EVENT_VECT_MMG1->Branch("Fint",&mmg1_integral);
    EVENT_VECT_MMG1->Branch("a0",&mmg1_a0);
    EVENT_VECT_MMG1->Branch("a1",&mmg1_a1);
    
     //-- uRWELL-TRD
    EVENT_VECT_URW = new TTree("urw_hits","uRWELL TTree with single track hit info");
    EVENT_VECT_URW->Branch("event_num",&event_num,"event_num/I");
    EVENT_VECT_URW->Branch("nxhit",&urw_nxhit,"urw_nxhit/I");
    EVENT_VECT_URW->Branch("nyhit",&urw_nyhit,"urw_nyhit/I");
    EVENT_VECT_URW->Branch("nclu",&urw_nclu,"urw_nclu/I");
    EVENT_VECT_URW->Branch("ntracks",&urw_ntracks,"urw_ntracks/I");
    EVENT_VECT_URW->Branch("xpos",&urw_xpos);
    EVENT_VECT_URW->Branch("zpos",&urw_zpos);
    EVENT_VECT_URW->Branch("dedx",&urw_dedx);
    EVENT_VECT_URW->Branch("parID",&urw_parID);
    EVENT_VECT_URW->Branch("zHist",&urw_zHist_vect);
    EVENT_VECT_URW->Branch("xposc",&urw_clu_xpos);
    EVENT_VECT_URW->Branch("zposc",&urw_clu_zpos);
    EVENT_VECT_URW->Branch("dedxc",&urw_clu_dedx);
    EVENT_VECT_URW->Branch("widthc",&urw_clu_width);
    EVENT_VECT_URW->Branch("xposc_max",&urw_clu_xpos_max,"urw_clu_xpos_max/f");
    EVENT_VECT_URW->Branch("zposc_max",&urw_clu_zpos_max,"urw_clu_zpos_max/f");
    EVENT_VECT_URW->Branch("dedxc_max",&urw_clu_dedx_max,"urw_clu_dedx_max/f");
    EVENT_VECT_URW->Branch("widthc_max",&urw_clu_width_max,"urw_clu_width_max/f");
    EVENT_VECT_URW->Branch("dedxc_tot",&urw_clu_dedx_tot,"urw_clu_dedx_tot/f");
    EVENT_VECT_URW->Branch("xch_max",&urw_xch_max,"urw_xch_max/f");
    EVENT_VECT_URW->Branch("amp_max",&urw_amp_max,"urw_amp_max/f");
    EVENT_VECT_URW->Branch("time_max",&urw_time_max,"urw_time_max/f");
    EVENT_VECT_URW->Branch("chi2",&urw_chi2cc);
    EVENT_VECT_URW->Branch("Fint",&urw_integral);
    EVENT_VECT_URW->Branch("a0",&urw_a0);
    EVENT_VECT_URW->Branch("a1",&urw_a1);
    
  #endif
  
  TStopwatch timer;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes=0, nb=0;
  if (MaxEvt>0) nentries=MaxEvt;  //-- limit number of events for test
  
  //==================================================================================================
  //                      E v e n t    L o o p
  //==================================================================================================
  printf("===============  Begin Event Loop - 1st evt=%lld, Last evt=%lld =============== \n",FirstEvt,MaxEvt);
  timer.Start();
  Long64_t jentry=0;
  //int el_count=0, pi_count=0, atlas_trigger_count=0;
  
  for (jentry=FirstEvt; jentry<nentries; jentry++) { //-- Event Loop --
    Count("EVT");
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    if (!(jentry%NPRT))
      printf("------- evt=%llu  f125_raw_count=%llu f125_pulse_count=%llu f250_wraw_count=%llu, srs_peak_count=%llu \n",jentry,f125_wraw_count, f125_pulse_count, f250_wraw_count, gem_peak_count);
    event_num=jentry;
    
    bool USE_MAX_TRD_HITS = false;
    bool USE_MAX_TRACKER_HITS = false;
    bool USE_TRD_XCORR = false;
    bool match=false, match_mmg1=false, match_urw=false;
    
    tgem_nhit=0;
    mmg1_nhit=0;
    urw_nxhit=0;
    urw_nyhit=0;
    tgem_nclu=0;
    mmg1_nclu=0;
    
    //-- Triple GEM-TRD
    tgem_xpos.clear();
    tgem_zpos.clear();
    tgem_dedx.clear();
    tgem_parID.clear();
    tgem_zHist->Reset();
    tgem_zHist_vect.clear();
    clu_xpos.clear();
    clu_zpos.clear();
    clu_dedx.clear();
    clu_width.clear();
    clu_xpos_max=0;
    clu_zpos_max=0;
    clu_dedx_max=0;
    clu_width_max=0;
    clu_dedx_tot=0;
    tgem_amp_max=0;
    tgem_xch_max=0;
    tgem_time_max=0;
    tgem_chi2cc.clear();
    tgem_integral.clear();
    
    //-- MMG1-TRD
    mmg1_xpos.clear();
    mmg1_zpos.clear();
    mmg1_dedx.clear();
    mmg1_parID.clear();
    mmg1_zHist->Reset();
    mmg1_zHist_vect.clear();
    mmg1_clu_xpos.clear();
    mmg1_clu_zpos.clear();
    mmg1_clu_dedx.clear();
    mmg1_clu_width.clear();
    mmg1_clu_xpos_max=0;
    mmg1_clu_zpos_max=0;
    mmg1_clu_dedx_max=0;
    mmg1_clu_width_max=0;
    mmg1_clu_dedx_tot=0;
    mmg1_amp_max=0;
    mmg1_xch_max=0;
    mmg1_time_max=0;
    mmg1_chi2cc.clear();
    mmg1_integral.clear();
    
    //-- uRWELL-TRD
    urw_xpos.clear();
    urw_zpos.clear();
    urw_dedx.clear();
    urw_parID.clear();
    urw_zHist->Reset();
    urw_zHist_vect.clear();
    urw_clu_xpos.clear();
    urw_clu_zpos.clear();
    urw_clu_dedx.clear();
    urw_clu_width.clear();
    urw_clu_xpos_max=0;
    urw_clu_zpos_max=0;
    urw_clu_dedx_max=0;
    urw_clu_width_max=0;
    urw_clu_dedx_tot=0;
    urw_amp_max=0;
    urw_xch_max=0;
    urw_time_max=0;
    urw_chi2cc.clear();
    urw_integral.clear();
    
    //==================================================================================================
    //                    Process SRS data
    //==================================================================================================
    
    //===========================================================
    //  GEMTracker (SRS) Correlations with TRD Prototypes
    //===========================================================
    
    ULong64_t gt_1_idx_x = 0, gt_1_idx_y=0, gt_2_idx_x = 0, gt_2_idx_y = 0, gt_3_idx_x = 0, gt_3_idx_y = 0, mmg1_idx_y=0, tgem_idx_y=0;
    int gt1_nhit=0, gt2_nhit=0, gt3_nhit=0;
    
    double gemtrkr1_xamp_max=-1., gemtrkr1_xch_max=-1000.;
    double gemtrkr1_yamp_max=-1., gemtrkr1_ych_max=-1000.;
    double gemtrkr2_xamp_max=-1., gemtrkr2_xch_max=-1000.;
    double gemtrkr2_yamp_max=-1., gemtrkr2_ych_max=-1000.;
    double gemtrkr3_xamp_max=-1., gemtrkr3_xch_max=-1000.;
    double gemtrkr3_yamp_max=-1., gemtrkr3_ych_max=-1000.;
    double tgem_yamp_max=-1., tgem_ych_max=-1000.;
    double mmg1_yamp_max=-1., mmg1_ych_max=-1000.;
    
    double gemtrkr_1_peak_pos_y[gem_peak_count];
    double gemtrkr_1_peak_pos_x[gem_peak_count];
    double gemtrkr_1_peak_x_height[gem_peak_count];
    double gemtrkr_1_peak_y_height[gem_peak_count];
    double gemtrkr_2_peak_pos_y[gem_peak_count];
    double gemtrkr_2_peak_pos_x[gem_peak_count];
    double gemtrkr_2_peak_x_height[gem_peak_count];
    double gemtrkr_2_peak_y_height[gem_peak_count];
    double gemtrkr_3_peak_pos_y[gem_peak_count];
    double gemtrkr_3_peak_pos_x[gem_peak_count];
    double gemtrkr_3_peak_x_height[gem_peak_count];
    double gemtrkr_3_peak_y_height[gem_peak_count];
    double mmg1_peak_pos_y[gem_peak_count];
    double mmg1_peak_y_height[gem_peak_count];
    double tgem_peak_pos_y[gem_peak_count];
    double tgem_peak_y_height[gem_peak_count];
    
    for (ULong64_t i=0; i<gem_peak_count; i++) {
      gemtrkr_1_peak_pos_y[gem_peak_count] = -1000;
      gemtrkr_1_peak_pos_x[gem_peak_count] = -1000;
      gemtrkr_1_peak_x_height[gem_peak_count] = -1000;
      gemtrkr_1_peak_y_height[gem_peak_count] = -1000;
      gemtrkr_2_peak_pos_y[gem_peak_count] = -1000;
      gemtrkr_2_peak_pos_x[gem_peak_count] = -1000;
      gemtrkr_2_peak_x_height[gem_peak_count] = -1000;
      gemtrkr_2_peak_y_height[gem_peak_count] = -1000;
      gemtrkr_3_peak_pos_y[gem_peak_count] = -1000;
      gemtrkr_3_peak_pos_x[gem_peak_count] = -1000;
      gemtrkr_3_peak_x_height[gem_peak_count] = -1000;
      gemtrkr_3_peak_y_height[gem_peak_count] = -1000;
      mmg1_peak_pos_y[gem_peak_count] = -1000;
      mmg1_peak_y_height[gem_peak_count] = -1000;
      tgem_peak_pos_y[gem_peak_count] = -1000;
      tgem_peak_y_height[gem_peak_count] = -1000;
    }
   
    for (ULong64_t i=0; i<gem_peak_count; i++) { //-- SRS Peaks Loop
      
      if (gem_peak_plane_name->at(i) == "GEMTR1X") {
        gemtrkr_1_peak_x_height[gt_1_idx_x] = gem_peak_height->at(i);
        if (gemtrkr_1_peak_x_height[gt_1_idx_x]>TRKR_THRESH) {
          gemtrkr_1_peak_pos_x[gt_1_idx_x] = gem_peak_real_pos->at(i);
          if (gemtrkr_1_peak_pos_x[gt_1_idx_x]<0) gemtrkr_1_peak_pos_x[gt_1_idx_x]+=51.; else if (gemtrkr_1_peak_pos_x[gt_1_idx_x]>0) gemtrkr_1_peak_pos_x[gt_1_idx_x]-=51.; gemtrkr_1_peak_pos_x[gt_1_idx_x]*=-1.; gemtrkr_1_peak_pos_x[gt_1_idx_x]+=51.;
          gt_1_idx_x++; Count("gt1_x");
        }
        //gt_1_idx_x++; Count("gt1_x");
      } if (gem_peak_plane_name->at(i) == "GEMTR1Y") {
          gemtrkr_1_peak_y_height[gt_1_idx_y] = gem_peak_height->at(i);
          if (gemtrkr_1_peak_y_height[gt_1_idx_y]>TRKR_THRESH) {
            gemtrkr_1_peak_pos_y[gt_1_idx_y] = gem_peak_real_pos->at(i);
            if (gemtrkr_1_peak_pos_y[gt_1_idx_y]<0) gemtrkr_1_peak_pos_y[gt_1_idx_y]+=51.; else if (gemtrkr_1_peak_pos_y[gt_1_idx_y]>0) gemtrkr_1_peak_pos_y[gt_1_idx_y]-=51.;  gemtrkr_1_peak_pos_y[gt_1_idx_y]*=-1.;  gemtrkr_1_peak_pos_y[gt_1_idx_y]+=51.;
            gt_1_idx_y++; Count("gt1_y");
          }
          //gt_1_idx_y++; Count("gt1_y");
      } if (gem_peak_plane_name->at(i) == "GEMTR2X") {
        gemtrkr_2_peak_x_height[gt_2_idx_x] = gem_peak_height->at(i);
        if (gemtrkr_2_peak_x_height[gt_2_idx_x]>TRKR_THRESH) {
          gemtrkr_2_peak_pos_x[gt_2_idx_x] = gem_peak_real_pos->at(i);
          if (gemtrkr_2_peak_pos_x[gt_2_idx_x]<0) gemtrkr_2_peak_pos_x[gt_2_idx_x]+=51.; else if (gemtrkr_2_peak_pos_x[gt_2_idx_x]>0) gemtrkr_2_peak_pos_x[gt_2_idx_x]-=51.;  gemtrkr_2_peak_pos_x[gt_2_idx_x]*=-1.;  gemtrkr_2_peak_pos_x[gt_2_idx_x]+=51.;
          gt_2_idx_x++; Count("gt2_x");
        }
        //gt_2_idx_x++; Count("gt2_x");
      } if (gem_peak_plane_name->at(i) == "GEMTR2Y") {
          gemtrkr_2_peak_y_height[gt_2_idx_y] = gem_peak_height->at(i);
          if (gemtrkr_2_peak_y_height[gt_2_idx_y]>TRKR_THRESH) {
            gemtrkr_2_peak_pos_y[gt_2_idx_y] = gem_peak_real_pos->at(i);
            if (gemtrkr_2_peak_pos_y[gt_2_idx_y]<0) gemtrkr_2_peak_pos_y[gt_2_idx_y]+=51.; else if (gemtrkr_2_peak_pos_y[gt_2_idx_y]>0) gemtrkr_2_peak_pos_y[gt_2_idx_y]-=51.;  gemtrkr_2_peak_pos_y[gt_2_idx_y]*=-1.;  gemtrkr_2_peak_pos_y[gt_2_idx_y]+=51.;
            gt_2_idx_y++; Count("gt2_y");
        }  //gt_2_idx_y++; Count("gt2_y");
      } if (gem_peak_plane_name->at(i) == "GEMTR3X") {
        gemtrkr_3_peak_x_height[gt_3_idx_x] = gem_peak_height->at(i);
        if (gemtrkr_3_peak_x_height[gt_3_idx_x]>TRKR_THRESH) {
          gemtrkr_3_peak_pos_x[gt_3_idx_x] = gem_peak_real_pos->at(i);
          if (gemtrkr_3_peak_pos_x[gt_3_idx_x]<0) gemtrkr_3_peak_pos_x[gt_3_idx_x]+=51.; else if (gemtrkr_3_peak_pos_x[gt_3_idx_x]>0) gemtrkr_3_peak_pos_x[gt_3_idx_x]-=51.;  gemtrkr_3_peak_pos_x[gt_3_idx_x]*=-1.;  gemtrkr_3_peak_pos_x[gt_3_idx_x]+=51.;
          gt_3_idx_x++; Count("gt3_x");
        }
        //gt_3_idx_x++; Count("gt3_x");
      } if (gem_peak_plane_name->at(i) == "GEMTR3Y") {
          gemtrkr_3_peak_y_height[gt_3_idx_y] = gem_peak_height->at(i);
          if (gemtrkr_3_peak_y_height[gt_3_idx_y]>TRKR_THRESH) {
            gemtrkr_3_peak_pos_y[gt_3_idx_y] = gem_peak_real_pos->at(i);
            if (gemtrkr_3_peak_pos_y[gt_3_idx_y]<0) gemtrkr_3_peak_pos_y[gt_3_idx_y]+=51.; else if (gemtrkr_3_peak_pos_y[gt_3_idx_y]>0) gemtrkr_3_peak_pos_y[gt_3_idx_y]-=51.;  gemtrkr_3_peak_pos_y[gt_3_idx_y]*=-1.;  gemtrkr_3_peak_pos_y[gt_3_idx_y]+=51.;
            gt_3_idx_y++; Count("gt3_y");
          }
          //gt_3_idx_y++; Count("gt3_y");
      } if (gem_peak_plane_name->at(i) == "MMG1TRDY" && RunNum>runSwitch) {
        mmg1_peak_y_height[mmg1_idx_y] = gem_peak_height->at(i);
        if (mmg1_peak_y_height[mmg1_idx_y]>1300.) {
          mmg1_peak_pos_y[mmg1_idx_y] = gem_peak_real_pos->at(i);
          if (mmg1_peak_pos_y[mmg1_idx_y]<0) mmg1_peak_pos_y[mmg1_idx_y]+=51.; else if (mmg1_peak_pos_y[mmg1_idx_y]>0) mmg1_peak_pos_y[mmg1_idx_y]-=51.;  mmg1_peak_pos_y[mmg1_idx_y]*=-1.;  mmg1_peak_pos_y[mmg1_idx_y]+=51.;
          mmg1_idx_y++; Count("mmg1_y");
        }
        //mmg1_idx_y++; Count("mmg1_y");
      } if (gem_peak_plane_name->at(i) == "VU_GEMTRDY") {
        tgem_peak_y_height[tgem_idx_y] = gem_peak_height->at(i);
        if (tgem_peak_y_height[tgem_idx_y]>1300.) {
          tgem_peak_pos_y[tgem_idx_y] = gem_peak_real_pos->at(i);
          if (tgem_peak_pos_y[tgem_idx_y]<0) tgem_peak_pos_y[tgem_idx_y]+=51.; else if (tgem_peak_pos_y[tgem_idx_y]>0) tgem_peak_pos_y[tgem_idx_y]-=51.;  tgem_peak_pos_y[tgem_idx_y]*=-1.;  tgem_peak_pos_y[tgem_idx_y]+=51.;
          tgem_idx_y++; Count("tgem_y");
        }
        //gem_idx_y++; Count("tgem_y");
      }
    }
    
    for (ULong64_t j=0; j<gt_1_idx_y; j++) {
      if (gemtrkr_1_peak_pos_y[j]>-1.) hgemtrkr_1_peak_y->Fill(gemtrkr_1_peak_pos_y[j]);
      hgemtrkr_1_peak_y_height->Fill(gemtrkr_1_peak_y_height[j]);
      if (gemtrkr_1_peak_y_height[j]>gemtrkr1_yamp_max) {
        gemtrkr1_yamp_max=gemtrkr_1_peak_y_height[j];
        gemtrkr1_ych_max=gemtrkr_1_peak_pos_y[j];
      }
      for (ULong64_t k=0; k<gt_1_idx_x; k++) {
        hgemtrkr_1_peak_xy->Fill(gemtrkr_1_peak_pos_x[k], gemtrkr_1_peak_pos_y[j]);
        gt1_nhit++;
      }
    }
    for (ULong64_t j=0; j<gt_2_idx_y; j++) {
      if (gemtrkr_2_peak_pos_y[j]>-1.) hgemtrkr_2_peak_y->Fill(gemtrkr_2_peak_pos_y[j]);
      hgemtrkr_2_peak_y_height->Fill(gemtrkr_2_peak_y_height[j]);
      if (gemtrkr_2_peak_y_height[j]>gemtrkr2_yamp_max) {
        gemtrkr2_yamp_max=gemtrkr_2_peak_y_height[j];
        gemtrkr2_ych_max=gemtrkr_2_peak_pos_y[j];
      }
      for (ULong64_t k=0; k<gt_2_idx_x; k++) {
        hgemtrkr_2_peak_xy->Fill(gemtrkr_2_peak_pos_x[k], gemtrkr_2_peak_pos_y[j]);
        gt2_nhit++;
      }
    }
    for (ULong64_t j=0; j<gt_3_idx_y; j++) {
      if (gemtrkr_3_peak_pos_y[j]>-1.) hgemtrkr_3_peak_y->Fill(gemtrkr_3_peak_pos_y[j]);
      hgemtrkr_3_peak_y_height->Fill(gemtrkr_3_peak_y_height[j]);
      if (gemtrkr_3_peak_y_height[j]>gemtrkr3_yamp_max) {
        gemtrkr3_yamp_max=gemtrkr_3_peak_y_height[j];
        gemtrkr3_ych_max=gemtrkr_3_peak_pos_y[j];
      }
      for (ULong64_t k=0; k<gt_3_idx_x; k++) {
        hgemtrkr_3_peak_xy->Fill(gemtrkr_3_peak_pos_x[k], gemtrkr_3_peak_pos_y[j]);
        gt3_nhit++;
      }
    }
    for (ULong64_t k=0; k<gt_1_idx_x; k++) {
      if (gemtrkr_1_peak_pos_x[k]>-1.) hgemtrkr_1_peak_x->Fill(gemtrkr_1_peak_pos_x[k]);
      hgemtrkr_1_peak_x_height->Fill(gemtrkr_1_peak_x_height[k]);
      if (gemtrkr_1_peak_x_height[k]>gemtrkr1_xamp_max) {
        gemtrkr1_xamp_max=gemtrkr_1_peak_x_height[k];
        gemtrkr1_xch_max=gemtrkr_1_peak_pos_x[k];
      }
    }
    for (ULong64_t k=0; k<gt_2_idx_x; k++) {
      if (gemtrkr_2_peak_pos_x[k]>-1.) hgemtrkr_2_peak_x->Fill(gemtrkr_2_peak_pos_x[k]);
      hgemtrkr_2_peak_x_height->Fill(gemtrkr_2_peak_x_height[k]);
      if (gemtrkr_2_peak_x_height[k]>gemtrkr2_xamp_max) {
        gemtrkr2_xamp_max=gemtrkr_2_peak_x_height[k];
        gemtrkr2_xch_max=gemtrkr_2_peak_pos_x[k];
      }
    }
    for (ULong64_t k=0; k<gt_3_idx_x; k++) {
      if (gemtrkr_3_peak_pos_x[k]>-1.) hgemtrkr_3_peak_x->Fill(gemtrkr_3_peak_pos_x[k]);
      hgemtrkr_3_peak_x_height->Fill(gemtrkr_3_peak_x_height[k]);
      if (gemtrkr_3_peak_x_height[k]>gemtrkr3_xamp_max) {
        gemtrkr3_xamp_max=gemtrkr_3_peak_x_height[k];
        gemtrkr3_xch_max=gemtrkr_3_peak_pos_x[k];
      }
    }
    for (ULong64_t k=0; k<mmg1_idx_y; k++) {
      if (mmg1_peak_pos_y[k]>-1.) mmg1_peak_y->Fill(mmg1_peak_pos_y[k]);
      hmmg1_peak_y_height->Fill(mmg1_peak_y_height[k]);
      if (mmg1_peak_y_height[k]>mmg1_yamp_max) {
        mmg1_yamp_max=mmg1_peak_y_height[k];
        mmg1_ych_max=mmg1_peak_pos_y[k];
      }
    }
    for (ULong64_t k=0; k<tgem_idx_y; k++) {
      if (tgem_peak_pos_y[k]>-1.) tgem_peak_y->Fill(tgem_peak_pos_y[k]);
      htgem_peak_y_height->Fill(tgem_peak_y_height[k]);
      if (tgem_peak_y_height[k]>tgem_yamp_max) {
        tgem_yamp_max=tgem_peak_y_height[k];
        tgem_ych_max=tgem_peak_pos_y[k];
      }
    }  
    
    if (tgem_ych_max>0. && gemtrkr1_ych_max>0.) hgemtrkr_1_tgem->Fill(tgem_ych_max, gemtrkr1_ych_max);
    if (mmg1_ych_max>0. && gemtrkr1_ych_max>0.) hgemtrkr_1_mmg1->Fill(mmg1_ych_max, gemtrkr1_ych_max);
    if (gemtrkr1_xch_max>0.) hgemtrkr_1_max_xch->Fill(gemtrkr1_xch_max);
    if (gemtrkr1_xamp_max>0.) hgemtrkr_1_max_xamp->Fill(gemtrkr1_xamp_max);
    if (gemtrkr2_xch_max>0.) hgemtrkr_2_max_xch->Fill(gemtrkr2_xch_max);
    if (gemtrkr2_xamp_max>0.) hgemtrkr_2_max_xamp->Fill(gemtrkr2_xamp_max);
    if (gemtrkr3_xch_max>0.) hgemtrkr_3_max_xch->Fill(gemtrkr3_xch_max);
    if (gemtrkr3_xamp_max>0.) hgemtrkr_3_max_xamp->Fill(gemtrkr3_xamp_max);
    if (gemtrkr1_xch_max>0. && gemtrkr1_ych_max>0.) hgemtrkr_1_max_xy->Fill(gemtrkr1_xch_max, gemtrkr1_ych_max);
    if (gemtrkr2_xch_max>0. && gemtrkr2_ych_max>0.) hgemtrkr_2_max_xy->Fill(gemtrkr2_xch_max, gemtrkr2_ych_max);
    if (gemtrkr3_xch_max>0. && gemtrkr3_ych_max>0.) hgemtrkr_3_max_xy->Fill(gemtrkr3_xch_max, gemtrkr3_ych_max);
    
    //=============== END SRS Data Processing & Correlations ==============
    
    //==================================================================================================
    //                    Process Fa125  Pulse  data
    //==================================================================================================
    
    #if 1
      
      f125_fit->Reset();
      mmg1_f125_fit->Reset();
      urw_f125_fit->Reset();
      #if (USE_PULSE>0)
        hevtk->Reset();
        hevtck->Reset();
      #endif
      
      double chi2cc_tgem=-999., chi2cc_mmg1=-999.,  chi2cc_urw=-999;
      double integral_tgem=0., integral_mmg1=0., integral_urw=0.;
      double tgem_ampmax=-1., mmg1_ampmax=-1., urw_xampmax=-1., urw_yampmax=-1., tgem_xchanmax=-1, mmg1_xchanmax=-1, urw_xchanmax=-1.,  urw_ychanmax=-1.;
      int tgem_timemax=0, mmg1_timemax=0, urw_xtimemax=0, urw_ytimemax=0;
      ULong64_t tgem_idx_x=0, mmg1_idx_x=0, urw_idx_x=0, urw_idx_y=0;
      double tgem_pos_x[f125_pulse_count], mmg1_pos_x[f125_pulse_count], urw_pos_x[f125_pulse_count], urw_pos_y[f125_pulse_count], tgem_amp_x[f125_pulse_count], mmg1_amp_x[f125_pulse_count], urw_amp_x[f125_pulse_count], urw_amp_y[f125_pulse_count];
      int gem_trk_hit=0, mmg1_trk_hit=0, urw_trk_hit=0;
      
      if (gt_3_idx_x>0 && gt_1_idx_x>0) {//(gt_1_idx_x>0 || gt_2_idx_x>0)) { //--External tracking condition
        
        Count("trk_hit");
        float x1=gemtrkr1_xch_max, x2=gemtrkr2_xch_max, x3=gemtrkr3_xch_max;
        float a=(x3-x1)/(z3-z1);
        float b=((x1)*z3-(x3)*z1)/(z3-z1);
        float gem_extr = a*zgem+b;
        float mmg1_extr = a*zmmg1+b;
        float urw_extr = a*zurw+b;
        f125_el_tracker_hits->Fill(gem_extr);
        mmg1_f125_el_tracker_hits->Fill(mmg1_extr);
        urw_f125_x_tracker_hits->Fill(urw_extr);
        
        for (ULong64_t i=0; i<f125_pulse_count; i++) { //--- Fadc125 Pulse Loop
          
        	float peak_amp = f125_pulse_peak_amp->at(i);
        	float ped = f125_pulse_pedestal->at(i);
        	if (0 > ped || ped > 200 ) ped = 100;
        	float amp=peak_amp-ped;
        	if (amp<0) amp=0;
        	float time=f125_pulse_peak_time->at(i);
        	int fADCSlot = f125_pulse_slot->at(i);
        	int fADCChan = f125_pulse_channel->at(i);
        	
        	int tripGemChan = Get3GEMChan(fADCChan, fADCSlot, RunNum);
        	int mmg1Chan = GetMMG1Chan(fADCChan, fADCSlot, RunNum);
          int urwXChan = GetRWELLXChan(fADCChan, fADCSlot, RunNum);
        	
        	if (tripGemChan>-1) {
            //amp = Get3GEMCalib(amp, tripGemChan, RunNum);
            if (amp>TGEM_THRESH) {
            float gemChan_x = tripGemChan*0.4+3.2; // to [mm]
            float gem_correction = 0.;
            if (RunNum < runAluminum) {
              gem_correction = -0.657655 + (gemChan_x)*0.00164153; //= -0.570621 + (gemChan_x)*0.00385134;
            } else {
              gem_correction = -0.657016 + (gemChan_x)*0.00255499;
            }
            tgem_residuals->Fill((gemChan_x-gem_extr));
            tgem_residualscorr->Fill((gemChan_x-gem_extr)-gem_correction);
            tgem_residual_ch->Fill(gemChan_x, (gemChan_x-gem_extr));
            tgem_residual_chcorr->Fill(gemChan_x, (gemChan_x-gem_extr-gem_correction));
            gem_trk_hit=0;
            if (abs(gemChan_x-gem_extr-gem_correction)<=5.) { //within 10 mm
              Count ("gem_trk_hit");
              if (!match) {
                f125_el_tracker_eff->Fill(gem_extr);
                match = true;
              }
              gem_trk_hit++;
            } //-- END within 10mm
            }
      	  }
      	  if (amp>MMG1_THRESH && mmg1Chan>-1) {
            //amp = GetMMGCalib(amp, mmg1Chan, RunNum);
            if (amp>MMG1_THRESH) {
            float mmg1Chan_x = mmg1Chan*0.4+3.2; //-- to [mm]
            float mmg1_correction = 0.;
            if (RunNum < runAluminum) {
              mmg1_correction = -0.182357 + (mmg1Chan_x)*0.00608014; //= -0.22918 + (mmg1Chan_x)*0.00544838;
            } else {
              mmg1_correction = -0.276239 + (mmg1Chan_x)*0.00523602;
            }
            mmg1_residuals->Fill((mmg1Chan_x-mmg1_extr));
            mmg1_residualscorr->Fill((mmg1Chan_x-mmg1_extr)-mmg1_correction);
            mmg1_residual_ch->Fill(mmg1Chan_x, (mmg1Chan_x-mmg1_extr));
            mmg1_residual_chcorr->Fill(mmg1Chan_x, (mmg1Chan_x-mmg1_extr)-mmg1_correction);
            mmg1_trk_hit=0;
            if (abs(mmg1Chan_x-mmg1_extr-mmg1_correction)<=5.) { //within 10 mm
              Count ("mmg1_trk_hit");
              if (!match_mmg1) {
                mmg1_f125_el_tracker_eff->Fill(mmg1_extr);
                match_mmg1 = true;
              }
              mmg1_trk_hit++;
            } //-- END within 10mm
            }
      	  }
          if (amp>URW_THRESH && urwXChan>-1) {
            float urwChan_x = urwXChan*0.8+3.2; //-- to [mm]
            float urw_correction = 0.;
            if (RunNum < runAluminum) {
              urw_correction = 0.688852 + (urwChan_x)*0.00359395; //= 0.973761 + (urwChan_x)*0.00708953;
            } else {
              urw_correction = 0.489856 + (urwChan_x)*0.00431468;
            }
            urw_x_residuals->Fill((urwChan_x-urw_extr));
            urw_x_residualscorr->Fill((urwChan_x-urw_extr)-urw_correction);
            urw_x_residual_ch->Fill(urwChan_x, (urwChan_x-urw_extr));
            urw_x_residual_chcorr->Fill(urwChan_x, (urwChan_x-urw_extr-urw_correction));
            urw_trk_hit=0;
            if (abs(urwChan_x-urw_extr-urw_correction)<=5.) { //within 10 mm
              Count ("urw_trk_hit");
              if (!match_urw) {
                urw_f125_x_tracker_eff->Fill(urw_extr);
                match_urw = true;
              }
              urw_trk_hit++;
            } //-- END within 10mm
      	  }
          
      	} //--- end Fa125 Pulse Loop ---
      } //-- End external tracker condition !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      
      //===================================================================
      //                    Chi^2 Fit Calculation
      //===================================================================
      
      char f125Title[80]; sprintf(f125Title,"GEM-TRD: Event=%lld Run=%d; z pos [time *8ns]; y pos [ch #]",jentry,RunNum);
      f125_fit->SetTitle(f125Title);
      if (f125_fit->GetEntries()!=0) {
        std::pair<Double_t, Double_t> fitResult  = TrkFit(f125_fit,fx1,"fx1",1);
        chi2cc_tgem = fitResult.first;
        integral_tgem = fitResult.second;
      }
      double a0 = fx1.GetParameter(0);
      double a1 = fx1.GetParameter(1);
      
      char mmg1f125Title[80]; sprintf(mmg1f125Title,"MMG1-TRD: Event=%lld Run=%d; z pos [time *8ns]; y pos [ch #]",jentry,RunNum);
      mmg1_f125_fit->SetTitle(mmg1f125Title);
      if (mmg1_f125_fit->GetEntries()!=0) {
        std::pair<Double_t, Double_t> fitResult  = TrkFit(mmg1_f125_fit,fx2,"fx2",1);
        chi2cc_mmg1 = fitResult.first;
        integral_mmg1 = fitResult.second;
      }
      double a0_mmg1 = fx2.GetParameter(0);
      double a1_mmg1 = fx2.GetParameter(1);
      
      char urwf125Title[80]; sprintf(urwf125Title,"uRWELL-TRD:  Event=%lld Run=%d; z pos [time *8ns]; y pos [ch #]",jentry,RunNum);
      urw_f125_fit->SetTitle(urwf125Title);
      if (urw_f125_fit->GetEntries()!=0) {
        std::pair<Double_t, Double_t> fitResult  = TrkFit(urw_f125_fit,fx2,"fx2",1);
        chi2cc_urw = fitResult.first;
        integral_urw = fitResult.second;
      }
      double a0_urw = fx2.GetParameter(0);
      double a1_urw = fx2.GetParameter(1);
      
      //==============================================
      //        First fa125 pulse loop
      //==============================================
      
      for (ULong64_t i=0; i<f125_pulse_count; i++) { //--- Fadc125 Pulse Loop
        
        float peak_amp = f125_pulse_peak_amp->at(i);
        float ped = f125_pulse_pedestal->at(i);
       	if (0 > ped || ped > 200 ) ped = 100;
       	float amp=peak_amp-ped;
       	if (amp<0) amp=0;
       	float time=f125_pulse_peak_time->at(i);
       	int fADCSlot = f125_pulse_slot->at(i);
       	int fADCChan = f125_pulse_channel->at(i);
       	
        if (/*RunNum<6154 && */time>177.) continue;
        
        int tripGemChan = Get3GEMChan(fADCChan, fADCSlot, RunNum);
        float tripGemChan_x = tripGemChan*0.4+3.2; // to [mm]
       	int mmg1Chan = GetMMG1Chan(fADCChan, fADCSlot, RunNum);
       	float mmg1Chan_x = mmg1Chan*0.4+3.2; // to [mm]
        int urwXChan = GetRWELLXChan(fADCChan, fADCSlot, RunNum);
       	float urwChan_x = urwXChan*0.8+3.2; // to [mm]
        int urwYChan = GetRWELLYChan(fADCChan, fADCSlot, RunNum);
       	float urwChan_y = urwYChan*0.8+3.2; // to [mm]
        
       	if (tripGemChan>-1) {
          //amp = Get3GEMCalib(amp, tripGemChan, RunNum);
          if (amp>TGEM_THRESH) {
          f125_el_amp2ds->Fill(time,tripGemChan,amp);
          if (tgem_ampmax<amp) {
            tgem_ampmax=amp;
            tgem_timemax=time;
            tgem_xchanmax=tripGemChan_x;
          }
          }
        }
        if (mmg1Chan>-1) {
          //amp = GetMMGCalib(amp, mmg1Chan, RunNum);
          if (amp>MMG1_THRESH) {
          mmg1_f125_el_amp2ds->Fill(time,mmg1Chan,amp);
          if (mmg1_ampmax<amp) {
            mmg1_ampmax=amp;
            mmg1_timemax=time;
            mmg1_xchanmax=mmg1Chan_x;
          }
          }
        }
        if (amp>URW_THRESH && urwXChan>-1) {
          urw_f125_x_amp2ds->Fill(time,urwXChan,amp);
          if (urw_xampmax<amp) {
            urw_xampmax=amp;
            urw_xtimemax=time;
            urw_xchanmax=urwChan_x;
          }
        }
        if (amp>URW_THRESH && urwYChan>-1) {
          urw_f125_y_amp2ds->Fill(time,urwYChan,amp);
          if (urw_yampmax<amp) {
            urw_yampmax=amp;
            urw_ytimemax=time;
            urw_ychanmax=urwChan_y;
          }
        }
      } //--END first f125 pulse loop
      
      if (tgem_ampmax > TGEM_THRESH) { f125_el_amp2d_max->Fill(tgem_timemax,tgem_xchanmax,tgem_ampmax); f125_xVSamp_max->Fill(tgem_xchanmax,tgem_ampmax); }
      if (mmg1_ampmax > MMG1_THRESH) { mmg1_f125_el_amp2d_max->Fill(mmg1_timemax,mmg1_xchanmax,mmg1_ampmax); mmg1_f125_xVSamp_max->Fill(mmg1_xchanmax,mmg1_ampmax); };
      if (urw_xampmax > URW_THRESH) { urw_f125_x_amp2d_max->Fill(urw_xtimemax,urw_xchanmax,urw_xampmax); urw_f125_xVSamp_max->Fill(urw_xchanmax,urw_xampmax); }
      if (urw_yampmax > URW_THRESH) urw_f125_y_amp2d_max->Fill(urw_ytimemax,urw_ychanmax,urw_yampmax);
      tgem_mmg1_max_xcorr->Fill(tgem_xchanmax,mmg1_xchanmax);
      tgem_urw_max_xcorr->Fill(urw_xchanmax,tgem_xchanmax);
      urw_mmg1_max_xcorr->Fill(urw_xchanmax,mmg1_xchanmax);
      htgem_max_xy->Fill(tgem_xchanmax,tgem_ych_max);
      hmmg1_max_xy->Fill(mmg1_xchanmax,mmg1_ych_max);
      if (urw_xampmax > URW_THRESH && urw_yampmax > URW_THRESH) hurw_max_xy->Fill(urw_xchanmax,urw_ychanmax);
      
      
      //==============================================
      //        Second fa125 pulse loop
      //==============================================
      
   /*   if (tgem_ampmax>TGEM_THRESH && (mmg1_ampmax>MMG1_THRESH || qgem_ampmax>QGEM_THRESH) && (urw_xampmax>URW_THRESH && urw_yampmax>URW_THRESH)) USE_MAX_TRD_HITS = true;
      if (((gemtrkr1_xamp_max>TRKR_THRESH && gemtrkr1_yamp_max>TRKR_THRESH) && (gemtrkr2_xamp_max>TRKR_THRESH && gemtrkr2_yamp_max>TRKR_THRESH)) || ((gemtrkr1_xamp_max>TRKR_THRESH && gemtrkr1_yamp_max>TRKR_THRESH) && (gemtrkr3_xamp_max>TRKR_THRESH && gemtrkr3_yamp_max>TRKR_THRESH)) || ((gemtrkr2_xamp_max>TRKR_THRESH && gemtrkr2_yamp_max>TRKR_THRESH) && (gemtrkr3_xamp_max>TRKR_THRESH && gemtrkr3_yamp_max>TRKR_THRESH))) USE_MAX_TRACKER_HITS = true;
      /////////if () USE_TRD_XCORR = true;
      */
      //-- SET EXTERNAL TRACKING REQUIREMENT --
      //if (USE_MAX_TRD_HITS) {
      //if (USE_MAX_TRACKER_HITS) {
      //if (USE_TRD_XCORR) {
      
      
      for (ULong64_t i=0; i<f125_pulse_count; i++) {
        
      	float peak_amp = f125_pulse_peak_amp->at(i);
      	float ped = f125_pulse_pedestal->at(i);
      	if (0 > ped || ped > 200 ) ped = 100;
      	float amp=peak_amp-ped;
      	float time=f125_pulse_peak_time->at(i);
      	int fADCSlot = f125_pulse_slot->at(i);
      	int fADCChan = f125_pulse_channel->at(i);
      	
        if ((!match && !match_mmg1) || time>177.) continue;
        
      	int tripGemChan = Get3GEMChan(fADCChan, fADCSlot, RunNum);
      	double tripGemChan_x = tripGemChan*0.4 + 3.2;
      	int mmg1Chan = GetMMG1Chan(fADCChan, fADCSlot, RunNum);
      	double mmg1Chan_x = mmg1Chan*0.4 + 3.2;
        int urwXChan = GetRWELLXChan(fADCChan, fADCSlot, RunNum);
       	double urwChan_x = urwXChan*0.8+3.2; // to [mm]
        int urwYChan = GetRWELLYChan(fADCChan, fADCSlot, RunNum);
       	double urwChan_y = urwYChan*0.8+3.2; // to [mm]
        
      	if (amp<0) amp=0;
      	
      	if (tripGemChan>-1) {
          //amp = Get3GEMCalib(amp, tripGemChan, RunNum);
          if (amp>TGEM_THRESH) {
            f125_fit->Fill(time,tripGemChan,amp);
            if (67.<=time && time<=147. && (abs(mmg1_xchanmax-(tripGemChan_x*0.908525-9.22087))<5. && abs(urw_xchanmax*0.881493-tripGemChan_x-7.886)<5.)) {
              tgem_pos_x[tgem_idx_x] = tripGemChan_x;
              tgem_amp_x[tgem_idx_x] = amp;
              tgem_idx_x++;
            }
            //if (abs(mmg1_xchanmax-(tripGemChan_x*0.908525-9.12579))<5. && abs(urw_xchanmax*0.881493-tripGemChan_x-7.886)<5.) { //-- Within 5mm
            if (match) { //tgem hit matched with external track
            hgemPulseDiff_el->Fill(mmg1_xchanmax-((tripGemChan_x*0.908525)-9.12579));
      	    f125_el->Fill(amp);
            f125_el_amp2d->Fill(time,tripGemChan,amp);
            f125_xVSamp->Fill(tripGemChan,amp);
            f125_timeVSamp->Fill(time,amp);
      	    tgem_xpos.push_back(tripGemChan);
      	    tgem_dedx.push_back(amp);
      	    tgem_zpos.push_back(time);
      	    tgem_nhit++;
      	    tgem_zHist->Fill(time,amp);
            bool foundRun = false;
            for (int j=0; j<listSize; j++) {
              if (RunNum==noRadList[j]) {foundRun=true; break;}
            }
            if (foundRun==true) {tgem_parID.push_back(0);}
            else {tgem_parID.push_back(1);}
            }
          }
    	  }
    	  if (mmg1Chan>-1) {
          //amp = GetMMGCalib(amp, mmg1Chan, RunNum);
          if (amp>MMG1_THRESH) {
            mmg1_f125_fit->Fill(time,mmg1Chan,amp);
            if (52.<=time && time<=152. && (abs(tgem_xchanmax*0.908525-mmg1Chan_x-9.12579)<5. && abs(urw_xchanmax*0.829002-mmg1Chan_x-18.4591))<5.) {
              mmg1_pos_x[mmg1_idx_x] = mmg1Chan_x;
              mmg1_amp_x[mmg1_idx_x] = amp;
              mmg1_idx_x++;
          }
            //if (abs(tgem_xchanmax*0.908525-mmg1Chan_x-9.12579)<5. && abs(urw_xchanmax*0.829002-mmg1Chan_x-18.4591)<5.) { //-- Within 5mm
            if (match_mmg1) { //mmg1 hit matched with external track
            hmmg1PulseDiff_el->Fill(tgem_xchanmax*0.908525-mmg1Chan_x-9.12579);
      	    mmg1_f125_el_amp2d->Fill(time,mmg1Chan,amp);
            mmg1_f125_xVSamp->Fill(mmg1Chan,amp);
            mmg1_f125_timeVSamp->Fill(time,amp);
      	    mmg1_f125_el->Fill(amp);
      	    mmg1_xpos.push_back(mmg1Chan);
      	    mmg1_dedx.push_back(amp);
      	    mmg1_zpos.push_back(time);
      	    mmg1_nhit++;
      	    mmg1_zHist->Fill(time,amp);
            bool foundRun = false;
            for (int j=0; j<listSize; j++) {
              if (RunNum==noRadList[j]) {foundRun=true; break;}
            }
            if (foundRun==true) {mmg1_parID.push_back(0);}
            else {mmg1_parID.push_back(1);}
            }
          }
        }
        if (amp>URW_THRESH && urwXChan>-1) {
            urw_f125_fit->Fill(time,urwXChan,amp);
            if (50.<=time && time<=140. && (abs(urwChan_x*0.881493-7.886-tgem_xchanmax)<5. && abs(urwChan_x*0.829002-mmg1_xchanmax-18.4591))<5.) {
              urw_pos_x[urw_idx_x] = urwChan_x;
              urw_amp_x[urw_idx_x] = amp; 
              urw_idx_x++;
            }
            if (abs(urwChan_x*0.881493-tgem_xchanmax-7.886)<5. && abs(urwChan_x*0.829002-mmg1_xchanmax-18.4591)<5.) {
            hurwPulseDiff_el->Fill(urwChan_x*0.881493-7.886-tgem_xchanmax);
            hurwPulseDiff_mmg->Fill(urwChan_x*0.829002-mmg1_xchanmax-18.4385);
            urw_f125_x_amp2d->Fill(time,urwXChan,amp);
            urw_f125_xVSamp->Fill(urwXChan,amp);
            urw_f125_x_timeVSamp->Fill(time,amp);
            urw_f125_el_x->Fill(amp);
            urw_xpos.push_back(urwXChan);
            urw_dedx.push_back(amp);
            urw_zpos.push_back(time);
            urw_nxhit++;
            urw_zHist->Fill(time,amp);
            bool foundRun = false;
            for (int j=0; j<listSize; j++) {
              if (RunNum==noRadList[j]) {foundRun=true; break;}
            }
            if (foundRun==true) {urw_parID.push_back(0);}
              else {urw_parID.push_back(1);}
            }
        }
        if (amp>URW_THRESH && urwYChan>-1) {
            if (50.<=time && time<=140. && (abs(urwChan_x*0.881493-tgem_xchanmax-7.886)<5. || abs(urwChan_x*0.829002-mmg1_xchanmax-18.4591)<5.)) {
              urw_pos_y[urw_idx_y] = urwChan_y;
              urw_amp_y[urw_idx_y] = amp;
              urw_idx_y++;
            }
            if (abs(urwChan_x*0.881493-tgem_xchanmax-7.886)<5. || abs(urwChan_x*0.829002-mmg1_xchanmax-18.4591)<5.) {
            urw_f125_y_amp2d->Fill(time,urwYChan,amp);
            urw_f125_yVSamp->Fill(urwYChan,amp);
            urw_f125_el_y->Fill(amp);
            //urw_ypos.push_back(urwYChan);
            //urw_dedx.push_back(amp);
            //urw_zpos.push_back(time);
            urw_nyhit++;
            //urw_zHist->Fill(time,amp);
            }
        }
    	} //--- end Fa125 Pulse Loop ---
      
      //} //END Selection
      
      //==================== Max Amplitude histos ============================
        if (tgem_ampmax>TGEM_THRESH) {
          f125_el_max->Fill(tgem_ampmax);
          f125_timeVSamp_max->Fill(tgem_timemax,tgem_ampmax);
          if (tgem_timemax>133) f125_el_max_late->Fill(tgem_ampmax);
        }
        if (mmg1_ampmax>MMG1_THRESH) {
          mmg1_f125_el_max->Fill(mmg1_ampmax);
          mmg1_f125_timeVSamp_max->Fill(mmg1_timemax,mmg1_ampmax);
          if (mmg1_timemax>140) mmg1_f125_el_max_late->Fill(mmg1_ampmax);
        }
        if (urw_xampmax>URW_THRESH) {
          urw_f125_el_xmax->Fill(urw_xampmax);
          urw_f125_x_timeVSamp_max->Fill(urw_xtimemax,urw_xampmax);
          if (urw_xtimemax>131) urw_f125_el_xmax_late->Fill(urw_xampmax);
        }
        if (urw_yampmax>URW_THRESH) {
          urw_f125_el_ymax->Fill(urw_yampmax);
        }
      
      tgem_amp_max=tgem_ampmax;
      tgem_time_max=tgem_timemax;
      tgem_xch_max=tgem_xchanmax;
      tgem_chi2cc.push_back(chi2cc_tgem);
      tgem_integral.push_back(integral_tgem);
      mmg1_amp_max=mmg1_ampmax;
      mmg1_time_max=mmg1_timemax;
      mmg1_xch_max=mmg1_xchanmax;
      mmg1_chi2cc.push_back(chi2cc_mmg1);
      mmg1_integral.push_back(integral_mmg1);
      urw_amp_max=urw_xampmax;
      urw_time_max=urw_xtimemax;
      urw_xch_max=urw_xchanmax;
      urw_chi2cc.push_back(chi2cc_urw);
      urw_integral.push_back(integral_urw);
      
      //if (abs(tgem_el_chan_max-(mmg1_xchanmax-1.35))<100. || abs(mmg1_el_chan_max-(tgem_xchanmax-1.35))<100.) gem_mmg1_max_xcorr->Fill(tgem_xchanmax, mmg1_xchanmax);
      
      for (int i=1; i<21; i++) {
        tgem_zHist_vect.push_back(tgem_zHist->GetBinContent(i));
        mmg1_zHist_vect.push_back(mmg1_zHist->GetBinContent(i));
        urw_zHist_vect.push_back(urw_zHist->GetBinContent(i));
      }
      
      if (tgem_nhit>0.) htgem_nhits->Fill(tgem_nhit);
      if (mmg1_nhit>0.) hmmg1_nhits->Fill(mmg1_nhit);
      if (urw_nxhit>0.) hurw_nxhits->Fill(urw_nxhit);
      if (urw_nyhit>0.) hurw_nyhits->Fill(urw_nyhit);
      if (gt1_nhit>0.) hgt1_nhits->Fill(gt1_nhit);
      if (gt2_nhit>0.) hgt2_nhits->Fill(gt2_nhit);
      if (gt3_nhit>0.) hgt3_nhits->Fill(gt3_nhit);
      
        if (gemtrkr1_xamp_max>TRKR_THRESH && tgem_ampmax>TGEM_THRESH) tgem_gt1_xcorr->Fill(tgem_xchanmax, gemtrkr1_xch_max);
        if (gemtrkr1_xamp_max>TRKR_THRESH && urw_xampmax>URW_THRESH) urw_gt1_xcorr->Fill(urw_xchanmax, gemtrkr1_xch_max);
        if (gemtrkr1_xamp_max>TRKR_THRESH && mmg1_ampmax>MMG1_THRESH) mmg1_gt1_xcorr->Fill(mmg1_xchanmax, gemtrkr1_xch_max);
        if (gemtrkr2_xamp_max>TRKR_THRESH && tgem_ampmax>TGEM_THRESH) tgem_gt2_xcorr->Fill(tgem_xchanmax, gemtrkr2_xch_max);
        if (gemtrkr2_xamp_max>TRKR_THRESH && urw_xampmax>URW_THRESH) urw_gt2_xcorr->Fill(urw_xchanmax, gemtrkr2_xch_max);
        if (gemtrkr2_xamp_max>TRKR_THRESH && mmg1_ampmax>MMG1_THRESH) mmg1_gt2_xcorr->Fill(mmg1_xchanmax, gemtrkr2_xch_max);
        if (gemtrkr3_xamp_max>TRKR_THRESH && tgem_ampmax>TGEM_THRESH) tgem_gt3_xcorr->Fill(tgem_xchanmax, gemtrkr3_xch_max);
        if (gemtrkr3_xamp_max>TRKR_THRESH && urw_xampmax>URW_THRESH) urw_gt3_xcorr->Fill(urw_xchanmax, gemtrkr3_xch_max);
        if (gemtrkr3_xamp_max>TRKR_THRESH && mmg1_ampmax>MMG1_THRESH) mmg1_gt3_xcorr->Fill(mmg1_xchanmax, gemtrkr3_xch_max);
        if (gemtrkr1_xamp_max>TRKR_THRESH && gemtrkr3_xamp_max>TRKR_THRESH) hTrackDiff->Fill(gemtrkr1_xch_max - gemtrkr3_xch_max);
      
      for (ULong64_t j=0; j<urw_idx_x; j++) {
          for (ULong64_t i=0; i<urw_idx_y; i++) {
            if (urw_amp_y[i]>URW_THRESH && urw_amp_x[j]>URW_THRESH) hurw_xy->Fill(urw_pos_x[j], urw_pos_y[i]);
          }
          for (ULong64_t i=0; i<tgem_idx_x; i++) {
            if (tgem_amp_x[i]>TGEM_THRESH && urw_amp_x[j]>URW_THRESH) urw_tgem_xcorr->Fill(urw_pos_x[j], tgem_pos_x[i]);
          }
          for (ULong64_t i=0; i<mmg1_idx_x; i++) {
            if (mmg1_amp_x[i]>MMG1_THRESH && urw_amp_x[j]>URW_THRESH) urw_mmg1_xcorr->Fill(urw_pos_x[j], mmg1_pos_x[i]);
          }
      }
      for (ULong64_t j=0; j<tgem_idx_x; j++) {
          for (ULong64_t i=0; i<tgem_idx_y; i++) {
            if (tgem_peak_y_height[i]>TRKR_THRESH+600. && tgem_amp_x[j]>TGEM_THRESH) htgem_xy->Fill(tgem_pos_x[j], tgem_peak_pos_y[i]);
          }
          for (ULong64_t i=0; i<mmg1_idx_x; i++) {
            if (mmg1_amp_x[i]>MMG1_THRESH && tgem_amp_x[j]>TGEM_THRESH) tgem_mmg1_xcorr->Fill(tgem_pos_x[j], mmg1_pos_x[i]);
          }
      }
      for (ULong64_t i=0; i<tgem_idx_y; i++) {
          for (ULong64_t j=0; j<mmg1_idx_y; j++) {
            if (mmg1_peak_y_height[j]>TRKR_THRESH+600. && tgem_peak_y_height[i]>TRKR_THRESH+600.) tgem_mmg1_ycorr->Fill(tgem_peak_pos_y[i], mmg1_peak_pos_y[j]);
          }
      }
      for (ULong64_t j=0; j<mmg1_idx_y; j++) {
          for (ULong64_t i=0; i<mmg1_idx_x; i++) {
            if (mmg1_amp_x[i]>MMG1_THRESH && mmg1_peak_y_height[j]>TRKR_THRESH+600.) hmmg1_xy->Fill(mmg1_pos_x[i], mmg1_peak_pos_y[j]);
          }
      }
      
      
    #endif
    //======================= End Process Fa125 Pulse data ================================
    
    //==================================================================================================
    //                    Process Fa125  RAW data
    //==================================================================================================
    
    #ifdef USE_125_RAW
      #ifdef VERBOSE
        if (jentry<MAX_PRINT) printf("------------------ Fadc125  wraw_count = %llu ---------\n", f125_wraw_count);
      #endif
      mhevt->Reset();
      mhevtc->Reset();
      mhevtf->Reset();
      hevt->Reset();
      hevtc->Reset();
      hevtf->Reset();
      
      for (ULong64_t i=0;i<f125_wraw_count; i++) { // --- fadc125 channels loop
        
        if (!match && !match_mmg1) continue;
        //cout<<"******** EVT# "<<jentry<<" RAW LOOP COND.: match=="<<match<<",  match_mmg1=="<<match_mmg1<<" **********"<<endl;
        int fadc_window = f125_wraw_samples_count->at(i);
        int fADCSlot = f125_wraw_slot->at(i);
        int fADCChan = f125_wraw_channel->at(i);
        int tripGemChan = Get3GEMChan(fADCChan, fADCSlot, RunNum);
        int mmg1Chan = GetMMG1Chan(fADCChan, fADCSlot, RunNum);
        int urwXChan = GetRWELLXChan(fADCChan, fADCSlot, RunNum);
        double DEDX_THR = TGEM_THRESH, mDEDX_THR = MMG1_THRESH, uDEDX_THR = URW_THRESH;
        int TimeWindowStart = 60;
        int TimeWindowStart_m = 50;
        int TimeMin = 0;
        int TimeMax = 140;
        
        //--ped calculation for f125 raw data
        int nped = 0, ped = 100;
        double ped_sum = 0.;
        for (int si=TimeWindowStart-15; si<TimeWindowStart; si++) {
          int ped_samp = f125_wraw_samples->at(f125_wraw_samples_index->at(i)+si);
          ped_sum += ped_samp;
          nped++;
        }
        ped = ped_sum / nped;
        if (0. > ped || ped > 200 ) ped = 100;
        //--FOR MMG
        int nped_m = 0, ped_m = 100;
        double ped_m_sum = 0.;
        for (int si=TimeWindowStart_m-15; si<TimeWindowStart_m; si++) {
          int ped_m_samp = f125_wraw_samples->at(f125_wraw_samples_index->at(i)+si);
          ped_m_sum += ped_m_samp;
          nped_m++;
        }
        ped_m = ped_m_sum / nped_m;
        if (0. > ped_m || ped_m > 200 ) ped_m = 100;
        
        for (int si=0; si<fadc_window; si++) {
          int time=si;
          int adc = f125_wraw_samples->at(f125_wraw_samples_index->at(i)+si);
          if (tripGemChan>-1) {
            adc = adc - ped;
            //adc = Get3GEMCalib(adc, tripGemChan, RunNum);
            if (adc>4090) printf("!!!!!!!!!!!!!!!!!!!!!! ADC 125 overflow: %d \n",adc);
            if (adc>DEDX_THR) {
              time-=TimeWindowStart;
              ///////////////if ( TimeMin > time || time > TimeMax ) continue; // --- drop early and late hits ---
              hevtc->SetBinContent(100-time,tripGemChan+1,adc/100.);
              hevt->SetBinContent(100-time,tripGemChan+1,adc/100.);
            }
          }
          if (mmg1Chan>-1) {
            adc = adc - ped_m;
            //adc = GetMMGCalib(adc, mmg1Chan, RunNum);
            if (adc > mDEDX_THR) {
              time-=(TimeWindowStart_m); //+40
              mhevtc->SetBinContent(120-time,mmg1Chan+1,adc/100.);
              mhevt->SetBinContent(120-time,mmg1Chan+1,adc/100.);
            }
          }
        } // --  end of samples loop
      } // -- end of fadc125 raw channels loop
      
      #ifdef SHOW_EVT_DISPLAY
          #if (USE_PULSE>0)
            c2->cd(1); hevtk->Draw("box");
          #else
            c2->cd(1); hevt->Draw("colz");
            c2->cd(6); mhevt->Draw("colz");
          #endif
          c2->cd(2);   hevtf->Draw("text");
          c2->cd(7);   mhevtf->Draw("text");
          c2->Modified(); c2->Update();
      #endif
      
      //==================================================================================================
      //            Begin NN Clustering & Track Fitting
      //==================================================================================================
      #if (USE_CLUST>0)
        // -------------------------------   hist dist clustering         ------------------------
        //--GEM-TRD
        float clust_Xmax[MAX_CLUST];
        float clust_Zmax[MAX_CLUST];
        float clust_Emax[MAX_CLUST];
        
        float clust_Xpos[MAX_CLUST];
        float clust_Zpos[MAX_CLUST];
        float clust_dEdx[MAX_CLUST];
        float clust_Size[MAX_CLUST];
        float clust_Width[MAX_CLUST][3];  // y1, y2, dy ; strips
        float clust_Length[MAX_CLUST][3]; // x1, x2, dx ; time
        float hits_Xpos[500];
        float hits_Zpos[500];
        float hits_dEdx[500];
        float hits_Size[MAX_CLUST];
        float hits_Width[MAX_CLUST];  // y1, y2, dy ; strips
        float hits_Length[MAX_CLUST]; // x1, x2, dx ; time
        //--MMG1TRD
        float mmg1_clust_Xmax[MAX_CLUST];
        float mmg1_clust_Zmax[MAX_CLUST];
        float mmg1_clust_Emax[MAX_CLUST];
        
        float mmg1_clust_Xpos[MAX_CLUST];
        float mmg1_clust_Zpos[MAX_CLUST];
        float mmg1_clust_dEdx[MAX_CLUST];
        float mmg1_clust_Size[MAX_CLUST];
        float mmg1_clust_Width[MAX_CLUST][3];  // y1, y2, dy ; strips
        float mmg1_clust_Length[MAX_CLUST][3]; // x1, x2, dx ; time
        float mmg1_hits_Xpos[500];
        float mmg1_hits_Zpos[500];
        float mmg1_hits_dEdx[500];
        float mmg1_hits_Size[MAX_CLUST];
        float mmg1_hits_Width[MAX_CLUST];  // y1, y2, dy ; strips
        float mmg1_hits_Length[MAX_CLUST]; // x1, x2, dx ; time
        
        for (int k=0; k<MAX_CLUST; k++) {
	        clust_Xpos[k]=0; clust_Zpos[k]=0; clust_dEdx[k]=0;  clust_Size[k]=0;
	        clust_Xmax[k]=0; clust_Zmax[k]=0; clust_Emax[k]=0;
          clust_Width[k][0]=999999;   	clust_Width[k][1]=-999999;   	clust_Width[k][2]=0;
          clust_Length[k][0]=999999;  	clust_Length[k][1]=-999999;  	clust_Length[k][2]=0;
          
          mmg1_clust_Xpos[k]=0; mmg1_clust_Zpos[k]=0; mmg1_clust_dEdx[k]=0;  mmg1_clust_Size[k]=0;
          mmg1_clust_Xmax[k]=0; mmg1_clust_Zmax[k]=0; mmg1_clust_Emax[k]=0;
          mmg1_clust_Width[k][0]=999999;     mmg1_clust_Width[k][1]=-999999;    mmg1_clust_Width[k][2]=0;
          mmg1_clust_Length[k][0]=999999;    mmg1_clust_Length[k][1]=-999999;   mmg1_clust_Length[k][2]=0;
        }
        int nclust=0, mmg1_nclust=0;
        #if (USE_PULSE>0)
          TH2F* hp = hevtk; // -- hevtk and hevtck should be same bin size
          TH2F* hpc = hevtck;
        #else
          TH2F* hp = hevt; // -- hevt and hevtc should be same bin size
          TH2F* hpc = hevtc;
          TH2F* hmp = mhevt;
          TH2F* hmpc = mhevtc;
        #endif
        //--GEM
        int nx=hp->GetNbinsX();    int ny=hp->GetNbinsY();
        double xmi=hp->GetXaxis()->GetBinLowEdge(1);     double xma=hp->GetXaxis()->GetBinUpEdge(nx);
        double ymi=hp->GetYaxis()->GetBinLowEdge(1);     double yma=hp->GetYaxis()->GetBinUpEdge(ny);
        double binx = (xma-xmi)/nx;      double biny = (yma-ymi)/ny;
        //--MMG1
        int nmx=hmp->GetNbinsX();    int nmy=hmp->GetNbinsY();
        double xmmi=hmp->GetXaxis()->GetBinLowEdge(1);   double xmma=hmp->GetXaxis()->GetBinUpEdge(nmx);
        double ymmi=hmp->GetYaxis()->GetBinLowEdge(1);   double ymma=hmp->GetYaxis()->GetBinUpEdge(nmy);
        double binmx = (xmma-xmmi)/nmx;      double binmy = (ymma-ymmi)/nmy;
        #ifdef VERBOSE
          printf("nx=%d,ny=%d,xmi=%f,xma=%f,ymi=%f,yma=%f\n",nx,ny,xmi,xma,ymi,yma);
        #endif
        #if (USE_PULSE>0)
          float CL_DIST=3.3; // mm
          double THR2 = 0.01;
        #else
          float CL_DIST=3.; //2.9; // mm
          double THR2 = 2.; //1.2
        #endif
        
        for (int iy=0; iy<ny; iy++) {  //-------------------- Clustering Loop (GEMTRD) ------------------------------------
          for (int ix=0; ix<nx; ix++) {
            double c1 = hpc->GetBinContent(ix+1,iy+1);                    // energy
            //if (c1>0) cout<<"      *** GEM-TRD CL ENERGY="<<c1<<" FOR ix,iy="<<ix<<endl;
            double x1=double(ix)/double(nx)*(xma-xmi)+xmi+binx/2.;    // drift time
            //if (c1>THR2) cout<<"      *** GEM-TRD CL DRIFT POSITION="<<x1<<" FOR ix,iy="<<ix<<endl;
            //cout<<"***GEM*** nx="<<nx<<" xma="<<xma<<" xmi="<<xmi<<" binx="<<binx<<endl;
            double y1=double(iy)/double(ny)*(yma-ymi)+ymi+biny/2.;    // X strip
            if (c1<THR2) continue;
            if (nclust==0) {
	            clust_Xpos[nclust]=y1; clust_Zpos[nclust]=x1;  clust_dEdx[nclust]=c1;  clust_Size[nclust]=1;
	            clust_Xmax[nclust]=y1; clust_Zmax[nclust]=x1;  clust_Emax[nclust]=c1;
              clust_Width[nclust][0]=y1;   	clust_Width[nclust][1]=y1;   	clust_Width[nclust][2]=0;
              clust_Length[nclust][0]=x1;  	clust_Length[nclust][1]=x1;  	clust_Length[nclust][2]=0;
              nclust++; continue;
            }
            int added=0;
            for (int k=0; k<nclust; k++) {
              double dist=sqrt(pow((y1-clust_Xpos[k]),2.)+pow((x1-clust_Zpos[k]),2.)); //--- dist hit to clusters
              //if (dist>0) cout<<"      *** GEM-TRD DIST="<<dist<<" FOR NCLUST="<<nclust<<" k="<<k<<endl;
              if (dist<CL_DIST) {
                clust_Xpos[k]=(y1*c1+clust_Xpos[k]*clust_dEdx[k])/(c1+clust_dEdx[k]);  //--  new X pos
                clust_Zpos[k]=(x1*c1+clust_Zpos[k]*clust_dEdx[k])/(c1+clust_dEdx[k]);  //--  new Z pos
	              if (c1>clust_Emax[k]) {
		              clust_Xmax[k]=y1;
		              clust_Zmax[k]=x1;
		              clust_Emax[k]=c1;
	              }
                clust_dEdx[k]=c1+clust_dEdx[k];  // new dEdx
                clust_Size[k]=1+clust_Size[k];  // clust size in pixels
                if (y1<clust_Width[k][0]) clust_Width[k][0]=y1; if (y1>clust_Width[k][1]) clust_Width[k][1]=y1; clust_Width[k][2]=clust_Width[k][1]-clust_Width[k][0];
                if (x1<clust_Length[k][0]) clust_Length[k][0]=x1;if (x1>clust_Length[k][1]) clust_Length[k][1]=x1;clust_Length[k][2]=clust_Length[k][1]-clust_Length[k][0];
                hpc->SetBinContent(ix,iy,k+1.);
                added=1; break;
              }
            }
            if (added==0) {
              if (nclust+1>=MAX_CLUST) continue;
	            clust_Xpos[nclust]=y1; clust_Zpos[nclust]=x1;  clust_dEdx[nclust]=c1;  clust_Size[nclust]=1;
	            clust_Xmax[nclust]=y1; clust_Zmax[nclust]=x1;  clust_Emax[nclust]=c1;
              clust_Width[nclust][0]=y1;   	clust_Width[nclust][1]=y1;   	clust_Width[nclust][2]=0;
              clust_Length[nclust][0]=x1;  	clust_Length[nclust][1]=x1;  	clust_Length[nclust][2]=0;
              nclust++;
            }
          }
        } //---------------------- End Clustering Loop (GEMTRD) ------------------------------
        
        for (int iy=0; iy<nmy; iy++) {  //-------------------- Clustering Loop (MMG1TRD) ------------------------------------
          for (int ix=0; ix<nmx; ix++) {
            double c1 = hmpc->GetBinContent(ix+1,iy+1);                       // energy
            //if (c1>0) cout<<"      *** MMG-TRD CL ENERGY="<<c1<<" FOR ix,iy="<<ix<<endl;
            double x1=double(ix)/double(nmx)*(xmma-xmmi)+xmmi+binmx/2.;   // drift time
            //cout<<"***MMG1*** nmx="<<nmx<<" xmma="<<xmma<<" xmmi="<<xmmi<<" binmx="<<binmx<<endl;
            double y1=double(iy)/double(nmy)*(ymma-ymmi)+ymmi+binmy/2.;        // X strip
            if (c1<THR2) continue;
            if (mmg1_nclust==0) {
              mmg1_clust_Xpos[mmg1_nclust]=y1; mmg1_clust_Zpos[mmg1_nclust]=x1;  mmg1_clust_dEdx[mmg1_nclust]=c1;  mmg1_clust_Size[mmg1_nclust]=1;
              mmg1_clust_Xmax[mmg1_nclust]=y1; mmg1_clust_Zmax[mmg1_nclust]=x1;  mmg1_clust_Emax[mmg1_nclust]=c1;
              mmg1_clust_Width[mmg1_nclust][0]=y1;    mmg1_clust_Width[mmg1_nclust][1]=y1;    mmg1_clust_Width[mmg1_nclust][2]=0;
              mmg1_clust_Length[mmg1_nclust][0]=x1;   mmg1_clust_Length[mmg1_nclust][1]=x1;   mmg1_clust_Length[mmg1_nclust][2]=0;
              mmg1_nclust++; continue;
            }
            int mmg1_added=0;
            for (int k=0; k<mmg1_nclust; k++) {
              double dist=sqrt(pow((y1-mmg1_clust_Xpos[k]),2.)+pow((x1-mmg1_clust_Zpos[k]),2.)); //--- dist hit to clusters
              //if (dist>0) cout<<"      *** MMG-TRD DIST="<<dist<<" FOR NCLUST="<<mmg1_nclust<<" k="<<k<<endl;
              if (dist<CL_DIST) {
                mmg1_clust_Xpos[k]=(y1*c1+mmg1_clust_Xpos[k]*mmg1_clust_dEdx[k])/(c1+mmg1_clust_dEdx[k]);  //--  new X pos
                mmg1_clust_Zpos[k]=(x1*c1+mmg1_clust_Zpos[k]*mmg1_clust_dEdx[k])/(c1+mmg1_clust_dEdx[k]);  //--  new Z pos
                if (c1>clust_Emax[k]) {
		              mmg1_clust_Xmax[k]=y1;
		              mmg1_clust_Zmax[k]=x1;
		              mmg1_clust_Emax[k]=c1;
	              }
                mmg1_clust_dEdx[k]=c1+mmg1_clust_dEdx[k];  // new dEdx
                mmg1_clust_Size[k]=1+mmg1_clust_Size[k];  // clust size in pixels
                if (y1<mmg1_clust_Width[k][0]) mmg1_clust_Width[k][0]=y1; if (y1>mmg1_clust_Width[k][1]) mmg1_clust_Width[k][1]=y1; mmg1_clust_Width[k][2]=mmg1_clust_Width[k][1]-mmg1_clust_Width[k][0];
                if (x1<mmg1_clust_Length[k][0]) mmg1_clust_Length[k][0]=x1;if (x1>mmg1_clust_Length[k][1]) mmg1_clust_Length[k][1]=x1; mmg1_clust_Length[k][2]=mmg1_clust_Length[k][1]-mmg1_clust_Length[k][0];
                hmpc->SetBinContent(ix,iy,k+1.);
                mmg1_added=1; break;
              }
            }
            if (mmg1_added==0) {
              if (mmg1_nclust+1>=MAX_CLUST) continue;
              mmg1_clust_Xpos[mmg1_nclust]=y1; mmg1_clust_Zpos[mmg1_nclust]=x1;  mmg1_clust_dEdx[mmg1_nclust]=c1;  mmg1_clust_Size[mmg1_nclust]=1;
              mmg1_clust_Xmax[mmg1_nclust]=y1; mmg1_clust_Zmax[mmg1_nclust]=x1;  mmg1_clust_Emax[mmg1_nclust]=c1;
              mmg1_clust_Width[mmg1_nclust][0]=y1;    mmg1_clust_Width[mmg1_nclust][1]=y1;    mmg1_clust_Width[mmg1_nclust][2]=0;
              mmg1_clust_Length[mmg1_nclust][0]=x1;   mmg1_clust_Length[mmg1_nclust][1]=x1;   mmg1_clust_Length[mmg1_nclust][2]=0;
              mmg1_nclust++;
            }
          }
        } //---------------------- End Clustering Loop (MMG1TRD) ------------------------------
        #if (USE_PULSE>0)
          int MinClustSize=1;
          double MinClustWidth=0.00;
          double MinClustLength=0.0;
          double MaxClustLength=5.; //5
          double zStart =  0.; // mm
          double zEnd   = 29.; // mm
        #else
          int MinClustSize=5;//1;
          double MinClustWidth=0.05;//0.001;
          double MinClustLength=0.01;
          double MaxClustLength=3.5; //4
          double zStart =  0.; // mm
          double zEnd   = 30.; // mm
          //double dEmin  = 2.; //
        #endif
        
        double maxClust_dEdx=0., totalClust_dEdx=0.;
        int ii=0;
        #ifdef VERBOSE
          printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
          printf("                Xpos   Ypos   Zpos       E    Width  Length   Size \n");
        #endif
        for (int k=0; k<nclust; k++) {
          #ifdef VERBOSE
            if (k<30) printf("%2d Clust(%2d): %6.1f %6.1f %8.1f %6.2f %6.2f %8.1f  ",k,k+1,clust_Xpos[k],clust_Zpos[k],clust_dEdx[k],clust_Width[k][2],clust_Length[k][2],clust_Size[k]);
          #endif
          //-------------  Cluster Filter (GEMTRD) -----------------
          if ((clust_Size[k]>=MinClustSize && zStart<clust_Zpos[k] && clust_Zpos[k]<zEnd && clust_Width[k][2]>MinClustWidth) || clust_Length[k][2]<MaxClustLength) {
            //if (abs(mmg1_xchanmax-((clust_Xpos[k]+3.2)*0.908525-9.22087))<5.) {
            //cout<<"EVENT="<<event_num<<" GEM_External_Track_Match="<<match<<" MMG_External_Track_Match="<<match_mmg1<<endl;
            if (match) {
        	    #if (USE_MAXPOS>0)
	              hgemClusterDiff_el->Fill(mmg1_xchanmax-((clust_Xmax[k]+3.2)*0.908525-9.22087));
                hits_Xpos[ii]=clust_Xmax[k];
	              hits_Zpos[ii]=clust_Zmax[k];
                //cout<<"  *** GEM-TRD HITS_ZPOSITION="<<clust_Zmax[k]<<" FOR k="<<k<<endl;
              #else
                hgemClusterDiff_el->Fill(mmg1_xchanmax-((clust_Xpos[k]+3.2)*0.908525-9.22087));
                hits_Xpos[ii]=clust_Xpos[k];
        	      hits_Zpos[ii]=clust_Zpos[k];
                //cout<<"  *** GEM-TRD HITS_ZPOSITION="<<clust_Zpos[k]<<" FOR k="<<k<<endl;
              #endif
        	    hits_dEdx[ii]=clust_dEdx[k];
              hits_Width[ii]=clust_Width[k][2];
              hits_Length[ii]=clust_Length[k][2];
        	    ii++;
        	    if (clust_dEdx[k]>maxClust_dEdx) maxClust_dEdx=clust_dEdx[k];
              totalClust_dEdx+=clust_dEdx[k];
              #ifdef VERBOSE
                if (k<30) printf("\n");
              #endif
            }
        	} else {
          #ifdef VERBOSE
            if (k<30) printf(" <--- skip \n");
          #endif
          }
        }
        int nhits=ii;
        clu_dedx_tot=totalClust_dEdx;
        
        double maxClust_m_dEdx=0., totalClust_m_dEdx=0.;
        int mmg1_ii=0;
        for (int k=0; k<mmg1_nclust; k++) {
          //-------------  Cluster Filter (MMG1TRD) -----------------
          if ((mmg1_clust_Size[k]>=MinClustSize && zStart<mmg1_clust_Zpos[k] && mmg1_clust_Zpos[k]<zEnd && mmg1_clust_Width[k][2]>MinClustWidth) || mmg1_clust_Length[k][2]<MaxClustLength+2.) {
            //if (abs(tgem_xchanmax*0.908525-(mmg1_clust_Xpos[k]+3.2)-9.22087)<5.) {
            if (match_mmg1) {
              #if (USE_MAXPOS>0)
                hmmg1ClusterDiff_el->Fill(tgem_xchanmax*0.908525-(mmg1_clust_Xpos[k]+3.2)-9.22087);
                mmg1_hits_Xpos[mmg1_ii]=mmg1_clust_Xmax[k];
                mmg1_hits_Zpos[mmg1_ii]=mmg1_clust_Zmax[k];
              #else
                hmmg1ClusterDiff_el->Fill(tgem_xchanmax*0.908525-(mmg1_clust_Xmax[k]+3.2)-9.22087);
                mmg1_hits_Xpos[mmg1_ii]=mmg1_clust_Xpos[k];
                mmg1_hits_Zpos[mmg1_ii]=mmg1_clust_Zpos[k];
              #endif
              mmg1_hits_dEdx[mmg1_ii]=mmg1_clust_dEdx[k];
              mmg1_hits_Width[mmg1_ii]=mmg1_clust_Width[k][2];
              mmg1_hits_Length[mmg1_ii]=mmg1_clust_Length[k][2];
              mmg1_ii++;
              if (mmg1_clust_dEdx[k]>maxClust_m_dEdx) maxClust_m_dEdx=mmg1_clust_dEdx[k];
              totalClust_m_dEdx+=mmg1_clust_dEdx[k];
            }
          }
        }
        int mmg1_nhits=mmg1_ii;
        mmg1_clu_dedx_tot=totalClust_m_dEdx;
        // ----------------------- end hist dist clustering ---------------------------------
        
        //=================================== Draw HITS and CLUST  ============================================
        #ifdef SHOW_EVTbyEVT
          char hevtTitle[80]; sprintf(hevtTitle,"GEM-TRD: Event=%lld Run=%d; z pos [mm]; y pos [mm]",jentry,RunNum);
          hevt->SetTitle(hevtTitle);
          #if (USE_PULSE>0)
            hevtk->SetTitle(hevtTitle);
          #endif
          char mhevtTitle[80]; sprintf(mhevtTitle,"MMG1-TRD: Event=%lld Run=%d; z pos [mm]; y pos [mm]",jentry,RunNum);
          mhevt->SetTitle(mhevtTitle);
          #ifdef VERBOSE
            printf("hits_SIZE=%d  Clust size = %d \n",nhits,nclust);
          #endif
          	c2->cd(1); gPad->Modified(); gPad->Update();
          	int COLMAP[]={1,2,3,4,6,5};
          	int pmt=22, pmt0=20; // PM type
            int max2draw=nclust;
            for (int i=0; i<max2draw; i++) {
              #if (USE_MAXPOS>0)
	              TMarker m = TMarker(clust_Zmax[i], clust_Xmax[i], pmt);
              #else
          	    TMarker m = TMarker(clust_Zpos[i], clust_Xpos[i], pmt);
              #endif
          	  int tcol=2;
          	  if (clust_Size[i]<MinClustSize) pmt=22; else pmt=pmt0;
          	  int mcol = COLMAP[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(pmt);
          	  m.SetMarkerSize(0.7+clust_dEdx[i]/300);
          	  m.DrawClone();
          	}
            gPad->Modified(); gPad->Update();
            
            //--MMG1TRD
            c2->cd(6); gPad->Modified(); gPad->Update();
            int mCOLMAP[]={1,2,3,4,6,5};
            int mpmt=22, mpmt0=20; // PM type
            int mmax2draw=mmg1_nclust;
            for (int i=0; i<mmax2draw; i++) {
              #if (USE_MAXPOS>0)
	              TMarker m = TMarker(mmg1_clust_Zmax[i], mmg1_clust_Xmax[i], mpmt);
              #else
                TMarker m = TMarker(mmg1_clust_Zpos[i], mmg1_clust_Xpos[i], mpmt);
              #endif
              int tcol=2;
              if (mmg1_clust_Size[i]<MinClustSize) mpmt=22; else mpmt=mpmt0;
              int mcol = mCOLMAP[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(mpmt);
              m.SetMarkerSize(0.7+mmg1_clust_dEdx[i]/300);
              m.DrawClone();
            }
            gPad->Modified(); gPad->Update();
        #endif
        
        #if (USE_GNN==1)   // GNN MC
          //----------------------------------------------------------
          //--   Send to Model simulation
          //----------------------------------------------------------
          
          #ifdef VERBOSE
            printf("**> Start Model simulation nclust=%d nhits=%d \n",nclust,nhits);
          #endif
          std::vector<int> tracks(nhits, 0);
          std::vector<float> Xcl;
          std::vector<float> Zcl;
          Xcl.clear();
          Zcl.clear();
          for (int n=0; n<nhits; n++) {
          	Xcl.push_back(hits_Xpos[n]);
          	Zcl.push_back(hits_Zpos[n]);
          }
          doPattern(Xcl, Zcl, tracks);  //---- call GNN ---
          //-- MMG1
          std::vector<int> mmg1_tracks(mmg1_nhits, 0);
          std::vector<float> mmg1_Xcl;
          std::vector<float> mmg1_Zcl;
          mmg1_Xcl.clear();
          mmg1_Zcl.clear();
          for (int n=0; n<mmg1_nhits; n++) {
            mmg1_Xcl.push_back(mmg1_hits_Xpos[n]);
            mmg1_Zcl.push_back(mmg1_hits_Zpos[n]);
          }
          doPattern(mmg1_Xcl, mmg1_Zcl, mmg1_tracks);  //---- call GNN ---
          #ifdef VERBOSE
            printf("**> End Model simulation \n"); //===================================================
          #endif
          #ifdef SHOW_EVTbyEVT
              c2->cd(2); gPad->Modified(); gPad->Update();
              int COLMAP2[]={1,2,3,4,6,5};
              for(ULong64_t i=0; i<tracks.size(); i++) {
                #ifdef VERBOSE
                  if (i<30) printf("i=%d trk=%d |  %8.2f,%8.2f\n",i, tracks[i], Xcl[i], Zcl[i]);
              	#endif
                TMarker m = TMarker(hits_Zpos[i], hits_Xpos[i], 24);
              	int tcol = min(tracks[i], 6);
              	int mcol = COLMAP2[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(41);     m.SetMarkerSize(1.5);
            	  m.DrawClone();  gPad->Modified(); gPad->Update();
              }
              //-- MMG1
              c2->cd(7); gPad->Modified(); gPad->Update();
              int mCOLMAP2[]={1,2,3,4,6,5};
              for(ULong64_t i=0; i<mmg1_tracks.size(); i++) {
                TMarker m = TMarker(mmg1_hits_Zpos[i], mmg1_hits_Xpos[i], 24);
                int tcol = min(mmg1_tracks[i], 6);
                int mcol = mCOLMAP2[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(41);     m.SetMarkerSize(1.5);
                m.DrawClone();  gPad->Modified(); gPad->Update();
              }
              #ifdef VERBOSE
                printf("\n\n");
                printf("**> End Cluster Plot \n");
              #endif
          #endif
          //--------------------------------------------------
          //----           Track fitting                 -----
          //--------------------------------------------------
          
          #ifdef VERBOSE
            printf("==> GNN: tracks sort  : trk_siz=%ld \r\n", tracks.size());
          #endif
          //-----------------   tracks sorting -------------
          std::vector<std::vector<float>> TRACKS;
          TRACKS.resize(nhits);
          //std::vector<float> hit_coord(2, 0);
          std::vector<int>  TRACKS_N(nhits, 0);
          for (int i=0; i<nhits; i++)  { TRACKS_N[i] = 0;  }
          //std::vector<float> xz(2,0);
          for (int i2=0; i2<nhits; i2++) {
          	int num =  tracks[i2];
          	int num2 = std::max(0, std::min(num, nhits - 1));
            #ifdef VERBOSE
              if (i2<20) printf("==> lstm3:track sort i=%d  : num=%d(%d) x=%f z=%f \n", i2, num, num2,  Xcl[i2],Zcl[i2]);
          	#endif
            //xz[0]=Xcl[i2];
            //xz[1]=Zcl[i2];
          	TRACKS[num2].push_back(Xcl[i2]);
            TRACKS[num2].push_back(Zcl[i2]);
          	TRACKS_N[num2]++;
          }
          //-- MMG1
          std::vector<std::vector<float>> mmg1_TRACKS;
          mmg1_TRACKS.resize(mmg1_nhits);
          //std::vector<float> mmg1_hit_coord(2, 0);
          std::vector<int>  mmg1_TRACKS_N(mmg1_nhits, 0);
          for (int i=0; i<mmg1_nhits; i++)  { mmg1_TRACKS_N[i] = 0;  }
          //std::vector<float> mmg1_xz(2, 0);
          for (int i2=0; i2<mmg1_nhits; i2++) {
            int num =  mmg1_tracks[i2];
            int num2 = std::max(0, std::min(num, mmg1_nhits - 1));
            //mmg1_xz[0]=mmg1_Xcl[i2];
            //mmg1_xz[1]=mmg1_Zcl[i2];
            mmg1_TRACKS[num2].push_back(mmg1_Xcl[i2]);
            mmg1_TRACKS[num2].push_back(mmg1_Zcl[i2]);
            mmg1_TRACKS_N[num2]++;
          }
          #if (DEBUG > 1)
            for (int i2 = 0; i2 < nhits; i2++) {
              printf(" trdID=%d n_hits=%d v_size=%d \n",i2,TRACKS_N[i2],TRACKS[i2].size());
              for (ULong64_t i3 = 0; i3 < TRACKS[i2].size(); i3+=2) {
                printf(" trkID=%d  hit=%d x=%f z=%f \n",i2,i3/2,TRACKS[i2].at(i3),TRACKS[i2].at(i3+1));
              }
              if ( TRACKS_N[i2]>0) printf("\n");
            }
          #endif
          //------------------ end tracks sorting --------------------
          
          #if (USE_FIT==1)
            //-----------------------------------
            //---       linear fitting        ---
            //-----------------------------------
            static TMultiGraph *mg;
            if (mg != NULL ) delete mg;
            mg = new TMultiGraph();
            int NTRACKS=0;
            int MIN_HITS=2;
            Double_t p0, p1;
            
            for (int i2=1; i2<nhits; i2++) {  //-- GEM tracks loop; zero track -> noise
              
             if (TRACKS_N[i2]<MIN_HITS) continue;   //---- select 2 (x,z) and more hits on track ----
            	#ifdef VERBOSE
                printf("==> fit: start trk: %d \r\n", i2);
              #endif
            	std::vector<Double_t> x;
            	std::vector<Double_t> y;
            	for (int i3=0; i3<(int)TRACKS[i2].size(); i3+=2) {
            	  #ifdef VERBOSE
                  printf(" trkID=%d  hit=%d x=%f z=%f \n",i2,i3/2,TRACKS[i2].at(i3),TRACKS[i2].at(i3+1));
            	  #endif
                x.push_back(TRACKS[i2].at(i3+1));
            	  y.push_back(TRACKS[i2].at(i3));
            	}
              #ifdef SHOW_EVTbyEVT
              	gErrorIgnoreLevel = kBreak; // Suppress warning messages from empty fit data
              	TGraph *g = new TGraph(TRACKS_N[i2], &x[0], &y[0]);  g->SetMarkerStyle(21); g->SetMarkerColor(i2);
              	TF1 *f = new TF1("f", "[1] * x + [0]");
              	g->Fit(f,"Q");
                //  --- get fit parameters ---
                TF1 *ffunc=g->GetFunction("f");
                p0=ffunc->GetParameter(0);
                p1=ffunc->GetParameter(1);
                Double_t chi2x_nn = ffunc->GetChisquare();
                Double_t Ndfx_nn = ffunc->GetNDF();
                double chi2nn=chi2x_nn/Ndfx_nn;
                #ifdef VERBOSE
                  printf("+++++>>  Track = %d fit: p0=%f p1=%f (%f deg) ff(15)=%f chi2nn=%f \n",i2,p0,p1,p1/3.1415*180.,ffunc->Eval(15.),chi2nn);
                #endif
            	  mg->Add(g,"p");
              #endif
              NTRACKS++;
            }  //-- end GEM tracks loop --
            
            //-- MMG1
            static TMultiGraph *mmg1_mg;
            if (mmg1_mg != NULL ) delete mmg1_mg;
            mmg1_mg = new TMultiGraph();
            int mmg1_NTRACKS=0;
            int mmg1_MIN_HITS=2;
            Double_t mp0, mp1;
            
            for (int i2=1; i2<mmg1_nhits; i2++) {  //-- MMG1 tracks loop; zero track -> noise
              
              if (mmg1_TRACKS_N[i2]<mmg1_MIN_HITS) continue;   //---- select 2 (x,z) or more hits on track ----
              std::vector<Double_t> mmg1_x;
              std::vector<Double_t> mmg1_y;
              for (int i3=0; i3<(int)mmg1_TRACKS[i2].size(); i3+=2) {
                mmg1_x.push_back(mmg1_TRACKS[i2].at(i3+1));
                mmg1_y.push_back(mmg1_TRACKS[i2].at(i3));
              }
              #ifdef SHOW_EVTbyEVT
                gErrorIgnoreLevel = kBreak; // Suppress warning messages from empty fit data
                TGraph *mmg1_g = new TGraph(mmg1_TRACKS_N[i2], &mmg1_x[0], &mmg1_y[0]);  mmg1_g->SetMarkerStyle(21); mmg1_g->SetMarkerColor(i2);
                TF1 *mmg1_f = new TF1("mmg1_f", "[1] * x + [0]");
                mmg1_g->Fit(mmg1_f,"Q");
                //  --- get fit parameters ---
                TF1 *mmg1_ffunc = mmg1_g->GetFunction("mmg1_f");
                mp0 = mmg1_ffunc->GetParameter(0);
                mp1 = mmg1_ffunc->GetParameter(1);
                Double_t mmg1_chi2x_nn = mmg1_ffunc->GetChisquare();
                Double_t mmg1_Ndfx_nn = mmg1_ffunc->GetNDF();
                double mmg1_chi2nn = mmg1_chi2x_nn/mmg1_Ndfx_nn;
                mmg1_mg->Add(mmg1_g, "p");
              #endif
              mmg1_NTRACKS++;
            }  //-- end MMG1 tracks loop --
            
            #ifdef SHOW_EVTbyEVT
              if (NTRACKS<1 && mmg1_NTRACKS<1) continue;  // --- skip event ----
              //if (nhits<3)    continue;  // --- skip event ----
              ///////////if (gem_trk_hit<1) continue;
                char mgTitle[80]; sprintf(mgTitle,"GEM ML-FPGA response, #Tracks=%d; z pos [mm]; y pos [mm]",NTRACKS);
                mg->SetTitle(mgTitle);
                c2->cd(3); mg->Draw("AP");
                mg->GetXaxis()->SetLimits(Xmin,Xmax);
                mg->SetMinimum(Ymin);
                mg->SetMaximum(Ymax);
                gPad->Modified(); gPad->Update();
                //-- MMG1
                char mmg1_mgTitle[80]; sprintf(mmg1_mgTitle,"MMG1 ML-FPGA response, #Tracks=%d; z pos [mm]; y pos [mm]",mmg1_NTRACKS);
                mmg1_mg->SetTitle(mmg1_mgTitle);
                c2->cd(8); mmg1_mg->Draw("AP");
                mmg1_mg->GetXaxis()->SetLimits(Xmin,Xmax);
                mmg1_mg->SetMinimum(Ymin);
                mmg1_mg->SetMaximum(Ymax);
                gPad->Modified(); gPad->Update();
            #endif
          #endif // USE_FIT
        #endif // USE_GNN MC
        
        //******************************************************************************
        #ifdef SHOW_EVTbyEVT
            cout<<"Event#="<<event_num<<" #ofGEMTracks="<<NTRACKS<<" #ofMMGTracks="<<mmg1_NTRACKS<<endl;
            cout<<"GEM_External_Track_Match="<<match<<" MMG_External_Track_Match="<<match_mmg1<<endl;
            #ifdef WRITE_CSV
              WriteToCSV(csvFile,event_num,NTRACKS,chi2cc_gem);
            #endif
            //c3->cd(1);  hCal_sum->Draw();         gPad->Modified();   gPad->Update();
            //c3->cd(2);  hCal_pulse->Draw("hist"); gPad->Modified();   gPad->Update();
            //c3->cd(3);  hCher_pulse->Draw("hist"); gPad->Modified();  gPad->Update();
            //c3->cd(4);  hPresh_pulse->Draw("hist"); gPad->Modified(); gPad->Update();
            //c3->cd(5);  hMult_pulse->Draw("hist"); gPad->Modified(); gPad->Update();
            
            c2->cd(4);  f125_fit->Draw("box");    gPad->Modified(); gPad->Update();
            c2->cd(5);  hevt->Draw("colz");       gPad->Modified();   gPad->Update();
            c2->cd(9);  mmg1_f125_fit->Draw("box"); gPad->Modified(); gPad->Update();
            c2->cd(10);  mhevt->Draw("colz");       gPad->Modified(); gPad->Update();
            printf("All done, click middle of canvas ...\n");
            #ifdef VERBOSE
              printf(" a0=%f a1=%f (%f deg)  fx1(150)=%f chi2cc_gem=%f  \n",a0,a1,a1/3.1415*180.,fx1.Eval(150.),chi2cc_gem);
            #endif
            if (match || match_mmg1/*NTRACKS>0 && mmg1_NTRACKS>0*/) c2->cd(1); gPad->WaitPrimitive();
        #endif
        if (maxClust_dEdx!=0.) hClusterMaxdEdx_el->Fill(maxClust_dEdx);
        if (totalClust_dEdx!=0.) hClusterTotaldEdx_el->Fill(totalClust_dEdx);
        if (maxClust_m_dEdx!=0.) hmmg1ClusterMaxdEdx_el->Fill(maxClust_m_dEdx);
        if (totalClust_m_dEdx!=0.) hmmg1ClusterTotaldEdx_el->Fill(totalClust_m_dEdx);
      #endif   // --- End if USE_CLUST>0 ---
    #endif   //=======================  End Fa125 RAW process Loop  =====================================
    
    //============ END GEMTRD Pattern Recognition Tracking ==================
    
    //=====================================================================================
    //===                Fill Root TTree Hits                                            ===
    //=====================================================================================
    
    #ifdef SAVE_TRACK_HITS
      #if USE_CLUST
      tgem_nclu=nhits;
      tgem_ntracks=NTRACKS;
      for (int n=0; n<nhits; n++) {
        clu_xpos.push_back(hits_Xpos[n]);
        clu_zpos.push_back(hits_Zpos[n]);
        clu_dedx.push_back(hits_dEdx[n]);
        clu_width.push_back(hits_Width[n]);
        if (hits_dEdx[n] > clu_dedx_max) {
          clu_dedx_max=hits_dEdx[n];
          clu_xpos_max=hits_Xpos[n];
          clu_zpos_max=hits_Zpos[n];
          clu_width_max=hits_Width[n];
        }
        //clu_length.push_back(hits_Length[n]);
        //f125_el_clu2d->Fill(hits_Zpos[n],hits_Xpos[n],hits_dEdx[n]);
      }
      
      mmg1_nclu=mmg1_nhits;
      mmg1_ntracks=mmg1_NTRACKS;
      for (int n=0; n<mmg1_nhits; n++) {
        mmg1_clu_xpos.push_back(mmg1_hits_Xpos[n]);
        mmg1_clu_zpos.push_back(mmg1_hits_Zpos[n]);
        mmg1_clu_dedx.push_back(mmg1_hits_dEdx[n]);
        mmg1_clu_width.push_back(mmg1_hits_Width[n]);
        if (mmg1_hits_dEdx[n] > mmg1_clu_dedx_max) {
          mmg1_clu_dedx_max=mmg1_hits_dEdx[n];
          mmg1_clu_xpos_max=mmg1_hits_Xpos[n];
          mmg1_clu_zpos_max=mmg1_hits_Zpos[n];
          mmg1_clu_width_max=mmg1_hits_Width[n];
        }
      }
      #endif
      if (tgem_nhit>0) EVENT_VECT_GEM->Fill();
      if (mmg1_nhit>0) EVENT_VECT_MMG1->Fill();
      if (urw_nxhit>0 /*|| urw_nyhit>0.*/) EVENT_VECT_URW->Fill();
    #endif
  } // ------------------------ END of event loop  ------------------------------
  
  timer.Stop();
  cout<<"***>>> End Event Loop, Elapsed Time:"<<endl; timer.Print();
  #ifdef WRITE_CSV
    csvFile.close();
  #endif
  cout<<" Total events = "<<(MaxEvt - FirstEvt)<<endl;
  //cout<<" # of electrons="<<el_count<<" # of pions="<<pi_count<<" # of atlas triggers="<<atlas_trigger_count<<endl;
  
  //---Drift time distribution plot ---
  TH1D *f125_drift = f125_el_amp2d->ProjectionX("f125_drift",100,135);
  TH1D *f125_drift_c = new TH1D(*f125_drift); f125_drift_c->SetStats(0);
  double tgemDriftScale = 1./f125_drift_c->GetEntries();
  f125_drift_c->Scale(tgemDriftScale);
  TH1D *mmg1_f125_drift = mmg1_f125_el_amp2d->ProjectionX("mmg1_f125_drift",95,130);
  TH1D *mmg1_f125_drift_c = new TH1D(*mmg1_f125_drift); mmg1_f125_drift_c->SetStats(0);
  double mmg1DriftScale = 1./mmg1_f125_drift_c->GetEntries();
  mmg1_f125_drift_c->Scale(mmg1DriftScale);
  TH1D *urw_f125_xdrift = urw_f125_x_amp2d->ProjectionX("urw_f125_xdrift",25,50);
  TH1D *urw_f125_xdrift_c = new TH1D(*urw_f125_xdrift); urw_f125_xdrift_c->SetStats(0);
  double urwxDriftScale = 1./urw_f125_xdrift_c->GetEntries();
  urw_f125_xdrift_c->Scale(urwxDriftScale);
  f125_drift_c->SetLineColor(4);  f125_drift_c->SetLineWidth(2);  HistList->Add(f125_drift_c);
  mmg1_f125_drift_c->SetLineColor(2); mmg1_f125_drift_c->SetLineWidth(2);  HistList->Add(mmg1_f125_drift_c);
  urw_f125_xdrift_c->SetLineColor(3); urw_f125_xdrift_c->SetLineWidth(2);  HistList->Add(urw_f125_xdrift_c);
  f125_drift_c->GetXaxis()->SetTitle("Time Response (*8ns)");
  f125_drift_c->GetYaxis()->SetTitle("1. / nEntries");
  f125_drift_c->SetTitle("Drift Time Distribution");
  TLegend *l1 = new TLegend(0.75,0.65,0.9,0.9);
  l1->SetNColumns(2);
  l1->SetTextSize(0.025);
  l1->AddEntry(f125_drift_c, "TRIP-GEM", "l");
  l1->AddEntry(mmg1_f125_drift_c, "MMG1", "l");
  l1->AddEntry(urw_f125_xdrift_c, "uRW X", "l");
  
  TCanvas *c4 = new TCanvas("c4","Drift Time Distribution", 1200, 800);
  c4->cd();
  f125_drift_c->Draw();
  mmg1_f125_drift_c->Draw("same");
  urw_f125_xdrift_c->Draw("same");
  l1->Draw();
  
  //=====================================================================================
  //===                 S A V E   H I S T O G R A M S                                ====
  //=====================================================================================
  TFile* fOut;
  #if ANALYZE_MERGED
    char rootFileName[256]; sprintf(rootFileName, "RootOutput/ps25/merged/Run_%06d_%0dEntries_Output.root",RunNum,nEntries);
  #else
    char rootFileName[256]; sprintf(rootFileName, "RootOutput/ps25/Run_%06d_Output.root",RunNum);
  #endif
  fOut = new TFile(rootFileName, "RECREATE");
  fOut->cd();
  cout<<"Writing Output File: "<<rootFileName<<endl;
  HistList->Write("HistDQM", TObject::kSingleKey);
  c4->Write();
  fOut->Close();
  delete fOut;
  
  //=====================================================================================
  //===                 S A V E   T R A C K   H I T   T T R E E S                    ====
  //=====================================================================================
  #ifdef SAVE_TRACK_HITS
    printf("Writing TTree Hit Info File... \n");
    fHits->cd();
    EVENT_VECT_GEM->Write();
    EVENT_VECT_MMG1->Write();
    EVENT_VECT_URW->Write();
    fHits->Close();
    printf("TTree File Written & Closed OK \n");
  #endif

  //=====================================================================================
  //===                 P L O T     H I S T O G R A M S                               ===
  //=====================================================================================
  #if ANALYZE_MERGED
    const char *OutputDir="RootOutput/ps25/merged";
  #else
    const char *OutputDir="RootOutput/ps25";
  #endif
  #ifdef SAVE_PDF
    char ctit[120];
    #if ANALYZE_MERGED
      sprintf(G_DIR,"%s/Run_%06d_%06dEntries",OutputDir,RunNum,nEntries);
    #else
      sprintf(G_DIR,"%s/Run_%06d",OutputDir,RunNum);
    #endif
    sprintf(ctit,"File=%s",G_DIR);
    bool COMPACT=false;
    TCanvas *cc;
    int nxd=3;
    int nyd=5;
    char pdfname[120];  sprintf(pdfname,"%s_evdisp.pdf",G_DIR);  //c0->Print(pdfname);
    
    //---------------------  page 1 --------------------
    htitle(" Count ");   // if (!COMPACT) cc=NextPlot(0,0);
    nxd=2; nyd=4;
    cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hcount->Draw();
    cc=NextPlot(nxd,nyd);  htgem_nhits->Draw();
    cc=NextPlot(nxd,nyd);  hmmg1_nhits->Draw();
    cc=NextPlot(nxd,nyd);  hurw_nxhits->Draw();
    cc=NextPlot(nxd,nyd);  hurw_nyhits->Draw();
    cc=NextPlot(nxd,nyd);  hgt1_nhits->Draw();
    cc=NextPlot(nxd,nyd);  hgt2_nhits->Draw();
    cc=NextPlot(nxd,nyd);  hgt3_nhits->Draw();
    
    //---------------------  page 2a --------------------
    htitle(" TRD Correlations  ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);  hgemtrkr_1_tgem->Draw("colz");
    cc=NextPlot(nxd,nyd);  hgemtrkr_1_mmg1->Draw("colz");
    cc=NextPlot(nxd,nyd);  tgem_mmg1_ycorr->Draw("colz");
    
    htitle(" TRD Correlations  ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);  tgem_mmg1_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  tgem_mmg1_max_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  urw_mmg1_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  urw_mmg1_max_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  urw_tgem_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  tgem_urw_max_xcorr->Draw("colz");
    
    
    //---------------------  page 2a --------------------
    htitle(" TRD Correlations  ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);  htgem_xy->Draw("colz");
    cc=NextPlot(nxd,nyd);  htgem_max_xy->Draw("colz");
    cc=NextPlot(nxd,nyd);  hmmg1_xy->Draw("colz");
    cc=NextPlot(nxd,nyd);  hmmg1_max_xy->Draw("colz");
    cc=NextPlot(nxd,nyd);  hurw_xy->Draw("colz");
    cc=NextPlot(nxd,nyd);  hurw_max_xy->Draw("colz");
    
    
    htitle(" X Correlations  ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);  tgem_gt1_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  tgem_gt2_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  tgem_gt3_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  urw_gt1_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  urw_gt2_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  urw_gt3_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  mmg1_gt1_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  mmg1_gt2_xcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);  mmg1_gt3_xcorr->Draw("colz");
    
    //---------------------  page 2a --------------------
    htitle(" SRS GEM-TRKR 1  ");   if (!COMPACT) cc=NextPlot(0,0);
    //nxd=2; nyd=4;
    cc=NextPlot(nxd,nyd); hgemtrkr_1_peak_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_1_max_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_1_peak_x->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_1_peak_y->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_1_peak_x_height->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_1_peak_y_height->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_1_max_xch->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_1_max_xamp->Draw();
    
    htitle(" SRS GEM-TRKR 2  ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd); hgemtrkr_2_peak_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_2_max_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_2_peak_x->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_2_peak_y->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_2_peak_x_height->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_2_peak_y_height->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_2_max_xch->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_2_max_xamp->Draw();
    
    htitle(" SRS GEM-TRKR 3  ");   if (!COMPACT) cc=NextPlot(0,0);
    //nxd=2; nyd=3;
    cc=NextPlot(nxd,nyd); hgemtrkr_3_peak_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_3_max_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_3_peak_x->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_3_peak_y->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_3_peak_x_height->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_3_peak_y_height->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_3_max_xch->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr_3_max_xamp->Draw();
    
    htitle(" SRS GEM / MMG1  ");   if (!COMPACT) cc=NextPlot(0,0);
    //nxd=2; nyd=3;
    cc=NextPlot(nxd,nyd); mmg1_peak_y->Draw();
    cc=NextPlot(nxd,nyd); hmmg1_peak_y_height->Draw();
    cc=NextPlot(nxd,nyd); tgem_peak_y->Draw();
    cc=NextPlot(nxd,nyd); htgem_peak_y_height->Draw();
    //TBox fbox(xbc1,ybc1,xbc2,ybc2);  //---- draw box cut ---
    //fbox.Draw("same");
    //fbox.SetLineColor(kRed);
    //fbox.SetFillStyle(0);
    //fbox.SetLineWidth(1);
    
    //--------------------- new page --------------------
    htitle(" fADC125 Pulse (Hit) Track Differences ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);  hgemPulseDiff_el->Draw();
    cc=NextPlot(nxd,nyd);  hmmg1PulseDiff_el->Draw();
    cc=NextPlot(nxd,nyd);  hurwPulseDiff_el->Draw();
    cc=NextPlot(nxd,nyd);  hurwPulseDiff_mmg->Draw();
    cc=NextPlot(nxd,nyd);  hTrackDiff->Draw();
    
    htitle(" fADC125 Raw (Clustering) Track Differences ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);  hgemClusterDiff_el->Draw();
    cc=NextPlot(nxd,nyd);  hmmg1ClusterDiff_el->Draw();
    cc=NextPlot(nxd,nyd);  hClusterMaxdEdx_el->Draw();
    cc=NextPlot(nxd,nyd);  hmmg1ClusterMaxdEdx_el->Draw();
    cc=NextPlot(nxd,nyd);  hClusterTotaldEdx_el->Draw();
    cc=NextPlot(nxd,nyd);  hmmg1ClusterTotaldEdx_el->Draw();
    
   //---------------------  page 3 --------------------
    htitle("  TRD (fa125) Amp Distributions ");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  f125_el->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  mmg1_f125_el->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  urw_f125_el_x->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  urw_f125_el_y->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  f125_el_max->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  mmg1_f125_el_max->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  urw_f125_el_xmax->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  urw_f125_el_ymax->Draw();
    
   //---------------------  page 3a --------------------
    htitle("  GEM-TRD (fa125) Amp 2D");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   f125_el_amp2ds->Draw("colz");
    cc=NextPlot(nxd,nyd);   f125_el_amp2d->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_el_amp2ds->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_el_amp2d->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_x_amp2ds->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_x_amp2d->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_y_amp2ds->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_y_amp2d->Draw("colz");
    
    htitle("  Amplitude vs Channel, 2D");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   f125_xVSamp->Draw("colz");
    cc=NextPlot(nxd,nyd);   f125_xVSamp_max->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_xVSamp->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_xVSamp_max->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_xVSamp->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_xVSamp_max->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_yVSamp->Draw("colz");
    
    htitle("  Amplitude vs Time, 2D");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   f125_timeVSamp->Draw("colz");
    cc=NextPlot(nxd,nyd);   f125_timeVSamp_max->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_timeVSamp->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_timeVSamp_max->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_x_timeVSamp->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_x_timeVSamp_max->Draw("colz");
    
    htitle("  fa125 Max Amp 2D");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   f125_el_amp2d_max->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_el_amp2d_max->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_x_amp2d_max->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_f125_y_amp2d_max->Draw("colz");
    
    
    htitle("  External Tracking");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   f125_el_tracker_hits->Draw();
    cc=NextPlot(nxd,nyd);   f125_el_tracker_eff->Draw();
    cc=NextPlot(nxd,nyd);   mmg1_f125_el_tracker_hits->Draw();
    cc=NextPlot(nxd,nyd);   mmg1_f125_el_tracker_eff->Draw();
    cc=NextPlot(nxd,nyd);   urw_f125_x_tracker_hits->Draw();
    cc=NextPlot(nxd,nyd);   urw_f125_x_tracker_eff->Draw();
    
    htitle("  External Tracking");    if (!COMPACT) cc=NextPlot(0,0);
    nxd=2; nyd=3;
    
    auto AddFitText = [](TF1* fitFunc) -> TPaveText* {
      TPaveText* pt = new TPaveText(0.15,0.75,0.5,0.85,"NDC");
      pt->SetFillColor(0);
      pt->SetTextColor(kRed);
      pt->SetTextSize(0.035);
      pt->SetBorderSize(1);
      pt->AddText(Form("pol0 = %.3f #pm %.3f",fitFunc->GetParameter(0), fitFunc->GetParError(0)));
      return pt;
    };
    
    cc=NextPlot(nxd,nyd);   if (f125_el_tracker_hits->GetEntries()!=0) { f125_el_tracker_hits->Sumw2(); f125_el_tracker_eff->Sumw2();  f125_el_tracker_eff->Divide(f125_el_tracker_hits);  f125_el_tracker_eff->SetMaximum(2.);  f125_el_tracker_eff->Draw();  f125_el_tracker_eff->Fit("pol0","Q","",33.,70.);  TF1* fitFunc=f125_el_tracker_eff->GetFunction("pol0");  fitFunc->SetLineColor(kRed);  fitFunc->Draw("same");  AddFitText(fitFunc)->Draw("same"); }
    cc=NextPlot(nxd,nyd);   if (mmg1_f125_el_tracker_hits->GetEntries()!=0) { mmg1_f125_el_tracker_hits->Sumw2(); mmg1_f125_el_tracker_eff->Sumw2(); mmg1_f125_el_tracker_eff->Divide(mmg1_f125_el_tracker_hits);  mmg1_f125_el_tracker_eff->SetMaximum(2.);  mmg1_f125_el_tracker_eff->Draw();  mmg1_f125_el_tracker_eff->Fit("pol0","Q","",15.,54.);  TF1* fitFunc=mmg1_f125_el_tracker_eff->GetFunction("pol0");  fitFunc->SetLineColor(kRed);  fitFunc->Draw("same");  AddFitText(fitFunc)->Draw("same"); }
    cc=NextPlot(nxd,nyd);   if (urw_f125_x_tracker_hits->GetEntries()!=0) { urw_f125_x_tracker_hits->Sumw2(); urw_f125_x_tracker_eff->Sumw2();  urw_f125_x_tracker_eff->Divide(urw_f125_x_tracker_hits);  urw_f125_x_tracker_eff->SetMaximum(2.);  urw_f125_x_tracker_eff->Draw();  urw_f125_x_tracker_eff->Fit("pol0","Q","",42.,85.);  TF1* fitFunc=urw_f125_x_tracker_eff->GetFunction("pol0");  fitFunc->SetLineColor(kRed);  fitFunc->Draw("same");  AddFitText(fitFunc)->Draw("same"); }
    
    htitle("  External Tracking");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   tgem_residuals->Draw();
    cc=NextPlot(nxd,nyd);   tgem_residualscorr->Draw();
    cc=NextPlot(nxd,nyd);   mmg1_residuals->Draw();
    cc=NextPlot(nxd,nyd);   mmg1_residualscorr->Draw();
    cc=NextPlot(nxd,nyd);   urw_x_residuals->Draw();
    cc=NextPlot(nxd,nyd);   urw_x_residualscorr->Draw();
    
    htitle("  External Tracking");    if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);   tgem_residual_ch->Draw("colz");
    cc=NextPlot(nxd,nyd);   tgem_residual_chcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_residual_ch->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_residual_chcorr->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_x_residual_ch->Draw("colz");
    cc=NextPlot(nxd,nyd);   urw_x_residual_chcorr->Draw("colz");
    
   //------------- MAX COMPARISONS ---------------
   
    htitle("  TRD (fa125) Max Amp Comparisons ");    if (!COMPACT) cc=NextPlot(0,0);
    f125_el_max_late->Scale(1./(f125_el_max_late->GetEntries()));  f125_el_max_late->SetLineColor(2); f125_el_max_late->SetStats(0);
    cc=NextPlot(nxd,nyd);  f125_el_max_late->Draw("same");
    mmg1_f125_el_max_late->Scale(1./(mmg1_f125_el_max_late->GetEntries()));  mmg1_f125_el_max_late->SetLineColor(2); mmg1_f125_el_max_late->SetStats(0);
    cc=NextPlot(nxd,nyd);  mmg1_f125_el_max_late->Draw("same");  
    urw_f125_el_xmax_late->Scale(1./(urw_f125_el_xmax_late->GetEntries()));  urw_f125_el_xmax_late->SetLineColor(2); urw_f125_el_xmax_late->SetStats(0);
    cc=NextPlot(nxd,nyd);  urw_f125_el_xmax_late->Draw("same");
    
    //--- close PDF file ----
    cc=NextPlot(-1,-1);
  #endif
  cout<<"========== END OF RUN "<<RunNum<<" =========="<<endl;
}
