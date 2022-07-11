TFile* analysisFile;
TTree* analysisTree;

const int fTrackSource = 0;
const int fClusterAlgo = 0;

TString analysisFileName;

const int fMaxNclusters = 100;
const int fMaxNtracks = 100;
const int fMaxNhepmcp = 1000;
const int fMaxNmcpart = 1000;

const int fMaxNmatch = fMaxNclusters*fMaxNtracks;

const int fNActiveCalo = 6;
int fActiveCaloID[fNActiveCalo] = {3, 10, 6, 7, 1, 8};
TString fCaloName[12] = {"FHCAL", "FEMC", "DRCALO", "EEMC", "CEMC", "EHCAL", "HCALIN", "HCALOUT", "LFHCAL", "EEMCG", "BECAL", "FOCAL"};

const int fMaxNtrackProj = fMaxNtracks*fNActiveCalo;


// Arrays for tree branches

int fNtracks;
int fTrackLepton;
int fTrackID[fMaxNtracks];
int fTrackTrueID[fMaxNtracks];
int fTrackCharge[fMaxNtracks];
double fTrackPx[fMaxNtracks];
double fTrackPy[fMaxNtracks];
double fTrackPz[fMaxNtracks];
double fTrackPionLL[fMaxNtracks];
double fTrackKaonLL[fMaxNtracks];
double fTrackProtonLL[fMaxNtracks];

int fNtrackProj;
int fTrackProjID[fMaxNtrackProj];
int fTrackProjLayer[fMaxNtrackProj];
double fTrackProjX[fMaxNtrackProj];
double fTrackProjY[fMaxNtrackProj];
double fTrackProjZ[fMaxNtrackProj];
	 
int fNclusters;
int fClusterLepton;
int fClusterCalo[fMaxNclusters];
int fClusterID[fMaxNclusters];
double fClusterE[fMaxNclusters];
double fClusterX[fMaxNclusters];
double fClusterY[fMaxNclusters];
double fClusterZ[fMaxNclusters];
double fClusterEta[fMaxNclusters];
double fClusterPhi[fMaxNclusters];
double fClusterM02[fMaxNclusters];
double fClusterM20[fMaxNclusters];

int fNmatch;
int fMatchTrackID[fMaxNmatch];
int fMatchCaloID[fMaxNmatch];
int fMatchClustID[fMaxNmatch];

int fNhepmcp;
int fHepmcpLepton;
double fHepmcpxB;
double fHepmcpQ2;
int fHepmcpPDG[fMaxNhepmcp];
int fHepmcpBCID[fMaxNhepmcp];
double fHepmcpPx[fMaxNhepmcp];
double fHepmcpPy[fMaxNhepmcp];
double fHepmcpPz[fMaxNhepmcp];
double fHepmcpE[fMaxNhepmcp];

int fNmcpart;
int fMcpartLepton;
int fMcpartPDG[fMaxNmcpart];
int fMcpartBCID[fMaxNmcpart];
double fMcpartPx[fMaxNmcpart];
double fMcpartPy[fMaxNmcpart];
double fMcpartPz[fMaxNmcpart];
double fMcpartE[fMaxNmcpart];

// Functions

void BuildFactory(TString);
void RunFactory();
void CloseFactory();

void InitializeAnalysisTree();

void SetTracks();
void SetTrackProjections();
void SetClusters();
void SetMatch();
void SetMCEvent();


