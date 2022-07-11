#include "analysistreefactory.h"

void RunFactory() {

	SetTracks();
	SetClusters();
	SetMCEvent();
		
	analysisTree->Fill();

}



void BuildFactory(TString fFileName) {

	analysisFileName = fFileName;
	InitializeAnalysisTree();

}



void CloseFactory() {

	analysisFile->cd();
	analysisTree->Write();
	analysisFile->Close();

}



void SetTracks() {

	fNtracks = 0;
	fTrackLepton = trackID_ScatteredElectron;
	for(int itr = 0; itr < _nTracks; itr++) {

		if(_track_source[itr] == fTrackSource) {

			fTrackID[fNtracks] = _track_ID[itr];
			fTrackTrueID[fNtracks] = _track_trueID[itr];
			fTrackCharge[fNtracks] = _track_charge[itr];
			fTrackPx[fNtracks] = _track_px[itr];
			fTrackPy[fNtracks] = _track_py[itr];
			fTrackPz[fNtracks] = _track_pz[itr];
			fTrackPionLL[fNtracks] = _track_pion_LL[itr];		
			fTrackKaonLL[fNtracks] = _track_kaon_LL[itr];		
			fTrackProtonLL[fNtracks] = _track_proton_LL[itr];		
	
			// This is used as array index above, so
			// only increment at the END!
			fNtracks++;

		}
	}

	// Call this in SetTracks to make sure that track
	// arrays are filled before track projections
	SetTrackProjections();
	SetMatch();
}



void SetTrackProjections() {

	fNtrackProj = 0;
	for(int itr = 0; itr < fNtracks; itr++) {

		int trID = fTrackID[itr];

		for(int icalo = 0; icalo < fNActiveCalo; icalo++) {

			int fCaloID = fActiveCaloID[icalo];
			int caloLayer = ReturnProjectionIndexForCalorimeter(fCaloID, false);
			
			for(int iproj = 0; iproj < _nProjections; iproj++) {
				
				if(_track_ProjTrackID[iproj] == trID && _track_ProjLayer[iproj] == caloLayer) {
					
					fTrackProjID[fNtrackProj] = trID;
					fTrackProjLayer[fNtrackProj] = fCaloID;
					fTrackProjX[fNtrackProj] = _track_Proj_x[iproj];
					fTrackProjY[fNtrackProj] = _track_Proj_y[iproj];
					fTrackProjZ[fNtrackProj] = _track_Proj_z[iproj];
			
					// This is used as array index above, so
					// only increment at the END!
					fNtrackProj++;
				}
			}
		} 
	}	
}



void SetClusters() {

	fClusterLepton = -1;
	int lepCaloID = -1; 
        int lepClustID = -1;
	if(trackID_ScatteredElectron !=-1) {
		lepCaloID      = _track_matchECal[trackID_ScatteredElectron].caloid;
		lepClustID     = _track_matchECal[trackID_ScatteredElectron].id;
	}

	fNclusters = 0;
	for(int icalo = 0; icalo < fNActiveCalo; icalo++) {
	
		int fCaloID = fActiveCaloID[icalo];

		for(int iclust = 0; iclust < _clusters_calo[fClusterAlgo][fCaloID].size(); iclust++) {	
		
			fClusterCalo[fNclusters] = fCaloID;
			fClusterID[fNclusters] = iclust;
			fClusterE[fNclusters] = (_clusters_calo[fClusterAlgo][fCaloID].at(iclust)).cluster_E;
			fClusterX[fNclusters] = (_clusters_calo[fClusterAlgo][fCaloID].at(iclust)).cluster_X;
			fClusterY[fNclusters] = (_clusters_calo[fClusterAlgo][fCaloID].at(iclust)).cluster_Y;
			fClusterZ[fNclusters] = (_clusters_calo[fClusterAlgo][fCaloID].at(iclust)).cluster_Z;
			fClusterEta[fNclusters] = (_clusters_calo[fClusterAlgo][fCaloID].at(iclust)).cluster_Eta;
			fClusterPhi[fNclusters] = (_clusters_calo[fClusterAlgo][fCaloID].at(iclust)).cluster_Phi;
			fClusterM02[fNclusters] = (_clusters_calo[fClusterAlgo][fCaloID].at(iclust)).cluster_M02;
			fClusterM20[fNclusters] = (_clusters_calo[fClusterAlgo][fCaloID].at(iclust)).cluster_M20;
	
			if(fCaloID == lepCaloID && iclust == lepClustID) {
				fClusterLepton = fNclusters;
			}
	
			// This is used as array index above, so
			// only increment at the END!
			fNclusters++;

		}
	}
}



void SetMatch() {

	fNmatch = 0;

	int caloID, clustID;

	for(int itr = 0; itr < fNtracks; itr++) {

		int trID = fTrackID[itr];
			
		// First look for EMCAL matches
		caloID = _track_matchECal[trID].caloid;
		clustID = _track_matchECal[trID].id;
		if(caloID != -1 && clustID != -1) {
			fMatchTrackID[fNmatch] = trID;
			fMatchCaloID[fNmatch] = caloID;
			fMatchClustID[fNmatch] = clustID;
			fNmatch++;
		}

		// Now look for HCAL matches
		caloID = _track_matchHCal[trID].caloid;
		clustID = _track_matchHCal[trID].id;
		if(caloID != -1 && clustID != -1) {
			fMatchTrackID[fNmatch] = trID;
			fMatchCaloID[fNmatch] = caloID;
			fMatchClustID[fNmatch] = clustID;
			fNmatch++;
		}


	}
}



void SetMCEvent() {

	fHepmcpxB = _hepmcp_x2;
	fHepmcpQ2 = _hepmcp_Q2;

	fNhepmcp = _nHepmcp;
	fHepmcpLepton = -1;
	for(int imcp = 0; imcp < _nHepmcp; imcp++) {

		// Check if this is the scattered lepton
		if(_hepmcp_status[imcp] == 1 && _hepmcp_PDG[imcp] == 11 && _hepmcp_m1[imcp] == 10001) {
			fHepmcpLepton = imcp;
		}

		fHepmcpPDG[imcp] = _hepmcp_PDG[imcp];
		fHepmcpBCID[imcp] = _hepmcp_BCID[imcp];
		fHepmcpPx[imcp] = _hepmcp_px[imcp];	
		fHepmcpPy[imcp] = _hepmcp_py[imcp];	
		fHepmcpPz[imcp] = _hepmcp_pz[imcp];	
		fHepmcpE[imcp] = _hepmcp_E[imcp];			
	}
	
	fNmcpart = _nMCPart;
	fMcpartLepton = -1;
	for(int imcp = 0; imcp < _nMCPart; imcp++) {

		// Check if this is the scattered lepton
		if(_mcpart_BCID[imcp] == _hepmcp_BCID[fHepmcpLepton] && _mcpart_PDG[imcp] == 11) {
			fMcpartLepton = imcp;
		}

		fMcpartPDG[imcp] = _mcpart_PDG[imcp];
		fMcpartBCID[imcp] = _mcpart_BCID[imcp];
		fMcpartPx[imcp] = _mcpart_px[imcp];	
		fMcpartPy[imcp] = _mcpart_py[imcp];	
		fMcpartPz[imcp] = _mcpart_pz[imcp];	
		fMcpartE[imcp] = _mcpart_E[imcp];
	}

}



void InitializeAnalysisTree() {

	analysisFile = new TFile(analysisFileName, "RECREATE");
	analysisTree = new TTree("T", "DIS analysis tree");	

	// These are branches pointing to values/arrays explicitly set in this macro
	// (want full control over how these are filled)
	analysisTree->Branch("track.n",		&fNtracks,	"track.n/I");
	analysisTree->Branch("track.lepton",	&fTrackLepton,	"track.lepton/I");
	analysisTree->Branch("track.ID",	fTrackID,	"track.ID[track.n]/I");
	analysisTree->Branch("track.trueID",	fTrackTrueID,	"track.trueID[track.n]/I");
	analysisTree->Branch("track.charge",	fTrackCharge,	"track.charge[track.n]/I");
	analysisTree->Branch("track.px",	fTrackPx,	"track.px[track.n]/D");
	analysisTree->Branch("track.py",	fTrackPy,	"track.py[track.n]/D");
	analysisTree->Branch("track.pz",	fTrackPz,	"track.pz[track.n]/D");
	analysisTree->Branch("track.pionLL", 	fTrackPionLL,	"track.pionLL[track.n]/D");
	analysisTree->Branch("track.kaonLL", 	fTrackKaonLL,	"track.kaonLL[track.n]/D");
	analysisTree->Branch("track.protonLL", 	fTrackProtonLL,	"track.protonLL[track.n]/D");

	analysisTree->Branch("trackProj.n",	&fNtrackProj,		"trackProj.n/I");
	analysisTree->Branch("trackProj.ID",	fTrackProjID,		"trackProj.ID[trackProj.n]/I");
	analysisTree->Branch("trackProj.layer",	fTrackProjLayer,	"trackProj.layer[trackProj.n]/I");
	analysisTree->Branch("trackProj.x",	fTrackProjX,		"trackProj.x[trackProj.n]/D");
	analysisTree->Branch("trackProj.y",	fTrackProjY,		"trackProj.y[trackProj.n]/D");
	analysisTree->Branch("trackProj.z",	fTrackProjZ,		"trackProj.z[trackProj.n]/D");
		
	analysisTree->Branch("cluster.n",	&fNclusters,	"cluster.n/I");
	analysisTree->Branch("cluster.lepton",	&fClusterLepton,"cluster.lepton/I");
	analysisTree->Branch("cluster.calo",	fClusterCalo,	"cluster.calo[cluster.n]/I");
	analysisTree->Branch("cluster.ID",	fClusterID,	"cluster.ID[cluster.n]/I");
	analysisTree->Branch("cluster.E",	fClusterE,	"cluster.E[cluster.n]/D");
	analysisTree->Branch("cluster.x",	fClusterX,	"cluster.x[cluster.n]/D");
	analysisTree->Branch("cluster.y",	fClusterY,	"cluster.y[cluster.n]/D");
	analysisTree->Branch("cluster.z",	fClusterZ,	"cluster.z[cluster.n]/D");
	analysisTree->Branch("cluster.eta",	fClusterEta,	"cluster.eta[cluster.n]/D");
	analysisTree->Branch("cluster.phi",	fClusterPhi,	"cluster.phi[cluster.n]/D");
	analysisTree->Branch("cluster.M02",	fClusterM02,	"cluster.M02[cluster.n]/D");
	analysisTree->Branch("cluster.M20",	fClusterM20,	"cluster.M20[cluster.n]/D");

	analysisTree->Branch("match.n",		&fNmatch,	"match.n/I");
	analysisTree->Branch("match.trackID",	fMatchTrackID,	"match.trackID[match.n]/I");
	analysisTree->Branch("match.caloID",	fMatchCaloID,	"match.caloID[match.n]/I");
	analysisTree->Branch("match.clustID",	fMatchClustID,	"match.clustID[match.n]/I");

	analysisTree->Branch("hepmcp.n",	&fNhepmcp,	"hepmcp.n/I");
	analysisTree->Branch("hepmcp.lepton",	&fHepmcpLepton,	"hepmcp.lepton/I");
	analysisTree->Branch("hepmcp.xB",	&fHepmcpxB,	"hepmcp.xB/D");
	analysisTree->Branch("hepmcp.Q2",	&fHepmcpQ2,	"hepmcp.Q2/D");
	analysisTree->Branch("hepmcp.PDG",	fHepmcpPDG,	"hepmcp.PDG[hepmcp.n]/I");	
	analysisTree->Branch("hepmcp.BCID",	fHepmcpBCID,	"hepmcp.BCID[hepmcp.n]/I");	
	analysisTree->Branch("hepmcp.px",	fHepmcpPx,	"hepmcp.px[hepmcp.n]/D");	
	analysisTree->Branch("hepmcp.py",	fHepmcpPy,	"hepmcp.py[hepmcp.n]/D");	
	analysisTree->Branch("hepmcp.pz",	fHepmcpPz,	"hepmcp.pz[hepmcp.n]/D");	
	analysisTree->Branch("hepmcp.E",	fHepmcpE,	"hepmcp.E[hepmcp.n]/D");	

	analysisTree->Branch("mcpart.n",	&fNmcpart,	"mcpart.n/I");
	analysisTree->Branch("mcpart.lepton",	&fMcpartLepton,	"mcpart.lepton/I");
	analysisTree->Branch("mcpart.PDG",	fMcpartPDG,	"mcpart.PDG[mcpart.n]/I");	
	analysisTree->Branch("mcpart.BCID",	fMcpartBCID,	"mcpart.BCID[mcpart.n]/I");	
	analysisTree->Branch("mcpart.px",	fMcpartPx,	"mcpart.px[mcpart.n]/D");	
	analysisTree->Branch("mcpart.py",	fMcpartPy,	"mcpart.py[mcpart.n]/D");	
	analysisTree->Branch("mcpart.pz",	fMcpartPz,	"mcpart.pz[mcpart.n]/D");	
	analysisTree->Branch("mcpart.E",	fMcpartE,	"mcpart.E[mcpart.n]/D");	

	// These are branches pointing to arrays/values set by treeProcessing
	analysisTree->Branch("vertex.x",	&_vertex_x,		"vertex.x/F");
	analysisTree->Branch("vertex.y",	&_vertex_y,		"vertex.y/F");
	analysisTree->Branch("vertex.z",	&_vertex_z,		"vertex.z/F");
	analysisTree->Branch("vertex.truex",	&_vertex_true_x,	"vertex.truex/F");
	analysisTree->Branch("vertex.truey",	&_vertex_true_y,	"vertex.truey/F");
	analysisTree->Branch("vertex.truez",	&_vertex_true_z,	"vertex.truez/F");

}


