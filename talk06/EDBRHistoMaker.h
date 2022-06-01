#include "TGraphAsymmErrors.h"
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TChain.h"
#include "math.h"
#include "TLorentzVector.h"
#include "roccor/RoccoR.cc"
#include "TRandom.h"
//#include "get_rochester_scale.C"
#define Pi 3.1415926
using namespace std;
/// The large arrays that were here are now GONE.
/// Instead, we have this helper that holds the
/// information of all our histograms.

class HistoFactory {
	public:
		std::vector<std::string> vars;
		std::vector<int> nBins;
		std::vector<float> minBin;
		std::vector<float> maxBin;
		void setHisto(std::string s, int n, float min, float max) {
			vars.push_back(s);
			nBins.push_back(n);
			minBin.push_back(min);
			maxBin.push_back(max);
		}
};

class EDBRHistoMaker {
	public:
		EDBRHistoMaker(TTree *tree = 0, TFile *fileTMP = 0, TH1F* hR1 = 0, std::string out_name_="", RoccoR* rc_in=0,TString isChannel="");
		virtual ~EDBRHistoMaker();
		/// This is the tree structure. This comes directly from MakeClass
		TTree* fChain;   //!pointer to the analyzed TTree or TChain
		TFile* fileTMP_;
		TH1F* hR1_;
		Int_t fCurrent; //!current Tree number in a TChain
		bool setUnitaryWeights_;
		// Declaration of leaf types
		Float_t HLT_SF;
		Float_t         mT_pf;
		Float_t         mT2_pf;
		Float_t         mT_puppi;
		Float_t         mT2_puppi;
		Float_t         ml1g;
		Float_t         ml2g;
                Bool_t          photon_flag;
                Bool_t          lepton_flag;
		UInt_t          lumi;
		Int_t           channel;
		Int_t           lep1_pid;
		Int_t           lep2_pid;
		Float_t         lep1pt;
		Float_t         lep2pt;
		Float_t         lep1eta;
		Float_t         lep2eta;
		Float_t         lep1phi;
		Float_t         lep2phi;
		Int_t           lep1_isprompt;
		Int_t           lep2_isprompt;
		Int_t           lep1_is_tight;
		Int_t           lep2_is_tight;
		Int_t           lep1_charge;
		Int_t           lep2_charge;
		Int_t           photon_selection;
		Int_t           n_loose_mu;
		Int_t           n_loose_ele;
		Int_t           n_photon;
		Float_t         photonet;
		Float_t         photoneta;
		Float_t         photonphi;
		Int_t           photon_isprompt;
		Int_t           photon_gen_matching;
		Float_t         mll;
		Float_t         mllg;
		Float_t         ptll;
		Float_t         mt;
		Float_t         met;
		Float_t         metup;
		Float_t         puppimet;
		Float_t         PuppiMET_pt;
		Float_t         PuppiMET_phi;
		Float_t         PuppiMET_phiJERDown;
		Float_t         PuppiMET_phiJERUp;
		Float_t         PuppiMET_phiJESDown;
		Float_t         PuppiMET_phiJESUp;
		Float_t         PuppiMET_ptJERDown;
		Float_t         PuppiMET_ptJERUp;
		Float_t         PuppiMET_ptJESDown;
		Float_t         PuppiMET_ptJESUp;
		Float_t         MET_T1_pt;
		Float_t         MET_T1_phi;
		Float_t         MET_T1Smear_pt;
		Float_t         MET_T1Smear_phi;
		Float_t         MET_T1_pt_jer1Up;
		Float_t         MET_T1_phi_jer1Up;
		Float_t         MET_T1Smear_pt_jer1Up;
		Float_t         MET_T1Smear_phi_jer1Up;
		Float_t         MET_T1_pt_jesTotalUp;
		Float_t         MET_T1_phi_jesTotalUp;
		Float_t         MET_T1Smear_pt_jesTotalUp;
		Float_t         MET_T1Smear_phi_jesTotalUp;
		Float_t         MET_T1_pt_jer1Down;
		Float_t         MET_T1_phi_jer1Down;
		Float_t         MET_T1Smear_pt_jer1Down;
		Float_t         MET_T1Smear_phi_jer1Down;
		Float_t         MET_T1_pt_jesTotalDown;
		Float_t         MET_T1_phi_jesTotalDown;
		Float_t         MET_T1Smear_pt_jesTotalDown;
		Float_t         MET_T1Smear_phi_jesTotalDown;
		Float_t         puppimetphi;
		Float_t         rawmet;
		Float_t         rawmetphi;
		Float_t         metphi;
		Float_t         gen_weight;
		Int_t           npu;
		Float_t         ntruepu;
		Int_t           n_pos;
		Int_t           n_minus;
		Int_t           n_num;
		Int_t           MET_pass;
		Int_t           npvs;
                Float_t         btag_weight_medium;
                Float_t         btag_weight_medium_up;
                Float_t         btag_weight_medium_down;
                Float_t         btag_weight_loose;
                Float_t         btag_weight_loose_up;
                Float_t         btag_weight_loose_down;
		Float_t         puWeight;
		Float_t         puWeightUp;
		Float_t         puWeightDown;
                Float_t         LHEScaleWeight[9];
                Float_t         LHEPdfWeight[103];
		Float_t         L1PreFiringWeight_Muon_Nom;
		Float_t         L1PreFiringWeight_Muon_StatDn;
		Float_t         L1PreFiringWeight_Muon_StatUp;
		Float_t         L1PreFiringWeight_Muon_SystDn;
		Float_t         L1PreFiringWeight_Muon_SystUp;
		Float_t         L1PreFiringWeight_Nom;
		Float_t         L1PreFiringWeight_Up;
		Float_t         L1PreFiringWeight_Dn;
                Bool_t          Flag_goodVertices;
                Bool_t          Flag_globalSuperTightHalo2016Filter;
                Bool_t          Flag_HBHENoiseFilter;
                Bool_t          Flag_HBHENoiseIsoFilter;
                Bool_t          Flag_EcalDeadCellTriggerPrimitiveFilter;
                Bool_t          Flag_BadPFMuonFilter;
                Bool_t          Flag_BadPFMuonDzFilter;
                Bool_t          Flag_eeBadScFilter;
                Bool_t          Flag_ecalBadCalibFilter;
		Bool_t          EcalDeadCellTriggerPrimitiveFilter;
		Float_t         PSWeight[4]; 
		Float_t         scalef;
                Float_t         ele_id_scale;
                Float_t         ele_id_scale_Up;
                Float_t         ele_id_scale_Down;
                Float_t         ele_reco_scale;
                Float_t         ele_reco_scale_Up;
                Float_t         ele_reco_scale_Down;
                Float_t         muon_id_scale;
                Float_t         muon_id_scale_Up;
                Float_t         muon_id_scale_Down;
                Float_t         muon_iso_scale;
                Float_t         muon_iso_scale_Up;
                Float_t         muon_iso_scale_Down;
                Float_t         photon_id_scale;
                Float_t         photon_id_scale_Up;
                Float_t         photon_id_scale_Down;
                Float_t         photon_veto_scale;
                Float_t         photon_veto_scale_Up;
                Float_t         photon_veto_scale_Down;
		Float_t         actualWeight;
                Float_t         drll;
                Float_t         drl1a;
                Float_t         drl2a;
                Float_t         yVlep;
                Float_t         phiVlep;
                Float_t         ptVlep;
		Char_t          HLT_ee;
		Char_t          HLT_mm;
		Char_t          HLT_emu;
		Int_t           n_bjets_loose_deepcsv;
		Int_t           n_bjets_loose_deepFlavB;
		Int_t           n_bjets20_loose_deepcsv;
		Int_t           n_bjets20_loose_deepFlavB;
		Int_t           n_bjets_medium_deepcsv;
		Int_t           n_bjets_medium_deepFlavB;
		Int_t           n_bjets20_medium_deepcsv;
		Int_t           n_bjets20_medium_deepFlavB;
		Int_t           njets50;
		Int_t           njets40;
		Int_t           njets30;
		Int_t           njets20;
		Int_t           njets50_pc;
		Int_t           njets40_pc;
		Int_t           njets30_pc;
		Int_t           njets20_pc;
		Int_t           n_bjets_loose_deepcsv_pc;
		Int_t           n_bjets_loose_deepFlavB_pc;
		Int_t           n_bjets20_loose_deepcsv_pc;
		Int_t           n_bjets20_loose_deepFlavB_pc;
		Int_t           n_bjets_medium_deepcsv_pc;
		Int_t           n_bjets_medium_deepFlavB_pc;
		Int_t           n_bjets20_medium_deepcsv_pc;
		Int_t           n_bjets20_medium_deepFlavB_pc;

		// List of branches
		TBranch        *b_mT_pf;   //!
		TBranch        *b_mT2_pf;   //!
		TBranch        *b_mT_puppi;   //!
		TBranch        *b_mT2_puppi;   //!
		TBranch        *b_ml1g;   //!
		TBranch        *b_ml2g;   //!
		TBranch        *b_lumi;   //!
		TBranch        *b_channel;   //!
		TBranch        *b_lep1_pid;   //!
		TBranch        *b_lep2_pid;   //!
		TBranch        *b_lep1pt;   //!
		TBranch        *b_lep2pt;   //!
		TBranch        *b_lep1eta;   //!
		TBranch        *b_lep2eta;   //!
		TBranch        *b_lep1phi;   //!
		TBranch        *b_lep2phi;   //!
		TBranch        *b_lep1_isprompt;   //!
		TBranch        *b_lep2_isprompt;   //!
		TBranch        *b_lep1_is_tight;   //!
		TBranch        *b_lep2_is_tight;   //!
		TBranch        *b_lep1_charge;   //!
		TBranch        *b_lep2_charge;   //!
		TBranch        *b_photon_selection;   //!
		TBranch        *b_n_loose_mu;   //!
		TBranch        *b_n_loose_ele;   //!
		TBranch        *b_n_photon;   //!
		TBranch        *b_photonet;   //!
		TBranch        *b_photoneta;   //!
		TBranch        *b_photonphi;   //!
		TBranch        *b_photon_isprompt;   //!
		TBranch        *b_photon_gen_matching;   //!
		TBranch        *b_drll;   //!
		TBranch        *b_drl1a;   //!
		TBranch        *b_drl2a;   //!
		TBranch        *b_mllg;   //!
		TBranch        *b_mll;   //!
		TBranch        *b_ptll;   //!
		TBranch        *b_mt;   //!
		TBranch        *b_met;   //!
		TBranch        *b_metup;   //!
		TBranch        *b_puppimet;   //!
		TBranch        *b_puppimetphi;   //!
		TBranch        *b_PuppiMET_pt;   //!
		TBranch        *b_PuppiMET_phi;   //!
		TBranch        *b_PuppiMET_T1_pt;   //!
		TBranch        *b_PuppiMET_phiJERDown;   //!
		TBranch        *b_PuppiMET_phiJERUp;   //!
		TBranch        *b_PuppiMET_phiJESDown;   //!
		TBranch        *b_PuppiMET_phiJESUp;   //!
		TBranch        *b_PuppiMET_ptJERDown;   //!
		TBranch        *b_PuppiMET_ptJERUp;   //!
		TBranch        *b_PuppiMET_ptJESDown;   //!
		TBranch        *b_PuppiMET_ptJESUp;   //!
		TBranch        *b_MET_T1_pt;   //!
		TBranch        *b_MET_T1_phi;   //!
		TBranch        *b_MET_T1Smear_pt;   //!
		TBranch        *b_MET_T1Smear_phi;   //!
		TBranch        *b_MET_T1_pt_jer1Up;   //!
		TBranch        *b_MET_T1_phi_jer1Up;   //!
		TBranch        *b_MET_T1Smear_pt_jer1Up;   //!
		TBranch        *b_MET_T1Smear_phi_jer1Up;   //!
		TBranch        *b_MET_T1_pt_jesTotalUp;   //!
		TBranch        *b_MET_T1_phi_jesTotalUp;   //!
		TBranch        *b_MET_T1Smear_pt_jesTotalUp;   //!
		TBranch        *b_MET_T1Smear_phi_jesTotalUp;   //!
		TBranch        *b_MET_T1_pt_jer1Down;   //!
		TBranch        *b_MET_T1_phi_jer1Down;   //!
		TBranch        *b_MET_T1Smear_pt_jer1Down;   //!
		TBranch        *b_MET_T1Smear_phi_jer1Down;   //!
		TBranch        *b_MET_T1_pt_jesTotalDown;   //!
		TBranch        *b_MET_T1_phi_jesTotalDown;   //!
		TBranch        *b_MET_T1Smear_pt_jesTotalDown;   //!
		TBranch        *b_MET_T1Smear_phi_jesTotalDown;   //!
		TBranch        *b_rawmet;   //!
		TBranch        *b_rawmetphi;   //!
		TBranch        *b_metphi;   //!
		TBranch        *b_gen_weight;   //!
		TBranch        *b_npu;   //!
		TBranch        *b_ntruepu;   //!
		TBranch        *b_n_pos;   //!
		TBranch        *b_n_minus;   //!
		TBranch        *b_n_num;   //!
		TBranch        *b_MET_pass;   //!
		TBranch        *b_npvs;   //!
		TBranch        *b_btag_weight_medium;   //!
		TBranch        *b_btag_weight_medium_up;   //!
		TBranch        *b_btag_weight_medium_down;   //!
		TBranch        *b_btag_weight_loose;   //!
		TBranch        *b_btag_weight_loose_up;   //!
		TBranch        *b_btag_weight_loose_down;   //!
		TBranch        *b_puWeight;   //!
		TBranch        *b_puWeightUp;   //!
		TBranch        *b_puWeightDown;   //!
                TBranch        *b_LHEScaleWeight;
                TBranch        *b_LHEPdfWeight;
		TBranch        *b_PSWeight;   //!
		TBranch        *b_L1PreFiringWeight_Muon_Nom;   //!
		TBranch        *b_L1PreFiringWeight_Muon_StatDn;   //!
		TBranch        *b_L1PreFiringWeight_Muon_StatUp;   //!
		TBranch        *b_L1PreFiringWeight_Muon_SystDn;   //!
		TBranch        *b_L1PreFiringWeight_Muon_SystUp;   //!
		TBranch        *b_L1PreFiringWeight_Nom;   //!
		TBranch        *b_L1PreFiringWeight_Up;   //!
		TBranch        *b_L1PreFiringWeight_Dn;   //!
                TBranch        *b_Flag_goodVertices;
                TBranch        *b_Flag_globalSuperTightHalo2016Filter;
                TBranch        *b_Flag_HBHENoiseFilter;
                TBranch        *b_Flag_HBHENoiseIsoFilter;
                TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;
                TBranch        *b_Flag_BadPFMuonFilter;
                TBranch        *b_Flag_BadPFMuonDzFilter;
                TBranch        *b_Flag_eeBadScFilter;
                TBranch        *b_Flag_ecalBadCalibFilter;
                TBranch        *b_EcalDeadCellTriggerPrimitiveFilter;
		TBranch        *b_scalef;   //!
                TBranch        *b_ele_id_scale;   //!
                TBranch        *b_ele_id_scale_Up;   //!
                TBranch        *b_ele_id_scale_Down;   //!
                TBranch        *b_ele_reco_scale;   //!
                TBranch        *b_ele_reco_scale_Up;   //!
                TBranch        *b_ele_reco_scale_Down;   //!
                TBranch        *b_muon_id_scale;   //!
                TBranch        *b_muon_id_scale_Up;   //!
                TBranch        *b_muon_id_scale_Down;   //!
                TBranch        *b_muon_iso_scale;   //!
                TBranch        *b_muon_iso_scale_Up;   //!
                TBranch        *b_muon_iso_scale_Down;   //!
                TBranch        *b_photon_id_scale;   //!
                TBranch        *b_photon_id_scale_Up;   //!
                TBranch        *b_photon_id_scale_Down;   //!
                TBranch        *b_photon_veto_scale;   //!
                TBranch        *b_photon_veto_scale_Up;   //!
                TBranch        *b_photon_veto_scale_Down;   //!
		TBranch        *b_n_bjets_loose_deepcsv;   //!
		TBranch        *b_n_bjets_loose_deepFlavB;   //!
		TBranch        *b_n_bjets20_loose_deepcsv;   //!
		TBranch        *b_n_bjets20_loose_deepFlavB;   //!
		TBranch        *b_n_bjets_medium_deepcsv;   //!
		TBranch        *b_n_bjets_medium_deepFlavB;   //!
		TBranch        *b_n_bjets20_medium_deepcsv;   //!
		TBranch        *b_n_bjets20_medium_deepFlavB;   //!
		TBranch        *b_njets50;   //!
		TBranch        *b_njets40;   //!
		TBranch        *b_njets30;   //!
		TBranch        *b_njets20;   //!
		TBranch        *b_njets50_pc;   //!
		TBranch        *b_njets40_pc;   //!
		TBranch        *b_njets30_pc;   //!
		TBranch        *b_njets20_pc;   //!
		TBranch        *b_n_bjets_loose_deepcsv_pc;   //!
		TBranch        *b_n_bjets_loose_deepFlavB_pc;   //!
		TBranch        *b_n_bjets20_loose_deepcsv_pc;   //!
		TBranch        *b_n_bjets20_loose_deepFlavB_pc;   //!
		TBranch        *b_n_bjets_medium_deepcsv_pc;   //!
		TBranch        *b_n_bjets_medium_deepFlavB_pc;   //!
		TBranch        *b_n_bjets20_medium_deepcsv_pc;   //!
		TBranch        *b_n_bjets20_medium_deepFlavB_pc;   //!
		TBranch        *b_HLT_ee;//!
		TBranch        *b_HLT_mm;//!
		TBranch        *b_HLT_emu;//!

		// Basic functions directly from MakeClass
		Int_t GetEntry(Long64_t entry);
		Long64_t LoadTree(Long64_t entry);
		void Init(TTree *tree);
		void Loop(std::string outFileName,float luminosity,int isBarrel,TString isChannel);
		void Loop_SFs_mc(std::string outFileNamei,float luminosity,int isBarrel,TString isChannel,TString year);
		// Our added functions
		void createAllHistos(TString isChannel);
		void printAllHistos();
		void saveAllHistos(std::string outFileName);
		void setUnitaryWeights(bool setuniw = false) {
			setUnitaryWeights_ = setuniw;
		}

		int check(float pt, std::vector<float> * ptZ) {
			int goodw = 1;
			for (unsigned int i = 0; i < ptZ->size(); i++) {
				if (pt == ptZ->at(i)) {
					goodw = 0;
					break;
				}
			}

			return goodw;
		}
		// Our added variables
		int nVars;
		bool debug_;

		///meng
		TTree *treename;
		TFile *fout;
		std::string out_name;
		void endjob(){
			fout->cd();
			treename->Write();
			fout->Write();
			fout->Close();
			delete fout;
		}
		///lu
		// fro rochester correction
		RoccoR rc;
		float get_rochester_scale(bool isdata, float charge_temp, float pt, float eta, float phi, int nl, float r1);


		// The histograms
		HistoFactory hs;
		std::map<std::string, TH1D*> theHistograms;
		TH2D *hmjmzz;
		TH1D *hmzzNEW;
};

void EDBRHistoMaker::Init(TTree *tree) {
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).

	// Set branch addresses and branch pointers
	if (!tree)
		return;
	fChain = tree;
	fCurrent = -1;
	fChain->SetMakeClass(1);

	fout = new TFile(out_name.c_str(), "RECREATE");
	treename = new TTree("outtree","outtree");
	cout<<"begin make outfile tree"<<endl;
	treename->Branch("channel", &channel,   "channel/I");
	treename->Branch("ptVlep", &ptVlep,   "ptVlep/F");
	treename->Branch("phiVlep", &phiVlep,   "phiVlep/F");
	treename->Branch("yVlep", &yVlep,   "yVlep/F");
	treename->Branch("lep1_pid", &lep1_pid, "lep1_pid/I");
	treename->Branch("lep2_pid", &lep2_pid, "lep1_pid/I");
	treename->Branch("lep1pt", &lep1pt,     "lep1pt/F");
	treename->Branch("lep2pt", &lep2pt,     "lep2pt/F");
	treename->Branch("lep1eta", &lep1eta,   "lep1eta/F");
	treename->Branch("lep2eta", &lep2eta,   "lep2eta/F");
	treename->Branch("lep1phi", &lep1phi,   "lep1phi/F");
	treename->Branch("lep2phi", &lep2phi,   "lep1phi/F");
	treename->Branch("lep1_isprompt", &lep1_isprompt, "lep1_isprompt/I");
	treename->Branch("lep2_isprompt", &lep2_isprompt, "lep2_isprompt/I");
	treename->Branch("lep1_is_tight", &lep1_is_tight, "lep1_is_tight/I");
	treename->Branch("lep2_is_tight", &lep2_is_tight, "lep2_is_tight/I");
	treename->Branch("lep1_charge", &lep1_charge, "lep1_charge/I");
	treename->Branch("lep2_charge", &lep2_charge, "lep2_charge/I");
	treename->Branch("photon_selection", &photon_selection, "photon_selection/I");
	treename->Branch("n_loose_mu", &n_loose_mu,   "n_loose_mu/I");
	treename->Branch("n_loose_ele", &n_loose_ele, "n_loose_ele/I");
	treename->Branch("n_photon", &n_photon, "n_photon/I");
	treename->Branch("photonet", &photonet,   "photonet/F");
	treename->Branch("photoneta", &photoneta, "photoneta/F");
	treename->Branch("photonphi", &photonphi, "photonphi/F");
	treename->Branch("photon_isprompt", &photon_isprompt, "photon_isprompt/I");
	treename->Branch("photon_gen_matching", &photon_gen_matching, "photon_gen_matching/I");
	treename->Branch("drll", &drll, "drll/F");
	treename->Branch("drl1a", &drl1a, "drl1a/F");
	treename->Branch("drl2a", &drl2a, "drl2a/F");
	treename->Branch("mll", &mll, "mll/F");
	treename->Branch("mllg", &mllg, "mllg/F");
        treename->Branch("ml1g", &ml1g, "ml1g/F");
        treename->Branch("ml2g", &ml2g, "ml2g/F");
	treename->Branch("ptll", &ptll, "ptll/F");
	treename->Branch("mt", &mt, "mt/F");
	treename->Branch("met", &met, "met/F");
	treename->Branch("metup", &metup, "metup/F");
	treename->Branch("puppimet", &puppimet, "puppimet/F");
	treename->Branch("puppimetphi", &puppimetphi, "puppimetphi/F");
	treename->Branch("rawmet", &rawmet, "rawmet/F");
	treename->Branch("rawmetphi", &rawmetphi, "rawmetphi/F");
	treename->Branch("metphi", &metphi, "metphi/F");
	treename->Branch("gen_weight", &gen_weight, "gen_weight/F");
	treename->Branch("npu", &npu, "npu/I");
	treename->Branch("ntruepu", &ntruepu, "ntruepu/F");
	treename->Branch("n_pos", &n_pos, "n_pos/I");
	treename->Branch("n_minus", &n_minus, "n_minus/I");
	treename->Branch("n_num", &n_num, "n_num/I");
	treename->Branch("MET_pass", &MET_pass, "MET_pass/I");
	treename->Branch("npvs", &npvs, "npvs/I");
	treename->Branch("btag_weight_medium",&btag_weight_medium,"btag_weight_medium/F");
	treename->Branch("btag_weight_medium_up",&btag_weight_medium_up,"btag_weight_medium_up/F");
	treename->Branch("btag_weight_medium_down",&btag_weight_medium_down,"btag_weight_medium_down/F");
	treename->Branch("btag_weight_loose",&btag_weight_loose,"btag_weight_loose/F");
	treename->Branch("btag_weight_loose_up",&btag_weight_loose_up,"btag_weight_loose_up/F");
	treename->Branch("btag_weight_loose_down",&btag_weight_loose_down,"btag_weight_loose_down/F");
	treename->Branch("puWeight", &puWeight, "puWeight/F");
	treename->Branch("puWeightUp", &puWeightUp, "puWeightUp/F");
	treename->Branch("puWeightDown", &puWeightDown, "puWeightDown/F");
	treename->Branch("scalef", &scalef, "scalef/F");
        treename->Branch("LHEScaleWeight", LHEScaleWeight, "LHEScaleWeight[9]/F");
        treename->Branch("LHEPdfWeight", LHEPdfWeight, "LHEPdfWeight[103]/F");
	treename->Branch("actualWeight", &actualWeight, "actualWeight/F");
        treename->Branch("ele_id_scale",&ele_id_scale,"ele_id_scale/F");
        treename->Branch("ele_id_scale_Up",&ele_id_scale_Up,"ele_id_scale_Up/F");
        treename->Branch("ele_id_scale_Down",&ele_id_scale_Down,"ele_id_scale_Down/F");
        treename->Branch("ele_reco_scale",&ele_reco_scale,"ele_reco_scale/F");
        treename->Branch("ele_reco_scale_Up",&ele_reco_scale_Up,"ele_reco_scale_Up/F");
        treename->Branch("ele_reco_scale_Down",&ele_reco_scale_Down,"ele_reco_scale_Down/F");
        treename->Branch("muon_id_scale",&muon_id_scale,"muon_id_scale/F");
        treename->Branch("muon_id_scale_Up",&muon_id_scale_Up,"muon_id_scale_Up/F");
        treename->Branch("muon_id_scale_Down",&muon_id_scale_Down,"muon_id_scale_Down/F");
        treename->Branch("muon_iso_scale",&muon_iso_scale,"muon_iso_scale/F");
        treename->Branch("muon_iso_scale_Up",&muon_iso_scale_Up,"muon_iso_scale_Up/F");
        treename->Branch("muon_iso_scale_Down",&muon_iso_scale_Down,"muon_iso_scale_Down/F");
        treename->Branch("photon_id_scale",&photon_id_scale,"photon_id_scale/F");
        treename->Branch("photon_id_scale_Up",&photon_id_scale_Up,"photon_id_scale_Up/F");
        treename->Branch("photon_id_scale_Down",&photon_id_scale_Down,"photon_id_scale_Down/F");
        treename->Branch("photon_veto_scale",&photon_veto_scale,"photon_veto_scale/F");
        treename->Branch("photon_veto_scale_Up",&photon_veto_scale_Up,"photon_veto_scale_Up/F");
        treename->Branch("photon_veto_scale_Down",&photon_veto_scale_Down,"photon_veto_scale_Down/F");
        treename->Branch("HLT_ee",&HLT_ee,"HLT_ee/B");
        treename->Branch("HLT_mm",&HLT_mm,"HLT_mm/B");
        treename->Branch("HLT_emu",&HLT_emu,"HLT_emu/B");
        treename->Branch("mT_pf", &mT_pf, "mT_pf/F");
        treename->Branch("HLT_SF", &HLT_SF, "HLT_SF/F");
        treename->Branch("mT2_pf", &mT2_pf, "mT2_pf/F");
        treename->Branch("mT_puppi", &mT_puppi, "mT_puppi/F");
        treename->Branch("mT2_puppi", &mT2_puppi, "mT2_puppi/F");
        treename->Branch("ml1g", &ml1g, "ml1g/F");
        treename->Branch("ml2g", &ml2g, "ml2g/F");
	treename->Branch("L1PreFiringWeight_Muon_Nom", &L1PreFiringWeight_Muon_Nom, "L1PreFiringWeight_Muon_Nom/F");
	treename->Branch("L1PreFiringWeight_Muon_StatDn", &L1PreFiringWeight_Muon_StatDn, "L1PreFiringWeight_Muon_StatDn/F");
	treename->Branch("L1PreFiringWeight_Muon_StatUp", &L1PreFiringWeight_Muon_StatUp, "L1PreFiringWeight_Muon_StatUp/F");
	treename->Branch("L1PreFiringWeight_Muon_SystDn", &L1PreFiringWeight_Muon_SystDn, "L1PreFiringWeight_Muon_SystDn/F");
	treename->Branch("L1PreFiringWeight_Muon_SystUp", &L1PreFiringWeight_Muon_SystUp, "L1PreFiringWeight_Muon_SystUp/F");
	treename->Branch("L1PreFiringWeight_Nom", &L1PreFiringWeight_Nom, "L1PreFiringWeight_Nom/F");
	treename->Branch("L1PreFiringWeight_Up", &L1PreFiringWeight_Up, "L1PreFiringWeight_Up/F");
	treename->Branch("L1PreFiringWeight_Dn", &L1PreFiringWeight_Dn, "L1PreFiringWeight_Dn/F");
        treename->Branch("Flag_goodVertices", &Flag_goodVertices, "Flag_goodVertices/B");
        treename->Branch("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter,"Flag_globalSuperTightHalo2016Filter/B");
        treename->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, "Flag_HBHENoiseFilter/B");
        treename->Branch("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, "Flag_HBHENoiseIsoFilter/B");
        treename->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, "Flag_EcalDeadCellTriggerPrimitiveFilter/B");
        treename->Branch("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, "Flag_BadPFMuonFilter/B");
        treename->Branch("Flag_BadPFMuonDzFilter", &Flag_BadPFMuonDzFilter, "Flag_BadPFMuonDzFilter/B");
        treename->Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter, "Flag_eeBadScFilter/B");
        treename->Branch("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, "Flag_ecalBadCalibFilter/B");
        treename->Branch("EcalDeadCellTriggerPrimitiveFilter", &EcalDeadCellTriggerPrimitiveFilter, "EcalDeadCellTriggerPrimitiveFilter/B");
	treename->Branch("PuppiMET_phi", &PuppiMET_phi, "PuppiMET_phi/F");
	treename->Branch("PuppiMET_phiJERDown", &PuppiMET_phiJERDown, "PuppiMET_phiJERDown/F");
	treename->Branch("PuppiMET_phiJERUp", &PuppiMET_phiJERUp, "PuppiMET_phiJERUp/F");
	treename->Branch("PuppiMET_phiJESDown", &PuppiMET_phiJESDown, "PuppiMET_phiJESDown/F");
	treename->Branch("PuppiMET_phiJESUp", &PuppiMET_phiJESUp, "PuppiMET_phiJESUp/F");
	treename->Branch("PuppiMET_pt", &PuppiMET_pt, "PuppiMET_pt/F");
	treename->Branch("PuppiMET_ptJERDown", &PuppiMET_ptJERDown, "PuppiMET_ptJERDown/F");
	treename->Branch("PuppiMET_ptJERUp", &PuppiMET_ptJERUp, "PuppiMET_ptJERUp/F");
	treename->Branch("PuppiMET_ptJESDown", &PuppiMET_ptJESDown, "PuppiMET_ptJESDown/F");
	treename->Branch("PuppiMET_ptJESUp", &PuppiMET_ptJESUp, "PuppiMET_ptJESUp/F");
	treename->Branch("MET_T1_pt", &MET_T1_pt, "MET_T1_pt/F");
	treename->Branch("MET_T1_phi", &MET_T1_phi, "MET_T1_phi/F");
	treename->Branch("MET_T1Smear_pt", &MET_T1Smear_pt, "MET_T1Smear_pt/F");
	treename->Branch("MET_T1Smear_phi", &MET_T1Smear_phi, "MET_T1Smear_phi/F");
	treename->Branch("MET_T1_pt_jer1Up", &MET_T1_pt_jer1Up, "MET_T1_pt_jer1Up/F");
	treename->Branch("MET_T1_phi_jer1Up", &MET_T1_phi_jer1Up, "MET_T1_phi_jer1Up/F");
	treename->Branch("MET_T1Smear_pt_jer1Up", &MET_T1Smear_pt_jer1Up, "MET_T1Smear_pt_jer1Up/F");
	treename->Branch("MET_T1Smear_phi_jer1Up", &MET_T1Smear_phi_jer1Up, "MET_T1Smear_phi_jer1Up/F");
	treename->Branch("MET_T1_pt_jesTotalUp", &MET_T1_pt_jesTotalUp, "MET_T1_pt_jesTotalUp/F");
	treename->Branch("MET_T1_phi_jesTotalUp", &MET_T1_phi_jesTotalUp, "MET_T1_phi_jesTotalUp/F");
	treename->Branch("MET_T1Smear_pt_jesTotalUp", &MET_T1Smear_pt_jesTotalUp, "MET_T1Smear_pt_jesTotalUp/F");
	treename->Branch("MET_T1Smear_phi_jesTotalUp", &MET_T1Smear_phi_jesTotalUp, "MET_T1Smear_phi_jesTotalUp/F");
	treename->Branch("MET_T1_pt_jer1Down", &MET_T1_pt_jer1Down, "MET_T1_pt_jer1Down/F");
	treename->Branch("MET_T1_phi_jer1Down", &MET_T1_phi_jer1Down, "MET_T1_phi_jer1Down/F");
	treename->Branch("MET_T1Smear_pt_jer1Down", &MET_T1Smear_pt_jer1Down, "MET_T1Smear_pt_jer1Down/F");
	treename->Branch("MET_T1Smear_phi_jer1Down", &MET_T1Smear_phi_jer1Down, "MET_T1Smear_phi_jer1Down/F");
	treename->Branch("MET_T1_pt_jesTotalDown", &MET_T1_pt_jesTotalDown, "MET_T1_pt_jesTotalDown/F");
	treename->Branch("MET_T1_phi_jesTotalDown", &MET_T1_phi_jesTotalDown, "MET_T1_phi_jesTotalDown/F");
	treename->Branch("MET_T1Smear_pt_jesTotalDown", &MET_T1Smear_pt_jesTotalDown, "MET_T1Smear_pt_jesTotalDown/F");
	treename->Branch("MET_T1Smear_phi_jesTotalDown", &MET_T1Smear_phi_jesTotalDown, "MET_T1Smear_phi_jesTotalDown/F");
	treename->Branch("n_bjets_loose_deepcsv", &n_bjets_loose_deepcsv, "n_bjets_loose_deepcsv/I");
	treename->Branch("n_bjets_loose_deepFlavB", &n_bjets_loose_deepFlavB, "n_bjets_loose_deepFlavB/I");
	treename->Branch("n_bjets20_loose_deepcsv", &n_bjets20_loose_deepcsv, "n_bjets20_loose_deepcsv/I");
	treename->Branch("n_bjets20_loose_deepFlavB", &n_bjets20_loose_deepFlavB, "n_bjets20_loose_deepFlavB/I");
	treename->Branch("n_bjets_medium_deepcsv", &n_bjets_medium_deepcsv, "n_bjets_medium_deepcsv/I");
	treename->Branch("n_bjets_medium_deepFlavB", &n_bjets_medium_deepFlavB, "n_bjets_medium_deepFlavB/I");
	treename->Branch("n_bjets20_medium_deepcsv", &n_bjets20_medium_deepcsv, "n_bjets20_medium_deepcsv/I");
	treename->Branch("n_bjets20_medium_deepFlavB", &n_bjets20_medium_deepFlavB, "n_bjets20_medium_deepFlavB/I");
	treename->Branch("njets50", &njets50, "njets50/I");
	treename->Branch("njets40", &njets40, "njets40/I");
	treename->Branch("njets30", &njets30, "njets30/I");
	treename->Branch("njets20", &njets20, "njets20/I");
	treename->Branch("njets50_pc", &njets50_pc, "njets50_pc/I");
	treename->Branch("njets40_pc", &njets40_pc, "njets40_pc/I");
	treename->Branch("njets30_pc", &njets30_pc, "njets30_pc/I");
	treename->Branch("njets20_pc", &njets20_pc, "njets20_pc/I");
	treename->Branch("n_bjets_loose_deepcsv_pc", &n_bjets_loose_deepcsv_pc, "n_bjets_loose_deepcsv_pc/I");
	treename->Branch("n_bjets_loose_deepFlavB_pc", &n_bjets_loose_deepFlavB_pc, "n_bjets_loose_deepFlavB_pc/I");
	treename->Branch("n_bjets20_loose_deepcsv_pc", &n_bjets20_loose_deepcsv_pc, "n_bjets20_loose_deepcsv_pc/I");
	treename->Branch("n_bjets20_loose_deepFlavB_pc", &n_bjets20_loose_deepFlavB_pc, "n_bjets20_loose_deepFlavB_pc/I");
	treename->Branch("n_bjets_medium_deepcsv_pc", &n_bjets_medium_deepcsv_pc, "n_bjets_medium_deepcsv_pc/I");
	treename->Branch("n_bjets_medium_deepFlavB_pc", &n_bjets_medium_deepFlavB_pc, "n_bjets_medium_deepFlavB_pc/I");
	treename->Branch("n_bjets20_medium_deepcsv_pc", &n_bjets20_medium_deepcsv_pc, "n_bjets20_medium_deepcsv_pc/I");
	treename->Branch("n_bjets20_medium_deepFlavB_pc", &n_bjets20_medium_deepFlavB_pc, "n_bjets20_medium_deepFlavB_pc/I");
	treename->Branch("PSWeight", PSWeight, "PSWeight/F");
	cout<<"make outfile tree end"<<endl;

	fChain->SetBranchAddress("mT_pf", &mT_pf, &b_mT_pf);
	fChain->SetBranchAddress("mT2_pf", &mT2_pf, &b_mT2_pf);
	fChain->SetBranchAddress("mT_puppi", &mT_puppi, &b_mT_puppi);
	fChain->SetBranchAddress("mT2_puppi", &mT2_puppi, &b_mT2_puppi);
	fChain->SetBranchAddress("ml1g", &ml1g, &b_ml1g);
	fChain->SetBranchAddress("ml2g", &ml2g, &b_ml2g);
	fChain->SetBranchAddress("channel", &channel, &b_channel);
	fChain->SetBranchAddress("lep1_pid", &lep1_pid, &b_lep1_pid);
	fChain->SetBranchAddress("lep2_pid", &lep2_pid, &b_lep2_pid);
	fChain->SetBranchAddress("lep1pt", &lep1pt, &b_lep1pt);
	fChain->SetBranchAddress("lep2pt", &lep2pt, &b_lep2pt);
	fChain->SetBranchAddress("lep1eta", &lep1eta, &b_lep1eta);
	fChain->SetBranchAddress("lep2eta", &lep2eta, &b_lep2eta);
	fChain->SetBranchAddress("lep1phi", &lep1phi, &b_lep1phi);
	fChain->SetBranchAddress("lep2phi", &lep2phi, &b_lep2phi);
	fChain->SetBranchAddress("lep1_isprompt", &lep1_isprompt, &b_lep1_isprompt);
	fChain->SetBranchAddress("lep2_isprompt", &lep2_isprompt, &b_lep2_isprompt);
	fChain->SetBranchAddress("lep1_is_tight", &lep1_is_tight, &b_lep1_is_tight);
	fChain->SetBranchAddress("lep2_is_tight", &lep2_is_tight, &b_lep2_is_tight);
	fChain->SetBranchAddress("lep1_charge", &lep1_charge, &b_lep1_charge);
	fChain->SetBranchAddress("lep2_charge", &lep2_charge, &b_lep2_charge);
	fChain->SetBranchAddress("photon_selection", &photon_selection, &b_photon_selection);
	fChain->SetBranchAddress("n_loose_mu", &n_loose_mu, &b_n_loose_mu);
	fChain->SetBranchAddress("n_loose_ele", &n_loose_ele, &b_n_loose_ele);
	fChain->SetBranchAddress("n_photon", &n_photon, &b_n_photon);
	fChain->SetBranchAddress("photonet", &photonet, &b_photonet);
	fChain->SetBranchAddress("photoneta", &photoneta, &b_photoneta);
	fChain->SetBranchAddress("photonphi", &photonphi, &b_photonphi);
	fChain->SetBranchAddress("photon_isprompt", &photon_isprompt, &b_photon_isprompt);
	fChain->SetBranchAddress("photon_gen_matching", &photon_gen_matching, &b_photon_gen_matching);
	fChain->SetBranchAddress("drll", &drll, &b_drll);
	fChain->SetBranchAddress("drl1a", &drl1a, &b_drl1a);
	fChain->SetBranchAddress("drl2a", &drl2a, &b_drl2a);
	fChain->SetBranchAddress("mll", &mll, &b_mll);
	fChain->SetBranchAddress("mllg", &mllg, &b_mllg);
	fChain->SetBranchAddress("ptll", &ptll, &b_ptll);
	fChain->SetBranchAddress("mt", &mt, &b_mt);
	fChain->SetBranchAddress("met", &met, &b_met);
	fChain->SetBranchAddress("metup", &metup, &b_metup);
	fChain->SetBranchAddress("puppimet", &puppimet, &b_puppimet);
	fChain->SetBranchAddress("puppimetphi", &puppimetphi, &b_puppimetphi);
	fChain->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi, &b_PuppiMET_phi);
	fChain->SetBranchAddress("PuppiMET_phiJERDown", &PuppiMET_phiJERDown, &b_PuppiMET_phiJERDown);
	fChain->SetBranchAddress("PuppiMET_phiJERUp", &PuppiMET_phiJERUp, &b_PuppiMET_phiJERUp);
	fChain->SetBranchAddress("PuppiMET_phiJESDown", &PuppiMET_phiJESDown, &b_PuppiMET_phiJESDown);
	fChain->SetBranchAddress("PuppiMET_phiJESUp", &PuppiMET_phiJESUp, &b_PuppiMET_phiJESUp);
	fChain->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt, &b_PuppiMET_pt);
	fChain->SetBranchAddress("PuppiMET_ptJERDown", &PuppiMET_ptJERDown, &b_PuppiMET_ptJERDown);
	fChain->SetBranchAddress("PuppiMET_ptJERUp", &PuppiMET_ptJERUp, &b_PuppiMET_ptJERUp);
	fChain->SetBranchAddress("PuppiMET_ptJESDown", &PuppiMET_ptJESDown, &b_PuppiMET_ptJESDown);
	fChain->SetBranchAddress("PuppiMET_ptJESUp", &PuppiMET_ptJESUp, &b_PuppiMET_ptJESUp);
	fChain->SetBranchAddress("MET_T1_pt", &MET_T1_pt, &b_MET_T1_pt);
	fChain->SetBranchAddress("MET_T1_phi", &MET_T1_phi, &b_MET_T1_phi);
	fChain->SetBranchAddress("MET_T1Smear_pt", &MET_T1Smear_pt, &b_MET_T1Smear_pt);
	fChain->SetBranchAddress("MET_T1Smear_phi", &MET_T1Smear_phi, &b_MET_T1Smear_phi);
	fChain->SetBranchAddress("MET_T1_pt_jer1Up", &MET_T1_pt_jer1Up, &b_MET_T1_pt_jer1Up);
	fChain->SetBranchAddress("MET_T1_phi_jer1Up", &MET_T1_phi_jer1Up, &b_MET_T1_phi_jer1Up);
	fChain->SetBranchAddress("MET_T1Smear_pt_jer1Up", &MET_T1Smear_pt_jer1Up, &b_MET_T1Smear_pt_jer1Up);
	fChain->SetBranchAddress("MET_T1Smear_phi_jer1Up", &MET_T1Smear_phi_jer1Up, &b_MET_T1Smear_phi_jer1Up);
	fChain->SetBranchAddress("MET_T1_pt_jesTotalUp", &MET_T1_pt_jesTotalUp, &b_MET_T1_pt_jesTotalUp);
	fChain->SetBranchAddress("MET_T1_phi_jesTotalUp", &MET_T1_phi_jesTotalUp, &b_MET_T1_phi_jesTotalUp);
	fChain->SetBranchAddress("MET_T1Smear_pt_jesTotalUp", &MET_T1Smear_pt_jesTotalUp, &b_MET_T1Smear_pt_jesTotalUp);
	fChain->SetBranchAddress("MET_T1Smear_phi_jesTotalUp", &MET_T1Smear_phi_jesTotalUp, &b_MET_T1Smear_phi_jesTotalUp);
	fChain->SetBranchAddress("MET_T1_pt_jer1Down", &MET_T1_pt_jer1Down, &b_MET_T1_pt_jer1Down);
	fChain->SetBranchAddress("MET_T1_phi_jer1Down", &MET_T1_phi_jer1Down, &b_MET_T1_phi_jer1Down);
	fChain->SetBranchAddress("MET_T1Smear_pt_jer1Down", &MET_T1Smear_pt_jer1Down, &b_MET_T1Smear_pt_jer1Down);
	fChain->SetBranchAddress("MET_T1Smear_phi_jer1Down", &MET_T1Smear_phi_jer1Down, &b_MET_T1Smear_phi_jer1Down);
	fChain->SetBranchAddress("MET_T1_pt_jesTotalDown", &MET_T1_pt_jesTotalDown, &b_MET_T1_pt_jesTotalDown);
	fChain->SetBranchAddress("MET_T1_phi_jesTotalDown", &MET_T1_phi_jesTotalDown, &b_MET_T1_phi_jesTotalDown);
	fChain->SetBranchAddress("MET_T1Smear_pt_jesTotalDown", &MET_T1Smear_pt_jesTotalDown, &b_MET_T1Smear_pt_jesTotalDown);
	fChain->SetBranchAddress("MET_T1Smear_phi_jesTotalDown", &MET_T1Smear_phi_jesTotalDown, &b_MET_T1Smear_phi_jesTotalDown);
	fChain->SetBranchAddress("rawmet", &rawmet, &b_rawmet);
	fChain->SetBranchAddress("rawmetphi", &rawmetphi, &b_rawmetphi);
	fChain->SetBranchAddress("metphi", &metphi, &b_metphi);
	fChain->SetBranchAddress("gen_weight", &gen_weight, &b_gen_weight);
	fChain->SetBranchAddress("npu", &npu, &b_npu);
	fChain->SetBranchAddress("ntruepu", &ntruepu, &b_ntruepu);
	fChain->SetBranchAddress("n_pos", &n_pos, &b_n_pos);
	fChain->SetBranchAddress("n_minus", &n_minus, &b_n_minus);
	fChain->SetBranchAddress("n_num", &n_num, &b_n_num);
	fChain->SetBranchAddress("MET_pass", &MET_pass, &b_MET_pass);
	fChain->SetBranchAddress("npvs", &npvs, &b_npvs);
	fChain->SetBranchAddress("btag_weight_medium", &btag_weight_medium, &b_btag_weight_medium);
	fChain->SetBranchAddress("btag_weight_medium_up", &btag_weight_medium_up, &b_btag_weight_medium_up);
	fChain->SetBranchAddress("btag_weight_medium_down", &btag_weight_medium_down, &b_btag_weight_medium_down);
	fChain->SetBranchAddress("btag_weight_loose", &btag_weight_loose, &b_btag_weight_loose);
	fChain->SetBranchAddress("btag_weight_loose_up", &btag_weight_loose_up, &b_btag_weight_loose_up);
	fChain->SetBranchAddress("btag_weight_loose_down", &btag_weight_loose_down, &b_btag_weight_loose_down); 
	fChain->SetBranchAddress("n_bjets_loose_deepcsv", &n_bjets_loose_deepcsv, &b_n_bjets_loose_deepcsv);
	fChain->SetBranchAddress("n_bjets_loose_deepFlavB", &n_bjets_loose_deepFlavB, &b_n_bjets_loose_deepFlavB);
	fChain->SetBranchAddress("n_bjets20_loose_deepcsv", &n_bjets20_loose_deepcsv, &b_n_bjets20_loose_deepcsv);
	fChain->SetBranchAddress("n_bjets20_loose_deepFlavB", &n_bjets20_loose_deepFlavB, &b_n_bjets20_loose_deepFlavB);
	fChain->SetBranchAddress("n_bjets_medium_deepcsv", &n_bjets_medium_deepcsv, &b_n_bjets_medium_deepcsv);
	fChain->SetBranchAddress("n_bjets_medium_deepFlavB", &n_bjets_medium_deepFlavB, &b_n_bjets_medium_deepFlavB);
	fChain->SetBranchAddress("n_bjets20_medium_deepcsv", &n_bjets20_medium_deepcsv, &b_n_bjets20_medium_deepcsv);
	fChain->SetBranchAddress("n_bjets20_medium_deepFlavB", &n_bjets20_medium_deepFlavB, &b_n_bjets20_medium_deepFlavB);
	fChain->SetBranchAddress("njets50", &njets50, &b_njets50);
	fChain->SetBranchAddress("njets40", &njets40, &b_njets40);
	fChain->SetBranchAddress("njets30", &njets30, &b_njets30);
	fChain->SetBranchAddress("njets20", &njets20, &b_njets20);
	fChain->SetBranchAddress("njets50_pc", &njets50_pc, &b_njets50_pc);
	fChain->SetBranchAddress("njets40_pc", &njets40_pc, &b_njets40_pc);
	fChain->SetBranchAddress("njets30_pc", &njets30_pc, &b_njets30_pc);
	fChain->SetBranchAddress("njets20_pc", &njets20_pc, &b_njets20_pc);
	fChain->SetBranchAddress("n_bjets_loose_deepcsv_pc", &n_bjets_loose_deepcsv_pc, &b_n_bjets_loose_deepcsv_pc);
	fChain->SetBranchAddress("n_bjets_loose_deepFlavB_pc", &n_bjets_loose_deepFlavB_pc, &b_n_bjets_loose_deepFlavB_pc);
	fChain->SetBranchAddress("n_bjets20_loose_deepcsv_pc", &n_bjets20_loose_deepcsv_pc, &b_n_bjets20_loose_deepcsv_pc);
	fChain->SetBranchAddress("n_bjets20_loose_deepFlavB_pc", &n_bjets20_loose_deepFlavB_pc, &b_n_bjets20_loose_deepFlavB_pc);
	fChain->SetBranchAddress("n_bjets_medium_deepcsv_pc", &n_bjets_medium_deepcsv_pc, &b_n_bjets_medium_deepcsv_pc);
	fChain->SetBranchAddress("n_bjets_medium_deepFlavB_pc", &n_bjets_medium_deepFlavB_pc, &b_n_bjets_medium_deepFlavB_pc);
	fChain->SetBranchAddress("n_bjets20_medium_deepcsv_pc", &n_bjets20_medium_deepcsv_pc, &b_n_bjets20_medium_deepcsv_pc);
	fChain->SetBranchAddress("n_bjets20_medium_deepFlavB_pc", &n_bjets20_medium_deepFlavB_pc, &b_n_bjets20_medium_deepFlavB_pc);
	fChain->SetBranchAddress("puWeight", &puWeight, &b_puWeight);
	fChain->SetBranchAddress("puWeightUp", &puWeightUp, &b_puWeightUp);
	fChain->SetBranchAddress("puWeightDown", &puWeightDown, &b_puWeightDown);
	fChain->SetBranchAddress("scalef", &scalef, &b_scalef);
	fChain->SetBranchAddress("L1PreFiringWeight_Muon_Nom", &L1PreFiringWeight_Muon_Nom, &b_L1PreFiringWeight_Muon_Nom);
	fChain->SetBranchAddress("L1PreFiringWeight_Muon_StatDn", &L1PreFiringWeight_Muon_StatDn, &b_L1PreFiringWeight_Muon_StatDn);
	fChain->SetBranchAddress("L1PreFiringWeight_Muon_StatUp", &L1PreFiringWeight_Muon_StatUp, &b_L1PreFiringWeight_Muon_StatUp);
	fChain->SetBranchAddress("L1PreFiringWeight_Muon_SystDn", &L1PreFiringWeight_Muon_SystDn, &b_L1PreFiringWeight_Muon_SystDn);
	fChain->SetBranchAddress("L1PreFiringWeight_Muon_SystUp", &L1PreFiringWeight_Muon_SystUp, &b_L1PreFiringWeight_Muon_SystUp);
	fChain->SetBranchAddress("L1PreFiringWeight_Nom", &L1PreFiringWeight_Nom, &b_L1PreFiringWeight_Nom);
	fChain->SetBranchAddress("L1PreFiringWeight_Up", &L1PreFiringWeight_Up, &b_L1PreFiringWeight_Up);
	fChain->SetBranchAddress("L1PreFiringWeight_Dn", &L1PreFiringWeight_Dn, &b_L1PreFiringWeight_Dn);
	fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
	fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
	fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
	fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
	fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
	fChain->SetBranchAddress("EcalDeadCellTriggerPrimitiveFilter", &EcalDeadCellTriggerPrimitiveFilter, &b_EcalDeadCellTriggerPrimitiveFilter);
	fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
	fChain->SetBranchAddress("Flag_BadPFMuonDzFilter", &Flag_BadPFMuonDzFilter, &b_Flag_BadPFMuonDzFilter);
	fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
	fChain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
	fChain->SetBranchAddress("LHEScaleWeight", LHEScaleWeight, &b_LHEScaleWeight);
        fChain->SetBranchAddress("LHEPdfWeight", LHEPdfWeight, &b_LHEPdfWeight);
	fChain->SetBranchAddress("PSWeight", PSWeight, &b_PSWeight);
        fChain->SetBranchAddress("ele_id_scale", &ele_id_scale, &b_ele_id_scale);
        fChain->SetBranchAddress("ele_id_scale_Up", &ele_id_scale_Up, &b_ele_id_scale_Up);
        fChain->SetBranchAddress("ele_id_scale_Down", &ele_id_scale_Down, &b_ele_id_scale_Down);
        fChain->SetBranchAddress("ele_reco_scale", &ele_reco_scale, &b_ele_reco_scale);
        fChain->SetBranchAddress("ele_reco_scale_Up", &ele_reco_scale_Up, &b_ele_reco_scale_Up);
        fChain->SetBranchAddress("ele_reco_scale_Down", &ele_reco_scale_Down, &b_ele_reco_scale_Down);
        fChain->SetBranchAddress("muon_id_scale", &muon_id_scale, &b_muon_id_scale);
        fChain->SetBranchAddress("muon_id_scale_Up", &muon_id_scale_Up, &b_muon_id_scale_Up);
        fChain->SetBranchAddress("muon_id_scale_Down", &muon_id_scale_Down, &b_muon_id_scale_Down);
        fChain->SetBranchAddress("muon_iso_scale", &muon_iso_scale, &b_muon_iso_scale);
        fChain->SetBranchAddress("muon_iso_scale_Up", &muon_iso_scale_Up, &b_muon_iso_scale_Up);
        fChain->SetBranchAddress("muon_iso_scale_Down", &muon_iso_scale_Down, &b_muon_iso_scale_Down);
        fChain->SetBranchAddress("photon_id_scale", &photon_id_scale, &b_photon_id_scale);
        fChain->SetBranchAddress("photon_id_scale_Up", &photon_id_scale_Up, &b_photon_id_scale_Up);
        fChain->SetBranchAddress("photon_id_scale_Down", &photon_id_scale_Down, &b_photon_id_scale_Down);
        fChain->SetBranchAddress("photon_veto_scale", &photon_veto_scale, &b_photon_veto_scale);
        fChain->SetBranchAddress("photon_veto_scale_Up", &photon_veto_scale_Up, &b_photon_veto_scale_Up);
        fChain->SetBranchAddress("photon_veto_scale_Down", &photon_veto_scale_Down, &b_photon_veto_scale_Down);
        fChain->SetBranchAddress("HLT_ee", &HLT_ee, &b_HLT_ee);
        fChain->SetBranchAddress("HLT_mm", &HLT_mm, &b_HLT_mm);
        fChain->SetBranchAddress("HLT_emu", &HLT_emu, &b_HLT_emu);

}

EDBRHistoMaker::EDBRHistoMaker(TTree* tree, TFile* fileTMP, TH1F* hR1, std::string out_name_, RoccoR* rc_in,TString isChannel) {
	fChain = 0;

	// Which category do we want to analyze?
	fileTMP_ = fileTMP;
	hR1_ = hR1;
	out_name = out_name_;
	rc = *rc_in;
	debug_ = true;
	Init(tree);
	createAllHistos(isChannel);
	printAllHistos();
}

float EDBRHistoMaker::get_rochester_scale(bool isdata, float charge_temp, float pt, float eta, float phi, int nl, float r1){
	int charge = int(charge_temp/fabs(charge_temp));
	// data correction
        double corr;

	if(isdata) 
           corr = rc.kScaleDT(charge, pt, eta, phi, 0, 0);

	// MC without genPT avalible
	if((!isdata))
		corr = rc.kSmearMC(charge, pt, eta, phi, nl, r1, 0, 0);
        return corr;
}
EDBRHistoMaker::~EDBRHistoMaker() {
	if (!fChain)
		return;
	delete fChain->GetCurrentFile();
}

Int_t EDBRHistoMaker::GetEntry(Long64_t entry) {
	// Read contents of entry.
	if (!fChain)
		return 0;
	return fChain->GetEntry(entry);
}

Long64_t EDBRHistoMaker::LoadTree(Long64_t entry) {
	// Set the environment to read one entry
	if (!fChain)
		return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0)
		return centry;
	if (fChain->GetTreeNumber() != fCurrent) {
		fCurrent = fChain->GetTreeNumber();
	}
	return centry;
}

//-------------------------
// Infrastructure functions
//-------------------------

void EDBRHistoMaker::createAllHistos(TString isChannel) {

	/// This part substitutes the big arrays that used to be 
	/// in the beginning of this file.
	/// Much simpler to create histos now: just add them to
	/// hs with hs.setHisto(name,nbins,min,max);
	hs.setHisto("mll", 15, 20, 300);
	hs.setHisto("photonet", 8, 20, 105);
	hs.setHisto("photoneta", 12, -1.5, 1.5);
	hs.setHisto("photonphi", 16, -3.14, 3.14);
	hs.setHisto("lep1pt", 15, 20, 180);
	hs.setHisto("lep1eta", 12, -2.4, 2.4);
	hs.setHisto("lep1phi", 16, -3.14, 3.14);
	hs.setHisto("lep2pt", 15, 20, 180);
	hs.setHisto("lep2eta", 12, -2.4, 2.4);
	hs.setHisto("lep2phi", 16, -3.14, 3.14);
        hs.setHisto("mllg", 12, 40, 500);
        hs.setHisto("ml1g", 10, 5, 200);
        hs.setHisto("ml2g", 10, 5, 300);
	hs.setHisto("ptll", 12, 30, 250);
	hs.setHisto("phiVlep", 16, -4, 4);
	hs.setHisto("yVlep", 20, -5, 5);
	hs.setHisto("MET_T1_pt", 12,20, 260);
	hs.setHisto("MET_T1_phi", 12,-3.2,3.2);
	hs.setHisto("npvs", 15, 5, 65);
	hs.setHisto("mT_pf", 20, 0, 250);
	hs.setHisto("mT2_pf",20, 0, 200);
	hs.setHisto("mT_puppi", 20, 0, 250);
	hs.setHisto("mT2_puppi", 20, 0, 200);
	hs.setHisto("PuppiMET_pt", 12,20, 260);
	hs.setHisto("PuppiMET_phi", 12,-3.2, 3.2);
        hs.setHisto("njets20_pc", 14,0, 14);
        hs.setHisto("njets30_pc", 14,0, 14);
        hs.setHisto("njets40_pc", 14,0, 14);
        hs.setHisto("njets50_pc", 14,0, 14);
	hs.setHisto("n_bjets20_deepcsv", 6, 1, 7);
	hs.setHisto("n_bjets20_deepFlavB", 6, 1, 7);
        hs.setHisto("drl1a", 20, 0, 10);
        hs.setHisto("drl2a", 20, 0, 10);
        hs.setHisto("drll", 20, 0, 10);
	char buffer[256];
	char buffer2[256];

	nVars = hs.vars.size();

	for (int i = 0; i != nVars; ++i) {
		sprintf(buffer, "%s_"+isChannel, hs.vars[i].c_str());
		sprintf(buffer2, "%s;%s;Number of events;", hs.vars[i].c_str(),hs.vars[i].c_str());
		TH1D* histogram = new TH1D(buffer, buffer2, hs.nBins[i], hs.minBin[i],hs.maxBin[i]);
		histogram->SetStats(kFALSE);
		histogram->SetDirectory(0);
		histogram->Sumw2();
		theHistograms[hs.vars[i]] = histogram;
	}

}

void EDBRHistoMaker::printAllHistos() {
	printf("We have %i histograms \n", int(theHistograms.size()));
	typedef std::map<std::string, TH1D*>::iterator it_type;
	for (it_type iterator = theHistograms.begin();
			iterator != theHistograms.end(); iterator++) {
	}
}

void EDBRHistoMaker::saveAllHistos(std::string outFileName) {

	TFile* outFile = TFile::Open(outFileName.c_str(), "RECREATE");

	for (int i = 0; i != nVars; ++i) {
		std::string name = hs.vars[i];
		const TH1D* thisHisto = this->theHistograms[name];
		thisHisto->Write();
	}
	outFile->Close();
}

///----------------------------------------------------------------
/// This is the important function, the loop over all events.
/// Here we fill the histograms according to cuts, weights,
/// and can also filter out events on an individual basis.
///----------------------------------------------------------------
void EDBRHistoMaker::Loop(std::string outFileName,float luminosity,int isBarrel,TString isChannel) {

	if (fChain == 0)
		return;
	int numbe_out = 0;
	Long64_t nentries = fChain->GetEntries();
//        nentries=100;
	std::cout <<"run loop "<< "nentries " << nentries << std::endl;
	Float_t nn;
	Long64_t npp =0;
	Long64_t nmm =0;
	//	std::cout << "numberofnp:" << npp << "  numberofnm:" << nmm << std::endl;
	Long64_t nbytes = 0, nb = 0;
	TLorentzVector Zp4, photonp4, jet1p4, jet2p4,lep1p4,lep2p4;
	float sum=0;
	for (Long64_t jentry = 0; jentry < nentries; jentry++) {
		float r1=gRandom->Rndm(jentry);
		float r2=gRandom->Rndm(jentry*2);
		drll=-1e2;
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0)
			break;

		nb = fChain->GetEntry(jentry);
		nbytes += nb;

		lep1p4.SetPtEtaPhiM(lep1pt, lep1eta, lep1phi, 0.105666);
		lep2p4.SetPtEtaPhiM(lep2pt, lep2eta, lep2phi, 0.000511);
		mll=(lep1p4+lep2p4).M();
		yVlep=(lep1p4+lep2p4).Eta();
		phiVlep=(lep1p4+lep2p4).Phi();
		ptVlep=(lep1p4+lep2p4).Pt();

		Zp4.SetPtEtaPhiM(ptVlep, yVlep, phiVlep, mll);
                photonp4.SetPtEtaPhiM(photonet, photoneta, photonphi,0);               
                mllg=(Zp4+photonp4).M();

		if (gen_weight > 0)
			nn = 1;
		else
			nn = -1;

		actualWeight = 1;//data
		if (fabs(lep1phi-lep2phi)>Pi) drll = sqrt((lep1eta-lep2eta)*(lep1eta-lep2eta)+(2*Pi-fabs(lep1phi-lep2phi))*(2*Pi-fabs(lep1phi-lep2phi)));
		else drll = sqrt((lep1eta-lep2eta)*(lep1eta-lep2eta)+(fabs(lep1phi-lep2phi))*(fabs(lep1phi-lep2phi)));

		if (setUnitaryWeights_) {
			if (jentry == 0)
				printf("Unitary weights set!\n");
			actualWeight = 1;
		}

		TString filename = fileTMP_->GetName();
		//data
		actualWeight=1;
                bool photon_channel=false,lepton_channel=false;
                Bool_t HLT=0,Filter=0;
                if(filename.Contains("MuonEG"))
			HLT=HLT_emu;
		else if(filename.Contains("Muon"))
			HLT= !(HLT_emu) && HLT_mm;
		else if(filename.Contains("Ele"))
			HLT= !(HLT_emu) && !(HLT_mm) && HLT_ee;
     
		if(filename.Contains("17")||filename.Contains("18"))
			Filter = Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter  && Flag_eeBadScFilter && Flag_ecalBadCalibFilter;      
		else
			Filter = Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_BadPFMuonDzFilter && Flag_eeBadScFilter;      
 
                if(isChannel=="emu") 
			lepton_channel= ( HLT && channel==1 && fabs(lep1_pid)==13 && fabs(lep2_pid)==11 && lep1pt>20 && lep2pt>25 && fabs(lep1eta) < 2.4 && fabs(lep1eta) < 2.5) && lep1_charge*lep2_charge>0 && drll>0.5 && Filter;

                if(isBarrel==1)  photon_channel=( (fabs(photoneta) < 1.4442) );
                else if(isBarrel==0) photon_channel= ( fabs(photoneta) < 2.5 && fabs(photoneta)>1.566 );
                else photon_channel=( (fabs(photoneta) < 1.4442) || ( fabs(photoneta) < 2.5 && fabs(photoneta)>1.566 ));

	        if( lepton_channel && n_loose_ele==1 && n_loose_mu==1 && lep1_is_tight==1 && lep2_is_tight==1 && mll > 20 && ptll > 30 && PuppiMET_pt > 20 /* MET_T1_pt > 20*/ && /*n_bjets20_loose_deepcsv_pc>=1*/ n_bjets20_medium_deepFlavB_pc == 0 && mT_puppi>0 ){
//			cout<<jentry<<" pass"<<endl;
			sum = sum + actualWeight;
			numbe_out++;
			treename->Fill();
		}
		else{
//			cout<<jentry<<" fail"<<endl;
			continue;
		}

		(theHistograms["mll"])->Fill(mll, actualWeight);
		(theHistograms["photonet"])->Fill(photonet, actualWeight);
		(theHistograms["photoneta"])->Fill(photoneta, actualWeight);
		(theHistograms["photonphi"])->Fill(photonphi, actualWeight);
		(theHistograms["lep1pt"])->Fill(lep1pt, actualWeight);
		(theHistograms["lep1eta"])->Fill(lep1eta, actualWeight);
		(theHistograms["lep1phi"])->Fill(lep1phi, actualWeight);
		(theHistograms["lep2pt"])->Fill(lep2pt, actualWeight);
		(theHistograms["lep2eta"])->Fill(lep2eta, actualWeight);
		(theHistograms["lep2phi"])->Fill(lep2phi, actualWeight);
		(theHistograms["mllg"])->Fill(mllg, actualWeight);
                (theHistograms["ml1g"])->Fill(ml1g, actualWeight);
                (theHistograms["ml2g"])->Fill(ml2g, actualWeight);
		(theHistograms["ptll"])->Fill(ptll, actualWeight);
		(theHistograms["phiVlep"])->Fill(phiVlep, actualWeight);
		(theHistograms["yVlep"])->Fill(yVlep, actualWeight);
		(theHistograms["MET_T1_pt"])->Fill(MET_T1_pt, actualWeight);
		(theHistograms["MET_T1_phi"])->Fill(MET_T1_phi, actualWeight);
		(theHistograms["npvs"])->Fill(npvs, actualWeight);
                (theHistograms["mT_pf"])->Fill(mT_pf, actualWeight);
                (theHistograms["mT2_pf"])->Fill(mT2_pf, actualWeight);
                (theHistograms["mT_puppi"])->Fill(mT_puppi, actualWeight);
                (theHistograms["mT2_puppi"])->Fill(mT2_puppi, actualWeight);
		(theHistograms["PuppiMET_pt"])->Fill(PuppiMET_pt, actualWeight);
		(theHistograms["PuppiMET_phi"])->Fill(PuppiMET_phi, actualWeight);
                (theHistograms["njets20_pc"])->Fill(njets20_pc, actualWeight);
                (theHistograms["njets30_pc"])->Fill(njets30_pc, actualWeight);
                (theHistograms["njets40_pc"])->Fill(njets40_pc, actualWeight);
                (theHistograms["njets50_pc"])->Fill(njets50_pc, actualWeight);
		(theHistograms["n_bjets20_deepcsv"])->Fill(n_bjets20_medium_deepcsv_pc, actualWeight);
		(theHistograms["n_bjets20_deepFlavB"])->Fill(n_bjets20_medium_deepFlavB_pc, actualWeight);
                (theHistograms["drl1a"])->Fill(drl1a, actualWeight);
                (theHistograms["drl2a"])->Fill(drl2a, actualWeight);
                (theHistograms["drll"])->Fill(drll, actualWeight);
	}     //end loop over entries
	cout << "after cut: " << numbe_out << "*actualweight " << actualWeight
		<< " result " << sum <<"; yields "<<sum<< endl;
	this->saveAllHistos(outFileName);
}

void EDBRHistoMaker::Loop_SFs_mc(std::string outFileName,float luminosity,int isBarrel,TString isChannel,TString year){
	if (fChain == 0)
		return;
	int numbe_out = 0;
	Long64_t nentries = fChain->GetEntriesFast();
//	nentries=100;
	std::cout << "nentries " << nentries << std::endl;
	Long64_t npp=0;
	Long64_t nmm=0;
	Float_t nn;
	Long64_t nbytes = 0, nb = 0;

	TLorentzVector Zp4, photonp4, jet1p4, jet2p4, lep1p4, lep2p4;
	float sum=0;

	TString y;
	if(year=="16") y="16post";
	else y=year;
	TFile*f1 = new TFile("/home/pku/guanz/HLT/HLT_mva_less/results/HLTSF_20"+y+".root");
	TH1D*h1 = (TH1D*)f1->Get("HLTSF_20"+y);
	TH1D*h2 = (TH1D*)f1->Get("HLTSF_UP_ERR_20"+y);
	TH1D*h3 = (TH1D*)f1->Get("HLTSF_DOWN_ERR_20"+y);

	for (Long64_t jentry = 0; jentry < nentries; jentry++) {
		float r1=gRandom->Rndm(jentry);
		float r2=gRandom->Rndm(jentry*2);
		drll=-1e2;
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0)
			break;

		nb = fChain->GetEntry(jentry);
		nbytes += nb;

		//rochester correction
		TString filename = fileTMP_->GetName();

                lep1p4.SetPtEtaPhiM(lep1pt, lep1eta, lep1phi, 0.105666);
                lep2p4.SetPtEtaPhiM(lep2pt, lep2eta, lep2phi, 0.000511);
                mll=(lep1p4+lep2p4).M();
                yVlep=(lep1p4+lep2p4).Eta();
                phiVlep=(lep1p4+lep2p4).Phi();
                ptVlep=(lep1p4+lep2p4).Pt();

                Zp4.SetPtEtaPhiM(ptVlep, yVlep, phiVlep, mll);
                photonp4.SetPtEtaPhiM(photonet, photoneta, photonphi,0);        
                mllg=(Zp4+photonp4).M();

		if(filename.Contains("plj") || filename.Contains("fake")){
                        MET_T1Smear_pt = MET_T1_pt; 
                        MET_T1Smear_phi= MET_T1_phi; 
		}

		if (jentry % 100000 == 0)
			std::cout << "Entry num " << jentry << std::endl;

		if (gen_weight > 0)
			nn = 1;
		else
			nn = -1;

		if (fabs(lep1phi-lep2phi)>Pi) drll = sqrt((lep1eta-lep2eta)*(lep1eta-lep2eta)+(2*Pi-fabs(lep1phi-lep2phi))*(2*Pi-fabs(lep1phi-lep2phi)));
		else drll = sqrt((lep1eta-lep2eta)*(lep1eta-lep2eta)+(fabs(lep1phi-lep2phi))*(fabs(lep1phi-lep2phi)));

		if (setUnitaryWeights_) {
			if (jentry == 0)
				printf("Unitary weights set!\n");
			actualWeight = 1;
		}
                photon_flag=false;lepton_flag=false;
                if(filename.Contains("plj")){//fake photon
                       if( (photon_selection==2 || photon_selection==3 || photon_selection==4 || photon_selection==5) )
			       photon_flag=1;
                       if(lep1_is_tight==1 && lep2_is_tight==1)
			       lepton_flag=1;
		}
		else if(filename.Contains("fake")){//fake lepton
                        if(photon_selection==1)
				photon_flag=1;
			if( !(lep1_is_tight==1 && lep2_is_tight==1)  )
			      lepton_flag=1;
		}
                else if(filename.Contains("WWG_emu_tot")){ //signal
                        if(photon_selection==1 && photon_isprompt==1)
                                photon_flag=1;
                        else if( (photon_selection==2 || photon_selection==3 || photon_selection==4 || photon_selection==5) && photon_isprompt==1 )
                                photon_flag=1;
                        if(lep1_is_tight==1 && lep2_is_tight==1 && lep1_isprompt==1 && lep2_isprompt==1)
                                lepton_flag=1;
                        else if( ! (lep1_is_tight==1 && lep2_is_tight==1) && lep1_isprompt==1 && lep2_isprompt==1 )
                                lepton_flag=1;
                }
		else if(filename.Contains("WJets")){
//			if( ! (lep1_is_tight==1 && lep2_is_tight==1) )
				lepton_flag=1;
		}
		else{//MC
                     if(photon_selection==1 && photon_isprompt==1) 
			     photon_flag=1;
                     if(lep1_is_tight==1 && lep2_is_tight==1 && lep1_isprompt==1 && lep2_isprompt==1) 
			     lepton_flag=1;
		}
		float l1,l2;
		if(lep1pt>lep2pt){ l1=lep1pt;l2=lep2pt;}
		else {l1=lep2pt;l2=lep1pt;}
                HLT_SF=h1->GetBinContent(h1->FindBin(l1,l2));

                if(filename.Contains("18")){L1PreFiringWeight_Nom=1;L1PreFiringWeight_Muon_Nom=1;}
                if(filename.Contains("16") )
                        L1PreFiringWeight_Muon_Nom=1;
		actualWeight = scalef*ele_id_scale*ele_reco_scale*muon_id_scale*muon_iso_scale*photon_id_scale*photon_veto_scale*puWeight*btag_weight_medium*L1PreFiringWeight_Muon_Nom*L1PreFiringWeight_Nom*HLT_SF;//mc
		if(filename.Contains("plj") || filename.Contains("fake")){
			actualWeight = scalef;luminosity=1;
		    //if (filename.Contains("fake") && lep1pt<40) actualWeight = scalef*0.90;
        }// mc

		bool photon_channel=false,lepton_channel=false;
                Bool_t HLT=0,Filter=0;
                if(filename.Contains("plj") || filename.Contains("fake")) HLT=1;
                else HLT = (HLT_ee || HLT_mm || HLT_emu);

		if(filename.Contains("17")||filename.Contains("18"))
			Filter = Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter  && Flag_eeBadScFilter && Flag_ecalBadCalibFilter;      
		else
			Filter = Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_eeBadScFilter;      

//                cout<<jentry<<" "<<Flag_goodVertices<<" "<<Flag_globalSuperTightHalo2016Filter<<" "<<Flag_HBHENoiseFilter<<" "<<Flag_HBHENoiseIsoFilter<<" "<<Flag_EcalDeadCellTriggerPrimitiveFilter<<" "<<Flag_BadPFMuonFilter<<" "<<Flag_BadPFMuonDzFilter<<" "<<Flag_eeBadScFilter<<" "<<Flag_ecalBadCalibFilter<<endl; 
                if(isChannel=="emu") 
			lepton_channel= ( HLT && channel==1 && fabs(lep1_pid)==13 && fabs(lep2_pid)==11 && lep1pt>20 && lep2pt>25 && fabs(lep1eta) < 2.4 && fabs(lep1eta) < 2.5) && lep1_charge*lep2_charge>0 && drll>0.5 && Filter ;

                if(isBarrel==1)  photon_channel=( (fabs(photoneta) < 1.4442) );
                else if(isBarrel==0) photon_channel= ( fabs(photoneta) < 2.5 && fabs(photoneta)>1.566 );
                else photon_channel=( (fabs(photoneta) < 1.4442) || ( fabs(photoneta) < 2.5 && fabs(photoneta)>1.566 ));

		if( lepton_channel && n_loose_ele==1 && n_loose_mu==1 && mll > 20 && ptll > 30 && lepton_flag && PuppiMET_pt > 20 &&/* MET_T1Smear_pt > 20 &&*/ /*n_bjets20_loose_deepcsv_pc >= 1*/ n_bjets20_medium_deepFlavB_pc == 0 && mT_puppi>0 ){
			if(gen_weight>0) npp++;
			if(gen_weight<0) nmm++;
			numbe_out++;
			treename->Fill();
		}
		else
			continue;

		sum = (sum + actualWeight);

		(theHistograms["mll"])->Fill(mll, actualWeight);
		(theHistograms["photonet"])->Fill(photonet, actualWeight);
		(theHistograms["photoneta"])->Fill(photoneta, actualWeight);
		(theHistograms["photonphi"])->Fill(photonphi, actualWeight);
		(theHistograms["lep1pt"])->Fill(lep1pt, actualWeight);
		(theHistograms["lep1eta"])->Fill(lep1eta, actualWeight);
		(theHistograms["lep1phi"])->Fill(lep1phi, actualWeight);
		(theHistograms["lep2pt"])->Fill(lep2pt, actualWeight);
		(theHistograms["lep2eta"])->Fill(lep2eta, actualWeight);
		(theHistograms["lep2phi"])->Fill(lep2phi, actualWeight);
		(theHistograms["mllg"])->Fill(mllg, actualWeight);
                (theHistograms["ml1g"])->Fill(ml1g, actualWeight);
                (theHistograms["ml2g"])->Fill(ml2g, actualWeight);
		(theHistograms["ptll"])->Fill(ptll, actualWeight);
		(theHistograms["phiVlep"])->Fill(phiVlep, actualWeight);
		(theHistograms["yVlep"])->Fill(yVlep, actualWeight);
		(theHistograms["MET_T1_pt"])->Fill(MET_T1Smear_pt, actualWeight);
		(theHistograms["MET_T1_phi"])->Fill(MET_T1Smear_phi, actualWeight);
		(theHistograms["npvs"])->Fill(npvs, actualWeight);
                (theHistograms["mT_pf"])->Fill(mT_pf, actualWeight);
                (theHistograms["mT2_pf"])->Fill(mT2_pf, actualWeight);
                (theHistograms["mT_puppi"])->Fill(mT_puppi, actualWeight);
                (theHistograms["mT2_puppi"])->Fill(mT2_puppi, actualWeight);
		(theHistograms["PuppiMET_pt"])->Fill(PuppiMET_pt, actualWeight);
		(theHistograms["PuppiMET_phi"])->Fill(PuppiMET_phi, actualWeight);
                (theHistograms["njets20_pc"])->Fill(njets20_pc, actualWeight);
                (theHistograms["njets30_pc"])->Fill(njets30_pc, actualWeight);
                (theHistograms["njets40_pc"])->Fill(njets40_pc, actualWeight);
                (theHistograms["njets50_pc"])->Fill(njets50_pc, actualWeight);
		(theHistograms["n_bjets20_deepcsv"])->Fill(n_bjets20_medium_deepcsv_pc, actualWeight);
		(theHistograms["n_bjets20_deepFlavB"])->Fill(n_bjets20_medium_deepFlavB_pc, actualWeight);
                (theHistograms["drl1a"])->Fill(drl1a, actualWeight);
                (theHistograms["drl2a"])->Fill(drl2a, actualWeight);
                (theHistograms["drll"])->Fill(drll, actualWeight);


	}
	f1->Close();
	cout << "after cut: " << numbe_out << "; actualweight" << actualWeight<<endl;
	cout<< " total events: " << sum <<"; yields "<<sum*luminosity<<endl;
	this->saveAllHistos(outFileName);
}
