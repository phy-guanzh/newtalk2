#include "TGraphAsymmErrors.h"
#include "Riostream.h"
#include <vector>
#include <string>
#include "TROOT.h"
#include "TError.h"
#include "TFile.h"
#include "TCollection.h"
#include "TKey.h"
#include "EDBRHistoMaker.h"
#include "EDBRHistoPlotter_all.h"
#include "CMSTDRStyle.h"
#include "roccor/RoccoR.cc"
#include "TH2.h"
void loopPlot(int isBarrel,TString isChannel) {
	gErrorIgnoreLevel = kFatal; //suppresses all info messages

	setTDRStyle(); //TDR style
	//#####################EDIT THE OPTIONS##############################

	TString channel,lep_ch;
	if(isBarrel==1)channel="b";else if(isBarrel==0)channel="e"; else channel="a";
	double lumiValue1,lumiValue2,lumiValue3,lumiValue4,lumiValue; 
        lumiValue4=59.7;
	lumiValue3=41.5;
	lumiValue2=16.81;
	lumiValue1=19.52;
	/// Should we scale the histograms to data?
	bool scaleToData = false;
	// Should we scale only wjets to make total MC = DATA?
	bool scaleOnlyWJets = false;
	/// Should we plot the Data/Bkg and Data-Bkg/Error ratios?
	bool makeRatio = true;
	/// Should we REDO histograms?
	bool redoHistograms = true;
	/// Should we put the signal MC stacked on top of the background (or just plot the signal alone)?
	bool isSignalStackOnBkg = false;
	bool dopileupreweight = false;

	// Path to wherever the files with the trees are. 
	vector<TString> pathToTrees = {"../Common/rootfiles/"};
	std::string outputDir = "./fig-output_16all/";
	RoccoR  rc;
        /// file for scale factors
        std::vector<std::string> bkgLegend; 

// Setup names of data files for trees.
	const int nDATA = 3;
	std::cout << "set data imformation, we have " << nDATA << "data file"
			<< std::endl;
	std::string data[nDATA] = { "MuonEG","Muon","Ele"};
	std::vector < TString > fData;
	std::vector<std::string> dataLabels;
	for (int ii = 0; ii < nDATA; ii++) {
		for (int j = 0; j< pathToTrees.size(); j++) {
			if(j==0){
				dataLabels.push_back(data[ii]);
				dataLabels.push_back(data[ii]);
				fData.push_back(pathToTrees[j] +"optimal_emua_"+ data[ii] + "16.root");
				fData.push_back(pathToTrees[j] +"optimal_emua_"+ data[ii] + "16pre.root");
			}
			else if(j==1){ 
				dataLabels.push_back(data[ii]);
				fData.push_back(pathToTrees[j] +"optimal_emua_"+ data[ii] + "17.root");
			}
			else if(j==2){ 
				dataLabels.push_back(data[ii]);
				fData.push_back(pathToTrees[j] +"optimal_emua_"+ data[ii] + "18.root");
			}
		}
	}
	cout<<"Data size "<<fData.size()<<endl;
	// set mc imformation
	const int nMC = 7;//
	std::cout << "set data imformation, we have " << nMC << "mc file"
		<< std::endl;
        std::string mc[nMC] ={"fakeL","WW","VV","ST","DYJets","TTJets","tZq"};
	std::vector< TString > fMC;
	std::vector<std::string> mcLabels;
	for (int ii = 0; ii < nMC; ii++) {
		for (int j = 0; j< pathToTrees.size(); j++) {
			if(j==0){
				mcLabels.push_back(mc[ii]);
				mcLabels.push_back(mc[ii]);
				fMC.push_back(pathToTrees[j] +"optimal_emua_"+ mc[ii] + "16.root");
				fMC.push_back(pathToTrees[j] +"optimal_emua_"+ mc[ii] + "16pre.root");
			}
			if(j==1){
				mcLabels.push_back(mc[ii]);
				fMC.push_back(pathToTrees[j] +"optimal_emua_"+ mc[ii] + "17.root");
			}
			if(j==2){
				mcLabels.push_back(mc[ii]);
				fMC.push_back(pathToTrees[j] +"optimal_emua_"+ mc[ii] + "18.root");
			}
			cout<<pathToTrees[j]<<"optimal_emua_"<<mc[ii]<<".root"<<endl;
		}
	}
	const int nmc=nMC;
	std::string bkgLabels[nmc] ={"fakeL","WW","VV","tW","DY","t#bar{t}","tZq"};

	cout<<"MC size "<<fMC.size()<<endl;
	std::vector<double> kFactorsMC;
	for (int index = 0; index < nmc; index++) {
                bkgLegend.push_back(bkgLabels[index]);
	}
	for(int i=0;i<fMC.size();i++){
		if(fMC[i].Contains("fake") || fMC[i].Contains("plj"))
			kFactorsMC.push_back(1);
		else{
			if(fMC[i].Contains("18")) kFactorsMC.push_back(59.7);
			else if(fMC[i].Contains("17")) kFactorsMC.push_back(41.52);
			else if(fMC[i].Contains("16pre")) kFactorsMC.push_back(19.52);
			else if(fMC[i].Contains("16")) kFactorsMC.push_back(16.81);
		}
	}
	// set mcsig information
	const int nMCSig = 1;
	std::cout << "set data imformation, we have " << nMCSig << "mcsig file"
		<< std::endl;
	std::string mcSig[nMCSig] = { "WWG_emu_tot"};
	std::vector < TString > fMCSig;
	std::vector<std::string> mcLabelsSig;
	for (int ii = 0; ii < nMCSig; ii++) {
		for (int j = 0; j< pathToTrees.size(); j++) {
			if(j==0){ 
				mcLabelsSig.push_back(mcSig[ii]);
				mcLabelsSig.push_back(mcSig[ii]);
				fMCSig.push_back(pathToTrees[j] + "optimal_emua_" + mcSig[ii] +"16.root");
				fMCSig.push_back(pathToTrees[j] + "optimal_emua_" + mcSig[ii] +"16pre.root");
			}
			else if(j==1){
				mcLabelsSig.push_back(mcSig[ii]);
				fMCSig.push_back(pathToTrees[j] + "optimal_emua_" + mcSig[ii] + "17.root");
			}
			else if(j==2){
				mcLabelsSig.push_back(mcSig[ii]);
				fMCSig.push_back(pathToTrees[j] + "optimal_emua_" + mcSig[ii] + "18.root");
			}
		}
	}
	cout<<"MC Sig size "<<fMCSig.size()<<endl;

	std::vector<double> kFactorsMCSig;
	for (int i = 0; i < fMCSig.size(); i++) {
		if(fMCSig[i].Contains("18")) kFactorsMCSig.push_back(59.7);
		else if(fMCSig[i].Contains("17")) kFactorsMCSig.push_back(41.52);
		else if(fMCSig[i].Contains("16pre")) kFactorsMCSig.push_back(19.52);
		else if(fMCSig[i].Contains("16")) kFactorsMCSig.push_back(16.81);
	}

	std::cout << "set data mc mcsig Jet information done" << std::endl;

	std::vector < std::string > fHistosData;
	std::vector < std::string > fHistosMC;
	std::vector < std::string > fHistosMCSig;
        std::vector < TString > histName={"m_{#font[12]{ll}}","p_{T}^{#gamma}","#eta_{#gamma}","#phi_{#gamma}","p_{T}^{#font[12]{l}_{#mu}}","#eta_{#font[12]{l}_{#mu}}","#phi_{#font[12]{l}_{#mu}}","p_{T}^{#font[12]{l}_{e}}","#eta_{#font[12]{l}_{e}}","#phi_{#font[12]{l}_{e}}","m_{#font[12]{ll}#gamma}","m_{#font[12]{l_{1}}#gamma}","m_{#font[12]{l_{2}}#gamma}","p_{T}^{#font[12]{ll}}","#phi_{#font[12]{ll}}","#eta_{#font[12]{ll}}","PFMET","#phi_{PFMET}","npvs","m^{PF}_{T_{WW}}","m^{PF}_{T_{W_{2}}}","m^{Puppi}_{T_{WW}}","m^{Puppi}_{T_{W_{2}}}","PuppiMET","#phi_{PuppiMET}","N_{jets20}","N_{jets30}","N_{jets40}","N_{jets50}","N_{bjets20_{DeepCSV}}","N_{bjets20_{DeepFlavB}}","#Delta R(#font[12]{l}_{1},#gamma)","#Delta R(#font[12]{l}_{2},#gamma)","#Delta R(#font[12]{ll})"};
	char buffer[256], out_buffer[256];
	printf("All strings set\n");

	/// ----------------------------------------------------------------
	/// This first part is the loop over trees to create histogram files 
	/// ----------------------------------------------------------------


	printf("\nStart making histograms\n\n");

	//loop over data files and make histograms individually for each of them
	TH1F* hisRatio = 0;

	for (int i = 0; i < fData.size(); i++) {
//		continue;
		std::cout << "\n-------\nRunning over " << dataLabels[i].c_str()
			<< std::endl;
		std::cout << "The file is " << fData.at(i) << std::endl; 
		TString tag,sample;
		if(fData.at(i).Contains("16")) {tag="16";rc.init("roccor/RoccoR2016bUL.txt"); lumiValue=16.8;}
		if(fData.at(i).Contains("pre")){tag="16pre";rc.init("roccor/RoccoR2016aUL.txt"); lumiValue=19.5;}
		if(fData.at(i).Contains("17")) {tag="17";rc.init("roccor/RoccoR2017UL.txt");lumiValue=41.5; }
		if(fData.at(i).Contains("18")) {tag="18";rc.init("roccor/RoccoR2018UL.txt");lumiValue=59.7; }
		sprintf(buffer, "./output-slimmed-rootfiles/histos_"+isChannel+channel+"_%s"+tag+".root",dataLabels[i].c_str());
		sprintf(out_buffer, "./output-slimmed-rootfiles/optimal_"+isChannel+channel+"_%s"+tag+".root",dataLabels[i].c_str());
		fHistosData.push_back(buffer);

		std::cout << "retrieve "<<i+1<<"th data file" << std::endl;
		TFile *fileData = TFile::Open(fData.at(i));
		TTree *treeData = (TTree*) fileData->Get("outtree");
		std::cout << "retrieve ith mc file" << std::endl;
		if (redoHistograms) {
			EDBRHistoMaker* maker = new EDBRHistoMaker(treeData, fileData,
					hisRatio, out_buffer, &rc,isChannel);
			maker->setUnitaryWeights(true);
			maker->Loop(buffer,lumiValue,isBarrel,isChannel);
			maker->endjob();
			fileData->Close();
		}

	}  //end loop on data files
	printf("Loop over data done\n");



	//loop over MC files and make histograms individually for each of them
	for (int i = 0; i < fMC.size(); i++) {
//		continue;
		std::cout << "\n-------\nRunning over " << mcLabels[i].c_str()
			<< std::endl;
		std::cout << "The file is " << fMC.at(i) << std::endl;
		TString tag;
		if(fMC.at(i).Contains("16")) {tag="16";rc.init("roccor/RoccoR2016bUL.txt"); lumiValue=16.8;}
		if(fMC.at(i).Contains("pre")){tag="16pre";rc.init("roccor/RoccoR2016aUL.txt"); lumiValue=19.5;}
		if(fMC.at(i).Contains("17")) {tag="17";rc.init("roccor/RoccoR2017UL.txt"); lumiValue=41.5;}
		if(fMC.at(i).Contains("18")) {tag="18";rc.init("roccor/RoccoR2018UL.txt"); lumiValue=59.7; }
		sprintf(buffer, "./output-slimmed-rootfiles/histos_"+isChannel+channel+"_%s"+tag+".root", mcLabels[i].c_str());
		sprintf(out_buffer, "./output-slimmed-rootfiles/optimal_"+isChannel+channel+"_%s"+tag+".root", mcLabels[i].c_str());
		fHistosMC.push_back(buffer);
		std::cout << "test" << std::endl;

		if (redoHistograms) {
			std::cout << "retrieve ith mc file" << std::endl;
			TFile *fileMC = TFile::Open(fMC.at(i));
			std::cout << "retrieve tree of mc file" << std::endl;
			TTree *treeMC = (TTree*) fileMC->Get("outtree");
			EDBRHistoMaker* maker = new EDBRHistoMaker(treeMC, fileMC,
					hisRatio, out_buffer, &rc,isChannel);
			maker->setUnitaryWeights(false);
			maker->Loop_SFs_mc(buffer,lumiValue,isBarrel,isChannel,tag);
			maker->endjob();
			fileMC->Close();
		}

	}  //end loop on MC files

	printf("Loop over MC done\n");

	//loop over MC signal files and make histograms individually for each of them
	for (int  i = 0; i < fMCSig.size(); i++) {
//		continue;
		std::cout << "\n-------\nRunning over " << mcLabelsSig[i].c_str()
			<< std::endl;
		std::cout << "The file is " << fMCSig.at(i) << std::endl;
		TString tag;
		if(fMCSig.at(i).Contains("16")) {tag="16";rc.init("roccor/RoccoR2016bUL.txt"); lumiValue= 16.8;}
		if(fMCSig.at(i).Contains("pre")){tag="16pre";rc.init("roccor/RoccoR2016aUL.txt"); lumiValue=19.5;}
		if(fMCSig.at(i).Contains("17")) {tag="17";rc.init("roccor/RoccoR2017UL.txt"); lumiValue=41.5;}
		if(fMCSig.at(i).Contains("18")) {tag="18";rc.init("roccor/RoccoR2018UL.txt"); lumiValue=59.7; }
		sprintf(buffer, "./output-slimmed-rootfiles/histos_"+isChannel+channel+"_%s"+tag+".root",mcLabelsSig[i].c_str());
		sprintf(out_buffer, "./output-slimmed-rootfiles/optimal_"+isChannel+channel+"_%s"+tag+".root",mcLabelsSig[i].c_str());
		fHistosMCSig.push_back(buffer);

		if (redoHistograms) {
			std::cout << "retrieve ith mcsig file" << std::endl;
			TFile *fileMCSig = TFile::Open(fMCSig.at(i));
			std::cout << "retrieve tree of mcsig file" << std::endl;
			TTree *treeMCSig = (TTree*) fileMCSig->Get("outtree");
			std::cout<<"open MC sig file OK"<<endl;
			EDBRHistoMaker* maker = new EDBRHistoMaker(treeMCSig, fileMCSig,
					hisRatio, out_buffer,&rc,isChannel);
			maker->setUnitaryWeights(false);
			maker->Loop_SFs_mc(buffer,lumiValue,isBarrel,isChannel,tag);
			maker->endjob();
			fileMCSig->Close();
		}
	}  //end loop on MC files

	// ------------------------------------------------------------------
	// This second part is the loop over histograms to create stack plots
	// ------------------------------------------------------------------  
	printf("\nStart looping over histograms\n\n");
	//make nice plots
	std::vector < std::string > listOfHistos;
	if (nMC > 0) {
		// Open one of the histogram files just to get the list of histograms
		// produced, then loop over all the histograms inheriting 
		// from TH1 contained in the file.
		TString tag;
		if(fMCSig.at(0).Contains("16")) {tag="16";rc.init("roccor/RoccoR2016bUL.txt"); lumiValue=16.8;}
		if(fMCSig.at(0).Contains("pre")){tag="16pre";rc.init("roccor/RoccoR2016aUL.txt"); lumiValue=19.5;}
		if(fMCSig.at(0).Contains("17")) {tag="17";rc.init("roccor/RoccoR2017UL.txt"); lumiValue=41.5;}
		if(fMCSig.at(0).Contains("18")) {tag="18";rc.init("roccor/RoccoR2018UL.txt"); lumiValue=59.7; }
		sprintf(buffer, "./output-slimmed-rootfiles/histos_"+isChannel+channel+"_%s"+tag+".root", mcLabels[0].c_str());
		std::cout << "Opening " << buffer << std::endl;
		TFile* oneFile = TFile::Open(buffer);
		TIter next(oneFile->GetListOfKeys());
		TKey *key;
		while ((key = (TKey*) next())) {
			TClass *cl = gROOT->GetClass(key->GetClassName());
			if (!cl->InheritsFrom("TH1"))
				continue;
			TH1 *hTMP = (TH1*) key->ReadObj();
			std::string hName = hTMP->GetName();
			printf("Histogram found: %s\n", hName.c_str());
			listOfHistos.push_back(hName);
		}      //end while loop
		oneFile->Close();
	}      //end if fmc size >0

	std::cout << "Creating plotter" << std::endl;
	EDBRHistoPlotter *plotter = new EDBRHistoPlotter("./", fHistosData,
			fHistosMC, fHistosMCSig, 35.9, scaleToData, scaleOnlyWJets,
			makeRatio, isSignalStackOnBkg, kFactorsMC, kFactorsMCSig);
	std::cout << "Set output dir" << std::endl;
	plotter->setOutDir(outputDir,isBarrel);
	plotter->setDebug(false);

	//colors are assigned in the same order of mcLabels
	// For ZZ
	////// {DYJetsToLL_HT-200to400,DYJetsToLL_HT-200to400,DYJetsToLL_HT-600toInf}
	std::vector<int> fColorsMC;

	fColorsMC.push_back(kGreen+2);
	fColorsMC.push_back(kGreen+2);
	fColorsMC.push_back(kGreen-3);
	fColorsMC.push_back(kGreen-3);
	fColorsMC.push_back(kCyan);
	fColorsMC.push_back(kCyan);
	fColorsMC.push_back(40);
	fColorsMC.push_back(40);
	fColorsMC.push_back(kYellow-7);
	fColorsMC.push_back(kYellow-7);
	fColorsMC.push_back(kBlue-6);
	fColorsMC.push_back(kBlue-6);
	fColorsMC.push_back(kBlue);
	fColorsMC.push_back(kBlue);
	fColorsMC.push_back(kGreen);
	fColorsMC.push_back(kGreen);


	std::vector<int> fColorsMCSig;
	fColorsMCSig.push_back(kRed-7);
	fColorsMCSig.push_back(kRed-7);

	plotter->setFillColor(fColorsMC);
	plotter->setLineColor(fColorsMCSig);

	int numOfHistos = listOfHistos.size();
	for (int i = 0; i != numOfHistos; ++i){
		plotter->makeStackPlots(listOfHistos.at(i),histName[i],isBarrel,bkgLegend);
	}
	printf("Plotting done\n");
	delete plotter;
}

int main() {
//	loopPlot(2,"emu","17");
//	loopPlot(2,"emu","18");
	loopPlot(2,"emu");
//	loopPlot(2,"emu","16pre");
	return 0;
}

