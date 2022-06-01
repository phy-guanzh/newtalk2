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
#include "EDBRHistoPlotter.h"
#include "CMSTDRStyle.h"
#include "roccor/RoccoR.cc"
#include "TH2.h"
void loopPlot(int isBarrel,TString isChannel,std::string year) {
	gErrorIgnoreLevel = kFatal; //suppresses all info messages

	setTDRStyle(); //TDR style
	//#####################EDIT THE OPTIONS##############################

	TString channel,lep_ch;
	if(isBarrel==1)channel="b";else if(isBarrel==0)channel="e"; else channel="a";
	double lumiValue1 ;
	if(year=="18")
		lumiValue1 = 59.7;
	else if(year=="17")
		lumiValue1 = 41.5;
        else if(year=="16"){lumiValue1=16.81;}
        else {lumiValue1=19.52;}
	double lumiValue;
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
        TString y;
        if(year=="16pre") y="16";
        else y=year;
//	vector<TString> pathToTrees = {"/home/pku/anying/cms/rootfiles/WWg/20"+y+"/"};
	vector<TString> pathToTrees = {"../Common/rootfiles/"};
	std::string outputDir = "./fig-output_"+year+"/";
	RoccoR  rc;
        /// file for scale factors
        std::vector<std::string> bkgLegend; 

// Setup names of data files for trees.
	const int nDATA = 3;
	std::cout << "set data imformation, we have " << nDATA << "data file"
			<< std::endl;
	std::string data[nDATA] = { "MuonEG","Muon","Ele"};
	std::vector < TString > fData;
	for (int j = 0; j< pathToTrees.size(); j++) {
		for (int ii = 0; ii < nDATA; ii++) {
			fData.push_back(pathToTrees[j] +"optimal_emua_"+ data[ii] +year+ ".root");
		}
	}
	const int kk=fData.size();
	std::string dataLabels[kk] = { "MuonEG","Muon","Ele"};
	cout<<"Data size "<<fData.size()<<endl;
	// set mc imformation
	const int nMC = 7;//
	std::cout << "set data imformation, we have " << nMC << "mc file"
		<< std::endl;
        std::string mc[nMC] ={"fakeL","WW","VV","ST","DYJets","TTJets","tZq"};
	std::vector< TString > fMC;
	for (int ii = 0; ii < nMC; ii++) {
		fMC.push_back(pathToTrees[0] +"optimal_emua_"+ mc[ii] +year+ ".root");
		cout<<pathToTrees[0]<<"optimal_emua_"<<mc[ii] <<year<<".root"<<endl;
	}
	const int nmc=fMC.size();
	//std::string bkgLabels[nmc]={"fakeL","VV","tW","DY","t#bar{t}","tZq"};
	std::string bkgLabels[nmc] ={"fakeL","WW","VV","tW","DY","t#bar{t}","tZq"};
	//double kFactorsMC_array[nmc] ={1,lumiValue1,lumiValue1,lumiValue1,lumiValue1,lumiValue1,lumiValue1};
	double kFactorsMC_array[nmc] ={1,lumiValue1,lumiValue1,lumiValue1,lumiValue1,lumiValue1,lumiValue1};

	cout<<"MC size "<<fMC.size()<<endl;
	std::vector<double> kFactorsMC;
	for (int index = 0; index < nmc; index++) {
                bkgLegend.push_back(bkgLabels[index]);
		kFactorsMC.push_back(kFactorsMC_array[index]);
	}
	// set mcsig information
	const int nMCSig = 1;
	std::cout << "set data imformation, we have " << nMCSig << "mcsig file"
		<< std::endl;
	std::string mcSig[nMCSig] = { "WWG_emu_tot"};
	std::vector < TString > fMCSig;
	for (int ii = 0; ii < nMCSig; ii++) {
		fMCSig.push_back(pathToTrees[0] + "optimal_emua_" + mcSig[ii] + year+".root");
	}
	cout<<"MC Sig size "<<fMCSig.size()<<endl;
	const int nmcsig=fMCSig.size();
	std::string mcLabelsSig[nmcsig] = { "WWG_emu_tot"};
	double kFactorsSig_array[nmcsig] = {lumiValue1};
	std::vector<double> kFactorsMCSig;
	for (int index = 0; index < nmcsig; index++) {
		kFactorsMCSig.push_back(kFactorsSig_array[index]);
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
		if(fData.at(i).Contains("17")) {tag="17";rc.init("roccor/RoccoR2017UL.txt");lumiValue=41.52; }
		if(fData.at(i).Contains("18")) {tag="18";rc.init("roccor/RoccoR2018UL.txt");lumiValue=59.7; }
		sprintf(buffer, "./output-slimmed-rootfiles/histos_"+isChannel+channel+"_%s.root",dataLabels[i].c_str());
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
		std::cout << "\n-------\nRunning over " << mc[i].c_str()
			<< std::endl;
		std::cout << "The file is " << fMC.at(i) << std::endl;
		TString tag;
                if(fMC.at(i).Contains("16")) {tag="16";rc.init("roccor/RoccoR2016bUL.txt"); lumiValue=16.8;}
                if(fMC.at(i).Contains("pre")){tag="16pre";rc.init("roccor/RoccoR2016aUL.txt"); lumiValue=19.5;}
		if(fMC.at(i).Contains("17")) {tag="17";rc.init("roccor/RoccoR2017UL.txt"); lumiValue=41.52;}
		if(fMC.at(i).Contains("18")) {tag="18";rc.init("roccor/RoccoR2018UL.txt"); lumiValue=59.7; }
		sprintf(buffer, "./output-slimmed-rootfiles/histos_"+isChannel+channel+"_%s.root", mc[i].c_str());
		sprintf(out_buffer, "./output-slimmed-rootfiles/optimal_"+isChannel+channel+"_%s"+tag+".root", mc[i].c_str());
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
                if(fMCSig.at(i).Contains("16")) {tag="16";rc.init("roccor/RoccoR2016bUL.txt"); lumiValue=16.8;}
                if(fMCSig.at(i).Contains("pre")){tag="16pre";rc.init("roccor/RoccoR2016aUL.txt"); lumiValue=19.5;}
		if(fMCSig.at(i).Contains("17")) {tag="17";rc.init("roccor/RoccoR2017UL.txt"); lumiValue=41.52;}
		if(fMCSig.at(i).Contains("18")) {tag="18";rc.init("roccor/RoccoR2018UL.txt"); lumiValue=59.7; }
		sprintf(buffer, "./output-slimmed-rootfiles/histos_"+isChannel+channel+"_%s.root",mcLabelsSig[i].c_str());
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
		sprintf(buffer, "./output-slimmed-rootfiles/histos_"+isChannel+channel+"_%s.root", mc[0].c_str());
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
			fHistosMC, fHistosMCSig, lumiValue, scaleToData, scaleOnlyWJets,
			makeRatio, isSignalStackOnBkg, kFactorsMC, kFactorsMCSig);
	std::cout << "Set output dir" << std::endl;
	plotter->setOutDir(outputDir,isBarrel);
	plotter->setDebug(false);

	//colors are assigned in the same order of mcLabels
	// For ZZ
	////// {DYJetsToLL_HT-200to400,DYJetsToLL_HT-200to400,DYJetsToLL_HT-600toInf}
	std::vector<int> fColorsMC;

	fColorsMC.push_back(kGreen+2);
	fColorsMC.push_back(kViolet);
	fColorsMC.push_back(kCyan);
	fColorsMC.push_back(40);
	fColorsMC.push_back(kYellow-7);
	fColorsMC.push_back(kBlue-6);
	fColorsMC.push_back(kBlue);
	fColorsMC.push_back(kGreen);
	fColorsMC.push_back(kGreen-6);
	fColorsMC.push_back(kGreen-7);
	fColorsMC.push_back(kGreen-9);
	fColorsMC.push_back(kGreen-10);


	std::vector<int> fColorsMCSig;
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
//	loopPlot(1,"mm");
///	loopPlot(1,"emu");
	loopPlot(2,"emu","18");
	loopPlot(2,"emu","17");
	loopPlot(2,"emu","16");
	loopPlot(2,"emu","16pre");
	return 0;
}

