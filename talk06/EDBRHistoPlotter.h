#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TVectorD.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TLatex.h" 
#include "CMS_lumi.C"

class EDBRHistoPlotter {
public:

	EDBRHistoPlotter(std::string nameInDir,
			std::vector<std::string> nameFileDATA,
			std::vector<std::string> nameFileMC,
			std::vector<std::string> nameFileMCSig, double targetLumi,
			bool scaleToData, bool scaleOnlyWJets, bool makeRatio, bool isSignalStackOnBkg,
			std::vector<double> kFactorsMC, std::vector<double> kFactorsMCSig) {
		std::cout << "Plotter constructor" << std::endl;
		nameInDir_ = nameInDir;
		fileNamesMC = nameFileMC;
		fileNamesMCSig = nameFileMCSig;
		fileNamesDATA = nameFileDATA;
		kFactorsMC_ = kFactorsMC;
		kFactorsSig_ = kFactorsMCSig;
		targetLumi_ = targetLumi;
		scaleToData_ = scaleToData;
		scaleOnlyWJets_ = scaleOnlyWJets;
		makeRatio_ = makeRatio;
		mcTotalLabels_=mcTotalLabels;
		isSignalStackOnBkg_ = isSignalStackOnBkg;
		debug_ = true;
		if (fileNamesDATA.size() != 0)
			isDataPresent_ = true;
		else
			isDataPresent_ = false;
		std::cout << "Check" << std::endl;
		EDBRColors.resize(20, kWhite);
		EDBRLineColors.resize(20, kWhite);
		std::cout << "Check" << std::endl;
		labels.resize(0);
		labelsSig.resize(0);
//		mcTotalLabels.resize(0);
		std::cout << "Check" << std::endl;
		makeLabels();
		std::cout << "Check" << std::endl;
		if (fileNamesMCSig.size() != kFactorsSig_.size()) {
			cout << "======> Mismatch in size of input MC Sig arrays !!! "
					<< fileNamesMCSig.size() << "  " << kFactorsSig_.size()
					<< endl;
		}

		printf("Target lumi is %g fb-1\n", targetLumi);
		std::cout << "k factors for MC backgrounds are: " << std::endl;
		int myKindex = 0;
		for (std::vector<double>::iterator it = kFactorsMC_.begin();
				it != kFactorsMC_.end(); ++it) {
			std::cout << *it << " for " << fileNamesMC.at(myKindex)
					<< std::endl;
			myKindex++;
		}
		myKindex = 0;
		for (std::vector<double>::iterator it = kFactorsSig_.begin();
				it != kFactorsSig_.end(); ++it) {
			std::cout << *it << " for " << fileNamesMCSig.at(myKindex)
					<< std::endl;
			myKindex++;
		}

		std::cout << std::endl;
	} //end constructor

	virtual ~EDBRHistoPlotter() {
	} //end destructor

	/// Members
	std::vector<std::string> fileNamesMC;
	std::vector<std::string> fileNamesMCSig;
	std::vector<std::string> fileNamesDATA;
	std::vector<std::string> labels;
	std::vector<std::string> labelsSig;
	std::vector<std::string> mcTotalLabels;
	std::vector<TFile*> filesMC;
	std::vector<TFile*> filesMCSig;
	std::vector<TFile*> filesDATA;
	std::vector<TH1D*> histosMC;
	std::vector<TH1D*> histosMCSig;
	std::vector<TH1D*> histosMCSigOrig;
	std::vector<TH1D*> histosDATA;
	std::vector<int> EDBRColors;
	std::vector<int> EDBRLineColors;
	std::vector<double> kFactorsMC_;
	std::vector<double> kFactorsSig_;

	std::string nameInDir_;
	std::string nameOutDir_;
	std::vector<std::string> mcTotalLabels_;
	double dataIntegral_;
	double targetLumi_;
	bool scaleToData_;
	bool scaleOnlyWJets_;
	bool makeRatio_;
	bool isSignalStackOnBkg_;
	bool isDataPresent_;
	bool debug_;

	/// Functions
	void cleanupMC();
	void cleanupMCSig();
	void cleanupDATA();
	void makeLabels();
	void makeStackPlots(std::string histoName,TString histTitle,int isBarrel,vector<std::string> bkgLegend);
	void setOutDir(std::string outDirNew,int isBarrel);

	/// set debug mode
	void setDebug(bool debug) {
		debug_ = debug;
	}

	/// get reasonable colors for stacks.
	int getFillColor(int index) {
		if (index < 20)
			return EDBRColors[index];
		return kWhite;
	}

	/// set reasonable colors for stacks.
	void setFillColor(std::vector<int> colorList) {
		unsigned int ind = 0;
		while (ind < 20 && ind < colorList.size()) {
			EDBRColors.at(ind) = colorList.at(ind);
			ind++;
		}
	}

	/// get reasonable colors for stacks.
	int getLineColor(int index) {
		if (index < 20)
			return EDBRLineColors[index];
		return kWhite;
	}

	/// set reasonable colors for stacks.
	void setLineColor(std::vector<int> colorList) {
		unsigned int ind = 0;
		while (ind < 20 && ind < colorList.size()) {
			EDBRLineColors.at(ind) = colorList.at(ind);
			ind++;
		}
	}

};

void EDBRHistoPlotter::cleanupMC() {
	for (size_t i = 0; i != filesMC.size(); ++i) {
		filesMC.at(i)->Close();
	}
	filesMC.clear();

	for (size_t i = 0; i != histosMC.size(); ++i) {
		histosMC.at(i)->Delete();
	}
	histosMC.clear();
}

void EDBRHistoPlotter::cleanupMCSig() {
	for (size_t i = 0; i != filesMCSig.size(); ++i) {
		filesMCSig.at(i)->Close();
	}
	filesMCSig.clear();

	for (size_t i = 0; i != histosMCSig.size(); ++i) {
		histosMCSig.at(i)->Delete();
		histosMCSigOrig.at(i)->Delete();
	}
	histosMCSig.clear();
	histosMCSigOrig.clear();
}

void EDBRHistoPlotter::cleanupDATA() {
	for (size_t i = 0; i != filesDATA.size(); ++i) {
		filesDATA.at(i)->Close();
	}
	filesDATA.clear();

	for (size_t i = 0; i != histosDATA.size(); ++i) {
		histosDATA.at(i)->Delete();
	}
	histosDATA.clear();
}

void EDBRHistoPlotter::makeLabels() {
	for (size_t i = 0; i != fileNamesMC.size(); i++) {
		TString s1 = fileNamesMC.at(i);
		TString s2 = "_.";
		TObjArray* tokens = s1.Tokenize(s2);
		std::string aLabel = ((TObjString*) (tokens->At(1)))->String().Data();
		std::string aLabel2 = ((TObjString*) (tokens->At(2)))->String().Data();
		labels.push_back(aLabel);
		//    labels.push_back(aLabel + aLabel2);
	}
	std::cout << "Labels MC done" << std::endl;

	for (size_t i = 0; i != fileNamesMCSig.size(); i++) {
		TString s1 = fileNamesMCSig.at(i);
		TString s2 = "_.";
		TObjArray* tokens = s1.Tokenize(s2);

		std::string aLabelType =
				((TObjString*) (tokens->At(1)))->String().Data();

		std::string aLabelCoupling =
				((TObjString*) (tokens->At(2)))->String().Data();

		//  std::string aLabelMass = ((TObjString*)(tokens->At(3)))->String().Data();
		//  std::cout << "Right3?" << std::endl;

		std::cout << s1.Data();

		std::string aLabel = aLabelType + aLabelCoupling;

		labelsSig.push_back(aLabel);

	}
	std::cout << "Labels MC signal done" << std::endl;

}

///set output directories for plots.
void EDBRHistoPlotter::setOutDir(std::string outDirNew,int isBarrel) {
	char buffer[256];
	nameOutDir_ = outDirNew;

        if(isBarrel==1)
		sprintf(buffer, "%s/pdf_b", nameOutDir_.c_str());
        else if(isBarrel==0)    sprintf(buffer, "%s/pdf_e", nameOutDir_.c_str());
	else  sprintf(buffer, "%s/pdf_a", nameOutDir_.c_str());
	printf("%s\n", buffer);
	gSystem->mkdir(buffer, true);

        if(isBarrel==1)
		sprintf(buffer, "%s/root_b", nameOutDir_.c_str());
	else if(isBarrel==0)    sprintf(buffer, "%s/root_e", nameOutDir_.c_str());
	else sprintf(buffer, "%s/root_a", nameOutDir_.c_str());

	printf("%s\n", buffer);
	gSystem->mkdir(buffer, true);

}

void EDBRHistoPlotter::makeStackPlots(std::string histoName,TString histTitle,int isBarrel,vector<std::string> bkgLegend) {

	cleanupMC();
	cleanupMCSig();
	cleanupDATA();

	//printf("Making histo %s\n",histoName.c_str());
	std::cout << "\rMaking histo " << histoName.c_str() << std::endl;

	TCanvas* cv = new TCanvas(("cv_" + histoName).c_str(),
			("cv_" + histoName).c_str(), 600, 600);

	//create 3 pads in the canvas
	TPad* fPads1 = NULL;
	TPad* fPads2 = NULL;
    //  TPad* fPads3 = NULL;

	if (makeRatio_ && isDataPresent_) {
		fPads1 = new TPad("pad1", "", 0.00, 0.30, 0.99, 0.99);
		fPads2 = new TPad("pad2", "", 0.00, 0.00, 0.99, 0.30);
		//fPads2 = new TPad("pad2", "", 0.00, 0.20, 0.99, 0.40);
		//fPads3 = new TPad("pad3", "", 0.00, 0.00, 0.99, 0.20);
		fPads1->SetFillColor(0);
		fPads1->SetLineColor(0);
		fPads2->SetFillColor(0);
		fPads2->SetLineColor(0);
		fPads1->SetBottomMargin(0.012);
		fPads2->SetTopMargin(0.005);
		fPads2->SetBottomMargin(0.5);

		//fPads3->SetFillColor(0);
		//fPads3->SetLineColor(0);
		fPads1->Draw();
		fPads2->Draw();
		//fPads3->Draw();
	}

	//============ Data vs MC plots ==============

	if (makeRatio_ && isDataPresent_) {
		fPads1->cd();
	}
	// fPads1->cd();
	///--------------------
	/// Make the DATA stack
	///--------------------

	TH1D* sumDATA = NULL;

	for (size_t i = 0; i != fileNamesDATA.size(); ++i) {
		filesDATA.push_back(
				TFile::Open((nameInDir_ + fileNamesDATA.at(i)).c_str()));
	}

	for (size_t i = 0; i != filesDATA.size(); ++i) {
		TH1D* histo = (TH1D*) (filesDATA.at(i)->Get(histoName.c_str())->Clone("clone"));
		histo->SetDirectory(0);
		histosDATA.push_back(histo);
	}
	if (histosDATA.size() != 0) {
		sumDATA = (TH1D*) (histosDATA.at(0)->Clone("masterDATA"));
		sumDATA->Reset();
		sumDATA->SetDirectory(0);
	}
	//sumDATA->SetBinErrorOption(TH1::kPoisson);

	for (size_t i = 0; i != histosDATA.size(); ++i) {
		sumDATA->Add(histosDATA.at(i));
	}
	sumDATA->SetBinErrorOption(TH1::kPoisson);
	double sumDataIntegral = 0;
	if (isDataPresent_)
		sumDataIntegral = sumDATA->Integral();

	///------------------
	/// Make the MC stack
	///------------------

	TH1D* sumMC = NULL;
	TH1D* hdiff = NULL;
	double sumBkgAtTargetLumi = 0;

	for (size_t i = 0; i != fileNamesMC.size(); ++i) {
		filesMC.push_back(
				TFile::Open((nameInDir_ + fileNamesMC.at(i)).c_str()));
	}

	for (size_t i = 0; i != filesMC.size(); ++i) {
		TH1D* histo = (TH1D*) (filesMC.at(i)->Get(histoName.c_str())->Clone(
				labels.at(i).c_str()));
		histo->SetDirectory(0);
		histo->SetFillColor(getFillColor(i));
		histo->SetLineColor(getFillColor(i));

		TString filename = filesMC.at(i)->GetName();
		histo->Scale(kFactorsMC_.at(i));
		histosMC.push_back(histo);
	}

	if (histosMC.size() != 0) {
		sumMC = (TH1D*) (histosMC.at(0)->Clone("masterMC"));
		sumMC->Reset();
		sumMC->SetDirectory(0);
		hdiff = (TH1D*) (histosMC.at(0)->Clone("masterMC"));
		hdiff->Reset();
		hdiff->SetDirectory(0);
	}

	/// Do we normalize to data or to lumi?
	/// NOTICE THAT THIS DEPENDS ON THE HISTOGRAMS HAVING BEING
	/// CORRECTLY FILLED WITH PUweight*LumiWeight*GenWeight
	for (size_t is = 0; is != histosMC.size(); is++) {
		if (scaleToData_ && isDataPresent_) {
			histosMC.at(is)->Scale(
					targetLumi_ * sumDataIntegral / sumBkgAtTargetLumi);
		} 
	      /*else 
		{
			histosMC.at(is)->Scale(targetLumi_);
		}*/
	}

	THStack* hs = new THStack("hs", "");

     // Make a histogram just for the sum
	for (size_t i = 0; i != histosMC.size(); ++i) {
		histosMC.at(i)->SetFillColor(getFillColor(i));
		sumMC->Add(histosMC.at(i));

		hs->Add(histosMC.at(i));
	}

	sumMC->SetFillStyle(3005);
	sumMC->SetFillColor(1);
        sumMC->SetMarkerSize(0);
      //sumMC->SetLineColor(kBlack);
        sumMC->SetLineColor(0);
	sumMC->SetLineWidth(2);

	if (scaleToData_ && isDataPresent_) {
		std::cout << "===> Residual DATA/MC Scale Factor is: "
				<< sumDataIntegral / sumBkgAtTargetLumi << std::endl;
	}

	///-------------------------------
	/// Add the MC signal to the stack
	///-------------------------------

	for (size_t i = 0; i != fileNamesMCSig.size(); ++i) {
		filesMCSig.push_back(
				TFile::Open((nameInDir_ + fileNamesMCSig.at(i)).c_str()));
	}

	for (size_t i = 0; i != filesMCSig.size(); ++i) {

		TH1D* histo = (TH1D*) (filesMCSig.at(i)->Get(histoName.c_str())->Clone(
				labelsSig.at(i).c_str()));
		TH1D* histoOrig =
				(TH1D*) (filesMCSig.at(i)->Get(histoName.c_str())->Clone(
						labelsSig.at(i).c_str()));
		histo->SetDirectory(0);
		histo->SetLineColor(kRed-7);
		histo->SetFillColor(kRed-7);

		histoOrig->SetDirectory(0);
		histoOrig->SetLineColor(getLineColor(i));
		histoOrig->SetFillColor(getLineColor(i));
		if (i % 2 == 0)
			histoOrig->SetFillStyle(3004);
		else
			histoOrig->SetFillStyle(3005);
		//histo->Scale(kFactor_); //============= SCALE FACTORS FOR SIGNAL? ==== FIXME

		histosMCSig.push_back(histo);
		histosMCSigOrig.push_back(histoOrig);
	}

	//scale the MC signal histogram
	if (histosMCSig.size() != kFactorsSig_.size())
		cout << "+++++++++++++++++ Mismatch in size of input MC Sig arrays !!!"
				<< endl;

	for (size_t is = 0; is != histosMCSig.size(); is++) {

		histosMCSig.at(is)->Scale(/*targetLumi_ **/ kFactorsSig_.at(is));
		histosMCSigOrig.at(is)->Scale(/*targetLumi_ **/ kFactorsSig_.at(is));

		hs->Add(histosMCSig.at(is));
		sumMC->Add(histosMCSig.at(is));
		//add the signal to the total background
		// histosMCSig.at(is)->Add(sumMC);
	}

	///-----------------------------------
	/// Draw both MC and DATA in the stack
	///-----------------------------------
	//hs->SetBinErrorOption(TH1::kPoisson);
	double max=hs->GetMaximum();
	hs->SetMaximum(max*50);
	hs->Draw("HIST ");
	hs->GetYaxis()->SetTitle("Events/bin");  //40.24pb-1");
//	hs->GetYaxis()->SetNdivisions(404);
	hs->GetYaxis()->SetLabelSize(0.06);  //40.24pb-1");
	hs->GetYaxis()->SetTitleOffset(1.17);
//  hs->GetYaxis()->CenterTitle();

	///------------------------------------
	///
	///for error band
	///
	///------------------------------------
	TH1D* htemp = (TH1D*)(hs->GetStack()->Last());
	htemp->SetBinErrorOption(TH1::kPoisson);
	const int binsize = htemp->GetSize()-2;
	Double_t x[binsize];
	Double_t y[binsize];
	Double_t xerror_l[binsize];
	Double_t xerror_r[binsize];
        Double_t yerror_u[binsize];
        Double_t yerror_d[binsize];
	for(int i=0;i<binsize;i++)
	{
		x[i]=htemp->GetBinCenter(i+1);
		y[i]=htemp->GetBinContent(i+1);
		xerror_l[i]=0.5*htemp->GetBinWidth(i+1);
		xerror_r[i]=0.5*htemp->GetBinWidth(i+1);
		yerror_u[i]=htemp->GetBinErrorUp(i+1);
		yerror_d[i]=htemp->GetBinErrorLow(i+1);
		if(htemp->GetBinContent(i+1)==0)
		{yerror_u[i]=0.;yerror_d[i]=0.;}
	}
	TGraphAsymmErrors* gr = new TGraphAsymmErrors(binsize, x, y, xerror_l,xerror_r, yerror_d, yerror_u); 
	gr->SetFillColor(1);
	gr->SetFillStyle(3005);
	gr->Draw("SAME 2");

	double maximumMC = 3.5 * sumMC->GetMaximum();
	double maximumDATA = -100;
	if (isDataPresent_)
		maximumDATA = 1.15 * sumDATA->GetMaximum();
	double maximumForStack = -100;
	if (isDataPresent_)
		maximumForStack = (maximumMC > maximumDATA ? maximumMC : maximumDATA);
	else
		maximumForStack = maximumMC;
	hs->SetMaximum(maximumForStack);
	// Some hacks for better aestetics
	// Extra vertical space in eta plots
	hs->SetMinimum(0.1);
        TLatex latex;
        latex.SetTextSize(0.1);
	latex.SetTextFont(42);
        double poi_X=htemp->GetBinLowEdge(htemp->GetNbinsX()/3);
        latex.DrawLatex(poi_X,maximumMC*1.5,"#font[132]{e}#mu#gamma_{#font[132]{barrel}}");
	if (isDataPresent_) {
		sumDATA->SetMarkerColor(1);
		sumDATA->SetMarkerStyle(20);
                TGraphAsymmErrors * g = new TGraphAsymmErrors(sumDATA);
		for (int i = 0; i < g->GetN(); ++i) {
			double x=sumDATA->GetXaxis()->GetBinCenter(i+1);
			double y=sumDATA->GetBinContent(i+1);
			double ey=sumDATA->GetBinError(i+1);
			g->SetPoint(i,x,y);
			g->SetPointEXlow(i, 0);
			g->SetPointEXhigh(i, 0);
			g->SetPointEYlow(i, ey);
			g->SetPointEYhigh(i, ey);
		}
		g->SetFillColor(kBlack);
		g->SetLineColor(1);
		g->SetMarkerStyle(20);
		g->SetMarkerSize(1.0);
		g->Draw("SAME E P");
	}

	// histosMCSig.at(0)->Draw("SAME HIST");

	// For the legend, we have to tokenize the name "histos_XXX.root"
        TLegend* leg1 = new TLegend(0.52, 0.5, 0.82, 0.9);
        TLegend* leg2 = new TLegend(0.2, 0.5, 0.5, 0.9);
        TLegend* leg3 = new TLegend(0.2, 0.63,0.57, 0.7);
        leg1->SetTextSize(0.04);
        leg2->SetTextSize(0.04);
        ofstream ftxt("./yields.txt");
	if (isDataPresent_){
		double yieldsDataErr=0;
	      //double yieldsData = sumDATA->IntegralAndError(0,sumDATA->GetNbinsX(),yieldsDataErr);;
		double yieldsData = sumDATA->IntegralAndError(0,sumDATA->GetNbinsX(),yieldsDataErr);
		char yData[100];sprintf(yData,"%0.f",yieldsData);
		char yDataErr[100];sprintf(yDataErr,"%0.f",yieldsDataErr);
		TString samplesData = "Data";
  	        TString LabelData = samplesData +" ["+ yData+"#pm"+yDataErr+"]";
//		TString LabelData = samplesData; 
		ftxt<<samplesData.Data()<<" "<<yData<<"$pm$"<<yDataErr<<""<<endl;

		leg1->AddEntry(sumDATA, LabelData, "ep");
	      //double yieldsMC=sumMC->GetSumOfWeights();
		double yieldsMCerr;
		double yieldsMC=sumMC->IntegralAndError(0,sumMC->GetNbinsX(),yieldsMCerr);
		char yMC[100];sprintf(yMC,"%0.f",yieldsMC);
		char yMCerr[100];sprintf(yMCerr,"%0.f",yieldsMCerr);
		TString samplesMC = "Pred. #oplus stat.";
	      TString LabelMC = samplesMC +" ["+ yMC+ "#pm"+yMCerr+" ]";
//		TString LabelMC = samplesMC ;
		ftxt<<samplesMC<<" "<<yMC<< "$pm$"<<yMCerr<<""<<endl;
		leg1->AddEntry(sumMC, LabelMC);

		if (histosMCSig.size() > 0) {
			char rescalingLabel[64];
			for (size_t i = 0; i != histosMCSig.size(); ++i) {
				sprintf(rescalingLabel, " (x%g)", kFactorsSig_.at(i));
				std::string rescalingStr(rescalingLabel);
				TH1D* h1;
				if (kFactorsSig_.at(i) != 2.0) {
					if (i == 0){
						h1=(TH1D*)histosMCSig.at(i)->Clone();
						//h1->Add(histosMCSig.at(i+1),1);
						//h1->Add(histosMCSig.at(i+2),1);
						double yieldsMCSigErr;
						double yieldsMCSig= h1->IntegralAndError(0,h1->GetNbinsX(),yieldsMCSigErr);
						yieldsMCSig = h1->GetSumOfWeights();
						char ySig[100];sprintf(ySig,"%.1f",yieldsMCSig);
						char ySigErr[100];sprintf(ySigErr,"%.1f",yieldsMCSigErr);
						TString samplesMCSig = "WW#gamma";
					        TString LabelSig = samplesMCSig +" ["+ ySig + "#pm"+ySigErr +"]";
//						TString LabelSig = samplesMCSig;
						ftxt<<samplesMCSig<<" "<<ySig<< "$pm$"<<ySigErr<<""<<endl;
						leg1->AddEntry(histosMCSig.at(i), LabelSig, "lf");
					}
				} else
					leg1->AddEntry(histosMCSig.at(i), (labelsSig.at(i)).c_str(),
							"lf");
			}
		}
	}
	cout<<"test"<<endl;
	int num=0;
	for(size_t i = 0; i != histosMC.size(); ++i){
		mcTotalLabels.push_back(bkgLegend.at(i));
                double y_tmp=histosMC.at(i)->GetSumOfWeights();
                if(y_tmp>0) num++;
	}
	num=(num+3)/2-3;
	cout<<"num "<<num<<endl;
	for (size_t i = 0; i != histosMC.size(); ++i) {
		TH1D* h1=(TH1D*)histosMC.at(i)->Clone();
		double yerr;
		double yields=h1->IntegralAndError(0,histosMC.at(i)->GetNbinsX(),yerr);
		yields = h1->GetSumOfWeights();
                if(yields<0){yields=0;yerr=0;}
		char y[100];sprintf(y,"%.1f",yields);
		char ye[100];sprintf(ye,"%.1f",yerr);
		TString samples ;
		if(mcTotalLabels.at(i).find("plj")!=string::npos)
			samples="Nonprompt #gamma";
		else if(mcTotalLabels.at(i).find("fakeL")!=string::npos)
			samples="Nonprompt #font[12]{l}";
		else    samples= mcTotalLabels.at(i).c_str();
		TString LabelMC = samples +" ["+ y+ "#pm"+ye+"]";
//		TString LabelMC = samples ;
                if(yields>0){
			if(i<num)leg1->AddEntry(h1, LabelMC, "f");
			else   leg2->AddEntry(h1, LabelMC, "f");
		}
		ftxt<<samples<<" "<<y<< "$pm$"<<ye<<""<<endl;
		cout<<LabelMC<<endl;
	}

	leg1->SetFillColor(kWhite);
	leg1->Draw();
	leg2->Draw();
      //leg2->SetFillColor(kWhite);
      //leg2->Draw("same");
      //leg3->Draw("same");
      // Nice labels
      //TMathText* l = makeCMSPreliminaryTop(13, 0.50, 0.935);
      //TMathText* l = makeCMSLumi(13,551.7,0.6,0.935);
      //l->Draw();
	//============ Data-MC/Error ==============

	TLine* lineAtZero = NULL;
	TLine* lineAtPlusTwo = NULL;
	TLine* lineAtMinusTwo = NULL;
	if (makeRatio_ && isDataPresent_) {
		fPads2->cd();

		fPads2->SetGridx();
		fPads2->SetGridy();

		double thisYmin = 0.5;
		double thisYmax = 1.8;

		TVectorD nsigma_x(sumDATA->GetNbinsX());
		TVectorD nsigma_y(sumDATA->GetNbinsX());
		TVectorD nsigma_ex(sumDATA->GetNbinsX());
		TVectorD nsigma_ey(sumDATA->GetNbinsX());
		TH1D*nominal=(TH1D*)sumMC->Clone("nominal");
		TH1D*nomNoErr=(TH1D*)nominal->Clone("nomNoErr");
		nominal->SetMarkerSize(0);
		nominal->SetFillColor(1);
		nominal->SetLineColor(1);
		nominal->SetFillStyle(3005);
		for (int i = 1; i<= nomNoErr->GetNbinsX(); ++i){nomNoErr->SetBinError(i,0);}
		nominal->Divide(nomNoErr);
		nominal->SetTitle("");
		nominal->Draw("E2");
		for (int ibin = 0; ibin != sumDATA->GetNbinsX(); ++ibin) {

			double Data = sumDATA->GetBinContent(ibin + 1);
			double Bkg = sumMC->GetBinContent(ibin + 1);
			double eData = sumDATA->GetBinError(ibin + 1);
			double eBkg = sumMC->GetBinError(ibin + 1);
			double x = sumDATA->GetBinCenter(ibin + 1);

			if (Data !=0 && Bkg !=0 ){
				nsigma_x[ibin] = x;
				nsigma_y[ibin] = Data / Bkg;
				nsigma_ex[ibin] = 0;
				nsigma_ey[ibin] = sqrt((eData*eData/(Bkg*Bkg)) + (Data*Data*eBkg*eBkg/(Bkg*Bkg*Bkg*Bkg)));
			}
			else {
				nsigma_x[ibin] = +999999;
				nsigma_y[ibin] = 0;
				nsigma_ex[ibin] = 0;
				nsigma_ey[ibin] = 0;
			}


		}

		if (nsigma_x.GetNoElements() != 0) {
			TGraph *nsigmaGraph = new TGraphErrors(nsigma_x, nsigma_y,nsigma_ex,nsigma_ey);
			nominal->SetTitle("");
			nominal->GetYaxis()->SetRangeUser(thisYmin, thisYmax);
			nominal->GetYaxis()->SetNdivisions(404);
			nominal->GetYaxis()->SetTitle("Data/Pred.");
//			nominal->GetYaxis()->CenterTitle();
			nominal->GetYaxis()->SetTitleOffset(0.42);
			nominal->GetXaxis()->SetLabelOffset(0.04);
			nominal->GetYaxis()->SetTitleSize(0.14);
			nominal->GetYaxis()->SetLabelSize(0.13);
			nominal->GetXaxis()->SetLabelSize(0.13);
			nominal->GetXaxis()->SetTitleSize(0.20);
			TString title=histoName.c_str();
		      //nominal->GetXaxis()->SetTitle(histTitle+" [GeV]");
			nominal->GetXaxis()->SetTitle(histTitle);
			nominal->GetXaxis()->SetLimits(sumMC->GetXaxis()->GetXmin(),
					sumMC->GetXaxis()->GetXmax());
			nominal->GetXaxis()->SetRangeUser(sumMC->GetXaxis()->GetXmin(),
					sumMC->GetXaxis()->GetXmax());
			nominal->GetXaxis()->SetTitleOffset(0.9);
			nsigmaGraph->SetMarkerSize(1.0);
			nsigmaGraph->Draw("same pe");
		}

		fPads2->Update();

		lineAtZero = new TLine(sumMC->GetXaxis()->GetXmin(), 1,
				sumMC->GetXaxis()->GetXmax(), 1);
		lineAtZero->SetLineColor(2);
		lineAtZero->Draw();
		lineAtPlusTwo = new TLine(sumMC->GetXaxis()->GetXmin(), 1,
				sumMC->GetXaxis()->GetXmax(), 1);
		lineAtPlusTwo->SetLineColor(2);
		lineAtPlusTwo->SetLineStyle(2);
		lineAtPlusTwo->Draw();
		lineAtMinusTwo = new TLine(sumMC->GetXaxis()->GetXmin(), 1,
				sumMC->GetXaxis()->GetXmax(), 1);
		lineAtMinusTwo->SetLineColor(2);
		lineAtMinusTwo->SetLineStyle(2);
		lineAtMinusTwo->Draw();
	}

      //char bufferstr[20];
      //sprintf(bufferstr, "%.3f",targetLumi_);
      //CMS_lumi(fPads1,4,0, bufferstr);
	std::ostringstream strs;
	strs << targetLumi_;
	std::string lumivalue = strs.str();
	CMS_lumi(fPads1, 4, 0, lumivalue);

	// Save the picture
	char buffer[256];
	cv->SetLogy(false);
        if(isBarrel==1)	
		                 sprintf(buffer, "%s/root_b/can_%s.root", nameOutDir_.c_str(),histoName.c_str());
	else if(isBarrel==0)     sprintf(buffer, "%s/root_e/can_%s.root", nameOutDir_.c_str(),histoName.c_str());
	else                     sprintf(buffer, "%s/root_a/can_%s.root", nameOutDir_.c_str(),histoName.c_str());
	cv->SaveAs(buffer);

	if(isBarrel==1)
		               sprintf(buffer, "%s/pdf_b/can_%s.pdf", nameOutDir_.c_str(),histoName.c_str());
	else if(isBarrel==0)   sprintf(buffer, "%s/pdf_e/can_%s.pdf", nameOutDir_.c_str(),histoName.c_str());
	else                   sprintf(buffer, "%s/pdf_a/can_%s.pdf", nameOutDir_.c_str(),histoName.c_str());
	cv->SaveAs(buffer);

	if (makeRatio_ && isDataPresent_) {
		fPads1->cd();
		fPads1->SetLogy(true);
	} else {
		cv->SetLogy(true);
		fPads1->SetLogy(true);
	}
	//-- resize y axis --
	hs->SetMaximum(70 * maximumForStack);
	//
	cv->SetLogy(true);
	cv->Update();
	cv->Modified(true);
	cv->Draw();
        if(isBarrel==1)
		sprintf(buffer, "%s/pdf_b/LOG_can_%s.pdf", nameOutDir_.c_str(),histoName.c_str());
        else if(isBarrel==0)    sprintf(buffer, "%s/pdf_e/LOG_can_%s.pdf", nameOutDir_.c_str(),histoName.c_str());
	else  sprintf(buffer, "%s/pdf_a/LOG_can_%s.pdf", nameOutDir_.c_str(),histoName.c_str());
	cv->SaveAs(buffer);

}

