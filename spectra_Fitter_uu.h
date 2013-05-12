//#include <iostream>
#include <TObjArray.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TString.h>
#include <TH1F.h>
#include <TSystem.h>
#include <TPad.h>
#include <TMath.h>
#include <TF1.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TSpectrum.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TStyle.h>

using namespace std;

//Create TGraphErrors to hold the spectra
TGraphErrors *spectraGraph[2][1][6];

//Create TGraphErrors to hold the fit results
TGraphErrors *fitResultGraph[6][4];


//Define Particle Masses
Double_t mPion = .139570;
Double_t mKaon = .493667;
Double_t mProton = .938272;


//General Naming Properties
TString collision = "UU";
TString energy = "_193GeV_";
TString extension = ".txt";
TString gif = ".gif";
TString png = ".png";
TString CentBin = "_CentBin";
//TString nCentBin[9] = {"0","1","2","3","4",
//	"5","6","7","8"}; 
TString particle[6] =
{"piPlus",
	"piMinus",
	"kPlus",
	"kMinus",
	"pPlus",
	"pMinus"};


TString *histoNameArray[6];
TFile *histoFile;

//TPC_TOF = 0 for TPC
//TPC_TOF = 1 for TOF

void setHistoNameArray(Int_t CENTBIN, Int_t TPC_TOF){
	// **removed CENTBIN dependence

	if (TPC_TOF == 1){
		histoNameArray[0]=piPlus_TOF;
		histoNameArray[1]=piMinus_TOF;
		histoNameArray[2]=kPlus_TOF;
		histoNameArray[3]=kMinus_TOF;
		histoNameArray[4]=pPlus_TOF;
		histoNameArray[5]=pMinus_TOF;
	}

	//Now for TPC Histograms
	if (TPC_TOF == 0){
		histoNameArray[0]=piPlus_TPC;
		histoNameArray[1]=piMinus_TPC;
		histoNameArray[2]=kPlus_TPC;
		histoNameArray[3]=kMinus_TPC;
		histoNameArray[4]=pPlus_TPC;
		histoNameArray[5]=pMinus_TPC;
	}

}//End setHistoNameArray


TFile *treeFile;
TTree *eventTree;
//Double_t nCentEvents[9];
Double_t nGoodEvents;

void countNGoodEvents(){
/*
	Float_t vpdVz,Vz,refMult,vertexMostTracks_f;
	nGoodEvents=0;

	//Load the Event Tree
	eventTree = (TTree*)treeFile->Get("mEventInfo");
	eventTree->SetMakeClass(1);

	//Get the total number of events in event tree
	Long64_t nTotalEvents = eventTree->GetEntries();
	cout <<"Total Number of Events: " <<nTotalEvents <<endl;

	//Load data
	eventTree->SetBranchAddress("vpdVz",&vpdVz);
	eventTree->SetBranchAddress("Vz",&Vz);
	eventTree->SetBranchAddress("refMult",&refMult);
	eventTree->SetBranchAddress("vertexMostTracks",&vertexMostTracks_f);

	//loop over the event tree entries
	for(Long64_t i=0; i<nTotalEvents; i++){
		eventTree->GetEntry(i);

//		if (refMult <= 200) { continue; }
		if (vpdVz >=30 || vpdVz <=-30) { continue; }
		if (vpdVz-Vz >=5 || vpdVz-Vz <=-5) { continue; }
//		if ( (int)vertexMostTracks != 0 ) { continue; }

		nGoodEvents++;
	}

	gSystem->Sleep(100);
*/

	nGoodEvents=6723459; //for completeish data set - all days, with some missing jobs
}


//Define the Student T Distribution for Fitting
Double_t studentT(Double_t *x, Double_t *par){

	//Parameter Definitions
	//par[0] is the Amplitude
	//par[1] is the "width" (NDF of Student T Dist.)
	//par[2] is the shift

	Float_t xx = x[0];
	Double_t f = TMath::Gamma((par[1]+1.)/2.)/(sqrt(TMath::Pi()*par[1])*TMath::Gamma(par[1]/2.));
	f = f*pow((1+pow(xx-par[2],2)/par[1]),-(par[1]+1.)/2.);
	f = f*par[0];

	return f;
}


void drawStudentT(){
	TF1 *f1 = new TF1("studentT",studentT,-1,5,3);
	f1->SetParameters(100,.01,2.);
	f1->SetParNames("Amplitude","Width","Shift");
	f1->Draw();

}

//Make Legends
TLegend *legend, *legend1;
TLegendEntry *entry;
/* void makeLegends(Int_t TPC_TOF){ */

//Set the draw options of all the Spectra and Draws them
void drawSpectra(Int_t TPC_TOF){
	cout <<"Drawing the Spectra!" <<endl;

	//Define Variables
	Int_t markerStyle;

	//Create a Canvas
	TCanvas *spectraCanvas[2];
	if (TPC_TOF == 1){
		spectraCanvas[TPC_TOF] = new TCanvas("tofCanvas","tofCanvas",0,0,1600,900);
		spectraCanvas[TPC_TOF]->Divide(3,2);
		markerStyle = 24;
	}
	else if (TPC_TOF == 0){
		spectraCanvas[TPC_TOF] = new TCanvas("tpcCanvas","tpcCanvas",0,0,1600,900);
		spectraCanvas[TPC_TOF]->Divide(3,2);
		markerStyle = 20;
	}

	//Define the Centrality Colors
	/* CentBin     |  Color
0: 70-80%   | kCyan+2
1: 60-70%   | kSpring+5
2: 50-60%   | kMagenta+4
3: 40-50%   | kOrange+4
4: 30-40%   | kOrange-2
5: 20-30%   | kMagenta+2
6: 10-20%   | kBlue+1
7: 05-10%   | kGreen+2
8: 00-05%   | kRed
	 */

	for (Int_t i=0; i<6; i++){
		spectraGraph[TPC_TOF][0][i]->SetMarkerStyle(markerStyle);
		spectraGraph[TPC_TOF][0][i]->SetMarkerColor(kCyan+2);
		spectraGraph[TPC_TOF][0][i]->SetLineColor(kCyan+2);
	}

	//Set Properties of the First Spectra for Drawing purposes 
	TString title[2] = {"dE/dx Spectra","Time of Flight Spectra"};
	spectraGraph[TPC_TOF][0][0]->SetTitle("#pi^{+} "+title[TPC_TOF]);
	spectraGraph[TPC_TOF][0][1]->SetTitle("#pi^{-} "+title[TPC_TOF]);
	spectraGraph[TPC_TOF][0][2]->SetTitle("k^{+} "+title[TPC_TOF]);
	spectraGraph[TPC_TOF][0][3]->SetTitle("k^{-} "+title[TPC_TOF]);
	spectraGraph[TPC_TOF][0][4]->SetTitle("p "+title[TPC_TOF]);
	spectraGraph[TPC_TOF][0][5]->SetTitle("#bar{p} "+title[TPC_TOF]);

	spectraGraph[TPC_TOF][0][0]->GetYaxis()->SetRangeUser(1e-1,1e3);
	spectraGraph[TPC_TOF][0][1]->GetYaxis()->SetRangeUser(1e-1,1e3);
	spectraGraph[TPC_TOF][0][2]->GetYaxis()->SetRangeUser(1e-3,1e1);
	spectraGraph[TPC_TOF][0][3]->GetYaxis()->SetRangeUser(1e-3,1e1);
	spectraGraph[TPC_TOF][0][4]->GetYaxis()->SetRangeUser(1e-3,1e1);
	spectraGraph[TPC_TOF][0][5]->GetYaxis()->SetRangeUser(1e-3,1e0);


	for (Int_t i=0; i<6; i++){
		spectraGraph[TPC_TOF][0][i]->GetYaxis()->SetTitle("#frac{1}{2#pi m_{T}} #frac{d^{2}N}{dm_{T}dy}");
		spectraGraph[TPC_TOF][0][i]->GetXaxis()->SetTitle("p_{T} (GeV)");
//		spectraGraph[TPC_TOF][0][i]->GetYaxis()->SetTitleOffset(1.75);
	}

	//Create a Legend
	legend = new TLegend(0.5,0.67,0.88,0.88);
	legend->SetFillColor(kWhite);

	//Create the Legend Entries
	entry = new TLegendEntry();
	entry = legend->AddEntry("","0 - 1% Central","lp");
	entry->SetLineColor(spectraGraph[TPC_TOF][0][1]->GetLineColor());
	entry->SetMarkerStyle(spectraGraph[TPC_TOF][0][1]->GetMarkerStyle());
	entry->SetMarkerColor(spectraGraph[TPC_TOF][0][1]->GetMarkerColor());


	//gSystem->Sleep(10000);

	//Draw the First Spectra for each particle so the frame will be set
	spectraCanvas[TPC_TOF]->cd(1);
	gPad->SetLogy();
	spectraGraph[TPC_TOF][0][0]->Draw("ALP");

	spectraCanvas[TPC_TOF]->cd(2);
	gPad->SetLogy();
	spectraGraph[TPC_TOF][0][1]->Draw("ALP");

	spectraCanvas[TPC_TOF]->cd(3);
	gPad->SetLogy();
	spectraGraph[TPC_TOF][0][2]->Draw("ALP");

	spectraCanvas[TPC_TOF]->cd(4);
	gPad->SetLogy();
	spectraGraph[TPC_TOF][0][3]->Draw("ALP");

	spectraCanvas[TPC_TOF]->cd(5);
	gPad->SetLogy();
	spectraGraph[TPC_TOF][0][4]->Draw("ALP");

	spectraCanvas[TPC_TOF]->cd(6);
	gPad->SetLogy();
	spectraGraph[TPC_TOF][0][5]->Draw("ALP");

/*	for(Int_t y=0; y<8; y++){
		spectraCanvas[TPC_TOF]->cd(1);
		spectraGraph[TPC_TOF][y][0]->Draw("LP");

		spectraCanvas[TPC_TOF]->cd(2);
		spectraGraph[TPC_TOF][y][1]->Draw("LP");

		spectraCanvas[TPC_TOF]->cd(3);
		spectraGraph[TPC_TOF][y][2]->Draw("LP");

		spectraCanvas[TPC_TOF]->cd(4);
		spectraGraph[TPC_TOF][y][3]->Draw("LP");

		spectraCanvas[TPC_TOF]->cd(5);
		spectraGraph[TPC_TOF][y][4]->Draw("LP");

		spectraCanvas[TPC_TOF]->cd(6);
		spectraGraph[TPC_TOF][y][5]->Draw("LP");
	}
*/
	spectraCanvas[TPC_TOF]->cd(1);
	legend->Draw();

}



//Draw the TOF and TPC Spectra Together
void drawTPCandTOF(){

	//Create a Canvas
	TCanvas *combinedCanvas = new TCanvas("combinedCanvas","combinedCanvas",0,0,1600,900);
	combinedCanvas->Divide(3,2);
	combinedCanvas->SetLeftMargin(.15);
	combinedCanvas->SetRightMargin(.15);

	//Create a Frame and Draw it on Each of the Pads
	TH1F *frame[6];
	frame[0] = combinedCanvas->cd(1)->DrawFrame(0,1e-2,1,1e3);
	frame[1] = combinedCanvas->cd(2)->DrawFrame(0,1e-2,1,1e3);
	frame[2] = combinedCanvas->cd(3)->DrawFrame(0,1e-3,1,5e1);
	frame[3] = combinedCanvas->cd(4)->DrawFrame(0,1e-3,1,5e1);
	frame[4] = combinedCanvas->cd(5)->DrawFrame(0,1e-3,1,1e2);
	frame[5] = combinedCanvas->cd(6)->DrawFrame(0,.3e-3,1,5e0);

	//Set the Titles of the Frames
	frame[0]->SetTitle("#pi^{+} TPC and TOF Spectra");
	frame[1]->SetTitle("#pi^{-} TPC and TOF Spectra");
	frame[2]->SetTitle("K^{+} TPC and TOF Spectra");
	frame[3]->SetTitle("K^{-} TPC and TOF Spectra");
	frame[4]->SetTitle("p TPC and TOF Spectra");
	frame[5]->SetTitle("#bar{p} TPC and TOF Spectra");

	//Set X Axis Titles
	for (Int_t i=0; i<6; i++){
		frame[i]->SetYTitle("#frac{1}{2#pi m_{T}} #frac{d^{2}N}{dm_{T}dy}");
		frame[i]->SetXTitle("p_{T} (GeV)");
		frame[i]->GetYaxis()->SetTitleOffset(1.75);
	}

	//Set the Pads to Log Scale
	for (Int_t i=1; i<7; i++){
		combinedCanvas->cd(i)->SetLogy();
		combinedCanvas->GetPad(i)->SetLeftMargin(.15);
		combinedCanvas->GetPad(i)->SetRightMargin(0.001);
	}

	//Draw All the Others
	for (Int_t i=0; i<2; i++){
		for (Int_t j=0; j<1; j++){ //changed to only do one cent bin
			for (Int_t k=0; k<6; k++){

				combinedCanvas->cd(k+1);
				spectraGraph[i][j][k]->Draw("LP");

			}
		}  
	}

	//Draw Legends
	combinedCanvas->cd(1);
	legend->Draw();
//	combinedCanvas->cd(2);
//	legend1->Draw();

}//End DrawTPCandTOF


//Do Efficiency Corrections
void efficiencyCorrection(){

	//Load the Files
//	TFile *file[6];
//	file[0] = new TFile("EfficiencyDocs/efficiency_allCents/piPlus_effFit19.root");
//	file[1] = new TFile("EfficiencyDocs/efficiency_allCents/piMinus_effFit19.root");
//	file[2] = new TFile("EfficiencyDocs/efficiency_allCents/kPlus_effFit19.root");
//	file[3] = new TFile("EfficiencyDocs/efficiency_allCents/kMinus_effFit19.root");
//	file[4] = new TFile("EfficiencyDocs/efficiency_allCents/pPlus_effFit19.root");
//	file[5] = new TFile("EfficiencyDocs/efficiency_allCents/pMinus_effFit19.root");

	//Load the Efficiency Functions
//	TF1 *effFunc[6];
//	effFunc[0] = (TF1 *)file[0]->Get("pipfunct_rp");
//	effFunc[1] = (TF1 *)file[1]->Get("pimfunct_rp");
//	effFunc[2] = (TF1 *)file[2]->Get("kapfunct_rp");
//	effFunc[3] = (TF1 *)file[3]->Get("kamfunct_rp");
//	effFunc[4] = (TF1 *)file[4]->Get("prpfunct_rp");
//	effFunc[5] = (TF1 *)file[5]->Get("prmfunct_rp");


	//Loop Over the Centrality Bins
/*	for (Int_t i=0; i<9; i++){

		//First we need to efficiency correct the TPC Spectra

		//Loop Over the Particle Species
		for (Int_t j=0; j<6; j++){

			//Temp Values to hold the pre-correct Point Values
			Double_t yieldTPC, mTm0TPC, yieldTPCErr, mTm0TPCErr, effCorrection;
			Int_t nPoints = spectraGraph[0][i][j]->GetN();

			//Get the Current Values of the Points and their errors
			//Then scale them according to the efficiency curve
			for (Int_t k=0; k<nPoints; k++){
				spectraGraph[0][i][j]->GetPoint(k,mTm0TPC,yieldTPC);
				yieldTPCErr = spectraGraph[0][i][j]->GetErrorY(k);
				mTm0TPCErr = spectraGraph[0][i][j]->GetErrorX(k); 

				//Get the Efficiency Correction from the Eff Curve
				effCorrection = effFunc[j]->Eval(mTm0TPC);

				//Set the Scaled point Values
				spectraGraph[0][i][j]->SetPoint(k,mTm0TPC,yieldTPC/effCorrection);
				spectraGraph[0][i][j]->SetPointError(k,mTm0TPCErr,yieldTPCErr/effCorrection);

			}

		}
	}
*/	

	//Set the TOF Scaled Value
	/*************************************************************************
	 **As a temporary solution we are scaling the TOF spectra by an additional
	 **value of .6
	 *************************************************************************/
/*

	for (Int_t i=0; i<9; i++){
		for (Int_t j=0; j<6; j++){
			Double_t yieldTPC, yieldTOF, mTm0TOF, mTm0TPC;
			Double_t yieldTPCErr, yieldTOFErr, mTm0TOFErr, mTm0TPCErr;
			Double_t effTPC, effTOF;

			//Get the Number of Points in the TOF Graph
			Int_t nPointsTOF = spectraGraph[1][i][j]->GetN();

			//Loop Over the points of the Graph
			for (Int_t k=0; k<nPointsTOF; k++){
				spectraGraph[1][i][j]->GetPoint(k,mTm0TOF,yieldTOF);
				mTm0TOFErr = spectraGraph[1][i][j]->GetErrorX(k);
				yieldTOFErr = spectraGraph[1][i][j]->GetErrorY(k);

				//Get the TPC Efficiency Correction
//				effTPC = 0.6*effFunc[j]->Eval(mTm0TOF);

				//Set the TOF Efficiency
				//effTOF = effTPC+.6;

				//Set the New Value of the TOF spectra
				spectraGraph[1][i][j]->SetPoint(k,mTm0TOF,yieldTOF/effTPC);
				spectraGraph[1][i][j]->SetPointError(k,mTm0TOFErr,yieldTOFErr/effTPC);  
			}

		}
	}
*/

	//Now the Offset Efficiency between the TPC curve and the TOF Curve 
	//needs to be computed and the TOF curve scaled up

	//Loop Over the Centrality Bins
	//  for (Int_t i=0; i<9; i++){

	/*   Double_t tpc_tof_Offset[5]; //One for each overlap bin - to be averaged later */


	//Loop Over the Particles
	//   for (Int_t j=0; j<6; j++){

	/*     Double_t offsetEff=0; */
	/*     Int_t nOverlapBins = 5; */
	/*     Double_t yieldTPC, yieldTOF, mTm0TOF, mTm0TPC; */
	/*     Double_t yieldTPCErr, yieldTOFErr, mTm0TOFErr, mTm0TPCErr; */
	/*     for (Int_t k=0; k<nOverlapBins; k++){ */
	/*   	spectraGraph[0][i][j]->GetPoint(k+11,mTm0TPC,yieldTPC); */
	/*   	spectraGraph[1][i][j]->GetPoint(k,mTm0TOF,yieldTOF); */

	/* 	//Initialize Offset Efficiency */
	/* 	tpc_tof_Offset[k] = 0.01; */

	/* 	//Find the efficiency by looping over a value until TOF matches TPC */
	/* 	do { */

	/* 	  //Scale TOF yield by Offset Efficiency */
	/* 	  //yieldTOF = yieldTOF*tpc_tof_Offset[k]; */

	/* 	  //Increment Offset Efficiency */
	/* 	  tpc_tof_Offset[k] += .05; */

	/* 	}while (yieldTOF*tpc_tof_Offset[k] < yieldTPC); */

	/* 	//Add the Found Offset Efficiency to the average */
	/* 	offsetEff += tpc_tof_Offset[k]; */

	/*     } */

	/*     //Find the Average Offset */
	/*     offsetEff = offsetEff/(nOverlapBins+0.0); */
	/*     cout <<offset <<endl; */

	/*     //Print the Offset Efficiencies */
	/*     for (Int_t k=0; k<nOverlapBins; k++){ */

	/*     } */

	/* offsetEff = +.6; */

	/* //Get the Number of Points in the TOF Graph */
	/* Int_t nPointsTOF = spectraGraph[1][i][j]->GetN(); */

	/* //Scale up the TOF spectra by the Average Offset */
	/* for (Int_t k=0; k<nPointsTOF; k++){ */

	/* 	//Get the Current Value of the TOF yeild */
	/* 	spectraGraph[1][i][j]->GetPoint(k,mTm0TOF,yieldTOF); */
	/* 	yieldTOFErr = spectraGraph[1][i][j]->GetErrorY(k); */
	/* 	mTm0TOFErr  = spectraGraph[1][i][j]->GetErrorX(k); */

	/* 	//Set the new point scaled up by the average offset */
	/* 	spectraGraph[1][i][j]->SetPoint(k,mTm0TOF,yieldTOF*offsetEff); */
	/* 	spectraGraph[1][i][j]->SetPointError(k,mTm0TOFErr,yieldTOFErr*offsetEff); */
	/* } */


}


//Write Spectra Values out to text File in three columns
void writeTextFiles(){

	//Directory for Saving Text Files
	TString headDir = "TextFiles/";

	//Variables used in loop
	TString OUTPUTFILE;
	Double_t mTm0, yield, yieldErr;
	Int_t nPointsTOF,nPointsTPC;


	for (Int_t j=0; j<9; j++){
		for(Int_t k=0; k<6; k++){

			//Set Output File Name
//			OUTPUTFILE = headDir+collision+energy+particle[k]+CentBin+nCentBin[j]+extension;
			OUTPUTFILE = headDir+collision+energy+particle[k]+extension;

			//Open the output file
			ofstream outfile (OUTPUTFILE);
			if (!outfile.is_open()){
				cout <<"ERROR: writeTextFiles() reports: OUTPUT FILE FAILED TO OPEN!";
				exit (EXIT_FAILURE);
			}

			cout <<"Writing Text File: " <<OUTPUTFILE <<endl;

			//Get the number of points in the graph
			nPointsTPC = spectraGraph[0][j][k]->GetN();
			nPointsTOF = spectraGraph[1][j][k]->GetN();

			//Loop Over the Points in the TPC Graph
			for (Int_t l=0; l<nPointsTPC; l++){

				//Get the Values from the Graphs
				spectraGraph[0][j][k]->GetPoint(l,mTm0,yield);
				yieldErr = spectraGraph[0][j][k]->GetErrorY(l);

				outfile <<setw(11) <<setprecision(6) <<mTm0
					<<setw(14) <<setprecision(6) <<yield
					<<setw(14) <<setprecision(6) <<yieldErr <<"\n";
			}

			for (Int_t l=0; l<nPointsTOF; l++){

				//Get the Values from the Graphs
				spectraGraph[1][j][k]->GetPoint(l,mTm0,yield);
				yieldErr = spectraGraph[1][j][k]->GetErrorY(l);

				outfile <<setw(11) <<setprecision(6) <<mTm0
					<<setw(14) <<setprecision(6) <<yield
					<<setw(14) <<setprecision(6) <<yieldErr <<"\n";
			}

			//Close the OutputFile
			outfile.close();

		}
	}

}//End Write Text Files

//Draw Single Plots
void drawSinglePlots(){

	//Create a new Canvas
	TCanvas *singleCanvas = new TCanvas();

	//combinedCanvas->cd(1)->Copy(singleCanvas);


}
