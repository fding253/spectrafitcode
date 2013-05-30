//Loads the histograms from Sam's root file and fits them with a 
//Student T distribution.  Then computes the particle yield and
//finally draws the spectra

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;


#include "histogramList_TOF.C"
#include "histogramList_TPC.C"
#include "spectra_Fitter_uu.h"
//#include "UtilityFunctions.h"
//#include "histoNameArray.C"

void spectra_CentBin(Int_t CENTBIN,Int_t TPC_TOF){

	Double_t ptBinSize;
	if (TPC_TOF==0) { ptBinSize = ptBinSize_TPC; }
	else if (TPC_TOF==1) { ptBinSize = ptBinSize_TOF; }

	cout <<"Doing Centrality Bin: " <<CENTBIN <<endl;

	//Define the Histogram Name Array
	setHistoNameArray(CENTBIN,TPC_TOF);

	//Create a Pointer to hold the histogram
	TH1F *htemp = new TH1F();

	//Create Objects needed for Peak Finding
	TSpectrum *tempSpectrum;
	if (TPC_TOF == 1)
		tempSpectrum = new TSpectrum(3,1.25);
	if (TPC_TOF == 0)
		tempSpectrum = new TSpectrum(4,1.25);

	//Define Variables needed for Peak Fitting
	Int_t nFoundPeaks, nPeak;
	Int_t *peakIndex, *peakBin; 
	Float_t *tempPeakLocation, *peakLocation, *peakAmplitude;
	Double_t protonMean[40];  //This is used to fix the means of the antiprotons
	Double_t protonWidth[40]; //This is used to fix the means of the antiprotons
	Double_t kaonWidth[6]; //This is used to store the first six widths of the kaon (and antikaon) to fix the remainer of the widths

	// **added: expected 1/beta mean **
	Double_t betaExpected[6][N_beta];
	Double_t dedxExpected[6][N_beta];

	TFile *dedxFile = new TFile("../dEdx_points.root","READ");
	TF1 *siPion = (TF1*) dedxFile->Get("siPion");
	TF1 *siKaon = (TF1*) dedxFile->Get("siKaon");
	TF1 *siProt = (TF1*) dedxFile->Get("siProton");
	TF1 *siElec = (TF1*) dedxFile->Get("siElectron");

	//Create Objects needed for Yield Extraction
	TFitResultPtr fitResult;
	Double_t mParticle;
	Double_t rawYield, rawYieldError;
	Double_t invYield, invYieldError;
	Double_t pt, yieldNorm;

	Int_t fitStatus;

	//Used to seed the current mean for kaons and protons in the TPC Spectra
	Double_t prevAmp, prevMean, prevWidth;  

	//Create the Fit Function
	//TF1 *fitFunction = new TF1("studentT",studentT,-1,5,3);
	TF1 *fitFunction;
	TF1 *singleGaussian = new TF1("Gaussian","gaus(0)",0.3,5.0);
	TF1 *doubleGaussian = new TF1("DblGaussian","gaus(0)+gaus(3)",0.3,5.0);

	TF1 *drawPeak = new TF1("SglGaussian","gaus");


	//Set Starting Locations for Each particle for TPC
//	Double_t kaonStart = 1.75;

	Int_t startptBin, endptBin;

	//Loop over the Particle Species
	for (Int_t i=0; i<6; i++){
	
		cout <<"Doing Particle Number: " <<i <<" Histogram: " <<endl;//<<histoNameArray[i][j]<<endl; 

		//Set the Starting Bin Number
		if(i<=3){ //pions and kaons
			if (TPC_TOF == 1){
				startptBin = 3;
				endptBin = 29;
			}
			else if (TPC_TOF == 0){
				startptBin = 2; //start with histogram3 for pions/kaons
				endptBin = 19;
			}
		}
		else if(i>=4) { //p's
			if(TPC_TOF == 1){
				startptBin = 5;
				endptBin = 29;
			}
			else if (TPC_TOF == 0){
				startptBin = 3;
				endptBin = 19;
			}
		}

		//In the majority of cases a single Gaussian will be used as the fit function
		fitFunction = singleGaussian;
		fitFunction->SetParNames(particle[i]+" Amp",particle[i]+" Mean",particle[i]+" Width");

		//Define the Spectra Graph
		spectraGraph[TPC_TOF][CENTBIN][i] = new TGraphErrors();

		//Reset fitStatus
		fitStatus = 1; // TH1::Fit returns 0 if fit OK

		for (Int_t j=startptBin; j<endptBin; j++){

			//Get the Histogram from the File
			htemp = (TH1F *)histoFile->Get(histoNameArray[i][j]);

			//Set The Histogram Title
	//		htemp->SetTitle(histoNameArray[i][j]);

			//Label axes
			htemp->GetYaxis()->SetTitle("counts");
			if(TPC_TOF == 1){
				htemp->GetXaxis()->SetTitle("1/#beta");
			}else if(TPC_TOF == 0){
				htemp->GetXaxis()->SetTitle("ln(dE/dx)");
			}

			//Set the X axis Range
			if (TPC_TOF == 1) {
				if(i<=1) htemp->GetXaxis()->SetRangeUser(0.8,1.5);
				else if(i<=3) htemp->GetXaxis()->SetRangeUser(0.8,2.1);
				else if(i>3)  htemp->GetXaxis()->SetRangeUser(0.8,3.0);
			}
			else if (TPC_TOF ==0) {
				if(i<=1) htemp->GetXaxis()->SetRangeUser(0.4,3.0);
				else if(i>1) htemp->GetXaxis()->SetRangeUser(0.4,4.0);
			}

			//Find the Peaks in the Histogram
			if (TPC_TOF == 1)
				nFoundPeaks = tempSpectrum->Search(htemp,2.,"nobackground",0.050);
			else if (TPC_TOF == 0)
				nFoundPeaks = tempSpectrum->Search(htemp,1.,"nobackground",0.0050);

			cout <<"Found this many peaks: " <<nFoundPeaks <<endl;

			//Get X position of the peaks
			tempPeakLocation = tempSpectrum->GetPositionX();

			//The TSpectrum::GetPositionX() function ranks the peaks in order of 
			//decreasing amplitude, so they need to be reorded by X-location
			peakIndex = new Int_t[nFoundPeaks];
			peakLocation = new Float_t[nFoundPeaks];
			TMath::Sort(nFoundPeaks,tempPeakLocation,peakIndex,false);
			for (Int_t k=0; k<nFoundPeaks; k++){
				peakLocation[k] = tempPeakLocation[peakIndex[k]];
			}

			//Find the Bin Contents of the bin with the peaks
			peakBin = new Int_t[nFoundPeaks];
			peakAmplitude = new Float_t[nFoundPeaks];
			for (Int_t l=0; l<nFoundPeaks; l++){
				peakBin[l] = htemp->GetXaxis()->FindBin(peakLocation[l]);
				peakAmplitude[l] = htemp->GetBinContent(peakBin[l]);
			}

			//Print Peak Information
			if (1){
				cout <<"Peak 1: " << peakLocation[0] <<endl
					 <<"Peak 2: " << peakLocation[1] <<endl
					 <<"Peak 3: " << peakLocation[2] <<endl;
			}


			//Set the Particle Mass
			//0 = pion  |  1 = kaon  |  2 = proton	
			if (i == 0 || i == 1)
				mParticle = mPion;
			else if (i == 2 || i == 3)
				mParticle = mKaon;
			else if (i == 4 || i == 5)
				mParticle = mProton;

			//Set the pt location - used in the TGraphs
			Double_t pt = j*ptBinSize+.5*ptBinSize;
			cout << "Doing pt bin "<< j <<", pt=" << pt << endl;

			//Determine which peak should be fit and what the particle mass is
			//Int_t nPeak;

			//--------------------------------------------------------------------------------
			//TOF  -  SETTING UP THE FIT
			//--------------------------------------------------------------------------------
			if (TPC_TOF == 1){

				//Set expected 1/beta
				betaExpected[i][j]=TMath::Sqrt((mParticle*mParticle)/(pt*pt) + 1);
			
				//Determine Which Peak to Fit
				//0 = pion  |  1 = kaon  |  2 = proton	
				if (i == 0 || i == 1)
					nPeak = 0;
				if (i == 2 || i == 3)
					nPeak = 1;
				if (i == 4 || i == 5)
					nPeak = 2;

				//Initialize the Fit Parameters and set their bounds
				fitFunction->SetParameters(peakAmplitude[nPeak],peakLocation[nPeak],.01);
				//fitFunction->SetParLimits(0,peakAmplitude[nPeak]-500.,peakAmplitude[nPeak]+500.);
				fitFunction->SetParLimits(1,peakLocation[nPeak]-.03,peakLocation[nPeak]+.03);
				fitFunction->SetParLimits(2,0,.05);

				cout << "Peak to be fit found at: " << peakLocation[nPeak] << endl;

				//Check to make sure the found peak's mean is close to the previous peak's mean
				if (j != startptBin){

					//If the found peak's mean is too different from the previous mean the previous mean is used
					//as the seed for the current fit, in addition set the current peak properties to previous peak's
					if (fitFunction->GetParameter(1) > prevMean+.1 || fitFunction->GetParameter(1) < prevMean-.1){
						fitFunction->SetParameter(1,prevMean);
						//fitFunction->SetParameter(2,prevWidth);
						peakAmplitude[nPeak] = prevAmp;
						peakLocation[nPeak] = prevMean;	    
				cout << "Peak found is too different from previous mean. Use prev at:" << peakLocation[nPeak] << endl;
					}

				}

				//Check found peak against betaExpected
				//if(peakLocation[nPeak])

				if(fitFunction->GetParameter(1) > betaExpected[i][j]+0.08 || fitFunction->GetParameter(1) < betaExpected[i][j]-0.08) {
					fitFunction->SetParameter(1, betaExpected[i][j]);
					peakLocation[nPeak] = betaExpected[i][j];	    
				
					//Get estimate of amplitude: bin contents at expected beta
					int expbin = htemp->GetXaxis()->FindBin(betaExpected[i][j]);
					peakAmplitude[nPeak] = htemp->GetBinContent(expbin);
				cout << "Peak found is too different from expected. Use expected at:" << peakLocation[nPeak] << endl;
				}

				//Finally, set the Limits on the parameter based on what peak was chosen
				fitFunction->SetParLimits(1,peakLocation[nPeak]-.2,peakLocation[nPeak]+.2);
				fitFunction->SetParLimits(2,0.002,0.1);
				if( i<2 || ((i==2||i==3) && j>=10) || ((i==4||i==5)&&j>=15) ){
					fitFunction->SetParLimits(2,0.002,0.02);
				}

				//Set the Fit Range
				fitFunction->SetRange(peakLocation[nPeak]-.1,peakLocation[nPeak]+.1);
				
				//The Kaons sometimes need to be fit with a double gaussian to prevent them
				//merging into the pion peak
				if ((i == 2 || i == 3) && j >= 12){

					//Change the Fit Function to a double Gaussian
					fitFunction = doubleGaussian;

					//Set the Parameter Names
					fitFunction->SetParNames(particle[i-2]+" Amp.",particle[i-2]+" Mean",particle[i-2]+" Width",
							particle[i]+" Amp",particle[i]+" Mean",particle[i]+" Width");

					fitFunction->SetParameter(1,peakLocation[0]);
					fitFunction->SetParameter(2,peakAmplitude[0]);

					//Set the Parameter Limits of the Pion Peak
					fitFunction->SetParLimits(1,peakLocation[0]-0.01,peakLocation[0]+0.01);
					fitFunction->SetParLimits(2,0.006,0.011);

					//Set the Parameters of the Kaon Peak
					fitFunction->SetParameter(3,prevAmp);
					fitFunction->SetParameter(4,prevMean);
					fitFunction->FixParameter(5,prevWidth);

					//Set the Parameter Limits of the Kaon Peak
					fitFunction->SetParLimits(3,prevAmp-100,prevAmp+10);
					fitFunction->SetParLimits(4,prevMean-.1,prevMean+.1);
					//fitFunction->SetParLimits(5,prevWidth-.01,prevWidth+.01);

					//Set the Range of the Fit
					fitFunction->SetRange(fitFunction->GetParameter(1)-.05,fitFunction->GetParameter(4)+.04);

				}




			}

			//--------------------------------------------------------------------------------
			//TPC  -  SETTING UP THE FIT
			//--------------------------------------------------------------------------------
			else if (TPC_TOF == 0){

				//Set expected dedx
				if (i == 0 || i == 1)
					dedxExpected[i][j]=log(siPion->Eval(pt));
				else if (i == 2 || i == 3)
					dedxExpected[i][j]=log(siKaon->Eval(pt));
				else if (i == 4 || i == 5)
					dedxExpected[i][j]=log(siProt->Eval(pt));

				//If the Peak Finder finds 3 or four peaks
				if (nFoundPeaks == 3 || nFoundPeaks == 4){

					//Ocassionally an electron peak is found so it needs to be skipped
					if (nFoundPeaks == 4){

						//0 = pion  |  2 = kaon  |  3 = proton	
						if (i == 0 || i == 1)
							nPeak = 0;
						if (i == 2 || i == 3)
							nPeak = 2;
						if (i == 4 || i == 5)
							nPeak = 3;
					}

					if (nFoundPeaks == 3){

						//0 = pion  |  1 = kaon  |  2 = proton	
						if (i == 0 || i == 1)
							nPeak = 0;
						if (i == 2 || i == 3)
							nPeak = 1;
						if (i == 4 || i == 5)
							nPeak = 2;	  
					}
				
					//Set the seed Parameters for the Fit
					fitFunction->SetParameters(peakAmplitude[nPeak],peakLocation[nPeak],.01);

					//Check to make sure the found peak's mean is close to the previous peak's mean
					/*DELETED*/

					//Check found peak against dedxExpected
					if(fitFunction->GetParameter(1) > dedxExpected[i][j]+0.3 || fitFunction->GetParameter(1) < dedxExpected[i][j]-0.3) {
						fitFunction->SetParameter(1, dedxExpected[i][j]);
						peakLocation[nPeak] = dedxExpected[i][j];	    

						//Get estimate of amplitude: bin contents at expected beta
						int expbin = htemp->GetXaxis()->FindBin(betaExpected[i][j]);
						peakAmplitude[nPeak] = htemp->GetBinContent(expbin);
					}

					//Finally, set the Limits on the parameter based on what peak was chosen
					fitFunction->SetParLimits(1,peakLocation[nPeak]-.2,peakLocation[nPeak]+.2);
					fitFunction->SetParLimits(2,0,.16);

					cout << "Fitting peak found at: " << peakLocation[nPeak] << endl;
				}

				//if there are less than three found peaks just use the previous peak as the seed
				else if(nFoundPeaks < 3 && fitStatus==0){ //if there was a successful previous fit
					fitFunction->SetParameters(prevAmp,prevMean,prevWidth);
					fitFunction->SetParLimits(1,prevMean-.2,prevMean+.2);
					fitFunction->SetParLimits(2,0,.16);
				cout << "Using previous peak: Fitting peak found at: " << peakLocation[nPeak] << endl;
				}
				else if(nFoundPeaks < 3 && fitStatus!=0){ //if no previous fit, use Bichsel
					fitFunction->SetParameters(10e5,dedxExpected[i][j],0.1*dedxExpected[i][j]);
					fitFunction->SetParLimits(1,0.8*dedxExpected[i][j],1.1*dedxExpected[i][j]);
					fitFunction->SetParLimits(2,0.05,.16);
				cout << "Using Bichsel expected: " << dedxExpected[i][j] << "Fitting peak found at: " << peakLocation[nPeak] << endl;
				}

				//Set the Function Fit Range
				fitFunction->SetRange(fitFunction->GetParameter(1)-.5,fitFunction->GetParameter(1)+.5);

			}

			//The Kaons and Proton will Benifit from limiting the range of the fit based on local minimums
/*			
				//Look for a local Minimum in the range of the fit
*/
			//Compute the Average kaon width if j>10 (to be used to fix the mean)
/*			DELETED */

			//The Kaons sometimes need to be fit with a double gaussian to prevent them
			//merging into the pion peak
			if (TPC_TOF == 0 && (i == 2 || i == 3) && j >= 10){

				//Change the Fit Function to a double Gaussian
				fitFunction = doubleGaussian;

				//Set the Parameter Names
				fitFunction->SetParNames(particle[i-2]+" Amp.",particle[i-2]+" Mean",particle[i-2]+" Width",
						particle[i]+" Amp",particle[i]+" Mean",particle[i]+" Width");

				//Set the Parameters of the First peak to those of the pion
				/* DELETED */

				fitFunction->SetParameter(1,.9);
				fitFunction->SetParameter(2,.05);

				//Set the Parameter Limits of the Pion Peak
				fitFunction->SetParLimits(1,.8,1.);
				fitFunction->SetParLimits(2,0.01,.15);

				//Set the Parameters of the Kaon Peak
				fitFunction->SetParameter(3,prevAmp);
				fitFunction->SetParameter(4,prevMean);
				//fitFunction->SetParameter(5,.03);
				fitFunction->FixParameter(5,prevWidth);
				//fitFunction->FixParameter(5,kaonWidthAvg);
				cout <<"Using Average Kaon Width" <<endl;

				//Set the Parameter Limits of the Kaon Peak
				fitFunction->SetParLimits(3,prevAmp-100,prevAmp+10);
				fitFunction->SetParLimits(4,prevMean-.1,prevMean+.1);
				//fitFunction->SetParLimits(5,prevWidth-.01,prevWidth+.01);

				if (j >= 14)
					fitFunction->SetParLimits(3,prevAmp-50,prevAmp+5);

				//Set the Range of the Fit
				fitFunction->SetRange(fitFunction->GetParameter(1)-.4,fitFunction->GetParameter(4)+.3);

			}

			//The Anti-Protons can be fit using the means of the protons as a fixed parameter
			/* DELETED */

			//Draw the Histogram
			canvas->cd(1);
			htemp->Draw();

			//Wait if the nFoundPeaks is other than three
			//if (nFoundPeaks != 3)
			//gSystem->Sleep(10000);

			//Fit the Histogram
			fitResult = htemp->Fit(fitFunction,"RS");
			fitStatus = fitResult;

			//Set the Current Mean to the prevMean
			prevAmp  = fitFunction->GetParameter(0);
			prevMean = fitFunction->GetParameter(1);
			prevWidth= fitFunction->GetParameter(2);

			//If the Protons were fit save the proton mean
			if (i == 4){
				protonMean[j] = fitFunction->GetParameter(1);
				protonWidth[j] = fitFunction->GetParameter(2);
			}

			//If the Kaons or Antikaons were fit store their first six widths
//			if (TPC_TOF == 0 && (i == 2 || i == 3) && j < 10){
//				kaonWidth[j-4] = fitFunction->GetParameter(2);
//			}

			//If the Double Gaussian was used...
			if (TPC_TOF == 0 && (i == 2|| i == 3) && j >= 10){
				prevAmp = fitFunction->GetParameter(3);
				prevMean = fitFunction->GetParameter(4);
				prevWidth = fitFunction->GetParameter(5);
			}
			if (TPC_TOF == 1 && (i == 2|| i == 3) && j >= 12){
				prevAmp = fitFunction->GetParameter(3);
				prevMean = fitFunction->GetParameter(4);
				prevWidth = fitFunction->GetParameter(5);
			}

			// **added**
			//Draw the fit peak
			if ( (i == 2|| i == 3)) {
				if ( (TPC_TOF==0 && j >= 10) || (TPC_TOF==1 && j>=12)) {
					drawPeak->SetParameter(0,fitFunction->GetParameter(3));
					drawPeak->SetParameter(1,fitFunction->GetParameter(4));
					drawPeak->SetParameter(2,fitFunction->GetParameter(5));
				}
			}else {
				drawPeak->SetParameter(0,fitFunction->GetParameter(0));
				drawPeak->SetParameter(1,fitFunction->GetParameter(1));
				drawPeak->SetParameter(2,fitFunction->GetParameter(2));
			}
			drawPeak->SetLineColor(kRed);
			drawPeak->Draw("same");

			//Update the Canvas Pad
			gPad->Update();

			//Save Animated GIFS
			if (saveFits){
				
				//Save as an Animated GIF
				TString fitTitle = fitDir+collision+energy+particle[i]+gif+gifSpeed;
				canvas->SaveAs(fitTitle);

				//Save as .png
				fitTitle = fitDir+collision+energy+histoNameArray[i][j]+png;
				canvas->SaveAs(fitTitle);
			}



			//Look at a particular Histogram
			//if (TPC_TOF == 0 && /*(CENTBIN == 0|| CENTBIN == 1)*/CENTBIN ==5 && (i == 3) && j ==4 )
			//if (TPC_TOF == 1 && i == 5 && CENTBIN < 3)
			//if (TPC_TOF == 0 && (i == 2 || i==3))
			//gSystem->Sleep(10000);

			// //Create the TGraphErrors to hold the Fit Results
			// fitResultGraph[i][0] = new TGraphErrors();
			// fitResultGraph[i][1] = new TGraphErrors();
			// fitResultGraph[i][2] = new TGraphErrors();
			// fitResultGraph[i][3] = new TGraphErrors();

			// //Fill the Fit Result Graphs

			// fitResultGraph[i][0]->SetPoint(j,pointX,fitFunction->GetParameter(0));
			// fitResultGraph[i][0]->SetPointError(j,mTm0BinSize/2.,fitFunction->GetParError(0));
			// fitResultGraph[i][1]->SetPoint(j,pointX,fitFunction->GetParameter(1));
			// fitResultGraph[i][1]->SetPointError(j,mTm0BinSize/2.,fitFunction->GetParError(1));
			// fitResultGraph[i][2]->SetPoint(j,pointX,fitFunction->GetParameter(2));
			// fitResultGraph[i][2]->SetPointError(j,mTm0BinSize/2.,fitFunction->GetParError(2));
			// fitResultGraph[i][3]->SetPoint(j,pointX,fitFunction->GetChisquare()/fitFunction->GetNDF());

			//Draw the fit Result Graphs and Update the Pad
			// canvas1->cd(1);
			// fitResultGraph[i][0]->Draw("AL*");
			// canvas1->cd(2);
			// fitResultGraph[i][1]->Draw("AL*");
			// canvas1->cd(3);
			// fitResultGraph[i][2]->Draw("AL*");
			// canvas1->cd(4);
			// fitResultGraph[i][3]->Draw("AL*");

			// gPad->Update();


			//Extract Raw Yields from Histograms
			// rawYield = fitFunction->Integral(-10,10)/htemp->GetBinWidth(5);
			// rawYieldError = fitFunction->IntegralError(-10,10,fitResult->GetParams(),
			//       fitResult->GetCovarianceMatrix().GetMatrixArray())/htemp->GetBinWidth(5);

			Double_t Amplitude = fitFunction->GetParameter(0);
			Double_t Width     = fitFunction->GetParameter(2);
			Double_t AmplitudeErr = fitFunction->GetParError(0);
			Double_t WidthErr  = fitFunction->GetParError(2);

			//If the Double Gaussian is used the second peak's properties must be used
			if (TPC_TOF == 0 && (i == 2 || i == 3) && j >= 10){
				Amplitude = fitFunction->GetParameter(3);
				Width = fitFunction->GetParameter(5);
				AmplitudeErr = fitFunction->GetParError(3);
				WidthErr = fitFunction->GetParError(5);
			}

			rawYield = Amplitude*Width*TMath::Sqrt(2*TMath::Pi())/htemp->GetBinWidth(5);
			rawYieldError = rawYield*TMath::Sqrt(AmplitudeErr*AmplitudeErr/(Amplitude*Amplitude)+(WidthErr*WidthErr/(Width*Width)));

			//Correction Factor to bring Sams Data to parity with Expected Results
			Double_t CorrFact=1.;// = 100.;

			//Define Normalization
			yieldNorm = (1./pt)*(1./ptBinSize)*(1./yBinSize)*(1./nGoodEvents)*(1./(2.*TMath::Pi()))*(1./CorrFact);

			//Compute the Invariant Yields
			invYield = rawYield*yieldNorm;
			invYieldError = rawYieldError*yieldNorm;

			cout <<"-------rawYield=" <<rawYield <<"------------" <<endl;
			cout <<"-------invYield=" <<invYield <<"------------" <<endl;

			//Fill the Spectra Graph
			spectraGraph[TPC_TOF][CENTBIN][i]->SetPoint(j-startptBin, pt, invYield);
			spectraGraph[TPC_TOF][CENTBIN][i]->SetPointError(j-startptBin,ptBinSize/2.,invYieldError);

			//gSystem->Sleep(1000);

			//Delete the pointers used in the loop
			delete peakIndex;
			delete peakLocation;
			delete peakBin;
			delete peakAmplitude;


		}//End Loop Over mTm0Bins

		//Print Status
		cout <<"------- DONE FITTING !!! ---------" <<endl;

	}//End Loop Over Particle Species

}//End spectra_CentBin

void spectra_Fitter_uu(){

	//Load the Root File and the Necessary TTree
//	treeFile = new TFile("v4_ntuples/UU193_NTuples_cent1_completeish.root");
	countNGoodEvents();

	histoFile = new TFile("../v4_ntuples/UU193_histos_cent1_completeish.root");
	
	//Set the Y axis to Log Scale for the Histograms
	canvas->SetLogy();
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(112);

	//Loop over the Centrality Bins and do the TPC Spectra
//	for (Int_t c=8; c>=0; c--){
		spectra_CentBin(0,0);
//	}

	//Loop over the Centrality Bins and do the TOF Spectra
//	for (Int_t c=8; c>=0; c--){
		spectra_CentBin(0,1);
//	}


	//Apply TOF efficiency scaling
	efficiencyCorrection();

	//Write out spectra data to text file
	writeTextFiles();

	//Draw the TOF Spectra
	drawSpectra(1);
				
	//Draw the TPC Spectra
	drawSpectra(0);

	//Make Legends
	//makeLegends(1);

	//Draw Both the Spectra
	drawTPCandTOF();

	//Do the Efficiency Correction
//	efficiencyCorrection();

}
