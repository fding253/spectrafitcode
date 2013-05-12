//Loads the histograms from Sam's root file and fits them with a 
//Student T distribution.  Then computes the particle yield and
//finally draws the spectra

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;


#include "histogramList_TOF.C"
#include "histogramList_TPC.C"
#include "spectra_Fitter_196.h"
#include "UtilityFunctions.h"
//#include "histoNameArray.C"

//Create Canvas for Drawing the TOF Histograms
TCanvas *canvas = new TCanvas("canvas");


//Define Bin Property Values
Double_t mTm0BinSize = .025;
Double_t yBinSize = .1;

//Run Time Options
Bool_t saveFits = false;
TString fitDir = "FitGIFS/";
TString gifSpeed = "+50";

void spectra_CentBin(Int_t CENTBIN,Int_t TPC_TOF){

  cout <<"Doing Centrality Bin: " <<CENTBIN <<endl;

  //Define the Histogram Name Array
  setHistoNameArray(CENTBIN,TPC_TOF);

  //Create a Pointer to hold the histogram
  TH1F *htemp = new TH1F();

  //Create Objects needed for Peak Finding
  TSpectrum *tempSpectrum;
  if (TPC_TOF == 1)
    tempSpectrum = new TSpectrum(3,.25);
  if (TPC_TOF == 0)
    tempSpectrum = new TSpectrum(4,.25);

  //Define Variables needed for Peak Fitting
  Int_t nFoundPeaks, nPeak;
  Int_t *peakIndex, *peakBin; 
  Float_t *tempPeakLocation, *peakLocation, *peakAmplitude;
  Double_t protonMean[40];  //This is used to fix the means of the antiprotons
  Double_t protonWidth[40]; //This is used to fix the means of the antiprotons
  Double_t kaonWidth[6]; //This is used to store the first six widths of the kaon (and antikaon) to fix the remainer of the widths

  //Create Objects needed for Yeild Extraction
  TFitResultPtr fitResult;
  Double_t mParticle;
  Double_t rawYield, rawYieldError;
  Double_t invYield, invYieldError;
  Double_t mT, yieldNorm;

  //Used to seed the current mean for kaons and protons in the TPC Spectra
  Double_t prevAmp, prevMean, prevWidth;  
  
  //Create the Fit Function
  //TF1 *fitFunction = new TF1("studentT",studentT,-1,5,3);
  TF1 *fitFunction;
  TF1 *singleGaussian = new TF1("Gaussian","gaus(0)",-0.4,0.4);
  TF1 *doubleGaussian = new TF1("DblGaussian","gaus(0)+gaus(3)",-0.4,.4);



  //Set the Starting Bin Number
  Int_t startmTm0Bin, endmTm0Bin;
  if (TPC_TOF == 1){
    startmTm0Bin = 15;
    endmTm0Bin = 39;
  }
  else if (TPC_TOF == 0){
    startmTm0Bin = 4;
    endmTm0Bin = 20;
  }

  //Set Starting Locations for Each particle for TPC
  Double_t kaonStart = 1.75;


  //Loop over the Particle Species
  for (Int_t i=0; i<6; i++){

    cout <<"Doing Particle Number: " <<i <<" Histogram: " <<endl;//<<histoNameArray[i][j]<<endl; 

    //In the majority of cases a single Gaussian will be used as the fit function
    fitFunction = singleGaussian;
    fitFunction->SetParNames(particle[i]+" Amp",particle[i]+" Mean",particle[i]+" Width");

    //Define the Spectra Graph
    spectraGraph[TPC_TOF][CENTBIN][i] = new TGraphErrors();

    for (Int_t j=startmTm0Bin; j<endmTm0Bin; j++){
      
      //Get the Histogram from the File
      htemp = (TH1F *)inputFile->Get(histoNameArray[i][j]);
      
      //Set The Histogram Title
      htemp->SetTitle(histoNameArray[i][j]);
      
      //Set the X axis Range
      if (TPC_TOF == 1)
	htemp->GetXaxis()->SetRangeUser(-1.00,1.25);
      else if (TPC_TOF ==0)
	htemp->GetXaxis()->SetRangeUser(0,5.0);

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
      if (0){
	cout <<"Peak 1: " <<peakLocation[0] <<endl
	     <<"Peak 2: " <<peakLocation[1] <<endl
	     <<"Peak 3: " <<peakLocation[2] <<endl;
      }


      //Set the Particle Mass
      //0 = pion  |  1 = kaon  |  2 = proton	
      if (i == 0 || i == 1)
	mParticle = mPion;
      else if (i == 2 || i == 3)
	mParticle = mKaon;
      else if (i == 4 || i == 5)
	mParticle = mProton;
      

      //Determine which peak should be fit and what the particle mass is
      //Int_t nPeak;
      
      //--------------------------------------------------------------------------------
      //TOF  -  SETTING UP THE FIT
      //--------------------------------------------------------------------------------
      if (TPC_TOF == 1){
	
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

	//Check to make sure the found peak's mean is close to the previous peak's mean
	if (j != startmTm0Bin){
	  
	  //If the found peak's mean is too different from the previous mean the previous mean is used
	  //as the seed for the current fit, in addition set the current peak properties to previous peak's
	  if (fitFunction->GetParameter(1) > prevMean+.3 || fitFunction->GetParameter(1) < prevMean-.3){
	    fitFunction->SetParameter(1,prevMean);
	    //fitFunction->SetParameter(2,prevWidth);
	    peakAmplitude[nPeak] = prevAmp;
	    peakLocation[nPeak] = prevMean;	    
	  }
	}
	
	//Finally, set the Limits on the parameter based on what peak was chosen
	fitFunction->SetParLimits(1,peakLocation[nPeak]-.2,peakLocation[nPeak]+.2);
	fitFunction->SetParLimits(2,0,.1);
	
	
	//Set the Fit Range
	fitFunction->SetRange(peakLocation[nPeak]-.1,peakLocation[nPeak]+.1);
	
    
      }

      //--------------------------------------------------------------------------------
      //TPC  -  SETTING UP THE FIT
      //--------------------------------------------------------------------------------
      else if (TPC_TOF == 0){
	
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
	  if (j != startmTm0Bin){
	    
	    //If the found peak's mean is too different from the previous mean the previous mean is used
	    //as the seed for the current fit, in addition set the current peak properties to previous peak's
	    if (fitFunction->GetParameter(1) > prevMean+.3 || fitFunction->GetParameter(1) < prevMean-.3){
	      fitFunction->SetParameter(1,prevMean);
	      //fitFunction->SetParameter(2,prevWidth);
	      peakAmplitude[nPeak] = prevAmp;
	      peakLocation[nPeak] = prevMean;	    
	    }
	  }

	  else if (j == startmTm0Bin){

	    //Kaons
	    if (i == 2){
	      fitFunction->SetParameter(1,kaonStart);
	      peakLocation[nPeak] = kaonStart;
	    }	  
	    //Anti-Kaons
	    if (i == 3){
	      fitFunction->SetParameter(1,kaonStart+.1);
	      peakLocation[nPeak] = kaonStart+.1;
	    }	  


	  }

	  //Finally, set the Limits on the parameter based on what peak was chosen
	  fitFunction->SetParLimits(1,peakLocation[nPeak]-.2,peakLocation[nPeak]+.2);
	  fitFunction->SetParLimits(2,0,.1);
	}
	
	//if there are less than three found peaks just use the previous peak as the seed
	else if(nFoundPeaks < 3){
	  fitFunction->SetParameters(prevAmp,prevMean,prevWidth);
	  fitFunction->SetParLimits(1,prevMean-.2,prevMean+.2);
	  fitFunction->SetParLimits(2,0,.1);
	}

	//Set the Function Fit Range
	fitFunction->SetRange(fitFunction->GetParameter(1)-.5,fitFunction->GetParameter(1)+.5);

      }

      //The Kaons and Proton will Benifit from limiting the range of the fit based on local minimums
      if (i >= 2){

	//Look for a local Minimum in the range of the fit
	Double_t lowerBound, upperBound, minValue, lookLower;
	fitFunction->GetRange(lowerBound,upperBound);
	
	if (TPC_TOF == 1)
	  lookLower = .05;
	else if (TPC_TOF == 0)
	  lookLower = .3;

	Int_t minBin = findLocalMinBin(htemp,fitFunction->GetParameter(1)-lookLower,fitFunction->GetParameter(1));
	
	if (minBin >= 0){
	  minValue = htemp->GetBinCenter(minBin);
	  
	  cout <<"Found Local Min Bin at: " <<minValue <<endl;
	  
	  //If a local min is found limit the range of the fit accordingly
	  if (minValue >= fitFunction->GetParameter(1))
	    fitFunction->SetRange(lowerBound,minValue);
	  else if (minValue < fitFunction->GetParameter(1))
	    fitFunction->SetRange(minValue,upperBound);
	}
	
	Double_t temp1,temp2;
	fitFunction->GetRange(temp1,temp2);
	cout <<"Low Bound: " <<temp1 <<" High Bound: " <<temp2 <<endl;
      }

      //Compute the Average kaon width if j>10 (to be used to fix the mean)
      Double_t kaonWidthAvg = 0;
      if (j >= 10 && (i == 2 || i ==3)){
	
	for (Int_t p=0; p<6; p++)
	  kaonWidthAvg += kaonWidth[p];

	kaonWidthAvg = kaonWidthAvg/6;
	
      }

      //The Kaons sometimes need to be fit with a double gaussian to prevent them
      //merging into the pion peak
      if (TPC_TOF == 0 && (i == 2 || i == 3) && j >= 10){

      	//Change the Fit Function to a double Gaussian
	fitFunction = doubleGaussian;

	//Set the Parameter Names
	fitFunction->SetParNames(particle[i-2]+" Amp.",particle[i-2]+" Mean",particle[i-2]+" Width",
				 particle[i]+" Amp",particle[i]+" Mean",particle[i]+" Width");

      	//Set the Parameters of the First peak to those of the pion
	if (CENTBIN < 7)
	  fitFunction->SetParameter(0,10e4);
	else 
	  fitFunction->SetParameter(0,10e3);

	fitFunction->SetParameter(1,.9);
	fitFunction->SetParameter(2,.05);

	//Set the Parameter Limits of the Pion Peak
	fitFunction->SetParLimits(1,.8,1.);
	fitFunction->SetParLimits(2,0,.1);
	
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
      if (TPC_TOF == 1 && i == 5){

	fitFunction->FixParameter(1,protonMean[j]);

	if (j >= 25)
	  fitFunction->FixParameter(2,protonWidth[j]);

	if (CENTBIN == 0 && j>23)
	  fitFunction->SetParLimits(0,prevAmp-2,prevAmp+.5);

      }


      //Draw the Histogram
      canvas->cd(1);
      htemp->Draw();

      //Wait if the nFoundPeaks is other than three
      //if (nFoundPeaks != 3)
      //gSystem->Sleep(10000);

      //Fit the Histogram
      fitResult = htemp->Fit(fitFunction,"RS");

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
      if (TPC_TOF == 0 && (i == 2 || i == 3) && j < 10){
	kaonWidth[j-4] = fitFunction->GetParameter(2);
      }

      //If the Double Gaussian was used...
      if (TPC_TOF == 0 && (i == 2|| i == 3) && j >= 10){
	prevAmp = fitFunction->GetParameter(3);
	prevMean = fitFunction->GetParameter(4);
	prevWidth = fitFunction->GetParameter(5);
      }

      //Update the Canvas Pad
      gPad->Update();

      //Save Animated GIFS
      if (saveFits){
	
	//Create File Title
	TString fitTitle = fitDir+collision+energy+particle[i]+CentBin+nCentBin[CENTBIN]+gif+gifSpeed;

	//Save as an Animated GIF
	canvas->SaveAs(fitTitle);

      }


      //Set the mTm0 location - used in the TGraphs
      Double_t pointX = j*mTm0BinSize+.5*mTm0BinSize;

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

      rawYield = Amplitude*Width*sqrt(2*3.14159)/htemp->GetBinWidth(5);
      rawYieldError = rawYield*sqrt(pow(AmplitudeErr/Amplitude,2)+pow(WidthErr/Width,2));

      //Calculate mT
      mT = j*mTm0BinSize+.5*mTm0BinSize + mParticle;

      //Correction Factor to bring Sams Data to parity with Expected Results
      Double_t CorrFact = 100.;

      //Define Normalization
      yieldNorm = (1./mT)*(1./mTm0BinSize)*(1./yBinSize)*(1./nCentEvents[CENTBIN])*(1./(2.*3.14159))*(1./CorrFact);

      //Compute the Invariant Yields
      invYield = rawYield*yieldNorm;
      invYieldError = rawYieldError*yieldNorm;

      cout <<"-------" <<rawYield <<"------------" <<endl;
      cout <<"-------" <<invYield <<"------------" <<endl;

      //Fill the Spectra Graph
      spectraGraph[TPC_TOF][CENTBIN][i]->SetPoint(j-startmTm0Bin, pointX, invYield);
      spectraGraph[TPC_TOF][CENTBIN][i]->SetPointError(j-startmTm0Bin,mTm0BinSize/2.,invYieldError);
 
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

void spectra_Fitter_196(){

  //Load the Root File and the Necessary TTree
  inputFile = new TFile("AuAu_19_6_Complete.root");
  countNCentEvents();

  //Set the Y axis to Log Scale for the Histograms
  canvas->SetLogy();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(112);

  //Loop over the Centrality Bins and do the TPC Spectra
  for (Int_t c=8; c>=0; c--){
    spectra_CentBin(c,0);
  }

  //Loop over the Centrality Bins and do the TOF Spectra
  for (Int_t c=8; c>=0; c--){
    spectra_CentBin(c,1);
  }
  


  //Draw the TOF Spectra
  drawSpectra(1);

  //Draw the TPC Spectra
  drawSpectra(0);

  //Make Legends
  //makeLegends(1);

  //Draw Both the Spectra
  drawTPCandTOF();

  //Do the Efficiency Correction
  efficiencyCorrection();

}
