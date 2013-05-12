#include <iostream>
#include <histogramList.C>

using namespace std;

TString *histoNameArray[6];


void setHistoNameArray(Int_t centbin){

  cout <<"In setHistoNameArray() centbin is: " <<centbin <<endl;

  if (centbin  == 0){
    histoNameArray[0]=piPlus_TOF_CentBin0;
    histoNameArray[1]=piMinus_TOF_CentBin0;
    histoNameArray[2]=kPlus_TOF_CentBin0;
    histoNameArray[3]=kMinus_TOF_CentBin0;
    histoNameArray[4]=pPlus_TOF_CentBin0;
    histoNameArray[5]=pMinus_TOF_CentBin0;
    }


  cout <<piPlus_TOF_CentBin0 <<endl;
  cout <<histoNameArray[0] <<endl;

  for (Int_t i=0; i<40; i++)
    cout <<histoNameArray[0][i] <<endl;




}
