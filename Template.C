#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TString.h"
#include <TStopwatch.h>
#include <TComplex.h>
#include <vector>
#include "include/toyflowinputs.h"
#include <TSystem.h> 

using namespace std;

//Declaration of histograms and constants
TH1D *hDeltaPhiHM = new TH1D("hDeltaPhiHM","hDeltaPhiHM",200, 0.0, 2.0*TMath::Pi());
TH1D *hDeltaPhiLM = new TH1D("hDeltaPhiLM","hDeltaPhiLM",200, 0.0, 2.0*TMath::Pi());
const int NH=2;
const double Fit_min=1.0;
const double Fit_max=3.0;
//Member functions
void LoadData(TString inputfile);


//Main function that calculates the fit
void Template(){
	//getting histgrams from input files
	LoadData(inputfilename);

	//defining formula for G(fourier)
	TString strformula = "[0]*(1";
	for (Int_t ih=0; ih<NH; ih++){
		strformula += Form("+2*[%d]*TMath::Cos(%d*(x-[%d]))",ih+1,ih+2,NH+ih+1);
	}
	strformula+=")";
	cout<<strformula<<endl;

	// Y'=Y_HM-F*Y_LM
	hDeltaPhiHM->Add(hDeltaPhiLM,-1);// need to get F in there somehow

	//Fit Y'

	//chi^2

	//find best F parameters


}

void LoadData(TString inputfilename){
	TFile *fIn = TFile::Open(inputfilename,"read");

	hDeltaPhiHM = (TH1D*)fIn->Get("hDeltaPhiHM");
	hDeltaPhiLM = (TH1D*)fIn->Get("hDeltaPhiLM"); 	
}