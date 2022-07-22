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
TH1D *hYPrime[F_trials+1];
TF1 *fFitarray[F_trials+1];

const int NH=2;
const double F_min=1.0;
const double F_max=3.0;
const int F_trials = 100;

for( int iF=0; iF<=F_trials; iF++){
		hYPrime[iF] = new TH1D(Form("hYPrime%03d"),Form("F=%d",Farray[iF]),200,0.0, 2.0*TMath::Pi());
	}
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

	//defining fit function + setting parameters
	TString ParName[NH+1]={"const","v_2","v_3"};
	TF1 *fFit = new TF1("fFit", strformula, 0, 2.0*TMath::Pi());
	fFit->SetParameter(0,1E4);
	for (Int_t i=1; i<=NH; i++) fFit->SetParameter(i,0.06);
	for (Int_t i=0; i<=NH+1; i++) fFit->SetParName(i,ParName[i]);

	//defining the different values for F
	double F_stepsize = (F_max-F_min)/F_trials;
	double Farray[F_trials+1]={};
	double F = F_min;
	Farray[0]=F;
	for( int iF=1; iF<=F_trials; iF++){
		F=+F_stepsize;
		Farray[iF]=F;
	}

	for( int iF=0; iF<=F_trials; iF++){
		//Multiplying Y_LM with F so that Add function can be used later for substraction.
		hDeltaPhiLM->Scale(Farray[iF]);
		// Y'=Y_HM-F*Y_LM
		hYPrime[iF]=hDeltaPhiHM->Add(hDeltaPhiLM,-1);
		//Fitting Y'
		hYPrime[iF]->Fit("fFit");
		fFitarray[iF]=fFit;
	}
	
	//chi^2

	//find best F parameters


}

void LoadData(TString inputfilename){
	TFile *fIn = TFile::Open(inputfilename,"read");

	hDeltaPhiHM = (TH1D*)fIn->Get("hDeltaPhiHM");
	hDeltaPhiLM = (TH1D*)fIn->Get("hDeltaPhiLM"); 	
}