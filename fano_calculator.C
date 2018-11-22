#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include "TFile.h"
#include "TObject.h"
#include "TKey.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLegend.h"
#include <TROOT.h>
#include "TSystem.h"
#include "TGaxis.h"
#include <cstdlib>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjString.h"
#include "TH2D.h"
using namespace std;

void Enable_and_Set_Branches(TTree* & tree);

// Setting parameters //////////////////////////////////////////////////

  // range for the number of electrons per cluster
  int emin =1540 ; int emax = 1620 ; //1530 1610
  int ohdu_numer = 4;
  //number of bins to take into account for chi2
  int bines = 20;
  float minePix = 0; // will be clasified as 1e-

  ////////////////////////////////////////////////////////////////////////

  int Entries_mc = 1;
  int Entries_exp = 1;
  //int xmin = 0; // range for histograms
  //int xmax =0; // range for histograms
  //int xBary_min=0;int xBary_max=100;
  //int yBary_min=0;int yBary_max=100;

  int nbins = 10;

  const int maxClusterSize = 50000;

////////////////////////////////////////////////////////////////////////

  int runID; int ohdu; int expoStart;
  int nSat; int flag;
  int xMin; int xMax;
  int yMin; int yMax;
  Float_t e; Float_t n;
  Float_t xBary; Float_t yBary;
  Float_t xVar; Float_t yVar;
  int  nSavedPix;
  int xPix[maxClusterSize];
  int yPix[maxClusterSize];
  Float_t ePix[maxClusterSize];

using namespace std;

////////////////////////////////////////////////////////////////////////
//  Older compilers
////////////////////////////////////////////////////////////////////////
string itos(const int i){
	ostringstream ossAux;
	ossAux << i;
	return ossAux.str();
}



void fano_calculator(){
//int main(int argc, char* argv[]){

// Experimental Data ///////////////////////////////////////////////////
// Get input files//////////////////////////////////////////////////////

cout<<"min ePix: "<< minePix<<endl;

					ofstream myfile;



					//myfile.open ("/home/mariano/MEGAsync/images_from_mkids/analisis/fe55_ST136/fano.txt");
					myfile.open ("/home/mariano/MEGAsync/images_from_mkids/analisis/Con_LED/142K/fano.txt");


            //TFile * f_exp = TFile::Open("./55Fe_exp.root");
            TFile * f_exp = TFile::Open("/home/mariano/MEGAsync/images_from_mkids/analisis/Con_LED/142K/output_2.root");
            if (!f_exp->IsOpen()) {std::cerr << "ERROR: cannot open the root file with experimental data" << std::endl;}
            if (!f_exp->IsOpen()) {std::cerr << "ERROR: cannot open the root file with experimental data" << std::endl;}
            TTree * texp = (TTree*) f_exp->Get("hitSumm");



            //for the total of clusters
            TH1D * h_exp_e_total  =  new TH1D("h_exp_e_total","", bines, emin, emax);
            h_exp_e_total -> Sumw2();

            //for each n
            TH1D * h_exp_e[nbins+1];
            TString histname_tmp_exp;
            for(int ystarbin=0;ystarbin<nbins+1;++ystarbin){
            histname_tmp_exp  = "h_exp_e";
            histname_tmp_exp += ystarbin;
            h_exp_e[ystarbin] = new TH1D(histname_tmp_exp.Data(),"",bines,emin,emax);
            h_exp_e[ystarbin]->Sumw2();
            }

            // Get information from trees///////////////////////////////////////////
            int Entries_exp = texp -> GetEntries();
            cout<<"Entries in experimental data file: "<<Entries_exp<<endl;
            Enable_and_Set_Branches(texp);

          	vector<double> events_exp(50000);

    				vector<double> mean_exp_fit(nbins+1);
    				vector<double> sigma_exp_fit(nbins+1);
    				vector<double> fano_exp_fit(nbins+1);
            vector<double> events_exp_fit(50000);

    				vector<double> mean_exp_fit_error(nbins+1);
    				vector<double> sigma_exp_fit_error(nbins+1);
    				vector<double> fano_exp_fit_error(nbins+1);


				for (int ene=1; ene<nbins+1; ene++){
					int i = 0;
					for(int i_event=0;i_event<Entries_exp; i_event++){
					texp->GetEntry(i_event);

				//			if (ohdu == ohdu_numer) {
								if (e>emin && e<emax){  // number of electrons
                  if (n==ene){
    									h_exp_e[ene]->Fill(e);
    									h_exp_e_total->Fill(e);
    									events_exp[i]=e;
									    i++;
									}
								}
					//		}
						}



						// aca calculamos varianza y media

						h_exp_e[ene]->Fit("gaus");
						TF1 *fit1 = (TF1*)h_exp_e[ene]->GetFunction("gaus");

						TCanvas*  canvitas   = new TCanvas("clusters","clusters",800,500);
						canvitas->cd();
						fit1->Draw("HIST E1");
						h_exp_e[ene]->Draw("HIST E1 same");

						sigma_exp_fit[ene]= fit1->GetParameter(2);
						sigma_exp_fit_error[ene]= fit1->GetParError(2);
						mean_exp_fit[ene]= fit1->GetParameter(1);
						mean_exp_fit_error[ene]= fit1->GetParError(1);
						fano_exp_fit[ene]= pow(sigma_exp_fit[ene],2)/mean_exp_fit[ene];
						fano_exp_fit_error[ene]= pow(pow(2*sigma_exp_fit[ene]*sigma_exp_fit_error[ene]/mean_exp_fit[ene],2)+pow(mean_exp_fit_error[ene]*sigma_exp_fit[ene]*sigma_exp_fit[ene]/(mean_exp_fit[ene]*mean_exp_fit[ene]),2),0.5);
						events_exp_fit[ene]=h_exp_e[ene]->GetEntries();



					}
						// save to file
						myfile << "Mean" << "	" <<"Sigma" << "	" <<"Fano factor" << "	" <<"Fano error" << "	" <<", increasing n" <<endl;
						for(int i = 1; i < nbins+1; i ++) {
						myfile  << mean_exp_fit[i]  << "	" << sigma_exp_fit[i] << "	" << fano_exp_fit[i] << "	" << fano_exp_fit_error[i]  << "	" <<  events_exp_fit[i] << endl;
						}

            h_exp_e_total->Fit("gaus");
            TF1 *fit2 = (TF1*)h_exp_e_total->GetFunction("gaus");


            h_exp_e_total->Draw("HIST E2");
            fit2->Draw("HIST same");

            sigma_exp_fit[0]= fit2->GetParameter(2);
						sigma_exp_fit_error[0]= fit2->GetParError(2);
						mean_exp_fit[0]= fit2->GetParameter(1);
						mean_exp_fit_error[0]= fit2->GetParError(1);
						fano_exp_fit[0]= pow(sigma_exp_fit[0],2)/mean_exp_fit[0];
						fano_exp_fit_error[0]= pow(pow(2*sigma_exp_fit[0]*sigma_exp_fit_error[0]/mean_exp_fit[0],2)+pow(mean_exp_fit_error[0]*sigma_exp_fit[0]*sigma_exp_fit[0]/(mean_exp_fit[0]*mean_exp_fit[0]),2),0.5);
						events_exp_fit[0]=h_exp_e_total->GetEntries();

            myfile  << "Para todos los n" << endl;
            myfile  << mean_exp_fit[0]  << "	" << sigma_exp_fit[0] << "	" << fano_exp_fit[0] << "	" << fano_exp_fit_error[0]  << "	" <<  events_exp_fit[0] << endl;
	    
	cout<<"Fano: "<<fano_exp_fit[0]<<endl;
	cout<<"Error: "<<fano_exp_fit_error[0]<<endl;
	cout<<"Eventos: "<<events_exp_fit[0]<<endl;

	myfile.close();
//return 0;
}


////////////////////////////////////////////////////////////////////////

void Enable_and_Set_Branches(TTree* & tree){

  tree->SetBranchStatus("*",1); //enable all branches

  tree->SetBranchAddress ("runID",&runID);
  tree->SetBranchAddress ("ohdu",&ohdu);
  tree->SetBranchAddress ("expoStart",&expoStart);
  tree->SetBranchAddress ("nSat",&nSat);
  tree->SetBranchAddress ("flag",&flag);
  tree->SetBranchAddress ("xMin",&xMin);
  tree->SetBranchAddress ("xMax",&xMax);
  tree->SetBranchAddress ("n",&n);
  tree->SetBranchAddress ("yMin",&yMin);
  tree->SetBranchAddress ("yMax",&yMax);
  tree->SetBranchAddress ("e",&e);
  tree->SetBranchAddress ("xBary",&xBary);
  tree->SetBranchAddress ("yBary",&yBary);
  tree->SetBranchAddress ("xVar",&xVar);
  tree->SetBranchAddress ("yVar",&yVar);
  tree->SetBranchAddress ("nSavedPix",&nSavedPix);
  tree->SetBranchAddress ("xPix",&xPix);
  tree->SetBranchAddress ("yPix",&yPix);
  tree->SetBranchAddress ("ePix",&ePix);
}
