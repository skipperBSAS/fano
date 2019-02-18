/*

Este script es una modificacion de fano_calculator escrito por Mariano.

Domingo 17 de Febrero de 2019 - Dario Rodrigues

Este script fue escrito para estudiar como evoluciona el numero medio de el pico de 5.9 keV del 55Fe en funcion del tamanio de cluster
Grafica los ajustes hasta n=6 y luego un TGraph con la media en funcion de n, esto se guarda en un archivo.
Esta version toma los datos de un branch que Andre agrego a los rootfiles donde todo ya esta en cargas usando la calibracion no lineal
que se derivo usando un ajuste gaussiano iterativo de los picos correspondientes a #electrones.

Esta version permite quitar columnas malas
*/


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
#include <fstream>
using namespace std;

void Enable_and_Set_Branches(TTree* & tree);

// Setting parameters //////////////////////////////////////////////////

  // range for the number of electrons per cluster
  //int emin =1540 ; int emax = 1620 ; //1530 1610
  int emin =1450 ; int emax = 1850 ; //1530 1610
  int ohdu_numer = 4;
  //number of bins to take into account for chi2
  int bines = 100;

    //const int nBads=6;
  //int Bads[nBads]={8, 110, 158, 194, 230, 322};
  
  //const int nBads=6;
  //int Bads[nBads]={8, 110, 158, 194, 230, 322};

  //const int nBads=25;
  //int Bads[nBads]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 109, 110, 111, 157, 158, 159, 193, 194, 195, 229, 230, 231, 321, 322, 323};

// Calibration parameters //////////////////////////////////////////////

  double alfa = 0.0043;
  double beta = 0.0000000003;

////////////////////////////////////////////////////////////////////////

  int Entries_exp = 1;
  //int xmin = 0; // range for histograms
  //int xmax =0; // range for histograms
  //int xBary_min=0;int xBary_max=100;
  //int yBary_min=0;int yBary_max=100;

  int nbins = 6;

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
  
  Float_t nCal;
  Double_t aduPix[maxClusterSize];
  Float_t eMeanFix[maxClusterSize];

using namespace std;

////////////////////////////////////////////////////////////////////////
//  Older compilers
////////////////////////////////////////////////////////////////////////
string itos(const int i){
	ostringstream ossAux;
	ossAux << i;
	return ossAux.str();
}

void mean_evolution(){
//int main(int argc, char* argv[]){

// Get input files//////////////////////////////////////////////////////

			std::ofstream myfile;
			myfile.open ("fano.txt", std::ios_base::app);
									
            //TFile * f_exp = TFile::Open("output_2_142K.root");
            //TFile * f_exp = TFile::Open("meanFix__cal__hits_skp_55Fe_2000Samp_2_2.root");
			TFile * f_exp = TFile::Open("merge.root");
            if (!f_exp->IsOpen()) {std::cerr << "ERROR: cannot open the root file with experimental data" << std::endl;}
            TTree * texp = (TTree*) f_exp->Get("hitSumm");
            
            //for the total of clusters
            //TH1D * h_exp_e_total  =  new TH1D("h_exp_e_total","", bines, emin, emax);
            //h_exp_e_total -> Sumw2();

            //for each n
            TF1 * fit[nbins+1];
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
		
			double cal = 0;
			double ecal = 0;
			
////////////////////////////////////////////////////////////////////////

TCanvas * c1   = new TCanvas("Electrons","Electrons",800,500);
c1->Divide(3,2);


				for (int ene=1; ene<nbins+1; ene++){
					
					int i = 0;
					c1->cd(ene);
					
						for(int i_event=0;i_event<Entries_exp; i_event++){
						texp->GetEntry(i_event);
										
								//for (int p = 0; p < nSavedPix; ++p){		
										//cal=232.518-0.00744*ePix[p]+0.00000321*pow(ePix[p],2);
										//ecal+=ePix[p]/cal*228.873;
										//ecal+=alfa*ePix[p]+beta*ePix[p]*ePix[p];											
      							//}
      							
										if (e>emin && e<emax){  // number of electrons

											bool noBadPixInCluster = true;
											//cout<<"nSavedPix: "<<nSavedPix<<endl;
											for (int p = 0; p < nSavedPix; ++p){
												for (int j = 0; j < nBads; ++j){
													if(xPix[p]==Bads[j]){
														noBadPixInCluster = false;
														//cout<<"p: "<<p<<"   xPix= "<<xPix[p]<<endl;
														break;
													}
													//cout<<p<<endl;
												}
												if (noBadPixInCluster ==false) break;
												//cout<<"termino loop sobre j  p vale: "<<p<<endl;
											}


											if (noBadPixInCluster){
												
													if (n==ene){
														//h_exp_e[ene]->Fill(ecal);
														//h_exp_e_total->Fill(ecal);
														//events_exp[i]=ecal;

														h_exp_e[ene]->Fill(nCal);
														events_exp[i]=nCal;
														i++;
													}
												
											}	
										}
								
						ecal = 0;
						}
						cout<<"numero de eventos guardados en el histograma: "<<i<<endl;

////////////////////////////////////////////////////////////////////////

						TF1 fDG("fDG","gaus +gaus(3)",1500,1850);
						fDG.SetParameter(0,100);
						//fDG.SetParameter(1,1570);
						fDG.SetParameter(1,1620);
						fDG.SetParameter(2,15);
						fDG.SetParameter(3,10);
						//fDG.SetParameter(4,1740);
						fDG.SetParameter(4,1780);
						fDG.SetParameter(5,16);
						
						h_exp_e[ene]->Draw();
						//fDG.Draw("same");
						
						//c1->WaitPrimitive();
						
						h_exp_e[ene]->Fit("fDG","");
						fit[ene] = (TF1*)h_exp_e[ene]->GetFunction("fDG");
									
						mean_exp_fit[ene]= fit[ene]->GetParameter(1);
						mean_exp_fit_error[ene]= fit[ene]->GetParError(1);

						sigma_exp_fit[ene]= fit[ene]->GetParameter(2);
						sigma_exp_fit_error[ene]= fit[ene]->GetParError(2);
											
						fano_exp_fit[ene]= pow(sigma_exp_fit[ene],2)/mean_exp_fit[ene];
						fano_exp_fit_error[ene]= pow(pow(2*sigma_exp_fit[ene]*sigma_exp_fit_error[ene]/mean_exp_fit[ene],2)+pow(mean_exp_fit_error[ene]*sigma_exp_fit[ene]*sigma_exp_fit[ene]/(mean_exp_fit[ene]*mean_exp_fit[ene]),2),0.5);
						
						events_exp_fit[ene]=h_exp_e[ene]->GetEntries();
					
							
							//h_exp_e[ene]->Draw("HIST E1");
							//fit[ene]->Draw("APL same");
									
				}
							
						// save to file
						myfile << "N" << "	" << "Mean" << "	" << "Mean Error" << "	" <<"Sigma" << "	" <<"Fano factor" << "	" <<"Fano error" << "	" <<" increasing n" <<endl;
						for(int i = 1; i < nbins+1; i ++) {
							myfile << i << "	" << mean_exp_fit[i] << "	" << mean_exp_fit_error[i]  << "	" << sigma_exp_fit[i] << "	" << fano_exp_fit[i] << "	" << fano_exp_fit_error[i]  << "	" <<  events_exp_fit[i] << endl;
						}
						
////////////////////////////////////////////////////////////////////////

	/*
						h_exp_e_total->Fit("gaus");
						TF1 *fit2 = (TF1*)h_exp_e_total->GetFunction("gaus");

						h_exp_e_total->Draw("HIST E2");
						fit2->Draw("HIST same");
						
						mean_exp_fit[0]= fit2->GetParameter(1);
						mean_exp_fit_error[0]= fit2->GetParError(1);

						sigma_exp_fit[0]= fit2->GetParameter(2);
						sigma_exp_fit_error[0]= fit2->GetParError(2);
							
						fano_exp_fit[0]= pow(sigma_exp_fit[0],2)/mean_exp_fit[0];
						fano_exp_fit_error[0]= pow(pow(2*sigma_exp_fit[0]*sigma_exp_fit_error[0]/mean_exp_fit[0],2)+pow(mean_exp_fit_error[0]*sigma_exp_fit[0]*sigma_exp_fit[0]/(mean_exp_fit[0]*mean_exp_fit[0]),2),0.5);
						
						events_exp_fit[0]=h_exp_e_total->GetEntries();

						myfile  << "Para todos los n" << endl;
						myfile  << mean_exp_fit[0]  << "	" << sigma_exp_fit[0] << "	" << fano_exp_fit[0] << "	" << fano_exp_fit_error[0]  << "	" <<  events_exp_fit[0] << endl;
						
									
////////////////////////////////////////////////////////////////////////
	    
	cout<<"Fano: "<<fano_exp_fit[0]<<endl;
	cout<<"Error: "<<fano_exp_fit_error[0]<<endl;
	cout<<"Eventos: "<<events_exp_fit[0]<<endl;
	
	*/
	myfile.close();

////////////////////////////////////////////////////////////////////////
   TCanvas *c2 = new TCanvas("c2","Evolucion del valor medio",200,10,700,500);

   c2->SetFillColor(42);
   c2->SetGrid();

   cout<<endl<<endl;
   printf(" n    Media       error    Fano    error  porcentual ");
   cout<<endl;
   
   const Int_t n = 6;
   Double_t x[n], y[n],ex[n], ey[n], sigma[n], esigma[n], fano[n], efano[n];
   for (Int_t i=0;i<n;i++) {
	   x[i] = i+1;
	   ex[i] = 0;
	   y[i] = mean_exp_fit[i+1]; 
	   ey[i] = mean_exp_fit_error[i+1]; 

	   sigma[i]=sigma_exp_fit[i+1];
	   esigma[i]=sigma_exp_fit_error[i+1];

	   fano[i]=pow(sigma[i],2)/y[i];
	   efano[i]=fano[i]*sqrt(pow((2*esigma[i]/sigma[i]),2)+pow((ey[i]/y[i]),2));

	   printf(" %i %f %f %f %f %f\n",i+1,y[i], ey[i], fano[i], efano[i], 100*efano[i]/fano[i]);
   }
   TGraph *gr = new TGraphErrors(n,x,y,ex,ey);
   gr->SetLineColor(2);
   gr->SetLineWidth(4);
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->SetTitle("Evolucion del valor medio");
   gr->GetXaxis()->SetTitle("Tamanio del cluster");
   gr->GetYaxis()->SetTitle("Media (electrones)");
   gr->Draw("ACP");

   // TCanvas::Update() draws the frame, after which one can change it
   c2->Update();
   c2->GetFrame()->SetFillColor(21);
   c2->GetFrame()->SetBorderSize(12);
   c2->Modified();
	
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

// Los siguientes son los branches que agrego Andre y tienen la carga calibrada por Cluster
  tree->SetBranchAddress ("nCal",&nCal);
  tree->SetBranchAddress ("aduPix",aduPix);
  tree->SetBranchAddress ("eMeanFix",eMeanFix);
}