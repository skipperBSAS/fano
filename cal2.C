/*
Enero 2019 - Dario Rodrigues
Este script permite hacer una calibracion no lineal ajustando de forma iterativa los picos correspondientes a j electrones usando gaussianas.
Permite poner cortes de calidad en chi2, la distancia al pico anterior y la estadistica del pico.
Al finalizar muestra el histograma de las cargas y debajo los puntos que pudieron ajusterse con gaussianas.
Como resultado se puede obtener la calibracion, aunque en el ajuste final resta quitar los puntos vacios.
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

#include <chrono>
#include <thread>

using namespace std;

void Enable_and_Set_Branches(TTree* & tree);

////////////////////////////////////////////////////////////////////////

int Entries_exp = 1;

int runID; int ohdu; 
int x;int y;
double pix;
  
using namespace std;

////////////////////////////////////////////////////////////////////////
//  Older compilers
////////////////////////////////////////////////////////////////////////
string itos(const int i){
	ostringstream ossAux;
	ossAux << i;
	return ossAux.str();
}

void cal2(){
//int main(int argc, char* argv[]){

// Get input files//////////////////////////////////////////////////////
 			int npeaks = 1700;
 			int constantmin = 20;
 			double chi2max = 100;
 			
 			
			std::ofstream myfile;
			myfile.open ("calibration_142K.txt", std::ios_base::app);
									
            TFile * f_exp = TFile::Open("output_142K_Andre.root");
            if (!f_exp->IsOpen()) {std::cerr << "ERROR: cannot open the root file with experimental data" << std::endl;}
            
            TTree * texp = (TTree*) f_exp->Get("skPixTree");
            	
 			int emin=0; 
			int emax=230*npeaks;
			
			int bines;
			bines = (emax - emin)/20;
			
            // Get information from trees///////////////////////////////////////////
            int Entries_exp = texp -> GetEntries();
            cout<<"Entries in experimental data file: "<<Entries_exp<<endl;
            Enable_and_Set_Branches(texp);
			
////////////////////////////////////////////////////////////////////////

				TH1D * h  =  new TH1D("Electrons","", bines, emin, emax);
				TCanvas * c1   = new TCanvas("Spectrum","Spectrum",800,500);
				c1->Divide(1,2);
				c1->cd(1);
				
				for(int i_event=0;i_event<Entries_exp; i_event++){
					texp->GetEntry(i_event);
					h->Fill(pix);
				}
				
				h->Draw();
					
				int j = 1;
				
				double constant[npeaks];
				constant[1]=90000;
				
				double mean[npeaks];
				double mean_error[npeaks];
				
				mean[0]=0;	
				mean[1]=234;
				
				double delta[npeaks];
				delta[1]= 230;
				
				double sigma[npeaks];
				double sigma_error[npeaks];
				sigma[1] = 40;
				
				double chi2[npeaks];
				double x[npeaks];
				double y[npeaks];
				
				x[0]=0;	y[0]=0;

				myfile << "ADU" << "	" << "e" << "	" <<"Delta" << "	" <<"chi2" <<endl;
				
				for(int i = 1;i<npeaks; i++){
				
						cout<<"Ahora busca un pico mas adelante usando los siguientes parametros iniciales: "<<endl<<endl; 
								    
						cout<<"Constant:              "<<constant[i]<<endl;
						cout<<"Mean:                  "<<mean[i]<<endl;
						cout<<"Sigma:                 "<<sigma[i]<<endl;
						cout<<"Delta:                 "<<delta[i]<<endl<<endl;
											
						TF1 fDG("fDG","gaus",mean[i]-2*sigma[i],mean[i]+2*sigma[i]);
						fDG.SetParameter(0,constant[i]);
						fDG.SetParameter(1,mean[i]);
						fDG.SetParameter(2,sigma[i]);
															
						h->Draw("same");
						//h->SetMaximum(constant[i]*10);	
								
						/*
						if (mean[i]<2000){
							h->SetAxisRange(0,1.25*mean[i],"X");				
						}
						else
						{
							h->SetAxisRange(0.75*mean[i],1.25*mean[i],"X");
					    }
					    */
					    
					    if (i==npeaks-1){
							h->SetAxisRange(0,mean[i],"X");
							h->SetMaximum(140);	
					    }
						
						//c1->WaitPrimitive();
						//gSystem->Sleep (1000);
						
						cout<<"y obtiene los siguientes parametros del ajuste: "<<endl<<endl; 
						
						h->Fit("fDG","R");	// muestra los parametros obtenidos en pantalla											
									
						TF1 * fit = (TF1*)h->GetFunction("fDG");
						                   		
                       // Update parameters	////////////////////////////	
                       		
							constant[i] = fit->GetParameter(0);
							mean[i] = fit->GetParameter(1);
							mean_error[i]= fit->GetParError(1);
							
							sigma[i]= fit->GetParameter(2);		
							sigma_error[i]= fit->GetParError(2);
							
							delta[i+1] = mean[i]-mean[i-1];
							
							chi2[i]= fit->GetChisquare();	
								 					
							cout<<endl<<"i="<<i<<"   chi2[i]: "<<chi2[i]<<"    delta[i+1]: "<<delta[i+1]<<endl;
							
							if (chi2[i]<chi2max && constant[i]>constantmin && delta[i+1]>220 && delta[i+1]<240){
								
										cout<<endl<<"El ajuste fue bueno! y el nuevo pico esta a una distancia razonable del anterior :)"<<endl;				
														
											constant[i+1] = constant[i];		
												
											mean[i+1] = mean[i] + delta[i+1];
											sigma[i+1] = sigma[i];
										
											j = i;
											x[i]=mean[i];
											y[i]=i;
											
											//std::cout << std::setprecision(5) << mean[i]
											myfile << mean[i] << "	" << i << "  " << delta[i+1] << "  " << chi2[i] <<endl;
														
							}
							else
							{						
										cout<<endl<<"El ajuste fue malo o el nuevo pico no esta a una distancia razonable del anterior :/ "<<endl;
										cout<<"Entonces, toma como parametros iniciales del proximo ajuste los ultimos buenos. "<<endl<<endl; 
										
											constant[i+1] = constant[j];				
											delta[i+1] = delta[1];
											mean[i+1] = mean[j] + delta[1]*(i+1-j);
											sigma[i+1]= sigma[1];
											
											x[i]=0;
											y[i]=0;
															
							}						
					} 	
				
   c1->cd(2);
   c1->SetGrid();
 
   TF1 cali("cali","[0]*x+[1]*pow(x,2)",0,npeaks+1);
   
 
   TGraph* gr = new TGraph(npeaks,x,y);
   gr->Draw("A*");
   gr->SetTitle("Calibration");
   gr->GetYaxis()->SetTitle("Electrons");
   gr->GetXaxis()->SetTitle("ADUs");
   
   gr->Fit("cali","R");	
   cali.Draw("same");
   
   //TLine *line = new TLine(0,0,delta[1]*npeaks,npeaks);
   //line->SetLineColor(kRed);
   //line->Draw();	
   
							
}

////////////////////////////////////////////////////////////////////////

void Enable_and_Set_Branches(TTree* & tree){

  tree->SetBranchStatus("*",1); //enable all branches

  tree->SetBranchAddress ("runID",&runID);
  tree->SetBranchAddress ("ohdu",&ohdu);
  tree->SetBranchAddress ("x",&x);
  tree->SetBranchAddress ("y",&y);
  tree->SetBranchAddress ("pix",&pix);
}
