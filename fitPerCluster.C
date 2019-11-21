#include <stdio.h>
#include <iostream>
#include <TFile.h>

void fitPerCluster(std::string input){
        std::cout << input << std::endl;
        gROOT->Reset();
        TFile *f    = new TFile("input.root","RECREATE");
        TTree *tree = new TTree("T","Gaussian diff size of cluster");
        Int_t n;
        Double_t val;
        tree->Branch("n",&n,"n/I");
        tree->Branch("val",&val,"val/D");
        TRandom3 rndgen;

        for (Int_t k = 0; k < 5; ++k){
                for (Int_t i = 0; i<1000; i++){
                        n = k;
                        val = rndgen.Gaus(50.0,10);
                        std::cout << n << " " << val << std::endl;
                        tree->Fill();
                }
        }

        TH1D *h = new TH1D ("randomGauss","randomGauss",/*nbins=*/100,/*xmin=*/0,/*xmax=*/100);

        tree->Write();
        TFile *myfile = new TFile("input.root","READ");
        TTree *t = (TTree*) myfile->Get("T");
        t->SetBranchAddress("n",&n);
        t->SetBranchAddress("val",&val);
        Int_t entries = t->GetEntries();
        for (int j = 0 ; j < 5 ; ++j){
                for (Int_t i =0; i < entries ; ++i){
                        t->GetEntry(i);
                        if (n == j) h->Fill(val);
                }
                h->Fit("gaus","Q");
                TF1 *g = (TF1*)h->GetListOfFunctions()->FindObject("gaus");
                std::cout << "mean: " << g->GetParameter(1) << std::endl;
                gPad->WaitPrimitive();
                h->Reset("ICESM");
        }
}
