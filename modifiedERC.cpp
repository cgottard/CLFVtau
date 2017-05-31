#include <iostream>
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "nominal.h"
#include <iomanip>
#include "TLorentzVector.h"

using namespace std;

/*
C++ code to be compiled with g++ which uses ROOT libraries..
It plots the invariant mass of opposite sign lepton pairs
To compile
//g++ -std=c++11 -c ../macro/nominal.C -o nominal.o `root-config --cflags --glibs`
// g++ -std=c++11 EventReaderCompiled.cpp nominal.o -o ev `root-config --cflags --glibs`
*/

struct LEPTON
{
	float charge;
	float flavour;
	TLorentzVector v1;
};


struct JET
{
	float btag;
	TLorentzVector v2;
};

struct TOP
{
	int index [3];
	TLorentzVector v3;
};

//descending order sorting (for jet btag)
bool jetsort(JET& a1, JET& a2) {return (a1.btag > a2.btag);}

//descending order sorting (for top pt)
bool topsort(TOP& b1, TOP& b2) {return (b1.v3.Pt() > b2.v3.Pt());}

int main()
{
	TFile *ntuple = new TFile("~/Signal_v9ta.root", "READ"); //open the file in read mode
	TTree *mytree = (TTree*)ntuple->Get("nominal");
	nominal N;	
  	N.Init((TTree *)ntuple->Get("nominal"));
	//get the total number of events in the tree
	long entries = mytree->GetEntries();
	cout << "DEBUG: Number of entries " << entries << endl;	
    TH1F *hist = new TH1F("hist","top mass",80,0,400);
    TH1F *hist2 = new TH1F("hist2","W_pos mass",80,0,400);
	TH1F *hist3 = new TH1F("hist3","W_neg mass",80,0,400);
	TH1F *hist4 = new TH1F("hist4","SMtop_pos mass",80,0,400);
	TH1F *hist5 = new TH1F("hist5","SMtop_neg mass",80,0,400);

	cout << "DEBUG: Starting the loop "<< endl;

	//start of event loop	
	for(int i=0; i<entries; ++i)
	{
		float progress = 100.0 * ((float) i) / ((float) entries);       
    	if (!(i % 10)) cout << setprecision(3) << " [ " << progress << " % ] \r";
		
		//get i-th entry
		vector<LEPTON> leptons;
		vector <JET> jets; 
		vector<TOP> top;
		N.GetEntry(i);
		//cout << "DEBUG: Got entry " << i << endl;

		//selection
		for(int j=0; j<N.nJets; ++j)
		{
			jets.push_back(JET());	
			jets.back().btag = N.jet_mv2c10->at(j);
			jets.back().v2.SetPtEtaPhiE(N.jet_pt->at(j)/1000., N.jet_eta->at(j), N.jet_phi->at(j), N.jet_e->at(j)/1000.);
		}	
		//electrons
		for(int k=0; k<N.nEl; ++k)
		{			
			leptons.push_back(LEPTON());	
			leptons.back().charge  = N.el_charge->at(k);
			leptons.back().flavour = 11;
			leptons.back().v1.SetPtEtaPhiM(N.el_pt->at(k)/1000., N.el_eta->at(k), N.el_phi->at(k), 0.511e-3);		
		}	
		//muons	
		for(int l=0; l<N.nMu; ++l)
		{		
			leptons.push_back(LEPTON());		
			leptons.back().charge  = N.mu_charge->at(l);
			leptons.back().flavour = 13; 
			leptons.back().v1.SetPtEtaPhiM(N.mu_pt->at(l)/1000., N.mu_eta->at(l), N.mu_phi->at(l), 105.6e-3);
		}	
		//tauons
		for(int m=0; m<N.tau_pt->size(); ++m)
                {
                        leptons.push_back(LEPTON());
                        leptons.back().charge  = N.tau_charge->at(m);
                        leptons.back().flavour = 15;
                        leptons.back().v1.SetPtEtaPhiM(N.tau_pt->at(m)/1000., N.tau_eta->at(m), N.tau_phi->at(m), 1.776);
                }

		//sorting jets for max btag likelihood
		if(jets.size()>0) 
			sort(jets.begin(), jets.end(), jetsort);
	
		for(int j=1; j<jets.size(); ++j) //max btagged discarded
		{
			for(int k=0; k<leptons.size(); ++k)
			{
				for(int l=k+1; l<leptons.size(); ++l)
				{
					if(leptons[k].charge == leptons[l].charge) continue; //same charged leptons discarded
					if(leptons[k].flavour == leptons[l].flavour) continue; //same flavour leptons discarded
					top.push_back(TOP());		
					top.back().v3 = jets[j].v2 + leptons[k].v1 + leptons[l].v1;
					top.back().index[0] = j; // ljet               
					top.back().index[1] = k; // lep 1                
					top.back().index[2] = l; // lep 2					
				}
			}
		}
		if(top.size()==0) continue;
		//LFV top (highest pt)
		int A{99}, B{99}, C{99};
		if(top.size()>0)
		{
			sort(top.begin(), top.end(), topsort); 
			A = top[0].index[0];
			B = top[0].index[1];
			C = top[0].index[2];
			hist->Fill(top[0].v3.M());
		}
	

		//SM top - W masses
		double mu, t1, t2, t3, nu_Pz_pos, nu_Pz_neg, nu_E_pos, nu_E_neg;
		double Wmass = 80.4; 
		float met = N.met_met/1000.;		
		vector<TLorentzVector> nu_pos, nu_neg, W_pos, W_neg;	
		for(int j=0; j<leptons.size(); ++j)
		{
			
			if(j == B || j == C) continue; //discarding the leptons already used in the reconstruction of the LFV top 
			
			float Dphi = fabs(M_PI-fabs(N.met_phi - leptons[j].v1.Phi()));
			mu = (Wmass/2.) + leptons[j].v1.Pt()  * met * cos(Dphi);
			t1 = mu * leptons[j].v1.Pz() / pow(leptons[j].v1.Pt(),2.);
			t2 = pow(mu,2.) * pow(leptons[j].v1.Pz(),2.) / pow(leptons[j].v1.Pt(),4.);
			t3 = (pow(leptons[j].v1.E(),2.) * pow(met,2.) - pow(mu,2.)) / pow(leptons[j].v1.Pt(),2.);

			nu_Pz_pos = t1 + sqrt(t2 - t3);
			nu_Pz_neg = t1 - sqrt(t2 - t3);
				
			nu_E_pos = sqrt(pow(met,2.) + pow(nu_Pz_pos,2.));
			nu_E_neg = sqrt(pow(met,2.) + pow(nu_Pz_neg,2.));
			
			nu_pos.push_back(TLorentzVector());
			nu_pos.back().SetPxPyPzE(met* cos(N.met_phi), met * sin(N.met_phi), nu_Pz_pos, nu_E_pos);
			nu_neg.push_back(TLorentzVector());
			nu_neg.back().SetPxPyPzE(met * cos(N.met_phi), met * sin(N.met_phi), nu_Pz_neg, nu_E_neg);
		
			W_pos.push_back(TLorentzVector());
			W_pos.back() = leptons[j].v1 + nu_pos.back();
			W_neg.push_back(TLorentzVector());
			W_neg.back() = leptons[j].v1 + nu_neg.back();
		}
		
		if(W_pos.size()>0)
			hist2->Fill(W_pos[0].M());

		if(W_neg.size()>0)
			hist3->Fill(W_neg[0].M());

		//SM top mass
		vector<TLorentzVector> SMtop_pos, SMtop_neg;
		//pos
		for(int j=0; j<jets.size(); ++j)
		{
			if(j == A) continue;
			for(int k=0; k<W_pos.size(); ++k)
			{
				SMtop_pos.push_back(TLorentzVector());
				SMtop_pos.back() = jets[j].v2 + W_pos[k];
			}
		}
		//neg
		for(int j=0; j<jets.size(); ++j)
		{
			if(j == A) continue;
			for(int k=0; k<W_neg.size(); ++k)
			{
				SMtop_neg.push_back(TLorentzVector());
				SMtop_neg.back() = jets[j].v2 + W_neg[k];
			}
		}

		
		if(SMtop_pos.size()>0)
			hist4->Fill(SMtop_pos[0].M());
		
		if(SMtop_neg.size()>0)
			hist5->Fill(SMtop_neg[0].M());

  	} // end of event loop
	
	TCanvas c("c","c",800,600);
	hist->Draw();
	c.SaveAs("LFV_top_mass_modified.pdf");

	TCanvas c2("c2","c2",800,600);
	hist2->Draw();
	c2.SaveAs("W_pos_mass.pdf");	
	
	TCanvas c3("c3","c3",800,600);
	hist3->Draw();
	c3.SaveAs("W_neg_mass.pdf");

	TCanvas c4("c4","c4",800,600);
	hist4->Draw();
	c4.SaveAs("SMtop_pos_mass.pdf");

	TCanvas c5("c5","c5",800,600);
	hist5->Draw();
	c5.SaveAs("SMtop_neg_mass.pdf");
return 0;
}