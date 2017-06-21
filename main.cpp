
 
#include <iostream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TEntryList.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH1.h"
#include "nominal.h"
#include "mymap.h"
#include <map>
#include <utility>
#include <vector>
#include <iomanip>


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
	/* Only entered files and parameters will be evaluated.
	   Parameter and file labels have to match with mymap.h file.
	   Changes can be made in mymap.h as desired
	*/
	
	vector<string> file_label = {"Signal","2l2q","2l2v","3lv"}; 				//File Labels

	vector<string> param_label = {"hist","hist6","hist7","hist8","hist9","hist10","hist11","hist12"}; // Parameter labels
		
	//int entries = 0;
	int nFiles  = file_label.size();	
	string a, c, pdfname;
	
	for (int x=0; x<nFiles; ++x) 
	{
	
		ReadFile fe;
		a = fe.getFileName(file_label.at(x));	//contains the location of the file
		c = file_label.at(x);			//contains the filelable
		TString c_copy = file_label.at(x); 				
		
		TFile *ntuple = new TFile(a.c_str(), "READ"); 	//open the file in read mode
		TTree *mytree = (TTree*)ntuple->Get("nominal"); 
		nominal N;	
		N.Init(mytree);

		TFile *new_ntuple = new TFile(c_copy+"_new.root","RECREATE");// write to a new root file
		int entries=0;
		entries = mytree->GetEntries();	//get the total number of events in the tree	
		cout << "DEBUG: Number of entries " << entries << endl;	

		int SMcounter=0;
		int LFVcounter=0;
		int taus_in_LFV=0;
		int samesign=0;
		int oppsign=0;
		
		PlotHist h;
		
		//event selection 
		
		string myselection ="nEl>=1";// "(@tau_pt.size()==1)" //&& (nEl+nMu==2)";//nEl==2 || nMu==2"; //(mu_charge[0]==mu_charge[1] || el_charge[0]==el_charge[1])";
		N.fChain->Draw(">>eventlist", myselection.c_str(), "entrylist");
		TEntryList *elist = (TEntryList *) gDirectory->Get("eventlist");
		N.fChain->SetEntryList(elist);
		
		entries = elist->GetN();
	
		cout << "Number of entries after the selection is: " << entries << endl;
		
		cout << "DEBUG: Starting the loop "<< endl;
		//start of event loop	
		for(unsigned int i=0; i<entries; ++i)
		{
			
			float progress = 100.0 * ((float) i) / ((float) entries);      
			if (!(i % 10)) cout << setprecision(3) << " [ " << progress << " % ] \r";
			cout << "Working on File: [" << x+1 <<"/"<< nFiles << "] | " << file_label.at(x) << " | " << setprecision(3) << "[ " << progress << " % ] \r";

	
			N.fChain->GetEntry(elist->GetEntry(i));
    			//N.GetEntry(i);


			if(!(N.charge_1lep==N.charge_2lep || N.charge_2lep==N.charge_3lep || N.charge_3lep==N.charge_1lep)) {++oppsign;}
			if(N.charge_1lep==N.charge_2lep || N.charge_2lep==N.charge_3lep || N.charge_3lep==N.charge_1lep) ++samesign;

		
			//get i-th entry
			vector<LEPTON> leptons;
			vector<JET> jets; 
			vector<TOP> top;
		
			//cout << "DEBUG: Got entry " << i << endl;
		
			if(!(N.hasTrigMatchedEl || N.hasTrigMatchedMu)) continue;
		
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
			if(jets.size()>0) sort(jets.begin(), jets.end(), jetsort);
	
			//taus_in_LFV=0;
	
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
						//if(leptons[k].flavour==15 || leptons[l].flavour==15) ++taus_in_LFV;
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
				//hist->Fill(top[0].v3.M());
				h.getHist("hist")->Fill(top[0].v3.M());
				if(top[0].v3.Pt()>0.1) LFVcounter++;
			}
	
			h.getHist("hist")->GetYaxis()->SetTitleOffset(1.5);

			//SM top - W masses
			double mu, t1, t2, t3, nu_Pz_pos, nu_Pz_neg, nu_E_pos, nu_E_neg;
			double Wmass = 80.4; 
			float met = N.met_met/1000.;		
			vector<TLorentzVector> nu_pos, nu_neg, W_pos, W_neg;	
			for(int j=0; j<leptons.size(); ++j)
			{
			
				if(j == B || j == C) continue; //discarding the leptons already used in the reconstruction of the LFV top 
			
				float Dphi = fabs(M_PI-fabs(N.met_phi - leptons[j].v1.Phi()));
				mu = (Wmass/2.) + leptons[j].v1.Pt() * met * cos(Dphi);
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
			
			/*
			if(W_pos.size()>0)
				hist2->Fill(W_pos[0].M());

			if(W_neg.size()>0)
				hist3->Fill(W_neg[0].M());
			*/
		
			//SM top mass
			vector<TOP> SMtop;
		
			//pos
			for(int j=0; j<jets.size(); ++j)
			{
				if(j == A) continue;
				for(int k=0; k<W_pos.size(); ++k)
				{
					SMtop.push_back(TOP());
					SMtop.back().v3 = jets[j].v2 + W_pos[k];
				}	
			}
			//neg
			for(int j=0; j<jets.size(); ++j)
			{
				if(j == A) continue;
				for(int k=0; k<W_neg.size(); ++k)
				{
					SMtop.push_back(TOP());
					SMtop.back().v3 = jets[j].v2 + W_neg[k];
				}	
			}
		
			// combining the neg and pos tops to get one with the highest pT	
			
			if(SMtop.size()>0) {
				sort(SMtop.begin(), SMtop.end(), topsort); 
				//hist6->Fill(SMtop[0].v3.M());
				h.getHist("hist6")->Fill(SMtop[0].v3.M());
				if (SMtop[0].v3.Pt()>0.1) SMcounter++;
			}	
			
			/*
			if(SMtop_pos.size()>0)
				hist4->Fill(SMtop_pos[0].M()); 
		
			if(SMtop_neg.size()>0)
				hist5->Fill(SMtop_neg[0].M()); 
			*/
		
			for(int j=0; j<N.tau_pt->size(); ++j)
				h.getHist("hist7")->Fill(N.tau_pt->at(j)/1000);
				//hist7->Fill(N.tau_pt->at(j)/1000);

			for(int j=0; j<N.tau_phi->size(); ++j)
				h.getHist("hist8")->Fill(N.tau_phi->at(j));
				//hist8->Fill(N.tau_phi->at(j));
		
			for(int j=0; j<N.tau_eta->size(); ++j)
				h.getHist("hist9")->Fill(N.tau_eta->at(j));
				//hist9->Fill(N.tau_eta->at(j));
			
			h.getHist("hist10")->Fill(N.tau_pt->size());
			h.getHist("hist11")->Fill(N.nEl);
			h.getHist("hist12")->Fill(N.nMu);
			//hist10->Fill(N.tau_pt->size());
			//hist11->Fill(N.nEl);
			//hist12->Fill(N.nMu);
		
			for(int j=0; j<N.el_charge->size(); ++j)
				h.getHist("hist13")->Fill(N.el_charge->at(j));
				//hist13->Fill(N.el_charge->at(j));
	
			//for(int j=0; j<N.mu_charge->size(); ++j)
				//hist14->Fill(N.charge_1lep);

		} // end of event loop
		
		cout << "Number of reconstructed LFV tops (with pT>0.1) is: " << LFVcounter <<endl;
		cout << "Number of reconstructed SM tops (with pT>0.1) is: " << SMcounter <<endl;
		//cout << "Number of taus in LFV is: " << taus_in_LFV <<endl;
		cout << "Number of entries with same sign leptons is: " << samesign <<endl;
		cout << "Number of entries with opp sign leptons is: " << oppsign <<endl;
	
		for(int t=0; t<param_label.size(); ++t) {
			
			TCanvas *can = new TCanvas("c","c", 800,600);
			h.getHist(param_label.at(t))->Draw();
		  
			new_ntuple->cd();							 // tell program you are at new_ntuple
			h.getHist(param_label.at(t))->Write(); 		 //write to folder

			pdfname =c + "_" + param_label.at(t) + "_new" + ".pdf";
			can->SaveAs(pdfname.c_str());

			delete can;
		}
	}	       	
    return 0;
}
