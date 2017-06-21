#include "TH1F.h"
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <cmath>
#include "TH1F.h"
#include "TH1.h"
#include <vector>
#include "TFile.h"
#include <iostream>
#include <string>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "nominal.h" // Include nominal
#include <map>
#include <utility>

//Define pairs for root files associated with labels
class ReadFile{

private:
std::string filelabel= "";
std::string filename= "";
std::map<std::string, std::string> files = {
	std::make_pair<std::string, std::string>("Signal", "/lustre/user/aruina/rootfiles/Signal_v9ta.root"),
	std::make_pair<std::string, std::string>("2l2q",   "/lustre/user/aruina/rootfiles/2l2q_v9ta.root"),
	std::make_pair<std::string, std::string>("2l2v",   "/lustre/user/aruina/rootfiles/2l2v_v9ta.root"),
	std::make_pair<std::string, std::string>("3lv",    "/lustre/user/aruina/rootfiles/3lv_v9ta.root"),
};

public:
std::string getFileName(std::string filelabel);

};

//Define pairs for histograms associated with labels
class PlotHist{

private:
std::string histlabel= "";       
std::string histlabel2= "";          
std::map<std::string, TH1F> hists = {
	
	std::make_pair<std::string, TH1F>("hist",  {"hist",  "LFVtop mass;Energy [GeV];Events",80,0,400}),
	std::make_pair<std::string, TH1F>("hist6", {"hist6", "SMtop mass;Energy [GeV];Events",80,0,400}),
    std::make_pair<std::string, TH1F>("hist7", {"hist7", "tau_pt;[GeV]",80,-50,350}),
    std::make_pair<std::string, TH1F>("hist8", {"hist8", "tau_phi",80,-4,4}),
    std::make_pair<std::string, TH1F>("hist9", {"hist9", "tau_eta",80,-3,3}),
    std::make_pair<std::string, TH1F>("hist10",{"hist10","nTa",80,-1,5}),
	std::make_pair<std::string, TH1F>("hist11",{"hist11","nEl",80,-1,5}),
	std::make_pair<std::string, TH1F>("hist12",{"hist12","nMu",80,-1,5}),
};

public:
	TH1F* getHist(std::string histlabel){return &hists.at(histlabel);};
};

