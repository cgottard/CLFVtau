#include "mymap.h"

using namespace std;

string ReadFile::getFileName(string filelabel) // a function to get filename from label
{
	filename = files.at(filelabel);
	return filename;
}

/*

void PlotSB(string a){

  for (auto &histo : hists) //for each variable

  {

    string variable = histo.first;

    TCanvas *c = new TCanvas("plot", "plot",100,0,200);

    hedeh.getHist(variable)->Draw();

    string pdfname2 = a+"_"+variable+".pdf";

    c->SaveAs(pdfname2.c_str());

  

    delete pdfname2;

    delete c;

  }

}

*/

