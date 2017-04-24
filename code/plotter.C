#include "Lcore.h"

int Plot(const int opt);

int main(const int argc, const char * argv[]){

  if (argc < 2) return 0;
  MODEL::SetMODEL();
  GOLD::SetGOLD();

  int opt = atoi(argv[1]);

  Plot(opt);


  return 0;
}



int Plot(const int opt){

  if (opt == 1){//wave function
    double x[100], y[100];
    for (int i = 0; i < 100; i++){
      x[i] = i * 0.03;
      y[i] = MODEL::Ur(x[i] / Phys::hbar);
    }
    TGraph * g0 = new TGraph(100, x, y);
    g0->SetTitle("");
    g0->SetLineWidth(2);
    g0->SetLineColor(4);
    g0->GetXaxis()->SetTitle("r (fm)");
    g0->GetXaxis()->CenterTitle(true);
    g0->GetXaxis()->SetTitleOffset(1.15);
    g0->GetXaxis()->SetTitleSize(0.06);
    g0->GetXaxis()->SetLabelSize(0.06);
    g0->GetXaxis()->SetRangeUser(0.0, 3.0);
    g0->GetYaxis()->SetTitle("u(r) (GeV^{1/2})");
    g0->GetYaxis()->CenterTitle(true);
    g0->GetYaxis()->SetTitleOffset(1.15);
    g0->GetYaxis()->SetTitleSize(0.06);
    g0->GetYaxis()->SetLabelSize(0.06);
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    g0->DrawClone("AC");
    c0->Print("gallary/ur.pdf");
  }

  if (opt == 2){//potential
    double x[100], y[100];
    for (int i = 0; i < 100; i++){
      x[i] = i * 0.03;
      y[i] = MODEL::Veff(x[i] / Phys::hbar);
    }
    TGraph * g0 = new TGraph(100, x, y);
    g0->SetTitle("");
    g0->SetLineWidth(2);
    g0->SetLineColor(4);
    g0->GetXaxis()->SetTitle("r (fm)");
    g0->GetXaxis()->CenterTitle(true);
    g0->GetXaxis()->SetTitleOffset(1.15);
    g0->GetXaxis()->SetTitleSize(0.06);
    g0->GetXaxis()->SetLabelSize(0.06);
    g0->GetXaxis()->SetRangeUser(0.0, 3.0);
    g0->GetYaxis()->SetTitle("V_{eff}(r) (GeV)");
    g0->GetYaxis()->CenterTitle(true);
    g0->GetYaxis()->SetTitleOffset(1.15);
    g0->GetYaxis()->SetTitleSize(0.06);
    g0->GetYaxis()->SetLabelSize(0.06);
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    g0->DrawClone("AC");
    c0->Print("gallary/Vr.pdf");
  }

  if (opt == 3){//photoproduction cross section
    ifstream infile("photoproductionsigma.txt");
    double x[11], y[11];
    TString tmp;
    for (int i = 0; i < 11; i++){
      infile >> x[i] >> y[i] >> tmp;
    }
    infile.close();
    TGraph * g0 = new TGraph(11, x, y);
    g0->SetTitle("");
    g0->SetLineWidth(2);
    g0->SetLineColor(4);
    g0->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
    g0->GetXaxis()->CenterTitle(true);
    g0->GetXaxis()->SetTitleOffset(1.15);
    g0->GetXaxis()->SetTitleSize(0.06);
    g0->GetXaxis()->SetLabelSize(0.06);
    g0->GetXaxis()->SetRangeUser(6.0, 11.0);
    g0->GetYaxis()->SetTitle("#sigma (nb)");
    g0->GetYaxis()->CenterTitle(true);
    g0->GetYaxis()->SetTitleOffset(1.15);
    g0->GetYaxis()->SetTitleSize(0.06);
    g0->GetYaxis()->SetLabelSize(0.06);\
    g0->GetYaxis()->SetRangeUser(1.0e-9, 1.0e-6);
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLogy();
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    g0->DrawClone("AC");
    c0->Print("gallary/photosigma.pdf");
  }

  if (opt == 4){//formation probability
    ifstream infile("formprobability.txt");
    double x[46], y[46];
    TString tmp;
    for (int i = 0; i < 46; i++){
      infile >> x[i] >> y[i];
    }
    infile.close();
    TGraph * g0 = new TGraph(46, x, y);
    g0->SetTitle("");
    g0->SetLineWidth(2);
    g0->SetLineColor(4);
    g0->GetXaxis()->SetTitle("P_{J/#psi} (GeV)");
    g0->GetXaxis()->CenterTitle(true);
    g0->GetXaxis()->SetTitleOffset(1.15);
    g0->GetXaxis()->SetTitleSize(0.06);
    g0->GetXaxis()->SetLabelSize(0.06);
    g0->GetXaxis()->SetRangeUser(0.0, 9.0);
    g0->GetYaxis()->SetTitle("Probability");
    g0->GetYaxis()->CenterTitle(true);
    g0->GetYaxis()->SetTitleOffset(1.15);
    g0->GetYaxis()->SetTitleSize(0.06);
    g0->GetYaxis()->SetLabelSize(0.06);
    g0->GetYaxis()->SetNdivisions(6, 0, 0);
    g0->GetYaxis()->SetRangeUser(1.0e-11, 0.9);
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLogy();
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    g0->DrawClone("AC");
    c0->Print("gallary/formprobability.pdf");
  }




  return 0;
}
