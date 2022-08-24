#include "eff_bins.h"
using namespace std;

static int cno = 0;

TH1D* h1d_gen_npart[pbin][etabin] = {0};

TH2D* h2d_1to1_dp_vs_nhit[pbin][etabin] = {0};
TH2D* h2d_1to1_dphi_vs_nhit[pbin][etabin] = {0};
TH2D* h2d_1to1_dtheta_vs_nhit[pbin][etabin] = {0};
TH2D* h2d_1to1_chi2_vs_nhit[pbin][etabin] = {0};
TH1D* h1d_1to1_dp[pbin][etabin] = {0};
TH1D* h1d_1to1_dphi[pbin][etabin] = {0};
TH1D* h1d_1to1_dtheta[pbin][etabin] = {0};
TH1D* h1d_1to1_chi2[pbin][etabin] = {0};
TH1D* h1d_1to1_nhit[pbin][etabin] = {0};

TH2D* h2d_Nto1_dp_vs_nhit[pbin][etabin] = {0};
TH2D* h2d_Nto1_dphi_vs_nhit[pbin][etabin] = {0};
TH2D* h2d_Nto1_dtheta_vs_nhit[pbin][etabin] = {0};
TH2D* h2d_Nto1_chi2_vs_nhit[pbin][etabin] = {0};
TH1D* h1d_Nto1_dp[pbin][etabin] = {0};
TH1D* h1d_Nto1_dphi[pbin][etabin] = {0};
TH1D* h1d_Nto1_dtheta[pbin][etabin] = {0};
TH1D* h1d_Nto1_chi2[pbin][etabin] = {0};
TH1D* h1d_Nto1_nhit[pbin][etabin] = {0};

TH2D* h2d_0to1_dp_vs_nhit[pbin][etabin] = {0};
TH2D* h2d_0to1_dphi_vs_nhit[pbin][etabin] = {0};
TH2D* h2d_0to1_dtheta_vs_nhit[pbin][etabin] = {0};
TH2D* h2d_0to1_chi2_vs_nhit[pbin][etabin] = {0};
TH1D* h1d_0to1_dp[pbin][etabin] = {0};
TH1D* h1d_0to1_dphi[pbin][etabin] = {0};
TH1D* h1d_0to1_dtheta[pbin][etabin] = {0};
TH1D* h1d_0to1_chi2[pbin][etabin] = {0};
TH1D* h1d_0to1_nhit[pbin][etabin] = {0};

TH2D* h2d_dp_in_p_vs_eta = NULL;
TH2D* h2d_dphi_in_p_vs_eta = NULL;
TH2D* h2d_dtheta_in_p_vs_eta = NULL;
TH2D* h2d_nreco_in_p_vs_eta = NULL;
TH2D* h2d_ngen_in_p_vs_eta = NULL;

double eff_tracking[pbin][etabin] = {0};
double dup_tracking[pbin][etabin] = {0};
double fake_tracking[pbin][etabin] = {0};

const char* seeding_name[2] = {"Realistic seeding", "Truth seeding"};

void mcs(int number = 0, int xp = 0, int yp = 0, int x = 400, int y = 400, double pleft = 0.17, double pright = 0.1, double ptop = 0.1, double pbot = 0.13){
  char name[100];
  sprintf(name,"cc%d",number);
  TPad *mpad;
  if(gROOT->GetListOfCanvases()->FindObject(name)==NULL){
    TCanvas* c = new TCanvas(name,name, xp, yp, x, y);
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    c->cd();
    sprintf(name,"mpad%d",number);  
    mpad = new TPad(name,name,0.02,0.02,0.99,0.99,0,0,0);
    mpad->SetTickx();
    mpad->SetTicky();
    mpad->SetTopMargin(ptop);
    mpad->SetBottomMargin(pbot);
    mpad->SetLeftMargin(pleft);
    mpad->SetRightMargin(pright);
    mpad->Draw();
    mpad->cd();
  c->Modified();
  c->Update();
  }
  else {
    TCanvas *tmp = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(name);
    tmp->cd();
    mpad->cd();
  tmp->Modified();
  tmp->Update();  
  }
  gSystem->ProcessEvents();
}
void mclogz(int number = 0, int xp = 0, int yp = 0, int x = 400, int y = 400, double pleft = 0.17, double pright = 0.12, double ptop = 0.1, double pbot = 0.13){
  char name[100];
  sprintf(name,"cc%d",number);
  TPad *mpad;
  if(gROOT->GetListOfCanvases()->FindObject(name)==NULL){
    TCanvas* c = new TCanvas(name,name, xp, yp, x, y);
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    c->cd();
    sprintf(name,"mpad%d",number);  
    mpad = new TPad(name,name,0.02,0.02,0.97,0.97,0,0,0);
    mpad->SetTopMargin(ptop);
    mpad->SetBottomMargin(pbot);
    mpad->SetLeftMargin(pleft);
    mpad->SetRightMargin(pright);
    mpad->SetLogz(1);
    mpad->Draw();
    mpad->cd();
  c->Modified();
  c->Update();
  }
  else {
    TCanvas *tmp = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(name);
    tmp->cd();
    mpad->cd();
  tmp->Modified();
  tmp->Update();  
  }
  gSystem->ProcessEvents();
}

void fit_gaus(TH1D* temp_h1d, TF1* fit_gaus, double& peak_mean, double& peak_mean_err, double& peak_sigma, double& peak_sigma_err)
{
  // fit_gaus->SetLineColor(temp_h1d->GetLineColor());
  // fit_gaus->SetLineStyle(7);
  temp_h1d->Fit(fit_gaus,"0R");
  const int iterbin = 3;
  for (int iter = 0; iter < iterbin; ++iter)
  {
    if (fit_gaus!=NULL)
    {
      peak_mean = fit_gaus->GetParameter(1);
      peak_mean_err = fit_gaus->GetParError(1);
      peak_sigma = fit_gaus->GetParameter(2);
      peak_sigma_err = fit_gaus->GetParError(2);
      
      temp_h1d->Fit(fit_gaus,"0R","",peak_mean-3*peak_sigma,peak_mean+3*peak_sigma);
    }
  }
}

void plot_dp_vs_nhit_2D(const int ip = 0, const int ieta = 0, const int seed_option = 0)
{
  if (h2d_1to1_dp_vs_nhit[ip][ieta]->Integral()==0) return;

  mclogz(cno++);
  {
    float plot_yrange_lo = 0;
    float plot_yrange_hi = 10;
    float plot_xrange_lo = -0.2;
    float plot_xrange_hi = 0.2;

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("#Delta p/p");
    htemp.GetYaxis()->SetTitle("# of hits");
    myhset(&htemp,1.2,1.6,0.05,0.045);

    TLegend leg(0.50,0.71,0.84,0.83);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.033);
    leg.SetMargin(0.2);
    leg.SetFillStyle(0);

    h2d_1to1_dp_vs_nhit[ip][ieta]->Draw("colz");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.03);
    tl->DrawLatexNDC(0.2,0.92,Form("Single particle events, p [%.1f, %.1f) GeV/c, #eta [%.1f, %.1f)",p_lo[ip],p_hi[ip],eta_lo[ieta],eta_hi[ieta]));
    tl->SetTextSize(0.03);
    // tl->DrawLatexNDC(0.2,0.85,Form("# of reconstructed (generated) particles %.0f (%.0f)",h2d_1to1_dp_vs_nhit[ip][ieta]->Integral(),h1d_gen_npart[ip][ieta]->GetEntries()));
    tl->DrawLatexNDC(0.2,0.85,Form("%s",seeding_name[seed_option]));

    gROOT->ProcessLine( Form("cc%d->Print(\"figs/dp_vs_nhit/h2d_dp_vs_nhit_eta%d_p%d.pdf\")", cno-1, ieta, ip) );
  }
}

void plot_dp_1D(const int ip = 0, const int ieta = 0, const int seed_option = 0)
{
  // h2d_1to1_dp_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(3,10);
  h1d_1to1_dp[ip][ieta] = (TH1D*)h2d_1to1_dp_vs_nhit[ip][ieta]->ProjectionX( Form("h1d_1to1_dp_eta%d_p%d",ieta,ip) );
  h1d_1to1_dp[ip][ieta]->SetName( Form("h1d_1to1_dp_eta%d_p%d",ieta,ip) );
  h2d_1to1_dp_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(0,10);

  // h2d_Nto1_dp_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(3,10);
  h1d_Nto1_dp[ip][ieta] = (TH1D*)h2d_Nto1_dp_vs_nhit[ip][ieta]->ProjectionX( Form("h1d_Nto1_dp_eta%d_p%d",ieta,ip) );
  h1d_Nto1_dp[ip][ieta]->SetName( Form("h1d_Nto1_dp_eta%d_p%d",ieta,ip) );
  h2d_Nto1_dp_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(0,10);

  // h2d_0to1_dp_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(3,10);
  h1d_0to1_dp[ip][ieta] = (TH1D*)h2d_0to1_dp_vs_nhit[ip][ieta]->ProjectionX( Form("h1d_0to1_dp_eta%d_p%d",ieta,ip) );
  h1d_0to1_dp[ip][ieta]->SetName( Form("h1d_0to1_dp_eta%d_p%d",ieta,ip) );
  h2d_0to1_dp_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(0,10);

  // cout << ip << " " << ieta << " w/ integral " << h1d_1to1_dp[ip][ieta]->Integral() << endl;

  // h1d_1to1_dp[ip][ieta]->Rebin(2);
  // h1d_Nto1_dp[ip][ieta]->Rebin(2);
  // h1d_0to1_dp[ip][ieta]->Rebin(2);

  if (h1d_1to1_dp[ip][ieta]->Integral()==0) return;

  mcs(cno++);
  {
    float plot_yrange_lo = 0;
    float plot_yrange_hi = 1.7*h1d_1to1_dp[ip][ieta]->GetMaximum();
    float plot_xrange_lo = -0.2;
    float plot_xrange_hi = 0.2;

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("#Delta p/p");
    htemp.GetYaxis()->SetTitle("Counts");
    myhset(&htemp,1.2,1.6,0.05,0.045);

    TLegend leg(0.50,0.71,0.84,0.83);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.033);
    leg.SetMargin(0.2);
    leg.SetFillStyle(0);

    h1d_1to1_dp[ip][ieta]->SetLineWidth(2);
    h1d_1to1_dp[ip][ieta]->SetLineColor(kBlue);
    h1d_1to1_dp[ip][ieta]->SetMarkerColor(kBlue);
    h1d_1to1_dp[ip][ieta]->Draw("hsame");
    leg.AddEntry(h1d_1to1_dp[ip][ieta],"good (matched) tracks","l");

    h1d_Nto1_dp[ip][ieta]->SetLineWidth(2);
    h1d_Nto1_dp[ip][ieta]->SetLineColor(kGreen+1);
    h1d_Nto1_dp[ip][ieta]->SetMarkerColor(kGreen+1);
    h1d_Nto1_dp[ip][ieta]->Draw("hsame");
    leg.AddEntry(h1d_Nto1_dp[ip][ieta],"duplicate tracks","l");

    h1d_0to1_dp[ip][ieta]->SetLineWidth(2);
    h1d_0to1_dp[ip][ieta]->SetLineColor(kGray+1);
    h1d_0to1_dp[ip][ieta]->SetMarkerColor(kGray+1);
    h1d_0to1_dp[ip][ieta]->Draw("hsame");
    leg.AddEntry(h1d_0to1_dp[ip][ieta],"fake (unmatched) tracks","l");

    leg.Draw("same");

    TF1* fit_gaus_full = new TF1("fit_gaus_full","gaus",plot_xrange_lo,plot_xrange_hi);
    fit_gaus_full->SetParameter(1,0);
    fit_gaus_full->SetParameter(2,0.05);
    double temp_preco_mean = -9999; 
    double temp_preco_mean_err = 0;
    double temp_preco_sigma = -9999;
    double temp_preco_sigma_err = 0;
    if (h1d_1to1_dp[ip][ieta]) fit_gaus(h1d_1to1_dp[ip][ieta], fit_gaus_full, temp_preco_mean, temp_preco_mean_err, temp_preco_sigma, temp_preco_sigma_err);
    if (fit_gaus_full) fit_gaus_full->Draw("same");

    int good_mom_range_lo = h1d_1to1_dp[ip][ieta]->FindBin(temp_preco_mean[ip][ieta]-5*temp_preco_sigma[ip][ieta]);
    int good_mom_range_hi = h1d_1to1_dp[ip][ieta]->FindBin(temp_preco_mean[ip][ieta]+5*temp_preco_sigma[ip][ieta]);
    float npart_good_mom = h1d_1to1_dp[ip][ieta]->GetEntries();//h1d_1to1_dp[ip][ieta]->Integral(good_mom_range_lo,good_mom_range_hi);
    float npart_thrown = h1d_gen_npart[ip][ieta]->GetEntries();//h1d_gen_npart[ip][ieta]->Integral();
    eff_tracking[ip][ieta] = npart_good_mom/npart_thrown;

    float ntracks_dup = h1d_Nto1_dp[ip][ieta]->GetEntries();//h1d_Nto1_dp[ip][ieta]->Integral(good_mom_range_lo,good_mom_range_hi);
    float ntracks_fake = h1d_0to1_dp[ip][ieta]->GetEntries();//h1d_0to1_dp[ip][ieta]->Integral(good_mom_range_lo,good_mom_range_hi);
    dup_tracking[ip][ieta] = ntracks_dup/npart_thrown;//ntracks_dup/npart_good_mom;
    fake_tracking[ip][ieta] = ntracks_fake/npart_thrown;//ntracks_fake/npart_good_mom;

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.03);
    tl->DrawLatexNDC(0.2,0.92,Form("Single particle events, p [%.1f, %.1f) GeV/c, #eta [%.1f, %.1f)",p_lo[ip],p_hi[ip],eta_lo[ieta],eta_hi[ieta]));
    tl->SetTextSize(0.03);
    // tl->DrawLatexNDC(0.2,0.85,Form("# of reconstructed (generated) particles %.0f (%.0f)",h2d_1to1_dp_vs_nhit[ip][ieta]->Integral(),h1d_gen_npart[ip][ieta]->GetEntries()));
    tl->DrawLatexNDC(0.2,0.85,Form("%s",seeding_name[seed_option]));

    tl->SetTextSize(0.03);
    tl->SetTextColor(kRed);
    if (fit_gaus_full)
    {
      tl->DrawLatexNDC(0.2,0.80,Form("#mu = %.2f (#chi^{2}/ndf = %.1f)",temp_preco_mean,fit_gaus_full->GetChisquare()/fit_gaus_full->GetNDF()));
      tl->DrawLatexNDC(0.2,0.75,Form("#sigma = %.1f%%",temp_preco_sigma*100));
    }
    else
    {
      tl->DrawLatexNDC(0.2,0.80,"#mu = invalid");
      tl->DrawLatexNDC(0.2,0.75,"#sigma = invalid");
    }
    // tl->DrawLatexNDC(0.2,0.70,Form("effciency = %.2e/%.2e = %.1f%%",npart_good_mom,npart_thrown,100*eff_tracking[ip][ieta]));
    tl->DrawLatexNDC(0.2,0.70,Form("effciency = %.1f%%",100*eff_tracking[ip][ieta]));
    tl->DrawLatexNDC(0.2,0.65,Form("duplicate rate = %.1f%%",100*dup_tracking[ip][ieta]));
    tl->DrawLatexNDC(0.2,0.60,Form("fake rate = %.1f%%",100*fake_tracking[ip][ieta]));

    if (temp_preco_sigma>0) h2d_dp_in_p_vs_eta->SetBinContent(ip+1,ieta+1,temp_preco_sigma*100);
    else h2d_dp_in_p_vs_eta->SetBinContent(ip+1,ieta+1,0);

    int good_track_bin_lo = h1d_1to1_dp[ip][ieta]->FindBin(temp_preco_mean-5*temp_preco_sigma);
    int good_track_bin_hi = h1d_1to1_dp[ip][ieta]->FindBin(temp_preco_mean+5*temp_preco_sigma);
    if (temp_preco_sigma>0) h2d_nreco_in_p_vs_eta->SetBinContent(ip+1,ieta+1,h1d_1to1_dp[ip][ieta]->Integral(good_track_bin_lo,good_track_bin_hi));
    else h2d_nreco_in_p_vs_eta->SetBinContent(ip+1,ieta+1,0);

    gROOT->ProcessLine( Form("cc%d->Print(\"figs/dp_fit/h1d_dp_fit_eta%d_p%d.pdf\")", cno-1, ieta, ip) );
  }
}

void plot_chi2_1D(const int ip = 0, const int ieta = 0, const int seed_option = 0)
{
  // h2d_1to1_chi2_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(3,10);
  h1d_1to1_chi2[ip][ieta] = (TH1D*)h2d_1to1_chi2_vs_nhit[ip][ieta]->ProjectionX( Form("h1d_1to1_chi2_eta%d_p%d",ieta,ip) );
  h1d_1to1_chi2[ip][ieta]->SetName( Form("h1d_1to1_chi2_eta%d_p%d",ieta,ip) );
  h2d_1to1_chi2_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(0,10);

  // h2d_Nto1_chi2_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(3,10);
  h1d_Nto1_chi2[ip][ieta] = (TH1D*)h2d_Nto1_chi2_vs_nhit[ip][ieta]->ProjectionX( Form("h1d_Nto1_chi2_eta%d_p%d",ieta,ip) );
  h1d_Nto1_chi2[ip][ieta]->SetName( Form("h1d_Nto1_chi2_eta%d_p%d",ieta,ip) );
  h2d_Nto1_chi2_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(0,10);

  // h2d_0to1_chi2_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(3,10);
  h1d_0to1_chi2[ip][ieta] = (TH1D*)h2d_0to1_chi2_vs_nhit[ip][ieta]->ProjectionX( Form("h1d_0to1_chi2_eta%d_p%d",ieta,ip) );
  h1d_0to1_chi2[ip][ieta]->SetName( Form("h1d_0to1_chi2_eta%d_p%d",ieta,ip) );
  h2d_0to1_chi2_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(0,10);

  h1d_1to1_chi2[ip][ieta]->Rebin(2);
  h1d_Nto1_chi2[ip][ieta]->Rebin(2);
  h1d_0to1_chi2[ip][ieta]->Rebin(2);

  if (h1d_1to1_chi2[ip][ieta]->Integral()==0) return;

  mcs(cno++);
  {
    float plot_yrange_lo = 0;
    float plot_yrange_hi = 1.7*h1d_1to1_chi2[ip][ieta]->GetMaximum();
    float plot_xrange_lo = 0;
    float plot_xrange_hi = 50;

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("#chi^{2}");
    htemp.GetYaxis()->SetTitle("Counts");
    myhset(&htemp,1.2,1.6,0.05,0.045);

    TLegend leg(0.50,0.71,0.84,0.83);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.033);
    leg.SetMargin(0.2);
    leg.SetFillStyle(0);

    h1d_1to1_chi2[ip][ieta]->SetLineWidth(2);
    h1d_1to1_chi2[ip][ieta]->SetLineColor(kBlue);
    h1d_1to1_chi2[ip][ieta]->SetMarkerColor(kBlue);
    h1d_1to1_chi2[ip][ieta]->Draw("hsame");
    leg.AddEntry(h1d_1to1_chi2[ip][ieta],"good (matched) tracks","l");

    h1d_Nto1_chi2[ip][ieta]->SetLineWidth(2);
    h1d_Nto1_chi2[ip][ieta]->SetLineColor(kGreen+1);
    h1d_Nto1_chi2[ip][ieta]->SetMarkerColor(kGreen+1);
    h1d_Nto1_chi2[ip][ieta]->Draw("hsame");
    leg.AddEntry(h1d_Nto1_chi2[ip][ieta],"duplicate tracks","l");

    h1d_0to1_chi2[ip][ieta]->SetLineWidth(2);
    h1d_0to1_chi2[ip][ieta]->SetLineColor(kGray+1);
    h1d_0to1_chi2[ip][ieta]->SetMarkerColor(kGray+1);
    h1d_0to1_chi2[ip][ieta]->Draw("hsame");
    leg.AddEntry(h1d_0to1_chi2[ip][ieta],"fake (unmatched) tracks","l");

    leg.Draw("same");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.03);
    tl->DrawLatexNDC(0.2,0.92,Form("Single particle events, p [%.1f, %.1f) GeV/c, #eta [%.1f, %.1f)",p_lo[ip],p_hi[ip],eta_lo[ieta],eta_hi[ieta]));
    tl->DrawLatexNDC(0.2,0.85,Form("%s",seeding_name[seed_option]));

    gROOT->ProcessLine( Form("cc%d->Print(\"figs/dp_fit/h1d_chi2_eta%d_p%d.pdf\")", cno-1, ieta, ip) );
  }
}

void plot_nhit_1D(const int ip = 0, const int ieta = 0, const int seed_option = 0)
{
  // h2d_1to1_dp_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(3,10);
  h1d_1to1_nhit[ip][ieta] = (TH1D*)h2d_1to1_dp_vs_nhit[ip][ieta]->ProjectionY( Form("h1d_1to1_nhit_eta%d_p%d",ieta,ip) );
  h1d_1to1_nhit[ip][ieta]->SetName( Form("h1d_1to1_nhit_eta%d_p%d",ieta,ip) );
  h2d_1to1_dp_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(0,10);

  // h2d_Nto1_dp_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(3,10);
  h1d_Nto1_nhit[ip][ieta] = (TH1D*)h2d_Nto1_dp_vs_nhit[ip][ieta]->ProjectionY( Form("h1d_Nto1_nhit_eta%d_p%d",ieta,ip) );
  h1d_Nto1_nhit[ip][ieta]->SetName( Form("h1d_Nto1_nhit_eta%d_p%d",ieta,ip) );
  h2d_Nto1_dp_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(0,10);

  // h2d_0to1_dp_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(3,10);
  h1d_0to1_nhit[ip][ieta] = (TH1D*)h2d_0to1_dp_vs_nhit[ip][ieta]->ProjectionY( Form("h1d_0to1_nhit_eta%d_p%d",ieta,ip) );
  h1d_0to1_nhit[ip][ieta]->SetName( Form("h1d_0to1_nhit_eta%d_p%d",ieta,ip) );
  h2d_0to1_dp_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(0,10);

  if (h1d_1to1_nhit[ip][ieta]->Integral()==0) return;

  mcs(cno++);
  {
    float plot_yrange_lo = 0;
    float plot_yrange_hi = 1.7*h1d_1to1_nhit[ip][ieta]->GetMaximum();
    float plot_xrange_lo = -0.5;
    float plot_xrange_hi = 10.5;

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("N_{hits}");
    htemp.GetYaxis()->SetTitle("Counts");
    myhset(&htemp,1.2,1.6,0.05,0.045);

    TLegend leg(0.50,0.71,0.84,0.83);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.033);
    leg.SetMargin(0.2);
    leg.SetFillStyle(0);

    h1d_1to1_nhit[ip][ieta]->SetLineWidth(2);
    h1d_1to1_nhit[ip][ieta]->SetLineColor(kBlue);
    h1d_1to1_nhit[ip][ieta]->SetMarkerColor(kBlue);
    h1d_1to1_nhit[ip][ieta]->Draw("hsame");
    leg.AddEntry(h1d_1to1_nhit[ip][ieta],"good (matched) tracks","l");

    h1d_Nto1_nhit[ip][ieta]->SetLineWidth(2);
    h1d_Nto1_nhit[ip][ieta]->SetLineColor(kGreen+1);
    h1d_Nto1_nhit[ip][ieta]->SetMarkerColor(kGreen+1);
    h1d_Nto1_nhit[ip][ieta]->Draw("hsame");
    leg.AddEntry(h1d_Nto1_nhit[ip][ieta],"duplicate tracks","l");

    h1d_0to1_nhit[ip][ieta]->SetLineWidth(2);
    h1d_0to1_nhit[ip][ieta]->SetLineColor(kGray+1);
    h1d_0to1_nhit[ip][ieta]->SetMarkerColor(kGray+1);
    h1d_0to1_nhit[ip][ieta]->Draw("hsame");
    leg.AddEntry(h1d_0to1_nhit[ip][ieta],"fake (unmatched) tracks","l");

    leg.Draw("same");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.03);
    tl->DrawLatexNDC(0.2,0.92,Form("Single particle events, p [%.1f, %.1f) GeV/c, #eta [%.1f, %.1f)",p_lo[ip],p_hi[ip],eta_lo[ieta],eta_hi[ieta]));
    tl->DrawLatexNDC(0.2,0.85,Form("%s",seeding_name[seed_option]));

    gROOT->ProcessLine( Form("cc%d->Print(\"figs/dp_fit/h1d_nhit_eta%d_p%d.pdf\")", cno-1, ieta, ip) );
  }
}

void plot_dphi_1D(const int ip = 0, const int ieta = 0, const int seed_option = 0)
{
  // h2d_1to1_dphi_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(3,10);
  h1d_1to1_dphi[ip][ieta] = (TH1D*)h2d_1to1_dphi_vs_nhit[ip][ieta]->ProjectionX( Form("h1d_1to1_dphi_eta%d_p%d",ieta,ip) );
  h1d_1to1_dphi[ip][ieta]->SetName( Form("h1d_1to1_dphi_eta%d_p%d",ieta,ip) );
  h2d_1to1_dphi_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(0,10);

  // h2d_Nto1_dphi_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(3,10);
  h1d_Nto1_dphi[ip][ieta] = (TH1D*)h2d_Nto1_dphi_vs_nhit[ip][ieta]->ProjectionX( Form("h1d_Nto1_dphi_eta%d_p%d",ieta,ip) );
  h1d_Nto1_dphi[ip][ieta]->SetName( Form("h1d_Nto1_dphi_eta%d_p%d",ieta,ip) );
  h2d_Nto1_dphi_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(0,10);

  // h2d_0to1_dphi_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(3,10);
  h1d_0to1_dphi[ip][ieta] = (TH1D*)h2d_0to1_dphi_vs_nhit[ip][ieta]->ProjectionX( Form("h1d_0to1_dphi_eta%d_p%d",ieta,ip) );
  h1d_0to1_dphi[ip][ieta]->SetName( Form("h1d_0to1_dphi_eta%d_p%d",ieta,ip) );
  h2d_0to1_dphi_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(0,10);

  h1d_1to1_dphi[ip][ieta]->Rebin(5);
  h1d_Nto1_dphi[ip][ieta]->Rebin(5);
  h1d_0to1_dphi[ip][ieta]->Rebin(5);

  if (h1d_1to1_dphi[ip][ieta]->Integral()==0) return;

  mcs(cno++);
  {
    float plot_yrange_lo = 0;
    float plot_yrange_hi = 1.7*h1d_1to1_dphi[ip][ieta]->GetMaximum();
    float plot_xrange_lo = h1d_1to1_dphi[ip][ieta]->GetBinCenter(1);
    float plot_xrange_hi = h1d_1to1_dphi[ip][ieta]->GetBinCenter(h1d_1to1_dphi[ip][ieta]->GetNbinsX());

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("#Delta#phi [rad]");
    htemp.GetYaxis()->SetTitle("Counts");
    myhset(&htemp,1.2,1.6,0.05,0.045);

    TLegend leg(0.50,0.71,0.84,0.83);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.033);
    leg.SetMargin(0.2);
    leg.SetFillStyle(0);

    h1d_1to1_dphi[ip][ieta]->SetLineWidth(2);
    h1d_1to1_dphi[ip][ieta]->SetLineColor(kBlue);
    h1d_1to1_dphi[ip][ieta]->SetMarkerColor(kBlue);
    h1d_1to1_dphi[ip][ieta]->Draw("hsame");
    leg.AddEntry(h1d_1to1_dphi[ip][ieta],"good (matched) tracks","l");

    h1d_Nto1_dphi[ip][ieta]->SetLineWidth(2);
    h1d_Nto1_dphi[ip][ieta]->SetLineColor(kGreen+1);
    h1d_Nto1_dphi[ip][ieta]->SetMarkerColor(kGreen+1);
    h1d_Nto1_dphi[ip][ieta]->Draw("hsame");
    leg.AddEntry(h1d_Nto1_dphi[ip][ieta],"duplicate tracks","l");

    h1d_0to1_dphi[ip][ieta]->SetLineWidth(2);
    h1d_0to1_dphi[ip][ieta]->SetLineColor(kGray+1);
    h1d_0to1_dphi[ip][ieta]->SetMarkerColor(kGray+1);
    h1d_0to1_dphi[ip][ieta]->Draw("hsame");
    leg.AddEntry(h1d_0to1_dphi[ip][ieta],"fake (unmatched) tracks","l");

    leg.Draw("same");

    TF1* fit_gaus_full = new TF1("fit_gaus_full","gaus",plot_xrange_lo,plot_xrange_hi);
    fit_gaus_full->SetParameter(1,0);
    fit_gaus_full->SetParameter(2,0.05);
    double temp_phireco_mean = -9999; 
    double temp_phireco_mean_err = 0;
    double temp_phireco_sigma = -9999;
    double temp_phireco_sigma_err = 0;
    if (h1d_1to1_dphi[ip][ieta]) fit_gaus(h1d_1to1_dphi[ip][ieta], fit_gaus_full, temp_phireco_mean, temp_phireco_mean_err, temp_phireco_sigma, temp_phireco_sigma_err);
    if (fit_gaus_full) fit_gaus_full->Draw("same");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.03);
    tl->DrawLatexNDC(0.2,0.92,Form("Single particle events, p [%.1f, %.1f) GeV/c, #eta [%.1f, %.1f)",p_lo[ip],p_hi[ip],eta_lo[ieta],eta_hi[ieta]));
    tl->SetTextSize(0.03);
    // tl->DrawLatexNDC(0.2,0.85,Form("# of reconstructed (generated) particles %.0f (%.0f)",h2d_1to1_dphi_vs_nhit[ip][ieta]->Integral(),h1d_gen_npart[ip][ieta]->GetEntries()));
    tl->DrawLatexNDC(0.2,0.85,Form("%s",seeding_name[seed_option]));

    tl->SetTextSize(0.03);
    tl->SetTextColor(kRed);
    if (fit_gaus_full)
    {
      tl->DrawLatexNDC(0.2,0.80,Form("#mu = %.2f (#chi^{2}/ndf = %.1f)",temp_phireco_mean,fit_gaus_full->GetChisquare()/fit_gaus_full->GetNDF()));
      tl->DrawLatexNDC(0.2,0.75,Form("#sigma = %.1f mrad",temp_phireco_sigma*1000));
    }
    else
    {
      tl->DrawLatexNDC(0.2,0.80,"#mu = invalid");
      tl->DrawLatexNDC(0.2,0.75,"#sigma = invalid");
    }

    if (temp_phireco_sigma>0) h2d_dphi_in_p_vs_eta->SetBinContent(ip+1,ieta+1,temp_phireco_sigma*1000);
    else h2d_dphi_in_p_vs_eta->SetBinContent(ip+1,ieta+1,0);

    gROOT->ProcessLine( Form("cc%d->Print(\"figs/dphi_fit/h1d_dphi_fit_eta%d_p%d.pdf\")", cno-1, ieta, ip) );
  }
}

void plot_dtheta_1D(const int ip = 0, const int ieta = 0, const int seed_option = 0)
{
  // h2d_1to1_dtheta_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(3,10);
  h1d_1to1_dtheta[ip][ieta] = (TH1D*)h2d_1to1_dtheta_vs_nhit[ip][ieta]->ProjectionX( Form("h1d_1to1_dtheta_eta%d_p%d",ieta,ip) );
  h1d_1to1_dtheta[ip][ieta]->SetName( Form("h1d_1to1_dtheta_eta%d_p%d",ieta,ip) );
  h2d_1to1_dtheta_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(0,10);

  // h2d_Nto1_dtheta_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(3,10);
  h1d_Nto1_dtheta[ip][ieta] = (TH1D*)h2d_Nto1_dtheta_vs_nhit[ip][ieta]->ProjectionX( Form("h1d_Nto1_dtheta_eta%d_p%d",ieta,ip) );
  h1d_Nto1_dtheta[ip][ieta]->SetName( Form("h1d_Nto1_dtheta_eta%d_p%d",ieta,ip) );
  h2d_Nto1_dtheta_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(0,10);

  // h2d_0to1_dtheta_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(3,10);
  h1d_0to1_dtheta[ip][ieta] = (TH1D*)h2d_0to1_dtheta_vs_nhit[ip][ieta]->ProjectionX( Form("h1d_0to1_dtheta_eta%d_p%d",ieta,ip) );
  h1d_0to1_dtheta[ip][ieta]->SetName( Form("h1d_0to1_dtheta_eta%d_p%d",ieta,ip) );
  h2d_0to1_dtheta_vs_nhit[ip][ieta]->GetYaxis()->SetRangeUser(0,10);

  h1d_1to1_dtheta[ip][ieta]->Rebin(2);
  h1d_Nto1_dtheta[ip][ieta]->Rebin(2);
  h1d_0to1_dtheta[ip][ieta]->Rebin(2);

  if (h1d_1to1_dtheta[ip][ieta]->Integral()==0) return;

  mcs(cno++);
  {
    float plot_yrange_lo = 0;
    float plot_yrange_hi = 1.7*h1d_1to1_dtheta[ip][ieta]->GetMaximum();
    float plot_xrange_lo = h1d_1to1_dtheta[ip][ieta]->GetBinCenter(1);
    float plot_xrange_hi = h1d_1to1_dtheta[ip][ieta]->GetBinCenter(h1d_1to1_dtheta[ip][ieta]->GetNbinsX());

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("#Delta#theta [rad]");
    htemp.GetYaxis()->SetTitle("Counts");
    myhset(&htemp,1.2,1.6,0.05,0.045);

    TLegend leg(0.50,0.71,0.84,0.83);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.033);
    leg.SetMargin(0.2);
    leg.SetFillStyle(0);

    h1d_1to1_dtheta[ip][ieta]->SetLineWidth(2);
    h1d_1to1_dtheta[ip][ieta]->SetLineColor(kBlue);
    h1d_1to1_dtheta[ip][ieta]->SetMarkerColor(kBlue);
    h1d_1to1_dtheta[ip][ieta]->Draw("hsame");
    leg.AddEntry(h1d_1to1_dtheta[ip][ieta],"good (matched) tracks","l");

    h1d_Nto1_dtheta[ip][ieta]->SetLineWidth(2);
    h1d_Nto1_dtheta[ip][ieta]->SetLineColor(kGreen+1);
    h1d_Nto1_dtheta[ip][ieta]->SetMarkerColor(kGreen+1);
    h1d_Nto1_dtheta[ip][ieta]->Draw("hsame");
    leg.AddEntry(h1d_Nto1_dtheta[ip][ieta],"duplicate tracks","l");

    h1d_0to1_dtheta[ip][ieta]->SetLineWidth(2);
    h1d_0to1_dtheta[ip][ieta]->SetLineColor(kGray+1);
    h1d_0to1_dtheta[ip][ieta]->SetMarkerColor(kGray+1);
    h1d_0to1_dtheta[ip][ieta]->Draw("hsame");
    leg.AddEntry(h1d_0to1_dtheta[ip][ieta],"fake (unmatched) tracks","l");

    leg.Draw("same");

    TF1* fit_gaus_full = new TF1("fit_gaus_full","gaus",plot_xrange_lo,plot_xrange_hi);
    fit_gaus_full->SetParameter(1,0);
    fit_gaus_full->SetParameter(2,0.05);
    double temp_phireco_mean = -9999; 
    double temp_phireco_mean_err = 0;
    double temp_phireco_sigma = -9999;
    double temp_phireco_sigma_err = 0;
    if (h1d_1to1_dtheta[ip][ieta]) fit_gaus(h1d_1to1_dtheta[ip][ieta], fit_gaus_full, temp_phireco_mean, temp_phireco_mean_err, temp_phireco_sigma, temp_phireco_sigma_err);
    if (fit_gaus_full) fit_gaus_full->Draw("same");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.03);
    tl->DrawLatexNDC(0.2,0.92,Form("Single particle events, p [%.1f, %.1f) GeV/c, #eta [%.1f, %.1f)",p_lo[ip],p_hi[ip],eta_lo[ieta],eta_hi[ieta]));
    tl->SetTextSize(0.03);
    // tl->DrawLatexNDC(0.2,0.85,Form("# of reconstructed (generated) particles %.0f (%.0f)",h2d_1to1_dtheta_vs_nhit[ip][ieta]->Integral(),h1d_gen_npart[ip][ieta]->GetEntries()));
    tl->DrawLatexNDC(0.2,0.85,Form("%s",seeding_name[seed_option]));

    tl->SetTextSize(0.03);
    tl->SetTextColor(kRed);
    if (fit_gaus_full)
    {
      tl->DrawLatexNDC(0.2,0.80,Form("#mu = %.2f (#chi^{2}/ndf = %.1f)",temp_phireco_mean,fit_gaus_full->GetChisquare()/fit_gaus_full->GetNDF()));
      tl->DrawLatexNDC(0.2,0.75,Form("#sigma = %.1f mrad",temp_phireco_sigma*1000));
    }
    else
    {
      tl->DrawLatexNDC(0.2,0.80,"#mu = invalid");
      tl->DrawLatexNDC(0.2,0.75,"#sigma = invalid");
    }

    if (temp_phireco_sigma>0) h2d_dtheta_in_p_vs_eta->SetBinContent(ip+1,ieta+1,temp_phireco_sigma*1000);
    else h2d_dtheta_in_p_vs_eta->SetBinContent(ip+1,ieta+1,0);

    gROOT->ProcessLine( Form("cc%d->Print(\"figs/dtheta_fit/h1d_dtheta_fit_eta%d_p%d.pdf\")", cno-1, ieta, ip) );
  }
}

void dump_seeding_graph(const char* inFile = "particles_hists.root", const int seed_option = 0)
{
  mcs(-1);

  // input 
  TFile fin(inFile,"read");
  // fin.ls();

  for (int ip = 0; ip < pbin; ++ip)
  {
    for (int ieta = 0; ieta < etabin; ++ieta)
    {
      h1d_gen_npart[ip][ieta] = (TH1D*)fin.Get(Form("h1d_gen_npart_%d_%d",ip,ieta));

      h2d_1to1_dp_vs_nhit[ip][ieta] = (TH2D*)fin.Get(Form("h2d_1to1_dp_vs_nhit_%d_%d",ip,ieta));
      h2d_1to1_dphi_vs_nhit[ip][ieta] = (TH2D*)fin.Get(Form("h2d_1to1_dphi_vs_nhit_%d_%d",ip,ieta));
      h2d_1to1_dtheta_vs_nhit[ip][ieta] = (TH2D*)fin.Get(Form("h2d_1to1_dtheta_vs_nhit_%d_%d",ip,ieta));
      h2d_1to1_chi2_vs_nhit[ip][ieta] = (TH2D*)fin.Get(Form("h2d_1to1_chi2_vs_nhit_%d_%d",ip,ieta));

      h2d_Nto1_dp_vs_nhit[ip][ieta] = (TH2D*)fin.Get(Form("h2d_Nto1_dp_vs_nhit_%d_%d",ip,ieta));
      h2d_Nto1_dphi_vs_nhit[ip][ieta] = (TH2D*)fin.Get(Form("h2d_Nto1_dphi_vs_nhit_%d_%d",ip,ieta));
      h2d_Nto1_dtheta_vs_nhit[ip][ieta] = (TH2D*)fin.Get(Form("h2d_Nto1_dtheta_vs_nhit_%d_%d",ip,ieta));
      h2d_Nto1_chi2_vs_nhit[ip][ieta] = (TH2D*)fin.Get(Form("h2d_Nto1_chi2_vs_nhit_%d_%d",ip,ieta));

      h2d_0to1_dp_vs_nhit[ip][ieta] = (TH2D*)fin.Get(Form("h2d_0to1_dp_vs_nhit_%d_%d",ip,ieta));
      h2d_0to1_dphi_vs_nhit[ip][ieta] = (TH2D*)fin.Get(Form("h2d_0to1_dphi_vs_nhit_%d_%d",ip,ieta));
      h2d_0to1_dtheta_vs_nhit[ip][ieta] = (TH2D*)fin.Get(Form("h2d_0to1_dtheta_vs_nhit_%d_%d",ip,ieta));
      h2d_0to1_chi2_vs_nhit[ip][ieta] = (TH2D*)fin.Get(Form("h2d_0to1_chi2_vs_nhit_%d_%d",ip,ieta));
    }
  }

  cout << "p_edges = {";
  double p_edges[pbin+1] = {0};
  for (int ip = 0; ip < pbin+1; ++ip)
  {
    p_edges[ip] = p_lo[ip];
    if (ip == pbin) p_edges[ip] = p_hi[ip-1];
    cout << p_edges[ip];
    if (ip!=pbin) cout << ", ";
  }
  cout << "};" << endl;
  cout << "eta_edges = {";
  double eta_edges[etabin+1] = {0};
  for (int ieta = 0; ieta < etabin+1; ++ieta)
  {
    eta_edges[ieta] = eta_lo[ieta];
    if (ieta == etabin) eta_edges[ieta] = eta_hi[ieta-1];
    cout << eta_edges[ieta];
    if (ieta!=etabin) cout << ", ";
  }
  cout << "};" << endl;

  h2d_dp_in_p_vs_eta = new TH2D("h2d_dp_in_p_vs_eta","",pbin,p_edges,etabin,eta_edges); // z axis %
  h2d_dphi_in_p_vs_eta = new TH2D("h2d_dphi_in_p_vs_eta","",pbin,p_edges,etabin,eta_edges); // z axis mrad
  h2d_dtheta_in_p_vs_eta = new TH2D("h2d_dtheta_in_p_vs_eta","",pbin,p_edges,etabin,eta_edges); // z axis mrad

  h2d_nreco_in_p_vs_eta = new TH2D("h2d_nreco_in_p_vs_eta","",pbin,p_edges,etabin,eta_edges);
  h2d_ngen_in_p_vs_eta = new TH2D("h2d_ngen_in_p_vs_eta","",pbin,p_edges,etabin,eta_edges);

  // setup generation hists
  for (int ip = 0; ip < pbin; ++ip)
  {
    for (int ieta = 0; ieta < etabin; ++ieta)
    {
      h2d_ngen_in_p_vs_eta->SetBinContent(ip+1,ieta+1,h1d_gen_npart[ip][ieta]->Integral());
    }
  }

  // setup other hists
  for (int ip = 0; ip < pbin; ++ip)
  {
    for (int ieta = 0; ieta < etabin; ++ieta)
    {
      cout << ip << " " << ieta << endl;
      plot_dp_1D(ip, ieta, seed_option);
      plot_chi2_1D(ip, ieta, seed_option);
      plot_nhit_1D(ip, ieta, seed_option);

      plot_dp_vs_nhit_2D(ip, ieta, seed_option);

      plot_dphi_1D(ip, ieta, seed_option);
      plot_dtheta_1D(ip, ieta, seed_option);
    }
  }

  // plot_dp_2D();

  TFile* fout = new TFile("graph.root","recreate");
  h2d_dp_in_p_vs_eta->Write();
  h2d_dphi_in_p_vs_eta->Write();
  h2d_dtheta_in_p_vs_eta->Write();
  h2d_nreco_in_p_vs_eta->Write();
  h2d_ngen_in_p_vs_eta->Write();
  fout->Write();
  fout->Close();
}
