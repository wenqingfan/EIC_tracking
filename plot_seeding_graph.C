#include "eff_bins.h"
#include "DIS_kin.h"
using namespace std;

// kinematic hard constrains: 0<x<1, 0<y<1(0.01-0.98?) , Q2>0

static int cno = 0;

TH2D* h2d_ngen_in_p_vs_eta = NULL;
TH2D* h2d_nreco_in_p_vs_eta = NULL;
TH2D* h2d_eff_in_p_vs_eta = NULL;

TH2D* h2d_dp_in_p_vs_eta = NULL;
TH2D* h2d_dphi_in_p_vs_eta = NULL;
TH2D* h2d_dtheta_in_p_vs_eta = NULL;

TGraph* g_dp_in_p[etabin] = {0};
TGraph* g_dphi_in_p[etabin] = {0};
TGraph* g_dtheta_in_p[etabin] = {0};
TGraph* g_nreco_in_p[etabin] = {0};
TGraph* g_ngen_in_p[etabin] = {0};
TGraph* g_eff_in_p[etabin] = {0};

const char* seeding_name[2] = {"realistic seed", "truth seed"};

// ATHENA single track smearing
const int NETA_MAX = 500;
TH1F* Res_Handler = NULL;
TGraph *gmom_res[NETA_MAX];
TGraph *gdca_rphi_res[NETA_MAX];
TGraph *gdca_z_res[NETA_MAX];

void plot_eff_1D(const int option = 0)
{
  for (int ieta = 0; ieta < etabin; ++ieta)
  {
    double temp_p[pbin] = {0};
    double temp_eff[pbin] = {0};

    double temp_display_p[pbin] = {0};
    double temp_display_eff[pbin] = {0};
    int pbin_display = 0;
    bool flag_threshold = false;
    double p_threhold = 0, pt_threhold = 0;
    for (int ip = 0; ip < pbin; ++ip)
    {
      temp_p[ip] = h2d_eff_in_p_vs_eta->GetXaxis()->GetBinCenter(ip+1);
      temp_eff[ip] = h2d_eff_in_p_vs_eta->GetBinContent(ip+1,ieta+1);
      if (!flag_threshold && temp_eff[ip]>0.8)
      {
        if (ip>0 && temp_eff[ip-1]<0.8) flag_threshold = true;
        p_threhold = temp_p[ip];
      }

      if (temp_p[ip]>=1)
      {
        temp_display_p[pbin_display] = temp_p[ip];
        temp_display_eff[pbin_display] = temp_eff[ip];
        pbin_display++;
      }
    }
    if (flag_threshold) pt_threhold = p_threhold*sin(eta_to_theta(0.5*(eta_lo[ieta]+eta_hi[ieta])));

    g_eff_in_p[ieta] = new TGraph(pbin,temp_p,temp_eff);
    g_eff_in_p[ieta]->SetMarkerStyle(20);
    g_eff_in_p[ieta]->SetMarkerColor(kRed);
    g_eff_in_p[ieta]->SetLineColor(kRed);
    g_eff_in_p[ieta]->SetLineWidth(2);
    g_eff_in_p[ieta]->SetMarkerSize(0.5);

    mclogx(cno++);
    {
      float plot_xrange_lo = 0.1;
      float plot_xrange_hi = 60;
      float plot_yrange_lo = 0.4;
      float plot_yrange_hi = 1.1;

      float plot_xrange_hi = 60;
      if (0.5*(eta_lo[ieta]+eta_hi[ieta])<=2) plot_xrange_hi = 15;
      else if (0.5*(eta_lo[ieta]+eta_hi[ieta])<=3) plot_xrange_hi = 30;
      else plot_xrange_hi = 60;

      TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
      htemp.Draw();
      htemp.GetXaxis()->SetTitle("p [GeV/c]");
      htemp.GetYaxis()->SetTitle("efficiency");
      myhset(&htemp,1.2,1.6,0.05,0.045);

      g_eff_in_p[ieta]->Draw("pcsame");

      TLine l100(plot_xrange_lo,1,plot_xrange_hi,1);
      l100.SetLineWidth(2);
      l100.SetLineStyle(7);
      l100.Draw("same");

      TLine l80(plot_xrange_lo,0.8,plot_xrange_hi,0.8);
      l80.SetLineWidth(2);
      l80.SetLineStyle(7);
      l80.Draw("same");

      TLatex* tl = new TLatex();
      tl->SetTextAlign(11);
      tl->SetTextSize(0.030);
      tl->DrawLatexNDC(0.2,0.92,Form("Single particle events (%s), #eta [%.1f, %.1f)",seeding_name[option], eta_lo[ieta],eta_hi[ieta]));
      if (flag_threshold) 
      {
        tl->DrawLatexNDC(0.43,0.35,"Threshold (efficiency > 80%)");
        tl->DrawLatexNDC(0.43,0.30,Form("p = %.2f GeV/c (p_{T} = %.2fGeV/c)",p_threhold,pt_threhold));
      }
      else
      {
        tl->DrawLatexNDC(0.43,0.35,"Threshold (efficiency > 80%)");
        tl->DrawLatexNDC(0.43,0.30,"p < 0.1 GeV");
      }

      gROOT->ProcessLine( Form("cc%d->Print(\"figs/g_eff_vs_p_in_eta%d.pdf\")", cno-1,ieta) );
    }
  }
}

void plot_eff_2D(const int option = 0)
{
  TGaxis::SetMaxDigits(3);
  mclogxz(cno++);
  {
    float plot_xrange_lo = 0.1;
    float plot_xrange_hi = 60;
    float plot_yrange_lo = -4;
    float plot_yrange_hi = 4;
    float plot_zrange_lo = 1E2;
    float plot_zrange_hi = 5E3;
    h2d_ngen_in_p_vs_eta->GetZaxis()->SetRangeUser(plot_zrange_lo,plot_zrange_hi);

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("p [GeV/c]");
    htemp.GetYaxis()->SetTitle("#eta");
    myhset(&htemp,1.2,1.6,0.05,0.045);

    h2d_ngen_in_p_vs_eta->Draw("colzsame");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.030);
    tl->DrawLatexNDC(0.2,0.92,Form("Single particle events (%s), # of generated particles",seeding_name[option]));

    gROOT->ProcessLine( Form("cc%d->Print(\"figs/ngen_p_vs_eta.pdf\")", cno-1) );
  }

  mclogxz(cno++);
  {
    float plot_xrange_lo = 0.1;
    float plot_xrange_hi = 60;
    float plot_yrange_lo = -4;
    float plot_yrange_hi = 4;
    float plot_zrange_lo = 1E2;
    float plot_zrange_hi = 5E3;
    h2d_nreco_in_p_vs_eta->GetZaxis()->SetRangeUser(plot_zrange_lo,plot_zrange_hi);

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("p [GeV/c]");
    htemp.GetYaxis()->SetTitle("#eta");
    myhset(&htemp,1.2,1.6,0.05,0.045);

    h2d_nreco_in_p_vs_eta->Draw("colzsame");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.030);
    tl->DrawLatexNDC(0.2,0.92,Form("Single particle events (%s), # of reconstructed particles",seeding_name[option]));

    gROOT->ProcessLine( Form("cc%d->Print(\"figs/nreco_p_vs_eta.pdf\")", cno-1) );
  }

  mclogx(cno++);
  {
    float plot_xrange_lo = 0.1;
    float plot_xrange_hi = 60;
    float plot_yrange_lo = -4;
    float plot_yrange_hi = 4;
    float plot_zrange_lo = 0;
    float plot_zrange_hi = 1;
    h2d_eff_in_p_vs_eta->GetZaxis()->SetRangeUser(plot_zrange_lo,plot_zrange_hi);

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("p [GeV/c]");
    htemp.GetYaxis()->SetTitle("#eta");
    myhset(&htemp,1.2,1.6,0.05,0.045);

    TLegend leg(0.55,0.75,0.89,0.86);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.032);
    leg.SetFillStyle(0);

    h2d_eff_in_p_vs_eta->Draw("colzsame");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.030);
    tl->DrawLatexNDC(0.2,0.92,Form("Single particle events (%s), efficiency",seeding_name[option]));

    gROOT->ProcessLine( Form("cc%d->Print(\"figs/eff_p_vs_eta.pdf\")", cno-1) );
  }
}

void plot_dp_1D(const int option = 0)
{
  for (int ieta = 0; ieta < etabin; ++ieta)
  {
    double temp_p[pbin] = {0};
    double temp_dp[pbin] = {0};

    double temp_display_p[pbin] = {0};
    double temp_display_dp[pbin] = {0};
    int pbin_display = 0;
    for (int ip = 0; ip < pbin; ++ip)
    {
      temp_p[ip] = h2d_dp_in_p_vs_eta->GetXaxis()->GetBinCenter(ip+1);
      temp_dp[ip] = h2d_dp_in_p_vs_eta->GetBinContent(ip+1,ieta+1);

      if (temp_p[ip]>=1)
      {
        temp_display_p[pbin_display] = temp_p[ip];
        temp_display_dp[pbin_display] = temp_dp[ip];
        pbin_display++;
      }
    }

    g_dp_in_p[ieta] = new TGraph(pbin,temp_p,temp_dp);
    g_dp_in_p[ieta]->SetMarkerStyle(20);
    g_dp_in_p[ieta]->SetMarkerColor(kRed);
    g_dp_in_p[ieta]->SetMarkerSize(0.5);

    // int smear_graph_bin = Res_Handler->FindBin(0.5*(eta_lo[ieta]+eta_hi[ieta]));
    // if (smear_graph_bin>0 || smear_graph_bin>Res_Handler->GetNbinsX())
    // {
    //   gmom_res[smear_graph_bin-1]->SetMarkerStyle(20);
    //   gmom_res[smear_graph_bin-1]->SetMarkerColor(kBlack);
    //   gmom_res[smear_graph_bin-1]->SetMarkerSize(0.5);
    //   gmom_res[smear_graph_bin-1]->SetLineColor(kBlack);
    //   gmom_res[smear_graph_bin-1]->SetLineWidth(2);
    //   gmom_res[smear_graph_bin-1]->SetLineStyle(7);

    //   int n = gmom_res[smear_graph_bin-1]->GetN();
    //   double *y = gmom_res[smear_graph_bin-1]->GetY();
    //   for (int i=0;i<n;i++) { y[i] *= 100; }
    // }

    mcs(cno++);
    {
      float plot_xrange_lo = 0;//0.05;
      float plot_xrange_hi = 60;
      if (0.5*(eta_lo[ieta]+eta_hi[ieta])<=2) plot_xrange_hi = 15;
      else if (0.5*(eta_lo[ieta]+eta_hi[ieta])<=3) plot_xrange_hi = 30;
      else plot_xrange_hi = 60;

      float plot_yrange_lo = 0;
      float plot_yrange_hi = 2.5;//1.5*TMath::MaxElement(pbin_display,temp_display_dp);
      if (fabs(0.5*(eta_lo[ieta]+eta_hi[ieta]))>3) plot_yrange_hi = 5;

      TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
      htemp.Draw();
      htemp.GetXaxis()->SetTitle("p [GeV/c]");
      htemp.GetYaxis()->SetTitle("dp/p [%]");
      myhset(&htemp,1.2,1.6,0.05,0.045);

      TLegend leg(0.39,0.71,0.87,0.86);
      leg.SetBorderSize(0);
      leg.SetTextSize(0.035);
      leg.SetFillStyle(0);

      g_dp_in_p[ieta]->Draw("psame");
      leg.AddEntry(g_dp_in_p[ieta],"DD4HEP","p");
      // if (smear_graph_bin>0 || smear_graph_bin>Res_Handler->GetNbinsX())
      // {
      //   gmom_res[smear_graph_bin-1]->Draw("csame");
      //   leg.AddEntry(gmom_res[smear_graph_bin-1],"Fast sim (3T)","l");
      // }

      leg.Draw("same");

      TLatex* tl = new TLatex();
      tl->SetTextAlign(11);
      tl->SetTextSize(0.030);
      tl->DrawLatexNDC(0.2,0.92,Form("Single particle events (%s), #eta [%.1f, %.1f)",seeding_name[option], eta_lo[ieta],eta_hi[ieta]));

      gROOT->ProcessLine( Form("cc%d->Print(\"figs/g_dp_vs_p_in_eta%d.pdf\")", cno-1,ieta) );
    }
  }

  // for (int ip = 0; ip < pbin; ++ip)
  // {
  //   double temp_eta[etabin] = {0};
  //   double temp_dp[etabin] = {0};

  //   double temp_display_eta[etabin] = {0};
  //   double temp_display_dp[etabin] = {0};
  //   int etabin_display = 0;
  //   for (int ieta = 0; ieta < etabin; ++ieta)
  //   {
  //     temp_eta[ieta] = h2d_dp_in_p_vs_eta->GetYaxis()->GetBinCenter(ieta+1);
  //     temp_dp[ieta] = h2d_dp_in_p_vs_eta->GetBinContent(ip+1,ieta+1);

  //     if (fabs(temp_eta[ieta])<3.5)
  //     {
  //       temp_display_eta[etabin_display] = temp_eta[ieta];
  //       temp_display_dp[etabin_display] = temp_dp[ieta];
  //       etabin_display++;
  //     }
  //   }

  //   g_dp_in_eta[ip] = new TGraph(etabin_display,temp_display_eta,temp_display_dp);
  //   g_dp_in_eta[ip]->SetMarkerStyle(20);
  //   g_dp_in_eta[ip]->SetMarkerColor(kRed);
  //   g_dp_in_eta[ip]->SetMarkerSize(0.5);

  //   mcs(cno++);
  //   {
  //     float plot_xrange_lo = -4;
  //     float plot_xrange_hi = 4;

  //     float plot_yrange_lo = 0;
  //     float plot_yrange_hi = 5;//1.5*TMath::MaxElement(pbin_display,temp_display_dp);

  //     TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
  //     htemp.Draw();
  //     htemp.GetXaxis()->SetTitle("#eta");
  //     htemp.GetYaxis()->SetTitle("dp/p [%]");
  //     myhset(&htemp,1.2,1.6,0.05,0.045);

  //     TLegend leg(0.39,0.71,0.87,0.86);
  //     leg.SetBorderSize(0);
  //     leg.SetTextSize(0.035);
  //     leg.SetFillStyle(0);

  //     g_dp_in_eta[ip]->Draw("psame");
  //     leg.AddEntry(g_dp_in_eta[ip],"DD4HEP","p");

  //     leg.Draw("same");

  //     TLatex* tl = new TLatex();
  //     tl->SetTextAlign(11);
  //     tl->SetTextSize(0.030);
  //     tl->DrawLatexNDC(0.2,0.92,Form("Single particle events (%s), p [%.1f, %.1f)",p_lo[ip],p_hi[ip]));

  //     gROOT->ProcessLine( Form("cc%d->Print(\"figs/g_dp_vs_eta_in_p%d.pdf\")", cno-1,ip) );
  //   }
  // }
}

void plot_dp_2D(const int option = 0)
{
  mclogx(cno++);
  {
    float plot_xrange_lo = 0.1;
    float plot_xrange_hi = 60;
    float plot_yrange_lo = -4;
    float plot_yrange_hi = 4;
    float plot_zrange_lo = 0;
    float plot_zrange_hi = 5E0;//3E0;
    h2d_dp_in_p_vs_eta->GetZaxis()->SetRangeUser(plot_zrange_lo,plot_zrange_hi);

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("p [GeV/c]");
    htemp.GetYaxis()->SetTitle("#eta");
    myhset(&htemp,1.2,1.6,0.05,0.045);

    TLegend leg(0.55,0.75,0.89,0.86);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.032);
    leg.SetFillStyle(0);

    h2d_dp_in_p_vs_eta->Draw("colzsame");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.030);
    tl->DrawLatexNDC(0.2,0.92,Form("Single particle events (%s), dp/p resolution in %%",seeding_name[option]));

    gROOT->ProcessLine( Form("cc%d->Print(\"figs/dp_p_vs_eta.pdf\")", cno-1) );
  }
}

void plot_dphi_1D(const int option = 0)
{
  for (int ieta = 0; ieta < etabin; ++ieta)
  {
    double temp_p[pbin] = {0};
    double temp_dphi[pbin] = {0};

    double temp_display_p[pbin] = {0};
    double temp_display_dphi[pbin] = {0};
    int pbin_display = 0;
    for (int ip = 0; ip < pbin; ++ip)
    {
      temp_p[ip] = h2d_dphi_in_p_vs_eta->GetXaxis()->GetBinCenter(ip+1);
      temp_dphi[ip] = h2d_dphi_in_p_vs_eta->GetBinContent(ip+1,ieta+1);

      if (temp_p[ip]>=1)
      {
        temp_display_p[pbin_display] = temp_p[ip];
        temp_display_dphi[pbin_display] = temp_dphi[ip];
        pbin_display++;
      }
    }

    g_dphi_in_p[ieta] = new TGraph(pbin,temp_p,temp_dphi);
    g_dphi_in_p[ieta]->SetMarkerStyle(20);
    g_dphi_in_p[ieta]->SetMarkerColor(kRed);
    g_dphi_in_p[ieta]->SetMarkerSize(0.5);

    mcs(cno++);
    {
      float plot_xrange_lo = 0;
      float plot_xrange_hi = 60;
      float plot_yrange_lo = 0;
      float plot_yrange_hi = 1.5*TMath::MaxElement(pbin_display,temp_display_dphi);

      TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
      htemp.Draw();
      htemp.GetXaxis()->SetTitle("p [GeV/c]");
      htemp.GetYaxis()->SetTitle("#phi [mrad]");
      myhset(&htemp,1.2,1.6,0.05,0.045);

      g_dphi_in_p[ieta]->Draw("psame");

      TLatex* tl = new TLatex();
      tl->SetTextAlign(11);
      tl->SetTextSize(0.030);
      tl->DrawLatexNDC(0.2,0.92,Form("Single particle events (%s), #eta [%.1f, %.1f)",seeding_name[option], eta_lo[ieta],eta_hi[ieta]));

      gROOT->ProcessLine( Form("cc%d->Print(\"figs/g_dphi_vs_p_in_eta%d.pdf\")", cno-1,ieta) );
    }
  }
}

void plot_dphi_2D(const int option = 0)
{
  mclogxz(cno++);
  {
    float plot_xrange_lo = 0.1;
    float plot_xrange_hi = 60;
    float plot_yrange_lo = -4;
    float plot_yrange_hi = 4;
    float plot_zrange_lo = 5E-2;
    float plot_zrange_hi = 2E2;
    h2d_dphi_in_p_vs_eta->GetZaxis()->SetRangeUser(plot_zrange_lo,plot_zrange_hi);

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("p [GeV/c]");
    htemp.GetYaxis()->SetTitle("#eta");
    myhset(&htemp,1.2,1.6,0.05,0.045);

    h2d_dphi_in_p_vs_eta->Draw("colzsame");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.030);
    tl->DrawLatexNDC(0.2,0.92,Form("Single particle events (%s), #phi resolution in mrad",seeding_name[option]));

    gROOT->ProcessLine( Form("cc%d->Print(\"figs/dphi_p_vs_eta.pdf\")", cno-1) );
  }
}

void plot_dtheta_1D(const int option = 0)
{
  for (int ieta = 0; ieta < etabin; ++ieta)
  {
    double temp_p[pbin] = {0};
    double temp_dtheta[pbin] = {0};

    double temp_display_p[pbin] = {0};
    double temp_display_dtheta[pbin] = {0};
    int pbin_display = 0;
    for (int ip = 0; ip < pbin; ++ip)
    {
      temp_p[ip] = h2d_dtheta_in_p_vs_eta->GetXaxis()->GetBinCenter(ip+1);
      temp_dtheta[ip] = h2d_dtheta_in_p_vs_eta->GetBinContent(ip+1,ieta+1);

      if (temp_p[ip]>=1)
      {
        temp_display_p[pbin_display] = temp_p[ip];
        temp_display_dtheta[pbin_display] = temp_dtheta[ip];
        pbin_display++;
      }
    }

    g_dtheta_in_p[ieta] = new TGraph(pbin,temp_p,temp_dtheta);
    g_dtheta_in_p[ieta]->SetMarkerStyle(20);
    g_dtheta_in_p[ieta]->SetMarkerColor(kRed);
    g_dtheta_in_p[ieta]->SetMarkerSize(0.5);

    mcs(cno++);
    {
      float plot_xrange_lo = 0;
      float plot_xrange_hi = 60;
      float plot_yrange_lo = 0;
      float plot_yrange_hi = 1.5*TMath::MaxElement(pbin_display,temp_display_dtheta);

      TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
      htemp.Draw();
      htemp.GetXaxis()->SetTitle("p [GeV/c]");
      htemp.GetYaxis()->SetTitle("#theta [mrad]");
      myhset(&htemp,1.2,1.6,0.05,0.045);

      g_dtheta_in_p[ieta]->Draw("psame");

      TLatex* tl = new TLatex();
      tl->SetTextAlign(11);
      tl->SetTextSize(0.030);
      tl->DrawLatexNDC(0.2,0.92,Form("Single particle events (%s), #eta [%.1f, %.1f)",seeding_name[option], eta_lo[ieta],eta_hi[ieta]));

      gROOT->ProcessLine( Form("cc%d->Print(\"figs/g_dtheta_vs_p_in_eta%d.pdf\")", cno-1,ieta) );
    }
  }
}

void plot_dtheta_2D(const int option = 0)
{
  mclogxz(cno++);
  {
    float plot_xrange_lo = 0.1;
    float plot_xrange_hi = 60;
    float plot_yrange_lo = -4;
    float plot_yrange_hi = 4;
    // float plot_zrange_lo = 3E-2;
    // float plot_zrange_hi = 5E1;
    // h2d_dtheta_in_p_vs_eta->GetZaxis()->SetRangeUser(plot_zrange_lo,plot_zrange_hi);

    TH2F htemp("htemp","",10,plot_xrange_lo,plot_xrange_hi,10,plot_yrange_lo,plot_yrange_hi);
    htemp.Draw();
    htemp.GetXaxis()->SetTitle("p [GeV/c]");
    htemp.GetYaxis()->SetTitle("#eta");
    myhset(&htemp,1.2,1.6,0.05,0.045);

    h2d_dtheta_in_p_vs_eta->Draw("colzsame");

    TLatex* tl = new TLatex();
    tl->SetTextAlign(11);
    tl->SetTextSize(0.030);
    tl->DrawLatexNDC(0.2,0.92,Form("Single particle events (%s), #theta resolution in mrad",seeding_name[option]));

    gROOT->ProcessLine( Form("cc%d->Print(\"figs/dtheta_p_vs_eta.pdf\")", cno-1) );
  }
}

void setup_LDT(TFile* fin)
{
  cout << "Current max eta bin is " << NETA_MAX << ", adjust it inside fastsim.h when neccesary" << endl;

  if (!fin)
  {
    cout << "Cannot find the ATHENA setup file, abort..." << endl;
    exit(0);
  }

  Res_Handler = (TH1F*)fin->Get("Res_Handler");
  if (!Res_Handler)
  {
    cout << "No Res_Handler found, abort..." << endl;
    exit(0);
  }

  for(int ibin = 0; ibin < Res_Handler->GetNbinsX(); ibin++)
  {
    gmom_res[ibin] = (TGraph*)fin->Get(Form("gmom_res_%i",ibin));
    gdca_rphi_res[ibin] = (TGraph*)fin->Get(Form("gdca_rphi_res_%i",ibin));
    gdca_z_res[ibin] = (TGraph*)fin->Get(Form("gdca_z_res_%i",ibin));
  }
}

void plot_seeding_graph(const char* inFile = "graph.root", const int option = 0)
{
  mcs(-1);

  // input 
  TFile* fin = new TFile(inFile,"read");
  fin->ls();

  h2d_dp_in_p_vs_eta = (TH2D*)fin->Get("h2d_dp_in_p_vs_eta"); // z axis %
  h2d_dphi_in_p_vs_eta = (TH2D*)fin->Get("h2d_dphi_in_p_vs_eta"); // z axis mrad
  h2d_dtheta_in_p_vs_eta = (TH2D*)fin->Get("h2d_dtheta_in_p_vs_eta"); // z axis mrad

  h2d_ngen_in_p_vs_eta = (TH2D*)fin->Get("h2d_ngen_in_p_vs_eta");
  h2d_nreco_in_p_vs_eta = (TH2D*)fin->Get("h2d_nreco_in_p_vs_eta");
  h2d_eff_in_p_vs_eta = (TH2D*)h2d_nreco_in_p_vs_eta->Clone();
  h2d_eff_in_p_vs_eta->SetName("h2d_eff_in_p_vs_eta");
  h2d_eff_in_p_vs_eta->Divide(h2d_ngen_in_p_vs_eta);

  TFile* fin_LDT = new TFile("/Users/wenqingfan/EIC/fast_sim/ATHENA_Resolutions_r.root","read");
  setup_LDT(fin_LDT);

  plot_dp_1D(option);
  plot_dp_2D(option);
  plot_dphi_1D(option);
  plot_dphi_2D(option);
  plot_dtheta_1D(option);
  plot_dtheta_2D(option);
  plot_eff_1D(option);
  plot_eff_2D(option);
}
