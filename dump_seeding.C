#include "TCanvas.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TProfile.h"

#include <iostream>
#include <bitset>

#include "eff_bins.h"

TLeaf* trk_qoverp;
TLeaf* trk_theta;
TLeaf* trk_phi;
TLeaf* trk_hits;
TLeaf* trk_chi2;

TLeaf* mc_px;
TLeaf* mc_py;
TLeaf* mc_pz;
TLeaf* mc_pid;
TLeaf* mc_charge;
TLeaf* mc_status;
TLeaf* mc_gen_status;

class RealSeedReco
{
  private:
    unsigned int part_id;
    unsigned int scattered_electron_mcid;

    vector<int> matched_mc_index;
    vector<unsigned int> trk_index_0_match;
    vector<unsigned int> trk_index_1_match;
    vector<unsigned int> trk_index_N_match;

    // histograms
    TH2D* h2d_mc_gen_p_vs_eta;
    TH2D* h2d_mc_assc_p_vs_eta;
    TH1D* h1d_gen_npart[pbin][etabin];
    TH2D* h2d_1to1_dp_vs_nhit[pbin][etabin];
    TH2D* h2d_1to1_dphi_vs_nhit[pbin][etabin];
    TH2D* h2d_1to1_dtheta_vs_nhit[pbin][etabin];
    TH2D* h2d_1to1_chi2_vs_nhit[pbin][etabin];
    TH2D* h2d_Nto1_dp_vs_nhit[pbin][etabin];
    TH2D* h2d_Nto1_dphi_vs_nhit[pbin][etabin];
    TH2D* h2d_Nto1_dtheta_vs_nhit[pbin][etabin];
    TH2D* h2d_Nto1_chi2_vs_nhit[pbin][etabin];
    TH2D* h2d_0to1_dp_vs_nhit[pbin][etabin];
    TH2D* h2d_0to1_dphi_vs_nhit[pbin][etabin];
    TH2D* h2d_0to1_dtheta_vs_nhit[pbin][etabin];
    TH2D* h2d_0to1_chi2_vs_nhit[pbin][etabin];

  public:
    RealSeedReco(int _part_id)
    {
      cout << "Constructing analyzing module to study tracking performance of particle with ID " << _part_id << endl;
      part_id = abs(_part_id);

      scattered_electron_mcid = 0;

      InitBins();

      cout << "Initalize histograms..." << endl;
      h2d_mc_gen_p_vs_eta = new TH2D("h2d_mc_gen_p_vs_eta","p vs eta",600,0,60,40,-3.5,3.5);
      h2d_mc_gen_p_vs_eta->Sumw2();
      h2d_mc_assc_p_vs_eta = new TH2D("h2d_mc_assc_p_vs_eta","p vs eta",600,0,60,40,-3.5,3.5);
      h2d_mc_assc_p_vs_eta->Sumw2();

      for (int ip = 0; ip < pbin; ++ip)
      {
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          float p_mid = 0.5*(p_lo[ip]+p_hi[ip]);
          float eta_mid = 0.5*(eta_lo[ieta]+eta_hi[ieta]);
          float theta_mid = 2*atan(exp(-eta_mid));
          float pt_mid = p_mid*sin(theta_mid);

          h1d_gen_npart[ip][ieta] = new TH1D(Form("h1d_gen_npart_%d_%d",ip,ieta),"generated particle",1,0.5,1.5);
          h1d_gen_npart[ip][ieta]->Sumw2();

          h2d_1to1_dp_vs_nhit[ip][ieta] = new TH2D(Form("h2d_1to1_dp_vs_nhit_%d_%d",ip,ieta),"dp/p vs nhit",200,-0.2,0.2,11,-0.5,10.5); // FIX ME: if want to use to efficiecny def, need to make sure it's wide enough
          h2d_1to1_dp_vs_nhit[ip][ieta]->Sumw2();
          h2d_Nto1_dp_vs_nhit[ip][ieta] = new TH2D(Form("h2d_Nto1_dp_vs_nhit_%d_%d",ip,ieta),"dp/p vs nhit",200,-0.2,0.2,11,-0.5,10.5);
          h2d_Nto1_dp_vs_nhit[ip][ieta]->Sumw2();
          h2d_0to1_dp_vs_nhit[ip][ieta] = new TH2D(Form("h2d_0to1_dp_vs_nhit_%d_%d",ip,ieta),"dp/p vs nhit",200,-0.2,0.2,11,-0.5,10.5); 
          h2d_0to1_dp_vs_nhit[ip][ieta]->Sumw2();

          float phi_range_width = 5*(1.0/pt_mid+0.5)*1E-3;
          h2d_1to1_dphi_vs_nhit[ip][ieta] = new TH2D(Form("h2d_1to1_dphi_vs_nhit_%d_%d",ip,ieta),"dphi vs nhit",500,-phi_range_width,phi_range_width,11,-0.5,10.5);
          h2d_1to1_dphi_vs_nhit[ip][ieta]->Sumw2();
          h2d_Nto1_dphi_vs_nhit[ip][ieta] = new TH2D(Form("h2d_Nto1_dphi_vs_nhit_%d_%d",ip,ieta),"dphi vs nhit",500,-phi_range_width,phi_range_width,11,-0.5,10.5);
          h2d_Nto1_dphi_vs_nhit[ip][ieta]->Sumw2();
          h2d_0to1_dphi_vs_nhit[ip][ieta] = new TH2D(Form("h2d_0to1_dphi_vs_nhit_%d_%d",ip,ieta),"dphi vs nhit",500,-phi_range_width,phi_range_width,11,-0.5,10.5);
          h2d_0to1_dphi_vs_nhit[ip][ieta]->Sumw2();

          float theta_range_width = 5*(0.4/pt_mid+0.2)*1E-3;
          h2d_1to1_dtheta_vs_nhit[ip][ieta] = new TH2D(Form("h2d_1to1_dtheta_vs_nhit_%d_%d",ip,ieta),"dtheta vs nhit",300,-theta_range_width,theta_range_width,11,-0.5,10.5);
          h2d_1to1_dtheta_vs_nhit[ip][ieta]->Sumw2();
          h2d_Nto1_dtheta_vs_nhit[ip][ieta] = new TH2D(Form("h2d_Nto1_dtheta_vs_nhit_%d_%d",ip,ieta),"dtheta vs nhit",300,-theta_range_width,theta_range_width,11,-0.5,10.5);
          h2d_Nto1_dtheta_vs_nhit[ip][ieta]->Sumw2();
          h2d_0to1_dtheta_vs_nhit[ip][ieta] = new TH2D(Form("h2d_0to1_dtheta_vs_nhit_%d_%d",ip,ieta),"dtheta vs nhit",300,-theta_range_width,theta_range_width,11,-0.5,10.5);
          h2d_0to1_dtheta_vs_nhit[ip][ieta]->Sumw2();

          h2d_0to1_chi2_vs_nhit[ip][ieta] = new TH2D(Form("h2d_0to1_chi2_vs_nhit_%d_%d",ip,ieta),"chi2 vs nhit",200,0,50,11,-0.5,10.5); 
          h2d_0to1_chi2_vs_nhit[ip][ieta]->Sumw2();
          h2d_Nto1_chi2_vs_nhit[ip][ieta] = new TH2D(Form("h2d_Nto1_chi2_vs_nhit_%d_%d",ip,ieta),"chi2 vs nhit",200,0,50,11,-0.5,10.5); 
          h2d_Nto1_chi2_vs_nhit[ip][ieta]->Sumw2();
          h2d_1to1_chi2_vs_nhit[ip][ieta] = new TH2D(Form("h2d_1to1_chi2_vs_nhit_%d_%d",ip,ieta),"chi2 vs nhit",200,0,50,11,-0.5,10.5); 
          h2d_1to1_chi2_vs_nhit[ip][ieta]->Sumw2();
        }
      }
    }

    virtual ~RealSeedReco()
    {
      delete h2d_mc_gen_p_vs_eta;
      delete h2d_mc_assc_p_vs_eta;
      for (int ip = 0; ip < pbin; ++ip)
      {
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          delete h1d_gen_npart[ip][ieta];
          delete h2d_1to1_dp_vs_nhit[ip][ieta];
          delete h2d_1to1_dphi_vs_nhit[ip][ieta];
          delete h2d_1to1_dtheta_vs_nhit[ip][ieta];
          delete h2d_1to1_chi2_vs_nhit[ip][ieta];
          delete h2d_Nto1_dp_vs_nhit[ip][ieta];
          delete h2d_Nto1_dphi_vs_nhit[ip][ieta];
          delete h2d_Nto1_dtheta_vs_nhit[ip][ieta];
          delete h2d_Nto1_chi2_vs_nhit[ip][ieta];
          delete h2d_0to1_dp_vs_nhit[ip][ieta];
          delete h2d_0to1_dphi_vs_nhit[ip][ieta];
          delete h2d_0to1_dtheta_vs_nhit[ip][ieta];
          delete h2d_0to1_chi2_vs_nhit[ip][ieta];
        }
      }
    };

    void InitBins()
    {
      cout << "Initialization from header, print p and eta bins..." << endl;
      for (int ip = 0; ip < pbin; ++ip)
      {
        if (ip==0) cout << "p_lo[" << pbin << "] = {";
        cout << p_lo[ip];
        if (ip<pbin-1) cout << ", ";
        if (ip==pbin-1) cout << "};" << endl;
      }
      for (int ip = 0; ip < pbin; ++ip)
      {
        if (ip==0) cout << "p_hi[" << pbin << "] = {";
        cout << p_hi[ip];
        if (ip<pbin-1) cout << ", ";
        if (ip==pbin-1) cout << "};" << endl;
      }
      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        if (ieta==0) cout << "eta_lo[" << etabin << "] = {";
        cout << eta_lo[ieta];
        if (ieta<etabin-1) cout << ", ";
        if (ieta==etabin-1) cout << "};" << endl;
      }
      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        if (ieta==0) cout << "eta_hi[" << etabin << "] = {";
        cout << eta_hi[ieta];
        if (ieta<etabin-1) cout << ", ";
        if (ieta==etabin-1) cout << "};" << endl;
      }
    }

    void ClearMatchedParts()
    {
      matched_mc_index.clear();
      trk_index_0_match.clear();
      trk_index_1_match.clear();
      trk_index_N_match.clear();
    }

    void FillGenHists()
    {
      // cout << "mc length " << mc_px->GetLen() << endl;
      for (int imc = 0; imc < mc_px->GetLen(); ++imc)
      {
        if (mc_gen_status->GetValue(imc)!=1) continue; // only loop though stable particles
        if (abs(mc_pid->GetValue(imc))!=part_id) continue;

        TVector3 mc_vec(mc_px->GetValue(imc),mc_py->GetValue(imc),mc_pz->GetValue(imc));

        float mc_p = mc_vec.Mag();
        float mc_phi = mc_vec.Phi();
        float mc_theta = mc_vec.Theta();
        float mc_eta = -log(tan(mc_theta/2));

        h2d_mc_gen_p_vs_eta->Fill(mc_p,mc_eta);

        int ipbin = -9999;
        for (int ip = 0; ip < pbin; ++ip)
        {
          if (mc_p>=p_lo[ip] && mc_p<p_hi[ip]) ipbin = ip; 
        }
        int ietabin = -9999;
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          if (mc_eta>=eta_lo[ieta] && mc_eta<eta_hi[ieta]) ietabin = ieta; 
        }         

        if (ipbin>-1 && ietabin>-1) h1d_gen_npart[ipbin][ietabin]->Fill(1);
      }
    }

    void FillMatchedParts()
    {
      ClearMatchedParts();

      for (int itrk = 0; itrk < trk_qoverp->GetLen(); ++itrk)
      {
        TVector3 trk_vec(1,0,0);
        trk_vec.SetMag( fabs(1.0/trk_qoverp->GetValue(itrk)) );
        trk_vec.SetTheta( trk_theta->GetValue(itrk) );
        trk_vec.SetPhi( trk_phi->GetValue(itrk) ); 

        int index_match = -1;
        bool flag_match = true;
        for (int imc = 0; imc < mc_px->GetLen(); ++imc)
        {
          if (mc_gen_status->GetValue(imc)!=1) continue; // only loop though stable particles
          if (abs(mc_pid->GetValue(imc))!=part_id) continue;

          TVector3 mc_vec(mc_px->GetValue(imc),mc_py->GetValue(imc),mc_pz->GetValue(imc));

          flag_match = true;
          if (trk_qoverp->GetValue(itrk)*mc_charge->GetValue(imc)<0) flag_match = false; // opposite sign
          // if ( fabs(mc_vec.Mag()-trk_vec.Mag())/mc_vec.Mag()>0.1 ) flag_match = false; // dp/p < 10%
          // if ( fabs(mc_vec.Phi()-trk_vec.Phi())>0.05 ) flag_match = false; // dphi < 50mrad
          // if ( fabs(mc_vec.Theta()-trk_vec.Theta())>0.01 ) flag_match = false; // dtheta < 10mrad (5mrad)

          if (flag_match)
          {
            index_match = imc;
            break;
          }
        }

        matched_mc_index.push_back(index_match);        
      }

      // look for multiple match
      for (int itrk = 0; itrk < trk_qoverp->GetLen(); ++itrk)
      {
        if (matched_mc_index[itrk]<0)
        {
          trk_index_0_match.push_back(itrk);
        }
        else 
        {
          // cout << "matched mc " << itrk << ", " << matched_mc_index[itrk] << endl;
          trk_index_1_match.push_back(itrk); // just candidates, can contain N match cases

          bool flag_N_match = false;
          for (int i_1_match = 0; i_1_match < trk_index_1_match.size()-1; ++i_1_match)
          {
            int itrk_current = itrk;
            int itrk_saved = trk_index_1_match[i_1_match];

            int imc_current = matched_mc_index[itrk_current];
            int imc_saved = matched_mc_index[itrk_saved];
            
            if ( imc_current==imc_saved )
            { // if current track index has already been saved
              TVector3 mc_vec(mc_px->GetValue(imc_current),mc_py->GetValue(imc_current),mc_pz->GetValue(imc_current));

              TVector3 trk_vec_current(1,0,0);
              trk_vec_current.SetMag( fabs(1.0/trk_qoverp->GetValue(itrk_current)) );
              trk_vec_current.SetTheta( trk_theta->GetValue(itrk_current) );
              trk_vec_current.SetPhi( trk_phi->GetValue(itrk_current) ); 

              TVector3 trk_vec_saved(1,0,0);
              trk_vec_saved.SetMag( fabs(1.0/trk_qoverp->GetValue(itrk_saved)) );
              trk_vec_saved.SetTheta( trk_theta->GetValue(itrk_saved) );
              trk_vec_saved.SetPhi( trk_phi->GetValue(itrk_saved) ); 

              if ( fabs(trk_vec_current.Mag()-mc_vec.Mag())/mc_vec.Mag() < fabs(trk_vec_saved.Mag()-mc_vec.Mag())/mc_vec.Mag() )
              { // current track is a better match in term of reconstructed p, updated trk_index_1_match[i_1_match] to current trk index, put saved track to duplicate list
                trk_index_N_match.push_back(itrk_saved);
                trk_index_1_match[i_1_match] = itrk_current;
              }
              else
              { // saved track is a better match, put current track to duplicate list
                trk_index_N_match.push_back(itrk_current);
              }
              trk_index_1_match.pop_back(); // remove N match cases

              break;
            }
          }
        }
      }
    }

    void GetMatchedReso(unsigned int itrk_matched, unsigned int imc_matched, int& ipbin, int& ietabin, int& nhit, double& dpp, double& dphi, double& dtheta, double& chi2)
    {
      TVector3 trk_vec(1,0,0);
      trk_vec.SetMag( fabs(1.0/trk_qoverp->GetValue(itrk_matched)) );
      trk_vec.SetTheta( trk_theta->GetValue(itrk_matched) );
      trk_vec.SetPhi( trk_phi->GetValue(itrk_matched) );   

      float trk_p = trk_vec.Mag();
      float trk_phi = trk_vec.Phi();
      float trk_theta = trk_vec.Theta();
      float trk_eta = -log(tan(trk_theta/2));

      string full_str = std::bitset<12>(trk_hits->GetValue(itrk_matched)).to_string();
      std::bitset<3> b_vertex(full_str.substr(9,3));
      std::bitset<7> b_track(full_str.substr(2,7));
      std::bitset<10> b_full_track(full_str.substr(2,10));
      std::bitset<1> b_pid(full_str.substr(1,1));
      std::bitset<1> b_ecal(full_str.substr(0,1));
      nhit = b_full_track.count();

      TVector3 mc_vec(mc_px->GetValue(imc_matched),mc_py->GetValue(imc_matched),mc_pz->GetValue(imc_matched));

      float mc_p = mc_vec.Mag();
      float mc_phi = mc_vec.Phi();
      float mc_theta = mc_vec.Theta();
      float mc_eta = -log(tan(mc_theta/2));

      for (int ip = 0; ip < pbin; ++ip)
      {
        if (mc_p>=p_lo[ip] && mc_p<p_hi[ip]) ipbin = ip; 
      }
      for (int ieta = 0; ieta < etabin; ++ieta)
      {
        if (mc_eta>=eta_lo[ieta] && mc_eta<eta_hi[ieta]) ietabin = ieta; 
      }

      if (ipbin>-1 && ietabin>-1)
      {
        dpp = (trk_p-mc_p)/mc_p;
        dphi = trk_phi-mc_phi;
        dtheta = trk_theta-mc_theta;
        chi2 = trk_chi2->GetValue(itrk_matched);
      }
    }

    void FillTrackingHists()
    {
      if (trk_index_1_match.size()>1) cout << "WARNING: more than one match" << endl;

      for (int itrk = 0; itrk < trk_index_1_match.size(); ++itrk)
      {
        unsigned int itrk_matched = trk_index_1_match[itrk];
        unsigned int imc_matched = matched_mc_index[itrk_matched];
        if (trk_index_1_match.size()>1) cout << "1 to 1 match " << itrk_matched << ", " << imc_matched << endl;

        int ipbin = -9999, ietabin = -9999;
        int nhit_full_track = -9999;
        double dpp = -9999, dphi = -9999, dtheta = -9999, chi2 = -9999;
        GetMatchedReso(itrk_matched, imc_matched, ipbin, ietabin, nhit_full_track, dpp, dphi, dtheta, chi2);

        TVector3 mc_vec(mc_px->GetValue(imc_matched),mc_py->GetValue(imc_matched),mc_pz->GetValue(imc_matched));
        float mc_p = mc_vec.Mag();
        float mc_theta = mc_vec.Theta();
        float mc_eta = -log(tan(mc_theta/2));
        h2d_mc_assc_p_vs_eta->Fill(mc_p, mc_eta);

        if (ipbin>-1 && ietabin>-1)
        {
          h2d_1to1_dp_vs_nhit[ipbin][ietabin]->Fill( dpp, nhit_full_track );
          h2d_1to1_dphi_vs_nhit[ipbin][ietabin]->Fill( dphi, nhit_full_track );
          h2d_1to1_dtheta_vs_nhit[ipbin][ietabin]->Fill( dtheta, nhit_full_track );
          h2d_1to1_chi2_vs_nhit[ipbin][ietabin]->Fill( chi2, nhit_full_track );
        }
      }

      for (int itrk = 0; itrk < trk_index_N_match.size(); ++itrk)
      {
        unsigned int itrk_matched = trk_index_N_match[itrk];
        unsigned int imc_matched = matched_mc_index[itrk_matched];
        int ipbin = -9999, ietabin = -9999;
        int nhit_full_track = -9999;
        double dpp = -9999, dphi = -9999, dtheta = -9999, chi2 = -9999;
        GetMatchedReso(itrk_matched, imc_matched, ipbin, ietabin, nhit_full_track, dpp, dphi, dtheta, chi2);

        if (ipbin>-1 && ietabin>-1)
        {
          h2d_Nto1_dp_vs_nhit[ipbin][ietabin]->Fill( dpp, nhit_full_track );
          h2d_Nto1_dphi_vs_nhit[ipbin][ietabin]->Fill( dphi, nhit_full_track );
          h2d_Nto1_dtheta_vs_nhit[ipbin][ietabin]->Fill( dtheta, nhit_full_track );
          h2d_Nto1_chi2_vs_nhit[ipbin][ietabin]->Fill( chi2, nhit_full_track );
        }
      }

      for (int itrk = 0; itrk < trk_index_0_match.size(); ++itrk)
      {
        unsigned int itrk_matched = trk_index_0_match[itrk];
        unsigned int imc_matched = 2; // FIX ME: only works for single particle event, for DIS event not sure how to do it
        int ipbin = -9999, ietabin = -9999;
        int nhit_full_track = -9999;
        double dpp = -9999, dphi = -9999, dtheta = -9999, chi2 = -9999;
        GetMatchedReso(itrk_matched, imc_matched, ipbin, ietabin, nhit_full_track, dpp, dphi, dtheta, chi2);

        if (ipbin>-1 && ietabin>-1)
        {
          h2d_0to1_dp_vs_nhit[ipbin][ietabin]->Fill( dpp, nhit_full_track );
          h2d_0to1_dphi_vs_nhit[ipbin][ietabin]->Fill( dphi, nhit_full_track );
          h2d_0to1_dtheta_vs_nhit[ipbin][ietabin]->Fill( dtheta, nhit_full_track );
          h2d_0to1_chi2_vs_nhit[ipbin][ietabin]->Fill( chi2, nhit_full_track );
        }
      }
    }

    void Write()
    {
      h2d_mc_gen_p_vs_eta->Write();
      h2d_mc_assc_p_vs_eta->Write();
      for (int ip = 0; ip < pbin; ++ip)
      {
        for (int ieta = 0; ieta < etabin; ++ieta)
        {
          h1d_gen_npart[ip][ieta]->Write();
          h2d_1to1_dp_vs_nhit[ip][ieta]->Write();
          h2d_1to1_dphi_vs_nhit[ip][ieta]->Write();
          h2d_1to1_dtheta_vs_nhit[ip][ieta]->Write();
          h2d_1to1_chi2_vs_nhit[ip][ieta]->Write();
          h2d_Nto1_dp_vs_nhit[ip][ieta]->Write();
          h2d_Nto1_dphi_vs_nhit[ip][ieta]->Write();
          h2d_Nto1_dtheta_vs_nhit[ip][ieta]->Write();
          h2d_Nto1_chi2_vs_nhit[ip][ieta]->Write();
          h2d_0to1_dp_vs_nhit[ip][ieta]->Write();
          h2d_0to1_dphi_vs_nhit[ip][ieta]->Write();
          h2d_0to1_dtheta_vs_nhit[ip][ieta]->Write();
          h2d_0to1_chi2_vs_nhit[ip][ieta]->Write();
        }
      }
    }
};

void dump_seeding(const char* inFile = "TTree.root", const char* outFile = "hist.root", int nprocess = 0)
{
  TChain *tree = new TChain("events");
  tree->Add(inFile);

  trk_qoverp = (TLeaf*)tree->GetLeaf("outputTrackParameters","outputTrackParameters.qOverP");
  trk_phi = (TLeaf*)tree->GetLeaf("outputTrackParameters","outputTrackParameters.phi");
  trk_theta = (TLeaf*)tree->GetLeaf("outputTrackParameters","outputTrackParameters.theta");
  trk_hits = (TLeaf*)tree->GetLeaf("outputTrackParameters","outputTrackParameters.charge");
  trk_chi2 = (TLeaf*)tree->GetLeaf("outputTrackParameters","outputTrackParameters.timeError");

  mc_px = (TLeaf*)tree->GetLeaf("MCParticles","MCParticles.momentum.x");
  mc_py = (TLeaf*)tree->GetLeaf("MCParticles","MCParticles.momentum.y");
  mc_pz = (TLeaf*)tree->GetLeaf("MCParticles","MCParticles.momentum.z");
  mc_pid = (TLeaf*)tree->GetLeaf("MCParticles","MCParticles.PDG");
  mc_charge = (TLeaf*)tree->GetLeaf("MCParticles","MCParticles.charge");
  mc_status = (TLeaf*)tree->GetLeaf("MCParticles","MCParticles.simulationStatus");
  mc_gen_status = (TLeaf*)tree->GetLeaf("MCParticles","MCParticles.generatorStatus");

  RealSeedReco ana_pions(211);
  int nevent = tree->GetEntries();
  if (nprocess>0) nevent = nprocess;

  for(int ievent=0; ievent < nevent; ievent++) 
  {
    tree->GetEntry(ievent);
    if (ievent%1000==0) cout << "Event = " << ievent << "/" << tree->GetEntries() << endl;

    ana_pions.FillGenHists();
    ana_pions.FillMatchedParts();
    ana_pions.FillTrackingHists();
  }

  cout << "Writing results out to root files" << endl;
  TFile* fout = new TFile(outFile,"recreate");
  ana_pions.Write();
  fout->Write();
  fout->Close();
}
