#include <iostream>
#include <cmath>

#include "TChain.h"
#include "TFile.h"
#include "TH2D.h"
#include "TLorentzVector.h"

using namespace std;

void EPOS_analysisCode2()
{
	const int nopt = 2; // stat 0,8 / stat 6,8 / stat 9-80,81
	const int nch = 8;	// rho, K*, phi, Delta++, Sigma-, Sigma+, Lambda*, Xi0

	// Single System (pp) Multiplicity Limits
	const int Nchmin_FT0 = 0;
	const int Nchmax_FT0 = 250;
	const int Nchbins_FT0 = (Nchmax_FT0 - Nchmin_FT0);

	const int Nchmin_mid = 0;
	const int Nchmax_mid = 100;
	const int Nchbins_mid = (Nchmax_mid - Nchmin_mid);

	const float ptmin = 0.;
	const float ptmax = 50.;
	const float ptintervals = 0.1;
	const int ptbins = (ptmax - ptmin) / ptintervals;

	// Define Histograms
	TH2D *hNch[nopt];
	TH2D *hDec[nopt][nch];
	TH2D *hRes_ON[nopt][nch];
	TH2D *hRes_OFF[nopt][nch];

	for (int iopt = 0; iopt < nopt; iopt++)
	{
		hNch[iopt] = new TH2D(Form("hNch_opt%d", iopt), "", Nchbins_FT0, Nchmin_FT0, Nchmax_FT0, Nchbins_mid, Nchmin_mid, Nchmax_mid);

		for (int ich = 0; ich < nch; ich++)
		{
			hDec[iopt][ich] = new TH2D(Form("hDec_opt%d_ch%d", iopt, ich), "", Nchbins_FT0, Nchmin_FT0, Nchmax_FT0, ptbins, ptmin, ptmax);
			hRes_ON[iopt][ich] = new TH2D(Form("hRes_stat9_opt%d_ch%d_mode0", iopt, ich), "", Nchbins_FT0, Nchmin_FT0, Nchmax_FT0, ptbins, ptmin, ptmax);
			hRes_OFF[iopt][ich] = new TH2D(Form("hRes_opt%d_ch%d", iopt, ich), "", Nchbins_FT0, Nchmin_FT0, Nchmax_FT0, ptbins, ptmin, ptmax);
		} // ich
	} // iopt

	// EPOS particle IDs mapping
	const int pion_pid = 120;
	const int kaon_pid = 130;
	const int proton_pid = 1120;
	const int Lambda_pid = 2130;
	const int Xi_m_pid = 2330;

	const int rho0_pid = 111;
	const int Kstar_pid = 231;
	const int phi_pid = 331;
	const int Delta_pp_pid = 1221;
	const int Sigma_star_m_pid = 2131;
	const int Sigma_star_p_pid = 2231;
	const int Lambda_star_pid = 3124; 
	const int Xi_star_0_pid = 1331;

	const int res_pid[nch] = {rho0_pid, Kstar_pid, phi_pid, Delta_pp_pid, Sigma_star_m_pid, Sigma_star_p_pid, Lambda_star_pid, Xi_star_0_pid};
	const int dec_pid[nch] = {pion_pid, kaon_pid, kaon_pid, proton_pid, Lambda_pid, Lambda_pid, Lambda_pid, Xi_m_pid};
	const int dec_chain1[nch] = {pion_pid, pion_pid, kaon_pid, pion_pid, pion_pid, pion_pid, kaon_pid, pion_pid};
	const int dec_chain2[nch] = {pion_pid, kaon_pid, kaon_pid, proton_pid, Lambda_pid, Lambda_pid, proton_pid, Xi_m_pid};

	// Single input tree setup
	TChain chain("teposevent0");
	chain.Add("/home/sawan/Storage/EPOS_localOutputs/mergedEPOS_UrQMD2.root");

	Long64_t nentries = chain.GetEntries();
	cout << "Total events to process = " << nentries << endl;
    // nentries = 100000; // for QA test, comment out for full processing

	const int MAX_TRK = 10000;
	int np;
	float px[MAX_TRK], py[MAX_TRK], pz[MAX_TRK], e[MAX_TRK];
	int id[MAX_TRK], ist[MAX_TRK], ity[MAX_TRK], ior[MAX_TRK], jor[MAX_TRK];

	chain.SetBranchAddress("np", &np);
	chain.SetBranchAddress("px", px);
	chain.SetBranchAddress("py", py);
	chain.SetBranchAddress("pz", pz);
	chain.SetBranchAddress("e", e);
	chain.SetBranchAddress("id", id);
	chain.SetBranchAddress("ist", ist);
	chain.SetBranchAddress("ity", ity);
	chain.SetBranchAddress("ior", ior);
	chain.SetBranchAddress("jor", jor);

	// Event loop
	for (Long64_t ien = 0; ien < nentries; ien++)
	{
		chain.GetEntry(ien);
		if (ien % 10000 == 0) cout << "Processing Event " << ien << " / " << nentries << endl;

		bool INEL_0[nopt] = {false, false};
		int Nch_FT0[nopt] = {0, 0};
		int Nch_mid[nopt] = {0, 0};

		// 1. Calculate Multiplicity and INEL selection
		for (int itrk = 0; itrk < np; itrk++)
		{
			int abs_id = abs(id[itrk]);
			if (!(abs_id == pion_pid || abs_id == kaon_pid || abs_id == proton_pid))
				continue;

			TLorentzVector p(px[itrk], py[itrk], pz[itrk], e[itrk]);
			double eta = p.Eta();

			if (fabs(eta) < 1.0 && (!INEL_0[0] || !INEL_0[1]))
			{
				if (ist[itrk] == 0)
					INEL_0[0] = true;
				else if (ist[itrk] == 8)
					INEL_0[1] = true;
			}

			if ((eta > -3.3 && eta < -2.1) || (eta > 3.5 && eta < 4.9)) // FT0 Forward Region
			{
				if (ist[itrk] == 0)
					Nch_FT0[0]++;
				else if (ist[itrk] == 8)
					Nch_FT0[1]++;
			}

			if (fabs(eta) < 0.5) // Mid-rapidity Region
			{
				if (ist[itrk] == 0)
					Nch_mid[0]++;
				else if (ist[itrk] == 8)
					Nch_mid[1]++;
			}
		} // itrk

		if (!INEL_0[0] && !INEL_0[1])
			continue;
		if (INEL_0[0])
			hNch[0]->Fill(Nch_FT0[0], Nch_mid[0]);
		if (INEL_0[1])
			hNch[1]->Fill(Nch_FT0[1], Nch_mid[1]);

		// 2. Track Loop for Yields & Resonances
		for (int itrk = 0; itrk < np; itrk++)
		{
			TLorentzVector final_state(px[itrk], py[itrk], pz[itrk], e[itrk]);
			double pt = final_state.Pt();
			int abs_id = abs(id[itrk]);

			if (fabs(final_state.Rapidity()) < 0.5)
			{
				for (int ich = 0; ich < nch; ich++)
				{
					if (abs_id == dec_pid[ich])
					{
						if (ist[itrk] == 0 && INEL_0[0])
							hDec[0][ich]->Fill(Nch_FT0[0], pt);
						else if (ist[itrk] == 8 && INEL_0[1])
							hDec[1][ich]->Fill(Nch_FT0[1], pt);
					}
					if (abs_id == res_pid[ich])
					{
						if (ist[itrk] == 9 && INEL_0[0])
						{
							if (ity[itrk] == 80)
								hRes_ON[0][ich]->Fill(Nch_FT0[0], pt);
							else if (ity[itrk] == 81)
								hRes_ON[1][ich]->Fill(Nch_FT0[0], pt);
						}
						else if (ist[itrk] == 6 && INEL_0[1])
						{
							hRes_OFF[0][ich]->Fill(Nch_FT0[1], pt);
						}
					}
				} // ich
			} // rapidity

			// 3. Rescattering check using parent array index (ior)
			if (itrk + 1 >= np)
				continue;

			int mom_idx1 = ior[itrk];
			int mom_idx2 = ior[itrk + 1];

			if (mom_idx1 > 0 && mom_idx1 == mom_idx2 && (mom_idx1 - 1) < np)
			{
				int parent_array_idx = mom_idx1 - 1; // 1-based to 0-based conversion
				TLorentzVector Truth(px[parent_array_idx], py[parent_array_idx], pz[parent_array_idx], e[parent_array_idx]);

				if (fabs(Truth.Rapidity()) < 0.5)
				{
					int mom_id = abs(id[parent_array_idx]);
					for (int ich = 0; ich < nch; ich++)
					{
						if ((mom_id == res_pid[ich]) && ist[itrk] == 8 && INEL_0[1])
						{
							int d1_id = abs(id[itrk]);
							int d2_id = abs(id[itrk + 1]);

							if ((d1_id == dec_chain1[ich] && d2_id == dec_chain2[ich]) ||
								(d1_id == dec_chain2[ich] && d2_id == dec_chain1[ich]))
							{
								hRes_OFF[1][ich]->Fill(Nch_FT0[1], Truth.Pt());
							}
						}
					} // ich
				} // rapidity
			}
		} // itrk
	} // ien

	// Write output ROOT file
	TFile *outfile = new TFile("yield_outfile_pp13TeV.root", "recreate");

	for (int iopt = 0; iopt < nopt; iopt++)
	{
		hNch[iopt]->Write();

		for (int ich = 0; ich < nch; ich++)
		{
			hDec[iopt][ich]->Write();
			hRes_ON[iopt][ich]->Write();
			hRes_OFF[iopt][ich]->Write();
		} // ich
	} // iopt

	outfile->Close();
	cout << "Finished processing pp system!" << endl;
}