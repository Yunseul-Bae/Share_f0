#include <TCanvas.h>
#include <TLatex.h>
#include <TString.h>
#include <sstream>
#include <iomanip>
#include <string>
#include <iostream>

void drawInvMass() {
	TFile* fin = new TFile("$Pwork/results/HL_v2.root","read");
	// TFile* fin = new TFile("~/Downloads/AnalysisResults.root","read");

	THnSparse* hInvMassUS = (THnSparse*)fin->Get("lf-f0980analysis/hInvMass_f0980_US");
	THnSparse* hInvMassLSpp = (THnSparse*)fin->Get("lf-f0980analysis/hInvMass_f0980_LSpp");
	THnSparse* hInvMassLSmm = (THnSparse*)fin->Get("lf-f0980analysis/hInvMass_f0980_LSmm");

	const int nmult = 7;
	const int npt = 16;
	// const int nlpt = 2;
	// const int nrt = 3;

	//	0-1, 1-5, 5-10, 10-20, 20-30, 30-50, 50-100
	int m_min[nmult] = {0, 1, 5, 10, 20, 30, 50};
	int m_max[nmult] = {1, 5, 10, 20, 30, 50, 100};
	//	0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 13.0
	double p_min[npt] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0};
	double p_max[npt] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 13.0};

	// double lpt_min[nlpt] = {0, 5};

	// TH1D* hProjInvMassUS[nlpt][nrt][nmult][npt];
	// TH1D* hProjInvMassLSpp[nlpt][nrt][nmult][npt];
	// TH1D* hProjInvMassLSmm[nlpt][nrt][nmult][npt];
	// TH1D* hProjInvMassLS[nlpt][nrt][nmult][npt];
	// TH1D* hProjInvMassSub[nlpt][nrt][nmult][npt];
	// TH1D* hProjInvMassSub_InJet[nlpt][nmult][npt];
	TH1D* hProjInvMassUS[nmult][npt];
	TH1D* hProjInvMassLSpp[nmult][npt];
	TH1D* hProjInvMassLSmm[nmult][npt];
	TH1D* hProjInvMassLS[nmult][npt];
	TH1D* hProjInvMassSub[nmult][npt];

	// for (int i = 0 ; i < nlpt ; i++) { // 0, 1
		// hInvMassUS->GetAxis(4)->SetRangeUser( lpt_min[i], 1e4 );
		// hInvMassLSpp->GetAxis(4)->SetRangeUser( lpt_min[i], 1e4 );
		// hInvMassLSmm->GetAxis(4)->SetRangeUser( lpt_min[i], 1e4 );
		
		// for (int r = 0 ; r < nrt ; r++) { // 0, 1, 2
			// hInvMassUS->GetAxis(3)->SetRangeUser(r, r+0.99);
			// hInvMassLSpp->GetAxis(3)->SetRangeUser(r, r+0.99);
			// hInvMassLSmm->GetAxis(3)->SetRangeUser(r, r+0.99);

			//	centrality = j
			for (int j = 0 ; j < nmult ; j++) { // 0, 1, 2, 3
				hInvMassUS->GetAxis(2)->SetRangeUser( m_min[j], m_max[j]-0.001 );
				hInvMassLSpp->GetAxis(2)->SetRangeUser( m_min[j], m_max[j]-0.001 );
				hInvMassLSmm->GetAxis(2)->SetRangeUser( m_min[j], m_max[j]-0.001 );

				for (int k = 0 ; k < npt ; k++) { // 0 - 15
					hInvMassUS->GetAxis(1)->SetRangeUser( p_min[k], p_max[k]-0.001 );
					hInvMassLSpp->GetAxis(1)->SetRangeUser( p_min[k], p_max[k]-0.001 );
					hInvMassLSmm->GetAxis(1)->SetRangeUser( p_min[k], p_max[k]-0.001 );
					
					// hProjInvMassUS[i][r][j][k] = hInvMassUS->Projection(0,"e");
					// hProjInvMassUS[i][r][j][k]->SetName(Form("hInvMassUS_%d_%d_%d_%d",i,r,j,k));
					// hProjInvMassLSpp[i][r][j][k] = hInvMassLSpp->Projection(0,"e");
					// hProjInvMassLSpp[i][r][j][k]->SetName(Form("hProjInvMassLSpp_%d_%d_%d_%d",i,r,j,k));
					// hProjInvMassLSmm[i][r][j][k] = hInvMassLSmm->Projection(0,"e");
					// hProjInvMassLSmm[i][r][j][k]->SetName(Form("hProjInvMassLSmm_%d_%d_%d_%d",i,r,j,k));
					// hProjInvMassLS[i][r][j][k] = (TH1D*)hProjInvMassUS[i][r][j][k]->Clone();
					// hProjInvMassLS[i][r][j][k]->SetName(Form("hProjInvMassLS_%d_%d_%d_%d",i,r,j,k));
					// hProjInvMassLS[i][r][j][k]->Reset();
					hProjInvMassUS[j][k] = hInvMassUS->Projection(0,"e");
					hProjInvMassUS[j][k]->SetName(Form("hInvMassUS_%d_%d",j,k));
					hProjInvMassLSpp[j][k] = hInvMassLSpp->Projection(0,"e");
					hProjInvMassLSpp[j][k]->SetName(Form("hProjInvMassLSpp_%d_%d",j,k));
					hProjInvMassLSmm[j][k] = hInvMassLSmm->Projection(0,"e");
					hProjInvMassLSmm[j][k]->SetName(Form("hProjInvMassLSmm_%d_%d",j,k));
					hProjInvMassLS[j][k] = (TH1D*)hProjInvMassUS[j][k]->Clone();
					hProjInvMassLS[j][k]->SetName(Form("hProjInvMassLS_%d_%d",j,k));
					hProjInvMassLS[j][k]->Reset();
					
					// for (int p = 0 ; p < hProjInvMassUS[i][r][j][k]->GetNbinsX() ; p++) {
					// 	hProjInvMassLS[i][r][j][k]->SetBinContent( p+1, 2.0*sqrt(
					// 		hProjInvMassLSpp[i][r][j][k]->GetBinContent(p+1)*hProjInvMassLSmm[i][r][j][k]->GetBinContent(p+1) ) );
					// 	hProjInvMassLS[i][r][j][k]->SetBinError( p+1, sqrt( 
					// 		hProjInvMassLSpp[i][r][j][k]->GetBinContent(p+1)*hProjInvMassLSmm[i][r][j][k]->GetBinContent(p+1) ) * (
					// 		pow( hProjInvMassLSpp[i][r][j][k]->GetBinError(p+1)/hProjInvMassLSpp[i][r][j][k]->GetBinContent(p+1),2 ) +
					// 		pow( hProjInvMassLSmm[i][r][j][k]->GetBinError(p+1)/hProjInvMassLSmm[i][r][j][k]->GetBinContent(p+1),2 ) ) );
					// }

					for (int p = 0 ; p < hProjInvMassUS[j][k]->GetNbinsX() ; p++) {
						hProjInvMassLS[j][k]->SetBinContent( p+1, 2.0*sqrt(
							hProjInvMassLSpp[j][k]->GetBinContent(p+1)*hProjInvMassLSmm[j][k]->GetBinContent(p+1) ) );
						hProjInvMassLS[j][k]->SetBinError( p+1, sqrt( 
							hProjInvMassLSpp[j][k]->GetBinContent(p+1)*hProjInvMassLSmm[j][k]->GetBinContent(p+1) ) * (
							pow( hProjInvMassLSpp[j][k]->GetBinError(p+1)/hProjInvMassLSpp[j][k]->GetBinContent(p+1),2 ) +
							pow( hProjInvMassLSmm[j][k]->GetBinError(p+1)/hProjInvMassLSmm[j][k]->GetBinContent(p+1),2 ) ) );
					}					
					hProjInvMassSub[j][k] = (TH1D*)hProjInvMassUS[j][k]->Clone();
					hProjInvMassSub[j][k]->SetName(Form("hProjInvMassSub_%d_%d",j,k));
					hProjInvMassSub[j][k]->Add( hProjInvMassLS[j][k], -1.0 );
				}
			}
		// }
	// }

	// for(int i=0;i<nlpt;i++) {
	// 	for(int j=0;j<nmult;j++) {
	// 		for(int k=0;k<npt;k++) {
	// 			hProjInvMassSub_InJet[i][j][k] = (TH1D*)hProjInvMassSub[i][0][j][k]->Clone();
	// 			hProjInvMassSub_InJet[i][j][k]->SetName(Form("hProjInvMassSub_InJet_%d_%d_%d",i,j,k));
	// 			hProjInvMassSub_InJet[i][j][k]->Add( hProjInvMassSub[i][1][j][k], -1.0 );
	// 		}
	// 	}
	// }

	TFile* fout = new TFile("InvMassOut.root","recreate");

	TLatex latex;

	// for (int i = 0; i < nlpt; i++) {
	// 	for (int r = 0; r < nrt; r++) {
			for (int j = 0; j < nmult; j++) {
				// TCanvas *c = new TCanvas(Form("c_%d_%d_%d", i, r, j), Form("Canvas %d_%d_%d", i, r, j), 1000, 1000);
				TCanvas *c = new TCanvas(Form("c_%d", j), Form("Canvas %d", j), 1000, 1000);
				c->Divide(4,4,0.001,0.001);
    			c->SetGrid();

				for (int k = 0; k < npt; k++) {
					const int cIdx = k+1;
    				c->cd(cIdx);
					gStyle->SetOptTitle(0);
    				gStyle->SetOptStat(0);

					hProjInvMassSub[j][k]->Draw("hist");
					hProjInvMassSub[j][k]->GetXaxis()->SetTitle("M_{#pi#pi} (GeV/c^{2})");
					hProjInvMassSub[j][k]->GetXaxis()->SetRangeUser(0.4, 1.7);

					std::stringstream line1, line2, line3, line4;
					line1 << "LHC22o_apass6";
					line2 << "FT0M " << m_min[j] << " - " << m_max[j] << " %";
					line3 << "pp, 13.6 TeV";
					line4 << p_min[k] << " < p_{T} < " << p_max[k] << " GeV/c";

					std::string line1Str = line1.str();
					std::string line2Str = line2.str();
					std::string line3Str = line3.str();
					std::string line4Str = line4.str();

					latex.SetTextSize(0.03);
					latex.SetTextColor(kBlack);
					latex.DrawLatexNDC(0.6, 0.75, line1Str.c_str());
					latex.DrawLatexNDC(0.6, 0.70, line2Str.c_str());
					latex.DrawLatexNDC(0.6, 0.65, line3Str.c_str());
					latex.DrawLatexNDC(0.6, 0.60, line4Str.c_str());

				}
				TString fileName = Form("results_plot/202407_mult_%d_%d.pdf", m_min[j], m_max[j]);
				c->SaveAs(fileName);
				delete c;
			}
	// 	}
	// }
}
