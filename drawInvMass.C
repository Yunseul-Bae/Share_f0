#include <TCanvas.h>
#include <TLatex.h>
#include <TString.h>
#include <sstream>
#include <iomanip>
#include <string>
#include <iostream>

void drawInvMass() {
	TFile* fin = new TFile("/Volumes/Seul/HP_results/239395.root");
	// TFile* fin = new TFile("~/Downloads/AnalysisResults.root","read");

	THnSparse* hInvMassUS = (THnSparse*)fin->Get("lf-f0980analysis/hInvMass_f0980_US");
	THnSparse* hInvMassLSpp = (THnSparse*)fin->Get("lf-f0980analysis/hInvMass_f0980_LSpp");
	THnSparse* hInvMassLSmm = (THnSparse*)fin->Get("lf-f0980analysis/hInvMass_f0980_LSmm");

	const int nmult = 7;
	const int npt = 16;

	//	0-1, 1-5, 5-10, 10-20, 20-30, 30-50, 50-100
	int m_min[nmult] = {0, 1, 5, 10, 20, 30, 50};
	int m_max[nmult] = {1, 5, 10, 20, 30, 50, 100};
	//	0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 13.0
	double p_min[npt] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0};
	double p_max[npt] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 13.0};

	TH1D* hProjInvMassUS[nmult][npt];
	TH1D* hProjInvMassLSpp[nmult][npt];
	TH1D* hProjInvMassLSmm[nmult][npt];
	TH1D* hProjInvMassLS[nmult][npt];
	TH1D* hProjInvMassSub[nmult][npt];
	TH1D* hProjInvMassSub_D[nmult][npt];
	TH1D* hProjInvMassSub_S[nmult][npt];

	TCanvas *c1 = new TCanvas("c1", "Invariant mass distribution", 1400, 600);
	c1->Divide(2,1,0.001,0.001);

	//	centrality = j
	for (int j = 0 ; j < nmult ; j++) {
		hInvMassUS->GetAxis(2)->SetRangeUser( m_min[j], m_max[j]-0.001 );
		hInvMassLSpp->GetAxis(2)->SetRangeUser( m_min[j], m_max[j]-0.001 );
		hInvMassLSmm->GetAxis(2)->SetRangeUser( m_min[j], m_max[j]-0.001 );

		for (int k = 0 ; k < npt ; k++) {
			hInvMassUS->GetAxis(1)->SetRangeUser( p_min[k], p_max[k]-0.001 );
			hInvMassLSpp->GetAxis(1)->SetRangeUser( p_min[k], p_max[k]-0.001 );
			hInvMassLSmm->GetAxis(1)->SetRangeUser( p_min[k], p_max[k]-0.001 );

			hProjInvMassUS[j][k] = hInvMassUS->Projection(0,"e");
			hProjInvMassUS[j][k]->SetName(Form("hInvMassUS_%d_%d",j,k));
			hProjInvMassLSpp[j][k] = hInvMassLSpp->Projection(0,"e");
			hProjInvMassLSpp[j][k]->SetName(Form("hProjInvMassLSpp_%d_%d",j,k));
			hProjInvMassLSmm[j][k] = hInvMassLSmm->Projection(0,"e");
			hProjInvMassLSmm[j][k]->SetName(Form("hProjInvMassLSmm_%d_%d",j,k));
			hProjInvMassLS[j][k] = (TH1D*)hProjInvMassUS[j][k]->Clone();
			hProjInvMassLS[j][k]->SetName(Form("hProjInvMassLS_%d_%d",j,k));
			hProjInvMassLS[j][k]->Reset();

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
			
			//	Draw
			if ( j == 5 && k == 7) {
				c1->cd(1);
				gStyle->SetOptTitle(0);
				gStyle->SetOptStat(0);

				hProjInvMassSub_D[j][k] = (TH1D*)hProjInvMassSub[j][k]->Clone();

				hProjInvMassSub_D[j][k]->Draw("");
				hProjInvMassSub_D[j][k]->SetMarkerColor(kBlack);
				hProjInvMassSub_D[j][k]->SetMarkerStyle(4);
				hProjInvMassSub_D[j][k]->SetMarkerSize(0.5);
				hProjInvMassSub_D[j][k]->GetXaxis()->SetTitle("M_{#pi#pi} (GeV/c^{2})");
				hProjInvMassSub_D[j][k]->GetXaxis()->SetRangeUser(0.4, 1.7);
				hProjInvMassSub_D[j][k]->GetYaxis()->SetTitle("Counts");
				hProjInvMassLS[j][k]->Draw("same");
				hProjInvMassLS[j][k]->SetMarkerColor(kRed);
				hProjInvMassLS[j][k]->SetMarkerStyle(4);
				hProjInvMassLS[j][k]->SetMarkerSize(0.5);

				TLegend* leg = new TLegend(0.55,0.55,0.8,0.85);
				leg->SetTextSize(0.038);
				leg->SetLineWidth(0.0);
				leg->SetFillStyle(0);
				leg->AddEntry( (TObject*)0, "LHC22o apass7", "");
				leg->AddEntry( (TObject*)0, Form("FT0M %d#font[122]{-}%d %%",m_min[j],m_max[j]), "");
				leg->AddEntry( (TObject*)0, "pp 13.6 TeV, |#it{y}| < 0.5", "");
				// leg->AddEntry( (TObject*)0, Form("#it{p}_{T,lead} > 5 GeV/#it{c}, %s",trnsName[r]), "");
				leg->AddEntry( (TObject*)0, Form("%.1lf < #it{p}_{T} < %.1lf GeV/#it{c}",p_min[k],p_max[k]), "");
				leg->AddEntry( hProjInvMassSub_D[j][k], "Unlike-sign pairs","p");
				leg->AddEntry( hProjInvMassLS[j][k], "Like-sign pairs","p");
				leg->Draw();
			}			
			hProjInvMassSub[j][k]->Add( hProjInvMassLS[j][k], -1.0 );

			if ( j == 5 && k == 7) {
				c1->cd(2);
				hProjInvMassSub_S[j][k] = (TH1D*)hProjInvMassSub[j][k]->Clone();
				hProjInvMassSub_S[j][k]->Draw("l");
				
				hProjInvMassSub_S[j][k]->GetXaxis()->SetTitle("M_{#pi#pi} (GeV/c^{2})");
				hProjInvMassSub_S[j][k]->GetXaxis()->SetRangeUser(0.4, 1.7);
				hProjInvMassSub_S[j][k]->GetYaxis()->SetTitle("Counts");
				hProjInvMassSub_S[j][k]->SetMarkerColor(kBlue);
				hProjInvMassSub_S[j][k]->SetMarkerStyle(4);
				hProjInvMassSub_S[j][k]->SetMarkerSize(0.8);
				
				c1->SaveAs("plot/bg_subtraction.pdf");
			}
		}
	}

	TFile* fout = new TFile("/Users/seul/seul_workspace/f0ana/results/ana_results/InvMassOut.root","recreate");

	for (int j=0;j<nmult;j++) {
		for (int k=0;k<npt;k++) {
			hProjInvMassSub[j][k]->Write();
		}
	}

	// TLatex latex;

	for (int j = 0; j < nmult; j++) {
		// TCanvas *c = new TCanvas(Form("c_%d_%d_%d", i, r, j), Form("Canvas %d_%d_%d", i, r, j), 1000, 1000);
		TCanvas *c = new TCanvas(Form("c_%d", j), Form("Canvas %d", j), 1500, 700);
		c->Divide(6,3,0.001,0.001);
		c->SetGrid();

		for (int k = 0; k < npt; k++) {
			const int cIdx = k+1;
			c->cd(cIdx);
			gStyle->SetOptTitle(0);
			gStyle->SetOptStat(0);

			hProjInvMassSub[j][k]->Draw("hist");
			hProjInvMassSub[j][k]->GetXaxis()->SetTitle("M_{#pi#pi} (GeV/c^{2})");
			hProjInvMassSub[j][k]->GetXaxis()->SetRangeUser(0.4, 1.7);
			
			TLegend* leg = new TLegend(0.4,0.6,0.9,0.85);
			// leg->SetTextFont(43);
			leg->SetTextSize(0.045);
			leg->SetLineWidth(0.0);
			leg->SetFillStyle(0);
			leg->AddEntry( (TObject*)0, "LHC22o apass6", "");
			leg->AddEntry( (TObject*)0, Form("FT0M %d#font[122]{-}%d %%",m_min[j],m_max[j]), "");
			leg->AddEntry( (TObject*)0, "pp 13.6 TeV, |#it{y}| < 0.5", "");
			// leg->AddEntry( (TObject*)0, Form("#it{p}_{T,lead} > 5 GeV/#it{c}, %s",trnsName[r]), "");
			leg->AddEntry( (TObject*)0, Form("%.1lf < #it{p}_{T} < %.1lf GeV/#it{c}",p_min[k],p_max[k]), "");
			leg->Draw();

		}
		TString fileName = Form("plot/Invmass_mult_%d_%d.pdf", m_min[j], m_max[j]);
		c->SaveAs(fileName);
		delete c;
	}
}
