#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <iostream>
#include <cmath>
#include <cstring>

Double_t threshold = 0.139570*2.;

double PhaseSpaceFactor (double x, double T, double pT) {
	double F = ( x / ( sqrt( pow(x,2) + pow(pT,2) ) ) ) * exp( -sqrt( pow(x,2) + pow(pT,2) )/T );
 	return F;
}

double Width (double x, double M0, double Gam0, int Spin) {
	double w = pow( (x*x - threshold*threshold)/(M0*M0 - threshold*threshold),0.5+(double)Spin )*Gam0*M0/x;
 	return w;
}

double breitWigner (double x, double Amp, double M0, double Gam0, int Spin) {
 	double br = Amp*x*M0*Width(x, M0, Gam0, Spin);
 	br /= ( pow(M0*M0-x*x,2) + M0*M0*pow(Width(x, M0, Gam0, Spin),2) );
 	return br;
}

double background (double x, double ind, double b1, double b2) {
 	double bg = pow(x-threshold,ind)*exp( b1*(x-threshold) + b2*pow(x-threshold,2) );
 	return bg;
}

// double breitWigner_f0 (double *x, double* par){
//  	return breitWigner(x[0], par[1], par[0], par[2], 0);
// }

double sfit (double *x, double* par) {

	double bgA = par[0];
	double bgind = par[1];
	double bgb1 = par[2];
	double bgb2 = par[3];

	double f0m = par[4];
	double f0A = par[5];
	double f0g = par[6];
	int f0s = 0;

	double f2m = par[7];
	double f2A = par[8];
	double f2g = par[9];
	int f2s = 2;

	double rhom = par[10];
	double rhoA = par[11];
	double rhog = par[12];
	int rhos = 1;

	double T = par[13];
	double pT = par[14];

	double breitWigner_f0  = breitWigner(x[0],f0A,f0m,f0g,f0s);
	double breitWigner_f2  = breitWigner(x[0],f2A,f2m,f2g,f2s);
	double breitWigner_rho  = breitWigner(x[0],rhoA,rhom,rhog,rhos);
	double bgfunc = bgA*background(x[0],bgind,bgb1,bgb2);

	return ( breitWigner_f0 + breitWigner_f2 + breitWigner_rho )*PhaseSpaceFactor(x[0],T,pT) + bgfunc;
	// return breitWigner_f0 + breitWigner_f2 + breitWigner_rho + bgfunc;
}

//Main fitting function
void FitLP (int LP=0, int Rebinning=1) {
	TFile* fin = new TFile("/Users/seul/seul_workspace/f0ana/data/ana_results/InvMassOut.root","read");

	const int nmult = 7; //	4
	const int npt = 16;  //	4

	int m_min[nmult] = {0, 1, 5, 10, 20, 30, 50};
	int m_max[nmult] = {1, 5, 10, 20, 30, 50, 100};

	double p_min[npt] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0};
	double p_max[npt] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 13.0};

	TH1D* hMassSig[nmult][npt];
	TF1* fitBin[nmult][npt];
	TF1* fitBinComp[nmult][npt][4];

	TFitResultPtr fitResults[nmult][npt];

	TGraphErrors* gf0Width[nmult];
	TGraphErrors* gf0Mass[nmult];
	TGraphErrors* gChi2OverNDF[nmult];

	for(int i=0;i<nmult;i++){
		gf0Width[i] = new TGraphErrors();
		gf0Mass[i] = new TGraphErrors();
		gChi2OverNDF[i] = new TGraphErrors();
	}

	double Intgr_RawY_cntl[nmult][npt];
	double Intgr_RawY_stat[nmult][npt];
	TH1D* hPtStat[nmult];
	double PtBinnings[npt] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 13.0};

	int MagIndex[4] = {0,5,8,11};
	int lcolor[4] = {46, 30, 41, 28};

	double BGAmpMin = 1e2*2.0*1e1;
	double BGAmpMax = 3e7*2.0*1e1;

	double BGSlopeMin = -2.5;
	double BGSlopeMax = 3.0;

	double BGInd1Min = -5.0;
	double BGInd1Max = 5.0;

	double BGInd2Min = -0.5;
	double BGInd2Max = 5.0;

	double f0MassMinRange = 0.95;
	double f0MassMaxRange = 1.0;

	double f0AmpMin = 1.0*2.0;
	double f0AmpMax = 1e6*2.0;

	double f0WMin = 0.01;
	double f0WMax = 0.10;

	double f2MassMinRange = 1.2755-0.0012*0.002;
	double f2MassMaxRange = 1.2755+0.0012*0.002;

	double f2AmpMin = 1.0*2.0;
	double f2AmpMax = 1e6*2.0;

	double f2WMin = 0.1867-0.0024*0.02;
	double f2WMax = 0.1867+0.0029*0.02;

	double rhoAmpMin = 1.0*2.0;
	double rhoAmpMax = 1e6*2.0;

	double rhoMassMin = 0.7752;
	double rhoMassMax = 0.7754;

	double rhoWMin = 0.1491-0.00009;
	double rhoWMax = 0.1491+0.00009;

	for (int i=0;i<nmult;i++) {
		for (int j=0;j<npt;j++) {
			hMassSig[i][j] = (TH1D*)fin->Get(Form("hProjInvMassSub_%d_%d",i,j));
			// hMassSig[3][i][j] = (TH1D*)fin->Get(Form("hProjInvMassSub_InJet_%d_%d_%d",LP,i,j));
			// if (!hMassSig[i][j]) {
			// 	std::cerr << "Error: Histogram not found for LP=" << LP << ", i=" << i << ", j=" << j << std::endl;
			// continue;
			// }
		}
	}

	double nsigIntgr = 3.;

	// for(int r=0;r<nrt;r++){
	for (int i=0 ; i<nmult ; i++) {
		//	canvas
		// TCanvas* c = new TCanvas("c","c",800,600);
		TCanvas *c = new TCanvas(Form("c_%d", i), Form("Canvas %d", i), 1600, 1200);
			c->Divide(4,4,0.001,0.001);
			// c->Divide(4,4);
			// c->SetGrid();

		//	legend
		// TLegend* leg = new TLegend(0.45,0.4,0.95,0.95);	

		for (int j=0 ; j<npt ; j++) {
			//	canvas
			const int cIdx = j+1;
			c->cd(cIdx);
			// gPad->SetLeftMargin(0.14);
			// gPad->SetBottomMargin(0.14);
			// gPad->SetRightMargin(0.03);
			// gPad->SetTopMargin(0.05);
			gPad->SetTicks();
			gStyle->SetOptTitle(0);
			gStyle->SetOptStat(0);

			fitBin[i][j] = new TF1("ffit",sfit,0.7,1.8,15);

			for (int k=0;k<4;k++) {
				fitBinComp[i][j][k] = new TF1("ffit",sfit,0.7,1.8,15);
			}
			fitBin[i][j]->FixParameter(13, 0.160 );
			fitBin[i][j]->FixParameter(14, (p_min[j]+p_max[j])/2. );

			fitBin[i][j]->SetParameter(1, 2.0 );
			fitBin[i][j]->SetParameter(2, -4.0 );
			fitBin[i][j]->SetParameter(3, 0.0 );

			fitBin[i][j]->SetParLimits(0,hMassSig[i][j]->GetMaximum()*0.1,BGAmpMax);
			fitBin[i][j]->SetParLimits(1,BGSlopeMin,BGSlopeMax);
			fitBin[i][j]->SetParLimits(2,BGInd1Min,BGInd1Max);
			fitBin[i][j]->SetParLimits(3,BGInd2Min,BGInd2Max);

			fitBin[i][j]->SetParLimits(4,f0MassMinRange,f0MassMaxRange);
			fitBin[i][j]->SetParLimits(6,0.054,0.056);

			if ( LP==0 ) {
				fitBin[i][j]->SetParLimits(6, 0.01, 0.1);
			}
			fitBin[i][j]->SetParLimits(7,f2MassMinRange,f2MassMaxRange);
			fitBin[i][j]->SetParLimits(9,f2WMin,f2WMax);

			fitBin[i][j]->SetParLimits(10,rhoMassMin,rhoMassMax);
			fitBin[i][j]->SetParLimits(12,rhoWMin,rhoWMax);

			fitBin[i][j]->FixParameter(7, (f2MassMinRange+f2MassMaxRange)/2. );
			fitBin[i][j]->FixParameter(9, (f2WMin+f2WMax)/2. );

			fitBin[i][j]->FixParameter(10, (rhoMassMin+rhoMassMax)/2. );
			fitBin[i][j]->FixParameter(12, (rhoWMin+rhoWMax)/2. );

			if ( j==0 ) {
				f0AmpMin = 1.0*2.0*1e1;
				f2AmpMin = 1.0*2.0*1e1;
				rhoAmpMin = 1.0*2.0*1e1;

				rhoAmpMax = hMassSig[i][j]->GetMaximum()*1e4;
				f0AmpMax = hMassSig[i][j]->GetMaximum()*1e3;
				f2AmpMax = hMassSig[i][j]->GetMaximum()*1e2;

			} else if( j>0 ){
				f0AmpMin *= PhaseSpaceFactor(1, 0.160, (p_min[j-1]+p_max[j-1])/2.0 ) / PhaseSpaceFactor(1, 0.160, (p_min[j]+p_max[j])/2.0 );
				f0AmpMax *= PhaseSpaceFactor(1, 0.160, (p_min[j-1]+p_max[j-1])/2.0 ) / PhaseSpaceFactor(1, 0.160, (p_min[j]+p_max[j])/2.0 );
				f2AmpMin *= PhaseSpaceFactor(1, 0.160, (p_min[j-1]+p_max[j-1])/2.0 ) / PhaseSpaceFactor(1, 0.160, (p_min[j]+p_max[j])/2.0 );
				f2AmpMax *= PhaseSpaceFactor(1, 0.160, (p_min[j-1]+p_max[j-1])/2.0 ) / PhaseSpaceFactor(1, 0.160, (p_min[j]+p_max[j])/2.0 );
				rhoAmpMin *= PhaseSpaceFactor(1, 0.160, (p_min[j-1]+p_max[j-1])/2.0 ) / PhaseSpaceFactor(1, 0.160, (p_min[j]+p_max[j])/2.0 );
				rhoAmpMax *= PhaseSpaceFactor(1, 0.160, (p_min[j-1]+p_max[j-1])/2.0 ) / PhaseSpaceFactor(1, 0.160, (p_min[j]+p_max[j])/2.0 );
			}

			fitBin[i][j]->SetParLimits(5,f0AmpMin, f0AmpMax );
			fitBin[i][j]->SetParLimits(8,f2AmpMin, f2AmpMax);
			fitBin[i][j]->SetParLimits(11,rhoAmpMin, rhoAmpMax);

			if ( Rebinning>1 ) {
				hMassSig[i][j]->Rebin( Rebinning );
			}
			if ( j>8 ) {
				hMassSig[i][j]->Rebin( 2 );
			}
			fitResults[i][j] = (TFitResultPtr)hMassSig[i][j]->Fit(fitBin[i][j],"sbq0","",0.8,1.55);
			// if (!fitResults[i][j]) {
			// 	std::cerr << "Error: Fit failed for i=" << i << ", j=" << j << std::endl;
			// continue;
			// }

			for(int k=1;k<4;k++){
				if( k==7 || k==9 || k==10 || k==12 ) continue;
				fitBin[i][j]->FixParameter( k, fitBin[i][j]->GetParameter(k) );
			}
			fitResults[i][j] = (TFitResultPtr)hMassSig[i][j]->Fit(fitBin[i][j],"sb","",0.8,1.55);

			if (j<npt-1){
				gf0Width[i]->SetPoint( j, (p_min[j]+p_max[j])/2., fitBin[i][j]->GetParameter(6) );
				gf0Width[i]->SetPointError( j, 0, fitBin[i][j]->GetParError(6) );

				gf0Mass[i]->SetPoint( j, (p_min[j]+p_max[j])/2., fitBin[i][j]->GetParameter(4) );
				gf0Mass[i]->SetPointError( j, 0, fitBin[i][j]->GetParError(4) );

				gChi2OverNDF[i]->SetPoint( j, (p_min[j]+p_max[j])/2., fitBin[i][j]->GetChisquare()/fitBin[i][j]->GetNDF() );
			}
			
			for (int k=0;k<4;k++) {
				for (int p=0;p<15;p++) {
					fitBinComp[i][j][k]->SetParameter( p, fitBin[i][j]->GetParameter(p) );
				}
				for(int p=0;p<4;p++){
					if( k==p ) continue;
					fitBinComp[i][j][k]->SetParameter( MagIndex[p], 0 );
				}
			}

			hMassSig[i][j]->Draw("*");

			hMassSig[i][j]->GetXaxis()->SetRangeUser(0.8,1.85);
			hMassSig[i][j]->GetXaxis()->SetTitle("#it{M}_{#pi#pi} (GeV/#it{c}^{2})");
			hMassSig[i][j]->GetYaxis()->SetTitle(Form("Counts / %.1lf MeV",hMassSig[i][j]->GetBinWidth(1)*1e3));
			// hMassSig[i][j]->SetTitle("");

			// hMassSig[i][j]->GetXaxis()->SetTitleFont(43);
			// hMassSig[i][j]->GetXaxis()->SetLabelFont(43);
			// hMassSig[i][j]->GetYaxis()->SetTitleFont(43);
			// hMassSig[i][j]->GetYaxis()->SetLabelFont(43);

			// hMassSig[i][j]->GetXaxis()->SetTitleSize(32);
			// hMassSig[i][j]->GetYaxis()->SetTitleSize(32);
			// hMassSig[i][j]->GetXaxis()->SetLabelSize(28);
			// hMassSig[i][j]->GetYaxis()->SetLabelSize(28);

			// hMassSig[i][j]->SetMarkerStyle(24);
			// hMassSig[i][j]->SetMarkerStyle(7);
			// hMassSig[i][j]->SetMarkerColor(1);
			hMassSig[i][j]->SetLineColor(1);
			hMassSig[i][j]->GetFunction("ffit")->SetLineColor(kBlue);
			hMassSig[i][j]->GetFunction("ffit")->SetLineWidth(1);
			hMassSig[i][j]->SetMinimum(0);
			hMassSig[i][j]->SetMaximum( hMassSig[i][j]->GetBinContent( hMassSig[i][j]->GetXaxis()->FindBin(0.8) )*1.1 );

			for (int k=0;k<4;k++) {
				fitBinComp[i][j][k]->SetLineWidth(1);
				fitBinComp[i][j][k]->SetLineColor( lcolor[k] );
				if (k!=1) {
					fitBinComp[i][j][k]->SetLineStyle(3);
				}
				fitBinComp[i][j][k]->Draw("same");
			}

			// leg->Clear();
			TLegend* leg = new TLegend(0.45,0.4,0.8,0.8);

			leg->SetTextFont(43);
			// leg->SetTextSize(26);
			leg->SetLineWidth(0.0);
			leg->SetFillStyle(0);
			// leg->AddEntry( (TObject*)0, "ALICE work in progress", "");
			leg->AddEntry( (TObject*)0, "LHC22o apass6", "");
			leg->AddEntry( (TObject*)0, Form("FT0M %d#font[122]{-}%d %%",m_min[i],m_max[i]), "");
			leg->AddEntry( (TObject*)0, "pp 13.6 TeV, |#it{y}| < 0.8", "");
			// leg->AddEntry( (TObject*)0, Form("#it{p}_{T,lead} > 5 GeV/#it{c}, %s",trnsName[r]), "");
			leg->AddEntry( (TObject*)0, Form("%.1lf < #it{p}_{T} < %.1lf GeV/#it{c}",p_min[j],p_max[j]), "");
			leg->AddEntry( fitBinComp[i][j][0], "Residual bkg.", "l");
			leg->AddEntry( fitBinComp[i][j][1], "f_{0}(980)", "l");
			leg->AddEntry( fitBinComp[i][j][3], "#rho(770)^{0}", "l");
			leg->AddEntry( fitBinComp[i][j][2], "f_{2}(1270)", "l");

			leg->Draw();

			// c->SaveAs(Form("figs/Fit_%dLP/fitout_%d_%d.pdf",LP,i,j));
			// c->SaveAs(Form("plot/Fit_%dLP/fitout_%d_%d.pdf",LP,i,j));		

			fitBin[i][j]->SetParameter(0, 0);
			fitBin[i][j]->SetParameter(8, 0);
			fitBin[i][j]->SetParameter(11, 0);

			Intgr_RawY_cntl[i][j] = fitBin[i][j]->Integral(
				fitBin[i][j]->GetParameter(4) - fitBin[i][j]->GetParameter(6)*nsigIntgr,
				fitBin[i][j]->GetParameter(4) + fitBin[i][j]->GetParameter(6)*nsigIntgr );
			Intgr_RawY_cntl[i][j] /= hMassSig[i][j]->GetBinWidth(1);

			Intgr_RawY_stat[i][j] = fitBin[i][j]->IntegralError(
				fitBin[i][j]->GetParameter(4) - fitBin[i][j]->GetParameter(6)*nsigIntgr,
				fitBin[i][j]->GetParameter(4) + fitBin[i][j]->GetParameter(6)*nsigIntgr,
				fitBin[i][j]->GetParameters(), fitResults[i][j]->GetCovarianceMatrix().GetMatrixArray(), 1e-6 );
			Intgr_RawY_stat[i][j] /= hMassSig[i][j]->GetBinWidth(1);
		}
		TString fileName = Form("plot/Fit_%dLP/fitout_mult_%d_%d.pdf",LP, m_min[i], m_max[i]);
		c->SaveAs(fileName);
		delete c;

		// break;
	}
	// }

	// for(int r=0;r<3;r++){
	for (int i=0;i<nmult;i++) {
		hPtStat[i] = new TH1D(Form("hPtStat_%d",i),"",npt-1,PtBinnings);
		for (int j=0;j<npt-1;j++) {
			hPtStat[i]->SetBinContent( j+1, Intgr_RawY_cntl[i][j] );
			hPtStat[i]->SetBinError( j+1, Intgr_RawY_stat[i][j] );
		}
	}
	// }

	TFile* fout = new TFile(Form("/Users/seul/seul_workspace/f0ana/data/ana_results/F0IntgrOut_%dLP.root",LP),"recreate");
	// for(int r=0;r<3;r++){
	for (int i=0;i<nmult;i++) {
		hPtStat[i]->Write();
	}
	// }
}

void Fit() {
 	FitLP(0, 1);
 	// FitLP(1, 2);
}
