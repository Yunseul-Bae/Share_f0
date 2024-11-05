#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <cmath>
#include <ctime>
#include <sstream>
#include <iomanip>

double_t threshold = 0.139570*2.;   //  charged pion mass = 139.6 MeV

//  PS (Invariant mass, Kinetic freeze-out temperature = 160 MeV, pT)
double PhaseSpaceFactor(double x, double T, double pT) {
	double F = (x / (sqrt(pow(x,2) + pow(pT,2)))) * exp(-sqrt(pow(x,2) + pow(pT,2)) / T);
 	return F;
}

//  Width (Invariant mass, rest mass, Gamma0:rest width, spin)
double Width(double x, double M0, double Gam0, int Spin) {
	double w = pow((x*x - threshold*threshold) / (M0*M0-threshold*threshold),0.5+(double)Spin) * Gam0 * M0 / x;
 	return w;
}

//  rBW (Invariant mass, Amplitude, Rest mass, Width, Spin)
double breitWigner(double x, double Amp, double M0, double Gam0, int Spin) {
 	double br_u = Amp * x * Width(x, M0, Gam0, Spin) * M0;
    double br_d = (pow(x*x - M0*M0,2) + M0*M0*pow(Width(x, M0, Gam0, Spin),2));
    double br = br_u/br_d;
 	// br /= (pow(M0*M0 - x*x,2) + M0*M0*pow(Width(x, M0, Gam0, Spin),2));
 	// br /= (pow(x*x - M0*M0,2) + M0*M0*pow(Width(x, M0, Gam0, Spin),2));
 	return br;
}

//  Maxwell-Boltzmann-like distribution (Invariant mass, index)
double background(double x, double ind, double b1, double b2) {
 	// double bg = pow(x-threshold,ind) * exp( b1*(x-threshold) + b2*pow(x-threshold,2)); //  Ana note
 	double bg = pow(x-threshold,ind) * exp( b1*x + b2*x*x); //  Paper
 	return bg;
}

double fitting(double* x, double* par) {
    
    double bg_A = par[0];	//  Amplitude - B (free parameter)
	double bg_ind = par[1]; //  index
	double bg_b1 = par[2];
	double bg_b2 = par[3];

	double f0_m = par[4];	//  Invariant mass
	double f0_A = par[5];
	double f0_g = par[6]; 	//  Gamma0 - width
	int f0_s = 0; 			//  spin

	double f2_m = par[7];
	double f2_A = par[8];
	double f2_g = par[9];
	int f2_s = 2;

	double rho_m = par[10];
	double rho_A = par[11];
	double rho_g = par[12];
	int rho_s = 1;

	double T = par[13];		//  Kinetic freeze-out temperature
	double pT = par[14];

    double breitWigner_f0 = breitWigner(x[0],f0_A,f0_m,f0_g,f0_s);
	double breitWigner_f2 = breitWigner(x[0],f2_A,f2_m,f2_g,f2_s);
	double breitWigner_rho = breitWigner(x[0],rho_A,rho_m,rho_g,rho_s);
	double bgfunc = bg_A*background(x[0],bg_ind,bg_b1,bg_b2);

    // double total = (breitWigner_f0 + breitWigner_f2 + breitWigner_rho) * PhaseSpaceFactor(x[0], T, pT) + bgfunc;
    // return total;
	return (breitWigner_f0 + breitWigner_f2 + breitWigner_rho) * PhaseSpaceFactor(x[0], T, pT) + bgfunc;
}

void Fit_outPrm() {
    TFile* Fin = new TFile("/Users/seul/seul_workspace/f0ana/results/InvMassOut.root","read");

    //  Multiplicity binning
    const int nmult = 8;
    int mult_min[nmult] = {0, 1, 5, 10, 20, 30, 50, 0};
	int mult_max[nmult] = {1, 5, 10, 20, 30, 50, 100, 100};
    //  pT binning
    const int npt = 16;
    double pt_min[npt] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0};
	double pt_max[npt] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 13.0};

    TH1D* hMassSig[nmult][npt];             //  Get the mass signal
    TF1* fitBin[nmult][npt];                //  Fitting function
    TF1* fitBinComp[nmult][npt][4];         //  3 resonance + bg
	TFitResultPtr fitResults[nmult][npt];
	TGraphErrors* gf0Width_init[nmult];		//	Width of initial fitting function
	TGraphErrors* gf0Width[nmult];          //  Width of f0 fitting function
	TGraphErrors* gf0Mass[nmult];           //  f0 mass - peak of f0 fitting function
	TGraphErrors* gf0Chi2perNDF[nmult];     //  chi^2/ndf of fitting function

    // double nSigIntgr = 3.; //  nSigma - detector sigma?
    int MagIndex[4] = {0,5,8,11};
    //  int lcolor[4] = {46, 30, 41, 28};
    int lcolor[4] = {30, 2, 38, 28};

    // double Intgr_RawY_cntl[nmult][npt];
	// double Intgr_RawY_stat[nmult][npt];
	// TH1D* hPtStat[nmult];
	// double PtBinnings[npt] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 13.0};

    int i, j;
    for (i = 0; i < nmult; i++) {
		gf0Width_init[i] = new TGraphErrors();
        gf0Width[i] = new TGraphErrors();
        gf0Mass[i] = new TGraphErrors();
        gf0Chi2perNDF[i] = new TGraphErrors();
    }

	double BG_AmpMin = 1e2*2.0*1e1;
	double BG_AmpMax = 3e7*2.0*1e1;
	double BG_SlopeMin = -2.5;
	double BG_SlopeMax = 3.0;
	double BG_Ind1Min = -5.0;
	double BG_Ind1Max = 5.0;
	double BG_Ind2Min = -0.5;
	double BG_Ind2Max = 5.0;

	double f0_MassMin = 0.95;
	double f0_MassMax = 1.0;
	double f0_AmpMin = 1.0*2.0;
	double f0_AmpMax = 1e6*2.0;
	double f0_WMin = 0.01;
	double f0_WMax = 0.10;

	double f2_MassMin = 1.2755-0.0012*0.002;
	double f2_MassMax = 1.2755+0.0012*0.002;
	double f2_AmpMin = 1.0*2.0;
	double f2_AmpMax = 1e6*2.0;
	double f2_WMin = 0.1867-0.0024*0.02;
	double f2_WMax = 0.1867+0.0029*0.02;

	double rho_MassMin = 0.7752;
	double rho_MassMax = 0.7754;
	double rho_AmpMin = 1.0*2.0;
	double rho_AmpMax = 1e6*2.0;
	double rho_WMin = 0.1491-0.00009;
	double rho_WMax = 0.1491+0.00009;

    for (i = 0; i < nmult; i++) {
		for (j = 0; j < npt; j++) {
			hMassSig[i][j] = (TH1D*)Fin->Get(Form("hProjInvMassSub_%d_%d", i, j));
		}
	}

    for (i = 0; i < nmult; i++) {
        TCanvas *c = new TCanvas(Form("c_%d", i), Form("Canvas %d", i), 3200, 1800);
		// c->Divide(6,3,0.001,0.001);
		c->Divide(4,4,0.001,0.001);

        for (j = 0; j < npt; j++) {
            const int cIdx = j+1;
            c->cd(cIdx);
            gPad->SetTicks();
			gPad->SetRightMargin(0.3);
			gStyle->SetOptTitle(0);
			gStyle->SetOptStat(0);

            fitBin[i][j] = new TF1("ffit", fitting, 0.8, 1.76, 15); //  0.7 ~ 1.8, number of parameter
            for (int k = 0; k < 4; k++) {
                //  0: bg, 1: f0, 2: f2, 3: rho0
				fitBinComp[i][j][k] = new TF1("ffit", fitting, 0.8, 1.76, 15);
			}

			fitBin[i][j]->SetParLimits(0, hMassSig[i][j]->GetMaximum()*0.1, BG_AmpMax); //  bg Amplitude
			fitBin[i][j]->SetParameter(1, 2.0);
            fitBin[i][j]->SetParLimits(1, BG_SlopeMin, BG_SlopeMax);                    //  bg index
			fitBin[i][j]->SetParameter(2, -4.0);
            fitBin[i][j]->SetParLimits(2, BG_Ind1Min, BG_Ind1Max);
            fitBin[i][j]->SetParameter(3, 0.0);
			fitBin[i][j]->SetParLimits(3, BG_Ind2Min, BG_Ind2Max);
			
            fitBin[i][j]->SetParLimits(4,f0_MassMin,f0_MassMax);                        //  f0 mass
			// fitBin[i][j]->SetParLimits(6,0.054,0.056);
			fitBin[i][j]->SetParLimits(6,f0_WMin,f0_WMax);                              //  f0 width
            // fitBin[i][j]->SetParameter(6, 0.055);                                   	//  Initial width at 55 MeV
			// fitBin[i][j]->FixParameter(6, 0.055);
			
            // fitBin[i][j]->SetParLimits(7, f2_MassMin, f2_MassMax);
            fitBin[i][j]->FixParameter(7, (f2_MassMin+f2_MassMax)/2.);
			// fitBin[i][j]->SetParLimits(9, f2_WMin, f2_WMax);
			fitBin[i][j]->FixParameter(9, (f2_WMin+f2_WMax)/2.);

			// fitBin[i][j]->SetParLimits(10, rho_MassMin, rho_MassMax);
            fitBin[i][j]->FixParameter(10, (rho_MassMin+rho_MassMax)/2.);
			// fitBin[i][j]->SetParLimits(12, rho_WMin, rho_WMax);
			fitBin[i][j]->FixParameter(12, (rho_WMin+rho_WMax)/2.);

            fitBin[i][j]->FixParameter(13, 0.160);                                      //  temperature
			fitBin[i][j]->FixParameter(14, (pt_min[j]+pt_max[j])/2.);                   //  pT

            //  Amplitude
            if (j == 0) {
                f0_AmpMin = 1.0*2.0*1e1;
                f0_AmpMax = hMassSig[i][j]->GetMaximum()*1e3;
				f2_AmpMin = 1.0*2.0*1e1;
                f2_AmpMax = hMassSig[i][j]->GetMaximum()*1e2;
				rho_AmpMin = 1.0*2.0*1e1;
				rho_AmpMax = hMassSig[i][j]->GetMaximum()*1e4;
			} else if (j > 0) {
				f0_AmpMin *= PhaseSpaceFactor(1, 0.160, (pt_min[j-1]+pt_max[j-1])/2.) / PhaseSpaceFactor(1, 0.160, (pt_min[j]+pt_max[j])/2.);
				f0_AmpMax *= PhaseSpaceFactor(1, 0.160, (pt_min[j-1]+pt_max[j-1])/2.) / PhaseSpaceFactor(1, 0.160, (pt_min[j]+pt_max[j])/2.);
				f2_AmpMin *= PhaseSpaceFactor(1, 0.160, (pt_min[j-1]+pt_max[j-1])/2.) / PhaseSpaceFactor(1, 0.160, (pt_min[j]+pt_max[j])/2.);
				f2_AmpMax *= PhaseSpaceFactor(1, 0.160, (pt_min[j-1]+pt_max[j-1])/2.) / PhaseSpaceFactor(1, 0.160, (pt_min[j]+pt_max[j])/2.);
				rho_AmpMin *= PhaseSpaceFactor(1, 0.160, (pt_min[j-1]+pt_max[j-1])/2.) / PhaseSpaceFactor(1, 0.160, (pt_min[j]+pt_max[j])/2.);
				rho_AmpMax *= PhaseSpaceFactor(1, 0.160, (pt_min[j-1]+pt_max[j-1])/2.) / PhaseSpaceFactor(1, 0.160, (pt_min[j]+pt_max[j])/2.);
			}
			fitBin[i][j]->SetParLimits(5, f0_AmpMin, f0_AmpMax);
			fitBin[i][j]->SetParLimits(8, f2_AmpMin, f2_AmpMax);
			fitBin[i][j]->SetParLimits(11, rho_AmpMin, rho_AmpMax);

            if (j > 8) { // pT > 2.5
				hMassSig[i][j]->Rebin(4);
			}
			//	Initial fit - 1st step
			fitResults[i][j] = (TFitResultPtr)hMassSig[i][j]->Fit(fitBin[i][j], "sbq0", "", 0.8, 1.55);
			// if (j < npt) {
				gf0Width_init[i]->SetPoint(j, (pt_min[j]+pt_max[j])/2., fitBin[i][j]->GetParameter(6));
				gf0Width_init[i]->SetPointError(j, (pt_max[j]-pt_min[j])/2., fitBin[i][j]->GetParError(6));
			// }

            //  bg_ind, bg_b1, bg_b2
            //  bg 2nd step fit
			for (int k = 1; k < 4; k++) {
				// if ( k==7 || k==9 || k==10 || k==12 ) continue;
				fitBin[i][j]->FixParameter(k, fitBin[i][j]->GetParameter(k));
			}
			// fitResults[i][j] = (TFitResultPtr)hMassSig[i][j]->Fit(fitBin[i][j], "sb", "", 0.8, 1.55);
			fitResults[i][j] = (TFitResultPtr)hMassSig[i][j]->Fit(fitBin[i][j], "sbq", "", 0.8, 1.55);

            // if (j < npt) {
				gf0Width[i]->SetPoint(j, (pt_min[j]+pt_max[j])/2., fitBin[i][j]->GetParameter(6));
				gf0Width[i]->SetPointError(j, (pt_max[j]-pt_min[j])/2., fitBin[i][j]->GetParError(6));
				
				gf0Mass[i]->SetPoint(j, (pt_min[j]+pt_max[j])/2., fitBin[i][j]->GetParameter(4));
				gf0Mass[i]->SetPointError(j, 0, fitBin[i][j]->GetParError(4));

				gf0Chi2perNDF[i]->SetPoint(j, (pt_min[j]+pt_max[j])/2., fitBin[i][j]->GetChisquare()/fitBin[i][j]->GetNDF());
			// }

            // fitBinComp - setting parameters
            for (int k = 0; k < 4; k++) {
				for (int p = 0; p < 15; p++) {
					fitBinComp[i][j][k]->SetParameter( p, fitBin[i][j]->GetParameter(p) );
				}
                //  reset Amplitude
				for(int p = 0; p < 4; p++){
					if( k==p ) continue;
					fitBinComp[i][j][k]->SetParameter( MagIndex[p], 0 );
				}
			}
            //  Draw
            hMassSig[i][j]->Draw("p");
			hMassSig[i][j]->SetMarkerStyle(4);
			// hMassSig[i][j]->GetXaxis()->SetRangeUser(0.8,1.85);
			hMassSig[i][j]->GetXaxis()->SetRangeUser(0.8,1.5);
			hMassSig[i][j]->GetXaxis()->SetTitle("#it{M}_{#pi#pi} (GeV/#it{c}^{2})");
			hMassSig[i][j]->GetYaxis()->SetTitle(Form("Counts / %.1lf MeV/#it{c}^{2}", hMassSig[i][j]->GetBinWidth(1)*1e3));
			hMassSig[i][j]->SetLineColor(1);
			hMassSig[i][j]->GetFunction("ffit")->SetLineColor(kBlue);
			hMassSig[i][j]->GetFunction("ffit")->SetLineWidth(1);
			hMassSig[i][j]->SetMinimum(0);
			hMassSig[i][j]->SetMaximum(hMassSig[i][j]->GetBinContent(hMassSig[i][j]->GetXaxis()->FindBin(0.8))*1.1);

            for (int k = 0; k < 4; k++) {
				fitBinComp[i][j][k]->SetLineWidth(1);
				fitBinComp[i][j][k]->SetLineColor(lcolor[k]);
				
                if (k!=1) {
					fitBinComp[i][j][k]->SetLineStyle(3);
				}
				fitBinComp[i][j][k]->Draw("same");
			}

			// TLegend* leg = new TLegend(0.45,0.45,0.8,0.85);
			TLegend* leg = new TLegend(0.30,0.50,0.8,0.85);
			leg->SetTextFont(43);
			leg->SetLineWidth(0.0);
			leg->SetFillStyle(0);
			// leg->AddEntry( (TObject*)0, "ALICE work in progress", "");
			// leg->AddEntry( (TObject*)0, "LHC22o apass7", "");
			leg->AddEntry( (TObject*)0, Form("FT0M %d#font[122]{-}%d %%",mult_min[i],mult_max[i]), "");
			// leg->AddEntry( (TObject*)0, "pp 13.6 TeV, |#it{y}| < 0.8", "");
			// leg->AddEntry( (TObject*)0, Form("#it{p}_{T,lead} > 5 GeV/#it{c}, %s",trnsName[r]), "");
			leg->AddEntry( (TObject*)0, Form("%.1lf < #it{p}_{T} < %.1lf GeV/#it{c}",pt_min[j],pt_max[j]), "");
			leg->AddEntry( fitBinComp[i][j][0], "Residual bkg.", "l");
			leg->AddEntry( fitBinComp[i][j][1], "f_{0}(980)", "l");
            leg->AddEntry( fitBinComp[i][j][2], "f_{2}(1270)", "l");
			leg->AddEntry( fitBinComp[i][j][3], "#rho(770)^{0}", "l");
			leg->Draw();
			// delete leg;

			//	Latex
			TLatex latex;
			latex.SetTextSize(0.04);
    		latex.SetTextAlign(13);	//	Left align

			double x_ndc = 0.73;
			double y_ndc_start = 0.8;
			double y_ndc_step = 0.06;

			latex.DrawLatexNDC(x_ndc, y_ndc_start + 2*y_ndc_step, Form("#chi^{2}/NDF = %.2f/%d", fitBin[i][j]->GetChisquare(),fitBin[i][j]->GetNDF()));
			latex.DrawLatexNDC(x_ndc+0.07, y_ndc_start + y_ndc_step, Form("= %.2f", fitBin[i][j]->GetChisquare()/fitBin[i][j]->GetNDF()));

    		latex.SetTextColor(lcolor[3]);
    		latex.DrawLatexNDC(x_ndc, y_ndc_start, Form("M_{#rho^{0}} = %.3f (GeV/c^{2})", fitBin[i][j]->GetParameter(10)));
    		latex.DrawLatexNDC(x_ndc, y_ndc_start - y_ndc_step, Form("W_{#rho^{0}} = %.3f (GeV/c^{2})", fitBin[i][j]->GetParameter(12)));
    		latex.DrawLatexNDC(x_ndc, y_ndc_start - 2*y_ndc_step, Form("A_{#rho^{0}} = %.3f", fitBin[i][j]->GetParameter(11)));

			latex.SetTextColor(lcolor[1]);
			latex.DrawLatexNDC(x_ndc, y_ndc_start - 3*y_ndc_step, Form("M_{f_{0}} = %.3f (GeV/c^{2})", fitBin[i][j]->GetParameter(4)));
			latex.DrawLatexNDC(x_ndc, y_ndc_start - 4*y_ndc_step, Form("W_{f_{0}} = %.3f (GeV/c^{2})", fitBin[i][j]->GetParameter(6)));
			latex.DrawLatexNDC(x_ndc, y_ndc_start - 5*y_ndc_step, Form("A_{f_{0}} = %.3f", fitBin[i][j]->GetParameter(5)));

			latex.SetTextColor(lcolor[2]);
			latex.DrawLatexNDC(x_ndc, y_ndc_start - 6*y_ndc_step, Form("M_{f_{2}} = %.3f (GeV/c^{2})", fitBin[i][j]->GetParameter(7)));
			latex.DrawLatexNDC(x_ndc, y_ndc_start - 7*y_ndc_step, Form("W_{f_{2}} = %.3f (GeV/c^{2})", fitBin[i][j]->GetParameter(9)));
			latex.DrawLatexNDC(x_ndc, y_ndc_start - 8*y_ndc_step, Form("A_{f_{2}} = %.3f", fitBin[i][j]->GetParameter(8)));

			latex.SetTextColor(lcolor[0]);
			latex.DrawLatexNDC(x_ndc, y_ndc_start - 9*y_ndc_step, Form("B = %.3f", fitBin[i][j]->GetParameter(0)));
			latex.DrawLatexNDC(x_ndc, y_ndc_start - 10*y_ndc_step, Form("n = %.3f", fitBin[i][j]->GetParameter(1)));
			latex.DrawLatexNDC(x_ndc, y_ndc_start - 11*y_ndc_step, Form("c_{1} = %.3f", fitBin[i][j]->GetParameter(2)));
			latex.DrawLatexNDC(x_ndc, y_ndc_start - 12*y_ndc_step, Form("c_{2} = %.3f", fitBin[i][j]->GetParameter(3)));

            //  reset Amplitude
            fitBin[i][j]->SetParameter(0, 0);
			fitBin[i][j]->SetParameter(8, 0);
			fitBin[i][j]->SetParameter(11, 0);

			// Intgr_RawY_cntl[i][j] = fitBin[i][j]->Integral(
			// 	fitBin[i][j]->GetParameter(4) - fitBin[i][j]->GetParameter(6)*nSigIntgr,
			// 	fitBin[i][j]->GetParameter(4) + fitBin[i][j]->GetParameter(6)*nSigIntgr );
			// Intgr_RawY_cntl[i][j] /= hMassSig[i][j]->GetBinWidth(1);

			// Intgr_RawY_stat[i][j] = fitBin[i][j]->IntegralError(
			// 	fitBin[i][j]->GetParameter(4) - fitBin[i][j]->GetParameter(6)*nSigIntgr,
			// 	fitBin[i][j]->GetParameter(4) + fitBin[i][j]->GetParameter(6)*nSigIntgr,
			// 	fitBin[i][j]->GetParameters(), fitResults[i][j]->GetCovarianceMatrix().GetMatrixArray(), 1e-6 );
			// Intgr_RawY_stat[i][j] /= hMassSig[i][j]->GetBinWidth(1);
        }
		//  file name
		std::time_t t = std::time(nullptr);
		char dateStr[100];
		std::strftime(dateStr, sizeof(dateStr), "%Y%m%d", std::localtime(&t));
		TString fileName = Form("/Users/seul/seul_workspace/f0ana/O2P_ana/fit/plot/%s_fitout_mult_%d_%d.pdf", dateStr, mult_min[i], mult_max[i]);
		// fileName << "/Users/seul/seul_workspace/f0ana/O2P_ana/fit/plot" << dateStr << "_fitout_mult_%d_%d.pdf" << mult_min[i], mult_max[i];
        // TString fileName = Form("plot/fitout_mult_%d_%d.pdf", mult_min[i], mult_max[i]);
		c->SaveAs(fileName);
		delete c;

        // hPtStat[i] = new TH1D(Form("hPtStat_%d",i),"",npt-1,PtBinnings);
        // for (j = 0; j < npt-1; j++) {
        //     hPtStat[i]->SetBinContent(j+1, Intgr_RawY_cntl[i][j]);
		// 	hPtStat[i]->SetBinError(j+1, Intgr_RawY_stat[i][j]);
        // }
		// break;
    }

    // TFile* Fout = new TFile(Form("/Users/seul/seul_workspace/f0ana/results/F0IntgrOut.root"), "recreate");
	// for (i = 0; i < nmult; i++) {
    //     hPtStat[i]->Write();
	// }
    
    // std::time_t t = std::time(nullptr);
    // char dateStr[100];
    // std::strftime(dateStr, sizeof(dateStr), "%Y%m%d", std::localtime(&t));
    // std::stringstream fileName;
    // fileName << "/Users/seul/seul_workspace/f0ana/results/" << dateStr << "_Fit.root";
	// TFile* Ffitout = new TFile(fileName.str().c_str(), "recreate");

	// for (int i = 0; i < nmult; i++) {
	// 	gf0Width_init[i]->Write(Form("gf0Width_init_mult%d", i));
	// 	gf0Width[i]->Write(Form("gf0Width_mult%d", i));
    //     gf0Mass[i]->Write(Form("gf0Mass_mult%d", i));
    //     gf0Chi2perNDF[i]->Write(Form("gf0Chi2perNDF_mult%d", i));
    	
	// 	for (int j = 0; j < npt; j++) {
    //     	fitBin[i][j]->Write(Form("fitBin_%d_%d", i, j));
	// 		for (int k = 0; k < 4; k++) {
	// 			fitBinComp[i][j][k]->Write(Form("fitBinComp_%d_%d_%d", i, j, k));
	// 		}
	// 	}
	// }
	// Ffitout->Close();
}
