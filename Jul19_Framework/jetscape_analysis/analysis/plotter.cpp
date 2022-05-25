//run with root -l "plotter.cpp(\"PP.root\",\"PbPb.root\",\"plots.pdf\")"
int plotter(std::string pp_histFileName, std::string pbpb_histFileName , std::string plotFileName)
{
	std::cout << "Reading from " << pp_histFileName << " and " << pbpb_histFileName <<
		" and writing to " << plotFileName << std::endl;
	TFile * pp  = TFile :: Open ( pp_histFileName . c_str () ," READ ");
	TFile * pbpb  = TFile :: Open ( pbpb_histFileName . c_str () ," READ ");
	TH1F *pp_JetPt = (TH1F*) pp->Get("hJetPt_R0.4");
	TH1F *pbpb_JetPt = (TH1F*)pbpb->Get("hJetPt_R0.4");
	TH1F *pp_JetMultiplicity = (TH1F*) pp->Get("hNJet_R0.4");
	TH1F *pbpb_JetMultiplicity =(TH1F*) pbpb->Get("hNJet_R0.4");
	
	TH1F *pp_LepPt = (TH1F*) pp->Get("hLepPt");
	TH1F *pbpb_LepPt = (TH1F*)pbpb->Get("hLepPt");
	
	TCanvas canvas("canvas");
	
	// Jet multiplicity
	canvas.cd();
	pp_JetMultiplicity->SetLineColor(kBlack);
	pbpb_JetMultiplicity->SetLineColor(kRed);
	TLegend legend(0.7,0.6,0.85,0.75);
	legend.AddEntry(pp_JetMultiplicity,"PP");
	legend.AddEntry(pbpb_JetMultiplicity,"PbPb");
	pp_JetMultiplicity->GetXaxis()->SetTitle("Number of Jets");
	pp_JetMultiplicity->Draw();
	pbpb_JetMultiplicity->Draw("same");
	legend.Draw("same");
	canvas.Print("nJet.pdf");	
	canvas.Clear();
	
	//Jet Pt distribution
	canvas.cd();
	canvas.SetLogy();
	pp_JetPt->SetLineColor(kBlack);
	pbpb_JetPt->SetLineColor(kRed);
	pp_JetPt->GetXaxis()->SetTitle("p_{T}^{jet}(GeV)");
	pp_JetPt->GetYaxis()->SetTitle("#frac{1}{N_{ev}} #frac{dN}{dp^{jet}_{p_{T}}}");
	pp_JetPt->Draw();
	pbpb_JetPt->Draw("same");
	legend.Clear();
	legend.AddEntry(pp_JetPt,"PP");
	legend.AddEntry(pbpb_JetPt,"PbPb");
	legend.Draw("same");
	canvas.Print("jetPt.pdf");
	canvas.Clear();
	
	
	// Jet R_AA
	canvas.cd();
	canvas.SetLogy(false);
	TLine *line = new TLine(80.0, 1.0, 150.0, 1.0);
	line->SetLineColor(kBlue);
	line->SetLineStyle(2);
	//pp_JetPt->Rebin(2);
	//pbpb_JetPt->Rebin(2);
	pp_JetPt->Sumw2();
	pbpb_JetPt->Sumw2();
	TH1F *rAA = new TH1F(*pbpb_JetPt);
	rAA->Divide(pp_JetPt);
	rAA->GetXaxis()->SetTitle("p_{T}^{jet}(GeV)");
	rAA->GetYaxis()->SetTitle("R_{AA}");
	rAA->SetLineColor(kBlack);
	rAA->GetYaxis()->SetRangeUser(0.,1.5);
	rAA->GetXaxis()->SetRangeUser(80.,150.);
	rAA->SetStats(0);
	rAA->Draw("pe");
	line->Draw("same");
	canvas.Print("rAA.pdf");
	canvas.Clear();
	
	
	
	//Lepton Pt distribution
	canvas.cd();
	//canvas.SetLogy();
	pp_LepPt->SetLineColor(kBlack);
	pbpb_LepPt->SetLineColor(kRed);
	pp_LepPt->GetXaxis()->SetTitle("p_{T}^{lep}(GeV)");
	pp_LepPt->GetYaxis()->SetTitle("#frac{1}{N_{ev}} #frac{dN}{dp^{lep}_{p_{T}}}");
	pp_LepPt->Draw();
	pbpb_LepPt->Draw("same");
	legend.Clear();
	legend.AddEntry(pp_JetPt,"PP");
	legend.AddEntry(pbpb_JetPt,"PbPb");
	legend.Draw("same");
	canvas.Print("LepPt.pdf");
	canvas.Clear();
	
	
	// Lepton R_AA
	canvas.cd();
	canvas.SetLogy(false);
	TLine *line2 = new TLine(0.0, 1.0, 50.0, 1.0);
	line2->SetLineColor(kBlue);
	line2->SetLineStyle(2);
	//pp_JetPt->Rebin(2);
	//pbpb_JetPt->Rebin(2);
	pp_LepPt->Sumw2();
	pbpb_LepPt->Sumw2();
	TH1F *Lep_rAA = new TH1F(*pbpb_LepPt);
	Lep_rAA->Divide(pp_LepPt);
	Lep_rAA->GetXaxis()->SetTitle("p_{T}^{lep}(GeV)");
	Lep_rAA->GetYaxis()->SetTitle("R^{lep}_{AA}");
	Lep_rAA->SetLineColor(kBlack);
	Lep_rAA->GetYaxis()->SetRangeUser(-2.0,2.0);
	Lep_rAA->SetStats(0);
	Lep_rAA->Draw("pe");
	line2->Draw("same");
	canvas.Print("Lep_rAA.pdf");
	canvas.Clear();
	
	return 0;
}
