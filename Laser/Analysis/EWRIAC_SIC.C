{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Apr  8 13:33:49 2013) by ROOT version5.32/00
   TCanvas *c1 = new TCanvas("c1", "c1",246,44,700,500);
   gStyle->SetOptStat(0);
   c1->Range(-2375,-0.25,11375,2.25);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   TH1D *haxis = new TH1D("haxis","haxis",2200,-1000,10000);
   haxis->SetMinimum(0);
   haxis->SetMaximum(2);
   haxis->SetStats(0);
   haxis->GetXaxis()->SetTitle("Dose rate (rad/h)");
   haxis->GetYaxis()->SetTitle("EWRIAC (m^{-1})");
   haxis->Draw("axis");
   
   TGraph *graph = new TGraph(3);
   graph->SetName("EWRIAC_SIC");
   graph->SetTitle("EWRIAC_SIC");
   graph->SetFillColor(1);
   graph->SetMarkerStyle(20);
   graph->SetPoint(0,15,0.16);
   graph->SetPoint(1,400,0.69);
   graph->SetPoint(2,9000,1.43);
   
   TH1F *Graph_EWRIAC_SIC1 = new TH1F("Graph_EWRIAC_SIC1","EWRIAC_SIC",100,0,9898.5);
   Graph_EWRIAC_SIC1->SetMinimum(0.033);
   Graph_EWRIAC_SIC1->SetMaximum(1.557);
   Graph_EWRIAC_SIC1->SetDirectory(0);
   Graph_EWRIAC_SIC1->SetStats(0);
   graph->SetHistogram(Graph_EWRIAC_SIC1);
   
   
   TF1 *fit_log2 = new TF1("fit_log2","[0]*log([1]*x+1.)",-10,10000);
   fit_log2->SetFillColor(19);
   fit_log2->SetFillStyle(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#0000ff");
   fit_log2->SetLineColor(ci);
   fit_log2->SetLineWidth(3);
   fit_log2->SetChisquare(0.001280728);
   fit_log2->SetNDF(1);
   fit_log2->SetParameter(0,0.2320875);
   fit_log2->SetParError(0,0.01513256);
   fit_log2->SetParLimits(0,0,0);
   fit_log2->SetParameter(1,0.05090391);
   fit_log2->SetParError(1,0.0158747);
   fit_log2->SetParLimits(1,0,0);
   graph->GetListOfFunctions()->Add(fit_log2);
   graph->Draw("p");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
