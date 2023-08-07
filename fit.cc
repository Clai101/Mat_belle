#include<test.h>


int n_bins;
double min_val;
double max_val;

double logLikelihood;


const int n_par = 22;


std::vector<double> result;
std::vector<double> result1;
std::vector<double> data1;


std::vector<std::string> split(const std::string &s) {
  std::stringstream ss(s);
  std::string item;
  std::vector<std::string> res;
  while (std::getline(ss, item, ' ')) {
    res.push_back(item);
  }
  return res;
}


double CHE(double x, double *par) {
  double y, T_0, T_1, T_2, T_3, T_4, T_5, T_6;
  
  //double sqrt2pi = 2.50662827;
      
  y=(x-(min_val+max_val)/2)/((max_val - min_val)/2);
  T_0 = 1;
  T_1 = y;
  T_2 = 2*y*y - 1;
  T_3 = 4*y*y*y-3*y;
  T_4 = 8*y*y*y*y - 8*y*y + 1;
  T_5 = 16*y*y*y*y*y - 2*y*y*y + 5*y;
  T_6 = 32*y*y*y*y*y*y - 48*y*y*y*y + 18*y*y - 1;

  return (T_0 + par[15]*T_1 + par[16]*T_2 + par[17]*T_3);// + par[15]*T_4);
}



double g_g_g(double x, double *par) {

  double sqrt2pi = 2.50662827;

  double mean1   = par[4] + par[20];
  double mean2   = mean1 + par[7]*par[21];
  double mean_cb = mean1 + par[10]*par[21];
  double mean_4 = mean1 + par[7]*par[21];
  
  double sigmaG1 = par[5]*par[21];
  double sigmaG2 = par[8]*par[21];
  double sigmaCB = par[11]*par[21];
  double sigma2 = par[11]*par[21];
  
  double g1 = exp( -(x-mean1)*(x-mean1) / (2*sigmaG1*sigmaG1));
  g1 /= (sqrt2pi * sigmaG1);
  
  double g2 = exp( -(x-mean2)*(x-mean2) / (2*sigmaG2*sigmaG2));
  g2 /= (sqrt2pi * sigmaG2);
  
  double g3 = exp( -(x-mean_cb)*(x-mean_cb) / (2*sigmaCB*sigmaCB));
  g3 /= (sqrt2pi * sigmaCB);
  
  double g4 = exp( -(x-mean_4)*(x-mean_4) / (2*sigma2*sigma2));
  g4 /= (sqrt2pi * sigma2);


  return g1;
}


void fill_result(double *par) {
  double bin = (max_val - min_val) / n_bins;

  std::vector<double> signal(n_bins);
  std::vector<double> bg(n_bins);
  
  double sum_sig = 0;
  for (int j = 0; j<n_bins; j++) {
    double x = min_val + (j + 0.5) * bin;
    signal[j] = g_g_g(x,par);
    sum_sig += signal[j];
  }
  for (int j = 0; j<n_bins; j++)
    signal[j] /= sum_sig;
    
  for (int j = 0; j<n_bins; j++) {
    double x = min_val + (j + 0.5) * bin;
    result[j] =  par[0]*signal[j] + par[12]*CHE(x,par);
  }
}


void fcn1(int &npar, double *grad, double &fval, double *par, int iflag) {

  fill_result(par);

  fval = 0;
  for (int j = 0; j<n_bins; j++) {
    double mu = result[j];
    double nn = data1[j];
    if (mu > 0) {
      if (nn > 0)
	fval += 2.0*(mu-nn+nn*log(fabs(nn/mu)));
      else
	fval += 2.0*mu;
    }
    else fval += 1000*mu*mu;
    if (isnan(fval))
      printf("Fval is NaN\n");
  }
  logLikelihood = fval;
}


int main() {

  
  std::string  fname = "results/sqr_rm2l.root";

  TFile miofile(fname.c_str(),"read");
  TH1F *hist1;
  TCanvas *c2 = (TCanvas*)miofile.Get("c1");
  hist1 = (TH1F*)c2->GetPrimitive("hist");
  
  n_bins = hist1->GetNbinsX();
  min_val = hist1->GetBinLowEdge(1);
  max_val = hist1->GetBinLowEdge(n_bins+1);  

  TH1F *hist = new TH1F("hist", "", n_bins, min_val, max_val);

  result.resize(n_bins);
  result1.resize(n_bins);
  data1.resize(n_bins);

  TCanvas *c11 = new TCanvas("c11", "c11",72,101,700,500);
  c11->Range(1.897901,-48.82378,2.026446,256.9672);
  c11->SetFillColor(0);
  c11->SetBorderMode(0);
  c11->SetBorderSize(2);
  c11->SetLeftMargin(0.1719198);
  c11->SetRightMargin(0.05014326);
  c11->SetTopMargin(0.03991597);
  c11->SetBottomMargin(0.1596639);
  c11->SetFrameBorderMode(0);
  c11->SetFrameBorderMode(0);


  double par[n_par], err[n_par];
  std::vector<std::string> name;
  std::ofstream file_out("resol_fit_output.txt");	
	
  for (int j = 0; j<n_bins;j++) {
    data1[j] = hist1->GetBinContent(j+1);
  }
	
  TMinuit mn(n_par);
  mn.SetFCN(fcn1);
  int flag = 0;
	
  std::string s;
  std::ifstream file_in("resol_fit_input.txt");
  int k = 0;
	
  while (getline(file_in,s)) {
    std::vector<std::string> res;
    res = split(s);
    name.push_back(res[0]);		
    mn.mnparm(k, res[0], atof(res[1].c_str()), atof(res[2].c_str()), 0, 0, flag);
    k += 1;
  }
	
  //mn.mncomd("FIX 2 3 4 5 6 7 8 9 10 11 12", flag);
  //mn.mncomd("FIX 21 22", flag);
  mn.mncomd("MIGRAD", flag);
  mn.mncomd("MINOS 10000 1",flag);

	
  for (int i = 0; i < n_par; i++) {
    mn.GetParameter(i, par[i], err[i]);
    file_out << name[i] << " " << par[i] << " " << err[i] << "\n"; 
  }
	
  for (int i = 0; i < n_bins; ++i) {
    hist->SetBinContent(i+1, result[i]);
  }
  c11->cd();
  hist1->SetMinimum(0);	
  hist1->SetLabelSize(0.06f);
  hist1->GetYaxis()->SetLabelSize(0.06f);
  hist1->GetXaxis()->SetNdivisions(6,5,0);
  hist1->GetYaxis()->SetNdivisions(6,5,0);
  hist1->SetXTitle("M_{rec \\Lambda} (GeV/c^{2})");
  hist1->SetYTitle("Events / 20 MeV/c");
  hist1->GetXaxis()->SetTitleSize(0.06f);
  hist1->GetYaxis()->SetTitleSize(0.06f);
  hist1->SetMarkerColor(kViolet+4);
  hist1->SetMarkerStyle(20);
  hist1->SetMarkerSize(0.6);
  //hist1->SetStats(kFALSE);
  //hist->SetStats(kFALSE);

  hist1->DrawCopy("e");//->Rebin(4)->SetStats(kFALSE);

  hist->SetLineColor(kRed);


  hist->DrawCopy("same");//->Rebin(4)->SetStats(kFALSE);



  printf("fval is %f\n",logLikelihood);
  c11->SaveAs("./results/RM.png");
  return 0;
}
