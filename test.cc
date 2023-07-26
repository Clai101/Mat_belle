#include<test.h>

int main(int argc, char *argv[]) {
  std::string _name;
  std::cout << "_____________________________________________________-" << newl;
  std::cout << argc << newl;
  std::cout << argv[0] << newl;
  std::cout << argv[1] << newl;
  std::cout << argv[2] << newl;

  double down = 1.5;
  double up = 3.5;

  if(argc >= 1){
    _name = argv[1];
  }
  else{
    _name = argv[0];
  }

  std::string path_name = fs::current_path();
  std::string fold = "root";
  path_name = path_name + "/" + fold;
  std::vector<std::string> names;
  std::string ext = ".root";

  for (auto &entry : fs::directory_iterator(path_name)){
    names.push_back(entry.path().filename().c_str()); 
  }
  
  for (int i = 1; i <= 7 ; i++){


	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	TH1F *hist = new TH1F("hist", "", 100, down, up);
	
	for (auto iter {names.begin()}; iter != names.end(); ++iter) {
    std::string  fname = fold + "/" + (*iter).c_str();
    if (fname.find(ext) != std::string::npos){
      TFile *input = new TFile(fname.c_str(), "read");
      TTree *tree = (TTree*)input->Get("h1");
      
      TH1F *temp = new TH1F("temp", "", 100, down, up);
        if(argc >= 3){
          std::string _cut = argv[2];
          tree->Draw("en >> temp", _cut.c_str());
        }
       if(argc >= 4){
          std::string _cut = argv[2];
          std::string _ax = argv[3];
          tree->Draw((_ax + " >> temp").c_str(), _cut.c_str());
        }
        else{
          //tree->Draw("ml >> temp", " chrgu == 0 &&  chl == 3");
          tree->Draw("sqrt(rm2l) >> temp", ("chrgu == 0 &&  chl == 5 && rm2l > 0 && abs(mach) < 0.015 && abs(rm2n) < 0.1 && abs(ml - 2.28646) < 0.015 && chu ==" + std::to_string(i)).c_str());
          //tree->Draw("chu >> temp", "chrgl == 0");
          //tree->Draw("ml >> temp", "p < 0.8 && (ntr == 0 ||  ntr == 1) && chг >= 11 && chг <= 14" );
          //tree->Draw("en >> temp", "ntr == 0 && abs(m_sigm - 2.45397) < 0.015  && abs(m_lamc - 2.286) < 0.015 && p < 0.8");
          //tree->Draw("en >> temp", "p < 0.08 && en > 8 && chu == 5 && abs(ml - 2.28646) < 0.015 && abs(mach - 1.86484) < 0.015 && ntr == 0");
          //tree->Draw("en - ecm>> temp", "p < 0.5 && en > 8 && chu == 7 && abs(ml - 2.28646) < 0.015 && abs(mach - 1.86966) < 0.015 && ntr == 0");
          //tree->Draw("abs(mrec2_r) >> temp", " ntr == 0");
          //tree->Draw("en >> temp", "ntr == 0 && ch_u == 15 ");
        }
      hist->Add(temp, 1); // добавляем данные temp с коэффициентом 1
    }
	}
	
	//hist->GetYaxis()->SetRangeUser(0, 100);
  hist->Draw("[]");
  c1->Print(("./results/" + _name + "_chu" + std::to_string(i) + ".root").c_str());
  }
  
  return 0;
}

