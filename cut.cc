#include<test.h>

std::string rm_sub_str(std::string str, std::string substr) {
    size_t start_position_to_erase = str.find(substr);
    if (start_position_to_erase != std::string::npos) {
        str.erase(start_position_to_erase, substr.length());
    }
    return str;
}

void FindMinMaxFromTree(TTree* tree, const std::string& varName, float& down, float& up) {

    TCanvas* c1 = new TCanvas("c1", "c1", 1600, 1200);
    TH1F* hist = new TH1F("hist", "", 1000, -100, 100);
    tree->Draw((varName + " >> hist").c_str());
    TH1F* temp = (TH1F*)gDirectory->Get("hist"); // Retrieve the histogram 'hist'

    std::vector<int> binsWithMoreThanOneEvent;
    for (int bin = 1; bin <= 1000; ++bin) {
        double binContent = temp->GetBinContent(bin);
        if (binContent > 1) {
            down = temp->GetBinCenter(bin);
            break;
        }
    }
    for (int bin = 1000; bin > 1; --bin) {
        double binContent = temp->GetBinContent(bin);
        if (binContent > 1) {
            up = temp->GetBinCenter(bin);
            break;
        }
    }
}

std::vector<int> parseString(const std::string& input) {
    std::vector<int> result;
    std::istringstream iss(input);

    std::string segment;
    while (std::getline(iss, segment, ',')) {
        size_t dash_pos = segment.find('-');
        if (dash_pos == std::string::npos) {
            int number;
            std::istringstream(segment) >> number;
            result.push_back(number);
        } else {
            int start, end;
            std::istringstream(segment.substr(0, dash_pos)) >> start;
            std::istringstream(segment.substr(dash_pos + 1)) >> end;
            for (int i = start; i <= end; ++i) {
                result.push_back(i);
            }
        }
    }

    return result;
}


int main(int argc, char *argv[]) {
  std::cout << "_____________________________________________________" << newl;


  std::string _name = "noname";
  std::string out;
  std::string cut = " ";
  std::vector<int> chu;
  bool chub = false;
  int nbins = 100;
  float down, up;

  for (int i = 1; i < argc; ++i){
    
    std::string iter = "__" + std::string(argv[i]);
    std::cout << "iter = " << iter << newl;

    if (iter.find("__down = ") != std::string::npos) {
      down = std::stof(rm_sub_str(iter, "__down = ")); 
      std::cout << "down = " << down << newl;}
    if (iter.find("__up = ") != std::string::npos){
      up = std::stof(rm_sub_str(iter, "__up = ")); 
      std::cout << "up = " << up << newl;}
    if (iter.find("__out = ") != std::string::npos){
      out = rm_sub_str(iter, "__out = ").c_str();
      std::cout << "out: " << out << newl;}
    if (iter.find("__fname = ") != std::string::npos){
      _name = rm_sub_str(iter, "__fname = ").c_str(); 
      std::cout << "name: " << _name << newl;}
    if (iter.find("__cut = ") != std::string::npos){
      cut = rm_sub_str(iter, "__cut = ").c_str();
      std::cout << "cut: " << cut << newl;}
    if (iter.find("__chu = ") != std::string::npos){
      chu = parseString(rm_sub_str(iter, "__chu = ")); 
      chub = true;}
    if (iter.find("__nbins = ") != std::string::npos){
      nbins = std::stoi(rm_sub_str(iter, "__nbins = "));}
  }


  std::cout << "____________________________" << newl;
  std::string path_name = fs::current_path();
  std::string fold = "root";
  path_name = path_name + "/" + fold;
  std::vector<std::string> names;
  std::string ext = ".root";

  for (auto &entry : fs::directory_iterator(path_name))
    names.push_back(entry.path().filename().c_str()); 

  for (auto iter {names.begin()}; iter != names.end(); ++iter) {
    std::string  fname = fold + "/" + (*iter).c_str();
    if (fname.find(ext) != std::string::npos){
      TFile *input = new TFile(fname.c_str(), "read");
      TTree* tree = (TTree*)input->Get("h1");
      FindMinMaxFromTree(tree, out, down, up);
      input->Close();

      break;
  }}
  std::cout << "up, down: " << up << "   " << down << newl;

	
  if (chub){
  TCanvas *c2 = new TCanvas("c2", "c2", 1600, 1600);
  int sqr_l = ceil(sqrt(chu.size()));
  c2->Divide(sqr_l, sqr_l);
  int j = 1;
  for (auto i {chu.begin()}; i != chu.end(); ++i, j++){
    TCanvas *c1 = new TCanvas("c1", "c1", 1600, 1600);
    TH1F *hist = new TH1F("hist", "", nbins, down, up);
    TH1F *hist2 = new TH1F("hist", "", nbins, down, up);
    for (auto iter = names.begin(); iter != names.end(); ++iter) {
      std::string  fname = fold + "/" + (*iter).c_str();
      if (fname.find(ext) != std::string::npos){
        TFile *input = new TFile(fname.c_str(), "read");
        TTree *tree = (TTree*)input->Get("h1");
        TH1F *temp = new TH1F("temp", "", nbins, down, up);
        tree->Draw((out + " >> temp").c_str(), (cut + " && chu == " + std::to_string(*i)).c_str());
        hist->Add(temp, 1);
        input->Close();
      }}
      hist->Draw();
      c2->cd(j);
      hist2->Add(hist, 1);
      hist2->Draw();
      
      mkdir(("./results/" + _name).c_str(), 0777);
      mkdir(("./results/" + _name + "/" + "root").c_str(), 0777);
      mkdir(("./results/" + _name + "/" + "pdf").c_str(), 0777);
      c1->Print(("./results/" + _name + "/" + "root" + "/" + _name + "_chu" + std::to_string(*i) + ".root").c_str());
      c1->Print(("./results/" + _name + "/" + "pdf" + "/" + _name + "_chu" + std::to_string(*i) + ".pdf").c_str());
    }
    c2->Print(("./results/" + _name + "/" + _name + "_chu" + ".pdf").c_str());
    c2->Print(("./results/" + _name + "/" + _name + "_chu" + ".root").c_str());

  }
  else{
    TCanvas *c1 = new TCanvas("c1", "c1", 1600, 1200);
  	TH1F *hist = new TH1F("hist", "", nbins, down, up);
    for (auto iter {names.begin()}; iter != names.end(); ++iter) {
      std::string  fname = fold + "/" + (*iter).c_str();
      if (fname.find(ext) != std::string::npos){
        TFile *input = new TFile(fname.c_str(), "read");
        TTree *tree = (TTree*)input->Get("h1");
        std::cout << "up, down: " << up << "   " << down << newl;
        TH1F *temp = new TH1F("temp", "", nbins, down, up);
        tree->Draw((out + " >> temp").c_str(), (cut).c_str());
        hist->Add(temp, 1);
        input->Close();
    }}
    hist->Draw();
    c1->Print(("./results/" + _name + ".png").c_str());
    c1->Print(("./results/" + _name + ".root").c_str());
  }
  
  // ""

  return 0;
}

