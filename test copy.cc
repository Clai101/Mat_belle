#include<test.h>

std::string rm_sub_str(std::string str, std::string substr) {
    size_t start_position_to_erase = str.find(substr);
    if (start_position_to_erase != std::string::npos) {
        str.erase(start_position_to_erase, substr.length());
    }
    return str;
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
  std::cout << "_____________________________________________________-" << newl;
  std::cout << argc << newl;
  std::cout << argv[0] << newl;
  std::cout << argv[1] << newl;
  std::cout << argv[2] << newl;



  double down = 1.5;
  double up = 3.5;
  std::string _name = "noname";
  std::string cut = "";
  std::vector<int> chu;
  bool chub= true;

  for (auto iter {argv.begin()}; iter != argv.end(); ++iter){
    iter = iter + ("__").c_str();
    if (iter.find("__down = ") != std::string::npos){
      str::string tool = rm_sub_str(iter, "__down = ");
      down = str::stof(tool); 
    }
    if (iter.find("__up = ") != std::string::npos){
      str::string tool = rm_sub_str(iter, "__up = ");
      up = str::stof(tool); 
    }
    if (iter.find("__fname = ") != std::string::npos){
      str::string tool = rm_sub_str(iter, "__fname = ");
      _name = tool.c_str(); 
    }
    if (iter.find("__cut = ") != std::string::npos){
      str::string tool = rm_sub_str(iter, "__cut = ");
      cut = tool.c_str(); 
    }
    if (iter.find("__chu = ") != std::string::npos){
      str::string tool = rm_sub_str(iter, "__chu = ");
      chu = parseString(tool).c_str(); 
      chub = false;
    }
  }

  std::string path_name = fs::current_path();
  std::string fold = "root";
  path_name = path_name + "/" + fold;
  std::vector<std::string> names;
  std::string ext = ".root";

  for (auto &entry : fs::directory_iterator(path_name)){
    names.push_back(entry.path().filename().c_str()); 
  }
  
  //for (int i = 1; i <= 12 ; i++){


	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	TH1F *hist = new TH1F("hist", "", 100, down, up);
	
  if (chub){
    
	for (auto iter {names.begin()}; iter != names.end(); ++iter) {
  }
    std::string  fname = fold + "/" + (*iter).c_str();
    if (fname.find(ext) != std::string::npos){
      TFile *input = new TFile(fname.c_str(), "read");
      TTree *tree = (TTree*)input->Get("h1");
      
      TH1F *temp = new TH1F("temp", "", 100, down, up);

      tree->Draw("sqrt(rm2l) >> temp", "chl <= 4 && rm2l > 0 && abs(mach) < 0.015 && abs(ml - 2.28646) < 0.015" /*+ std::to_string(i)).c_str()*/);

      hist->Add(temp, 1); // добавляем данные temp с коэффициентом 1
    }
	}
	
	//hist->GetYaxis()->SetRangeUser(0, 100);
  hist->Draw("[]");
  c1->Print(("./results/" + _name + "/" + _name + ".root").c_str());
  
  
  return 0;
}

