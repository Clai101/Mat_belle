bool 
checkMother(Particle &mother, std::vector<Particle> &childern){
  const Gen_hepevt moth =  mother.pType().genHepevt();
  bool same_mother = true;
  for  (std::vector<Particle>::iterator chl = childern.begin(); chl!=childern.end(); ++chl){
    const Gen_hepevt chl_hep = chl.genHepevt();
    for  (std::vector<Particle>::iterator chl__ = chl+1; chl__!=childern.end(); ++chl__){
      const Gen_hepevt chl_hep__ = chl__.genHepevt();
      if (chl_hep == chl_hep__) return false;
    }
    same_mother = same_mother and searchMother(moth, chl_hep);
    if (not same_mother)
      return same_mother;
  }
  return same_mother;
}

bool
searchMother(const Gen_hepevt p, const Gen_hepevt sig){
  if (!p.mother()) return 0;
  else if (p.mother() == sig) return 1;
  else return checkMother(p.mother(), sig);
}


vector<Particle>
findchild (Particle &part){
  vector<Particle> childern;
  if (part.nChildren() == 0) return vector<Particle> vect{part};
  for (int i = 0; i < part.nChildren(), i++){
    vector<Particle> b = findchild(part.child(i));
    childern.insert(childern.end(), b.begin(), b.end());
  }
  return childern;
}

bool 
checkReal (Particle &part){
  vector<Particle> childern = findchild(part);
  return checkMother(part, childern);
}