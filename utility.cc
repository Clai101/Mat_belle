//
// $Id: utility.cc 10613 2008-09-03 11:36:24Z katayama $
//
#include <string>
#include <iostream>
#include "belle.h"
#include "particle/utility.h"
#include "particle/gammac.h"
#include "toolbox/FoxWolfr.h"
#include "toolbox/FuncPtr.h"
#include "toolbox/Thrust.h"
#include "mdst/mdst.h"
#include "kfitter/kmakemother.h"
#include "kid/atc_pid.h"
#include "eid/eid.h"
#include "mdst/Muid_mdst.h"
#include MDST_H
#include HEPEVT_H
#include <algorithm>
#include "belleCLHEP/Matrix/Vector.h"
#include "belleCLHEP/Matrix/SymMatrix.h"
#include "belleCLHEP/Matrix/Matrix.h"
#include "belleCLHEP/Vector/ThreeVector.h"
#include "belleCLHEP/Vector/LorentzVector.h"
#include "belleCLHEP/Geometry/Point3D.h"
#include "tables/evtcls.h"
#include "tables/evtvtx.h"
#include "tables/belletdf.h"
#include "helix/Helix.h"
#include "ip/IpProfile.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


void 
makeKPi(std::vector<Particle> &k_p, 
	std::vector<Particle> &k_m, 
	std::vector<Particle> &pi_p, 
	std::vector<Particle> &pi_m,
	const int flag){
  Mdst_charged_Manager &charged_mag = Mdst_charged_Manager::get_manager();
  
  //...Particle Type
  Ptype ptype_pion_plus("PI+");
  Ptype ptype_pion_minus("PI-");
  Ptype ptype_kaon_plus("K+");
  Ptype ptype_kaon_minus("K-");
 
  //...Fills pion and kaon lists with MDST_Charged Data Base
  for(std::vector<Mdst_charged>::iterator i = charged_mag.begin();
      i != charged_mag.end(); ++i){
#if defined(BELLE_DEBUG)
    try {
#endif
    if(flag && !good_charged(*i))continue;
    if((*i).charge() > 0.){
      Particle tmp1(*i, ptype_kaon_plus);
      k_p.push_back(tmp1);
      Particle tmp2(*i, ptype_pion_plus);
      pi_p.push_back(tmp2);
    }else{
      Particle tmp1(*i, ptype_kaon_minus);
      k_m.push_back(tmp1);
      Particle tmp2(*i, ptype_pion_minus);
      pi_m.push_back(tmp2);
    }
#if defined(BELLE_DEBUG)
    }
    catch(const std::string &e) {
      dout(Debugout::ERR,"utility") << "particle/functions/utility::makeKPi:" << e << std::endl;
	Mdst_charged_Manager::get_manager().dump();
	Mdst_trk_Manager::get_manager().dump();
      continue;
    }
#endif    
  }
}

void 
makeLepton(std::vector<Particle> &e_p,
	   std::vector<Particle> &e_m,
	   std::vector<Particle> &mu_p,
	   std::vector<Particle> &mu_m,
	   const int flag){
  Mdst_charged_Manager &charged_mag = Mdst_charged_Manager::get_manager();
  
  //...Particle Type
  Ptype ptype_elec_plus("E+");
  Ptype ptype_elec_minus("E-");
  Ptype ptype_muon_plus("MU+");
  Ptype ptype_muon_minus("MU-");
 
  //...Fills elec and muon lists with MDST_Charged Data Base
  for(std::vector<Mdst_charged>::iterator i = charged_mag.begin();
      i != charged_mag.end(); ++i){
#if defined(BELLE_DEBUG)
    try {
#endif
    if(flag && !good_charged(*i))continue;
    if((*i).charge() > 0.){
      Particle tmp1(*i, ptype_elec_plus);
      e_p.push_back(tmp1);
      Particle tmp2(*i, ptype_muon_plus);
      mu_p.push_back(tmp2);
    }else{
      Particle tmp1(*i, ptype_elec_minus);
      e_m.push_back(tmp1);
      Particle tmp2(*i, ptype_muon_minus);
      mu_m.push_back(tmp2);
    }
#if defined(BELLE_DEBUG)
    }
    catch(const std::string &e) {
      dout(Debugout::ERR,"utility") << "particle/functions/utility:makeLepton" << e << std::endl;
	Mdst_charged_Manager::get_manager().dump();
	Mdst_trk_Manager::get_manager().dump();
      continue;
    }
#endif    
  }
}

void 
withPCut(std::vector<Particle> &list, 
	 const double pL,
	 const double pR){
  for(int i=0;i<(int)list.size();++i){
    if(list[i].momentum().p().vect().mag() < pL ||
       list[i].momentum().p().vect().mag() > pR){
      list.erase(list.begin()+i);
      --i;
    }
  }
}

void 
withPCut(std::vector<Particle> &list, 
	 const double p){
  for(int i=0;i<(int)list.size();++i){
    if(list[i].momentum().p().vect().mag() < p){
      list.erase(list.begin()+i);
      --i;
    }
  }
}

HepLorentzVector
pStar(const Particle &p,
      const double elec, 
      const double posi,
      const double angle){  
  return pStar(p.momentum().p(),elec,posi,angle);
}

HepLorentzVector 
pStar(const Particle &p,
      const HepLorentzVector &el, 
      const HepLorentzVector &po){
  return pStar(p.momentum().p(),el,po);
}

HepLorentzVector
pStar(HepLorentzVector p,
      const double elec, 
      const double posi,
      const double angle){
  double pzP = sqrt(posi*posi-0.000511*0.000511);
  double pE  = sqrt(elec*elec-0.000511*0.000511);
  //HepLorentzVector boost_vector(0., 0., elec-posi, elec+posi);
  HepLorentzVector boost_vector(pE*sin(angle*0.001), 0., 
				pE*cos(angle*0.001)-pzP, elec+posi);
  p.boost(-boost_vector.boostVector());
  return p;
}

HepLorentzVector 
pStar(HepLorentzVector p,
      const HepLorentzVector &el, 
      const HepLorentzVector &po){
  HepLorentzVector boost_vector(el+po);
  p.boost(-boost_vector.boostVector());
  return p;
}

void 
withPSCut(std::vector<Particle> &list, 
	  const double pL,
	  const double pR){
  for(int i=0;i<(int)list.size();++i){
    HepLorentzVector P(pStar(list[i]));
    if(P.vect().mag() < pL ||
       P.vect().mag() > pR){
      list.erase(list.begin()+i);
      --i;
    }
  }
}

void 
withPSCut(std::vector<Particle> &list, 
	  const double p){
  for(int i=0;i<(int)list.size();++i){
    HepLorentzVector P(pStar(list[i]));
    if(P.vect().mag() < p){
      list.erase(list.begin()+i);
      --i;
    }
  }
}

void 
withSVD(std::vector<Particle> &list, 
	const unsigned n){ // z
  for(int i=0;i<(int)list.size();++i){
    if(list[i].mdstCharged()){
      int mhyp = 2;//pion
      if(abs(list[i].pType().lund()) == 321) mhyp = 3;
      if(abs(list[i].pType().lund()) == 11)  mhyp = 0;
      if(abs(list[i].pType().lund()) == 13)  mhyp = 1;
      if(abs(list[i].pType().lund()) == 2212)mhyp = 4;
      if(list[i].mdstCharged().trk().mhyp(mhyp).nhits(4) < (int)n){
	list.erase(list.begin()+i);
	--i;
      }
    }
  }
}

void 
withSVD2(std::vector<Particle> &list, 
	 const unsigned nR,
	 const unsigned nZ){
  for(int i=0;i<(int)list.size();++i){
    if(list[i].mdstCharged()){
      int mhyp = 2;//pion
      if(abs(list[i].pType().lund()) == 321) mhyp = 3;
      if(abs(list[i].pType().lund()) == 11)  mhyp = 0;
      if(abs(list[i].pType().lund()) == 13)  mhyp = 1;
      if(abs(list[i].pType().lund()) == 2212)mhyp = 4;
      if(list[i].mdstCharged().trk().mhyp(mhyp).nhits(3) < (int)nR ||
	 list[i].mdstCharged().trk().mhyp(mhyp).nhits(4) < (int)nZ){
	list.erase(list.begin()+i);
	--i;
      }
    }
  }
}

void 
withSVD(std::vector<Particle> &list, 
	const unsigned n,
	const unsigned childID){ // z
  for(int i=0;i<(int)list.size();++i){
    if(list[i].relation().child(childID).mdstCharged()){
      int mhyp = 2;//pion
      if(abs(list[i].relation().child(childID).pType().lund()) == 321) mhyp = 3;
      if(abs(list[i].relation().child(childID).pType().lund()) == 11)  mhyp = 0;
      if(abs(list[i].relation().child(childID).pType().lund()) == 13)  mhyp = 1;
      if(abs(list[i].relation().child(childID).pType().lund()) == 2212)mhyp = 4;
      if(list[i].relation().child(childID).mdstCharged().trk().mhyp(mhyp).nhits(4) < (int)n){
	list.erase(list.begin()+i);
	--i;
      }
    }
  }
}

void 
withSVD2(std::vector<Particle> &list, 
	 const unsigned nR,
	 const unsigned nZ,
	 const unsigned childID){
  for(int i=0;i<(int)list.size();++i){
    if(list[i].relation().child(childID).mdstCharged()){
      int mhyp = 2;//pion
      if(abs(list[i].relation().child(childID).pType().lund()) == 321) mhyp = 3;
      if(abs(list[i].relation().child(childID).pType().lund()) == 11)  mhyp = 0;
      if(abs(list[i].relation().child(childID).pType().lund()) == 13)  mhyp = 1;
      if(abs(list[i].relation().child(childID).pType().lund()) == 2212)mhyp = 4;
      if(list[i].relation().child(childID).mdstCharged().trk().mhyp(mhyp).nhits(3) < (int)nR ||
	 list[i].relation().child(childID).mdstCharged().trk().mhyp(mhyp).nhits(4) < (int)nZ){
	list.erase(list.begin()+i);
	--i;
      }
    }
  }
}

void 
withMassCut(std::vector<Particle> &list,
	    const double pL,
	    const double pR){
  for(int i=0;i<(int)list.size();++i){
    if(list[i].momentum().mass() < pL ||
       list[i].momentum().mass() > pR){
      list.erase(list.begin()+i);
      --i;
    }
  }
}

void 
withMassCut(std::vector<Particle> &list, 
	    std::vector<Particle> &nlist,
	    const double pL,
	    const double pR){
  for(int i=0;i<(int)list.size();++i){
    if(list[i].momentum().mass() < pL ||
       list[i].momentum().mass() > pR){
      nlist.push_back(list[i]);
      list.erase(list.begin()+i);
      --i;
    }
  }
}

void 
withMassDifCut(std::vector<Particle> &list, 
	       const double pL,
	       const double pR,
	       const unsigned child){
  for(int i=0;i<(int)list.size();++i){
    if(!withMassDifCut(list[i],pL,pR,child)){
      list.erase(list.begin()+i);
      --i;
    }
  }
}

unsigned
withMassDifCut(const Particle &p,
	       const double pL,
	       const double pR,
	       const unsigned child){
  double massdif(p.momentum().mass()-
		 p.relation().child(child).momentum().mass());
  if(massdif < pL || massdif > pR)return 0;
  return 1;
}

double
kaonId(const Particle &p,
       int accq0, int tofq0,int cdcq0,
       int ids0, int idb0){
  atc_pid kid(accq0,tofq0,cdcq0,ids0,idb0);
  return kid.prob(&(p.mdstCharged()));
}

void
withKaonId(std::vector<Particle> &list, const double prob,
	   int accq0, int tofq0, int cdcq0,
	   int ids0, int idb0){
  atc_pid kid(accq0,tofq0,cdcq0,ids0,idb0);
  for(int i=0;i<(int)list.size();++i){
    if(list[i].mdstCharged() && 
       kid.prob(&(list[i].mdstCharged())) >= prob){
      ;
    }else{
      list.erase(list.begin()+i);
      --i;
    }
  }
}

void
withPionId(std::vector<Particle> &list, const double prob,
	   int accq0, int tofq0, int cdcq0,
	   int ids0, int idb0){
  atc_pid kid(accq0,tofq0,cdcq0,ids0,idb0);
  for(int i=0;i<(int)list.size();++i){
    if(list[i].mdstCharged() && 
       kid.prob(&(list[i].mdstCharged())) < prob){
      ;
    }else{
      list.erase(list.begin()+i);
      --i;
    }
  }
}

void 
withMuId(std::vector<Particle> &list, 
	 const unsigned th,
	 const double prob,
	 const unsigned whichCriteria){
  for(int i=0;i<(int)list.size();++i){
    if(list[i].mdstCharged()){
      Muid_mdst muid(list[i].mdstCharged());
      if(whichCriteria == 2){
	if(muid.Muon_likelihood() < prob){
	  list.erase(list.begin()+i);
	  --i;
	}
      }else{
	if(muid.Level() < (int)th){
	  list.erase(list.begin()+i);
	  --i;
	}
      }
    }else{
      list.erase(list.begin()+i);
      --i;
    }
  }
}

void
withMuonId(std::vector<Particle> &list, const unsigned th){
  for(int i=0;i<(int)list.size();++i){
    if(list[i].mdstCharged() &&
       list[i].mdstCharged().muid() &&
       list[i].mdstCharged().muid().muon() >= (int)th){
      ;
    }else{
      list.erase(list.begin()+i);
      --i;
    }
  }
}
void
withEId(std::vector<Particle> &list, const double th){
  for(int i=0;i<(int)list.size();++i){
    if(list[i].mdstCharged()){
      eid post_elid(list[i].mdstCharged());
      if(post_elid.prob() < th){
	list.erase(list.begin()+i);
	--i;
      }
    }else{
      list.erase(list.begin()+i);
      --i;
    }
  }
}
void
makeGenHepKPi(std::vector<Particle> &k_p, 
	      std::vector<Particle> &k_m, 
	      std::vector<Particle> &pi_p, 
	      std::vector<Particle> &pi_m){
  Gen_hepevt_Manager &genMgr  = Gen_hepevt_Manager::get_manager();

  //...Particle Type
  Ptype ptype_pion_plus("PI+");
  Ptype ptype_pion_minus("PI-");
  Ptype ptype_kaon_plus("K+");
  Ptype ptype_kaon_minus("K-");

  for(std::vector<Gen_hepevt>::iterator i = genMgr.begin();
      i != genMgr.end(); ++i){
    if((*i).idhep() == ptype_kaon_plus.lund()){
      k_p.push_back(Particle(*i));
      continue;
    }
    if((*i).idhep() == ptype_kaon_minus.lund()){
      k_m.push_back(Particle(*i));
      continue;
    }
    if((*i).idhep() == ptype_pion_plus.lund()){
      pi_p.push_back(Particle(*i));
      continue;      
    }
    if((*i).idhep() == ptype_pion_minus.lund()){
      pi_m.push_back(Particle(*i));
      continue;
    }
  }
}

void
makeGenHepLepton(std::vector<Particle> &e_p,
		 std::vector<Particle> &e_m,
		 std::vector<Particle> &mu_p,
		 std::vector<Particle> &mu_m){
  Gen_hepevt_Manager &genMgr  = Gen_hepevt_Manager::get_manager();

  //...Particle Type
  Ptype ptype_elec_plus("E+");
  Ptype ptype_elec_minus("E-");
  Ptype ptype_muon_plus("MU+");
  Ptype ptype_muon_minus("MU-");

  for(std::vector<Gen_hepevt>::iterator i = genMgr.begin();
      i != genMgr.end(); ++i){
    if((*i).idhep() == ptype_elec_plus.lund()){
      e_p.push_back(Particle(*i));
      continue;
    }
    if((*i).idhep() == ptype_elec_minus.lund()){
      e_m.push_back(Particle(*i));
      continue;
    }
    if((*i).idhep() == ptype_muon_plus.lund()){
      mu_p.push_back(Particle(*i));
      continue;      
    }
    if((*i).idhep() == ptype_muon_minus.lund()){
      mu_m.push_back(Particle(*i));
      continue;
    }
  }
}

void 
setGenHepInfoF(std::vector<Particle> &list){
  //...final state particles
  if(list.size() == 0)return;

  for(std::vector<Particle>::iterator i = list.begin();
      i != list.end(); ++i){
    //...mdst_charged
    if(i->mdstCharged()){
      const Gen_hepevt & hep(get_hepevt(i->mdstCharged()));
      //      if (hep && i->pType().lund() == hep.idhep()){
      if (hep){
	i->relation().genHepevt(hep);
      }
    }
  }
}


VectorL 
child_p_sum(Particle part){
  if(part.nChildren() > 0){
    VectorL p = VectorL(0, 0, 0,0);
    for(int i = part.nChildren(); i >= 0; i--)
    {
      p += child_p_sum(part.child(i));
    }
    return p;
  }
  return part.p();
}

bool 
check_re(Particle part){
  float pm = part.pType().mass()*part.pType().mass();
  float m_child = child_p_sum(part).m2();
  return abs(pm - m_child) < 0.001;
}

bool
searchMother(const Gen_hepevt p, int mother_evt_id){
  if (!p.mother()) return false;
  if (p.mother().idhep() == mother_evt_id) return true;
  return searchMother(p.mother(), mother_evt_id);
}


bool 
checkMother(Particle &mother, std::vector<Particle> &childern){
  bool same_mother = true;
  VectorL b_mc =  VectorL(0, 0, 0, 0);
  for  (std::vector<Particle>::iterator chl = childern.begin(); chl!=childern.end(); ++chl){
    const Gen_hepevt chl_hep = chl->relation().genHepevt();
    for  (std::vector<Particle>::iterator chl__ = chl+1; chl__!=childern.end(); ++chl__){
      const Gen_hepevt chl_hep__ = chl__->relation().genHepevt();
      if (chl_hep == chl_hep__) return false;
    }
    same_mother = same_mother and searchMother(chl_hep, mother.pType().lund());
    if (not same_mother)
      return same_mother;
        b_mc = b_mc + VectorL(chl_hep.PX(), chl_hep.PY(), chl_hep.PZ(), chl_hep.E());
  }
  const Gen_hepevt m_hep = mother.relation().genHepevt();
  if (abs(b_mc.m2() - VectorL(m_hep.PX(), m_hep.PY(), m_hep.PZ(), m_hep.E()).m2()) < 0.0001 )
  return same_mother;
}

std::vector<Particle>
findchild (Particle &part){
  std::vector<Particle> childern;
  if (part.nChildren() == 0) { 
    std::vector<Particle> vect{part};
    return vect;
  }
  for (int i = 0; i < part.nChildren(); i++){
    std::vector<Particle> b = findchild(part.child(i));
    childern.insert(childern.end(), b.begin(), b.end());
  }
  return childern;
}

bool 
checkReal (Particle &part){
  std::vector<Particle> childern = findchild(part);
  return checkMother(part, childern);
}

void 
setUniqueGenHepInfoFBySvdAndDeltaP(std::vector<Particle> &list){
  int size = list.size();
  if(size == 0 || size == 1)return;
#if 0
  dout(Debugout::INFO,"utility") << "setUniqueGenHepInfoFBySvd--init--" << std::endl;
  for(int i=0;i<size;++i){
    if(list[i].genHepevt()){
      dout(Debugout::INFO,"utility") << i << ": hepID=" << list[i].genHepevt().get_ID() << std::endl;
    }
  }
#endif
  std::vector<int> list_checked;
  std::vector<int> list_best;
  for(int i=0;i<size;++i){
    if(!list[i].genHepevt())continue;    
    int best_index = i;
    Hep3Vector diff;
    diff.setX(list[i].p().x()-list[i].genHepevt().PX());
    diff.setY(list[i].p().y()-list[i].genHepevt().PY());
    diff.setZ(list[i].p().z()-list[i].genHepevt().PZ());
    double min_frac = diff.mag();
    int mhyp = 2;//pion
    if(abs(list[i].pType().lund()) == 321)mhyp = 3;
    if(abs(list[i].pType().lund()) == 11)mhyp = 0;
    if(abs(list[i].pType().lund()) == 13)mhyp = 1;
    if(abs(list[i].pType().lund()) == 2212)mhyp = 4;
    int svdHits = list[i].mdstCharged().trk().mhyp(mhyp).nhits(4);
    for(std::vector<int>::iterator it = list_checked.begin();
	it != list_checked.end();++it){
      if((int)(list[i].genHepevt().get_ID()) == *it)goto kokodayo;
    }
    for(int j=i+1;j<size;j++){
      if(!list[j].genHepevt())continue;
      if(list[i].genHepevt().get_ID() == list[j].genHepevt().get_ID()){
	Hep3Vector diff2;
	diff2.setX(list[j].p().x()-list[j].genHepevt().PX());
	diff2.setY(list[j].p().y()-list[j].genHepevt().PY());
	diff2.setZ(list[j].p().z()-list[j].genHepevt().PZ());
	int svdHits2 = list[j].mdstCharged().trk().mhyp(mhyp).nhits(4);
	if(svdHits2 > svdHits){
	  svdHits = svdHits2;
	  min_frac   = diff2.mag();
          best_index = j;
	}else if(svdHits2 == svdHits && diff2.mag() < min_frac){
	  min_frac   = diff2.mag();
	  best_index = j;
	}
      }
    }
    list_checked.push_back((int)(list[i].genHepevt().get_ID()));
    list_best.push_back(best_index);
    kokodayo:;
  }

  std::vector<int> remove_list;
  for(int i=0;i<size;i++){
    for(std::vector<int>::iterator it = list_best.begin();
	it != list_best.end();++it){
      if(i == *it)goto kokokana;
    }
    remove_list.push_back(i);
  kokokana:;
  }
  
  for(std::vector<int>::iterator it = remove_list.begin();
      it != remove_list.end();++it){
    list[*it].relation().resetGenHepevt();
  }
#if 0
  dout(Debugout::INFO,"utility") << "setUniqueGenHepInfoFBySvd--result--" << std::endl;
  for(int i=0;i<size;++i){
    if(list[i].genHepevt()){
      dout(Debugout::INFO,"utility") << i << ": hepID=" << list[i].genHepevt().get_ID() << std::endl;
    }
  }
#endif
}

void
addTrack2fit(kvertexfitter &kv, const std::vector<Particle> &list){
  for(unsigned i=0;i<list.size();++i){
    addTrack2fit(kv,list[i]);
  }
}



void
addTrack2fit(kvertexfitter &kv, const Particle &p){
  //No Correlation.
  kv.addTrack(p.momentum().p(),
	      p.momentum().x(),
	      p.momentum().dpx(),
	      p.pType().charge(),
	      p.momentum().mass());
  //p.pType().mass());
}

void
addBeam2fit(kvertexfitter &kv, 
	    const HepPoint3D &beam, 
	    const HepSymMatrix &errBeam){
  //No Correlation.
  kv.initialVertex(beam);
  kv.beamProfile(errBeam);
}

// added by T.H 2006/05/09
void
addBeam2fit(kvertexfitter &kv)
{
	const int flag_event_by_event_ip = 1;
	addBeam2fit(kv,
		IpProfile::position(flag_event_by_event_ip),
		IpProfile::position_err_b_life_smeared(flag_event_by_event_ip));
}

// added by T.H 2006/05/09, modified by T.H 2006/06/09
void
addTube2fit(kvertexfitter &kv)
{
	const int flag_event_by_event_ip = 1;
  kv.initialVertex(IpProfile::position(flag_event_by_event_ip));
	kv.tubeProfile(IpProfile::ip_tube_kfitterparticle(flag_event_by_event_ip));
}

void
addTrack2fit(kmassvertexfitter &kmv, const std::vector<Particle> &list){
  for(unsigned i=0;i<list.size();++i){
    addTrack2fit(kmv,list[i]);
  }
}

void
addTrack2fit(kmassvertexfitter &kmv, const Particle &p){
  //No Correlation.
  kmv.addTrack(p.momentum().p(),
	       p.momentum().x(),
	       p.momentum().dpx(),
	       p.pType().charge(),
	       p.momentum().mass());
	       //p.pType().mass());
}

void
addTrack2fit(kmassfitter &km, const std::vector<Particle> &list){
  for(unsigned i=0;i<list.size();++i){
    addTrack2fit(km,list[i]);
  }
}

void
addTrack2fit(kmassfitter &km, const Particle &p){
  //No Correlation.
  km.addTrack(p.momentum().p(),
	      p.momentum().x(),
	      p.momentum().dpx(),
	      p.pType().charge(),
	      p.momentum().mass());
	      //p.pType().mass());
}

unsigned
removeParticle(std::vector<Particle> &list, const Particle &p){
  unsigned count(0);
  for(unsigned i=0;i<list.size();++i){
#if 0
    if(list[i].relation().isIdenticalWith(p.relation())){
      list.erase(list.begin()+i);
      --i;
      ++count;
    }
#else
    if(list[i].relation().mdstCharged() && p.relation().mdstCharged()){
      if(list[i].relation().mdstCharged().get_ID() ==
	 p.relation().mdstCharged().get_ID()){
	list.erase(list.begin()+i);
	--i;
	++count;
      }
    }
#endif
  }
  return count;  
}

void 
calcuFoxWolfram(const std::vector<Particle> &list,
		double *r,
		const double e, 
		const double p,
		const double angle){
  std::vector<Hep3Vector> vec;
  //HepLorentzVector boost_vector(0., 0., e-p, e+p);
  double pzP = sqrt(p*p-0.000511*0.000511);
  double pE  = sqrt(e*e-0.000511*0.000511);
  HepLorentzVector boost_vector(pE*sin(angle*0.001), 0.,
				pE*cos(angle*0.001)-pzP, e+p);
  for(int i=0;i<(int)list.size();i++){
    HepLorentzVector b0(list[i].p());
    b0.boost(-(boost_vector.boostVector()));
    Hep3Vector tmp(b0.vect());
    vec.push_back(tmp);
  }
  
  FoxWolfram t;

  t = foxwolfram(vec.begin(), vec.end(), SelfFunc(Hep3Vector()));

  r[0] = t.R(0);
  r[1] = t.R(1);
  r[2] = t.R(2);
}

void 
calcuFoxWolfram(const std::vector<Particle> &list, 
		double *r,
		const HepLorentzVector &el,
		const HepLorentzVector &po){
  std::vector<Hep3Vector> vec;
  HepLorentzVector boost_vector(el+po);
  for(int i=0;i<(int)list.size();i++){
    HepLorentzVector b0(list[i].p());
    b0.boost(-(boost_vector.boostVector()));
    Hep3Vector tmp(b0.vect());
    vec.push_back(tmp);
  }
  
  FoxWolfram t;

  t = foxwolfram(vec.begin(), vec.end(), SelfFunc(Hep3Vector()));

  r[0] = t.R(0);
  r[1] = t.R(1);
  r[2] = t.R(2);
}

void 
massCut(const std::vector<Particle> &p1, 
	std::vector<Particle> &p2, 
	const double width){
  for(std::vector<Particle>::const_iterator i = p1.begin();
      i != p1.end(); ++i){
    double mass = (*i).momentum().mass();
    double nominal_mass = (*i).pType().mass();
    if(nominal_mass - width <= mass && mass <= nominal_mass + width){
      p2.push_back(*i);
    }
  }
}

void 
massDifCut(const std::vector<Particle> &p1, 
	   std::vector<Particle> &p2, 
	   const double width, 
	   const unsigned child){
  for(std::vector<Particle>::const_iterator i = p1.begin();
      i != p1.end(); ++i){
    double massDif = (*i).momentum().mass()-(*i).relation().child(child).momentum().mass();
    double nominal_massDif = (*i).pType().mass()-(*i).relation().child(child).pType().mass();
    if(nominal_massDif - width <= massDif && massDif <= nominal_massDif + width){
      p2.push_back(*i);
    }
  }
}

void 
deepCopy(std::vector<Particle> &p1, 
	 std::vector<Particle> &p2){
  for(std::vector<Particle>::iterator i = p1.begin();
      i != p1.end(); ++i){
    p2.push_back((*i).deepCopy());
  }
}

void 
deleteDeepCopiedObjects(std::vector<Particle> &p){
  for(std::vector<Particle>::iterator i = p.begin();
      i != p.end(); ++i){
    (*i).deepDelete();
  }
}

unsigned setGenHepInfoR(Particle &);
unsigned setGenHepInfoR_sub2(Particle &);
unsigned setGenHepInfoR_sub3(Particle &);
unsigned setGenHepInfoR_sub4(Particle &);
unsigned setGenHepInfoR_sub5(Particle &);

void 
setGenHepInfoR(std::vector<Particle> &list){
  //...reconstructed particles
  if(list.size() == 0)return;
  for(std::vector<Particle>::iterator i = list.begin();
      i != list.end(); ++i){
    setGenHepInfoR(*i);
  }
}

unsigned
setGenHepInfoR(Particle &p){
  //...supports 2, 3, 4, and 5
  
  //...returns = 0 : can NOT set.
  //...returns = 1 : can set.
  
  switch(p.relation().nChildren()){
  case 2:
    return setGenHepInfoR_sub2(p);
  case 3:
    return setGenHepInfoR_sub3(p);
  case 4:
    return setGenHepInfoR_sub4(p);
  case 5:
    return setGenHepInfoR_sub5(p);
  default:
    return 0;
  }
}

unsigned
setGenHepInfoR_sub2(Particle &p){
  if(!p.relation().child(0).genHepevt() || 
     !p.relation().child(1).genHepevt()){
    return 0;
  }

  if(p.relation().child(0).genHepevt().idhep() == 0 ||
     p.relation().child(1).genHepevt().idhep() == 0){
    return 0;
  }

  if(p.relation().child(0).genHepevt().get_ID() == 
     p.relation().child(1).genHepevt().get_ID()){
    return 0;
  }
  if(p.relation().child(0).genHepevt().mother() &&
     p.relation().child(1).genHepevt().mother()){
    if(p.relation().child(0).genHepevt().mother().get_ID() ==
       p.relation().child(1).genHepevt().mother().get_ID()){
      if(p.relation().child(0).genHepevt().mother().daLast() -
	 p.relation().child(0).genHepevt().mother().daFirst() + 1 == 2){
	if(p.pType().lund() == p.relation().child(0).genHepevt().mother().idhep()){
	  p.relation().genHepevt(p.relation().child(0).genHepevt().mother());
	  return 1;
	}
      }
    }
  }
  
  return 0;
}

unsigned
setGenHepInfoR_sub3(Particle &p){
  if(!p.relation().child(0).genHepevt()||
     !p.relation().child(1).genHepevt()||
     !p.relation().child(2).genHepevt()){
    return 0;
  }

  if(p.relation().child(0).genHepevt().idhep() == 0 ||
     p.relation().child(1).genHepevt().idhep() == 0 ||
     p.relation().child(2).genHepevt().idhep() == 0){
    return 0;
  }
  
  if(p.relation().child(0).genHepevt().get_ID() == 
     p.relation().child(1).genHepevt().get_ID() ||
     p.relation().child(0).genHepevt().get_ID() == 
     p.relation().child(2).genHepevt().get_ID() ||
     p.relation().child(1).genHepevt().get_ID() == 
     p.relation().child(2).genHepevt().get_ID()){
    return 0;
  }

  if(p.relation().child(0).genHepevt().mother() &&
     p.relation().child(1).genHepevt().mother() &&
     p.relation().child(2).genHepevt().mother()){
    if(p.relation().child(0).genHepevt().mother().get_ID() ==
       p.relation().child(1).genHepevt().mother().get_ID() &&
       p.relation().child(0).genHepevt().mother().get_ID() ==
       p.relation().child(2).genHepevt().mother().get_ID()){
      if(p.relation().child(0).genHepevt().mother().daLast() -
	 p.relation().child(0).genHepevt().mother().daFirst() + 1 == 3){
	if(p.pType().lund() == p.relation().child(0).genHepevt().mother().idhep()){
	  p.relation().genHepevt(p.relation().child(0).genHepevt().mother());
	  return 1;
	}
      }
    }
  }
  
  return 0;
}


unsigned
setGenHepInfoR_sub4(Particle &p){
  if(!p.relation().child(0).genHepevt() ||
     !p.relation().child(1).genHepevt() ||
     !p.relation().child(2).genHepevt() ||
     !p.relation().child(3).genHepevt()){
    return 0;
  }
  if(p.relation().child(0).genHepevt().idhep() == 0 ||
     p.relation().child(1).genHepevt().idhep() == 0 ||
     p.relation().child(2).genHepevt().idhep() == 0 ||
     p.relation().child(3).genHepevt().idhep() == 0){
    return 0;
  }
  if(p.relation().child(0).genHepevt().get_ID() == 
     p.relation().child(1).genHepevt().get_ID() ||
     p.relation().child(0).genHepevt().get_ID() == 
     p.relation().child(2).genHepevt().get_ID() ||
     p.relation().child(0).genHepevt().get_ID() == 
     p.relation().child(3).genHepevt().get_ID() ||
     p.relation().child(1).genHepevt().get_ID() == 
     p.relation().child(2).genHepevt().get_ID() ||
     p.relation().child(1).genHepevt().get_ID() == 
     p.relation().child(3).genHepevt().get_ID() ||
     p.relation().child(2).genHepevt().get_ID() == 
     p.relation().child(3).genHepevt().get_ID()){
    return 0;
  }
  if(p.relation().child(0).genHepevt().mother() &&
     p.relation().child(1).genHepevt().mother() &&
     p.relation().child(2).genHepevt().mother() &&
     p.relation().child(3).genHepevt().mother()){
    if(p.relation().child(0).genHepevt().mother().get_ID() ==
       p.relation().child(1).genHepevt().mother().get_ID() &&
       p.relation().child(0).genHepevt().mother().get_ID() ==
       p.relation().child(2).genHepevt().mother().get_ID() &&
       p.relation().child(2).genHepevt().mother().get_ID() ==
       p.relation().child(3).genHepevt().mother().get_ID()){
      if(p.relation().child(0).genHepevt().mother().daLast() -
	 p.relation().child(0).genHepevt().mother().daFirst() + 1 == 4){
	if(p.pType().lund() == p.relation().child(0).genHepevt().mother().idhep()){
	  p.relation().genHepevt(p.relation().child(0).genHepevt().mother());
	  return 1;
	}
      }
    }
  }
 
  return 0;
}

unsigned
setGenHepInfoR_sub5(Particle &p){
  if(!p.relation().child(0).genHepevt() ||
     !p.relation().child(1).genHepevt() ||
     !p.relation().child(2).genHepevt() ||
     !p.relation().child(3).genHepevt() ||
     !p.relation().child(4).genHepevt()){
    return 0;
  }

  if(p.relation().child(0).genHepevt().idhep() == 0 ||
     p.relation().child(1).genHepevt().idhep() == 0 ||
     p.relation().child(2).genHepevt().idhep() == 0 ||
     p.relation().child(3).genHepevt().idhep() == 0 ||
     p.relation().child(4).genHepevt().idhep() == 0){
    return 0;
  }

  if(p.relation().child(0).genHepevt().get_ID() == 
     p.relation().child(1).genHepevt().get_ID() ||
     p.relation().child(0).genHepevt().get_ID() == 
     p.relation().child(2).genHepevt().get_ID() ||
     p.relation().child(0).genHepevt().get_ID() == 
     p.relation().child(3).genHepevt().get_ID() ||
     p.relation().child(0).genHepevt().get_ID() == 
     p.relation().child(4).genHepevt().get_ID() ||
     p.relation().child(1).genHepevt().get_ID() == 
     p.relation().child(2).genHepevt().get_ID() ||
     p.relation().child(1).genHepevt().get_ID() == 
     p.relation().child(3).genHepevt().get_ID() ||
     p.relation().child(1).genHepevt().get_ID() == 
     p.relation().child(4).genHepevt().get_ID() ||
     p.relation().child(2).genHepevt().get_ID() == 
     p.relation().child(3).genHepevt().get_ID() ||
     p.relation().child(2).genHepevt().get_ID() == 
     p.relation().child(4).genHepevt().get_ID() ||
     p.relation().child(3).genHepevt().get_ID() == 
     p.relation().child(4).genHepevt().get_ID()){
    return 0;
  }

  if(p.relation().child(0).genHepevt().mother() &&
     p.relation().child(1).genHepevt().mother() &&
     p.relation().child(2).genHepevt().mother() &&
     p.relation().child(3).genHepevt().mother() &&
     p.relation().child(4).genHepevt().mother()){
    if(p.relation().child(0).genHepevt().mother().get_ID() ==
       p.relation().child(1).genHepevt().mother().get_ID() &&
       p.relation().child(2).genHepevt().mother().get_ID() ==
       p.relation().child(3).genHepevt().mother().get_ID() &&
       p.relation().child(0).genHepevt().mother().get_ID() ==
       p.relation().child(2).genHepevt().mother().get_ID() &&
       p.relation().child(0).genHepevt().mother().get_ID() ==
       p.relation().child(4).genHepevt().mother().get_ID()){
      if(p.relation().child(0).genHepevt().mother().daLast() -
	 p.relation().child(0).genHepevt().mother().daFirst() + 1 == 5){
	if(p.pType().lund() == p.relation().child(0).genHepevt().mother().idhep()){
	  p.relation().genHepevt(p.relation().child(0).genHepevt().mother());
	  return 1;
	}
      }
    }
  }
 
  return 0;
}

double
beamEnergyConstraint(const Particle &b, const double e, const double p,
		     const double angle){
  //b -- generally B0 or B0B
  //e -- electron beam energy
  //p -- positron beam enegry
  //angle --- crossing angle btw e & p (mrad)
  //HepLorentzVector boost_vector(0., 0., e-p, e+p);
  double pzP = sqrt(p*p-0.000511*0.000511);
  double pE  = sqrt(e*e-0.000511*0.000511);
  HepLorentzVector boost_vector(pE*sin(angle*0.001), 0.,
				pE*cos(angle*0.001)-pzP, e+p);
  HepLorentzVector b0(b.p());
  b0.boost(-(boost_vector.boostVector()));
  //double mass2 = e*p - b0.vect().mag2();
  double mass2 = boost_vector.mag2()*0.25 - b0.vect().mag2();
  double mass  = (mass2 > 0.) ? sqrt(mass2) :  -sqrt(-mass2);
  return mass;
}

double 
beamEnergyConstraint(const Particle &b,
		     const HepLorentzVector &el, 
		     const HepLorentzVector &po){
  //b -- generally B0 or B0B
  //e -- electron beam vector
  //p -- positron beam vector
  HepLorentzVector boost_vector(el+po);
  HepLorentzVector b0(b.p());
  b0.boost(-(boost_vector.boostVector()));
  double mass2 = boost_vector.mag2()*0.25 - b0.vect().mag2();
  double mass  = (mass2 > 0.) ? sqrt(mass2) :  -sqrt(-mass2);
  return mass;
}

Hep3Vector
calcuThrust(const std::vector<Particle> &list, const double e, const double p,
	    const double angle){
  std::vector<Hep3Vector> vec;
  //HepLorentzVector boost_vector(0., 0., e-p, e+p);
  double pzP = sqrt(p*p-0.000511*0.000511);
  double pE  = sqrt(e*e-0.000511*0.000511);
  HepLorentzVector boost_vector(pE*sin(angle*0.001), 0.,
				pE*cos(angle*0.001)-pzP, e+p);
  for(unsigned i=0;i<list.size();++i){
    HepLorentzVector b0(list[i].p());
    b0.boost(-(boost_vector.boostVector()));
    Hep3Vector tmp(b0.vect());
    vec.push_back(tmp);
  }
  
  return thrust(vec.begin(), vec.end(), SelfFunc(Hep3Vector()));
}

Hep3Vector 
calcuThrust(const std::vector<Particle> &list, 
	    const HepLorentzVector &el, 
	    const HepLorentzVector &po){
  std::vector<Hep3Vector> vec;
  HepLorentzVector boost_vector(el+po);
  for(unsigned i=0;i<list.size();++i){
    HepLorentzVector b0(list[i].p());
    b0.boost(-(boost_vector.boostVector()));
    Hep3Vector tmp(b0.vect());
    vec.push_back(tmp);
  }
  
  return thrust(vec.begin(), vec.end(), SelfFunc(Hep3Vector())); 
}

void 
makeGamma(std::vector<Particle> &gamma){
  Mdst_gamma_Manager &gamma_mag = Mdst_gamma_Manager::get_manager();
  for(std::vector<Mdst_gamma>::iterator i = gamma_mag.begin();
      i != gamma_mag.end(); ++i){
    gamma.push_back(Particle(*i));
  }
}

void 
makeEcl(std::vector<Particle> &ecl){
  Mdst_ecl_Manager &ecl_mag = Mdst_ecl_Manager::get_manager();
  for(std::vector<Mdst_ecl>::iterator i = ecl_mag.begin();
      i != ecl_mag.end(); ++i){
    ecl.push_back(Particle(*i));
  }
}

void 
makePi0(std::vector<Particle> &pi0){
  Mdst_pi0_Manager &pi0Mgr = Mdst_pi0_Manager::get_manager();  
  for(std::vector<Mdst_pi0>::iterator i = pi0Mgr.begin();
      i != pi0Mgr.end(); ++i){
    pi0.push_back(Particle(*i));
  }
}

unsigned
makeMother(kvertexfitter &kv,
	   Particle &mother){
  unsigned n = kv.tracks();
  kmakemother kmm;
  for(unsigned i=0;i<n;++i){
    kmm.addTrack(kv.momentum(i),
		 kv.position(i),
		 kv.error(i),
		 mother.relation().child(i).pType().charge());
    kmm.errVertexTrack(kv.errVertexTrack(i));
    for(unsigned j=i+1;j<n;++j){
      kmm.correlation(kv.correlation(i,j));
    }
  }
  kmm.vertex(kv.vertex());
  kmm.errVertex(kv.errVertex());
  unsigned err = kmm.make();
  if(err != 0)return 0;
  mother.momentum().momentumPosition(kmm.momentum(),
				     kmm.position(),
				     kmm.error());
  mother.momentum().decayVertex(kv.vertex(), kv.errVertex());
  return 1;
}

unsigned
makeMother(kmassvertexfitter &kmv,
	   Particle &mother){
  unsigned n = kmv.tracks();
  kmakemother kmm;
  for(unsigned i=0;i<n;++i){
    kmm.addTrack(kmv.momentum(i),
		 kmv.position(i),
		 kmv.error(i),
		 mother.relation().child(i).pType().charge());
    kmm.errVertexTrack(kmv.errVertexTrack(i));
    for(unsigned j=i+1;j<n;++j){
      kmm.correlation(kmv.correlation(i,j));
    }
  }
  kmm.vertex(kmv.vertex());
  kmm.errVertex(kmv.errVertex());
  unsigned err = kmm.make();
  if(err != 0)return 0;
  mother.momentum().momentumPosition(kmm.momentum(),
				     kmm.position(),
				     kmm.error());
  mother.momentum().decayVertex(kmv.vertex(), kmv.errVertex());
  return 1;
}

unsigned
makeMother(kmassfitter &km,
	   Particle &mother){
  unsigned n = km.tracks();
  kmakemother kmm;
  for(unsigned i=0;i<n;++i){
    kmm.addTrack(km.momentum(i),
		 km.position(i),
		 km.error(i),
		 mother.relation().child(i).pType().charge());
    if(km.fitWithVertex())kmm.errVertexTrack(km.errVertexTrack(i));
    for(unsigned j=i+1;j<n;++j){
      kmm.correlation(km.correlation(i,j));
    }
  }
  kmm.vertex(km.vertex());
  if(km.fitWithVertex()){
    kmm.errVertex(km.errVertex());
  }
  unsigned err = kmm.make();
  if(err != 0)return 0;
  mother.momentum().momentumPosition(kmm.momentum(),
				     kmm.position(),
				     kmm.error());
  if(km.fitWithVertex())
    mother.momentum().decayVertex(km.vertex(), km.errVertex());
  return 1;
}

unsigned
makeMother(kmakemother &kmm,
	   kvertexfitter &kv,
	   Particle &mother,
	   const unsigned notMake){
  unsigned n = kv.tracks();
  for(unsigned i=0;i<n;++i){
    kmm.addTrack(kv.momentum(i),
		 kv.position(i),
		 kv.error(i),
		 mother.relation().child(i).pType().charge());
    kmm.errVertexTrack(kv.errVertexTrack(i));
    for(unsigned j=i+1;j<n;++j){
      kmm.correlation(kv.correlation(i,j));
    }
  }
  kmm.vertex(kv.vertex());
  kmm.errVertex(kv.errVertex());
  if(notMake == 0){
    unsigned err = kmm.make();
    if(err != 0)return 0;
    mother.momentum().momentumPosition(kmm.momentum(),
				       kmm.position(),
				       kmm.error());
    mother.momentum().decayVertex(kv.vertex(), kv.errVertex());
  }
  return 1;
}

unsigned
makeMother(kmakemother &kmm,
	   kmassvertexfitter &kmv,
	   Particle &mother,
	   const unsigned notMake){
  unsigned n = kmv.tracks();
  for(unsigned i=0;i<n;++i){
    kmm.addTrack(kmv.momentum(i),
		 kmv.position(i),
		 kmv.error(i),
		 mother.relation().child(i).pType().charge());
    kmm.errVertexTrack(kmv.errVertexTrack(i));
    for(unsigned j=i+1;j<n;++j){
      kmm.correlation(kmv.correlation(i,j));
    }
  }
  kmm.vertex(kmv.vertex());
  kmm.errVertex(kmv.errVertex());
  if(notMake == 0){
    unsigned err = kmm.make();
    if(err != 0)return 0;
    mother.momentum().momentumPosition(kmm.momentum(),
				       kmm.position(),
				       kmm.error());
    mother.momentum().decayVertex(kmv.vertex(), kmv.errVertex());
  }
  return 1;
}

unsigned
makeMother(kmakemother &kmm,
	   kmassfitter &km,
	   Particle &mother,
	   const unsigned notMake){
  unsigned n = km.tracks();
  for(unsigned i=0;i<n;++i){
    kmm.addTrack(km.momentum(i),
		 km.position(i),
		 km.error(i),
		 mother.relation().child(i).pType().charge());
    if(km.fitWithVertex())kmm.errVertexTrack(km.errVertexTrack(i));
    for(unsigned j=i+1;j<n;++j){
      kmm.correlation(km.correlation(i,j));
    }
  }
  kmm.vertex(km.vertex());
  if(km.fitWithVertex()){
    kmm.errVertex(km.errVertex());
  }
  if(notMake == 0){
    unsigned err = kmm.make();
    if(err != 0)return 0;
    mother.momentum().momentumPosition(kmm.momentum(),
				       kmm.position(),
				       kmm.error());
    if(km.fitWithVertex())
      mother.momentum().decayVertex(km.vertex(), km.errVertex());
  }
  return 1;
}

void 
eraseVector(std::vector<Particle> &plist){
  plist.erase(plist.begin(),plist.end());
}

bool
findMotherParticle(int & motherIdhep,
                   int & now,
                   int & next,
                   int & motherID,
                   int &flag){
  next = 0;
  Gen_hepevt_Manager &genMgr = Gen_hepevt_Manager::get_manager();
  if(now <= 1)return false;
  if(!(genMgr[now-1]))return false;
  if(!(genMgr[now-1].mother()))return false;
  if(genMgr[now-1].mother().idhep() == motherIdhep){
    motherID = (int)(genMgr[now-1].mother().get_ID());
    flag = 1;
    return true;
  }else{
    next = (int)(genMgr[now-1].mother().get_ID());
    return false;
  }
}

void
fillChildId(int id, std::set<int,std::less<int> > &list){
  if(id < 1)return;
  Gen_hepevt_Manager &genMgr = Gen_hepevt_Manager::get_manager();
  if(genMgr[id-1]){
    if(genMgr[id-1].daFirst() == 0){
      list.insert(id);
    }else{
      int first = genMgr[id-1].daFirst();
      int last  = genMgr[id-1].daLast();
      for(int i=first;i<=last;++i){
        fillChildId(i,list);
      }
    }
  }
}

unsigned
setGenHepInfoR2(Particle &p){
  //...returns = 0 : can NOT set.
  //...returns = 1 : can set.
  unsigned nC = p.relation().nChildren();
  if(nC == 0)return 0;
  for(int i=0;i<(int)nC;++i)
    if(!(p.relation().child(i).genHepevt()))return 0;

  int motherIdhep = p.pType().lund();
  int now = (int)(p.relation().child(0).genHepevt().get_ID());
  int next, motherID;
  int flag = 0;
  while(!findMotherParticle(motherIdhep,now,next,motherID,flag)){
    if(next == 0)break;
    now = next;
  }
  if(flag == 0)return 0;
  
  // find a mother panther ID using child0!!!
  std::set<int,std::less<int> > childID;

  for(int i=0;i<(int)p.relation().nChildren();++i){
    fillChildId((int)(p.relation().child(i).genHepevt().get_ID()),childID);
  }
  std::set<int,std::less<int> > childIDofMother;
  fillChildId(motherID,childIDofMother);

  if(childID.size() > 0 &&
     childID.size() == childIDofMother.size()){

    std::set<int,std::less<int> >::iterator child1 = childID.begin();
    std::set<int,std::less<int> >::iterator child2 = childIDofMother.begin();

    int counter = 0;
    while(child1 != childID.end() &&
          *child1 == *child2){
      ++child1;
      ++child2;
      ++counter;
    }
    if((int)childID.size() != counter)return 0;
#if 0
    dout(Debugout::INFO,"utility") << "CHILD1 = " << childID.size() << std::endl;
    int index = 0;
    for(set<int,less<int> >::iterator i=childID.begin();
        i!=childID.end();++i){
      dout(Debugout::INFO,"utility") << index << ": " << *i << std::endl;
      ++index;
    }
    dout(Debugout::INFO,"utility") << "CHILD2 = " << childIDofMother.size() << std::endl;
    index = 0;
    for(set<int,less<int> >::iterator i=childIDofMother.begin();
        i!=childIDofMother.end();++i){
      dout(Debugout::INFO,"utility") << index << ": " << *i << std::endl;
      ++index;
    }
    Gen_hepevt_Manager &genMgr = Gen_hepevt_Manager::get_manager();
    for(int i=0;i<genMgr.count();++i){
      dout(Debugout::INFO,"utility") << genMgr[i].get_ID() << ": " 
           << genMgr[i].idhep() << "  "
           << genMgr[i].mo(0) << "  "
           << genMgr[i].mo(1) << "  "
           << genMgr[i].da(0) << "  "
           << genMgr[i].da(1) << std::endl;
    }
#endif
    p.relation().genHepevt(p.relation().child(0).genHepevt().mother());
    return 1;
  }
  return 0;
}

void
setGenHepInfoR2(std::vector<Particle> &list){
  if(list.size() == 0)return;
  for(std::vector<Particle>::iterator i = list.begin();
      i != list.end(); ++i){
    setGenHepInfoR2(*i);
  }
}

void 
setGenHepInfoG(std::vector<Particle> &list){
  //...gamma
  if(list.size() == 0)return;

  for(std::vector<Particle>::iterator i = list.begin();
      i != list.end(); ++i){
    //...mdst_gamma
    if(i->mdstGamma()){
      const Gen_hepevt & hep(gen_level(get_hepevt(i->mdstGamma())));
      if (hep && i->pType().lund() == hep.idhep()){
	i->relation().genHepevt(hep);
      }
    }
  }
}

void 
setUniqueGenHepInfoByDeltaP(std::vector<Particle> &list){
  //...this function can be used for MDST_CHARGED, MDST_GAMMA at least.
  int size = list.size();
  if(size == 0 || size == 1)return;
  std::vector<int> list_checked;
  std::vector<int> list_best;
  for(int i=0;i<size;++i){
    if(!list[i].genHepevt())continue;    

    for(std::vector<int>::iterator it = list_checked.begin();
        it != list_checked.end();++it){
      if((int)(list[i].genHepevt().get_ID()) == *it){
        //...already associated particle
        goto kokodayo2;
      }
    }
    {
      int best_index = i;
      Hep3Vector diff;
      diff.setX(list[i].p().x()-list[i].genHepevt().PX());
      diff.setY(list[i].p().y()-list[i].genHepevt().PY());
      diff.setZ(list[i].p().z()-list[i].genHepevt().PZ());
      double min_frac = diff.mag();
      //...search best particles with the same GehHepID.
      for(int j=i+1;j<size;j++){
	if(!list[j].genHepevt())continue;
	if(list[i].genHepevt().get_ID() == list[j].genHepevt().get_ID()){
	  Hep3Vector diff2;
	  diff2.setX(list[j].p().x()-list[j].genHepevt().PX());
	  diff2.setY(list[j].p().y()-list[j].genHepevt().PY());
	  diff2.setZ(list[j].p().z()-list[j].genHepevt().PZ());
	  if(diff2.mag() < min_frac){
	    min_frac   = diff2.mag();
          best_index = j;
	  }
	}
      }
      //...checked GenHepID
      list_checked.push_back((int)(list[i].genHepevt().get_ID()));
      //...index of the best particle
      list_best.push_back(best_index);
    }
  kokodayo2:;
  }

  //...search not-best particle.
  std::vector<int> remove_list;
  for(int i=0;i<size;i++){
    for(std::vector<int>::iterator it = list_best.begin();
        it != list_best.end();++it){
      if(i == *it)goto kokokana2;
    }
    //...index of not-best particle
    remove_list.push_back(i);
  kokokana2:;
  }
  
  //...reset GenHepInformation from not-best particle.
  for(std::vector<int>::iterator it = remove_list.begin();
      it != remove_list.end();++it){
    list[*it].relation().resetGenHepevt();
  }
  //dout(Debugout::INFO,"utility") << "reset Gamma = " << remove_list.size() << std::endl;
}

double
rc2wc(const double &xw, const double &yw,
      const double &xc, const double &yc,
      const double &xr, const double &yr)
{
  // c : xc = center
  // w : xw = (wire)position
  // r : xr = reference position
  const double crs = (xr-xc)*(yw-yc)-(yr-yc)*(xw-xc);
  const double dot = (xr-xc)*(xw-xc)+(yr-yc)*(yw-yc);
  return atan2(crs,dot);
}

bool
px2helix(const Hep3Vector &p,
	 const HepPoint3D &x, 	 
	 const double charge,
	 HepVector &a,
	 const double alpha)
{
  // pivot = (0,0,0)
  double pt = sqrt(p.x()*p.x()+p.y()*p.y());
  if(pt == 0.)return false;

  //HepVector a(5);
  a[2] = charge > 0. ? 1./pt : -1./pt;
  a[4] = p.z()/pt;

  double ux =  p.y()/pt;
  double uy = -p.x()/pt;
  double r  = alpha/a[2];

  double xc = x.x()+r*ux;
  double yc = x.y()+r*uy;

  double R = sqrt(xc*xc+yc*yc);

  a[0] = fabs(-r+R) < fabs(-r-R) ? -r+R : -r-R;

  a[1] = atan2((a[0]+r)*yc,(a[0]+r)*xc);
  if(a[1] < 0.)a[1] += 2.*M_PI;
  a[3] = x.z()+r*a[4]*rc2wc(x.x(),x.y(),xc,yc,0.,0.);

  return true;
}

bool
px2helix(const Hep3Vector &p,
	 const HepPoint3D &x,          
         const HepSymMatrix &dpx,
         const double charge,
         HepVector &a,
         HepSymMatrix &da,
         const double alpha)
{
  // charge = 1 or -1
  if(px2helix(p,x,charge,a,alpha)){
    HepMatrix dHdPX(5,6,0);
    
    // kappa
    dHdPX[2][0] = -p.x()*a[2]*a[2]*a[2]; // px
    dHdPX[2][1] = -p.y()*a[2]*a[2]*a[2]; // py

    // tanL
    dHdPX[4][0] = -a[2]*a[2]*a[2]*charge*p.z()*p.x(); // px
    dHdPX[4][1] = -a[2]*a[2]*a[2]*charge*p.z()*p.y(); // py
    dHdPX[4][2] =  a[2]*charge;  //pz

    double r = alpha/a[2];
    double X = x.x()+alpha*p.y()*charge;
    double Y = x.y()-alpha*p.x()*charge;
    double R = sqrt(X*X+Y*Y);
    if(R == 0. || X == 0.)return false;
    double s = fabs(-r+R) < fabs(-r-R) ? 1. : -1.; 
    double inK2 = p.x()*p.x()+p.y()*p.y(); // 1/kappa/kappa

    // drho
    dHdPX[0][0] = alpha*inK2*dHdPX[2][0]-s*Y*alpha*charge/R; // px
    dHdPX[0][1] = alpha*inK2*dHdPX[2][1]+s*X*alpha*charge/R; // py
    dHdPX[0][3] = s*X/R; // x
    dHdPX[0][4] = s*Y/R; // y

    // phi0
    double A = X*X/(X*X+Y*Y);
    dHdPX[1][0] = -A*alpha*charge/X; // px
    dHdPX[1][1] = -A*alpha*charge*Y/X/X; // py
    dHdPX[1][3] = -A*Y/X/X; // x
    dHdPX[1][4] =  A/X; // y

    double U = x.x()*p.x()+x.y()*p.y();
    double L = x.y()*p.x()-x.x()*p.y()-alpha*charge*inK2;
    if(L == 0. || U*U+L*L == 0.)return false;
    double B = L*L/(U*U+L*L);
    double phi = rc2wc(x.x(),x.y(),X,Y,0.,0.);

    // dz
    dHdPX[3][0] = -alpha*inK2*dHdPX[2][0]*a[4]*phi+
                   alpha/a[2]*dHdPX[4][0]*phi+
                   alpha/a[2]*a[4]*B*(x.x()/L-U/L/L*(x.y()-2.*alpha*charge*p.x())); // px
    dHdPX[3][1] = -alpha*inK2*dHdPX[2][1]*a[4]*phi+
                   alpha/a[2]*dHdPX[4][1]*phi+
                   alpha/a[2]*a[4]*B*(x.y()/L+U/L/L*(x.x()+2.*alpha*charge*p.y())); // py
    dHdPX[3][3] =  alpha/a[2]*a[4]*B*(p.x()/L+U/L/L*p.y()); // x
    dHdPX[3][4] =  alpha/a[2]*a[4]*B*(p.y()/L-U/L/L*p.x()); // y
    dHdPX[3][5] =  1.0;
    
    da = dpx.similarity(dHdPX);
    return true;
  }else{
    return false;
  }
}

double
calEcm(const double elec, 
       const double posi,
       const double angle)
{
  double pzP = sqrt(posi*posi-0.000511*0.000511);
  double pE  = sqrt(elec*elec-0.000511*0.000511);
  HepLorentzVector boost_vector(pE*sin(angle*0.001), 0., 
                                pE*cos(angle*0.001)-pzP, elec+posi);
  return boost_vector.mag();
}

bool
isHadronA(void)
{
  bool goodFlag = true;
  const double Ecm = calEcm();

  Evtcls_hadron_info_Manager &hadMgr = Evtcls_hadron_info_Manager::get_manager();
  if(hadMgr.count() != 0){
    if(hadMgr[0].Evis() < Ecm * 0.2)goodFlag = false;
    if(fabs(hadMgr[0].Pz()) > Ecm * 0.5)goodFlag = false;
    if(hadMgr[0].Esum() < Ecm * 0.025)goodFlag = false;
    if(hadMgr[0].Esum() > Ecm * 0.9)goodFlag = false;
    if(hadMgr[0].Ntrk() < 3)goodFlag = false;
  }else goodFlag = false;

  Evtvtx_primary_vertex_Manager &vtxMgr = Evtvtx_primary_vertex_Manager::get_manager();
  if(vtxMgr.count() != 0 && vtxMgr[0].quality() >= 2){
    if(fabs(vtxMgr[0].PV_z()) > 3.5)goodFlag = false;
    if(sqrt(vtxMgr[0].PV_x()*vtxMgr[0].PV_x()+
	    vtxMgr[0].PV_y()*vtxMgr[0].PV_y()) > 1.5)goodFlag = false;
  }

  return goodFlag;
}


void 
setGenHepInfoT(Particle &p){
  int nchildren = p.nChildren();
  if(!nchildren) return;

  // Check that genHep references exist;
  for(int i=0; i<nchildren; ++i)
    if(!p.relation().child(i).genHepevt()) return;

  // Check that child particles haven't same genHep reference;
  for(int i=0; i<nchildren-1; ++i)
    for(int j=i+1; j<nchildren; ++j)
      if(p.relation().child(i).genHepevt().get_ID() == 
         p.relation().child(j).genHepevt().get_ID() ) return;

  // Seek mother by the first daughter;
  const Gen_hepevt *mother(&(p.child(0).genHepevt()));
  while(mother->mother()){
    mother = &(mother->mother());
    if(mother->idhep() == p.pType().lund()) break;
  }
  if(mother->idhep() != p.pType().lund()) return;
  
  // Check for other children that have the same mother;
  for(int i=1; i<nchildren; ++i){
    const Gen_hepevt *tmp(&(p.child(i).genHepevt()));
    while(tmp->mother()){
      tmp = &tmp->mother();
      if(tmp == mother) break;
    }
    if(tmp != mother) return;
  }

  // Check that there are no another children from this mother;
  double e=0;
  for(int i=0; i<nchildren; ++i)
    e+=p.child(i).genHepevt().E();
  if(abs(e-mother->E())<.0001)
    p.relation().genHepevt(*mother);
}

void 
setGenHepInfoT(Particle &p, bool lost){
  int nchildren = p.nChildren();
  if(!nchildren) return;

  // Check that genHep references exist;
  for(int i=0; i<nchildren; ++i)
    if(!p.relation().child(i).genHepevt()) return;

  // Check that child particles haven't same genHep reference;
  for(int i=0; i<nchildren-1; ++i)
    for(int j=i+1; j<nchildren; ++j)
      if(p.relation().child(i).genHepevt().get_ID() == 
         p.relation().child(j).genHepevt().get_ID() ) return;

  // Seek mother by the first daughter;
  const Gen_hepevt *mother(&(p.child(0).genHepevt()));
  while(mother->mother()){
    mother = &(mother->mother());
    if(mother->idhep() == p.pType().lund()) break;
  }
  if(mother->idhep() != p.pType().lund()) return;
  
  // Check for other children that have the same mother;
  for(int i=1; i<nchildren; ++i){
    const Gen_hepevt *tmp(&(p.child(i).genHepevt()));
    while(tmp->mother()){
      tmp = &tmp->mother();
      if(tmp == mother) break;
    }
    if(tmp != mother) return;
  }

  // Check that there are no another children from this mother;
  double e=0;
  for(int i=0; i<nchildren; ++i)
    e+=p.child(i).genHepevt().E();
  if(abs(e-mother->E())<.0001 or lost)
    p.relation().genHepevt(*mother);
}

void 
setGenHepInfoT(std::vector<Particle> &p_list){
  for(unsigned int i=0; i<p_list.size(); ++i)
    setGenHepInfoT(p_list[i]);
}

void 
setGenHepInfoT(std::vector<Particle> &p_list, bool lost){
  for(unsigned int i=0; i<p_list.size(); ++i)
    setGenHepInfoT(p_list[i], lost);
}

/************ SetGenHepInfoP for PI0 ***********************************/
void 
setGenHepInfoP(Particle &p){
  if( !p.child(0) || !p.child(1) ) return;
  if( !p.child(0).mdstGamma() || !p.child(1).mdstGamma() ) return;
  if(  p.child(0).mdstGamma() ==  p.child(1).mdstGamma() ) return;

  const Gen_hepevt& hep0 = gen_level(get_hepevt(p.child(0).mdstGamma()));
  const Gen_hepevt& hep1 = gen_level(get_hepevt(p.child(1).mdstGamma()));
  if (hep0 && hep0.idhep() == 22) p.child(0).relation().genHepevt(hep0);
  if (hep1 && hep1.idhep() == 22) p.child(1).relation().genHepevt(hep1);
  if( p.child(0).genHepevt() && p.child(1).genHepevt() )
    if( p.child(0).genHepevt() !=  p.child(1).genHepevt() &&
	p.child(0).genHepevt().mother() == 
	p.child(1).genHepevt().mother() &&
	p.child(0).genHepevt().mother().idhep() == 111 )
      p.relation().genHepevt(p.child(0).genHepevt().mother());
}

void 
setGenHepInfoP(std::vector<Particle> &p_list){
  for(unsigned int i=0; i<p_list.size(); ++i)
    setGenHepInfoP(p_list[i]);
}

/*************** Make Mass Fit (used for pi0 fit) **********************/
void 
doMassFit(Particle &p){
  kmassfitter km;
  km.invariantMass(p.pType().mass());
  for(unsigned i=0; i<p.nChildren(); ++i)
    addTrack2fit(km,p.child(i));
  if(!km.fit())
    makeMother(km,p);
}

void 
doMassFit(std::vector<Particle> &p_list){
  for(unsigned int i=0; i<p_list.size(); ++i)
    doMassFit(p_list[i]);
}

/************ Set error matrix for gamma ***************************/
void 
setGammaError(Particle &p,
	      const HepPoint3D &gVertex,
	      const HepSymMatrix &errGVertex){
  if(!p.mdstGamma())return;

  HepSymMatrix errG(3,0);
  errG[0][0] = p.mdstGamma().ecl().error(0);
  errG[1][0] = p.mdstGamma().ecl().error(1);
  errG[1][1] = p.mdstGamma().ecl().error(2);
  errG[2][0] = p.mdstGamma().ecl().error(3);
  errG[2][1] = p.mdstGamma().ecl().error(4);
  errG[2][2] = p.mdstGamma().ecl().error(5);
  GammaParticle g(p.mdstGamma().ecl().energy(),
		  p.mdstGamma().ecl().phi(),
		  p.mdstGamma().ecl().theta(),
		  p.mdstGamma().ecl().r(),
		  errG);
  g.vertex(gVertex, errGVertex);
  
  HepSymMatrix errGnew(7,0);
  errGnew.sub(1,g.errorMomentumEnergy());
  errGnew.sub(5,errGVertex);
  
  p.momentum().momentumPosition(g.momentumEnergy(),
				gVertex,errGnew);
}

/************ Set error matrix for two gammas of pi0 ***************************/
void 
setGammasError(Particle &p,
	       const HepPoint3D &gVertex,
	       const HepSymMatrix &errGVertex){
  if( !p.child(0) || !p.child(1) ) return;
  if( !p.child(0).mdstGamma() || 
      !p.child(1).mdstGamma() ) return;
  for(int i=0; i<2; ++i){
    setGammaError(p.child(i),
		  gVertex,
		  errGVertex);
  }
}

void 
getEventInfo(int &expNo,
	     int &runNo,
	     int &evtNo,
	     bool &McFlag)
{
  expNo = runNo = evtNo = 0;
  McFlag = false;
  Belle_event_Manager &evtMgr = Belle_event_Manager::get_manager();
  if(evtMgr.count() != 0){
    if(evtMgr[0].ExpMC() == 2)McFlag = true;
    const int MASK28BIT = 0x0FFFFFFF;
    expNo = evtMgr[0].ExpNo();
    runNo = evtMgr[0].RunNo();
    evtNo = (int)(evtMgr[0].EvtNo() & MASK28BIT);
  }
}

unsigned 
getNRSvd(const Particle &p, const unsigned id)
{
  if(p.mdstCharged() && p.mdstCharged().trk()){
    unsigned tmpid = id;
    if(tmpid == TYPE_AUTO){
      switch(abs(p.lund())){
      case 11:
	tmpid = TYPE_E;
	break;
      case 13:
	tmpid = TYPE_MU;
	break;
      case 321:
	tmpid = TYPE_K;
	break;
      case 2212:
	tmpid = TYPE_P;
	break;
      default:
	tmpid = TYPE_PI;
      }
    }
    return p.mdstCharged().trk().mhyp(tmpid).nhits(3);
  }else if(p.mdstVee2()){
    int nsvdrp = getNRSvdVee(p, 0);
    int nsvdrm = getNRSvdVee(p, 1);
    return (nsvdrp < nsvdrm ? nsvdrp : nsvdrm);
  }else{
    return 0;
  }     
}

unsigned 
getNZSvd(const Particle &p, const unsigned id)
{
  if(p.mdstCharged() && p.mdstCharged().trk()){
    unsigned tmpid = id;
    if(tmpid == TYPE_AUTO){
      switch(abs(p.lund())){
      case 11:
	tmpid = TYPE_E;
	break;
      case 13:
	tmpid = TYPE_MU;
	break;
      case 321:
	tmpid = TYPE_K;
	break;
      case 2212:
	tmpid = TYPE_P;
	break;
      default:
	tmpid = TYPE_PI;
      }
    }
    return p.mdstCharged().trk().mhyp(tmpid).nhits(4);
  }else if(p.mdstVee2()){
    int nsvdzp = getNZSvdVee(p, 0);
    int nsvdzm = getNZSvdVee(p, 1);
    return (nsvdzp < nsvdzm ? nsvdzp : nsvdzm);
  }else{
    return 0;
  }
}

unsigned 
getNRSvdVee(const Particle &p, const unsigned id)
{
  if(p.mdstVee2()){
    Mdst_vee_daughters_add_Manager& vda_mgr = Mdst_vee_daughters_add_Manager::get_manager();
    if(vda_mgr.count() == 0 || !(p.mdstVee2().daut())){
      return getNRSvd(p.child(id));
    }else{
      Mdst_vee_daughters_add &vda = vda_mgr[p.mdstVee2().daut_ID()-1];
      return (id == 0? vda.nhits_p(0) : id == 1 ? vda.nhits_m(0) : 0);
    }
  }else{
    return 0;
  }
}

unsigned 
getNZSvdVee(const Particle &p, const unsigned id)
{
  if(p.mdstVee2()){
    Mdst_vee_daughters_add_Manager& vda_mgr = Mdst_vee_daughters_add_Manager::get_manager();
    if(vda_mgr.count() == 0 || !(p.mdstVee2().daut())){
      return getNZSvd(p.child(id));
    }else{
      Mdst_vee_daughters_add &vda = vda_mgr[p.mdstVee2().daut_ID()-1];
      return (id == 0? vda.nhits_p(1) : id == 1 ? vda.nhits_m(1) : 0);
    }
  }else{
    return 0;
  }
}

unsigned getNRSvdVeeP(const Particle &p){return getNRSvdVee(p, 0);}
unsigned getNRSvdVeeM(const Particle &p){return getNRSvdVee(p, 1);}
unsigned getNZSvdVeeP(const Particle &p){return getNZSvdVee(p, 0);}
unsigned getNZSvdVeeM(const Particle &p){return getNZSvdVee(p, 1);}

unsigned 
getHitSvd(const Particle &p, const unsigned id)
{
  if(p.mdstCharged() && p.mdstCharged().trk()){
    unsigned tmpid = id;
    if(tmpid == TYPE_AUTO){
      switch(abs(p.lund())){
      case 11:
	tmpid = TYPE_E;
	break;
      case 13:
	tmpid = TYPE_MU;
	break;
      case 321:
	tmpid = TYPE_K;
	break;
      case 2212:
	tmpid = TYPE_P;
	break;
      default:
	tmpid = TYPE_PI;
      }
    }
    return p.mdstCharged().trk().mhyp(tmpid).hit_svd();
  }else{
    return 0;
  }
}

unsigned 
getHitSvdVee(const Particle &p, const unsigned id)
{
  if(p.mdstVee2()){
    Mdst_vee_daughters_add_Manager& vda_mgr = Mdst_vee_daughters_add_Manager::get_manager();
    if(vda_mgr.count() == 0 || !(p.mdstVee2().daut())){
      return getHitSvd(p.child(id));
    }else{
      Mdst_vee_daughters_add &vda = vda_mgr[p.mdstVee2().daut_ID()-1];
      return (id == 0? vda.hit_p() : id == 1 ? vda.hit_m() : 0);
    }
  }else{
    return 0;
  }
}

unsigned getHitSvdVeeP(const Particle &p){return getHitSvdVee(p, 0);}
unsigned getHitSvdVeeM(const Particle &p){return getHitSvdVee(p, 1);}

Helix 
calMdstChargedHelix(const Particle &p, const unsigned id)
{
  if(p.mdstCharged() && p.mdstCharged().trk()){
    unsigned tmpid = id;
    if(tmpid == TYPE_AUTO){
      switch(abs(p.lund())){
      case 11:
	tmpid = TYPE_E;
	break;
      case 13:
	tmpid = TYPE_MU;
	break;
      case 321:
	tmpid = TYPE_K;
	break;
      case 2212:
	tmpid = TYPE_P;
	break;
      default:
	tmpid = TYPE_PI;
      }
    }
    HepPoint3D pivot(p.mdstCharged().trk().mhyp(tmpid).pivot(0),
		     p.mdstCharged().trk().mhyp(tmpid).pivot(1),
		     p.mdstCharged().trk().mhyp(tmpid).pivot(2));
    HepVector a(5);
    a[0] = p.mdstCharged().trk().mhyp(tmpid).helix(0);
    a[1] = p.mdstCharged().trk().mhyp(tmpid).helix(1);
    a[2] = p.mdstCharged().trk().mhyp(tmpid).helix(2);
    a[3] = p.mdstCharged().trk().mhyp(tmpid).helix(3);
    a[4] = p.mdstCharged().trk().mhyp(tmpid).helix(4);
    HepSymMatrix ea(5,0);
    ea[0][0] = p.mdstCharged().trk().mhyp(tmpid).error(0);
    ea[1][0] = p.mdstCharged().trk().mhyp(tmpid).error(1);
    ea[1][1] = p.mdstCharged().trk().mhyp(tmpid).error(2);
    ea[2][0] = p.mdstCharged().trk().mhyp(tmpid).error(3);
    ea[2][1] = p.mdstCharged().trk().mhyp(tmpid).error(4);
    ea[2][2] = p.mdstCharged().trk().mhyp(tmpid).error(5);
    ea[3][0] = p.mdstCharged().trk().mhyp(tmpid).error(6);
    ea[3][1] = p.mdstCharged().trk().mhyp(tmpid).error(7);
    ea[3][2] = p.mdstCharged().trk().mhyp(tmpid).error(8);
    ea[3][3] = p.mdstCharged().trk().mhyp(tmpid).error(9);
    ea[4][0] = p.mdstCharged().trk().mhyp(tmpid).error(10);
    ea[4][1] = p.mdstCharged().trk().mhyp(tmpid).error(11);
    ea[4][2] = p.mdstCharged().trk().mhyp(tmpid).error(12);
    ea[4][3] = p.mdstCharged().trk().mhyp(tmpid).error(13);
    ea[4][4] = p.mdstCharged().trk().mhyp(tmpid).error(14);
    return Helix(pivot,a,ea);
  }else{    
    HepPoint3D pivot;
    HepVector a(5);
    HepSymMatrix ea(5,0);
    return Helix(pivot, a, ea);
  }
}

void 
setGenHepInfoKs(Particle &p){
  if( !p.child(0) || !p.child(1) ) return;
  if( !p.child(0).mdstCharged() || !p.child(1).mdstCharged() ) return;
  if(  p.child(0).mdstCharged() ==  p.child(1).mdstCharged() ) return;

  const Gen_hepevt& hep0 = get_hepevt(p.child(0).mdstCharged());
  const Gen_hepevt& hep1 = get_hepevt(p.child(1).mdstCharged());
  if (hep0 && hep1 && (hep0.idhep()*hep1.idhep() > 0)) return;
  if (hep0 && abs(hep0.idhep()) == 211) p.child(0).relation().genHepevt(hep0);
  if (hep1 && abs(hep1.idhep()) == 211) p.child(1).relation().genHepevt(hep1);

  if( p.child(0).genHepevt() && p.child(1).genHepevt() )
    if( p.child(0).genHepevt() !=  p.child(1).genHepevt() &&
	p.child(0).genHepevt().mother() && p.child(1).genHepevt().mother() &&
	p.child(0).genHepevt().mother() == 
	p.child(1).genHepevt().mother() &&
	p.child(0).genHepevt().mother().idhep() == 310 &&
	p.child(0).genHepevt().mother().daLast() - 
	p.child(0).genHepevt().mother().daFirst() == 1){
      p.relation().genHepevt(p.child(0).genHepevt().mother());
    }else if( p.child(0).genHepevt() !=  p.child(1).genHepevt() &&
	      p.child(0).genHepevt().mother() && p.child(1).genHepevt().mother() &&
	      p.child(0).genHepevt().mother() == 
	      p.child(1).genHepevt().mother() &&
	      abs(p.child(0).genHepevt().mother().idhep()) == 311 &&
	      p.child(0).genHepevt().mother().daLast() - 
	      p.child(0).genHepevt().mother().daFirst() == 1 &&
	      p.child(0).genHepevt().mother().mother() &&
	      p.child(0).genHepevt().mother().mother().idhep() == 310 &&
	      p.child(0).genHepevt().mother().mother().daLast() - 
	      p.child(0).genHepevt().mother().mother().daFirst() == 0){
      p.relation().genHepevt(p.child(0).genHepevt().mother().mother());
    }
}

void 
setGenHepInfoKs(std::vector<Particle> &p_list){
  for(unsigned int i=0; i<p_list.size(); ++i)
    setGenHepInfoKs(p_list[i]);
}

void 
makeKs(std::vector<Particle> &ks, unsigned type){
  if(type == 2){
    Mdst_vee2_Manager &veeMgr = Mdst_vee2_Manager::get_manager();  
    for(std::vector<Mdst_vee2>::iterator i = veeMgr.begin();
        i != veeMgr.end(); ++i){
      Particle tmp(*i);
      if(tmp.lund() == 310){
        ks.push_back(Particle(*i));
      }
    }
  }
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
