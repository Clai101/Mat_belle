#include "my_belle.h"


#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


int code(int charg, int barion_num, int chanel){
    int code = 0;
    
    if (charg < 0)
    {
      code = code + 1 * 1e10;
    }
    code = code + abs(charg) * 1e9;
    
    if (barion_num < 0)
    {
      code = code + 1 * 1e8;
    }
    code = code + abs(barion_num) * 1e7;
    
    code = code + abs(chanel) * 1e4;
  
    return code;
}


void decode(int code, int& charg, int& barion_num, int& chanel) {
    
    int charg_sign = (code / static_cast<int>(1e10)) % 2; 
    charg = (code / static_cast<int>(1e9)) % static_cast<int>(1e1);

    int barion_sign = (code / static_cast<int>(1e8)) % 2;
    barion_num = (code / static_cast<int>(1e7)) % static_cast<int>(1e1);

    chanel = code % static_cast<int>(1e4);

    charg = (charg_sign == 1) ? -charg : charg;
    barion_num = (barion_sign == 1) ? -barion_num : barion_num;
}


template<typename T, typename... Args>
void make_comb(int charg, int barion_num, int chanel, std::vector<Particle> object, T type, Args... args)
{
    combination(new_part, object, args...);
    setUserInfo(new_part, code(charg, barion_num, chanel));
}

template<typename T, typename... Args>
void make_comb(int chanel, std::vector<Particle> object, T type, Args... args)
{
    combination(new_part, object, args...);
    setUserInfo(new_part, chanel);
}



using namespace std;
void User_reco::hist_def( void )
{ extern BelleTupleManager* BASF_Histogram;    
  t1 = BASF_Histogram->ntuple ("lmbda",
        "ml mach p chu chl chct en ecm ntr rm2n rm2l rm2nu chrgl chrgach chrgU");
  t2 = BASF_Histogram->ntuple ("sigma_check",
        "dm ml ms chl p chs chrgl chrgs abar ml");
};


int fill_tup(Particle lamc, /*vector<Particle> all,*/ double elec, double posi, double ecm, double r2, BelleTuple *t)
{    
  //int chb = dynamic_cast<UserInfo&>(B.userInfo()).channel();
  
  return 1;
};


void User_reco::event ( BelleEvent* evptr, int* status ) {
  
  *status=0;
  
  static int nevent=0;
  static int nwritt=0;
  if(++nevent<2 || !(nevent%1000)) cout << "Event number " << nevent
          << " selected " << nwritt << endl;
  
  Belle_runhead_Manager& rhdmgr = Belle_runhead_Manager::get_manager();
  Belle_runhead_Manager::const_iterator rhd = rhdmgr.begin();
  Evtcls_hadron_info_Manager&  ehimgr =
    Evtcls_hadron_info_Manager::get_manager();
  
  Evtcls_hadronic_flag_Manager&  ehadfl =
    Evtcls_hadronic_flag_Manager::get_manager();
  
  HepPoint3D ip_position = IpProfile::position();
  const HepSymMatrix& runIp_err = IpProfile::position_err();
  Mdst_vee2_Manager &vee2_mgr = Mdst_vee2_Manager::get_manager();
  Evtcls_hadron_info_Manager::iterator iti = ehimgr.begin();
  Evtcls_hadronic_flag_Manager::iterator eti = ehadfl.begin();
  
  double r2=0;
  int ntrk=0;
  double evis=0;
  double Pz=0;
  double hjmass=0;

  
  if (iti!=ehimgr.end()){
    r2 = (*iti).R2();
    //  ntrk = (*iti).Ntrk();
    //      evis = (*iti).Evis();
    //      Pz   = (*iti).Pz();
    //      hjmass = (*iti).HeavyJetMass();
  }
      
  double ecm = BeamEnergy::Ecm();
  double elec = BeamEnergy::E_HER();
  double posi = BeamEnergy::E_LER();
  VectorL beam =VectorL(elec*sin(0.022), 0.,elec*cos(0.022)-posi, elec+posi);
  
  /*************** Make particle lists ********************************/

  //Base particles

  std::vector<Particle> p, ap, k_p, k_m, pi_p, pi_m, pi0, gamma, all, e_p, e_m, mu_m, mu_p;

  //fill vectors
  makeProton(p, ap, 1);
  makeKPi(k_p, k_m, pi_p, pi_m, 1);
  makePi0(pi0);
  withEminCutPi0(pi0, 0.05);
  makeGamma(gamma);
  makeLepton(e_p, e_m, mu_p, mu_m, 1);

  //Cuts
  withDrDzCut(p, 1., 2.);
  withDrDzCut(e_m, 1., 2.);
  withDrDzCut(e_p, 1., 2.);
  withDrDzCut(mu_m, 1., 2.);
  withDrDzCut(mu_p, 1., 2.);
  withDrDzCut(k_m, 1., 2.);
  withDrDzCut(pi_p, 1., 2.);
  withDrDzCut(ap, 1., 2.);
  withDrDzCut(k_p, 1., 2.);
  withDrDzCut(pi_m, 1., 2.);
  withLeptonIdCut(e_p, e_m, mu_p, mu_m, 0.01, 0.1);
  withPSCut(mu_p, 1.);
  withPSCut(mu_m, 1.);
  withPSCut(e_p, 1.);
  withPSCut(e_m, 1.);

  withKaonIdCut(k_p, k_m, 0.6);
  withProtonIdCut(p, ap, 0.6);
  
  deepCopy(pi_p, all);
  deepCopy(pi_m, all);

  setGenHepInfoF(p);
  setGenHepInfoF(k_m);
  setGenHepInfoF(pi_p);
  
  //Undetected particles
  std::vector<Particle> lamc_p, lamc_m, lamct_p, lamct_m;
  std::vector<Particle> sigc_pp, sigc_mm, sigc0, asigc0;
  std::vector<Particle> lam, alam;
  std::vector<Particle> ups, rho, rho_2m, rho_2p, rho4, rho_ppm, rho_mmp;
  std::vector<Particle> D0, aD0, D_p, D_m;
  std::vector<Particle> k_s;
  
  combination(rho_mmp, m_ptypeRHO0, pi_m, pi_m, pi_p);
  setUserInfo(rho_mmp,  -1);
  combination(rho_ppm, m_ptypeRHO0, pi_p, pi_p, pi_m);
  setUserInfo(rho_ppm,  1);

  combination(rho, m_ptypeRHO0, pi_p, pi_m);
  setUserInfo(rho,  11);
  combination(rho4, m_ptypeRHO0, rho_2p, rho_2m);
  setUserInfo(rho4,  13);
  combination(rho_2p, m_ptypeRHO0, pi_p, pi_p);
  setUserInfo(rho_2p,  12);
  combination(rho_2m, m_ptypeRHO0, pi_m, pi_m);
  setUserInfo(rho_2m,  -12);

  makeKs(k_s);
  makeLambda(lam,alam);

  for(std::vector<Particle>::iterator l = k_s.begin(); l!=k_s.end(); ++l) {
    HepPoint3D V(l->mdstVee2().vx(),l->mdstVee2().vy(),0);
    Vector3 P(l->px(),l->py(),0);
    V=V-ip_position;
    V.setZ(0.);
    if (abs(l->mass()-0.4977)>0.03 || V.perp()<0.1 ||
    cos(V.angle(P))<0.99 || l->mdstVee2().z_dist()>1.) {
      k_s.erase(l); --l; continue;
    }
  }

  for(std::vector<Particle>::iterator l = lam.begin(); l!=lam.end(); ++l) {
      HepPoint3D V(l->mdstVee2().vx(),l->mdstVee2().vy(),0);
    Vector3 P(l->px(),l->py(),0);
    V=V-ip_position;
    V.setZ(0.);
    double p_id;
    if (l->child(0).pType().mass()>l->child(1).pType().mass())
      p_id=atc_pid(3,-1,5,4,2).prob(&(l->child(0).mdstCharged()));
    else p_id=atc_pid(3,-1,5,4,2).prob(&(l->child(1).mdstCharged()));
    if (abs(l->mass()-1.1157)>0.01 || l->mdstVee2().z_dist()>1. || p_id<0.6 ) {
      lam.erase(l); --l;
    }
  }


  for(std::vector<Particle>::iterator l = alam.begin(); l!=alam.end(); ++l) {
    HepPoint3D V(l->mdstVee2().vx(),l->mdstVee2().vy(),0);
    Vector3 P(l->px(),l->py(),0);
    V=V-ip_position;
    V.setZ(0.);
    double p_id;
    if (l->child(0).pType().mass()>l->child(1).pType().mass()) 
p_id=atc_pid(3,-1,5,4,2).prob(&(l->child(0).mdstCharged()));
    else p_id=atc_pid(3,-1,5,4,2).prob(&(l->child(1).mdstCharged()));                               
    if (abs(l->mass()-1.1157)>0.01 || l->mdstVee2().z_dist()>1. || p_id<0.6 ) {
      alam.erase(l); --l;
    }
  }

  /*D*/

  combination(D_p, m_ptypeDP, k_m, pi_p, pi_p, 0.05);
  setUserInfo(D_p,  1);
  combination(D_m, m_ptypeDM, k_p, pi_m, pi_m, 0.05);
  setUserInfo(D_m,  -1);

  combination(D_p, m_ptypeD0, k_s, pi_p, 0.05);
  setUserInfo(D_p, 2);
  combination(D_m, m_ptypeD0B, k_s, pi_m, 0.05);
  setUserInfo(D_m, -2);

  combination(D0, m_ptypeD0, k_m, pi_p, 0.05);
  combination(aD0, m_ptypeD0B, k_p, pi_m, 0.05);
  setUserInfo(D0,  11);
  setUserInfo(aD0,  11);

  combination(D0, m_ptypeD0, k_m, rho_2p, pi_m, 0.05);
  combination(aD0, m_ptypeD0B, k_p, rho_2m, pi_p, 0.05);
  setUserInfo(D0,  12);
  setUserInfo(aD0,  12);

  combination(D0, m_ptypeD0, k_m, pi0, pi_p, 0.05);
  combination(aD0, m_ptypeD0B, k_p, pi0, pi_m, 0.05);
  setUserInfo(D0, 13);
  setUserInfo(aD0, 13);

  combination(D0, m_ptypeD0, k_m, k_p, 0.05);
  combination(aD0, m_ptypeD0B, k_p, k_m, 0.05);
  setUserInfo(D0, 14);
  setUserInfo(aD0, 14);

  combination(D0, m_ptypeD0, k_s, pi_p, pi_m, 0.05);
  combination(aD0, m_ptypeD0B, k_s, pi_p, pi_m, 0.05);
  setUserInfo(D0, 15);
  setUserInfo(aD0, 15);

  combination(D0, m_ptypeD0, k_s, pi0, 0.05);
  combination(aD0, m_ptypeD0B, k_s, pi0, 0.05);
  setUserInfo(D0, 16);
  setUserInfo(aD0, 16);



  /*
  Lambda mesons
  */

  std::cout << "first negative usinf" << endl;

  combination(lamct_m, m_ptypeLAMC, ap, k_p, pi_m, 0.05);
  setUserInfo(lamct_m,  -1);

  combination(lamct_p, m_ptypeLAMC, p, k_m, pi_p, 0.05);
  setUserInfo(lamct_p,  1);

  combination(lamct_m, m_ptypeLAMC, alam, pi_m, 0.05);
  setUserInfo(lamct_m,  -2);

  combination(lamct_p, m_ptypeLAMC, lam, pi_p, 0.05);
  setUserInfo(lamct_p,  2);

  combination(lamct_m, m_ptypeLAMC, ap, k_s, 0.05);
  setUserInfo(lamct_m,  -3);

  combination(lamct_p, m_ptypeLAMC, p, k_s, rho, 0.05);
  setUserInfo(lamct_p,  3);

  combination(lamct_m, m_ptypeLAMC, ap, k_s, rho, 0.05);
  setUserInfo(lamct_m,  -4);

  combination(lamct_p, m_ptypeLAMC, p, k_s, 0.05);
  setUserInfo(lamct_p,  4);



  /*
  semileptonics mode
  */

  /*combination(lamc_p, m_ptypeLAMC, lam, e_p);
  combination(lamc_m, m_ptypeLAMC, alam, e_m);
  setUserInfo(lamc_p,  10);
  setUserInfo(lamc_m,  10);

  combination(lamc_p, m_ptypeLAMC, lam, mu_p);
  combination(lamc_m, m_ptypeLAMC, alam, mu_m);
  setUserInfo(lamc_p,  11);
  setUserInfo(lamc_m,  11);*/

  combination(lamc_p, m_ptypeLAMC, lam, pi_p, 0.05);
  setUserInfo(lamc_p,  2);

  combination(lamc_m, m_ptypeLAMC, alam, pi_m, 0.05);
  setUserInfo(lamc_m,  -2);    

  combination(lamc_p, m_ptypeLAMC, p, k_m, pi_p, 0.05);
  setUserInfo(lamc_p,  3);

  combination(lamc_m, m_ptypeLAMC, ap, k_p, pi_m, 0.05);
  setUserInfo(lamc_m,  -3);

  /*
    Sigma_c
  */
  
  combination(sigc_mm, m_ptypeSIGC0, lamc_m, pi_m);
  setUserInfo(sigc_mm,  -2);

  combination(sigc_pp, m_ptypeSIGC0, lamc_p, pi_p);
  setUserInfo(sigc_pp,  2);

  combination(sigc0, m_ptypeSIGC0, lamc_p, pi_m);
  setUserInfo(sigc0,  1);

  combination(asigc0, m_ptypeSIGC0, lamc_m, pi_p);
  setUserInfo(asigc0,  -1);

  for(std::vector<Particle>::iterator l = sigc_pp.begin(); l!=sigc_pp.end(); ++l) {
    if (l->mass()-l->child(0).mass()>0.18) {
      sigc_pp.erase(l); --l;
    }
  }
  for(std::vector<Particle>::iterator l = sigc_mm.begin(); l!=sigc_mm.end(); ++l) {
    if (l->mass()-l->child(0).mass()>0.18) {
      sigc_mm.erase(l); --l;
    }
  }
  for(std::vector<Particle>::iterator l = sigc0.begin(); l!=sigc0.end(); ++l) {
    if (l->mass()-l->child(0).mass()>0.18) {
      sigc0.erase(l); --l;
    }
  }
  for(std::vector<Particle>::iterator l = asigc0.begin(); l!=asigc0.end(); ++l) {
    if (l->mass()-l->child(0).mass()>0.18) {
      asigc0.erase(l); --l;
    }
  }


  /*
    All together    
  */
  
  combination(ups, m_ptypeUPS4, lamc_p, lamct_m);
  combination(ups, m_ptypeUPS4, lamc_m, lamct_p);
  setUserInfo(ups, 1);

  combination(ups, m_ptypeUPS4, lamc_p, lamct_m, rho);
  combination(ups, m_ptypeUPS4, lamc_m, lamct_p, rho);
  setUserInfo(ups, 2);

  combination(ups, m_ptypeUPS4, lamc_p, lamct_m, rho4);
  combination(ups, m_ptypeUPS4, lamc_m, lamct_p, rho4);
  setUserInfo(ups, 3);

  combination(ups, m_ptypeUPS4, lamc_m, D0, p);
  combination(ups, m_ptypeUPS4, lamc_p, aD0, ap);
  setUserInfo(ups, 4);

  combination(ups, m_ptypeUPS4, lamc_m, D0, p, rho);
  combination(ups, m_ptypeUPS4, lamc_p, aD0, ap, rho);
  setUserInfo(ups, 5);
  
  combination(ups, m_ptypeUPS4, lamc_m, D_p, p, pi_m);
  combination(ups, m_ptypeUPS4, lamc_p, D_m, ap, pi_p);
  setUserInfo(ups,  6);
  
  combination(ups, m_ptypeUPS4, lamc_m, D_p, p, rho_mmp);
  combination(ups, m_ptypeUPS4, lamc_p, D_m, ap, rho_ppm);
  setUserInfo(ups,  7);


  std::vector<Particle> sigc;

  deepCopy(sigc0, sigc);
  deepCopy(asigc0, sigc);
  deepCopy(sigc_mm, sigc);
  deepCopy(sigc_pp, sigc);



  for(int i=0; i<sigc.size(); ++i){

    Particle s = sigc[i];
    Particle l = s.child(0);
    int chs = dynamic_cast<UserInfo&>(s.userInfo()).channel();
    int chl = dynamic_cast<UserInfo&>(l.userInfo()).channel();

    t2->column("ml", l.mass());
    t2->column("ms", s.mass());
    t2->column("dm", s.mass() - l.mass() + l.pType().mass());     
    t2->column("chrgl", chl / abs(chl));     
    t2->column("chrgs", chs / abs(chs));     
    t2->column("chs", abs(chs));
    t2->column("chl", abs(chl));
    t2->column("p", pStar(s, elec, posi).vect().mag());
    t2->dumpData();

  }
    
  for(int j=0; j<ups.size(); ++j){
    Particle u=ups[j];

    short ntr=0;
    short n = u.nChildren();
    for(int jj=0; jj<all.size(); ++jj) 
    if (!checkSame(all[jj],u)) ntr++;

    Particle lam = u.child(0);
    Particle ach = u.child(1);
    int chl = dynamic_cast<UserInfo&>(lam.userInfo()).channel();
    int chu = dynamic_cast<UserInfo&>(u.userInfo()).channel();
    int chargU = 0;
    


    for (int k = 0; k < n; ++k){
      int p_charge = 0;
      if ((chu >= 4 & k == 2) || (chu == 6 & k == 3)){
        p_charge = u.child(k).charge();
      }
      else{
        p_charge = dynamic_cast<UserInfo&>(u.child(k).userInfo()).channel() / 100;
      }
      if (p_charge != 0){
        chargU = chargU + p_charge / abs(p_charge);
      }
    }

    t1->column("ml", lam.mass());
    t1->column("mach", ach.mass() - ach.pType().mass());
    t1->column("p", pStar(u, elec, posi).vect().mag());
    t1->column("chu", chu);     
    t1->column("chl", chl);
    t1->column("chct", dynamic_cast<UserInfo&>(ach.userInfo()).channel());
    t1->column("en", pStar(u, elec, posi).e());
    t1->column("ecm", ecm);
    t1->column("ntr", ntr);
    t1->column("chrgl", chl / abs(chl));     
    t1->column("chrgach", chu / abs(chu));
    t1->column("chrgU", chargU);
    t1->column("rm2n", (beam - u.p()).m2());
    t1->column("rm2l", (beam - (u.p() - lam.p())).m2());
    t1->column("rm2nu", (beam - (u.p() - lam.p() + lam.child(0).p() + lam.child(1).p())).m2());

    t1->dumpData();


*status = 1;
}

if (*status==1) nwritt++;

}

#if defined(BELLE_NAMESPACE)
}
#endif