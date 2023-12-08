#include "my_belle.h"


#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

using namespace std;
void User_reco::hist_def( void )
{ extern BelleTupleManager* BASF_Histogram;    
  t1 = BASF_Histogram->ntuple ("lmbda_lept",
    "ml mach p chxc chl __chxc ch_dot_ach chach en dm_dst ecm ntr rm2n fox rm2l rm2nu chrgl chrgach en_gam enc_gam pl");
  /*
  t2 = BASF_Histogram->ntuple ("lmbdat",
    "en ecm p ntr chxc chrgach chach mach rm2lc");
  */
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

  std::vector<Particle> p, ap, k_p, k_m, k0, pi_p, pi_m, pi0, gamma, gam, all, e_p, e_m, mu_m, mu_p, k_s, lam, alam;

  //fill vectors
  makeLambda(lam,alam);
  makeProton(p, ap, 1);

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
  withProtonIdCut(p, ap, 0.6);
  withDrDzCut(p, 1., 2.);
  withDrDzCut(ap, 1., 2.);

  if((p.size() + ap.size() + lam.size() + alam.size() < 1.5)) return;


  makeKs(k_s);
  makeKPi(k_p, k_m, pi_p, pi_m, 1);


  //Cuts
  withDrDzCut(e_m, 1., 2.);
  withDrDzCut(e_p, 1., 2.);
  withDrDzCut(mu_m, 1., 2.);
  withDrDzCut(mu_p, 1., 2.);
  withDrDzCut(k_m, 1., 2.);
  withDrDzCut(k_p, 1., 2.);
  withDrDzCut(pi_p, 1., 2.);
  withDrDzCut(pi_m, 1., 2.);

  withKaonIdCut(k_p, k_m, 0.6);

  
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


  //static int n_trash_track = 0, n_trash_str = 0, n_trash_bar = 0, n_trash_str_bar = 0;
  int n_str = k_s.size() + k_m.size() + k_p.size() + lam.size() + alam.size();
  if (n_str < 1.5) return;
  if (pi_p.size() + pi_m.size() > 12.5) return;
  



  makeLepton(e_p, e_m, mu_p, mu_m, 1);
  withLeptonIdCut(e_p, e_m, mu_p, mu_m, 0.01, 0.1);
  
  makePi0(pi0);
  //makeK0(k0, ak0);
  withEminCutPi0(pi0, 0.05);
  //withEminCutK0(k0, 0.05);
  makeGamma(gamma);

  setGenHepInfoF(p);
  setGenHepInfoF(k_m);
  setGenHepInfoF(pi_p);

  /*
  if (((pi_p.size() != pi_m.size()) && (k_s.size() + lam.size() + alam.size())) || (pi_p.size() + pi_m.size() > 12)) cout << "Happened" << endl;
  */
  //Undetected particles

  std::vector<Particle> ups;
  
  //Lambdac  

  std::vector<Particle> lamc_p, lamc_m;

  combination(lamc_p, m_ptypeLAMC, lam, e_p);
  combination(lamc_m, m_ptypeLAMC, alam, e_m);
  setUserInfo(lamc_p,  {{"chanel", 1}, {"charg", 1}, {"baryon_num", -1}});
  setUserInfo(lamc_m,  {{"chanel", 1}, {"charg", -1}, {"baryon_num", 1}});

  combination(lamc_p, m_ptypeLAMC, lam, mu_p);
  combination(lamc_m, m_ptypeLAMC, alam, mu_m);
  setUserInfo(lamc_p,  {{"chanel", 2}, {"charg", 1}, {"baryon_num", -1}});
  setUserInfo(lamc_m,  {{"chanel", 2}, {"charg", -1}, {"baryon_num", 1}});

  combination(lamc_m, m_ptypeLAMC, alam, pi_m, 0.05);
  combination(lamc_p, m_ptypeLAMC, lam, pi_p, 0.05);
  setUserInfo(lamc_m,  {{"chanel", 3}, {"charg", -1}, {"baryon_num", 1}});
  setUserInfo(lamc_p,  {{"chanel", 3}, {"charg", 1}, {"baryon_num", -1}});

  combination(lamc_m, m_ptypeLAMC, alam, pi_m, pi0, 0.05);
  combination(lamc_p, m_ptypeLAMC, lam, pi_p, pi0, 0.05);
  setUserInfo(lamc_m,  {{"chanel", 4}, {"charg", -1}, {"baryon_num", 1}});
  setUserInfo(lamc_p,  {{"chanel", 4}, {"charg", 1}, {"baryon_num", -1}});

  combination(lamc_m, m_ptypeLAMC, ap, k_p, pi_m, 0.05);
  combination(lamc_p, m_ptypeLAMC, p, k_m, pi_p, 0.05);
  setUserInfo(lamc_m, {{"chanel", 5}, {"charg", -1}, {"baryon_num", -1}});
  setUserInfo(lamc_p, {{"chanel", 5}, {"charg", 1}, {"baryon_num", 1}});

  combination(lamc_m, m_ptypeLAMC, ap, k_p, pi_m, 0.05);
  combination(lamc_p, m_ptypeLAMC, p, k_m, pi_p, 0.05);
  setUserInfo(lamc_m, {{"chanel", 5}, {"charg", -1}, {"baryon_num", -1}});
  setUserInfo(lamc_p, {{"chanel", 5}, {"charg", 1}, {"baryon_num", 1}});
  
  if (lamc_p.size()+lamc_m.size()==0) return;


  //D_pm

  std::vector<Particle> D_p, D_m;

  combination(D_p, m_ptypeDP, k_m, pi_p, pi_p, 0.05);
  setUserInfo(D_p,  {{"chanel", 1}, {"charg", 1}, {"baryon_num", 0}});
  combination(D_m, m_ptypeDM, k_p, pi_m, pi_m, 0.05);
  setUserInfo(D_m,  {{"chanel", 1}, {"charg", -1}, {"baryon_num", 0}});

  combination(D_p, m_ptypeDP, k_s, pi_p, 0.05);
  setUserInfo(D_p, {{"chanel", 2}, {"charg", 1}, {"baryon_num", 0}});
  combination(D_m, m_ptypeDM, k_s, pi_m, 0.05);
  setUserInfo(D_m, {{"chanel", 2}, {"charg", -1}, {"baryon_num", 0}});

  combination(D_p, m_ptypeDP, k_s, pi_p, pi_p, pi_m, 0.05);
  setUserInfo(D_p, {{"chanel", 3}, {"charg", 1}, {"baryon_num", 0}});
  combination(D_m, m_ptypeDM, k_s, pi_m, pi_m, pi_p, 0.05);
  setUserInfo(D_m, {{"chanel", 3}, {"charg", -1}, {"baryon_num", 0}});

  combination(D_p, m_ptypeDP, k_p, k_m, pi_p, 0.05);
  setUserInfo(D_p, {{"chanel", 4}, {"charg", 1}, {"baryon_num", 0}});
  combination(D_m, m_ptypeDM, k_p, k_m, pi_m, 0.05);
  setUserInfo(D_m, {{"chanel", 4}, {"charg", -1}, {"baryon_num", 0}});

  //D0

  std::vector<Particle> D0, aD0;

  combination(D0, m_ptypeD0, k_m, pi_p, 0.05);
  combination(aD0, m_ptypeD0B, k_p, pi_m, 0.05);
  setUserInfo(D0, {{"chanel", 1}, {"charg", 0}, {"baryon_num", 0}});
  setUserInfo(aD0, {{"chanel", 1}, {"charg", 0}, {"baryon_num", 0}});

  combination(D0, m_ptypeD0, k_m, pi_p, pi_p, pi_m, 0.05);
  combination(aD0, m_ptypeD0B, k_p, pi_m, pi_m, pi_p, 0.05);
  setUserInfo(D0,  {{"chanel", 2}, {"charg", 0}, {"baryon_num", 0}});
  setUserInfo(aD0,  {{"chanel", 2}, {"charg", 0}, {"baryon_num", 0}});

  combination(D0, m_ptypeD0, k_m, pi0, pi_p, 0.05);
  combination(aD0, m_ptypeD0B, k_p, pi0, pi_m, 0.05);
  setUserInfo(D0, {{"chanel", 3}, {"charg", 0}, {"baryon_num", 0}});
  setUserInfo(aD0, {{"chanel", 3}, {"charg", 0}, {"baryon_num", 0}});

  combination(D0, m_ptypeD0, k_m, k_p, 0.05);
  combination(aD0, m_ptypeD0B, k_p, k_m, 0.05);
  setUserInfo(D0, {{"chanel", 4}, {"charg", 0}, {"baryon_num", 0}});
  setUserInfo(aD0, {{"chanel", 4}, {"charg", 0}, {"baryon_num", 0}});

  combination(D0, m_ptypeD0, k_s, pi_p, pi_m, 0.05);
  combination(aD0, m_ptypeD0B, k_s, pi_p, pi_m, 0.05);
  setUserInfo(D0, {{"chanel", 5}, {"charg", 0}, {"baryon_num", 0}});
  setUserInfo(aD0, {{"chanel", 5}, {"charg", 0}, {"baryon_num", 0}});

  combination(D0, m_ptypeD0, k_s, pi0, 0.05);
  combination(aD0, m_ptypeD0B, k_s, pi0, 0.05);
  setUserInfo(D0, {{"chanel", 6}, {"charg", 0}, {"baryon_num", 0}});
  setUserInfo(aD0, {{"chanel", 6}, {"charg", 0}, {"baryon_num", 0}});

  combination(D0, m_ptypeD0, k_m, k_p, 0.05);
  combination(aD0, m_ptypeD0B,  k_m, k_p, 0.05);
  setUserInfo(D0, {{"chanel", 7}, {"charg", 0}, {"baryon_num", 0}});
  setUserInfo(aD0, {{"chanel", 7}, {"charg", 0}, {"baryon_num", 0}});

  //D star

  std::vector<Particle> D0_st, aD0_st;

  combination(D0_st, m_ptypeDstar0, D0, pi0, 0.05);
  combination(aD0_st, m_ptypeDstarB, aD0, pi0, 0.05);
  setUserInfo(D0_st, {{"chanel", 1}, {"charg", 0}, {"baryon_num", 0}});
  setUserInfo(aD0_st, {{"chanel", 1}, {"charg", 0}, {"baryon_num", 0}});

  combination(D0_st, m_ptypeDstar0, D0, gamma, 0.05);
  combination(aD0_st, m_ptypeDstarB, aD0, gamma, 0.05);
  setUserInfo(D0_st, {{"chanel", 2}, {"charg", 0}, {"baryon_num", 0}});
  setUserInfo(aD0_st, {{"chanel", 2}, {"charg", 0}, {"baryon_num", 0}});

  std::vector<Particle> D_p_st, D_m_st;

  combination(D_p_st, m_ptypeDstar0, D0, pi_p, 0.05);
  combination(D_m_st, m_ptypeDstarB, aD0, pi_m, 0.05);
  setUserInfo(D_p_st, {{"chanel", 1}, {"charg", 1}, {"baryon_num", 0}});
  setUserInfo(D_m_st, {{"chanel", 1}, {"charg", -1}, {"baryon_num", 0}});

  combination(D_p_st, m_ptypeDstar0, D_p, pi0, 0.05);
  combination(D_m_st, m_ptypeDstarB, D_m, pi0, 0.05);
  setUserInfo(D_p_st, {{"chanel", 2}, {"charg", 1}, {"baryon_num", 0}});
  setUserInfo(D_m_st, {{"chanel", 2}, {"charg", -1}, {"baryon_num", 0}});

  //cout << "Before " << D0_st.size() << '\n';

  for(std::vector<Particle>::iterator D = D0_st.begin(); D!=D0_st.end(); ++D) {
    if (dynamic_cast<UserInfo&>(D->userInfo()).channel().find("chanel")->second == 1){
      if (D->mass() - D->child(0).mass() > 0.16 or abs(D->child(0).mass() - D->child(0).pType().mass()) > 0.015){
      aD0_st.erase(D); --D; continue;}}
    if (dynamic_cast<UserInfo&>(D->userInfo()).channel().find("chanel")->second == 2){
      if (abs(D->mass() - D->child(0).mass() - 0.142014) > 0.030 or abs(D->child(0).mass() - D->child(0).pType().mass()) > 0.015){
        aD0_st.erase(D); --D; continue;}}}

  for(std::vector<Particle>::iterator D = aD0_st.begin(); D!=aD0_st.end(); ++D) {
    if (dynamic_cast<UserInfo&>(D->userInfo()).channel().find("chanel")->second == 1){
      if (D->mass() - D->child(0).mass() > 0.16 or abs(D->child(0).mass() - D->child(0).pType().mass()) > 0.015){
      aD0_st.erase(D); --D; continue;}}
    if (dynamic_cast<UserInfo&>(D->userInfo()).channel().find("chanel")->second == 2){
      if (abs(D->mass() - D->child(0).mass() - 0.142014) > 0.030 or abs(D->child(0).mass() - D->child(0).pType().mass()) > 0.015){
        aD0_st.erase(D); --D; continue;}}}

  for(std::vector<Particle>::iterator D = D_p_st.begin(); D!=D_p_st.end(); ++D) {
    if (D->mass() - D->child(0).mass() > 0.16 or abs(D->child(0).mass() - D->child(0).pType().mass()) > 0.015) {
      D_p_st.erase(D); --D; continue;}}

  for(std::vector<Particle>::iterator D = D_m_st.begin(); D!=D_m_st.end(); ++D) {
    if (D->mass() - D->child(0).mass() > 0.16 or abs(D->child(0).mass() - D->child(0).pType().mass()) > 0.015) {
      D_m_st.erase(D); --D; continue;}}

  //X_c


  std::vector<Particle> X_c, aX_c;

  combination(X_c, m_ptypeUPS4, D0, p);
  combination(aX_c, m_ptypeUPS4, aD0, ap);
  setUserInfo(X_c, {{"chanel", 1}, {"charg", 1}, {"baryon_num", 1}});
  setUserInfo(aX_c, {{"chanel", 1}, {"charg", -1},  {"baryon_num", -1}});


  combination(X_c, m_ptypeUPS4, D_p, p, pi_m);
  combination(aX_c, m_ptypeUPS4, D_m, ap, pi_p);
  setUserInfo(X_c, {{"chanel", 2}, {"charg", 1}, {"baryon_num", 1}});
  setUserInfo(aX_c, {{"chanel", 2}, {"charg", -1},  {"baryon_num", -1}});

/*
  combination(X_c, m_ptypeUPS4, D0, k_p, lam);
  combination(aX_c, m_ptypeUPS4, aD0, k_m, lam);
  setUserInfo(X_c, {{"chanel", 2}, {"charg", 1}, {"baryon_num", 1}});
  setUserInfo(aX_c, {{"chanel", 2}, {"charg", -1},  {"baryon_num", -1}});
*/

  combination(X_c, m_ptypeUPS4, D_p_st, p, pi_m);
  combination(aX_c, m_ptypeUPS4, D_m_st, ap, pi_p);
  setUserInfo(X_c,  {{"chanel", 3}, {"charg", 1}, {"baryon_num", 1}});
  setUserInfo(aX_c,  {{"chanel", 3}, {"charg", -1}, {"baryon_num", -1}});

  combination(X_c, m_ptypeUPS4, D0_st, p);
  combination(aX_c, m_ptypeUPS4, aD0_st, ap);
  setUserInfo(X_c,  {{"chanel", 4}, {"charg", 1}, {"baryon_num", 1}});
  setUserInfo(aX_c,  {{"chanel", 4}, {"charg", -1}, {"baryon_num", -1}});

  combination(ups, m_ptypeUPS4, lamc_m, X_c);
  combination(ups, m_ptypeUPS4, lamc_p, aX_c);
  

  for(int j=0; j<ups.size(); ++j){
    Particle u=ups[j];

    short ntr=0;
    for(int jj=0; jj<pi_p.size(); ++jj) 
    if (!checkSame(pi_p[jj],u)) ntr++;

    for(int jj=0; jj<pi_m.size(); ++jj) 
    if (!checkSame(pi_m[jj],u)) ntr++;

    

    float en_gam=0, enc_gam=0;
    for(int jj=0; jj<gamma.size(); ++jj){
    if (!checkSame(gamma[jj],u)){
      enc_gam = enc_gam + pStar(gamma[jj], elec, posi).e();
      en_gam = en_gam + (gamma[jj]).p().e();
    }} 

    Particle lam = u.child(0);
    Particle ach = u.child(1).child(0);
    Particle dot_ach = ach.child(0);
    Particle x_c = u.child(1);
    map <string, int> chl = dynamic_cast<UserInfo&>(lam.userInfo()).channel();
    map <string, int> chach = dynamic_cast<UserInfo&>(ach.userInfo()).channel();
    map <string, int> chxc = dynamic_cast<UserInfo&>(x_c.userInfo()).channel();
    float delta_dst = -400;
    int ch_dot_ach = 0;

    UserInfo __chxc = static_cast<UserInfo&>(x_c.userInfo());
    
    if (__chxc.channel().find("chanel")->second >= 3){
      map <string, int> _ = dynamic_cast<UserInfo&>(dot_ach.userInfo()).channel();
      ch_dot_ach = _["chanel"];
      delta_dst =  ach.mass() - ach.child(0).mass();
    }


    if ((beam - (x_c.p())).m2() > 4*4 or (ntr >= 1)) continue;    
    
    int chargU = 0;


    t1->column("en", pStar(u, elec, posi).e());
    t1->column("ecm", ecm);
    t1->column("p", pStar(u, elec, posi).vect().mag());
    t1->column("ntr", ntr);

    t1->column("en_gam", en_gam);
    t1->column("enc_gam", enc_gam);

    t1->column("chxc", chxc["chanel"]);     
    t1->column("__chxc", __chxc.channel().find("chanel")->second);     

    t1->column("chach", chach["chanel"]);
    t1->column("chrgach", chach["charg"]);
    t1->column("mach", ach.mass() - ach.pType().mass());
    t1->column("dm_dst", delta_dst);

    t1->column("ch_dot_ach", ch_dot_ach);

    t1->column("chrgl", chl["charg"]);     
    t1->column("chl", chl["chanel"]);
    t1->column("ml", lam.mass() - 2.28646);

    t1->column("rm2n", (beam - u.p()).m2());
    t1->column("rm2l", (beam - (x_c.p())).m2());
    t1->column("rm2nu", (beam - (x_c.p() + lam.child(0).p() + lam.child(1).p())).m2());

    t1->column("fox", r2);

    if(chl["chanel"] <= 2){
      t1->column("pl", pStar(lam.child(1), elec, posi).vect().mag());}
    else{
      t1->column("pl", 0);}
    

    t1->dumpData();;
    *status = 1; 
}

if (*status==1) {nwritt++;
  cout << "Chac " << pi_p.size() << " " <<  pi_m.size() << " ntrack " << pi_p.size() + pi_m.size() << " ks_size " << k_s.size() << " lam size " << lam.size() + alam.size() <<endl;
}
}

#if defined(BELLE_NAMESPACE)
}
#endif