#include "reco.h"
#include "userinfo.h"
#include "myutils.h"
// #include "genhep.h"
#include "helix/Helix.h"
#include EVTCLS_H //R2 distribution
#include "toolbox/Thrust.h"
#include "toolbox/FuncPtr.h"
#include "psi.h"
#include "cstdlib"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

using namespace std;

const double rand_m=abs(double(RAND_MAX+1));
const double dwid[] = 
  {0.015, 0.015, 0.01, 0.020, 0.015, 0.015, 0.02, 0.015, 0.015,
   0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015}; 
//  Kpi,   KK,   K3pi, K2pi, Ks2pi,  Ks2K,  K2pi,  Kspi,  KKpi
const double dwidst=0.003;

double dmass_best[]={10,10,10,10,10,10,10,10,10,10};
double dmass_best1[]={10,10,10,10,10,10,10,10,10,10};
double dmass_best2[]={10,10,10,10,10,10,10,10,10,10};
  
int idbest[]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
int idbest1[]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
int idbest2[]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

void User_reco::hist_def( void )
{ extern BelleTupleManager* BASF_Histogram;
 t1 = BASF_Histogram->ntuple (" Inclusive psi ", 
			      " mj bin pj rm lhe lhm "
			      " lv cth thr run ecm "
			      " ea r2 nt ne nl em f2 id ");
 t2 = BASF_Histogram->ntuple (" psi D ",
			      " mj bin pj rm lhe lhm "
			      " md pd vd rmx lvd0 bin1 chd "
			      " lv cth thr run ecm "
			      " mph mks "
			      " b1 b2 b3 mx lvds lvd mpi "
			      " rmxd rmxds rmfd rmfds rmd "
			      " mftd mftds "
			      " hepid hepidm "
			      " ea r2 nt ne nl em f2 id ");
 t3 = BASF_Histogram->ntuple (" psi D D ",
			      " mj bin pj rm lhe lhm "
			      " md1 pd1 bin1 chd1 "
			      " md2 pd2 bin2 chd2 chx "
			      " rmx mx mdd cthx px "
			      " ea r2 nt ne nl em f2 id ");
 t4 = BASF_Histogram->ntuple (" psi2 D ",
			      " mj bin pj rm md2 "
			      " lv cth dvz "
			      " md pd vd rmx bin1 chd "
			      " mx lvds lvd mpi "
			      " rmxd rmxds rmfd rmfds rmd "
			      " ea r2 nt ne nl em f2 id ");
 t5 = BASF_Histogram->ntuple (" psi2 Dst ",
			      " mp2 mj bin pj rm lhe lhm vp "
			      " lv cth dvz "
			      " md pd vd rmx bin1 bin2 chd "
			      " b1 b2 b3 mds pds chds mx "
			      " rmxd rmxds rmfd rmfds rmd mpi "
			      " hepid hepidm "
			      " mftd mftds "
			      " ea r2 nt ne nl em f2 ");

 t6 = BASF_Histogram->ntuple (" psi MC ",
			      " mj pj rm lv cth ndp " );

 t7 = BASF_Histogram->ntuple (" psi D MC ",
			      " mj pj rm lv cth " 
			      " id rmx mx pd " );

 t8 = BASF_Histogram->ntuple (" Inclusive psi2 ", 
			      " mp2 mj bin pj rm rm2 p2 lhe lhm"
			      " lv cth thr run ecm "
			      " ea r2 nt ne nl em f2 id ");
};

void init_best (void) {

  for(int im=0; im<10; ++im) {
    dmass_best[im]=10;
    dmass_best1[im]=10;
    dmass_best2[im]=10;
  
    idbest[im]=-1;
    idbest1[im]=-1;
    idbest2[im]=-1;
  }
}
void 
fitRM3(const Particle &p1, const Particle &p2, const Particle &p3, const VectorL UPS, const double m_fit, 
     VectorL &P_p1, VectorL &P_p2, VectorL &P_p3, double &chi2) {

  int ierr;
  double lambda=0;
  
  P_p1=p1.p();
  P_p2=p2.p();
  P_p3=p3.p();
  
  const double M_p1=p1.pType().mass();
  const double M_p2=p2.pType().mass();
  const double M_p3=p3.pType().mass();
  
  HepSymMatrix err_p1(3,0), err_p2(3,0), err_p3(3,0), unit(3,0);
  for(int im=0; im<3; ++im) {
    unit[im][im]=1;
    for(int in=0; in<3; ++in) { 
      err_p1[im][in]=p1.momentum().dp()[im][in];
      err_p2[im][in]=p2.momentum().dp()[im][in];
      err_p3[im][in]=p3.momentum().dp()[im][in];
      if (abs(err_p1[im][in])<1.e-10) err_p1[im][in]=1.e-10;
      if (abs(err_p2[im][in])<1.e-10) err_p2[im][in]=1.e-10;
      if (abs(err_p3[im][in])<1.e-10) err_p3[im][in]=1.e-10;
    }
  }
  
  HepVector v_p1(3,0), v_p2(3,0), v_p3(3,0);
  v_p1[0]=P_p1.x();
  v_p1[1]=P_p1.y();
  v_p1[2]=P_p1.z();
  
  v_p2[0]=P_p2.x();
  v_p2[1]=P_p2.y();
  v_p2[2]=P_p2.z();
  
  v_p3[0]=P_p3.x();
  v_p3[1]=P_p3.y();
  v_p3[2]=P_p3.z();
  
  HepVector v_p10(3,0), v_p20(3,0), v_p30(3,0);
  v_p10=v_p1;
  v_p20=v_p2;
  v_p30=v_p3;

  HepSymMatrix A(10,0);
  HepVector a(10,0), diff_a(10,0);

  HepVector diff_p1(3,0), diff_p2(3,0), diff_p3(3,0);
  
  HepSymMatrix diff2_p1(3,0), diff2_p2(3,0), diff2_p3(3,0);
  HepSymMatrix diff2_p1_p2(3,0), diff2_p1_p3(3,0), diff2_p2_p3(3,0);
  
  double diff;
  
  HepVector diff_chi_p1(3,0), diff_chi_p2(3,0), diff_chi_p3(3,0);
  HepVector diff_m_p1(3,0), diff_m_p2(3,0), diff_m_p3(3,0);
  HepVector diff_lam_p1(3,0), diff_lam_p2(3,0), diff_lam_p3(3,0);
    
  for(int itt=0; itt<30; ++itt) {
    
    double E_recoil=(UPS-P_p1-P_p2-P_p3).e();
    double M_recoil=(UPS-P_p1-P_p2-P_p3).m();

    double E_p1=P_p1.e();
    double E_p2=P_p2.e();
    double E_p3=P_p3.e();

    diff=M_recoil* M_recoil - m_fit* m_fit;
    if (abs(diff)<0.002) continue;

    diff_chi_p1 = 2. * err_p1.inverse(ierr) *  (v_p1-v_p10) ; 
    diff_m_p1 = -2. * (( E_recoil / E_p1 + 1. ) * v_p1 + v_p2 + v_p3);
    
    diff_chi_p2 = 2. * err_p2.inverse(ierr) *  (v_p2-v_p20) ; 
    diff_m_p2 = -2. * (( E_recoil / E_p2 + 1. ) * v_p2 + v_p1 + v_p3 );
    
    diff_chi_p3 = 2. * err_p3.inverse(ierr) *  (v_p3-v_p30) ; 
    diff_m_p3 = -2. * (( E_recoil / E_p3 + 1. ) * v_p3 + v_p1 + v_p2);
    
    diff2_p1 = 2. * err_p1.inverse(ierr) + 2. *
      lambda * ( E_recoil / E_p1 + 1. ) * unit;
    
    diff2_p2 = 2. * err_p2.inverse(ierr) + 2. *  
      lambda * ( E_recoil / E_p2 + 1. ) * unit;

    diff2_p3 = 2. * err_p3.inverse(ierr) + 2. *  
      lambda * ( E_recoil / E_p3 + 1. ) * unit;

    diff2_p1_p2 = + 2. * lambda * unit; 
    diff2_p1_p3 = + 2. * lambda * unit; 
    diff2_p2_p3 = + 2. * lambda * unit; 

    diff_lam_p1= -1. * diff_m_p1;
    diff_lam_p2= -1. * diff_m_p2;
    diff_lam_p3= -1. * diff_m_p3;
        
    for(int im=0; im<3; ++im) {
      a[im]  = diff_chi_p1[im] - lambda * diff_m_p1[im];
      a[im+3]= diff_chi_p2[im] - lambda * diff_m_p2[im];
      a[im+6]= diff_chi_p3[im] - lambda * diff_m_p3[im];

      A[im][9]   = diff_lam_p1[im]; 
      A[im+3][9] = diff_lam_p2[im]; 
      A[im+6][9] = diff_lam_p3[im]; 
      A[9][im]   = diff_lam_p1[im];
      A[9][im+3] = diff_lam_p2[im];
      A[9][im+6] = diff_lam_p3[im];

      for(int in=0; in<3; ++in) { 
	A[im][in]     = diff2_p1[im][in];
	A[im+3][in+3] = diff2_p2[im][in];
	A[im+6][in+6] = diff2_p3[im][in];

	A[im+3][in]   = diff2_p1_p2[im][in];
	A[im][in+3]   = diff2_p1_p2[im][in];

	A[im+6][in]   = diff2_p1_p3[im][in];
	A[im][in+6]   = diff2_p1_p3[im][in];

	A[im+6][in+3]   = diff2_p2_p3[im][in];
	A[im+3][in+6]   = diff2_p2_p3[im][in];
      }
    }

    A[9][9]=0;
    a[9]=-diff;
    
    diff_a = -0.5 * A.inverse(ierr) * a;

    for(int im=0; im<3; ++im) {
      diff_p1[im] = diff_a[im];
      diff_p2[im] = diff_a[im+3];
      diff_p3[im] = diff_a[im+6];
    }

    lambda+=diff_a[9];
    v_p1+=diff_p1;
    v_p2+=diff_p2;
    v_p3+=diff_p3;

    chi2 = ((v_p1-v_p10).T() * err_p1.inverse(ierr) * (v_p1-v_p10))[0] +
           ((v_p2-v_p20).T() * err_p2.inverse(ierr) * (v_p2-v_p20))[0] +
           ((v_p3-v_p30).T() * err_p3.inverse(ierr) * (v_p3-v_p30))[0];
    
    P_p1=VectorL(v_p1[0], v_p1[1], v_p1[2], 0);
    P_p2=VectorL(v_p2[0], v_p2[1], v_p2[2], 0);
    P_p3=VectorL(v_p3[0], v_p3[1], v_p3[2], 0);
    
    P_p1.setE(sqrt(P_p1.vect().mag2()+M_p1*M_p1));
    P_p2.setE(sqrt(P_p2.vect().mag2()+M_p2*M_p2));
    P_p3.setE(sqrt(P_p3.vect().mag2()+M_p3*M_p3));
    
  }
}  

void 
fitRM2(const Particle &p1, const Particle &p2, const VectorL UPS, const double m_fit, 
     VectorL &P_p1, VectorL &P_p2, double &chi2) {

  int ierr;
  double lambda=0;
  
  P_p1=p1.p();
  P_p2=p2.p();
  
  const double M_p1=p1.pType().mass();
  const double M_p2=p2.pType().mass();
  
  HepSymMatrix err_p1(3,0), err_p2(3,0), err_p3(3,0), unit(3,0);
  for(int im=0; im<3; ++im) {
    unit[im][im]=1;
    for(int in=0; in<3; ++in) { 
      err_p1[im][in]=p1.momentum().dp()[im][in];
      err_p2[im][in]=p2.momentum().dp()[im][in];
      if (abs(err_p1[im][in])<1.e-10) err_p1[im][in]=1.e-10;
      if (abs(err_p2[im][in])<1.e-10) err_p2[im][in]=1.e-10;
    }
  }
  
  HepVector v_p1(3,0), v_p2(3,0);
  v_p1[0]=P_p1.x();
  v_p1[1]=P_p1.y();
  v_p1[2]=P_p1.z();
  
  v_p2[0]=P_p2.x();
  v_p2[1]=P_p2.y();
  v_p2[2]=P_p2.z();
  
  HepVector v_p10(3,0), v_p20(3,0);
  v_p10=v_p1;
  v_p20=v_p2;

  HepSymMatrix A(7,0);
  HepVector a(7,0), diff_a(7,0);

  HepVector diff_p1(3,0), diff_p2(3,0);
  
  HepSymMatrix diff2_p1(3,0), diff2_p2(3,0);
  HepSymMatrix diff2_p1_p2(3,0), diff2_p1_p3(3,0);
  
  double diff;
  
  HepVector diff_chi_p1(3,0), diff_chi_p2(3,0);
  HepVector diff_m_p1(3,0), diff_m_p2(3,0);
  HepVector diff_lam_p1(3,0), diff_lam_p2(3,0);
    
  for(int itt=0; itt<30; ++itt) {
    
    double E_recoil=(UPS-P_p1-P_p2).e();
    double M_recoil=(UPS-P_p1-P_p2).m();

    double E_p1=P_p1.e();
    double E_p2=P_p2.e();

    diff=M_recoil* M_recoil - m_fit* m_fit;
    if (abs(diff)<0.002) continue;

    diff_chi_p1 = 2. * err_p1.inverse(ierr) *  (v_p1-v_p10) ; 
    diff_m_p1 = -2. * (( E_recoil / E_p1 + 1. ) * v_p1 + v_p2 );
    
    diff_chi_p2 = 2. * err_p2.inverse(ierr) *  (v_p2-v_p20) ; 
    diff_m_p2 = -2. * (( E_recoil / E_p2 + 1. ) * v_p2 + v_p1 );
        
    diff2_p1 = 2. * err_p1.inverse(ierr) + 2. *
      lambda * ( E_recoil / E_p1 + 1. ) * unit;
    
    diff2_p2 = 2. * err_p2.inverse(ierr) + 2. *  
      lambda * ( E_recoil / E_p2 + 1. ) * unit;

    diff2_p1_p2 = + 2. * lambda * unit; 
    diff2_p1_p3 = + 2. * lambda * unit; 

    diff_lam_p1= -1. * diff_m_p1;
    diff_lam_p2= -1. * diff_m_p2;
        
    for(int im=0; im<3; ++im) {
      a[im]  = diff_chi_p1[im] - lambda * diff_m_p1[im];
      a[im+3]= diff_chi_p2[im] - lambda * diff_m_p2[im];

      A[im][6]   = diff_lam_p1[im]; 
      A[im+3][6] = diff_lam_p2[im]; 
      
      A[6][im]   = diff_lam_p1[im];
      A[6][im+3] = diff_lam_p2[im];
      
      for(int in=0; in<3; ++in) { 
	A[im][in]     = diff2_p1[im][in];
	A[im+3][in+3] = diff2_p2[im][in];
	
	A[im+3][in]   = diff2_p1_p2[im][in];
	A[im][in+3]   = diff2_p1_p2[im][in];
      }
    }

    A[6][6]=0;
    a[6]=-diff;
    
    diff_a = -0.5 * A.inverse(ierr) * a;

    for(int im=0; im<3; ++im) {
      diff_p1[im] = diff_a[im];
      diff_p2[im] = diff_a[im+3];
    }

    lambda+=diff_a[6];
    v_p1+=diff_p1;
    v_p2+=diff_p2;

    chi2 = ((v_p1-v_p10).T() * err_p1.inverse(ierr) * (v_p1-v_p10))[0] +
      ((v_p2-v_p20).T() * err_p2.inverse(ierr) * (v_p2-v_p20))[0];
    
    P_p1=VectorL(v_p1[0], v_p1[1], v_p1[2], 0);
    P_p2=VectorL(v_p2[0], v_p2[1], v_p2[2], 0);
    
    P_p1.setE(sqrt(P_p1.vect().mag2()+M_p1*M_p1));
    P_p2.setE(sqrt(P_p2.vect().mag2()+M_p2*M_p2));
  }
}  

void User_reco::event ( BelleEvent* evptr, int* status )
{
  *status=0;
  static int nevent=0;
  if(++nevent<2 || !(nevent%100)) cout << "Event number " << nevent << endl;

  Belle_runhead_Manager& rhdmgr = Belle_runhead_Manager::get_manager();
  Belle_runhead_Manager::const_iterator rhd = rhdmgr.begin();
  Evtcls_hadron_info_Manager&  ehimgr =
    Evtcls_hadron_info_Manager::get_manager();

  Evtcls_hadronic_flag_Manager&  ehadfl =
    Evtcls_hadronic_flag_Manager::get_manager();
  
  HepPoint3D ip_position = IpProfile::position();
  const HepSymMatrix& runIp_err = IpProfile::position_err();

  Evtcls_hadron_info_Manager::iterator iti = ehimgr.begin();
  Evtcls_hadronic_flag_Manager::iterator eti = ehadfl.begin();
    
  double r2=0;
  int ntrk=0;
  double evis=0;
  double Pz=0;
  double hjmass=0;
  
  if (iti!=ehimgr.end()){
    r2 = (*iti).R2();
    ntrk = (*iti).Ntrk();
    evis = (*iti).Evis();
    Pz   = (*iti).Pz();
    hjmass = (*iti).HeavyJetMass();
  }  

  double elec=8.0, posi=3.5;

  if(rhd!=rhdmgr.end() && rhd->EHER() > 5.){
    elec = rhd->EHER();
    posi = rhd->ELER();
  }

  VectorL UPS=
    VectorL (elec*sin(0.022), 0.,elec*cos(0.022)-posi, elec+posi);
  Particle Ups(UPS, m_ptypeUPS4);
  
  double ecm=UPS.m();

  int hadflag=-1;
  if (eti!=ehadfl.end())
    hadflag = (*eti).hadronic_flag(2);
  
  double Run = IpProfile::RunNo();
  int Exp    = IpProfile::ExpNo();
  
  Gen_hepevt_Manager &genMgr  = Gen_hepevt_Manager::get_manager();
  int ndp = 2, njp = 0, flag = 0;
  
  for(std::vector<Gen_hepevt>::iterator i=genMgr.begin();
      i!=genMgr.end(); ++i)
    if ( i -> idhep() == 443) njp++;
  
  if (njp > 0) ndp=0;
  for(std::vector<Gen_hepevt>::iterator ie=genMgr.begin(); 
      ie!=genMgr.end(); ++ie){
    
    int idd=abs(ie->idhep());
    if (idd!= 411 && idd!= 421 &&
	idd!= 431 && idd!= 4122 && idd!= 4232 ) continue;
    ndp++;
  }
  
  //  if (njp>1.5) return; 
  //  if (njp > 1.5 || ndp < 1.5 ) return;
  
  for(std::vector<Gen_hepevt>::iterator i=genMgr.begin();
      i!=genMgr.end(); ++i){
    if ( i -> idhep() != 443) continue;
    //   if ( i -> E() < 11 ) continue;
    ndp=0;
    VectorL PJ(i->PX(),i->PY(),i->PZ(),i->E());
    VectorL PLP(0,0,0,0), PLM(0,0,0,0);

    for(std::vector<Gen_hepevt>::iterator il=genMgr.begin(); 
	il!=genMgr.end(); ++il) {
      if ((il->idhep() == 11 || il->idhep() == 13) &&
	  il->mother() == *i) 
	PLM = VectorL(il->PX(), il->PY(), il->PZ(), il->E());
      if ((il->idhep() == -11 || il->idhep() == -13) &&
	  il->mother() == *i) 
	PLP = VectorL(il->PX(), il->PY(), il->PZ(), il->E());
    }      
    PJ=PLP+PLM;

    double lv1d2=cos(boostT(PLP,PJ).vect().angle(boostT(UPS,PJ).vect()));
    double cth = pStar(PJ).vect().cosTheta();
    double pj = pStar(PJ).vect().mag();

    double rand_y = double(rand())/rand_m;
    if (rand_y > (1 - .1 * cth*cth) / 1.0 ) continue; //cc
    //    if (rand_y > (1 + 5.0 * cth*cth) / 6.0 ) continue; //uu

    double rand_x = double(rand())/rand_m;
    //    if (rand_x > (1 + 0.2 * lv1d2 * lv1d2 ) / 1.2 ) continue;
    rand_x = double(rand())/rand_m;
    //    if (rand_x > (1 + 0.6 * lv1d2 * lv1d2 ) / 1.6 ) continue; //uu

    double rand_p = double(rand())/rand_m;

    // cc ****************************** 
    if (pj < 3.5 && pj > 2.5 && rand_p > 0.6 ) continue;  
    if (pj < 2.5 && pj > 1.0 && rand_p > 0.7 ) continue;  
    if (pj < 1.0 &&             rand_p > 0.26 ) continue;
    
    // uu ****************************** 
    //  if (pj > 3.0 && pj < 3.5 && rand_p > 0.7 ) continue;  
    //  if (pj > 3.5 &&             rand_p > 0.3 ) continue;  
    
    flag=1;
    
    for(std::vector<Gen_hepevt>::iterator ie=genMgr.begin(); 
	ie!=genMgr.end(); ++ie){
      if (ie == i) continue;
      
      int idd=abs(ie->idhep());
      if (idd!= 411 && idd!= 421 &&
	  idd!= 431 && idd!= 4122 && idd!= 4232 ) continue;
      
      VectorL PD(ie->PX(),ie->PY(),ie->PZ(),ie->E());
      
      t7->column("mj", njp);
      t7->column("pj", pStar(PJ).vect().mag());
      t7->column("rm", (UPS-PJ).m()); 
      
      t7->column("lv", lv1d2);
      t7->column("cth", pStar(PJ).vect().cosTheta());
      t7->column("id", idd);
      t7->column("rmx", (UPS-PJ-PD).m());
      t7->column("mx", (PJ+PD).m());
      t7->column("pd", pStar(PD).vect().mag());
      
      t7->dumpData();
    }
  
    t6->column("mj", PJ.m());
    t6->column("pj", pStar(PJ).vect().mag());
    t6->column("rm", (UPS-PJ).m()); 
    t6->column("lv", lv1d2);
    t6->column("cth", cth);
    t6->column("ndp", ndp);
    t6->dumpData();
    
  }  
  
  //  if (njp == 1 && flag == 0) return;
  //  if (njp != 1 ) cout << "J/psi " << njp << endl;
  //  if ( flag == 0) return;
  //  if ( njp == 0) return;
  //   if (ntrk < 3.5) return;

  std::vector<Particle> e_p, e_m, mu_p, mu_m;
  std::vector<Particle> en_p, en_m, mun_p, mun_m;
  
  makeLepton(e_p, e_m, mu_p, mu_m, 1);
  withDrDzCut(mu_p, 2., 4.);
  withDrDzCut(mu_m, 2., 4.);
  withDrDzCut(e_p, 2., 4.);
  withDrDzCut(e_m, 2., 4.);
  
  deepCopy(e_p, en_p);
  deepCopy(e_m, en_m);
  deepCopy(mu_p, mun_p);
  deepCopy(mu_m, mun_m);

  withLeptonIdCut(e_p, e_m, mu_p, mu_m, 0.01, 0.1);

  if (e_p.size()*e_m.size()==0 && mu_m.size()*mu_p.size() == 0) return;
  
  int lep_size=e_p.size()+e_m.size()+mu_m.size()+mu_p.size();
  int e_size=e_p.size()+e_m.size();

  std::vector<Particle> ecl, psi, psi2;
  makeBrem(ecl);
  makePsi(psi, e_p, e_m, mu_p, mu_m, ecl, 2.6, 3.9);

  withPSCut(psi, 1.9);

  if (psi.size()<0.5) return;
  *status=1;

  /*************** Make particle lists ********************************/
  
  std::vector<Particle> k_p, k_m, pi_p, pi_m;
  std::vector<Particle> k_s, pr, apr, lam, alam;
  std::vector<Particle> pin_p, pin_m;
  std::vector<Particle> gam, gam_en, all, pi_0;
  
  makeKPi(k_p, k_m, pi_p, pi_m, 1);

  withDrDzCut(pi_p, 2., 4.);
  withDrDzCut(pi_m, 2., 4.);
  withDrDzCut(k_p,  2., 4.);
  withDrDzCut(k_m,  2., 4.);

  deepCopy(pi_p, pin_p);
  deepCopy(pi_m, pin_m);
  ntrk=pin_p.size()+pin_m.size();
  
  withKaonIdCut(k_p, k_m, 0.6);
  
  makeProton(pr, apr, 1);
  withDrDzCut(pr, 2., 4.);
  withDrDzCut(apr, 2., 4.);
  withProtonIdCut(pr, apr, 0.6);

  makeKs (k_s);
  for(std::vector<Particle>::iterator l = k_s.begin(); l!=k_s.end(); ++l) {
    HepPoint3D V(l->mdstVee2().vx(),l->mdstVee2().vy(),0);
    Vector3 P(l->px(),l->py(),0);
    V=V-ip_position;
    V.setZ(0.);
    if (abs(l->mass()-0.4977)>0.015 || V.perp()<0.1 || 
	cos(V.angle(P))<0.99 || l->mdstVee2().z_dist()>1. ) {
      k_s.erase(l); --l;
    }
  }

  makeLambda (lam,alam);

  for(std::vector<Particle>::iterator l = lam.begin(); l!=lam.end(); ++l) {
    HepPoint3D V(l->mdstVee2().vx(),l->mdstVee2().vy(),0);
    Vector3 P(l->px(),l->py(),0);
    V=V-ip_position;
    V.setZ(0.);
    if (abs(l->mass()-1.1157)>0.01 || V.perp()<0.1 || 
	cos(V.angle(P))<0.99 || l->mdstVee2().z_dist()>1. ) {
      lam.erase(l); --l;
    }
  }

  for(std::vector<Particle>::iterator l = alam.begin(); l!=alam.end(); ++l) {
    HepPoint3D V(l->mdstVee2().vx(),l->mdstVee2().vy(),0);
    Vector3 P(l->px(),l->py(),0);
    V=V-ip_position;
    V.setZ(0.);
    if (abs(l->mass()-1.1157)>0.01 || V.perp()<0.1 || 
	cos(V.angle(P))<0.99 || l->mdstVee2().z_dist()>1. ) {
      alam.erase(l); --l;
    }
  }

  makeGamma (gam);
  withEminCut(gam, 0.05);
  deepCopy (gam, gam_en);
  withPSCut(gam_en, 0.5);

  VectorL Gam, Pi0, Eta;
  double Emax_gam=0;
  for(std::vector<Particle>::iterator r = gam.begin(); 
      r!=gam.end(); ++r) {
    if (pStar(r->p(),elec,posi).e()>Emax_gam) {
      Emax_gam=pStar(r->p(),elec,posi).e();
      Gam=r->p();
    }
  }

  //  deepCopy(k_p, all);
  //  deepCopy(k_m, all);
  deepCopy(pi_p, all);
  deepCopy(pi_m, all);
  //  deepCopy(gam,  all);

  makePi0 (pi_0);
  withEminCutPi0(pi_0, 0.05);

  for(std::vector<Particle>::iterator l = pi_0.begin(); l!=pi_0.end(); ++l) 
    if (abs(l->mdstPi0().mass()-0.134)>0.015)
      {pi_0.erase(l); --l;}
  
  for(std::vector<Particle>::iterator l = pi_0.begin(); l!=pi_0.end(); ++l) 
    if (pStar(l->p(),elec,posi).e()>pStar(Pi0).e())
      Pi0=l->p();
  
  for(std::vector<Particle>::iterator l=gam.begin(); l!=gam.end(); ++l){
    for(std::vector<Particle>::iterator l1=gam.begin(); l1!=gam.end(); ++l1){ 
      if (abs((l->p()+l1->p()).m()-0.548 )<0.02 &&
	  pStar((l->p()+l1->p()),elec,posi).e()>pStar(Eta).e())
	Eta=(l->p()+l1->p());
    }
  }
  
  for(int p=0; p<pi_0.size(); ++p){
    setPi0Error(pi_0);
    doKmvFit(pi_0[p]);
  }

  for(int p=0; p<k_s.size(); ++p)
    doKmvFit(k_s[p]);

  for(int p=0; p<lam.size(); ++p)
    doKmvFit(lam[p]);

  for(int p=0; p<alam.size(); ++p)
    doKmvFit(alam[p]);

  setGenHepInfoF(e_p); 
  setGenHepInfoF(mu_p); 
  setGenHepInfoF(e_m);
  setGenHepInfoF(mu_m);

  setGenHepInfoF(pr);
  setGenHepInfoF(k_p);
  setGenHepInfoF(pi_p);
  setGenHepInfoF(apr);
  setGenHepInfoF(k_m);
  setGenHepInfoF(pi_m);
  setGenHepInfoP(pi_0);
  setGenHepInfoKs(k_s);

  std::vector<Particle> DA, D0, DB, Dp, Dm, Dsp, Dsm, Lamc, Alamc, DD;

  combination(D0,  m_ptypeD0 , k_m, pi_p, 0.1);
  setUserInfo(D0,  1);
  combination(DB,  m_ptypeD0B , k_p, pi_m, 0.1);
  setUserInfo(DB,  1);

  combination(DA,  m_ptypeD0 , k_p, k_m, 0.1);
  setUserInfo(DA,  2);

  combination(D0,  m_ptypeD0 , pi_p, pi_p, k_m, pi_m, 0.1);
  setUserInfo(D0,  3);
  combination(DB,  m_ptypeD0B , pi_m, pi_m, k_p, pi_p, 0.1);
  setUserInfo(DB,  3);

  combination(D0,  m_ptypeD0 , pi_0, k_m, pi_p, 0.1);
  setUserInfo(D0,  4);
  combination(DB,  m_ptypeD0B , pi_0, k_p, pi_m, 0.1);
  setUserInfo(DB,  4);

  combination(DA,  m_ptypeD0 , k_s, pi_m, pi_p, 0.1);
  setUserInfo(DA,  5);

  combination(DA,  m_ptypeD0 , k_s, k_m, k_p, 0.1);
  setUserInfo(DA,  6);

  combination(DA,  m_ptypeD0 , k_s,  pi_0, 0.1);
  setUserInfo(DA,  7);

  combination(Dm,  m_ptypeDM , pi_m, pi_m, k_p, 0.1);
  setUserInfo(Dm,  8);
  combination(Dp,  m_ptypeDP , pi_p, pi_p, k_m, 0.1);
  setUserInfo(Dp,  8);

  combination(Dm,  m_ptypeDM , pi_m, k_m, k_p, 0.1);
  setUserInfo(Dm,  9);
  combination(Dp,  m_ptypeDP , pi_p, k_p, k_m, 0.1);
  setUserInfo(Dp,  9);

  combination(Dm,  m_ptypeDM , k_s, pi_m, 0.1);
  setUserInfo(Dm,  10);
  combination(Dp,  m_ptypeDP , k_s, pi_p, 0.1);
  setUserInfo(Dp,  10);

  combination(Dm,  m_ptypeDM , pi_m, pi_m, k_s, pi_p, 0.1);
  setUserInfo(Dm,  11);
  combination(Dp,  m_ptypeDP , pi_p, pi_p, k_s, pi_m, 0.1);
  setUserInfo(Dp,  11);
  
  combination(Dm,  m_ptypeDM , k_s, k_m, 0.1);
  setUserInfo(Dm,  12);
  combination(Dp,  m_ptypeDP , k_s, k_p, 0.1);
  setUserInfo(Dp,  12);

  combination(Dsm,  m_ptypeDSM , pi_m, k_m, k_p, 0.1);
  setUserInfo(Dsm,  13);
  combination(Dsp,  m_ptypeDSP , pi_p, k_p, k_m, 0.1);
  setUserInfo(Dsp,  13);

  combination(Dsm,  m_ptypeDSM , k_m, k_s, 0.1);
  setUserInfo(Dsm,  14);
  combination(Dsp,  m_ptypeDSP , k_p, k_s, 0.1);
  setUserInfo(Dsp,  14);
	      
  combination(Lamc,  m_ptypeLAMC , pi_p, k_m, pr, 0.1);
  setUserInfo(Lamc,  15);
  combination(Alamc,  m_ptypeALAMC , pi_m, k_p, apr, 0.1);
  setUserInfo(Alamc,  15);

  combination(Lamc,  m_ptypeLAMC , k_s, pr, 0.1);
  setUserInfo(Lamc,  16);
  combination(Alamc,  m_ptypeALAMC , k_s, apr, 0.1);
  setUserInfo(Alamc,  16);

  combination(Lamc,  m_ptypeLAMC , lam, pi_p, 0.1);
  setUserInfo(Lamc,  17);
  combination(Alamc,  m_ptypeALAMC , alam, pi_m, 0.1);
  setUserInfo(Alamc,  17);

  float chisq;
  doKvFit(D0);
  doKvFit(DB);
  doKvFit(DA);
  doKvFit(Dp);
  doKvFit(Dm);
  doKvFit(Dsp);
  doKvFit(Dsm);
  doKvFit(Lamc);
  doKvFit(Alamc);

  deepCopy(DA,  DD);
  deepCopy(D0,  DD);
  deepCopy(DB,  DD);
  deepCopy(Dp,  DD);
  deepCopy(Dm,  DD);
  deepCopy(Dsp,  DD);
  deepCopy(Dsm,  DD);
  deepCopy(Lamc,  DD);
  deepCopy(Alamc,  DD);

  setGenHepInfoT(D0);
  setGenHepInfoT(DB);
  setGenHepInfoT(DA);
  setGenHepInfoT(DD);
  setGenHepInfoT(psi);

  for(int j=0; j<DD.size(); ++j) {
    Particle &D=DD[j];
    int chD = dynamic_cast<UserInfo&>(D.userInfo()).channel();
    float MD =dynamic_cast<UserInfo&>(D.userInfo()).vmass();
    float MD_nom=D.pType().mass();
    float dmassd = (MD-MD_nom)/2./dwid[chD-1];
    int bind = int(floor(dmassd+0.5));
    float mass_fit=MD_nom+bind*2.*dwid[chD-1];
    reFitDSb (D, mass_fit);
  }  
  
  doKmvFit(DA,chisq);
  doKmvFit(D0,chisq);
  doKmvFit(DB,chisq);
  doKmvFit(Dp,chisq);
  doKmvFit(Dm,chisq); 

  std::vector<Particle> Dstp, Dstm, Dst0, DstB, Dst, D2;
  combination(Dst,  m_ptypeDstarP , D0, pi_p, 0.025);
  setUserInfo(Dst,  1);
  combination(Dst,  m_ptypeDstarM , DB, pi_m, 0.025);
  setUserInfo(Dst,  1);
  combination(Dst,  m_ptypeDstarP , DA, pi_p, 0.025);
  setUserInfo(Dst,  1);
  combination(Dst,  m_ptypeDstarM , DA, pi_m, 0.025);
  setUserInfo(Dst,  1);
  combination(Dst,  m_ptypeDstarP , Dp, pi_0, 0.025);
  setUserInfo(Dst,  2);
  combination(Dst,  m_ptypeDstarM , Dm, pi_0, 0.025);
  setUserInfo(Dst,  2);
  combination(Dst,  m_ptypeDstar0 , DA, pi_0, 0.025);
  setUserInfo(Dst,  3);
  combination(Dst,  m_ptypeDstar0 , D0, pi_0, 0.025);
  setUserInfo(Dst,  3);
  combination(Dst,  m_ptypeDstarB , DB, pi_0, 0.025);
  setUserInfo(Dst,  3);

  combination(D2,  m_ptypeDstarB , DA, pi_m, 0.8);
  setUserInfo(D2,  1);
  combination(D2,  m_ptypeDstarB , DA, pi_p, 0.8);
  setUserInfo(D2,  1);
  combination(D2,  m_ptypeDstarB , DB, pi_m, 0.8);
  setUserInfo(D2,  1);
  combination(D2,  m_ptypeDstarB , D0, pi_p, 0.8);
  setUserInfo(D2,  1);
  combination(D2,  m_ptypeDstarB , Dp, pi_m, 0.8);
  setUserInfo(D2,  2);
  combination(D2,  m_ptypeDstarB , Dm, pi_p, 0.8);
  setUserInfo(D2,  2);

  setGenHepInfoT(Dst);
  doKvFit(Dst);
  doKmvFit(Dst,chisq);

  reFitPsiSb(psi, ecl); 

  combination (psi2, m_ptypePSI2, psi, pi_p, pi_m, 1.);
  int psi2_flag=0;
  
  for(std::vector<Particle>::iterator l = psi2.begin(); l!=psi2.end(); ++l) 
    if (abs(l->mass()-l->child(0).mass()+3.097-3.686)<0.012)
      psi2_flag++;
  
  doKvFit(psi2);
  doKmvFit(psi2,chisq);
  
  std::vector<Particle> B;
  combination (B, m_ptypeB0, psi, DD);
  
  std::vector<Particle> B2;
  combination (B2, m_ptypeB0, psi, D2);

  //  std::vector<Particle> B2;
  //  combination (B2, m_ptypeB0, psi2, DD);
  
  std::vector<Particle> Bds;
  combination (Bds, m_ptypeB0, psi, DA, DA);
  combination (Bds, m_ptypeB0, psi, DA, D0);
  combination (Bds, m_ptypeB0, psi, DA, DB);
  combination (Bds, m_ptypeB0, psi, DB, D0);
  combination (Bds, m_ptypeB0, psi, Dp, Dm);
  combination (Bds, m_ptypeB0, psi, Dsp, Dsm);

  std::vector<Particle> Bds2;
  combination (Bds2, m_ptypeB0, psi2, Dst);

  VectorL All(0,0,0,0);
  int nk=k_p.size()+k_m.size()+k_s.size();

  for(std::vector<Particle>::iterator l = all.begin(); l!=all.end(); ++l)
    All+=pStar(*l,elec,posi);
  
  for(int j=0; j<psi.size(); ++j) {
    Particle p=psi[j];
    float Mj =dynamic_cast<UserInfo&>(p.userInfo()).vmass();
    VectorL PP=p.p();
    double pj=pStar(PP,elec,posi).vect().mag();
    float dwd=0.03;
    float dmass = (Mj-3.097)/2./dwd;
    int bin = int(floor(dmass+0.5));

    Particle c0=p.child(0);
    Particle c1=p.child(1);
    Vector3 V = p.momentum().decayVertex()-ip_position;   
    
    double mulh1=Muid_mdst(c0.mdstCharged()).Muon_likelihood();
    double mulh2=Muid_mdst(c1.mdstCharged()).Muon_likelihood();
    double elh1= eid(c0.mdstCharged()).prob(3, -1, 5);
    double elh2= eid(c1.mdstCharged()).prob(3, -1, 5);

    VectorL P_apsi=pStar(UPS-PP,elec,posi);
    
    std::vector<Hep3Vector> vec;
    VectorL Rest;
    for(std::vector<Particle>::iterator r = all.begin(); r!=all.end(); ++r) {
      if (checkSame(*r,p)) continue;
      VectorL Ptmp=r->p();
      Rest+=Ptmp;
      HepLorentzVector b0(pStar(Ptmp,elec,posi));
      b0.boost(-(P_apsi.boostVector()));
      vec.push_back(b0.vect());
    }

    Vector3 thr= thrust(vec.begin(),vec.end(),SelfFunc(Hep3Vector()));
    double athr=cos(thr.angle(P_apsi.vect()));

    double rm=(UPS-PP).m();
    double cth=pStar(p,elec,posi).vect().cosTheta();
    double lv1d2=cos(boostT(c1,p).vect().angle(boostT(Ups,p).vect()));

    int id =0;
    if (p.relation().genHepevt())
      id=p.relation().genHepevt().idhep();
        
    if (c0.relation().genHepevt() && c1.relation().genHepevt()) {
      if (c0.relation().genHepevt().mother() &&
	  c1.relation().genHepevt().mother() ) {
	if (c0.relation().genHepevt().mother().E()>11 &&
	    c1.relation().genHepevt().mother().E()>11 ) {
	  id=443;  
	}
      }
    }
 
    t1->column("mj", Mj);
    t1->column("bin", bin);
    t1->column("pj", pj);
    t1->column("rm", rm); 
    t1->column("ecm", ecm); 

    t1->column("vp", V.perp());
    t1->column("lhm", min(mulh1, mulh2));
    t1->column("lhe", min(elh1, elh2));

    t1->column("lv", lv1d2);
    t1->column("cth", cth);
    t1->column("thr", athr);
    
    t1->column("r2", r2);
    t1->column("ea", All.e());
    t1->column("nt", ntrk);
    t1->column("ne", e_size);
    t1->column("nl", lep_size);
    t1->column("nk", nk);
    t1->column("em", Emax_gam);
    t1->column("f2", psi2_flag);
    t1->column("id", id);
    t1->column("run", 10000*Exp+Run);

    t1->dumpData();
  }


  for(int j=0; j<psi2.size(); ++j) {
    Particle p2=psi2[j];
    Particle p=p2.child(0);
    float Mj =dynamic_cast<UserInfo&>(p.userInfo()).vmass();
    float Mp2 =dynamic_cast<UserInfo&>(p2.userInfo()).vmass();
    if (abs(Mp2-3.686)>0.1) continue;
    VectorL PP=p.p();
    double pj=pStar(PP,elec,posi).vect().mag();
    float dwd=0.03;
    float dmass = (Mj-3.097)/2./dwd;
    int bin = int(floor(dmass+0.5));

    Particle c0=p.child(0);
    Particle c1=p.child(1);
    Vector3 V = p.momentum().decayVertex()-ip_position;   

    double mulh1=Muid_mdst(c0.mdstCharged()).Muon_likelihood();
    double mulh2=Muid_mdst(c1.mdstCharged()).Muon_likelihood();
    double elh1= eid(c0.mdstCharged()).prob(3, -1, 5);
    double elh2= eid(c1.mdstCharged()).prob(3, -1, 5);

    VectorL P_apsi=pStar(UPS-PP,elec,posi);
    
    std::vector<Hep3Vector> vec;
    VectorL Rest;
    for(std::vector<Particle>::iterator r = all.begin(); r!=all.end(); ++r) {
      if (checkSame(*r,p)) continue;
      VectorL Ptmp=r->p();
      Rest+=Ptmp;
      HepLorentzVector b0(pStar(Ptmp,elec,posi));
      b0.boost(-(P_apsi.boostVector()));
      vec.push_back(b0.vect());
    }

    Vector3 thr= thrust(vec.begin(),vec.end(),SelfFunc(Hep3Vector()));
    double athr=cos(thr.angle(P_apsi.vect()));

    double rm=(UPS-PP).m();
    double rm2=(UPS-p2.p()).m();
    double cth=pStar(p,elec,posi).vect().cosTheta();
    double lv1d2=cos(boostT(c1,p).vect().angle(boostT(Ups,p).vect()));

    int id =0;
    if (p.relation().genHepevt())
      id=p.relation().genHepevt().idhep();
   
    t8->column("mp2", Mp2);
    t8->column("mj", Mj);
    t8->column("bin", bin);
    t8->column("pj", pj);
    t8->column("p2", pStar(p2).vect().mag());
    t8->column("rm", rm); 
    t8->column("rm2", rm2); 
    t8->column("ecm", ecm); 

    t8->column("vp", V.perp());
    t8->column("lhm", min(mulh1, mulh2));
    t8->column("lhe", min(elh1, elh2));

    t8->column("lv", lv1d2);
    t8->column("cth", cth);
    t8->column("thr", athr);
    
    t8->column("r2", r2);
    t8->column("ea", All.e());
    t8->column("nt", ntrk);
    t8->column("ne", e_size);
    t8->column("nl", lep_size);
    t8->column("nk", nk);
    t8->column("em", Emax_gam);
    t8->column("f2", psi2_flag);
    t8->column("id", id);
    t8->column("run", 10000*Exp+Run);

    t8->dumpData();
  }

  init_best();
   
  for(int id=0; id<B.size(); ++id) {
    Particle p=B[id].child(0);
    Particle D=B[id].child(1);
    float Mj =dynamic_cast<UserInfo&>(p.userInfo()).vmass();
    Vector3 V = p.momentum().decayVertex()-ip_position;
    if (abs(Mj-3.097)>0.03 || V.perp()>0.1 ) continue;

    int chD = dynamic_cast<UserInfo&>(D.userInfo()).channel();
    float MD =dynamic_cast<UserInfo&>(D.userInfo()).vmass();
    float MD_nom =D.pType().mass();
    float dmass = (MD-MD_nom)/2./dwid[chD-1];
    int ibin = int(dmass+4.5);
    Vector3 VD = D.momentum().decayVertex()-ip_position;
    if (ibin<-0.5 || ibin>8.5 || VD.perp()>0.1 ) continue;
    
    double rm_p_D_fit=(UPS-B[id].p()).m();
    double Mdst=2.007, Md=1.8645;
    if (abs(D.pType().lund())==411) { Mdst=2.01; Md=1.869;}
    if (abs(D.pType().lund())==431) { Mdst=2.11; Md=1.969;}
    if (abs(D.pType().lund())==4122){ Mdst=2.285; Md=2.285;}

    if (rm_p_D_fit<0) continue; 
    if (abs(rm_p_D_fit-Mdst)<0.07 && 
	abs(dmass+4-ibin)<abs(dmass_best[ibin])){ 
      dmass_best[ibin]=dmass+4-ibin;
      idbest[ibin]=id;
    }
    if (abs(rm_p_D_fit-Md)<0.07 && 
	abs(dmass+4-ibin)<abs(dmass_best1[ibin])){ 
      dmass_best1[ibin]=dmass+4-ibin;
      idbest1[ibin]=id;
    }
    if (rm_p_D_fit<3. && rm_p_D_fit>1. &&
	abs(dmass+4-ibin)<abs(dmass_best2[ibin])){ 
      dmass_best2[ibin]=dmass+4-ibin;
      idbest2[ibin]=id;
    }
  }
  
  for(int id=0; id<B.size(); ++id) {
    Particle p=B[id].child(0);
    float Mj =dynamic_cast<UserInfo&>(p.userInfo()).vmass();
    VectorL PP=p.p();
    double pj=pStar(PP,elec,posi).vect().mag();
    float dwd=0.03;
    float dmass = (Mj-3.097)/2./dwd;
    int bin = int(floor(dmass+0.5));

    Particle c0=p.child(0);
    Particle c1=p.child(1);
    Vector3 V = p.momentum().decayVertex()-ip_position;   

    double mulh1=Muid_mdst(c0.mdstCharged()).Muon_likelihood();
    double mulh2=Muid_mdst(c1.mdstCharged()).Muon_likelihood();
    double elh1= eid(c0.mdstCharged()).prob(3, -1, 5);
    double elh2= eid(c1.mdstCharged()).prob(3, -1, 5);

    double cth=pStar(p,elec,posi).vect().cosTheta();
    double lv1d2=cos(boostT(c1,p).vect().angle(boostT(Ups,p).vect()));

    Particle D=B[id].child(1);

    int chD = dynamic_cast<UserInfo&>(D.userInfo()).channel();
    float MD =dynamic_cast<UserInfo&>(D.userInfo()).vmass();
    float MD_nom=D.pType().mass();
    float dmassd = (MD-MD_nom)/2./dwid[chD-1];
    int bind = int(floor(dmassd+0.5));
    int ibin = int(dmassd+4.5);
    if (abs(bind)>4.5) continue;

    VectorL PD=D.p();
    VectorL Rec=UPS-p.p();
    
    double rm=Rec.m();
    double rmx=(UPS-B[id].p()).m();
    //    if (rmx > 3. || rm > 6. ) continue;
    //    if (chD != 1 && chD != 3 && chD != 8 && chD !=10 &&
    //	chD != 13 && chD != 14 && chD != 15 ) continue;
    
    double mks=0, mph=0;
    if (chD == 13) {
      mph=(D.child(1).p()+D.child(2).p()).m();
      mks=(D.child(0).p()+D.child(2).p()).m();
    }
    
    Vector3 VD = D.momentum().decayVertex()-ip_position;   
    
    double min_pid=1;
    for(unsigned i=0; i<D.nChildren(); i++)
      if (abs(D.child(i).pType().lund())==211)
	if (min_pid>
	    atc_pid(3,1,5,2,3).prob(&(D.child(i).mdstCharged())))
	  min_pid=
	    atc_pid(3,1,5,2,3).prob(&(D.child(i).mdstCharged()));

    double lvd0=cos(boostT(PD,Rec).vect().angle(boostT(UPS,Rec)));
    double pdr=boostT(PD,Rec).vect().mag();
    
    t2->column("mj", Mj);
    t2->column("bin", bin);
    t2->column("pj", pj);
    t2->column("rm", rm); 
    t2->column("ecm", ecm); 

    t2->column("vp", V.perp());
    t2->column("lhm", min(mulh1, mulh2));
    t2->column("lhe", min(elh1, elh2));

    t2->column("lv", lv1d2);
    t2->column("cth", cth);
    
    t2->column("r2", r2);
    t2->column("ea", All.e());
    t2->column("nt", ntrk);
    t2->column("ne", e_size);
    t2->column("nl", lep_size);
    t2->column("nk", nk);
    t2->column("em", Emax_gam);
    t2->column("f2", psi2_flag);

    t2->column("md", MD);
    t2->column("vd", VD.perp());
    t2->column("dvz", (VD-V).z());
    t2->column("mpi", min_pid);
    t2->column("pd", pStar(D).vect().mag());
    t2->column("lvd0", lvd0);
    t2->column("chd", chD);
    t2->column("mph", mph);
    t2->column("mks", mks);
    t2->column("bin1", bind); 
    t2->column("b1", id-idbest[ibin]); 
    t2->column("b2", id-idbest1[ibin]); 
    t2->column("b3", id-idbest2[ibin]); 

    t2->column("rmx", rmx);
    t2->column("mx", (p.p()+D.p()).m());

    VectorL P_new_psi, P_new_D;
    double chi2_fit, chi2_fitD; 

    double m_fit=2.01;
    if (abs(D.pType().lund())== 421) m_fit=2.007;
    if (abs(D.pType().lund())== 431) m_fit=2.112;
    if (abs(D.pType().lund())== 4122) m_fit=2.285;
  
    fitRM2(p, D, UPS, m_fit, P_new_psi, P_new_D, chi2_fit);
    t2->column("mftds", m_fit);
    
    m_fit=1.869;
    if (abs(D.pType().lund())== 421) m_fit=1.8645;
    if (abs(D.pType().lund())== 431) m_fit=1.969;
    if (abs(D.pType().lund())== 4122) m_fit=2.285;
    t2->column("mftd", m_fit);
    
    VectorL P_new_psiD, P_new_DD;
    fitRM2(p, D, UPS, m_fit, P_new_psiD, P_new_DD, chi2_fit);

    Particle RecDs(UPS-P_new_psi, m_ptypeUPS4);
    double lv_ds=cos(boostT(D,RecDs).vect().angle(boostT(Ups,RecDs).vect()));

    Particle RecD(UPS-P_new_psiD, m_ptypeUPS4);
    double lv_d=cos(boostT(D,RecD).vect().angle(boostT(Ups,RecD).vect()));

    t2->column("rmxds", (UPS-P_new_psi-P_new_D).m());
    t2->column("rmfds", (UPS-P_new_psi).m());
    t2->column("lvds", lv_ds);

    t2->column("rmxd", (UPS-P_new_psiD-P_new_DD).m());
    t2->column("rmfd", (UPS-P_new_psiD).m());
    t2->column("lvd", lv_d);

    double bestrmd=10;
    VectorL Par=UPS-P_new_psi-P_new_D;

    if (D.lund()==411) {
      for(int ips=0; ips<pi_m.size(); ++ips) {
	double rmd= Par.m()-(Par-pi_m[ips].p()).m()-0.1454;
	if (bestrmd>rmd) bestrmd=rmd;
      }
      
      for(int ips=0; ips<pi_0.size(); ++ips) {
	double rmd= Par.m()-(Par-pi_0[ips].p()).m()-0.1406;
	if (bestrmd>rmd) bestrmd=rmd;
      }
    }


    if (D.lund()==-411) {
      for(int ips=0; ips<pi_p.size(); ++ips) {
	double rmd= Par.m()-(Par-pi_p[ips].p()).m()-0.1454;
	if (bestrmd>rmd) bestrmd=rmd;
      }
      
      for(int ips=0; ips<pi_0.size(); ++ips) {
	double rmd= Par.m()-(Par-pi_0[ips].p()).m()-0.1406;
	if (bestrmd>rmd) bestrmd=rmd;
      }
    }

    if (abs(D.lund())==421) {
      for(int ips=0; ips<pi_0.size(); ++ips) {
	double rmd= Par.m()-(Par-pi_0[ips].p()).m()-0.1421;
	if (bestrmd>rmd) bestrmd=rmd;
      }
    }

    t2->column("rmd", bestrmd);

    int hepid =0, hepid_mo=0, id=0;
    if (D.relation().genHepevt()){
      hepid=D.relation().genHepevt().idhep();
      if (p.relation().genHepevt()){
	Particle p_mc(p.relation().genHepevt());
	Particle D_mc(D.relation().genHepevt());
      } 
      if (D.relation().genHepevt().mother()) 
	hepid_mo=D.relation().genHepevt().mother().idhep();
    }
    if (p.relation().genHepevt())
      id=p.relation().genHepevt().idhep();
    
    t2->column("hepid", hepid); 
    t2->column("hepidm", hepid_mo); 
    t2->column("id", id);

    t2->dumpData();
  }

  init_best();
  
  for(int id=0; id<B2.size(); ++id) {
    Particle p=B2[id].child(0);
    Particle D=B2[id].child(1).child(0);
    Particle pi=B2[id].child(1).child(1);
    float Mj =dynamic_cast<UserInfo&>(p.userInfo()).vmass();
    Vector3 V = p.momentum().decayVertex()-ip_position;
    if (abs(Mj-3.097)>0.03 || V.perp()>0.1 ) continue;

    int chD = dynamic_cast<UserInfo&>(D.userInfo()).channel();
    float MD =dynamic_cast<UserInfo&>(D.userInfo()).vmass();
    float MD_nom =D.pType().mass();
    float dmass = (MD-MD_nom)/2./dwid[chD-1];
    int ibin = int(dmass+4.5);
    Vector3 VD = D.momentum().decayVertex()-ip_position;
    if (ibin<-0.5 || ibin>8.5 || VD.perp()>0.1 ) continue;
    
    double rm_p_D_fit=(UPS-B2[id].p()).m();
    double Mdst=2.007, Md=1.8645;

    if (rm_p_D_fit<0) continue; 
    if (abs(rm_p_D_fit-Mdst)<0.07 && 
	abs(dmass+4-ibin)<abs(dmass_best[ibin])){ 
      dmass_best[ibin]=dmass+4-ibin;
      idbest[ibin]=id;
    }
    if (abs(rm_p_D_fit-Md)<0.07 && 
	abs(dmass+4-ibin)<abs(dmass_best1[ibin])){ 
      dmass_best1[ibin]=dmass+4-ibin;
      idbest1[ibin]=id;
    }
    if (rm_p_D_fit<3. && rm_p_D_fit>1. &&
	abs(dmass+4-ibin)<abs(dmass_best2[ibin])){ 
      dmass_best2[ibin]=dmass+4-ibin;
      idbest2[ibin]=id;
    }
  }

  
  for(int id=0; id<B2.size(); ++id) {
    Particle p=B2[id].child(0);
    Particle D=B2[id].child(1).child(0);
    Particle pi=B2[id].child(1).child(1);
    float Mj =dynamic_cast<UserInfo&>(p.userInfo()).vmass();
    float dwd=0.03;
    float dmass = (Mj-3.097)/2./dwd;
    int bin = int(floor(dmass+0.5));
    Vector3 V = p.momentum().decayVertex()-ip_position;

    int chD = dynamic_cast<UserInfo&>(D.userInfo()).channel();
    float MD2 = (D.p()+pi.p()).m();
    float MD =dynamic_cast<UserInfo&>(D.userInfo()).vmass();
    float MD_nom =D.pType().mass();
    float dmassd = (MD-MD_nom)/2./dwid[chD-1];
    Vector3 VD = D.momentum().decayVertex()-ip_position;

    VectorL PP=p.p();
    double pj=pStar(PP,elec,posi).vect().mag();

    Particle c0=p.child(0);
    Particle c1=p.child(1);

    double mulh1=Muid_mdst(c0.mdstCharged()).Muon_likelihood();
    double mulh2=Muid_mdst(c1.mdstCharged()).Muon_likelihood();
    double elh1= eid(c0.mdstCharged()).prob(3, -1, 5);
    double elh2= eid(c1.mdstCharged()).prob(3, -1, 5);

    double cth=pStar(p,elec,posi).vect().cosTheta();
    double lv1d2=cos(boostT(c1,p).vect().angle(boostT(Ups,p).vect()));

    int bind = int(floor(dmassd+0.5));
    int ibin = int(dmassd+4.5);
    if (abs(bind)>4.5) continue;

    VectorL PD=D.p();
    VectorL Rec=UPS-p.p();
    
    double rm=Rec.m();
    double rmx=(UPS-B2[id].p()).m();

    double mks=0, mph=0;
    if (chD == 13) {
      mph=(D.child(1).p()+D.child(2).p()).m();
      mks=(D.child(0).p()+D.child(2).p()).m();
    }

    double min_pid=1;
    for(unsigned i=0; i<D.nChildren(); i++)
      if (abs(D.child(i).pType().lund())==211)
	if (min_pid>
	    atc_pid(3,1,5,2,3).prob(&(D.child(i).mdstCharged())))
	  min_pid=
	    atc_pid(3,1,5,2,3).prob(&(D.child(i).mdstCharged()));

    double lvd0=cos(boostT(PD,Rec).vect().angle(boostT(UPS,Rec)));
    double pdr=boostT(PD,Rec).vect().mag();

    t4->column("mj", Mj);
    t4->column("bin", bin);
    t4->column("pj", pj);
    t4->column("rm", rm); 
    t4->column("ecm", ecm); 

    t4->column("lhm", min(mulh1, mulh2));
    t4->column("lhe", min(elh1, elh2));

    t4->column("lv", lv1d2);
    t4->column("cth", cth);
    
    t4->column("r2", r2);
    t4->column("nt", ntrk);
    t4->column("ne", e_size);
    t4->column("nl", lep_size);
    t4->column("nk", nk);
    t4->column("em", Emax_gam);
    t4->column("f2", psi2_flag);

    t4->column("md", MD);
    t4->column("md2", MD2);
    t4->column("pd", pStar(D).vect().mag());
    t4->column("lvd0", lvd0);
    t4->column("chd", chD);
    t4->column("bin1", bind); 
    t4->column("b1", id-idbest[ibin]); 
    t4->column("b2", id-idbest1[ibin]); 
    t4->column("b3", id-idbest2[ibin]); 

    t4->column("rmx", rmx);
    t4->column("mx", (p.p()+D.p()).m());

    VectorL P_new_psi, P_new_D, P_new_pi;
    double chi2_fit, chi2_fitD; 

    double m_fit=2.01;
    fitRM3(p, D, pi, UPS, m_fit, P_new_psi, P_new_D, P_new_pi, chi2_fit);
    t4->column("mftds", m_fit);

    m_fit=1.869;
    t4->column("mftd", m_fit);

    VectorL P_new_psiD, P_new_DD, P_new_piD;
    fitRM3(p, D, pi, UPS, m_fit, P_new_psiD, P_new_DD, P_new_piD, chi2_fit);

    Particle RecDs(UPS-P_new_psi, m_ptypeUPS4);
    double lv_ds=cos(boostT(D,RecDs).vect().angle(boostT(Ups,RecDs).vect()));

    Particle RecD(UPS-P_new_psiD, m_ptypeUPS4);
    double lv_d=cos(boostT(D,RecD).vect().angle(boostT(Ups,RecD).vect()));

    t4->column("rmxds", (UPS-P_new_psi-P_new_D).m());
    t4->column("rmfds", (UPS-P_new_psi).m());
    t4->column("lvds", lv_ds);

    t4->column("rmxd", (UPS-P_new_psiD-P_new_DD).m());
    t4->column("rmfd", (UPS-P_new_psiD).m());
    t4->column("lvd", lv_d);

    double bestrmd=10;
    VectorL Par=UPS-P_new_psi-P_new_D;

    if (D.lund()==411) {
      for(int ips=0; ips<pi_m.size(); ++ips) {
	double rmd= Par.m()-(Par-pi_m[ips].p()).m()-0.1454;
	if (bestrmd>rmd) bestrmd=rmd;
      }
      
      for(int ips=0; ips<pi_0.size(); ++ips) {
	double rmd= Par.m()-(Par-pi_0[ips].p()).m()-0.1406;
	if (bestrmd>rmd) bestrmd=rmd;
      }
    }


    if (D.lund()==-411) {
      for(int ips=0; ips<pi_p.size(); ++ips) {
	double rmd= Par.m()-(Par-pi_p[ips].p()).m()-0.1454;
	if (bestrmd>rmd) bestrmd=rmd;
      }
      
      for(int ips=0; ips<pi_0.size(); ++ips) {
	double rmd= Par.m()-(Par-pi_0[ips].p()).m()-0.1406;
	if (bestrmd>rmd) bestrmd=rmd;
      }
    }

    if (abs(D.lund())==421) {
      for(int ips=0; ips<pi_0.size(); ++ips) {
	double rmd= Par.m()-(Par-pi_0[ips].p()).m()-0.1421;
	if (bestrmd>rmd) bestrmd=rmd;
      }
    }

    t4->column("rmd", bestrmd);

    int hepid =0, hepid_mo=0, id=0;
    if (D.relation().genHepevt()){
      hepid=D.relation().genHepevt().idhep();
      if (p.relation().genHepevt()){
	Particle p_mc(p.relation().genHepevt());
	Particle D_mc(D.relation().genHepevt());
      } 
      if (D.relation().genHepevt().mother()) 
	hepid_mo=D.relation().genHepevt().mother().idhep();
    }
    if (p.relation().genHepevt())
      id=p.relation().genHepevt().idhep();
    
    t4->column("hepid", hepid); 
    t4->column("hepidm", hepid_mo); 
    t4->column("id", id);

    t4->dumpData();
  }
  
  for(int id=0; id<Bds.size(); ++id){
    
    Particle p=Bds[id].child(0);
    if (p.lund()!=443) p=Bds[id].child(2);

    float Mj =dynamic_cast<UserInfo&>(p.userInfo()).vmass();
    VectorL PP=p.p();
    double pj=pStar(PP,elec,posi).vect().mag();
    float dwd=0.03;
    float dmass = (Mj-3.097)/2./dwd;
    int bin = int(floor(dmass+0.5));

    Particle c0=p.child(0);
    Particle c1=p.child(1);
    Vector3 V = p.momentum().decayVertex()-ip_position;   

    double mulh1=Muid_mdst(c0.mdstCharged()).Muon_likelihood();
    double mulh2=Muid_mdst(c1.mdstCharged()).Muon_likelihood();
    double elh1= eid(c0.mdstCharged()).prob(3, -1, 5);
    double elh2= eid(c1.mdstCharged()).prob(3, -1, 5);
    double cth=pStar(p,elec,posi).vect().cosTheta();
    double lv1d2=cos(boostT(c1,p).vect().angle(boostT(Ups,p).vect()));
    
    Particle D1=Bds[id].child(1);
    Particle D2=Bds[id].child(2);
    if (D2.lund()==443) D2=Bds[id].child(0);
    
    int chD1 = dynamic_cast<UserInfo&>(D1.userInfo()).channel();
    float MD1 =dynamic_cast<UserInfo&>(D1.userInfo()).vmass();
    float MD_nom=D1.pType().mass();
    float dmassd1 = (MD1-MD_nom)/2./dwid[chD1-1];
    int bind1 = int(floor(dmassd1+0.5));
    int ibin1 = int(dmassd1+4.5);

    int chD2 = dynamic_cast<UserInfo&>(D2.userInfo()).channel();
    float MD2 =dynamic_cast<UserInfo&>(D2.userInfo()).vmass();
    float dmassd2 = (MD2-MD_nom)/2./dwid[chD2-1];
    int bind2 = int(floor(dmassd2+0.5));
    int ibin2 = int(dmassd2+4.5);

    if (abs(bind1)>5.5) continue;
    if (abs(bind2)>5.5) continue;

    double rm=(UPS-p.p()).m();
    double rmx=pStar(UPS-Bds[id].p()).m();
    double cthx=pStar(Bds[id].p()).vect().cosTheta();
    double px=pStar(Bds[id].p()).vect().mag();
    double mx=Bds[id].mass();
   
    double Mdd = (D1.p()+D2.p()).m();

    Vector3 VD1 = D1.momentum().decayVertex()-ip_position;   
    Vector3 VD2 = D2.momentum().decayVertex()-ip_position;   

    double min_pid=1;
    for(unsigned i=0; i<D1.nChildren(); i++)
      if (abs(D1.child(i).pType().lund())==211)
	if (min_pid>
	    atc_pid(3,1,5,2,3).prob(&(D1.child(i).mdstCharged())))
	  min_pid=
	    atc_pid(3,1,5,2,3).prob(&(D1.child(i).mdstCharged()));

    for(unsigned i=0; i<D2.nChildren(); i++)
      if (abs(D2.child(i).pType().lund())==211)
	if (min_pid>
	    atc_pid(3,1,5,2,3).prob(&(D2.child(i).mdstCharged())))
	  min_pid=
	    atc_pid(3,1,5,2,3).prob(&(D2.child(i).mdstCharged()));


    t3->column("mj", Mj);
    t3->column("bin", bin);
    t3->column("pj", pj);
    t3->column("rm", rm); 
    
    t3->column("vp", V.perp());
    t3->column("lhm", min(mulh1, mulh2));
    t3->column("lhe", min(elh1, elh2));
    
    t3->column("lv", lv1d2);
    t3->column("cth", cth);
     
    t3->column("r2", r2);
    t3->column("ea", All.e());
    t3->column("nt", ntrk);
    t3->column("ne", e_size);
    t3->column("nl", lep_size);
    t3->column("nk", nk);
    t3->column("em", Emax_gam);
    t3->column("f2", psi2_flag);
    
    t3->column("md1", MD1);
    t3->column("md2", MD2);
    t3->column("mpi", min_pid);
    t3->column("pd1", pStar(D1).vect().mag());
    t3->column("pd2", pStar(D1).vect().mag());
    t3->column("chd1", chD1);
    t3->column("chd2", chD2);
    t3->column("bin1", bind1);    
    t3->column("bin2", bind2);    

    t3->column("rmx", rmx);
    t3->column("mx", mx);
    t3->column("px", px);
    t3->column("cthx", cthx);
    t3->column("mdd", Mdd);

    t3->dumpData();

  }

}
#if defined(BELLE_NAMESPACE)
}
#endif

