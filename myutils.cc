#include "belle.h"
#include "panther/panther.h"
#include "particle/utility.h"
#include "particle/combination.h"
#include "mdst/mdst.h"
#include "ip/IpProfile.h"
#include "myutils.h"
#include "userinfo.h"
#include HEPEVT_H
#include MDST_H
//
// P.Koppenburg - 2003 08 29
//
//   modified by S. Nishida (for syst study)       07/05/23
//   modified by S. Nishida (for namespace Belle)  07/04/28
//   modified by S. Nishida (for gcc 3.2.2)
#include <iostream>
#include <iomanip>
#include <cmath>
//#include "pi0eta_prob.h"

using namespace std;


#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

  /******************* Make Photon_converted ****************/
void
makeGamma_converted(vector<Particle> &gam_c)
{
  Mdst_vee2_Manager &veeMgr = Mdst_vee2_Manager::get_manager();
  for(std::vector<Mdst_vee2>::iterator i = veeMgr.begin();
      i != veeMgr.end(); ++i )
    {
      if ( i -> kind() == 4 )  gam_c.push_back( Particle(*i) ) ;
      //    if ( i -> kind() == 3 ) alam.push_back( Particle(*i) ) ;
    }
}

/*************** Make Lambda *****************/
void
makeLambda(vector<Particle> &lam, vector<Particle> &alam)
{
  Mdst_vee2_Manager &veeMgr = Mdst_vee2_Manager::get_manager();
  for(std::vector<Mdst_vee2>::iterator i = veeMgr.begin();
      i != veeMgr.end(); ++i )
    {
      if ( i -> kind() == 2 )  lam.push_back( Particle(*i) ) ;
      if ( i -> kind() == 3 ) alam.push_back( Particle(*i) ) ;
    }
}
//*************** Make Proton *****************

void 
makeProton(vector<Particle> &pr, 
	   vector<Particle> &apr, 
	   const int flag){
  Mdst_charged_Manager &charged_mag = Mdst_charged_Manager::get_manager();
  Ptype ptype_proton("P+");
  Ptype ptype_aproton("AP+");
  
  for(std::vector<Mdst_charged>::iterator i = charged_mag.begin();
      i != charged_mag.end(); ++i){
    if(flag && !good_charged(*i))continue;
    if((*i).charge() > 0.){
      Particle tmp1(*i, ptype_proton);
      pr.push_back(tmp1);
    }else{
      Particle tmp1(*i, ptype_aproton);
      apr.push_back(tmp1);
    }
  }
}
//*************** Make Vertex Fit *****************/
void
doKvFit(Particle &p, unsigned flag) {
  if(! &p.userInfo() ) setUserInfo(p);
  dynamic_cast<UserInfo&>(p.userInfo()).vmass(p.mass());
  double v_chi=-1;
  if( dynamic_cast<UserInfo&>(p.userInfo()).vchisq() != -1. ) return;
  
  kvertexfitter kv;
  int nfit=0;
  for(unsigned i=0; i<p.nChildren(); i++) 
      addTrack2fit(kv,p.child(i));
  
  int err=kv.fit();
  if(!err){
    v_chi=kv.chisq()/kv.dgf();
    dynamic_cast<UserInfo&>(p.userInfo()).vchisq(kv.chisq()/kv.dgf());
    if(flag==0) makeMother(kv,p);
  }
  dynamic_cast<UserInfo&>(p.userInfo()).vchisq(v_chi);
}

void
doKvFit_nogamma(Particle &p, unsigned flag) {
  if(! &p.userInfo() ) setUserInfo(p);
    
  dynamic_cast<UserInfo&>(p.userInfo()).vchisq(999.);

  kvertexfitter kv;
  int nfit=0;
  for(unsigned i=0; i<p.nChildren(); i++) 
    if(p.child(i).mdstCharged()) {
      addTrack2fit(kv,p.child(i));
      nfit++;
    }
  
  int err=kv.fit();
  if(!err){
    VectorL PP;
    for(unsigned i=0; i<nfit; i++)
      PP += kv.momentum(i);
    
    for(unsigned i=0; i<p.nChildren(); i++)
      if(!p.child(i).mdstCharged()) PP += p.child(i).p();
    
    dynamic_cast<UserInfo&>(p.userInfo()).vchisq(kv.chisq()/kv.dgf());
    dynamic_cast<UserInfo&>(p.userInfo()).vmass(PP.m());
  }
}

void
doKvFit(vector<Particle> &plist, unsigned flag) {
  for(unsigned i=0;i<plist.size();++i) doKvFit(plist[i], flag);
}

void
doKvFit_nogamma(vector<Particle> &plist, unsigned flag) {
  for(unsigned i=0;i<plist.size();++i) doKvFit_nogamma(plist[i], flag);
}

/*************** Make Mass & Vertex Fit ******************************/
float
doKmvFit(Particle &p){
  if(! &p.userInfo() ) setUserInfo(p);
  //  if( dynamic_cast<UserInfo&>(p.userInfo()).mchisq() != -1. ) return;
  dynamic_cast<UserInfo&>(p.userInfo()).vmass(p.mass());
  kmassvertexfitter kmv;
  kmv.invariantMass(p.pType().mass());
  for(unsigned i=0; i<p.nChildren(); i++)
    addTrack2fit(kmv,p.child(i));
  int err=kmv.fit();
  if(!err){
    dynamic_cast<UserInfo&>(p.userInfo()).mchisq(kmv.chisq()/kmv.dgf());
    makeMother(kmv,p);
    return kmv.chisq()/kmv.dgf();
  }
  else {
    return 999;
      }
}

void
doKmvFit(vector<Particle> &plist, float chisq) {
  for(unsigned i=0;i<plist.size();++i)
    float chisq = doKmvFit(plist[i]);
}

void
doKmFit(Particle &p){
  dynamic_cast<UserInfo&>(p.userInfo()).vmass(p.mass());
  kmassfitter kmv;
  kmv.invariantMass(p.pType().mass());
  for(unsigned i=0; i<p.nChildren(); i++)
    addTrack2fit(kmv,p.child(i));
  int err=kmv.fit();
  if(!err){
    //    dynamic_cast<UserInfo&>(p.userInfo()).mchisq(kmv.chisq()/kmv.dgf());
    makeMother(kmv,p);
  } else {
    //    dynamic_cast<UserInfo&>(p.userInfo()).mchisq(999);
  }
}

void
doKmFit(vector<Particle> &plist) {
  for(unsigned i=0;i<plist.size();++i) doKmFit(plist[i]);
}


float
doKmvFit(Particle &p, vector<Particle> &ecl){
  kmassvertexfitter kmv;
  kmv.invariantMass(p.pType().mass());
  for(unsigned i=0; i<p.nChildren(); i++)
    addTrack2fit(kmv,p.child(i));
  if (abs(p.child(0).lund())==11 && abs(p.child(1).lund())==11)  {
    for(vector<Particle>::iterator k = ecl.begin();
	k != ecl.end(); k++){
      if((*k).momentum().p().vect().mag() <  0.001 ) continue;
      if((*k).momentum().p().vect().mag() >  5.000 ) continue;
      if((*k).momentum().p().vect().angle(p.child(0).momentum().p().vect()) 
	 < 0.05 || 
	 (*k).momentum().p().vect().angle(p.child(1).momentum().p().vect()) 
	 < 0.05 )  {
	addTrack2fit(kmv,*k);
      }
    }
  }
  
  int err=kmv.fit();
  if(!err){
    //    dynamic_cast<UserInfo&>(p.userInfo()).mchisq(kmv.chisq()/kmv.dgf());
    makeMother(kmv,p);
    return kmv.chisq()/kmv.dgf();
  }
  else {
    return 999;
  }
}


void
doKmvFit(vector<Particle> &plist, vector<Particle> &ecl) {
  for(unsigned i=0;i<plist.size();++i) {
    doKmvFit(plist[i],ecl);
  }
}

/*************** Make Fit with Beam  *********************************/
void
doKbFit(Particle &p){
  if(! &p.userInfo() ) setUserInfo(p);
  if( dynamic_cast<UserInfo&>(p.userInfo()).bchisq() != -1. ) return;
  kvertexfitter kbf;
  addTrack2fit(kbf, p);
  kbf.initialVertex(IpProfile::position());
  kbf.beamProfile(IpProfile::position_err_b_life_smeared());
  int err=kbf.fit();
  if(!err){
    dynamic_cast<UserInfo&>(p.userInfo()).bchisq(kbf.chisq()/kbf.dgf());
    p.momentum().decayVertex(kbf.vertex(), kbf.errVertex());
  }
}

void
doKbFit(vector<Particle> &plist){
  for(unsigned i=0;i<plist.size();++i) doKbFit(plist[i]);
}

/*************** Make Vertex ReFit for D* ****************************/
void
doReFit(Particle &Dstar, int flag){
  Dstar.momentum(Dstar.child(0).p()+Dstar.child(1).p());
  if(flag) return;
  setUserInfo(Dstar.child(1));
  if(! &Dstar.child(0).userInfo() ) return;
  if(dynamic_cast<UserInfo&>(Dstar.child(0).userInfo()).bchisq()==-1.) return;
  kvertexfitter kbf;
  addTrack2fit(kbf,Dstar.child(1));
  kbf.initialVertex(Dstar.child(0).momentum().decayVertex());
  kbf.knownVertex();
  int err=kbf.fit();
  if(!err){
    dynamic_cast<UserInfo&>(Dstar.child(1).userInfo()).bchisq(kbf.chisq());
    Dstar.momentum(Dstar.child(0).p()+kbf.momentum(0));
    Dstar.momentum().decayVertex(kbf.vertex(), kbf.errVertex());
  }
}

void
doReFit(vector<Particle> &plist, int flag){
  for(unsigned i=0;i<plist.size();++i) doReFit(plist[i], flag);
}

/************ Set error matrix (dpx) for pi0 ***************************/
void setPi0Error(Particle &p){
  if( !p.child(0) || !p.child(1) ) return;
  if( !p.child(0).mdstGamma() || 
      !p.child(1).mdstGamma()    ) return;
  for(int i=0; i<2; i++){
    double E     = p.child(i).mdstGamma().ecl().energy();
    double phi   = p.child(i).mdstGamma().ecl().phi();
    double theta = p.child(i).mdstGamma().ecl().theta();

    double E2  = E*E;
    double E4  = E2*E2;
    double ct2 = cos(theta)*cos(theta);
    double st2 = sin(theta)*sin(theta);
    double st4 = st2*st2;

    const HepPoint3D pivot(0.,0.,0.);
    HepMatrix  tmp_a(5,1);
    tmp_a[0][0] = 0.;
    tmp_a[1][0] = phi-M_PI/2;
    tmp_a[2][0] = 1/E/sin(theta);
    tmp_a[3][0] = 0.;
    tmp_a[4][0] = tan(M_PI/2-theta);
    HepVector  a(tmp_a);

    double errE     = p.child(i).mdstGamma().ecl().error(0);
    double errPHI   = p.child(i).mdstGamma().ecl().error(2);
    double errTHETA = p.child(i).mdstGamma().ecl().error(5);

    HepSymMatrix Ea(5,0);
    Ea[0][0] = 1.;
    Ea[1][1] = errPHI;
    Ea[2][2] = errE/E4/st2 + errTHETA*ct2/E2/st4;
    Ea[3][3] = 1.;
    Ea[4][2] = errTHETA*cos(theta)/E/st4;
    Ea[4][4] = errTHETA/st4;

    Helix helix(pivot, a, Ea);

    HepLorentzVector momentum;
    HepPoint3D position;
    HepSymMatrix error(7,0);

    momentum = helix.momentum(0.,0.,position,error);
    p.child(i).momentum().momentumPosition(momentum,position,error);
  }
  doKmvFit(p);
}

void setPi0Error(vector<Particle> &p_list){
  for(int i=0; i<p_list.size(); ++i)
    setPi0Error(p_list[i]);
}

/*********** Make 2- & 3- body particles ***********************************/
Particle
makeParticle(const Ptype &ptype,
	     Particle &p1, Particle &p2,
	     const double &width){
  Particle new_p;
  if(checkSame(p1,p2)) return new_p;
  VectorL P = p1.momentum().p()+p2.momentum().p();
  double mass      = P.mag();
  double type_mass = ptype.mass();
  if( (type_mass - width <= mass && mass <= type_mass + width) || !width){
    new_p=Particle( P, ptype );
    new_p.relation().append(p1);
    new_p.relation().append(p2);
  }
  return new_p;
}

Particle
makeParticle(const Ptype &ptype,
	     Particle &p1, Particle &p2, Particle &p3,
	     const double &width){
  Particle new_p;
  if(checkSame(p1,p2) || checkSame(p2,p3) || checkSame(p1,p3))return new_p;
  VectorL P = p1.momentum().p()+p2.momentum().p()+p3.momentum().p();
  double mass      = P.mag();
  double type_mass = ptype.mass();
  if((type_mass - width <= mass && mass <= type_mass + width) || !width){
    new_p=Particle( P, ptype );
    new_p.relation().append(p1);
    new_p.relation().append(p2);
    new_p.relation().append(p3);
  }
  return new_p;
}

VectorL
boostT(Particle p, Particle p_boost){ //p --boost--> p_boost;
  VectorL tmp(p.p());
  tmp.boost(-p_boost.p().boostVector());
  return tmp;
}

VectorL
boostT(VectorL p, VectorL p_boost){ //p --boost--> p_boost;
  VectorL tmp(p);
  tmp.boost(-p_boost.boostVector());
  return tmp;
}

void
removeParticleT(vector<Particle> &p_list, Particle p){
  if(p.nChildren())
    for(int i=0; i<p.nChildren(); i++)
      removeParticleT(p_list, p.child(i));
  else removeParticle(p_list, p);
}

void SetGenHepInfoLam(Particle &p)
{
  if (!p.child(0) || !p.child(1))
    return;
  if (!p.child(0).mdstCharged() || !p.child(1).mdstCharged())
    return;
  if (p.child(0).mdstCharged() ==  p.child(1).mdstCharged())
    return;
  
  const Gen_hepevt& hep0 = get_hepevt(p.child(0).mdstCharged());
  const Gen_hepevt& hep1 = get_hepevt(p.child(1).mdstCharged());

  if (hep0 && hep1 && (hep0.idhep() * hep1.idhep() > 0)) //we need 1 particle and 1 anti-particle                                                                                   
    return;

  if (hep0 && hep0.idhep() == 2212) //proton                                             
    p.child(0).relation().genHepevt(hep0); //proton relation set                          

  if (hep1 && hep1.idhep() == -211) //pion                                               
    p.child(1).relation().genHepevt(hep1); //pion relation set                           
 
  if (p.child(0).genHepevt() && p.child(1).genHepevt()) //gened guys exist                
    {
      if ((p.child(0).genHepevt() !=  p.child(1).genHepevt()) && //gened p0 != gened p1 
	  (p.child(0).genHepevt().mother() && p.child(1).genHepevt().mother()) && 
	  //their moms exist                                                   
        
	  (p.child(0).genHepevt().mother() == p.child(1).genHepevt().mother()) && 
	  //their mom is the same                                                      
	  
	  (p.child(0).genHepevt().mother().idhep() == 3122) && //their mom is lambda 
	  
	  (p.child(0).genHepevt().mother().daLast() - p.child(0).genHepevt().mother().daFirst() == 1)) //need to ask P. N.                                                         
	{
	  p.relation().genHepevt(p.child(0).genHepevt().mother()); //lambda relation set 
	}
    }
}

void SetGenHepInfoLam(std::vector<Particle> &p_list)
{
  for(unsigned int i=0; i<p_list.size(); ++i)
    SetGenHepInfoLam(p_list[i]);
}

void SetGenHepInfoALam(Particle &p)
{
  if (!p.child(0) || !p.child(1))
    return;
  if (!p.child(0).mdstCharged() || !p.child(1).mdstCharged())
    return;
  if (p.child(0).mdstCharged() ==  p.child(1).mdstCharged())
    return;

  const Gen_hepevt& hep0 = get_hepevt(p.child(0).mdstCharged());
  const Gen_hepevt& hep1 = get_hepevt(p.child(1).mdstCharged());

  if (hep0 && hep1 && (hep0.idhep() * hep1.idhep() > 0)) 
    //we need 1 particle and 1 anti-particle
    return;
 
  if (hep1 && hep1.idhep() == -2212) //proton                                          
    p.child(1).relation().genHepevt(hep0); //proton relation set                          
  
  if (hep0 && hep0.idhep() == 211) //pion                                          
    p.child(0).relation().genHepevt(hep1); //pion relation set                          
  
  if (p.child(0).genHepevt() && p.child(1).genHepevt()) //gened guys exist             
    {
      if ((p.child(0).genHepevt() !=  p.child(1).genHepevt()) && //gened p0 != gened p1	
	  (p.child(0).genHepevt().mother() && p.child(1).genHepevt().mother()) && 
	  //their moms exist
	  (p.child(0).genHepevt().mother() == p.child(1).genHepevt().mother()) && 
	  //their mom is the same
	  (p.child(1).genHepevt().mother().idhep() == -3122) && //their mom is lambda 
	  (p.child(1).genHepevt().mother().daLast() - p.child(1).genHepevt().mother().daFirst() == 1)) //need to ask P. N. (maybe anti changes!!) 
	
	{
	  p.relation().genHepevt(p.child(1).genHepevt().mother()); //lambda relation set
	}
    }
}

  void SetGenHepInfoALam(std::vector<Particle> &p_list)
  {
    for(unsigned int i=0; i<p_list.size(); ++i)
      SetGenHepInfoALam(p_list[i]);
  }


  void Pi0Veto(std::vector<Particle> &gamma, double npi0sigma){
    Mdst_ecl_aux_Manager &recAuxMgr = Mdst_ecl_aux_Manager::get_manager();
    
    for(std::vector<Particle>::iterator i=gamma.begin(); i!=gamma.end()-1; )
      {
        Mdst_ecl& eclGam = (*i).mdstGamma().ecl();
        Mdst_ecl_aux& shower(recAuxMgr(Panther_ID(eclGam.get_ID())));
        Particle Gam1(*i);
        int pi0veto = -1;
        for(std::vector<Particle>::iterator j=i+1; j!=gamma.end(); )
          {
            Mdst_ecl& eclGam = (*i).mdstGamma().ecl();
            Mdst_ecl_aux& shower(recAuxMgr(Panther_ID(eclGam.get_ID())));
            Particle Gam2(*j);
	    
            HepLorentzVector Gam1P4=Gam1.p();
            HepLorentzVector Gam2P4=Gam2.p();
	    
            HepLorentzVector Pi0P4 = Gam1P4 + Gam2P4;
            if(abs(Pi0P4.m()-0.134977)<npi0sigma*0.0068)
              {
                pi0veto=1;
		//           cout << "pi0 mass: " << Pi0P4.m() << " caused veto!"<<endl;                        
              }
	    
            ++j;
          }
        if(pi0veto>0)
          {
            gamma.erase(i);
          }
        else
          {
            ++i;
          }
      }

  }
  /*
  //==================================================================================================
  // Pi0 and eta prob
  //=================================================================================================//
  double Pi0_Prob(double m, double p2, double theta){   // mass of pair, p2 of second gamma. First gamma assumed p1>1.5
    return Pi0_Eta_Prob(1, m, p2, theta);
  }
  //================================================================================================//
  double Eta_Prob(double m, double p2, double theta){   // mass of pair, p2 of second gamma. First gamma assumed p1>1.5
    return Pi0_Eta_Prob(2, m, p2, theta);
  }
  //================================================================================================//
  double Pi0_Eta_Prob(int what, double m, double p2, double theta){   // mass of pair, p2 of second gamma. First gamma assumed p1>1.5
    double p_min, p_max, m_min, m_max;
    int p_bins, m_bins ;
    double *prob ;
    if (what==1){
#include "pi0_prob.h"
    } else if (what==2){
#include "eta_prob.h"
    } else  return -1.;
    // get row
    if (p2<=0.) return -2;
    double logp = log(1000*p2)/log(10.);
    if ((logp<p_min) || (logp>p_max)) return 0;
    double drow = p_bins*(logp-p_min)/(p_max-p_min);
    int row = int(drow);
    // get column
    if ((m<m_min) || (m>m_max)) return 0;
    double dcol = m_bins*(m-m_min)/(m_max-m_min);
    int col = int(dcol);
    if ((row<0) || (row>=p_bins) || (col<0) || (col>=m_bins)){
      std::cout << "Pi0 Eta prob: " << what << " " << logp << " " << m << " get coords " << row << " " << col << std::endl ;
      return -99.;
    }
    // get first order prob
    int pos = m_bins*row + col ;
    if (pos>m_bins*p_bins) {
      std::cout << "Pi0 Eta prob: " << what << " " << logp << " " << m << " get coords " << row << " " << col << " i.e. " << pos << std::endl ;
      return -99.;
    }
    double fop = prob[pos] ;
    //double fop = prob[1] ;
    //    for (int i=0;i<10;++i) std::cout<<"i="<<i<<"; p[i]="<<prob[i]<<endl;
    //    if (fop>1.) std::cout<<"ppi0="<<fop<<"; coord="<< pos <<"; for test pos[10]="<<prob[10]<<"; theta="<<theta<< std::endl ;
    //     std::cout << "Pi0 Eta prob: " << what << " " << logp<< " "  << m
    //	       << " get coords " << row << " " << col << " i.e. " << pos << " -> prob " << fop << std::endl ;
    return fop ;
  }



  //==================================================================================================
  // Closest probability Interface
  //=================================================================================================//
  double Closest_Pi0_Probability(const Mdst_gamma& IN,Mdst_gamma& OUT,double& mm,Mdst_gamma* VETO){
    return Closest_Probability(1,IN,OUT,mm,VETO);
  }
  //================================================================================================//
  double Closest_Eta_Probability(const Mdst_gamma& IN,Mdst_gamma& OUT,double& mm,Mdst_gamma* VETO){
    return Closest_Probability(2,IN,OUT,mm,VETO);
  }
  //=================================================================================================//
  double Closest_Pi0_Probability(const Hep3Vector& p3,Mdst_gamma& OUT,double& mm,Mdst_gamma* VETO){
    return Closest_Probability(1,p3,OUT,mm,VETO);
  }
  //================================================================================================//
  double Closest_Eta_Probability(const Hep3Vector& p3,Mdst_gamma& OUT,double& mm,Mdst_gamma* VETO){
    return Closest_Probability(2,p3,OUT,mm,VETO);
  }
  //==================================================================================================
  // Closest probability
  //=================================================================================================//
  double Closest_Probability(int what, const Mdst_gamma& IN, Mdst_gamma& OUT, double& mm, Mdst_gamma* VETO){

    HepLorentzVector Reco1( IN.px(), IN.py(), IN.pz()  ) ;
    Reco1= HepLorentzVector( Reco1.px(), Reco1.py(), Reco1.pz(), Reco1.rho()) ;
    double bp = -100.;
    mm = -1;

    Mdst_gamma_Manager& GammaMgr = Mdst_gamma_Manager::get_manager();

    for (std::vector<Mdst_gamma>::iterator it1= GammaMgr.begin(); it1 !=GammaMgr.end(); ++it1){
      Mdst_gamma* D = &(*it1);
      if (VETO!=NULL){ if (D->get_ID()==VETO->get_ID()) continue ;}
      if (D->get_ID()==IN.get_ID()) continue ;

      HepLorentzVector PG( D->px(), D->py(), D->pz()  ) ;
      PG = HepLorentzVector( PG.px(), PG.py(), PG.pz(), PG.rho()) ;
      double m = (Reco1 + PG).m();
      double p = Pi0_Eta_Prob(what, m, PG.e(), PG.theta() );

      if (p>bp) {
	bp = p ;
	mm = m ;
	OUT = *D ;
      }
    }

    return bp;
  }
  
  //=================================================================================================//
  double Closest_Probability(int what, const Hep3Vector& p3, Mdst_gamma& OUT, double& mm, Mdst_gamma* VETO){

    HepLorentzVector Reco1= HepLorentzVector( p3, p3.mag() ) ;
    double bp = -100.;
    mm = -1;

    Mdst_gamma_Manager& GammaMgr = Mdst_gamma_Manager::get_manager();

    for (std::vector<Mdst_gamma>::iterator it1= GammaMgr.begin(); it1 !=GammaMgr.end(); ++it1){
      Mdst_gamma* D = &(*it1);
      if (VETO!=NULL){ if (D->get_ID()==VETO->get_ID()) continue ;}

      HepLorentzVector PG( D->px(), D->py(), D->pz()  ) ;
      PG = HepLorentzVector( PG.px(), PG.py(), PG.pz(), PG.rho()) ;
      double m = (Reco1 + PG).m();
      double p = Pi0_Eta_Prob(what, m, PG.e(), PG.theta() );

      if (p>bp) {
	bp = p ;
	mm = m ;
	OUT = *D ;
      }
    }
    return bp;
  }
*/







#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

