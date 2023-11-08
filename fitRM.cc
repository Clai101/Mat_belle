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
using namespace std;


#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
  
  void 
  fitRM(const Particle &p1, const Particle &p2, const VectorL UPS, const double M_fit, 
	VectorL &P_fit_1, VectorL &P_fit_2, double &chi2) {
    
    int ierr;
    double lambda=0;

    // cout << "*********************" << M_fit << endl;
    // first itteration : P_fit   <-- particle.p() 
    
    P_fit_1=p1.p();
    P_fit_2=p2.p();
    
    const double M_1=p1.pType().mass();
    const double M_2=p2.pType().mass();
    
    // fill momentum vector
    
    HepVector v_p1(3,0), v_p2(3,0);
    HepVector v_p10(3,0), v_p20(3,0);
    v_p1[0]=P_fit_1.x();
    v_p1[1]=P_fit_1.y();
    v_p1[2]=P_fit_1.z();
  
    v_p2[0]=P_fit_2.x();
    v_p2[1]=P_fit_2.y();
    v_p2[2]=P_fit_2.z();
  
    v_p10=v_p1;
    v_p20=v_p2;
  
    
    // fill momenutm error matrix
    
    HepSymMatrix err_p1(3,0), err_p2(3,0), unit(3,0);
    for(int im=0; im<3; ++im) {
      unit[im][im]=1;
      for(int in=0; in<3; ++in) { 
	err_p1[im][in]=p1.momentum().dp()[im][in];
	err_p2[im][in]=p2.momentum().dp()[im][in];
	if (abs(err_p1[im][in])<1.e-10) err_p1[im][in]=1.e-10;
	if (abs(err_p2[im][in])<1.e-10) err_p2[im][in]=1.e-10;
      }
    }
    
    // base a(0,1,2) 3-momentum of p1 ; a(4,5,6) - 3 momentum of p2 ; a(6) - lambda 
    HepSymMatrix A(7,0);
    HepVector a(7,0), diff_a(7,0);
    
    HepSymMatrix diff2_p1(3,0), diff2_p2(3,0);
    HepSymMatrix diff2_p1_p2(3,0), diff2_p2_p1(3,0);
    
    double diff;
    HepVector diff_p1(3,0), diff_p2(3,0);
    
    HepVector diff_chi_p1(3,0), diff_chi_p2(3,0);
    HepVector diff_m_p1(3,0), diff_m_p2(3,0);
    HepVector diff_lam_p1(3,0), diff_lam_p2(3,0);
    int n_fit=0;
    
    for(int itt=0; itt<100; ++itt) {
      
      n_fit++;
      double E_recoil=(UPS-P_fit_1-P_fit_2).e();
      double M_recoil=(UPS-P_fit_1-P_fit_2).m();

      double E_p1=P_fit_1.e();
      double E_p2=P_fit_2.e();

      diff=M_recoil* M_recoil - M_fit* M_fit;

      double diff_msdt= M_recoil - M_fit;
      if (abs(diff_msdt)<0.0001) continue;
      
      diff_chi_p1= 2. * err_p1.inverse(ierr) *  (v_p1-v_p10) ; 
      diff_m_p1= -2. * (( E_recoil / E_p1 + 1. ) * v_p1 + v_p2);
      
      diff_chi_p2= 2. * err_p2.inverse(ierr) *  (v_p2-v_p20) ; 
      diff_m_p2= -2. * (( E_recoil / E_p2 + 1. ) * v_p2 + v_p1);
      
      diff2_p1= 2. * err_p1.inverse(ierr) + 2. *
	lambda * ( E_recoil / E_p1 + 1. ) * unit;
    
      diff2_p2= 2. * err_p2.inverse(ierr) + 2. *  
	lambda * ( E_recoil / E_p2 + 1. ) * unit;
    
      diff2_p1_p2 = + 2. * lambda * unit; 
      diff2_p2_p1 = + 2. * lambda * unit; 
      
      diff_lam_p1= -1. * diff_m_p1;
      diff_lam_p2= -1. * diff_m_p2;
      
      
      for(int im=0; im<3; ++im) {
	a[im]=diff_chi_p1[im]-lambda*diff_m_p1[im];
	a[im+3]=diff_chi_p2[im]-lambda*diff_m_p2[im];
	A[im][6]=diff_lam_p1[im];
	A[6][im]=diff_lam_p1[im];
	for(int in=0; in<3; ++in) { 
	  A[im][in]=diff2_p1[im][in];
	  A[im+3][in+3]=diff2_p2[im][in];
	  A[im+3][in]=diff2_p1_p2[im][in];
	  A[im][in+3]=diff2_p2_p1[im][in];
	}
      }
      A[6][6]=0;
      a[6]=-diff;
      
      // use scale factor 0.3 for provide convergency (slow)
      diff_a = -0.3 * A.inverse(ierr) * a;
      
      for(int im=0; im<3; ++im) {
	diff_p1[im]=diff_a[im];
	diff_p2[im]=diff_a[im+3];
      }
      
      lambda+=diff_a[6];
      v_p1+=diff_p1;
      v_p2+=diff_p2;
      
      chi2 = ((v_p1-v_p10).T() * err_p1.inverse(ierr) * (v_p1-v_p10))[0] +
	((v_p2-v_p20).T() * err_p2.inverse(ierr) * (v_p2-v_p20))[0];

      //cout << itt << " diff=" << diff_msdt << " M_rec=" << M_recoil << " M_fit=" << M_fit << "; chi2=" << chi2 << endl;
      
      P_fit_1=VectorL(v_p1[0], v_p1[1], v_p1[2], 0);
      P_fit_2=VectorL(v_p2[0], v_p2[1], v_p2[2], 0);
      
      P_fit_1.setE(sqrt(P_fit_1.vect().mag2()+M_1*M_1));
      P_fit_2.setE(sqrt(P_fit_2.vect().mag2()+M_2*M_2));
      
    }

    // double M_recoil=(UPS-P_fit_1-P_fit_2).m();
    // cout<<"in fitrm: rmyf="<<M_recoil<< " n itter " << n_fit <<endl;

  }

#if defined(BELLE_NAMESPACE)
}
#endif
