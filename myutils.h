// P.Koppenburg - 2003 08 29
//     
//   modified by S.Nishida (accept 3-vector for syst. study) 07/05/23
//

#ifndef __PI0ETA_PROB__
#define __PI0ETA_PROB__

#include "belle.h"
#include "panther/panther.h"
#include MDST_H
#include "CLHEP/Vector/LorentzVector.h"


#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

using namespace std;
void makeGamma_converted  ( vector<Particle> &gam_c );
void makeLambda  ( vector<Particle> &lam, vector<Particle>  &alam );
void makeProton(vector<Particle> &pr, 
		vector<Particle> &apr, 
		const int = 1);

void doKvFit (Particle &p, unsigned flag=0);
void doKvFit_nogamma (Particle &p, unsigned flag=0);
void doKvFit (std::vector<Particle> &plist, unsigned flag=0);
void doKvFit_nogamma (std::vector<Particle> &plist, unsigned flag=0);

void doKmFit (Particle &p);
void doKmFit (std::vector<Particle> &plist);

float doKmvFit(Particle &p);
void doKmvFit(std::vector<Particle> &plist, float chisq);

float doKmvFit(Particle &p, std::vector<Particle> &ecl);
void doKmvFit(std::vector<Particle> &plist, std::vector<Particle> &ecl);

void doKbFit (Particle &p);
void doKbFit (std::vector<Particle> &plist);

void doReFit (Particle &p, int flag=0);
void doReFit (std::vector<Particle> &plist, int flag=0);

//void setGenHepInfoT(Particle &p);
//void setGenHepInfoT(std::vector<Particle> &p_list);

//void setGenHepInfoP(Particle &p);
//void setGenHepInfoP(std::vector<Particle> &p_list);

void setPi0Error(Particle &p);
void setPi0Error(std::vector<Particle> &p_list);

Particle makeParticle(const Ptype &ptype, 
		      Particle &p1, Particle &p2, 
		      const double &width = 0);
Particle makeParticle(const Ptype &ptype, 
		      Particle &p1, Particle &p2, Particle &p3, 
		      const double &width = 0);

VectorL boostT(Particle p, Particle p_boost);
VectorL boostT(VectorL p, VectorL p_boost);
 
void removeParticleT(std::vector<Particle> &p_list, Particle p);

void SetGenHepInfoLam(Particle &p);
void SetGenHepInfoLam(std::vector<Particle> &p_list);

void SetGenHepInfoALam(Particle &p);
void SetGenHepInfoALam(std::vector<Particle> &p_list);

void Pi0Veto(std::vector<Particle> &gamma, double npi0sigma);


  /*
    Probabilities given m     - the 2gamma mass, 
    E     - the 2nd gamma E, 
    theta - and the 2nd gamma polar angle
  */

  double Pi0_Prob(double m, double E, double theta);             
  double Eta_Prob(double m, double E, double theta);

  /*
    Returns highest Probability for IN   - high energy photon 
    VETO - This Mdst_gamma is excluded from the search
    returns:       - probability
    OUT  - low energy photon
    m    - mass for best combination
  */

  double Closest_Pi0_Probability(const Mdst_gamma& IN, Mdst_gamma& OUT, double& m, Mdst_gamma* VETO = NULL); 
  double Closest_Eta_Probability(const Mdst_gamma& IN, Mdst_gamma& OUT, double& m, Mdst_gamma* VETO = NULL); 
  
  double Closest_Pi0_Probability(const Hep3Vector& p3, Mdst_gamma& OUT, double& m, Mdst_gamma* VETO = NULL); 
  double Closest_Eta_Probability(const Hep3Vector& p3, Mdst_gamma& OUT, double& m, Mdst_gamma* VETO = NULL); 

  // internal usage
  double Pi0_Eta_Prob(int what, double m, double p2, double theta);
  double Closest_Probability(int,const Mdst_gamma&, Mdst_gamma&, double&, Mdst_gamma* ); 
  double Closest_Probability(int,const Hep3Vector &, Mdst_gamma&, double&, Mdst_gamma* ); 

  static const double ECL_F_THETA =  33 ; 
  static const double ECL_B_THETA = 128 ;
  static const double Degrees = 0.017453293 ;

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
