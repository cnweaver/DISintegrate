#include <array>
#include <iostream>
#include <cmath>
#include <cassert>
#include <cctype>
#include <deque>
#include <fstream>

#include "LHAPDF/LHAPDF.h"

#if USE_GSL == USE_CUBPACKPP
    #error Exactly one of USE_GSL and USE_CUBPACKPP must be set
#endif

#if USE_GSL
    #include <gsl/gsl_math.h>
    #include <gsl/gsl_monte.h>
    #include <gsl/gsl_monte_plain.h>
    #include <gsl/gsl_monte_vegas.h>
#endif

#if USE_CUBPACKPP
    #include <Cubpack++/cubpack.h>
#endif

#if USE_PHOTOSPLINE
	#include <photospline/splinetable.h>
#endif

#include "AdaptiveQuad.h"
#include "cl_options.h"

extern "C" void qninit_();
extern "C" void grqinp_(double&,int&);
extern "C" void grxinp_(double&,int&);
extern "C" void grxlim_(int&,double&);
extern "C" void grqlim_(int&,double&,double&);
extern "C" void grxdef_(int&,double&);
extern "C" void grqdef_(int&,double&,double&);
extern "C" double qstfxq_(char*, char*, double&, double&, int&, long int, long int);
extern "C" void qnlinc_(int&, char*, int&, double[], long int);
extern "C" void qnbook_(int&, char*, long int);
extern "C" {
	void openhelper_(int&, char*, long int);
	void closehelper_(int&);
	void writefile_(int&, char*, long int);
	void qnvershelper_(bool&, int&, int&);
	void f_init(); // for fortran I/O initialisation. Called in qninit
	void f_exit(); // opposite of f_init
	double xfromix_(int&);
	int iqfromq_(double&);
	int nflget_(int&);
	double qfromiq_(int&);
	void qnfilw_(int&, int&);
	void qnpset_(char*, int&, int&, double&, long int);
	double qpdfxq_(char*, double&, double&, int&, long int);
	//void evolsg_(int&, int&, int&);
	void grxdef(int i, double x);
	void grqdef(int i, double x, double y);
	void qnrset_(char*, double&, long int);
	//void evolnp_(char*, int&, int&, int&, long int);
	//void evolnm_(char*, int&, int&, int&, long int);
	void evplus_(char*, int&, int&, int&, long int); //wrapper to write
	void qnlset_(char*, bool&, long int);
	void qniset_(char*, int&, long int);
	void qaddsi_(char*, int&, double&, long int);
	void qnpnul_(char*, long int);  //wrapper to write
	void qthres_(double&, double&);
}

void qnlinc(int i, std::string s, int j, double x[]){
	qnlinc_(i, (char*)s.c_str(), j, x, s.size());
}


double qstfxq(std::string s,std::string t,double x ,double y,int i){
	return(qstfxq_((char*)s.c_str(), (char*)t.c_str(), x, y, i, s.size(), t.size()));
};


void qnbook(int i, std::string s){
	qnbook_(i, (char*)s.c_str(), s.size());
}

void qnfilw(int i, int j){
	qnfilw_(i, j);
}

void qnpset(std::string s, int i, int j, double x){
	qnpset_((char*)s.c_str(), i, j, x, s.size());
}

void qaddsi(std::string s, int i, double j){
	qaddsi_((char*)s.c_str(), i, j ,s.size());
}

double xfromix(int i){
	return(xfromix_(i));
}

int iqfromq(double x){
	return(iqfromq_(x));
}

int nflget(int x){
	return(nflget_(x));
}

double qfromiq(int i){
	return(qfromiq_(i));
}

double qpdfxq(std::string s, double x, double y, int & i){
	return(qpdfxq_((char*)s.c_str(), x, y, i, s.size()));
}

//void evolsg(int i, int j, int k){
//	evolsg_(i, j, k);
//}

void grxdef(int i, double x){
	grxdef_(i, x);
}

void grqdef(int i, double x, double y){
	grqdef_(i, x, y);
}

// sucessfully tested
void qnrset(std::string s, double x){
	qnrset_((char*)s.c_str(), x, s.size());
}

//void evolnp(std::string s, int i, int j, int k){
//	evolnp_((char*)s.c_str(), i, j, k, s.size());
//}

//void evolnm(std::string s, int i, int j, int k){
//	evolnm_((char*)s.c_str(), i, j, k, s.size());
//}

void evplus(std::string s, int i, int j, int k){
	evplus_((char*)s.c_str(), i, j, k, s.size());
}

void qnlset(std::string s, bool l){
	qnlset_((char*)s.c_str(), l, s.size());
}

void qniset(std::string s, int i){
	qniset_((char*)s.c_str(), i, s.size());
}

void grqinp(double x, int i){
	grqinp_(x, i);
}

void grxinp(double x, int i){
	grxinp_(x, i);
}

void qthres(double x, double y){
	qthres_(x, y);
}

enum class InteractionType{
	CC,
	NC
};

std::ostream& operator<<(std::ostream& os, const InteractionType it){
	switch(it){
		case InteractionType::CC:
			return(os << "CC");
		case InteractionType::NC:
			return(os << "NC");
	}
}

std::istream& operator>>(std::istream& is, InteractionType& it){
	std::string raw;
	is >> raw;
	if(raw=="CC")
		it=InteractionType::CC;
	else if(raw=="NC")
		it=InteractionType::NC;
	else
		is.setstate(std::ios_base::failbit);
	return is;
}

enum class LeptonFlavor{
	None,
	Electron,
	Muon,
	Tau
};

std::ostream& operator<<(std::ostream& os, const LeptonFlavor lf){
	switch(lf){
		case LeptonFlavor::None:
			return(os << "None");
		case LeptonFlavor::Electron:
			return(os << "Electron");
		case LeptonFlavor::Muon:
			return(os << "Muon");
		case LeptonFlavor::Tau:
			return(os << "Tau");
	}
}

std::istream& operator>>(std::istream& is, LeptonFlavor& lf){
	std::string raw;
	is >> raw;
	std::transform(raw.begin(),raw.end(),raw.begin(),[](char c)->char{return(std::tolower(c));});
	if(raw=="none" || raw=="")
		lf=LeptonFlavor::Electron;
	else if(raw=="electron" || raw =="e")
		lf=LeptonFlavor::Electron;
	else if(raw=="muon" || raw=="mu")
		lf=LeptonFlavor::Muon;
	else if(raw=="tau")
		lf=LeptonFlavor::Tau;
	else
		is.setstate(std::ios_base::failbit);
	return is;
}

template<typename T>
std::string to_string(const T& t){
    std::ostringstream oss;
    oss << t;
    return oss.str();
}

static const double pi=4*atan(1);

class CrossSectionComputer{
public:
    ///Fraction of the target material which is protons
    double fProton;
    ///Fraction of the target material which is neutrons
    double fNeutron;
    ///Whether the incoming lepton is charged or not (a neutrino)
    bool chargedLepton;
	///The scattering process to calculate
	InteractionType process;
	///The flavor of the out-going lepton, used to determine its mass
	LeptonFlavor flavor;
	///lepton sign
	double ql;
	///Mass of the final-state lepton in GeV
	///May be zero for neutral-current neutrino interactions.
	double outgoingLeptonMass;
	///Whether to ignore contributions from kinematically forbidden regions of phase space
	bool enforceKinematics;
	
	///
	static constexpr double AlphaEM=7.570988599e-03;
	///
	static constexpr double GFermi=1.16637e-5;
	///W^{+/-} mass in GeV
	static constexpr double WBosonMass=80.398;
	///Z^0 mass in GeV
	static constexpr double ZBosonMass=91.1876;
    ///bottom quark mass in GeV
    static constexpr double BottomMass=4.20;
    ///top quark mass in GeV
    static constexpr double TopMass=171.2;
	
	///electron mass in GeV
	static constexpr double ElectronMass=0.00051099891;
	///muon mass in GeV
	static constexpr double MuonMass=0.105658367;
	///tau mass in GeV
	static constexpr double TauMass=1.77682;
	
	static constexpr double LeptonPolarization=0.;
	
	static constexpr double Sin2ThetaC=0.05;
	
	static constexpr double Sin2ThetaW=0.22308;
	
	static constexpr double CouplingVu=0.196;//0.203;
	static constexpr double CouplingVd=-0.346;//-0.351;
	static constexpr double CouplingVe=-0.03759600;//-0.00538;
    static constexpr double CouplingVnu=0.5;
	static constexpr double CouplingAu=0.5;
	static constexpr double CouplingAd=-0.5;
	static constexpr double CouplingAe=-0.5;
    static constexpr double CouplingAnu=0.5;
	
	static double getChargedLeptonMass(LeptonFlavor flavor){
		switch(flavor){
			case LeptonFlavor::None: return 0;
			case LeptonFlavor::Electron: return ElectronMass;
			case LeptonFlavor::Muon: return MuonMass;
			case LeptonFlavor::Tau: return TauMass;
		}
	}
	
	static double getOutgoingLeptonMass(bool chargedLepton, InteractionType process, LeptonFlavor flavor){
		if(not chargedLepton and process==InteractionType::NC)
			return 0;
		//TODO: CC case for charged leptons?
		return getChargedLeptonMass(flavor);
	}
	
	CrossSectionComputer(InteractionType process, double ql, bool chargedLepton, LeptonFlavor flavor):
    fProton(0.5),fNeutron(0.5), //isoscalar by default
	process(process),ql(ql),chargedLepton(chargedLepton),flavor(flavor),
	outgoingLeptonMass(getOutgoingLeptonMass(chargedLepton,process,flavor)),enforceKinematics(true)
	{}
    
    ///Set the combination of nucleons in the target material.
    ///\param fractionProton the ratio of the proton number density to the total number density
    ///\param fractionNeutron the ratio of the proton number density to the total number density
    ///\pre fractionProton+fractionNeutron == 1 (within a reasonable tolerance)
    void setNucleonFractions(double fractionProton, double fractionNeutron){
        if(std::abs((fractionProton+fractionNeutron)-1)>2*std::numeric_limits<double>::epsilon())
            throw std::logic_error("Nucleon number fractions do not sum to one");
        fProton=fractionProton;
        fNeutron=fractionNeutron;
    }
	
	void setEnforceKinematics(bool enforce){
		enforceKinematics=enforce;
	}
	
	///Compute the per-nucleon target mass
	///\return the rest energy in GeV
	double targetMass() const{
		return(fProton*0.938272 + fNeutron*0.93956);
	}

	///Compute the center of mass energy for a lepton-nucleon collision
	///\param En nucleon (proton) total energy in GeV
	///\param El lepton (electron) total energy in GeV
	///\return Mandlestam s (the square of CoM energy) in GeV^2
	double CalculateCentreOfMassEnergy(double En, double El) const{
		//std::cout << "Calculating Center-of-Mass Energy"  << std::endl;
		const double Mn=targetMass(); //GeV
		//std::cout << "Mn=" << Mn << std::endl;
		const double Ml=(chargedLepton?getChargedLeptonMass(flavor):0);
		double s=pow((En+El),2) - pow(sqrt(En*En-Mn*Mn)-sqrt(El*El-Ml*Ml),2);
		//std::cout << "Lepton Beam Energy: " << El  << " GeV" << std::endl;
		//std::cout << "Nucleon Beam Energy: " << En  << " GeV" << std::endl;
		//std::cout << "Center-of-Mass Energy: " << sqrt(s)  << " GeV"<< std::endl;
		return(s);
	}
	
	///Compute the energy of the incoming lepton in the lab frame, assuming the 
	///nucleon is stationary. 
	///\param s Mandlestam s (the square of CoM energy) in GeV^2
	double CalculateLabLeptonEnergy(double s) const{
		const double Ml=(chargedLepton?getChargedLeptonMass(flavor):0); //GeV
		const double Mn=targetMass(); //GeV
		return (s + Ml*Ml - Mn*Mn)/(2*Mn);
	}
	
	///Check whether a given point in phase space is physically realizable.
	///Based on equations 6-8 of http://dx.doi.org/10.1103/PhysRevD.66.113007
	///S. Kretzer and M. H. Reno
	///"Tau neutrino deep inelastic charged current interactions"
	///Phys. Rev. D 66, 113007
	///\param x Bjorken x of the interaction
	///\param y Bjorken y of the interaction
	///\param E Incoming neutrino in energy in the lab frame ($E_\nu$)
	///\param M Mass of the target nucleon ($M_N$)
	///\param m Mass of the secondary lepton ($m_\tau$)
	static bool kinematicallyAllowed(double x, double y, double E, double M, double m){
		if(x>1) //Eq. 6 right inequality
			return(false);
		if(x<((m*m)/(2*M*(E-m)))) //Eq. 6 left inequality
			return(false);
		//denominator of a and b
		double d=2*(1+(M*x)/(2*E));
		//the numerator of a (or a*d)
		double ad=1-m*m*((1/(2*M*E*x))+(1/(2*E*E)));
		double term=1-((m*m)/(2*M*E*x));
		//the numerator of b (or b*d)
		double bd=sqrt(term*term-((m*m)/(E*E)));
		return((ad-bd)<=d*y && d*y<=(ad+bd)); //Eq. 7
	}
	
	double TotalCrossSection(double s, double Q2min){
		//std::cout << "Computing with s=" << s << ", Q2min=" << Q2min << std::endl;
//		struct paramsContainer{
//			const CrossSectionComputer& comp;
//			double s; //Center of mass energy
//			double Q2min;
//			size_t calls;
//			paramsContainer(const CrossSectionComputer& comp, double s, double Q2min):
//			comp(comp),s(s),Q2min(Q2min),calls(0){}
//		} params(*this,s,Q2min);
//		auto gslWrapper=[](double* k, size_t dim, void* paramsp)->double{
//			paramsContainer& params=*static_cast<paramsContainer*>(paramsp);
//			double lgQ2min = log10(params.Q2min);
//			double lgs     = log10(params.s);
//			double lgQ2    = lgQ2min+k[0]*(lgs-lgQ2min);
//			double Q2      = pow(10.,lgQ2);
//			double x       = pow(10.,(lgQ2-lgs)*(1-k[1]));
//			double jacob = log(10.) * log(10.) * Q2 * x * (lgs - lgQ2min) * (lgs - lgQ2min) * (1-k[0]);
//			double result= jacob*params.comp.DoublyDifferentialCrossSection(x,Q2,params.s);
//			//std::cout << params.s << ' ' << params.Q2min
//			//<< ' ' << k[0] << ' ' << k[1]
//			//<< ' ' << lgQ2min << ' ' << lgs << ' ' << lgQ2 << ' ' << Q2 << ' ' << x
//			//<< ' ' << jacob << ' ' << result << '\n';
//			params.calls++;
//			return result;
//		};
//		typedef double(*gslFuncType)(double*,size_t,void*);
//		gsl_rng_env_setup();
//		const gsl_rng_type* T=gsl_rng_default;
//		gsl_rng* r=gsl_rng_alloc(T);
//		size_t calls = 1e4;
//		double result, error;
//		gsl_monte_vegas_state* state = gsl_monte_vegas_alloc(2);
//		gsl_monte_function G = { (gslFuncType)gslWrapper, 2, &params};
//		double xl[2] = {0.,0.};
//		double xu[2] = {1.,1.};
//		gsl_monte_vegas_integrate (&G, xl, xu, 2, calls, r, state, &result, &error);
//		//while(std::abs(state->chisq - 1.0) > 0.5)
//		do{
//			gsl_monte_vegas_integrate (&G, xl, xu, 2, calls, r, state, &result, &error);
//		}while(std::abs(gsl_monte_vegas_chisq(state) - 1.0) > 0.5);
//		
//		std::cout << ' ' << params.calls << " evaluations\n";
//		gsl_monte_vegas_free(state);
//		return(result);

#if USE_GSL
        struct paramsContainer{
            const CrossSectionComputer& comp;
            double s; //Center of mass energy
            double Q2min;
            size_t calls;
            paramsContainer(const CrossSectionComputer& comp, double s, double Q2min):
            comp(comp),s(s),Q2min(Q2min),calls(0){}
        } params(*this,s,Q2min);
        auto gslWrapper=[](double* k, size_t dim, void* paramsp)->double{
            paramsContainer& params=*static_cast<paramsContainer*>(paramsp);
            double x=pow(10.,k[0]);
            double y=pow(10.,k[1]);
            double Q2=params.s*x*y;
            if(Q2<params.Q2min)
                return 0;
            double jacobian=params.s*x*x*y*log(10.)*log(10.);
            double result=jacobian*params.comp.DoublyDifferentialCrossSection(x,Q2,params.s);
            params.calls++;
            return result;
        };
        typedef double(*gslFuncType)(double*,size_t,void*);
        gsl_rng_env_setup();
        const gsl_rng_type* T=gsl_rng_default;
        gsl_rng* r=gsl_rng_alloc(T);
        size_t maxCalls = 10000;
        double result, error;
        gsl_monte_vegas_state* state = gsl_monte_vegas_alloc(2);
        gsl_monte_function G = { (gslFuncType)gslWrapper, 2, &params};
        double cMin=log10(Q2min)-log10(s);
        double xl[2] = {cMin,cMin};
        double xu[2] = {0.,0.};
        //gsl_monte_vegas_integrate (&G, xl, xu, 2, maxCalls, r, state, &result, &error);
        //while(std::abs(state->chisq - 1.0) > 0.5)
        do{
            gsl_monte_vegas_integrate(&G, xl, xu, 2, maxCalls, r, state, &result, &error);
        }while(std::abs(gsl_monte_vegas_chisq(state) - 1.0) > 0.5);
        
        //std::cout << ' ' << params.calls << " evaluations\n";
        gsl_monte_vegas_free(state);
        return(result);
#endif //USE_GSL
		
#if USE_CUBPACKPP
		double cMin=log10(Q2min)-log10(s);
		TRIANGLE xyDomain(Point(cMin,0),Point(0,cMin),Point(0,0));
		size_t calls=0;
		auto integrand=[=,&calls](double lx, double ly)->double{
			double x=pow(10.,lx);
			double y=pow(10.,ly);
			double Q2=s*x*y;
			double jacobian=s*x*x*y*log(10.)*log(10.);
			double result=jacobian*this->DoublyDifferentialCrossSection(x,Q2,s);
			//std::cout << "Q2=" << Q2 << " x=" << x << " y=" << y << std::endl;
			//std::cout << x << ' ' << y << ' ' << result << '\n';
			calls++;
			return result;
		};
		double tolerance=1e-3;
		std::size_t maxCalls=1e7;
		double totalXS=Integrate(integrand,xyDomain,0,tolerance,maxCalls);
		//std::cout << ' ' << calls << " evaluations\n";
        return totalXS;
#endif //USE_CUBPACKPP
	}

	double SinglyDifferentialCrossSection(double s, double y, double Q2min){
		double cMin=log10(Q2min)-log10(s)-log10(y);
		if(cMin>=0)
			return(0);
		//std::cout << "s=" << s << " y=" << y << " Q2min=" << Q2min << " cMin=" << cMin << '\n';
		//unsigned int calls=0;
		auto integrand=[=/*,&calls*/](double lx){
			//calls++;
			double x=pow(10.,lx);
			double Q2=s*x*y;
			double jacobian=s*x*x*log(10.);
			double result=jacobian*this->DoublyDifferentialCrossSection(x,Q2,s);
			return result;
		};
		double dsdE=AdaptiveQuad::integrate(integrand,cMin,0,1e-4);
		//std::cout << '\t' << calls << " cross section evaluations\n";
		return dsdE;
	}
    
    double AverageSinglyDifferentialCrossSection(double s, double yMin, double yMax, double Q2min){
        if(yMin>=yMax || yMin<0)
            return 0;
        double lyMax=log10(yMax);
        double lyMin=log10(yMin);
        
        double cMin=log10(Q2min)-log10(s);
        if(lyMax<cMin)
            return 0;
        if(lyMin<cMin)
            lyMin=cMin;
        
        auto integrand=[=/*,&calls*/](double ly){
            //calls++;
            double y=pow(10.,ly);
            double jacobian=y*log(10.);
            double result=jacobian*this->SinglyDifferentialCrossSection(s, y, Q2min);
            return result;
        };
        
        //std::cout << "computing average dsdy over log(y) in [" << lyMin << ',' << lyMax << "]\n"; 
        double total=AdaptiveQuad::integrate(integrand,lyMin,lyMax,1e-4);
        return total/(yMax-yMin);
    }
	
	///Compute the doubly-differential cross section in x and Q^2
	///\param x the value of x at which to compute
	///\param Q2 the value of Q^2 at which to compute
	///\return the cross section in pb (?)
	double DoublyDifferentialCrossSection(double x, double Q2, double s) const{
		double y=Q2/(s*x);
		if(enforceKinematics && not kinematicallyAllowed(x, y, CalculateLabLeptonEnergy(s), targetMass(), outgoingLeptonMass))
			return 0;
		double redSigma=(process==InteractionType::CC ?
		                 CalculateCCReducedCrossSection(x,Q2,s) :
		                 CalculateNCReducedCrossSection(x,Q2,s));
		double prop=(process==InteractionType::CC ?
		             CalculateCCPropagator(x,Q2) :
		             CalculateNCPropagator(x,Q2));
		//std::cout << "\tredSigma=" << redSigma << " prop=" << prop << '\n';
		return(redSigma*prop);
	}
	
	double CalculateCCPropagator(double x, double Q2) const{
		double conversion_factor= 0.3894e9;  // GeV^2 -> pb^2
		return GFermi*GFermi*conversion_factor*pow(WBosonMass*WBosonMass/(Q2+WBosonMass*WBosonMass),2)/(2*pi*x);
	}
	
	double CalculateNCPropagator(double x, double Q2) const{
		double conversion_factor= 0.3894e9;  // GeV^2 -> pb^2
		//return GFermi*GFermi*conversion_factor*pow(ZBosonMass*ZBosonMass/(Q2+ZBosonMass*ZBosonMass),2)/(2*pi*x);
		return conversion_factor*2*pi*AlphaEM*AlphaEM/(x*Q2*Q2);
	}
	
	double CalculateEMPropagator(double x, double Q2, double s) const{
		double y=Q2/(s*x);
		double Yplus=(1+pow((1-y),2));
		double conversion_factor= 0.3894e9;  // GeV^2 -> pb^2
		return(conversion_factor*2*pi*Yplus*AlphaEM*AlphaEM/(x*(Q2*Q2)));
	}
	
	double CalculateCCReducedCrossSection(double x, double Q2, double s) const{
		double y=Q2/(s*x);
		double Yplus=(1+pow((1-y),2));
		double Yminus=(1-pow((1-y),2));
		//std::cout << '\t' << x << ' ' << Q2 << ' ' << s << ' ' << y << ' ' << Yplus << ' ' << Yminus << '\n';
		
        double flterm, f2term, f3term;
        if(chargedLepton){
            flterm=-0.5*y*y*FLCC(x,Q2);
            f2term=0.5*Yplus*F2CC(x,Q2);
            f3term=-0.5*ql*Yminus*XF3CC(x,Q2);
        }
        else{
            flterm=-y*y*FLCCnu(x,Q2);
            f2term=Yplus*F2CCnu(x,Q2);
            f3term=-ql*Yminus*XF3CCnu(x,Q2);
        }
		//std::cout << "\tfl=" << flterm << " f2=" << f2term << " f3=" << f3term << '\n';
        return ((1+LeptonPolarization*ql)*(flterm+f2term+f3term));
	}
	
	double CalculateNCReducedCrossSection(double x, double Q2, double s) const{
		double y=Q2/(s*x);
		double Yplus=(1+pow((1-y),2));
		double Yminus=(1-pow((1-y),2));
		
        double flterm, f2term, f3term;
        if(chargedLepton){
            flterm=-y*y*(FLNC(x,Q2)+LeptonPolarization*FLP(x,Q2))/Yplus;
            f2term=Yplus*F2NC(x,Q2)-ql*LeptonPolarization*F2P(x,Q2);
            f3term=-1*ql*Yminus*(XF3(x,Q2)-ql*LeptonPolarization*XF3P(x,Q2))/Yplus;
        }
        else{
            flterm=-y*y*(FLNCnu(x,Q2));
            f2term=Yplus*F2NCnu(x,Q2);
            f3term=-1*ql*Yminus*(XF3NCnu(x,Q2));
        }
        return (f2term+f3term+flterm);
	}
	
	double FLCCnu(double x, double Q2) const{
		//std::cout << "\tFLCCnu(" << x << ',' << Q2 << ")";
		const double cos2thc= 1-Sin2ThetaC;
		const double sin2thc=Sin2ThetaC;
		double weights[20]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		if(ql<0){ //neutrino
			weights[0]= 1./3.;
			weights[1]= 0.5*fProton + 0.5*fNeutron*cos2thc;
			weights[2]= 0.5*fProton*cos2thc + 0.5*fNeutron;
			weights[3]= 0.5*sin2thc;
			weights[6]=-0.5*fProton + 0.5*fNeutron*cos2thc;
			weights[7]= 0.5*fProton*cos2thc - 0.5*fNeutron;
			qnlinc(12,"CCEP2",3,weights);
			weights[0]=0.5;
			weights[1]=0.;
			weights[2]=0.;
			weights[3]=0.;
			weights[6]=-0.5*fProton + 0.5*fNeutron;
			weights[7]= 0.5*fProton - 0.5*fNeutron;
			qnlinc(12,"CCEP2",4,weights);
			qnlinc(12,"CCEP2",5,weights);
			
			int iflag=0;
			double result= qstfxq("FL","CCEP2",x,Q2,iflag);
			//std::cout << result << '\n';
			return result;
		}
		else{ //anti-neutrino
			weights[0]= 1./3.;
			weights[1]= 0.5*fProton + 0.5*fNeutron*cos2thc;
			weights[2]= 0.5*fProton*cos2thc + 0.5*fNeutron;
			weights[3]= 0.5*sin2thc;
			weights[6]= 0.5*fProton - 0.5*fNeutron*cos2thc;
			weights[7]=-0.5*fProton*cos2thc + 0.5*fNeutron;
			qnlinc(13,"CCEM2",3,weights);
			weights[0]=0.5;
			weights[1]=0.;
			weights[2]=0.;
			weights[3]=0.;
			weights[6]= 0.5*fProton - 0.5*fNeutron;
			weights[7]=-0.5*fProton + 0.5*fNeutron;
			qnlinc(13,"CCEM2",4,weights);
			qnlinc(13,"CCEM2",5,weights);
			
			int iflag=0;
			return qstfxq("FL","CCEM2",x,Q2,iflag);
		}
	}
	
	double F2CCnu(double x, double Q2) const{
		const double cos2thc= 1-Sin2ThetaC;
		const double sin2thc=Sin2ThetaC;
		double weights[20]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		if(ql<0){ //neutrino
			weights[0]= 1./3.;
			weights[1]= 0.5*fProton + 0.5*fNeutron*cos2thc;
			weights[2]= 0.5*fProton*cos2thc + 0.5*fNeutron;
			weights[3]= 0.5*sin2thc;
			weights[6]=-0.5*fProton + 0.5*fNeutron*cos2thc;
			weights[7]= 0.5*fProton*cos2thc - 0.5*fNeutron;
			qnlinc(12,"CCEP2",3,weights);
			weights[0]=0.5;
			weights[1]=0.;
			weights[2]=0.;
			weights[3]=0.;
			weights[6]=-0.5*fProton + 0.5*fNeutron;
			weights[7]= 0.5*fProton - 0.5*fNeutron;
			qnlinc(12,"CCEP2",4,weights);
			qnlinc(12,"CCEP2",5,weights);
			
			int iflag=0;
			return qstfxq("F2","CCEP2",x,Q2,iflag);
		}
		else{ //anti-neutrino
			weights[0]= 1./3.;
			weights[1]= 0.5*fProton + 0.5*fNeutron*cos2thc;
			weights[2]= 0.5*fProton*cos2thc + 0.5*fNeutron;
			weights[3]= 0.5*sin2thc;
			weights[6]= 0.5*fProton - 0.5*fNeutron*cos2thc;
			weights[7]=-0.5*fProton*cos2thc + 0.5*fNeutron;
			qnlinc(13,"CCEM2",3,weights);
			weights[0]=0.5;
			weights[1]=0.;
			weights[2]=0.;
			weights[3]=0.;
			weights[6]= 0.5*fProton - 0.5*fNeutron;
			weights[7]=-0.5*fProton + 0.5*fNeutron;
			qnlinc(13,"CCEM2",4,weights);
			qnlinc(13,"CCEM2",5,weights);
			
			int iflag=0;
			return qstfxq("F2","CCEM2",x,Q2,iflag);
		}
	}
	
	double XF3CCnu(double x, double Q2) const{
		const double cos2thc= 1-Sin2ThetaC;
		const double sin2thc=Sin2ThetaC;
		double weights[20]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		if(ql<0){ //neutrino
			weights[1]=-0.5*fProton + 0.5*fNeutron*cos2thc;
			weights[2]= 0.5*fProton*cos2thc - 0.5*fNeutron;
			weights[3]= 0.5*sin2thc;
			weights[6]= 0.5*fProton + 0.5*fNeutron*cos2thc;
			weights[7]= 0.5*fProton*cos2thc + 0.5*fNeutron;
			qnlinc(14,"CCEP3",3,weights);
			weights[1]=-0.5*(fProton-fNeutron);
			weights[2]= 0.5*(fProton-fNeutron);
			weights[3]= 0.5;
			weights[4]=-0.5;
			weights[6]= 0.5;
			weights[7]= 0.5;
			qnlinc(14,"CCEP3",4,weights);
			// Can't leave b->t here!
			weights[5]= 0.5;
			weights[0]= 0.5*1/5.;
			qnlinc(14,"CCEP3",5,weights);
			
			int iflag=0;
			return qstfxq("XF3","CCEP3",x,Q2,iflag);
		}
		else{ //anti-neutrino
			weights[1]= 0.5*fProton - 0.5*fNeutron*cos2thc;
			weights[2]=-0.5*fProton*cos2thc + 0.5*fNeutron;
			weights[3]=-0.5*sin2thc;
			weights[6]= 0.5*fProton + 0.5*fNeutron*cos2thc;
			weights[7]= 0.5*fProton*cos2thc + 0.5*fNeutron;
			qnlinc(15,"CCEM3",3,weights);
			weights[1]= 0.5*(fProton-fNeutron);
			weights[2]=-0.5*(fProton-fNeutron);
			weights[3]=-0.5;
			weights[4]= 0.5;
			weights[6]= 0.5;
			weights[7]= 0.5;
			qnlinc(15,"CCEM3",4,weights);
			// Can't leave b->t here!
			weights[5]=-0.5;
			weights[0]=-0.5*1/5.;
			qnlinc(15,"CCEM3",5,weights);
			
			int iflag=0;
			return qstfxq("XF3","CCEM3",x,Q2,iflag);
		}
	}
	
	double FLNCnu(double x, double Q2) const{
		double pz=Q2/((ZBosonMass*ZBosonMass+Q2)*(4*(Sin2ThetaW)*(1-Sin2ThetaW)));
		double Aup=pz*pz*(CouplingVnu*CouplingVnu+CouplingAnu*CouplingAnu)*(CouplingVu*CouplingVu+CouplingAu*CouplingAu);
		double Adown=pz*pz*(CouplingVnu*CouplingVnu+CouplingAnu*CouplingAnu)*(CouplingVd*CouplingVd+CouplingAd*CouplingAd);
		double weights[20]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		weights[0]=(Aup+2*Adown)*2./3.0;
		weights[1]=2*fProton*Aup   + 2*fNeutron*Adown;
		weights[2]=2*fProton*Adown + 2*fNeutron*Aup  ;
		weights[3]=2*Adown;
		qnlinc(23,"NCF2",3,weights);
		weights[0]=(Aup+Adown);
		weights[4]=2*Aup;
		qnlinc(23,"NCF2",4,weights);
		weights[0]=(4*Aup+6*Adown)/5;
		weights[5]=2*Adown;
		qnlinc(23,"NCF2",5,weights);
		
		int iflag=0;
		return qstfxq("FL","NCF2",x,Q2,iflag);
	}
	
	double F2NCnu(double x, double Q2) const{
		double pz=Q2/((ZBosonMass*ZBosonMass+Q2)*(4*(Sin2ThetaW)*(1-Sin2ThetaW)));
		double Aup=pz*pz*(CouplingVnu*CouplingVnu+CouplingAnu*CouplingAnu)*(CouplingVu*CouplingVu+CouplingAu*CouplingAu);
		double Adown=pz*pz*(CouplingVnu*CouplingVnu+CouplingAnu*CouplingAnu)*(CouplingVd*CouplingVd+CouplingAd*CouplingAd);
		double weights[20]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		weights[0]=(Aup+2*Adown)*2./3.0;
		weights[1]=2*fProton*Aup   + 2*fNeutron*Adown;
		weights[2]=2*fProton*Adown + 2*fNeutron*Aup  ;
		weights[3]=2*Adown;
		qnlinc(23,"NCF2",3,weights);
		weights[0]=(Aup+Adown);
		weights[4]=2*Aup;
		qnlinc(23,"NCF2",4,weights);
		weights[0]=(4*Aup+6*Adown)/5;
		weights[5]=2*Adown;
		qnlinc(23,"NCF2",5,weights);
		
		int iflag=0;
		return qstfxq("F2","NCF2",x,Q2,iflag);
	}
	
	double XF3NCnu(double x, double Q2) const{
		double pz=Q2/((ZBosonMass*ZBosonMass+Q2)*(4*(Sin2ThetaW)*(1-Sin2ThetaW)));
		double Bup=  4*CouplingAu*CouplingAnu*CouplingVu*CouplingVnu*pz*pz;
		double Bdown=4*CouplingAd*CouplingAnu*CouplingVd*CouplingVnu*pz*pz;
		double weights[20]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		weights[6]=2*fProton*Bup   + 2*fNeutron*Bdown;
		weights[7]=2*fProton*Bdown + 2*fNeutron*Bup;
		qnlinc(25,"NCF3",3,weights);
		qnlinc(25,"NCF3",4,weights);
		qnlinc(25,"NCF3",5,weights);
		
		int iflag=0;
		return qstfxq("XF3","NCF3",x,Q2,iflag);
	}
	
	///Always computes the NLO version
	double FLCC(double x, double Q2) const{
		if (ql>0){
			double weights[20]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
			weights[0]=0.5;
			weights[6]=-fProton*0.5 + fNeutron*0.5;
			weights[7]= fProton*0.5 - fNeutron*0.5;
			qnlinc(12,"CCEP2",3,weights);
			qnlinc(12,"CCEP2",4,weights);
			// actually you cannot just leave 5 the same as 4 for ccep2 because the b
			// will get added in by the weights(1), you have to get rid of it explicitly
			//    /*
			weights[0]=0.4;
			weights[1]=fProton*0.5 + fNeutron*0.5; //always 0.5
			weights[2]=fProton*0.5 + fNeutron*0.5; //always 0.5
			weights[3]=0.5;
			weights[4]=0.5;
			//    */
			qnlinc(12,"CCEP2",5,weights);
			
			int iflag=0;
			double result= qstfxq("FL","CCEP2",x,Q2,iflag);
			return result;
		} else {
			double weights[20]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
			weights[0]=0.5;
			weights[6]= fProton*0.5 - fNeutron*0.5;
			weights[7]=-fProton*0.5 + fNeutron*0.5;
			qnlinc(13,"CCEM2",3,weights);
			qnlinc(13,"CCEM2",4,weights);
			weights[0]=0.4;
			weights[1]=fProton*0.5 + fNeutron*0.5; //always 0.5
			weights[2]=fProton*0.5 + fNeutron*0.5; //always 0.5
			weights[3]=0.5;
			weights[4]=0.5;
			qnlinc(13,"CCEM2",5,weights);
			
			int iflag=0;
			return qstfxq("FL","CCEM2",x,Q2,iflag);
		}
	}
	
	///Always computes the NLO version
	double F2CC(double x, double Q2) const{
		if (ql>0){
			double weights[20]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
			weights[0]=0.5;
			weights[6]=-fProton*0.5 + fNeutron*0.5;
			weights[7]= fProton*0.5 - fNeutron*0.5;
			qnlinc(12,"CCEP2",3,weights);
			qnlinc(12,"CCEP2",4,weights);
			// actually you cannot just leave 5 the same as 4 for ccep2 because the b
			// will get added in by the weights(1), you have to get rid of it explicitly
			//
			weights[0]=0.4;
			weights[1]=fProton*0.5 + fNeutron*0.5; //always 0.5
			weights[2]=fProton*0.5 + fNeutron*0.5; //always 0.5
			weights[3]=0.5;
			weights[4]=0.5;
			// end lines not to be used
			qnlinc(12,"CCEP2",5,weights);
			
			int iflag=0;
			return qstfxq("F2","CCEP2",x,Q2,iflag);
		} else {
			double weights[20]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
			weights[0]=0.5;
			weights[6]= fProton*0.5 - fNeutron*0.5;
			weights[7]=-fProton*0.5 + fNeutron*0.5;
			qnlinc(13,"CCEM2",3,weights);
			qnlinc(13,"CCEM2",4,weights);
			weights[0]=0.4;
			weights[1]=fProton*0.5 + fNeutron*0.5; //always 0.5
			weights[2]=fProton*0.5 + fNeutron*0.5; //always 0.5
			weights[3]=0.5;
			weights[4]=0.5;
			qnlinc(13,"CCEM2",5,weights);
			
			int iflag=0;
			return qstfxq("F2","CCEM2",x,Q2,iflag);
		}
	}
	
	///Always computes the NLO version
	double XF3CC(double x, double Q2) const{
		double cos2thc= 1-Sin2ThetaC;
		double sin2thc=Sin2ThetaC;
		if(ql>0){
			double weights[20]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
			weights[1]=-fProton*0.5 + fNeutron*0.5*cos2thc;
			weights[2]= fProton*0.5*cos2thc - fNeutron*0.5;
			weights[3]=0.5*sin2thc;
			weights[0]=0.;
			weights[6]=fProton*0.5 + fNeutron*0.5*cos2thc;
			weights[7]=fProton*0.5*cos2thc + fNeutron*0.5;
			qnlinc(14,"CCEP3",3,weights);
			weights[1]=-fProton*0.5 - fNeutron*0.5;
			weights[2]= fProton*0.5 + fNeutron*0.5;
			weights[3]=0.5;
			weights[0]=0.;
			weights[4]=-0.5;
			weights[6]=fProton*0.5 + fNeutron*0.5;
			weights[7]=fProton*0.5 + fNeutron*0.5;
			qnlinc(14,"CCEP3",4,weights);
			// Leave b->t since it is mostly mass suppressed
			qnlinc(14,"CCEP3",5,weights);
			int iflag=0;
			return qstfxq("XF3","CCEP3",x,Q2,iflag);
		} else {
			double weights[20]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
			weights[1]= fProton*0.5 - fNeutron*0.5*cos2thc;
			weights[2]=-fProton*0.5*cos2thc + fNeutron*0.5;
			weights[3]=-0.5*sin2thc;
			weights[0]=0.;
			weights[6]=fProton*0.5 + fNeutron*0.5*cos2thc;
			weights[7]=fProton*0.5*cos2thc + fNeutron*0.5;
			qnlinc(15,"CCEM3",3,weights);
			weights[1]= fProton*0.5 - fNeutron*0.5;
			weights[2]=-fProton*0.5 + fNeutron*0.5;
			weights[3]=-0.5;
			weights[0]=0.;
			weights[4]=0.5;
			weights[6]=fProton*0.5 + fNeutron*0.5;
			weights[7]=fProton*0.5 + fNeutron*0.5;
			qnlinc(15,"CCEM3",4,weights);
			// Leave b->t since it is mostly mass suppressed
			qnlinc(15,"CCEM3",5,weights);
			int iflag=0;
			return qstfxq("XF3","CCEM3",x,Q2,iflag);
		}
	}
	
	double FLNC(double x, double Q2) const{
		double pz=Q2/((ZBosonMass*ZBosonMass+Q2)*(4*Sin2ThetaW*(1-Sin2ThetaW)));
		
		double Aup=4.0/9.0 -2*pz*(2./3.)*CouplingVu*CouplingVe
			+ pz*pz*(CouplingVe*CouplingVe+CouplingAe*CouplingAe)*(CouplingVu*CouplingVu+CouplingAu*CouplingAu);
		double Adown=1.0/9.0 -2*pz*(-1./3.)*CouplingVd*CouplingVe
			+ pz*pz*(CouplingVe*CouplingVe+CouplingAe*CouplingAe)*(CouplingVd*CouplingVd+CouplingAd*CouplingAd);
		double weights[20]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		weights[0]=(Aup+2*Adown)/3.0;
		weights[1] = fProton*Aup + fNeutron*Adown;
		weights[2] = fProton*Adown + fNeutron*Aup;
		weights[3] = Adown;
		qnlinc(23,"NCF2",3,weights);
		weights[0]=(Aup+Adown)/2;
		weights[4]=Aup;
		qnlinc(23,"NCF2",4,weights);
		weights[0]=(2.*Aup+3.*Adown)/5.;
		weights[5]=Adown;
		qnlinc(23,"NCF2",5,weights);
		
		int iflag=0;
		double result=qstfxq("FL","NCF2",x,Q2,iflag);
		return result;
	}

	
	///Always computes the NLO version
	double FLP(double x, double Q2) const{
		double pz=Q2/((ZBosonMass*ZBosonMass+Q2)*(4*(Sin2ThetaW)*(1-Sin2ThetaW)));
		
		double Aup=2*pz*(2./3.)*CouplingAe*CouplingVu
			- 2*pz*pz*CouplingAe*CouplingVe*(CouplingVu*CouplingVu+CouplingAu*CouplingAu);
		double Adown=2*pz*(-1./3.)*CouplingAe*CouplingVd
			- 2*pz*pz*CouplingAe*CouplingVe*(CouplingVd*CouplingVd+CouplingAd*CouplingAd);
		
		double weights[20]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		weights[0]=(Aup+2*Adown)/3.0;
		weights[1] = fProton*Aup + fNeutron*Adown;
		weights[2] = fProton*Adown + fNeutron*Aup;
		weights[3] = Adown;
		qnlinc(22,"NCP2",3,weights);
		weights[0]=(Aup+Adown)/2;
		weights[4]=Aup;
		qnlinc(22,"NCP2",4,weights);
		weights[0]=(2.*Aup+3.*Adown)/5.;
		weights[5]=Adown;
		qnlinc(22,"NCP2",5,weights);
		int iflag=0;
		return qstfxq("FL","NCP2",x,Q2,iflag);
	}
	
	///Always computes the NLO version
	double F2NC(double x, double Q2) const{
		double pz=Q2/((ZBosonMass*ZBosonMass+Q2)*(4*(Sin2ThetaW)*(1-Sin2ThetaW)));
		
		double Aup=4.0/9.0 -2*pz*(2.0/3.0)*CouplingVu*CouplingVe
			+ pz*pz*(CouplingVe*CouplingVe+CouplingAe*CouplingAe)*(CouplingVu*CouplingVu+CouplingAu*CouplingAu);
		double Adown=1.0/9.0 -2*pz*(-1.0/3.0)*CouplingVd*CouplingVe
			+ pz*pz*(CouplingVe*CouplingVe+CouplingAe*CouplingAe)*(CouplingVd*CouplingVd+CouplingAd*CouplingAd);
		double weights[20]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		weights[0]=(Aup+2*Adown)/3.0;
		weights[1] = fProton*Aup + fNeutron*Adown;
		weights[2] = fProton*Adown + fNeutron*Aup;
		weights[3] = Adown;
		qnlinc(23,"NCF2",3,weights);
		weights[0]=(Aup+Adown)/2;
		weights[4]=Aup;
		qnlinc(23,"NCF2",4,weights);
		weights[0]=(2.*Aup+3.*Adown)/5.;
		weights[5]=Adown;
		qnlinc(23,"NCF2",5,weights);
		int iflag=0;
		return qstfxq("F2","NCF2",x,Q2,iflag);
	}
	
	///Always computes the NLO version
	double F2P(double x, double Q2) const{
		double pz= Q2/((ZBosonMass*ZBosonMass+Q2)*(4*Sin2ThetaW*(1-Sin2ThetaW)));
		
		double Aup=2*pz*(2.0/3.0)*CouplingAe*CouplingVu
			- 2*pz*pz*CouplingAe*CouplingVe*(CouplingVu*CouplingVu+CouplingAu*CouplingAu);
		double Adown=2*pz*(-1.0/3.0)*CouplingAe*CouplingVd
			- 2*pz*pz*CouplingAe*CouplingVe*(CouplingVd*CouplingVd+CouplingAd*CouplingAd);
		
		double weights[20]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		weights[0]=(Aup+2*Adown)/3.0;
		weights[1] = fProton*Aup + fNeutron*Adown;
		weights[2] = fProton*Adown + fNeutron*Aup;
		weights[3] = Adown;
		qnlinc(22,"NCP2",3,weights);
		weights[0]=(Aup+Adown)/2;
		weights[4]=Aup;
		qnlinc(22,"NCP2",4,weights);
		weights[0]=(2.*Aup+3.*Adown)/5.;
		weights[5]=Adown;
		qnlinc(22,"NCP2",5,weights);
		int iflag=0;
		return qstfxq("F2","NCP2",x,Q2,iflag);
	}
	
	///Always computes the NLO version
	double XF3(double x, double Q2) const{
		double pz=Q2/((ZBosonMass*ZBosonMass+Q2)*(4*Sin2ThetaW*(1-Sin2ThetaW)));
		
		double Bup=-2*(2./3.)*CouplingAu*CouplingAe*pz
			+ 4*CouplingAu*CouplingAe*CouplingVu*CouplingVe*pz*pz;
		double Bdown=-2*(-1./3.)*CouplingAd*CouplingAe*pz
			+ 4*CouplingAd*CouplingAe*CouplingVd*CouplingVe*pz*pz;
		
		double weights[20]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		weights[6]= fProton*Bup + fNeutron*Bdown;
		weights[7]= fProton*Bdown + fNeutron*Bup;
		
		qnlinc(25,"NCF3",3,weights);
		qnlinc(25,"NCF3",4,weights);
		qnlinc(25,"NCF3",5,weights);
		
		int iflag=0;
		return qstfxq("XF3","NCF3",x,Q2,iflag);
	}
	
	///Always computes the NLO version
	double XF3P(double x, double Q2) const{
		double pz=Q2/((ZBosonMass*ZBosonMass+Q2)*(4*Sin2ThetaW*(1-Sin2ThetaW)));
		
		double Bup=2*(2.0/3.0)*CouplingAu*CouplingVe*pz
			- 2*CouplingAu*CouplingVu*pz*pz*(CouplingVe*CouplingVe+CouplingAe*CouplingAe);
		double Bdown=2*(-1.0/3.0)*CouplingAd*CouplingVe*pz
			- 2*CouplingAd*CouplingVd*pz*pz*(CouplingVe*CouplingVe+CouplingAe*CouplingAe);
		
		double weights[20]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		weights[6]= fProton*Bup + fNeutron*Bdown;
		weights[7]= fProton*Bdown + fNeutron*Bup;
		
		qnlinc(24,"NCP3",3,weights);
		qnlinc(24,"NCP3",4,weights);
		qnlinc(24,"NCP3",5,weights);
		
		int iflag=0;
		return qstfxq("XF3","NCP3",x,Q2,iflag);
	}

};

void qcd_init(InteractionType process, std::string pdfName, int pdfSubset){
	LHAPDF::initPDFSet(pdfName,pdfSubset);
	
	qninit_();
	
	double mc2= pow(LHAPDF::getQMass(4),2);
	double mb2= pow(LHAPDF::getQMass(5),2);
	
	int Ix=409;
	int Iq=200; // was 61
//	double xlim2=97e-8;
//	double qlim1=29e-2;
//	double qlim2=200e3;
	double xlim2=5.e-13;
	double qlim1=5.e-1;
	double qlim2=2.e12;
	
	grxdef(Ix,xlim2);
	grqdef(Iq,qlim1,qlim2);
	double qmore=50.0;
	grqinp(qmore,1);
	grqinp(mc2,1);
	grqinp(mb2,1);
	qthres(mc2,mb2);
	//qthres(2.e12,4.e12); //effectively disable heavy quarks

	double Mz= CrossSectionComputer::ZBosonMass;
	double alfas=LHAPDF::alphasPDF(Mz);
	qnrset("ALFQ0",Mz*Mz);
	qnrset("ALFAS",alfas);
	qnrset("AAAR2",1.);
	qnrset("BBBR2",0.);
	qnrset("AAM2L",1.);
	qnrset("BBM2L",0.);
	qnrset("AAM2H",1.);
	qnrset("BBM2H",0.);
	qnrset("MCSTF",sqrt(mc2));
	qnrset("MBSTF",sqrt(mb2));
	qnrset("MCALF",sqrt(mc2));
	qnrset("MBALF",sqrt(mb2));
	
	//TODO why are these flags not set anymore?
	//qnlset("WTF2C",true);
	qnlset("W2STF",true);
	//qnlset("WTF2B",true) ;
	qnlset("CLOWQ",false) ;
	//qnlset("WTFLC",true) ;
	//qnlset("WTFLB",true) ;
	
	qniset("ORDER",2) ; //use NLO
	
	qnfilw(0,0);
	
	qnbook(2,"UPLUS");
	qnbook(3,"DPLUS");
	qnbook(4,"SPLUS");
	qnbook(5,"CPLUS");
	qnbook(6,"BPLUS");
	qnbook(7,"UPVAL"); //up^-
	qnbook(8,"DNVAL"); //down^-
	
	for(int my_iq=1; my_iq<=Iq; my_iq++){
		double nf=0.0; //TODO: should be int?
		double my_q=qfromiq(my_iq);
		for(int mif=1 ;mif<7; mif++){
			if(LHAPDF::xfx(0.01, sqrt(my_q), mif)>0.000000000001)
				nf=nf+1;
		}
		for(int mix=1; mix<=Ix; mix++){
			double mx=xfromix(mix);
			double singlet =0.0;
			//TODO: both branches same?
			if (process==InteractionType::CC){
				singlet=(LHAPDF::xfx(mx,sqrt(my_q), 1)+LHAPDF::xfx(mx,sqrt(my_q), 2)+LHAPDF::xfx(mx,sqrt(my_q), 3)+LHAPDF::xfx(mx,sqrt(my_q), 4)+LHAPDF::xfx(mx,sqrt(my_q), 5)+LHAPDF::xfx(mx,sqrt(my_q), -1)+LHAPDF::xfx(mx,sqrt(my_q), -2)+LHAPDF::xfx(mx,sqrt(my_q), -3)+LHAPDF::xfx(mx,sqrt(my_q), -4)+LHAPDF::xfx(mx,sqrt(my_q), -5));
			} else{
				singlet=(LHAPDF::xfx(mx,sqrt(my_q), 1)+LHAPDF::xfx(mx,sqrt(my_q), 2)+LHAPDF::xfx(mx,sqrt(my_q), 3)+LHAPDF::xfx(mx,sqrt(my_q), 4)+LHAPDF::xfx(mx,sqrt(my_q), 5)+LHAPDF::xfx(mx,sqrt(my_q), -1)+LHAPDF::xfx(mx,sqrt(my_q), -2)+LHAPDF::xfx(mx,sqrt(my_q), -3)+LHAPDF::xfx(mx,sqrt(my_q), -4)+LHAPDF::xfx(mx,sqrt(my_q), -5));
			}
			//TODO: here negative values are replaced by 0. Why is this not done above for the singlet?
			qnpset("SINGL",mix,my_iq,singlet);
			qnpset("UPLUS",mix,my_iq,std::max(LHAPDF::xfx(mx, sqrt(my_q), 2),0.0)+std::max(LHAPDF::xfx(mx, sqrt(my_q), -2),0.0)-(1.0/(nf))*singlet);
			qnpset("DPLUS",mix,my_iq,std::max(LHAPDF::xfx(mx, sqrt(my_q), 1),0.0)+std::max(LHAPDF::xfx(mx, sqrt(my_q), -1),0.0)-(1.0/(nf))*singlet);
			qnpset("SPLUS",mix,my_iq,std::max(LHAPDF::xfx(mx, sqrt(my_q), 3),0.0)+std::max(LHAPDF::xfx(mx, sqrt(my_q), -3),0.0)-(1.0/(nf))*singlet);
			if (nf>3.01){
				qnpset("CPLUS",mix,my_iq,std::max(LHAPDF::xfx(mx, sqrt(my_q), 4),0.0)+std::max(LHAPDF::xfx(mx, sqrt(my_q), -4),0.0)-(1.0/(nf))*singlet);
			} else{
				qnpset("CPLUS",mix,my_iq,0.0);
			}
			if (nf>4.01){
				qnpset("BPLUS",mix,my_iq,std::max(LHAPDF::xfx(mx, sqrt(my_q), 5),0.0)+std::max(LHAPDF::xfx(mx, sqrt(my_q), -5),0.0)-(1.0/(nf))*singlet);
			} else{
				qnpset("BPLUS",mix,my_iq,0.0);
			}
			qnpset("GLUON",mix,my_iq,std::max(LHAPDF::xfx(mx,sqrt(my_q), 0),0.0));
			qnpset("UPVAL",mix,my_iq,std::max(LHAPDF::xfx(mx, sqrt(my_q), 2),0.0)-std::max(LHAPDF::xfx(mx, sqrt(my_q), -2),0.0));
			qnpset("DNVAL",mix,my_iq,std::max(LHAPDF::xfx(mx, sqrt(my_q), 1),0.0)-std::max(LHAPDF::xfx(mx, sqrt(my_q), -1),0.0));
		}
	}
}

double extractUnits(std::string raw){
	if(raw=="picobarn")
		return 1.;
	if(raw=="millibarn")
		return 1.e-9;
	if(raw=="barn")
		return 1.e-12;
	if(raw=="m2")
		return 1.e-40;
	if(raw=="cm2")
		return 1.e-36;
	try{
		double value=std::stod(raw);
		return value;
	}catch(std::logic_error& ex){
		throw std::logic_error(std::string("Unrecognized unit: "+raw+": "+ex.what()));
	}
}

enum class OutputFormat{
	Unset, 
	ASCII,
	Spline
};

std::ostream& operator<<(std::ostream& os, const OutputFormat of){
	switch(of){
		case OutputFormat::Unset:
			return(os << "Unset");
		case OutputFormat::ASCII:
			return(os << "ASCII");
		case OutputFormat::Spline:
			return(os << "Spline");
	}
}
std::istream& operator>>(std::istream& is, OutputFormat& of){
	std::string raw;
	is >> raw;
	std::transform(raw.begin(),raw.end(),raw.begin(),[](char c)->char{return(std::tolower(c));});
	if(raw=="unset" || raw=="")
		of=OutputFormat::Unset;
	else if(raw=="ascii")
		of=OutputFormat::ASCII;
	else if(raw=="spline")
		of=OutputFormat::Spline;
	else
		is.setstate(std::ios_base::failbit);
	return is;
}

///generate a sequence of values of the form baseValue + index * step
///\param startIndex the first index (inclusive)
///\param endIndex the last index (inclusive)
///\param baseValue the value at index 0
///\param step the difference between values at adjacent indices
std::vector<double> generateSequence(int startIndex, int endIndex, double baseValue, double step){
	std::vector<double> seq;
	for(int idx=startIndex; idx<=endIndex; idx++)
		seq.push_back(baseValue+idx*step);
	return seq;
}

void padKnots(std::vector<double>& knots, const unsigned int order){
	assert(knots.size()>=2);
	double lowStep=knots[1]-knots[0];
	for(unsigned int i=0; i<order; i++)
		knots.insert(knots.begin(),knots.front()-lowStep);
	double highStep=knots[knots.size()-1]-knots[knots.size()-2];
	for(unsigned int i=0; i<order; i++)
		knots.push_back(knots.back()+highStep);
}

int main(int argc, char* argv[]){
	std::string pdfName="HERAPDF15NLO_EIG.LHgrid";
	int pdfSubset=0;
	auto interaction=InteractionType::CC;
	double ql=-1; //-1=neutrino/electron, 1=antineutrino/positron
	std::string neutrinoTypeString="nu"; //should be kept consistent with ql
	double Q2Min=1;
	bool computeTotal=false;
	std::string totalFile="Total_XS.dat";
	bool computeSingle=false;
	std::string singleFile="dsdE.dat";
	bool computeDouble=false;
	std::string doubleFile="dsdxdy.fits";
	double outputUnits=1.;
	double protonFraction=0.5;
	double ElMin=1e1;
	double ElMax=1e12;
	unsigned int ElSteps=220;
	unsigned int ElKnotSteps=75;
	unsigned int lxySteps=200;
	unsigned int xyKnotSteps=85;
	unsigned int zKnotSteps=100;
	OutputFormat format=OutputFormat::ASCII;
	bool singleDiffPeakHack=false;
	bool squareSingleTable=false;
    bool writeSplineSlices = true;
	double splineRegSingle=1e-10;
	LeptonFlavor leptonFlavor=LeptonFlavor::None; //original calculation always ignored this
	bool enforceKinematics=true; //original calculation did not enforce this
	
	OptionParser op;
	op.setBaseUsage("Compute cross sections for deep-inelastic scattering of neutrinos on stationary nuclear targets.\n"
	                "All energy units are GeV, and neutrino flavor differences are neglected.");
	op.addOption({"f","format"},format,"Output data format (ASCII/Spline)");
	op.addOption({"i","interaction"},interaction,"The interaction process (CC/NC)");
	op.addOption<std::string>({"p","particle"},
		[&](std::string raw){
			if(raw=="nu" || raw=="neutrino"){
				ql=-1;
				neutrinoTypeString="nu";
			}
			else if(raw=="anu" || raw=="nubar" || raw=="antineutrino"){
				ql=1;
				neutrinoTypeString="nubar";
			}
			else
				throw std::logic_error("Invalid incoming particle type: "+raw);
		},"Particle type (nu/nubar) to use as the incoming particle (default: "+neutrinoTypeString+")");
	op.addOption("pdf", pdfName, "Name of the LHAPDF PDFSet to use");
	op.addOption("pdfSubset", pdfSubset, "Subset of the selected PDFSet to use");
	op.addOption({"t","total"},[&]{computeTotal=true;},"Compute the total cross section (default: false)");
	op.addOption({"s","singly-differential"},[&]{computeSingle=true;},"Compute the singly-differential cross section (default: false)");
	op.addOption({"d","doubly-differential"},[&]{computeDouble=true;},"Compute the doubly-differential cross section (default: false)");
	op.addOption("totalXSOutput",totalFile,"The file to which the calculated total cross section data should be written");
	op.addOption("singleXSOutput",singleFile,"The file to which the calculated singly-differential cross section data should be written");
	op.addOption("doubleXSOutput",doubleFile,"The file to which the calculated doubly-differential cross section data should be written");
	op.addOption({"protonFraction"},protonFraction,"The number fraction of the target nuclei which consists of protons");
	op.addOption<std::string>({"u","unit"},[&](std::string s){outputUnits=extractUnits(s);},
	                           "Units to use for output. May be a string (picobarn,millibarn,barn,m2,cm2)\n"
	                           "or a numerical factor relative to picobarns (e.g. 1e-40 for m^2)\n"
	                           "(default: picobarns)");
	op.addOption("leptonFlavor",leptonFlavor,"The out-going lepton whose mass should be accounted for in the kinematic limits\n"
	             "of integration. May be electron, muon, tau, or empty/unset to neglect this and treat the out-going\n"
	             "lepton as massless.");
	op.addOption("enforceKinematics",enforceKinematics,
	             "Whether to ignore cross section contributions from kinematically-forbiden regions\n"
	             "of phase space");
	op.addOption("Q2min",Q2Min,"The minimum Q^2 value to consider in the calculation\n"
	             "(This implcitly determines the minimum x and y values considered)");
	op.addOption("Emin",ElMin,"The minimum incoming neutrino energy for which to compute cross sections (GeV)");
	op.addOption("Emax",ElMax,"The maximum incoming neutrino energy for which to compute cross sections (GeV)");
	op.addOption("Esteps",ElSteps,"The number of (logarithmic) steps in incoming neutrino energy for which to compute cross sections");
	op.addOption("Eknots",ElKnotSteps,"The number of knots to use in the incoming neutrino energy dimension for spline output.\n"
	             "Knots will be added for padding to ensure support over the full physical region");
	op.addOption("xysteps",lxySteps,"The number of (logarithmic) steps in x and y for which to compute cross sections");
	op.addOption("xyknots",xyKnotSteps,"The number of knots to use in the x and y dimensions for spline output.\n"
	             "Knots will be added for padding to ensure support over the full physical region");
	op.addOption("zknots",zKnotSteps,"The number of knots to use in the z dimensions for singly-differential square spline output.\n"
	             "Knots will be added for padding to ensure support over the full physical region");
	op.addOption("singleDiffPeakHack",singleDiffPeakHack,"At high energies, when E_out == E_l, replace the singly-differential cross section\n"
	             "with a value which will yield the correct average value with linear interpolation over the last two\n"
	             "tabulated E_out values.\n"
	             "See https://arxiv.org/abs/2112.13804 appendix D for details.");
	op.addOption("squareSingleTable",squareSingleTable,"Tabulate singly-differential cross section in E_l and\n"
	             "z = (Eout - Emin)/(Ein - Emin), to produce a fully populated square table.\n"
	             "See https://arxiv.org/abs/2112.13804 appendix D for details.");
	op.addOption({"rs","single-regularization"},splineRegSingle,"Regularization strength for singly-differential spline fit.");
	op.addConfigFileOption("config", "Read further configuration options from this file", "file");
	try{
		op.parseArgs(argc, argv);
	}catch(std::exception& ex){
		std::cerr << ex.what() << std::endl;
		return(1);
	}
	if(op.didPrintUsage())
		return 0;

	if(format==OutputFormat::Unset){
		std::cerr << "An output format must be selected" << std::endl;
		return 1;
	}	
#if ! USE_PHOTOSPLINE
	else if(format==OutputFormat::Spline){
		std::cerr << "Not compiled with spline output support" << std::endl;
		return 1;
	}
#endif

	qcd_init(interaction, pdfName, pdfSubset);
	CrossSectionComputer csc(interaction, ql, false, leptonFlavor);
	csc.setNucleonFractions(protonFraction,1.-protonFraction);
	csc.setEnforceKinematics(enforceKinematics);
	const double Etarget=csc.targetMass();
	
	const double lEStep=(log10(ElMax)-log10(ElMin))/(ElSteps-1);
	const double lEKnotStep=(log10(ElMax)-log10(ElMin))/(ElKnotSteps-1);
	double lxyMin=log10(Q2Min)-log10(csc.CalculateCentreOfMassEnergy(Etarget, ElMax));
	double lxyStep=(0-lxyMin)/(lxySteps-1);
	double lxyKnotStep=(0-lxyMin)/xyKnotSteps;

	std::vector<std::pair<double,double>> totalXS;
	//Total cross section
	if(computeTotal){
		std::cout << "Computing total cross section" << std::endl;
		for(unsigned int eiIdx=0; eiIdx!=ElSteps; eiIdx++){
			double El=ElMin*pow(10.,eiIdx*lEStep);
			double s=csc.CalculateCentreOfMassEnergy(Etarget, El);
			double sig=csc.TotalCrossSection(s, 1);
			totalXS.emplace_back(El,sig);
		}
		if(format==OutputFormat::ASCII){
			std::ofstream out(totalFile);
			for(auto entry : totalXS)
				out << entry.first << ' ' << outputUnits*entry.second << '\n';
		}
#if USE_PHOTOSPLINE
		else if(format==OutputFormat::Spline){
			std::vector<double> eCoords;
			photospline::ndsparse fitData(totalXS.size(),1);
			unsigned int idx=0;
			for(const auto& entry : totalXS){
				eCoords.push_back(log10(entry.first));
				fitData.insertEntry(log10(outputUnits*entry.second),&idx);
				idx++;
			}
			std::vector<double> weights(totalXS.size(),1.);
			std::vector<double> eKnots=generateSequence(-2,(int)(ElKnotSteps+1),log10(ElMin),lEKnotStep);
			photospline::splinetable<> spline;
			spline.fit(fitData,weights,std::vector<std::vector<double>>{eCoords},std::vector<uint32_t>{2},
				{eKnots},{1e-10},std::vector<uint32_t>{2},photospline::splinetable<>::no_monodim,true);
			//keys are SHOUTY because FITS
			spline.write_key("TARGETMASS",Etarget);
			spline.write_key("Q2MIN",Q2Min);
			spline.write_key("UNITS",1/outputUnits);
			spline.write_key("NUTYPE",neutrinoTypeString);
			spline.write_key("CURRENT",to_string(interaction));
			spline.write_fits(totalFile);
		}
#endif
	}

	//Doubly-Differential cross section
	if(computeDouble){
		std::cout << "Computing doubly-differential cross section" << std::endl;
		std::deque<std::pair<double,std::array<unsigned int,3>>> baseData;
		for(unsigned int eIdx=0; eIdx<ElSteps; eIdx++){
			double lE=log10(ElMin)+eIdx*lEStep;
			double El=pow(10.,lE);
			std::cout << "  E=" << El; std::cout.flush();
			double s=csc.CalculateCentreOfMassEnergy(Etarget, El);
			for(unsigned int xIdx=0; xIdx<lxySteps; xIdx++){
				double lx=lxyMin+xIdx*lxyStep;
				double x=pow(10.,lx);
				for(unsigned int yIdx=0; yIdx<lxySteps; yIdx++){
					double ly=lxyMin+yIdx*lxyStep;
					double y=pow(10.,ly);
					double Q2=s*x*y;
					if(Q2<Q2Min) //kinematic limit
						continue;
					double xs=outputUnits*s*x*csc.DoublyDifferentialCrossSection(x,Q2,s);
					if(xs!=0)
						baseData.push_back(std::make_pair(xs,std::array<unsigned int,3>{eIdx,xIdx,yIdx}));
					if(!std::isfinite(xs))
						std::cout << "DING! " << x << ' ' << y << ' ' << xs << '\n';;
					//std::cout << x << ' ' << y << ' ' << xs << '\n';
				}
			}
			std::cout << ' ' << baseData.size() << std::endl;
		}
		std::cout << "Sampled " << baseData.size() << " d2s/dxdy values" << std::endl;
		
		std::vector<double> eCoords, xyCoords;
		eCoords=generateSequence(0,ElSteps,log10(ElMin),lEStep);
		xyCoords=generateSequence(0,lxySteps,lxyMin,lxyStep);
		
		if(format==OutputFormat::ASCII){ //TODO: this mode is largely untested
			std::ofstream out(doubleFile);
			for(const auto& entry : baseData){
				double xs=entry.first;
				auto& coords=entry.second;
				out << pow(10.,eCoords[coords[0]]) // E
				 << '\t' << pow(10.,xyCoords[coords[1]]) //x
				 << '\t' << pow(10.,xyCoords[coords[2]]) //y
				 << '\t' << xs << '\n';
			}
		}
#if USE_PHOTOSPLINE
		if(format==OutputFormat::Spline){
			photospline::ndsparse fitData(baseData.size(),3);
			for(auto& entry : baseData)
				fitData.insertEntry(log10(entry.first),&entry.second[0]);
			std::vector<double> weights(baseData.size(),1.);
			std::vector<double> eKnots, xyKnots;
			eKnots=generateSequence(-2,(int)(ElKnotSteps+1),log10(ElMin),lEKnotStep);
			xyKnots=generateSequence(-2,(int)(xyKnotSteps+2),lxyMin,lxyKnotStep);
			
			photospline::splinetable<> spline;
			spline.fit(fitData,weights,std::vector<std::vector<double>>{eCoords,xyCoords,xyCoords},std::vector<uint32_t>{2,2,2},
					   {eKnots,xyKnots,xyKnots},{1e-10,1e-10,1e-10},std::vector<uint32_t>{2,2,2},
					   photospline::splinetable<>::no_monodim,true);
			
			if(std::isnan(*spline.get_coefficients())){
				std::cerr << "Doubly-differential cross section spline fit failed" << std::endl;
				return 1;
			}

			//keys are SHOUTY because FITS
			spline.write_key("TARGETMASS",Etarget);
			spline.write_key("Q2MIN",Q2Min);
			spline.write_key("UNITS",1/outputUnits);
			spline.write_key("NUTYPE",neutrinoTypeString);
			spline.write_key("CURRENT",to_string(interaction));
			spline.write_fits(doubleFile);
			
			double maxErr=0;
			double avgErr=0;
			decltype(baseData)::value_type maxErrPoint;
			for(auto entry : baseData){
				double xs=entry.first;
				double lcoords[3]={
					log10(ElMin)+entry.second[0]*lEStep,
					lxyMin+entry.second[1]*lxyStep,
					lxyMin+entry.second[2]*lxyStep
				};
				double sv=pow(10.,spline(lcoords));
				assert(sv!=1.0);
				double err=std::abs((xs-sv)/xs);
				if(err>maxErr){
					maxErr=err;
					maxErrPoint=entry;
				}
				avgErr+=err;
			}
			avgErr/=baseData.size();
			std::cout << "Maximum relative error was " << maxErr << std::endl;
			std::cout << " At E=" << pow(10.,log10(ElMin)+maxErrPoint.second[0]*lEStep)
				<< " x=" << pow(10.,lxyMin+maxErrPoint.second[1]*lxyStep)
				<< " y=" << pow(10.,lxyMin+maxErrPoint.second[2]*lxyStep) << std::endl;
			std::cout << " Coordinates: " << maxErrPoint.second[0] << ' ' << maxErrPoint.second[1] << ' ' << maxErrPoint.second[2] << std::endl;
			std::cout << "Average relative error was " << avgErr << std::endl;

			//TODO: possibly revive this code for producing diagnostic plots for spline fits
			/*{
				std::ofstream out("E_worst.txt");
				double lxWorst=lxyMin+maxErrPoint.second[1]*lxyStep;
				double lyWorst=lxyMin+maxErrPoint.second[2]*lxyStep;
				double x=pow(10.,lxWorst);
				double y=pow(10.,lyWorst);
				out << "# x=" << x << " y=" << y << '\n';
				double lcoords[3]={0,lxWorst,lyWorst};
				for(unsigned int eIdx=0; eIdx<2*ElSteps; eIdx++){
					double lE=log10(ElMin)+eIdx*(lEStep/2);
					double El=pow(10.,lE);
					double s=csc.CalculateCentreOfMassEnergy(Etarget, El);
					double Q2=s*x*y;
					if(Q2<Q2Min)
						continue;
					double xs=outputUnits*csc.DoublyDifferentialCrossSection(x,Q2,s);
					lcoords[0]=lE;
					double sxs=pow(10.,spline(lcoords));
					out << El << ' ' << xs << ' ' << sxs << '\n';
				}
			}
			{
				std::ofstream out("x_worst.txt");
				double lEWorst=log10(ElMin)+maxErrPoint.second[0]*lEStep;
				double lyWorst=lxyMin+maxErrPoint.second[2]*lxyStep;
				double El=pow(10.,lEWorst);
				double y=pow(10.,lyWorst);
				out << "# E=" << El << " y=" << y << '\n';
				double lcoords[3]={lEWorst,0,lyWorst};
				double s=csc.CalculateCentreOfMassEnergy(Etarget, El);
				for(unsigned int xIdx=0; xIdx<2*lxySteps; xIdx++){
					double lx=lxyMin+xIdx*(lxyStep/2);
					double x=pow(10.,lx);
					double Q2=s*x*y;
					if(Q2<Q2Min)
						continue;
					double xs=outputUnits*csc.DoublyDifferentialCrossSection(x,Q2,s);
					lcoords[1]=lx;
					double sxs=pow(10.,spline(lcoords));
					out << x << ' ' << xs << ' ' << sxs << '\n';
				}
			}
			{
				std::ofstream out("y_worst.txt");
				double lEWorst=log10(ElMin)+maxErrPoint.second[0]*lEStep;
				double lxWorst=lxyMin+maxErrPoint.second[1]*lxyStep;
				double El=pow(10.,lEWorst);
				double x=pow(10.,lxWorst);
				out << "# E=" << El << " x=" << x << '\n';
				double lcoords[3]={lEWorst,lxWorst,0};
				double s=csc.CalculateCentreOfMassEnergy(Etarget, El);
				for(unsigned int yIdx=0; yIdx<4*lxySteps; yIdx++){
					double ly=lxyMin+yIdx*(lxyStep/4);
					double y=pow(10.,ly);
					double Q2=s*x*y;
					if(Q2<Q2Min)
						continue;
					double xs=outputUnits*csc.DoublyDifferentialCrossSection(x,Q2,s);
					lcoords[2]=ly;
					double sxs=pow(10.,spline(lcoords));
					out << y << ' ' << xs << ' ' << sxs << '\n';
				}
			}*/

#if USE_CUBPACKPP
			//Total cross section check
			//Ensure that integrating the fitted spline produces a result similar to the
			//originally inegrated total cross section
			if(computeTotal){
				//TODO: this filename should be configurable
				std::ofstream out("Total_XS_s.dat");
				auto preciseIt=totalXS.begin();
				double factor=pow(10.,lEStep);
				for(unsigned int eiIdx=0; eiIdx!=ElSteps; eiIdx++){
					double El=ElMin*pow(10.,eiIdx*lEStep);
					double s=csc.CalculateCentreOfMassEnergy(Etarget, El);
					double cMin=log10(Q2Min)-log10(s);
					double lE=log10(El);
					TRIANGLE xyDomain(Point(cMin,0),Point(0,cMin),Point(0,0));
					auto integrand=[=,&spline](double lx, double ly)->double{
						double x=pow(10.,lx);
						double y=pow(10.,ly);
						//double jacobian=s*x*x*y*log(10.)*log(10.);
						double jacobian=x*y*log(10.)*log(10.);
						double lcoords[3]={lE,lx,ly};
						double result=jacobian*pow(10.,spline(lcoords));
						return result;
					};
					double sig=Integrate(integrand,xyDomain,0,1e-3,1e7);
					assert(std::abs(El-preciseIt->first)/El<1e-5);
					double relerr=std::abs(sig-preciseIt->second)/preciseIt->second;
					out << El << ' ' << sig << ' ' << preciseIt->first << ' ' << preciseIt->second << ' ' << relerr << std::endl;
					preciseIt++;
				}
			}
#endif //USE_CUBPACKPP
		}
#endif //USE_PHOTOSPLINE
	} //computeDouble

	if(computeSingle){ //singly-differential cross section
		std::cout << "Computing singly-differential cross section" << std::endl;
		std::deque<std::pair<double,std::array<unsigned int,2>>> baseData;
		std::vector<double> lenergies=generateSequence(0,ElSteps-1,log10(ElMin),lEStep);
		
		std::ofstream out;
		if(format==OutputFormat::ASCII)
			out.open(singleFile);
		else if(format==OutputFormat::Spline)
			out.open(singleFile+".raw");

        auto preciseIt=totalXS.begin();
        
		for(unsigned int eiIdx=0; eiIdx!=ElSteps; eiIdx++){
			double El=pow(10.,lenergies[eiIdx]);
			double s=csc.CalculateCentreOfMassEnergy(Etarget, El);

			auto computeEout=[El,ElMin,ElSteps,squareSingleTable,&lenergies](unsigned int outIdx){
				if(!squareSingleTable)
					return pow(10.,lenergies[outIdx]);
				double z=outIdx/(ElSteps-1.);
				return z*(El - ElMin)+ElMin;
			};

            double totalEstimate=0;
            
			for(unsigned int eoIdx=0; eoIdx!=ElSteps; eoIdx++){
				double Eout=computeEout(eoIdx);
				if(Eout>El || eiIdx==0){
					if(format==OutputFormat::ASCII)
						out << El << ' ' << Eout << ' ' << 0 << '\n';
					//don't break; this isn't terribly efficient, but helps keep the ASCII output uniform
					continue;
				}
				double y=1-(Eout/El);
				double dsdy=0;
				
				//At high energies the cross section has a peak and then drops very rapidly to zero as y->0+
				//At some point this peak becomes sharp enough to be hidden between tabulated dsdE values,
				//causing interpolation schemes to badly underestimate the cross section in this region.
				//Since the peak tends to represent a more important contribution than the drop to zero
				//(which will anyway be enforced once Eout>El) we cheat and replace the value with one 
				//computed from the average since the previous Eout to push interpolations in the right direction.
				if(singleDiffPeakHack && Eout==El){
					double Eprev=computeEout(eoIdx-1);
					double yPrev=1-(Eprev/El);
					double dsdyPrev=csc.SinglyDifferentialCrossSection(s, yPrev, Q2Min);
					
					double dsdyAvg=csc.AverageSinglyDifferentialCrossSection(s, 0, yPrev, Q2Min);
					//want to report dsdy such that dsdyAvg=(dsdyPrev+dsdy)/2
					std::cout << ' ' << El << "->" << Eout << " computed dsdyAvg=" << outputUnits*dsdyAvg << '\n';
					dsdy=2*dsdyAvg-dsdyPrev;
				}
				else
					dsdy=csc.SinglyDifferentialCrossSection(s, y, Q2Min);
				
				//accumulate estimate of the total cross section
				{
					double dy=0;
					if(eoIdx)
					dy=(1-(computeEout(eoIdx-1)/El))-y;
					totalEstimate+=(Eout<El?2:1)*dsdy*dy;
				}
                
                dsdy*=outputUnits;
                if(format==OutputFormat::ASCII){
					//out << El << ' ' << Eout << ' ' << dsdy << '\n';
					out << El << ' ' << Eout << ' ' << dsdy;
					/*if(eoIdx){
						double dsdyAvg=csc.AverageSinglyDifferentialCrossSection(s, y, 1-(computeEout(eoIdx-1)/El), Q2Min);
						out << ' ' << outputUnits*dsdyAvg;
					}*/
					out << '\n';
                }
				else if(format==OutputFormat::Spline){
					baseData.push_back(std::make_pair(dsdy,std::array<unsigned int,2>{eiIdx,eoIdx}));
					out.write((const char*)&dsdy,sizeof(dsdy));
					out.write((const char*)&eiIdx,sizeof(eiIdx));
					out.write((const char*)&eoIdx,sizeof(eoIdx));
				}
			}
            
			//Total cross section check
			if(computeTotal){
				std::cout << "Check: E1=" << El << ": " << outputUnits*preciseIt->second << ' ' << outputUnits*totalEstimate/2 << '\n';
				preciseIt++;
			}
		}
		out.close();
#if USE_PHOTOSPLINE
		if(format==OutputFormat::Spline){
			{
				std::cout << "Attempting to read " << singleFile+".raw" << std::endl;
				std::ifstream raw(singleFile+".raw");
				double dsdy;
				unsigned int eiIdx,eoIdx;
				if(!raw){
					std::cerr << "Failed to open single.raw" << std::endl;
					return 1;
				}
				while(raw.read((char*)&dsdy,sizeof(dsdy))){
					raw.read((char*)&eiIdx,sizeof(eiIdx));
					raw.read((char*)&eoIdx,sizeof(eoIdx));
					baseData.push_back(std::make_pair(dsdy,std::array<unsigned int,2>{eiIdx,eoIdx}));
				}
				std::cout << raw.good() << ' ' << raw.bad() << ' ' << raw.fail() << ' ' << raw.eof() << std::endl;
			}
			std::cout << "Sampled " << baseData.size() << " dsdy values" << std::endl;
			//TODO correct below to say that it is dsdy, or add factors to make it dsdE
			photospline::ndsparse fitData(baseData.size(),2);
			for(auto entry : baseData){
				if(entry.first<=0 || !std::isfinite(log10(entry.first))){
					//std::cout << "Bad dsdy value: " << entry.first << " at (" << entry.second[0] << ',' << entry.second[1] << ")\n";
					//if(entry.first<=0)
					//	fitData.insertEntry(-50,&entry.second[0]);
					continue;
				}
				fitData.insertEntry(log10(entry.first),&entry.second[0]);
			}
			std::vector<double> weights(baseData.size(),1.);
            
            std::vector<double> loenergies;
            if(squareSingleTable)
                loenergies=generateSequence(0,ElSteps-1,0,1./(ElSteps-1));
            else
                loenergies=lenergies;
            
			std::vector<double> eKnots;
			eKnots=generateSequence(-2,(int)(ElKnotSteps+1),log10(ElMin),lEKnotStep);
            
            std::vector<double> eoKnots;
			if(squareSingleTable){
				//linear knots
                //eoKnots=generateSequence(-2,(int)(zKnotSteps+1),0,1./(zKnotSteps-1));
				for(unsigned int i=0; i<zKnotSteps; i++){
					double z=1-pow(i/(zKnotSteps-1.)-1,2.);
					eoKnots.push_back(z);
				}
				padKnots(eoKnots,2);
			}
            else
                eoKnots=eKnots;
            
            std::cout << "Energy abcissas:";
            for(auto len : lenergies)
                std::cout << ' ' << len;
            std::cout << std::endl;
            
            if(squareSingleTable){
                std::cout << "z abcissas:";
                for(auto z : loenergies)
                    std::cout << ' ' << z;
                std::cout << std::endl;
            }
            
			std::cout << "Energy knots:";
			for(auto eknot : eKnots)
				std::cout << ' ' << eknot;
			std::cout << std::endl;
            
            if(squareSingleTable){
                std::cout << "z knots:";
                for(auto zknot : eoKnots)
                    std::cout << ' ' << zknot;
                std::cout << std::endl;
				if(writeSplineSlices){
					std::ofstream knotfile("single_differential_slices/zknots.dat");
					for(auto zknot : eoKnots)
						knotfile << zknot << '\n';
				}
            }
            
			photospline::splinetable<> spline;
			spline.fit(fitData,weights,std::vector<std::vector<double>>{lenergies,loenergies},
			           std::vector<uint32_t>{2,2},
			           {eKnots,eoKnots},{splineRegSingle*20,splineRegSingle},std::vector<uint32_t>{2,2},
			           photospline::splinetable<>::no_monodim,true);
			if(std::isnan(*spline.get_coefficients())){
				std::cerr << "Singly-differential cross section spline fit failed" << std::endl;
				return 1;
			}
			//keys are SHOUTY because FITS
			spline.write_key("TARGETMASS",Etarget);
			spline.write_key("Q2MIN",Q2Min);
			spline.write_key("UNITS",1/outputUnits);
			spline.write_key("NUTYPE",neutrinoTypeString);
			spline.write_key("CURRENT",to_string(interaction));
			spline.write_fits(singleFile);

			double maxErr=0;
			double avgErr=0;
			std::ofstream sliceFile;
			unsigned int lastEiIdx=0;
			auto write_dense_spline = [&]()->void{
				std::ofstream denseSpline("single_differential_slices/"+to_string(lastEiIdx)+"_spline.dat");
				double ein=pow(10.,log10(ElMin)+lastEiIdx*lEStep);
				denseSpline << "# E_in = " << ein << '\n';
				denseSpline << "# E_out z dsdE_spline\n";
				unsigned int zSteps=5*(ElSteps-1)+1;
				double zStep=1./(zSteps-1);
				for(unsigned int zIdx=0; zIdx<zSteps; zIdx++){
					double z=zIdx*zStep;
					double eout=z*(ein-ElMin)+ElMin;
					double lcoords[2]={
						lenergies[lastEiIdx],
						(squareSingleTable?z:log10(eout)),
					};
					double sv=pow(10.,spline(lcoords));
					denseSpline << eout << '\t' << z << '\t' << sv << '\n';
				}
			};
			decltype(baseData)::value_type maxErrPoint;
			for(auto entry : baseData){
				double dsdE=entry.first;
				if(!dsdE)
					continue;
                
				double lcoords[2]={
					lenergies[entry.second[0]],
					loenergies[entry.second[1]]
				};
                    
				double sv=pow(10.,spline(lcoords));
				//assert(sv!=1.0);
				double err=std::abs((dsdE-sv)/dsdE);
				double ein=pow(10.,log10(ElMin)+entry.second[0]*lEStep);
				double eout=(squareSingleTable?lcoords[1]*(ein-ElMin)+ElMin:pow(10.,lcoords[1]));
				if(writeSplineSlices){
					if(!sliceFile.is_open() || entry.second[0]!=lastEiIdx){
						if(sliceFile.is_open()){
							sliceFile.close();
							write_dense_spline();
						}
						sliceFile.open("single_differential_slices/"+to_string(entry.second[0])+".dat");
						sliceFile << "# E_in = " << ein << '\n';
						sliceFile << "# E_out z dsdE_computed dsdE_spline relative_error\n";
					}
					sliceFile << eout << '\t' << (squareSingleTable?lcoords[1]:(eout-ElMin)/(ein-ElMin))
						<< '\t' << dsdE << '\t' << sv << '\t' << err << '\n';
					lastEiIdx=entry.second[0];
				}
				if(err>maxErr){
					std::cout << "New worst error: " << err << " for " << dsdE << " vs. " << sv
					<< " at El=" << ein << " Eout=" << eout << '\n';
					maxErr=err;
					maxErrPoint=entry;
				}
				else if(err>.01){
					std::cout << "Large error: " << err << " for " << dsdE << " vs. " << sv 
					<< " at El=" << ein << " Eout=" << eout << '\n';
                    //std::cout << "  logical coordinates: " << entry.second[0] << ',' << entry.second[1] << '\n';
                    //std::cout << "  physical coordinates: " << lcoords[0] << ',' << lcoords[1] << '\n';
				}
				avgErr+=err;
			}
			if(sliceFile.is_open()){
				sliceFile.close();
				write_dense_spline();
			}
			avgErr/=baseData.size();
			{
				double ein=pow(10.,log10(ElMin)+maxErrPoint.second[0]*lEStep);
				double eout=(squareSingleTable?loenergies[maxErrPoint.second[1]]*(ein-ElMin)+ElMin:pow(10.,loenergies[maxErrPoint.second[1]]));
				std::cout << "Maximum relative error was " << maxErr << std::endl;
				std::cout << " At El=" << pow(10.,log10(ElMin)+maxErrPoint.second[0]*lEStep)
				          << " Eout=" << eout << std::endl;
				std::cout << " Coordinates: " << maxErrPoint.second[0] << ' ' << maxErrPoint.second[1] << std::endl;
				std::cout << "Average relative error was " << avgErr << std::endl;
			}
		}
#endif //USE_PHOTOSPLINE
	}
}
