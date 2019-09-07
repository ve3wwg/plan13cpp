//////////////////////////////////////////////////////////////////////
// Plan13.hpp -- Borrowed (see README) and cleaned up to eliminate
// Arduino code and made Linux friendly, by Warren Gay VE3WWG
///////////////////////////////////////////////////////////////////////

#ifndef PLAN13_HPP
#define PLAN13_HPP

#include <time.h>

class Plan13 {
	char		name[10];

	int 		aYear;			// Provided by setTime()
	int 		aMonth;
	int 		aMday;
	int 		aHour;
	int 		aMin;
	int 		aSec;

	double		EL;			// Elevaton
	double		AZ;			// Azimuth
	double		SLON;			// Lon, + East
	double		SLAT;			// Lat, + North
	double		RR;			// Range rate, km/s

	double 		rx, tx;
	double 		observer_lon;
	double 		observer_lat;
	int 		observer_height;

	unsigned long 	rxOutLong;
	unsigned long 	txOutLong;

	unsigned long 	rxFrequencyLong;
	unsigned long 	txFrequencyLong;
	double 		dopplerFactor;

	double		TN;                           /*                                    */

	double		E;
	double		N;

	double		CL;
	double		CS;
	double		SL;
	double		CO;
	double		SO;
	double		RE;
	double		FL;
	double		RP;
	double		XX;
	double		ZZ;
	double		D;
	double		R;
	double		Rx;
	double		Ry;
	double		Rz;
	double		Ex;
	double		Ey;
	double		Ez;
	double		Ny;
	double		Nx;
	double		Nz;
	double		Ox;		// Observer's XYZ coords at Earth's surface
	double		Oy;
	double		Oz;
	double		U;
	double		Ux;		// Observer's unit vectors UP EAST and NORTH in GEOCENTRIC coords.
	double		Uy;
	double		Uz;
	double		WW;		// Earth's rotation rate, rads/whole day
	double		WE;		//       ditto            radians/day
	double		W0;		//       ditto            radians/day
	double		VOx;		// Observer's velocity, GEOCENTRIC coords. (VOz=0)
	double		VOy;
	double		VOz;
	double		DE;		// Convert satellite Epoch to Day No. and Fraction of day
	double		GM;		// Earth's gravitational constant
	double		J2;		// 2nd Zonal coeff, earth's gravity field
	double		N0;		// Mean motion (rad/s)
	double		A0;		// Semi major axis (km)
	double		b0;		// Semi minor axis (km)
	double		SI;		// sin(IN)
	double		CI;		// cos(IN)
	double		PC;		// Precession const, rad/Day
	double		QD;		// Node precession rate, rad/day
	double		WD;		// Perigee precession rate, rad/day
	double		DC;		// Drag coeff. (Angular momentum rate)/(Ang mom)  s^-1
	double		YG;
	double		G0;
	double		MAS0;
	double		MASD;
	double		INS;
	double		CNS;
	double		SNS;		// Sun's Equation of centre terms
	double		EQC1;
	double		EQC2;
	double		TEG;		// Elapsed Time: Epoch - YG
	double		GHAE;		// GHA Aries, epoch
	double		MRSE;		// Mean RA Sun at Sat epoch
	double		MASE;		// Mean MA Sun  ..
	double		ax;
	double		ay;
	double		az;
	int		OLDRN;

	double		T;
	double		DT;
	double		KD;
	double		KDP;
	double		M;
	int		DR;
	long		RN;
	double		EA;
	double		C;
	double		S;
	double		DNOM;
	double		A;
	double		B;
	double		RS;
	double		Sx;
	double		Sy;
	//double   Sz;
	double		Vx;
	double		Vy;
	double		Vz;
	double		AP;
	double		CWw;
	double		SW;
	double		RAAN;
	double		CQ;
	double		SQ;
	double		CXx;
	double		CXy;
	double		CXz;
	double		CYx;
	double		CYy;
	double		CYz;
	double		CZx;
	double		CZy;
	double		CZz;
	double		SATx;
	double		SATy;
	double		SATz;
	double		ANTx;
	double		ANTy;
	double		ANTz;
	double		VELx;
	double		VELy;
	double		VELz;
	double		Ax;
	double		Ay;
	double		Az;
	double		Sz;
	double		GHAA;

	double		DS;
	double		DF;
	
	// Keplerians

	char		SAT[20];
	long		SATNO;
	double		YE;		// Epoch year
	double		TE;		// Epoch time (days)
	double		IN;		// Inclination (deg)
	double		RA;		// R.A.A.N (deg)
	double		EC;		// Eccentricity
	double		WP;		// Arg perifee (deg)
	double		MA;		// Mean anomaly (rev/d)
	double		MM;		// Mean motion (rev/d)
	double		M2;		// Decay rate (rev/d/d)
	long		RV;		// Orbit number
	double		ALON;		// Sat attitude (deg) 180 = nominal
	double		ALAT;		// Sat attitude (deg) 0 = nominal
	double		rxOut;		
	double		txOut;

	// Location

	char		LOC[20];
	double		LA;
	double		LO;
	double		HT;

	double		HR; // Hours
	double		DN;

	void		*arg = nullptr;
	void		(*write_cb)(const void *buf,unsigned bytes,void *arg) = nullptr;

	void satvec(void);
	double rad(double deg);
	double deg(double rad);
	double FNatn(double y, double x);
	double FNday(int year, int month, int day);
	double myFNday(int year, int month, int day, int uh, int um, int us);
	void rangevec(void);

	double get_element(const char *gstr,int gstart,int gstop);
	static void write(const void *buf,unsigned bytes,void *arg);
public:
	//////////////////////////////////////////////////////////////
	// Defaults to writing to stdout when wrfunc is supplied as
	// nullptr.
	//////////////////////////////////////////////////////////////

	Plan13(
	  void (*wrfunc)(const void *buf,unsigned bytes,void *arg)=nullptr,
	  void *arg=nullptr
	);

	void calculate(void);
	void print_header();
	void print();

	void init(void);
	void set_frequency(unsigned long rx_frequency,unsigned long tx_frequency);
	void set_location(double lon,double lat,int height);
	void set_time(time_t time_date);
	void set_time(int yearIn,int monthIn,int mDayIn,int hourIn,int minIn,int secIn);

	void set_elements(
		double YE_in,
		double TE_in,
		double IN_in,
		double RA_in,
		double EC_in,
		double WP_in,
		double MA_in,
		double MM_in, 
		double M2_in,
		double RV_in,
		double ALON_in
	);
	bool load_elements(const char *kep0,const char *kep1,const char *kep2);

	int get_doppler(unsigned long freq);
	int get_doppler64(unsigned long freq);

	double get_EL() { return EL; }
	double get_AZ() { return AZ; }
	double get_LON() { return SLON; }
	double get_LAT() { return SLAT; }
	double get_RR() { return RR; }

	double get_RE() { return RE; }
	double get_RS() { return RS; }
};

#endif // PLAN13_HPP

// End plan13.hpp

