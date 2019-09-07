//////////////////////////////////////////////////////////////////////
// plan13.cpp -- Borrowed from https://github.com/BackupGGCode/qrptracker.git
///////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>

#include "plan13.hpp"

#define DEBUG 0
#define ONEPPM 1.0e-6

#include <stdio.h>
#include <math.h>
#include <stdint.h>

static const double YM = 365.25;	// Days in a year
static const double YT = 365.2421970;	// Tropical year, days

//////////////////////////////////////////////////////////////////////
// Conversions degrees / radians
//////////////////////////////////////////////////////////////////////

double
Plan13::rad(double deg) {
	return M_PI / 180.0 * deg;
}

double
Plan13::deg(double rad) {
	return rad * 180.0 / M_PI;
}

double
Plan13::FNatn(double y,double x) {
	double a;

	if ( x != 0.0 )
		a = atan(y / x);
	else	a = M_PI / 2.0 * sin(y);

	if ( x < 0.0 ) {
		a = a + M_PI;
	}
	if ( a < 0.0 ) {
		a = a + 2.0 * M_PI;
	}
	return a;
}

//////////////////////////////////////////////////////////////////////
// Convert date to day number
//
// Function returns a general day number from year, month, and day.
// Value is (JulianDay - 1721409.5) or (AmsatDay + 722100)
//////////////////////////////////////////////////////////////////////

double
Plan13::FNday(int year,int month,int day) {
	double JulianDate;

#if DEBUG
	printf("## FNDay:  ");
	printf("Year: %d",year);
	printf(" Month: %d",month);
	printf(" Day: %d",day);
#endif
	if ( month <= 2 ) {
		year -= 1;
		month += 12;
	}

	JulianDate = (long)(year * YM) + (int)((month + 1) * 30.6) + (day - 428);

#if DEBUG
	printf(" JD: %lf\n",JulianDate);
#endif
	return JulianDate;
}

Plan13::Plan13(void (*wrfunc)(const void *buf,unsigned bytes,void *arg),void *arg) {
	if ( !wrfunc )
		write_cb = Plan13::write;	// Writs to stdout
	else	write_cb = wrfunc;		// User supplied write()
	this->arg = arg;			// User argument
}

void
Plan13::init() {

	// Observer's location
#if DEBUG
	printf("Start initSat()\n");
#endif
	LA = rad(observer_lat);
	LO = rad(observer_lon);
	HT = ((double) observer_height)/1000.0; // this needs to be in km
	CL = cos(LA);
	SL = sin(LA);
	CO = cos(LO);
	SO = sin(LO);

	// WGS-84 Earth Ellipsoid
	RE = 6378.137;
	FL = 1.0 / 298.257224;
	
	// IAU-76 Earth Ellipsoid
	// RE = 6378.140;
	// FL = 1.0 / 298.257;
	
	RP = RE * (1.0 - FL);
	XX = RE * RE;
	ZZ = RP * RP;
	
	D = sqrt(XX * CL * CL + ZZ * SL * SL);
	
	Rx = XX / D + HT;
	Rz = ZZ / D + HT;
	
	// Observer's unit vectors Up EAST and NORTH in geocentric coordinates
	Ux = CL * CO;
	Ex = -SO;
	Nx = -SL * CO;
	
	Uy = CL * SO;
	Ey = CO;
	Ny = -SL * SO;
	
	Uz = SL;
	Ez = 0;
	Nz = CL;
	
	// Observer's XYZ coordinates at earth's surface
	Ox = Rx * Ux;
	Oy = Rx * Uy;
	Oz = Rz * Uz;
	
	// Convert angles to radians, etc.
	RA = rad(RA);
	IN = rad(IN);
	WP = rad(WP);
	MA = rad(MA);
	MM = MM * 2.0 * M_PI;
	M2 = M2 * 2.0 * M_PI;
	
	// YM = 365.25;		Mean year, days
	// YT = 365.2421970;	Tropical year, days
	WW = 2.0 * M_PI / YT;	// Earth's rotation rate, rads/whole day 
	WE = 2.0 * M_PI + WW;	// Earth's rotation rate, rads/day
	W0 = WE / 86400;	// Earth's rotation rate, rads/sec
	
	// Observer's velocity, geocentric coordinates
	VOx = -Oy * W0;
	VOy = Ox * W0;
	
	// Convert satellite epoch to Day No. and fraction of a day
	DE = FNday(YE, 1, 0) + (int)TE;
	
	TE = TE - (int)TE;
#if DEBUG
	printf("DE: %f TE: %lf\n",DE,TE);
#endif
	// Average Precession rates
	GM = 3.986E5;				// Earth's gravitational constant km^3/s^2
	J2 = 1.08263E-3;			// 2nd Zonal coeff, Earth's gravity Field
	N0 = MM / 86400.0;			// Mean motion rads/s
	A0 = pow(GM / N0 / N0, 1.0 / 3.0);	// Semi major axis km
	b0 = A0 * sqrt(1.0 - EC * EC);		// Semi minor axis km
	SI = sin(IN);
	CI = cos(IN);
	PC = RE * A0 / (b0 * b0);
	PC = 1.5 * J2 * PC * PC * MM;		// Precession const, rad/day
	QD = -PC * CI;				// Node Precession rate, rad/day
	WD = PC *(5.0 * CI*CI - 1.0) / 2.0;	// Perigee Precession rate, rad/day
	DC = -2.0 * M2 / MM / 3.0;		// Drag coeff
	
	// Sideral and solar data. Never needs changing. Valid to year 2000+
	
	// GHAA, Year YG, Jan 0.0 
	YG = 2010;
	G0 = 99.5578;
	// MA Sun and rate, deg, deg/day
	MAS0 = 356.4485;
	MASD = 0.98560028;
	// Sun's inclination
	INS = rad(23.4380);
	CNS = cos(INS);
	SNS = sin(INS);
	// Sun's equation of center terms
	EQC1 = 0.03341;
	EQC2 = 0.00035;
	
	// Bring Sun data to satellite epoch
	TEG = (DE - FNday(YG, 1, 0)) + TE;	// Elapsed Time: Epoch - YG
	GHAE = rad(G0) + TEG * WE;		// GHA Aries, epoch
	MRSE = rad(G0) + (TEG * WW) + M_PI;	// Mean RA Sun at Sat Epoch
	MASE = rad(MAS0 + MASD * TEG);		// Mean MA Sun
	
	// Antenna unit vector in orbit plane coordinates
	CO = cos(rad(ALON));
	SO = sin(rad(ALON));
	CL = cos(rad(ALAT));
	SL = sin(rad(ALAT));
	ax = -CL * CO;
	ay = -CL * SO;
	az = -SL;
	
	// Miscellaneous
	OLDRN = -99999;
}

//////////////////////////////////////////////////////////////////////
// Calculate satellite position at DN, TN
//////////////////////////////////////////////////////////////////////

void
Plan13::satvec() {

	T = (DN - DE) + (TN - TE);//83.848;	// Elapsed T since epoch
#if DEBUG
	printf("T: %f\n",T);
#endif
	DT = DC * T / 2.0;			// Linear drag terms
	KD = 1.0 + 4.0 * DT;
	KDP = 1.0 - 7.0 * DT;
	M = MA + MM * T * (1.0 - 3.0 * DT); 	// Mean anomaly at YR,/ TN
	DR = (int)(M / (2.0 * M_PI));		// Strip out whole no of revs
	M = M - DR * 2.0 * M_PI;              	// M now in range 0 - 2PI
	RN = RV + DR + 1;                   	// Current orbit number
	
	// Solve M = EA - EC * sin(EA) for EA given M, by Newton's method 
	EA = M;					// Initail solution
	do	{
		C = cos(EA);
		S = sin(EA);
		DNOM = 1.0 - EC * C;
		D = (EA - EC * S - M) / DNOM;	// Change EA to better resolution 
		EA = EA - D;			// by this amount until converged 
	} while ( fabs(D) > 1.0E-5 );
	
	// Distances
	A = A0 * KD;
	B = b0 * KD;
	RS = A * DNOM;
	
	// Calculate satellite position and velocity in plane of ellipse
	Sx = A * (C - EC);
	Vx = -A * S / DNOM * N0;
	Sy = B * S;
	Vy = B * C / DNOM * N0;
	
	AP = WP + WD * T * KDP;
	CWw = cos(AP);
	SW = sin(AP);
	RAAN = RA + QD * T * KDP;
	CQ = cos(RAAN);
	SQ = sin(RAAN);
	
	// Plane -> celestial coordinate transformation, [C] = [RAAN]*[IN]*[AP]
	CXx = CWw * CQ - SW * CI * SQ;
	CXy = -SW * CQ - CWw * CI * SQ;
	CXz = SI * SQ;
	CYx = CWw * SQ + SW * CI * CQ;
	CYy = -SW * SQ + CWw * CI * CQ;
	CYz = -SI * CQ;
	CZx = SW * SI;
	CZy = CWw * SI;
	CZz = CI;
	
	// Compute satellite's position vector, ANTenna axis unit vector 
	// and velocity  in celestial coordinates. (Note: Sz = 0, Vz = 0)
	SATx = Sx * CXx + Sy * CXy;
	ANTx = ax * CXx + ay * CXy + az * CXz;
	VELx = Vx * CXx + Vy * CXy;
	SATy = Sx * CYx + Sy * CYy;
	ANTy = ax * CYx + ay * CYy + az * CYz;
	VELy = Vx * CYx + Vy * CYy;
	SATz = Sx * CZx + Sy * CZy;
	ANTz = ax * CZx + ay * CZy + az * CZz;
	VELz = Vx * CZx + Vy * CZy;
	
	// Also express SAT, ANT, and VEL in geocentric coordinates
	GHAA = GHAE + WE * T;		// GHA Aries at elaprsed time T
	C = cos(-GHAA);
	S = sin(-GHAA);
	Sx = SATx * C - SATy * S;
	Ax = ANTx * C - ANTy * S;
	Vx = VELx * C - VELy * S;
	Sy = SATx * S + SATy * C;
	Ay = ANTx * S + ANTy * C;
	Vy = VELx * S + VELy * C;
	Sz = SATz;
	Az = ANTz;
	Vz = VELz;
}

//////////////////////////////////////////////////////////////////////
// Compute and manipulate range/velocity/antenna vectors
//////////////////////////////////////////////////////////////////////

void
Plan13::rangevec() {

	// Range vector = sat vector - observer vector
	Rx = Sx - Ox;
	Ry = Sy - Oy;
	Rz = Sz - Oz;

	R = sqrt(Rx * Rx + Ry * Ry + Rz * Rz);    /* Range Magnitute */

	// Normalize range vector
	Rx = Rx / R;
	Ry = Ry / R;
	Rz = Rz / R;
	U = Rx * Ux + Ry * Uy + Rz * Uz;
	E = Rx * Ex + Ry * Ey;
	N = Rx * Nx + Ry * Ny + Rz * Nz;

	AZ = deg(FNatn(E, N));
	EL = deg(asin(U));
	
	// Solve antenna vector along unit range vector, -r.a = cos(SQ)
	// SQ = deg(acos(-(Ax * Rx + Ay * Ry + Az * Rz)));

	// Calculate sub-satellite Lat/Lon 
	SLON = deg(FNatn(Sy, Sx));		// Lon, + East
	SLAT = deg(asin(Sz / RS));		// Lat, + North
	
	if ( SLON > 180.0 )
		SLON -= 360.0;			// -ve is degrees West

	// Resolve Sat-Obs velocity vector along unit range vector. (VOz = 0)
	RR = (Vx - VOx) * Rx + (Vy - VOy) * Ry + Vz * Rz; // Range rate, km/sec
	//FR = rxFrequency * (1 - RR / 299792);
	dopplerFactor = RR / 299792.0;
	int rxDoppler = get_doppler(rxFrequencyLong);
	int txDoppler = get_doppler(txFrequencyLong);
	rxOutLong = rxFrequencyLong - rxDoppler;
	txOutLong = txFrequencyLong + txDoppler;
}

int
Plan13::get_doppler(unsigned long freq) {

	freq = (freq + 50000L) / 100000L;
	long factor = dopplerFactor * 1E11;
	int digit;
	double tally = 0.0;
	for (int x = 4; x > -1; x--) {
		digit = freq/pow(10,x);
		long bare = digit * pow(10,x);
		freq = freq - bare;
		double inBetween =  factor * (double(bare) / 1E6);
		tally += inBetween;
	}
	return int( tally + .5); //round
}

int
Plan13::get_doppler64(unsigned long freq) {

	//UNUSED: long factor = dopplerFactor * 1E11;
	uint64_t doppler_sixfour = freq * dopplerFactor;
	return (int) doppler_sixfour/1E11;
}

//////////////////////////////////////////////////////////////////////
// Setter method for uplink (tx) and downlink (rx) frequencies in Hz.
// These data need not be set if the doppler shifted frequencies are not needed.
//////////////////////////////////////////////////////////////////////

void
Plan13::set_frequency(
  unsigned long rxFrequency_in,		// downlink freq
  unsigned long txFrequency_in		// uplink freq
) {

	rxFrequencyLong = rxFrequency_in;
	txFrequencyLong = txFrequency_in;
}

//////////////////////////////////////////////////////////////////////
// Setter method for indicating the location of the ground station. 
// This and setElements() must be done before the calculate method is applied for the first time.
// Thereafter, however, it doesn't matter unless the groundstation changes position.
// 
// \param observer_lon_in the longitude of the observer, with east positive and west negative
// \param observer_lat_in the latitude of the observer with north positive and south negative
// \param height height of the observer, in meters
//////////////////////////////////////////////////////////////////////

void
Plan13::set_location(double observer_lon_in,double observer_lat_in,int height) {

	observer_lon = observer_lon_in;//-64.375; //0.06; // lon east is positive, west is negative
	observer_lat = observer_lat_in;//45.8958; //52.21; //Cambridge UK
	observer_height = height; //60m height in meters
}

//////////////////////////////////////////////////////////////////////
// Set date/time using Unix epoch time:
//////////////////////////////////////////////////////////////////////

void
Plan13::set_time(time_t time_date) {
	struct tm td;

	::gmtime_r(&time_date,&td);
	set_time(td.tm_year+1900,td.tm_mon+1,td.tm_mday,td.tm_hour,td.tm_min,td.tm_sec);
}

//////////////////////////////////////////////////////////////////////
// Setter method for UTC time at which the satellite is to be observed. 
// This is usually the current tim.
//////////////////////////////////////////////////////////////////////

void
Plan13::set_time(int yearIn,int monthIn,int mDayIn,int hourIn,int minIn,int secIn) {

	aYear = yearIn;
	aMonth = monthIn;
	aMday = mDayIn;
	aHour = hourIn;
	aMin  = minIn;
	aSec  = secIn;

	DN = FNday(aYear,aMonth,aMday);
	TN = ((double)aHour + ((double)aMin + ((double)aSec/60.0)) /60.0)/24.0;
	DN = (long)DN;
}

//////////////////////////////////////////////////////////////////////
// Sets the keplerian elements for the following calculations. 
//////////////////////////////////////////////////////////////////////

void
Plan13::set_elements(
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
) {
	
	YE = YE_in;
	TE = TE_in;
	IN = IN_in;
	RA = RA_in;
	EC = EC_in;
	WP = WP_in;
	MA = MA_in;
	MM = MM_in;
	M2 = M2_in;
	RV = RV_in;
	ALON = ALON_in;
}

//////////////////////////////////////////////////////////////////////
// A function that joins together the necessary functions for calculating 
// the satellite position. You must set the keplerian elements, time and observer lat/long
// before using this.
//////////////////////////////////////////////////////////////////////

void
Plan13::calculate() {
	init();
	satvec();
	rangevec();
}

void
Plan13::print_header() {
	write_cb("  Date / Time UTC    Sat   Azimuth   Elevation   Latitude   Longitude   RR\n",75,arg);
	write_cb("------------------- ----- ---------- ---------- ---------- ---------- ------\n",77,arg);
}

void
Plan13::print() {
	char buf[120];
	int n;

	n = snprintf(buf,sizeof buf,"%04d-%02d-%02d %02d:%02d:%02d "
		"%-5.5s %10.6lf %10.6lf %10.5lf %10.5lf %6.3f\n",
		aYear,aMonth,aMday,aHour,aMin,aSec,
		name,AZ,EL,SLAT,SLON,RR);
	write_cb(buf,strlen(buf),arg);
}

double
Plan13::get_element(const char *gstr,int gstart,int gstop) {
	int k, glength;
	char gestr[40];
	
	glength = gstop - gstart + 1;
	for ( k = 0; k <= glength; k++ )
		gestr[k] = gstr[gstart+k-1];
	gestr[glength] = 0;

	return atof(gestr);
}

bool
Plan13::load_elements(const char *kep0,const char *kep1,const char *kep2) {
	// for example ...
	// char line1[] = "1 28375U 04025K   09232.55636497 -.00000001  00000-0 12469-4 0   4653";
	// char line2[] = "2 28375 098.0531 238.4104 0083652 290.6047 068.6188 14.40649734270229";

	if ( !kep0 || !*kep0 )
		kep0 = "NoName";
	if ( !kep1 || strlen(kep1) < 69 || kep1[0] != '1' )
		return false;
	if ( !kep2 || strlen(kep2) < 69 || kep2[0] != '2' )
		return false;

	strncpy(name,kep0,sizeof name-1);
	name[sizeof name-1] = 0;

        set_elements(
		get_element(kep1,19,20) + 2000,		// Year
		get_element(kep1,21,32),		// TE: Elapsed time (Epoch - YG)
		get_element(kep2,9,16), 		// IN: Inclination (deg)
         	get_element(kep2,18,25), 		// RA: R.A.A.N (deg)
		get_element(kep2,27,33) * 1.0e-7,	// EC: Eccentricity
		get_element(kep2,35,42),		// WP: Arg perifee (deg)
		get_element(kep2,44,51),		// MA: Mean motion (rev/d)
		get_element(kep2,53,63), 		// MM: Mean motion (rev/d)
         	get_element(kep1,34,43),		// M2: Decay rate (rev/d/d)
		(get_element(kep2,64,68) + ONEPPM),	// RV: Orbit number
		0					// ALON: Sat attitude (deg)
	);
	return true;
}

//////////////////////////////////////////////////////////////////////
// Static default routine for output:
//////////////////////////////////////////////////////////////////////

void
Plan13::write(const void *buf,unsigned bytes,void *arg) {

	::write(1,buf,bytes);
}

//////////////////////////////////////////////////////////////////////
// 
// FOOTPRINT
// 
// In some programs there is a requirement to draw a circle around the sub-
// satellite point to indicate the field of view of the satellite. The
// following routine indicates how to do this. It is coded for clarity, not
// for speed. It would be better to store SIN(A) and COS(A) in a table
// rather than compute them every call. Since the footprint is left-right
// symmetric, only one half of the circle's points need be computed, and
// the other half can be inferred logically.
// 
// The sub-routine's output is a unit vector {X,Y,Z} in geocentric
// coordinates for the I'th point on the Earth's surface. This will then
// need transforming to map coordinates (mercator, spherical, linear or
// whatever), and then screen coordinates to suit the computer.
// 
//////////////////////////////////////////////////////////////////////

#if 0
void
Plan13::footprintOctagon(
  double points[32],
  double SLATin,
  double SLONin,
  double REin,
  double RSin
) {
	printf("SLAT: %lf, SLON: %lf, RE: %lf, RS: %lf\n",
		SLATin,SLONin,REin,RSin);

	double srad = acos(REin/RSin); // Beta in Davidoff diag. 13.2, this is in rad
	printf("srad: %lf\n",srad);

	double cla= cos(rad(SLATin));
	double sla = sin(rad(SLATin));
	double clo = cos(rad(SLONin));
	double slo = sin(rad(SLONin));
	double sra = sin(srad);
	double cra = cos(srad);

	for ( int i = 0; i < 16; i = i+2 ) {
		double a = 2 * M_PI * i / 16;
		printf("\ta: %f\n",a);
		double X = cra;
		double Y = sra*sin(a);
		double Z = sra*cos(a);
		double x = X*cla - Z*sla;
		double y = Y;
		double z = X*sla + Z*cla;
		X = x*clo - y*slo;
		printf("\tX: %f\n",X);
		Y = x*slo + y*clo;
		printf("\tY: %f\n",Y);
		Z = z; 
		points[i] = deg(FNatn(Y,X));
		points[i+1] = deg(asin(Z));

		printf("\t first X: %f, Y: %f, Z: %f, Long: %lf, Lat: %lf\n",
			X,Y,Z,points[i],points[i+1]);
	}
}
#endif

// End plan13.cpp
