/***************************************************************************\
*          PREDICT: A satellite tracking/orbital prediction program         *
*          Project started 26-May-1991 by John A. Magliacane, KD2BD         *
*                        Last update: 14-May-2006                           *
*****************************************************************************
*                                                                           *
* This program is free software; you can redistribute it and/or modify it   *
* under the terms of the GNU General Public License as published by the     *
* Free Software Foundation; either version 2 of the License or any later    *
* version.                                                                  *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU         *
* General Public License for more details.                                  *
*                                                                           *
*****************************************************************************
*          See the "CREDITS" file for the names of those who have           *
*   generously contributed their time, talent, and effort to this project.  *
\***************************************************************************/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <fcntl.h>
#include <termios.h>

#define VERSION "2.2.3"

/* Constants used by SGP4/SDP4 code */

#define km2mi       0.621371        /* km to miles */
#define deg2rad     1.745329251994330E-2    /* Degrees to radians */
#define pi      3.14159265358979323846  /* Pi */
#define pio2        1.57079632679489656 /* Pi/2 */
#define x3pio2      4.71238898038468967 /* 3*Pi/2 */
#define twopi       6.28318530717958623 /* 2*Pi  */
#define e6a     1.0E-6
#define tothrd      6.6666666666666666E-1   /* 2/3 */
#define xj2     1.0826158E-3        /* J2 Harmonic (WGS '72) */
#define xj3     -2.53881E-6     /* J3 Harmonic (WGS '72) */
#define xj4     -1.65597E-6     /* J4 Harmonic (WGS '72) */
#define xke     7.43669161E-2
#define xkmper      6.378137E3      /* WGS 84 Earth radius km */
#define xmnpda      1.44E3          /* Minutes per day */
#define ae      1.0
#define ck2     5.413079E-4
#define ck4     6.209887E-7
#define f       3.35281066474748E-3 /* Flattening factor */
#define ge      3.986008E5  /* Earth gravitational constant (WGS '72) */
#define s       1.012229
#define qoms2t      1.880279E-09
#define secday      8.6400E4    /* Seconds per day */
#define omega_E     1.00273790934   /* Earth rotations/siderial day */
#define omega_ER    6.3003879   /* Earth rotations, rads/siderial day */
#define zns     1.19459E-5
#define c1ss        2.9864797E-6
#define zes     1.675E-2
#define znl     1.5835218E-4
#define c1l     4.7968065E-7
#define zel     5.490E-2
#define zcosis      9.1744867E-1
#define zsinis      3.9785416E-1
#define zsings      -9.8088458E-1
#define zcosgs      1.945905E-1
#define zcoshs      1
#define zsinhs      0
#define q22     1.7891679E-6
#define q31     2.1460748E-6
#define q33     2.2123015E-7
#define g22     5.7686396
#define g32     9.5240898E-1
#define g44     1.8014998
#define g52     1.0508330
#define g54     4.4108898
#define root22      1.7891679E-6
#define root32      3.7393792E-7
#define root44      7.3636953E-9
#define root52      1.1428639E-7
#define root54      2.1765803E-9
#define thdt        4.3752691E-3
#define rho     1.5696615E-1
#define mfactor     7.292115E-5
#define sr      6.96000E5   /* Solar radius - km (IAU 76) */
#define AU      1.49597870691E8 /* Astronomical unit - km (IAU 76) */

/* Entry points of Deep() */

#define dpinit   1 /* Deep-space initialization code */
#define dpsec    2 /* Deep-space secular code        */
#define dpper    3 /* Deep-space periodic code       */

/* Flow control flag definitions */

#define ALL_FLAGS              -1
#define SGP_INITIALIZED_FLAG   0x000001 /* not used */
#define SGP4_INITIALIZED_FLAG  0x000002
#define SDP4_INITIALIZED_FLAG  0x000004
#define SGP8_INITIALIZED_FLAG  0x000008 /* not used */
#define SDP8_INITIALIZED_FLAG  0x000010 /* not used */
#define SIMPLE_FLAG            0x000020
#define DEEP_SPACE_EPHEM_FLAG  0x000040
#define LUNAR_TERMS_DONE_FLAG  0x000080
#define NEW_EPHEMERIS_FLAG     0x000100 /* not used */
#define DO_LOOP_FLAG           0x000200
#define RESONANCE_FLAG         0x000400
#define SYNCHRONOUS_FLAG       0x000800
#define EPOCH_RESTART_FLAG     0x001000
#define VISIBLE_FLAG           0x002000
#define SAT_ECLIPSED_FLAG      0x004000

/*
  TODO: This is a refactoring hack to ensure we get consistent before/after
        tests.  Freezes the clock to a specified time, making calculations
        deterministic. (remove when finished along with (2) call sites)
*/
extern char debug_freeze_time;
extern struct tm debug_frozen_tm;

// This struct represents an observation of a particular satellite
// from a particular reference point (on earth) at a particular time
// and is used primarily in the PyPredict code.
typedef struct observation {
    double epoch;
    char orbital_model[5];
    long norad_id;
    char name[25];
    double latitude;
    double longitude;
    double altitude;
    double orbital_velocity;
    double footprint;
    double eclipse_depth;
    double orbital_phase;
    char sunlit;
    long orbit;
    char geostationary;
    double azimuth;
    double elevation;
    double slant_range;
    char visibility;
    char has_aos;
    char decayed;
    double doppler;
} observation;

typedef struct sat_struct{
       char line1[70];       // First line of TLE
       char line2[70];       // Second line of TLE
       char name[25];        // Spacecraft Name
       long catnum;          // Catalog Number (a.k.a NORAD id)
       long setnum;          // Element Set No.
       char designator[10];  // Designator
       int year;             // Reference Epoch (year?) (?)
       double refepoch;      // Reference Epoch (?)
       double incl;          // Inclination
       double raan;          // RAAN
       double eccn;          // Eccentricity
       double argper;        // Arg of Perigee
       double meanan;        // Mean Anomaly
       double meanmo;        // Mean Motion
       double drag;          // Decay Rate
       double nddot6;        // Nddot/6
       double bstar;         // Bstar Drag Term
       long orbitnum;        // Orbit Number
} sat_struct;

extern sat_struct sat;

typedef struct qth_struct {
       char callsign[17];    // Observation Position Call Sign
       double stnlat;        // Observation Position Latitude
       double stnlong;       // Observation Position Longitude
       int stnalt;           // Observation Position Altitude
} qth_struct;

extern qth_struct qth;

// TODO: sat_db doesn't seem to be used much.  Verify it functions and annotate.
typedef struct sat_db_struct {  
       char name[25];
       long catnum;
       char squintflag;
       double alat;
       double alon;
       unsigned char transponders;
       char transponder_name[10][80];
       double uplink_start[10];
       double uplink_end[10];
       double downlink_start[10];
       double downlink_end[10];
       unsigned char dayofweek[10];
       int phase_start[10];
       int phase_end[10];
} sat_db_struct;

extern sat_db_struct sat_db;

/* Global variables for sharing data among functions... */

extern double  tsince, jul_epoch, jul_utc, eclipse_depth,
    sat_azi, sat_ele, sat_range, sat_range_rate,
    sat_lat, sat_lon, sat_alt, sat_vel, phase,
    sun_azi, sun_ele, daynum, fm, fk, age, aostime,
    lostime, ax, ay, az, rx, ry, rz, squint, alat, alon,
    sun_ra, sun_dec, sun_lat, sun_lon, sun_range, sun_range_rate,
    moon_az, moon_el, moon_dx, moon_ra, moon_dec, moon_gha, moon_dv;

extern char    qthfile[50], tlefile[50], dbfile[50], temp[80], output[25],
    ephem[5], sat_sun_status, findsun, calc_squint;

extern int indx, iaz, iel, ma256, isplat, isplong, Flags;

extern long    rv, irk;

extern unsigned char val[256];

/** Type definitions **/

/* Two-line-element satellite orbital data
   structure used directly by the SGP4/SDP4 code. */

typedef struct  {
           double  epoch, xndt2o, xndd6o, bstar, xincl,
               xnodeo, eo, omegao, xmo, xno;
           int     catnr, elset, revnum;
           char    sat_name[25], idesg[9];
        }  tle_t;

/* Geodetic position structure used by SGP4/SDP4 code. */

typedef struct  {
           double lat, lon, alt, theta;
        }  geodetic_t;

/* General three-dimensional vector structure used by SGP4/SDP4 code. */

typedef struct  {
           double x, y, z, w;
        }  vector_t;

/* Common arguments between deep-space functions used by SGP4/SDP4 code. */

typedef struct  {
               /* Used by dpinit part of Deep() */
           double  eosq, sinio, cosio, betao, aodp, theta2,
               sing, cosg, betao2, xmdot, omgdot, xnodot, xnodp;

               /* Used by dpsec and dpper parts of Deep() */
           double  xll, omgadf, xnode, em, xinc, xn, t;

               /* Used by thetg and Deep() */
           double  ds50;
        }  deep_arg_t;

/* Global structure used by SGP4/SDP4 code. */

extern geodetic_t obs_geodetic;

/* Two-line Orbital Elements for the satellite used by SGP4/SDP4 code. */

extern tle_t tle;

/* Functions for testing and setting/clearing flags used in SGP4/SDP4 code */

int isFlagSet(int flag);
int isFlagClear(int flag);
void SetFlag(int flag);
void ClearFlag(int flag);

/* Remaining SGP4/SDP4 code follows... */
int Sign(double arg);
double Sqr(double arg);
double Cube(double arg);
double Radians(double arg);
double Degrees(double arg);
double ArcSin(double arg);
double ArcCos(double arg);
void Magnitude(vector_t *v);
void Vec_Add(vector_t *v1, vector_t *v2, vector_t *v3);
void Vec_Sub(vector_t *v1, vector_t *v2, vector_t *v3);
void Scalar_Multiply(double k, vector_t *v1, vector_t *v2);
void Scale_Vector(double k, vector_t *v);
double Dot(vector_t *v1, vector_t *v2);
double Angle(vector_t *v1, vector_t *v2);
void Cross(vector_t *v1, vector_t *v2 ,vector_t *v3);
void Normalize(vector_t *v);
double AcTan(double sinx, double cosx);
double FMod2p(double x);
double Modulus(double arg1, double arg2);
double Frac(double arg);
int Round(double arg);
double Int(double arg);


void Convert_Sat_State(vector_t *pos, vector_t *vel);
double Julian_Date_of_Year(double year);
double Julian_Date_of_Epoch(double epoch);
int DOY (int yr, int mo, int dy);
double Fraction_of_Day(int hr, int mi, double se);
double Julian_Date(struct tm *cdate);
void Date_Time(double julian_date, struct tm *cdate);
double Delta_ET(double year);
double ThetaG(double epoch, deep_arg_t *deep_arg);
double ThetaG_JD(double jd);
void Calculate_Solar_Position(double time, vector_t *solar_vector);
int Sat_Eclipsed(vector_t *pos, vector_t *sol, double *depth);
void select_ephemeris(tle_t *tle);
void SGP4(double tsince, tle_t * tle, vector_t * pos, vector_t * vel);
void Deep(int ientry, tle_t * tle, deep_arg_t * deep_arg);

void SDP4(double tsince, tle_t * tle, vector_t * pos, vector_t * vel);
void Calculate_User_PosVel(double time, geodetic_t *geodetic, vector_t *obs_pos, vector_t *obs_vel);


void Calculate_Obs(double time, vector_t *pos, vector_t *vel, geodetic_t *geodetic, vector_t *obs_set);
void Calculate_RADec(double time, vector_t *pos, vector_t *vel, geodetic_t *geodetic, vector_t *obs_set);


/* .... SGP4/SDP4 functions end .... */

double FixAngle(double x);
double PrimeAngle(double x);

char *SubString(char *string, char start, char end);
void CopyString(char *source, char *destination, char start, char end);
char *Abbreviate(char *string, int n);

char KepCheck(char *line1, char *line2);
void InternalUpdate(int x);
char *noradEvalue(double value);

void Data2TLE(int x);
double ReadBearing(char *input);
char ReadTLE(char *line0, char *line1, char *line2);
char ReadQTH();
char ReadQTHFile();
char CopyFile(char *source, char *destination);
void SaveQTH();
void SaveTLE();

long DayNum(int m, int d, int y);
double CurrentDaynum();
char *Daynum2String(double daynum);

void FindMoon(double daynum);
void FindSun(double daynum);
void PreCalc(int x);

void Calc();
char AosHappens(int x);
char Decayed(int x, double time);
char Geostationary(int x);
double FindAOS();
double FindLOS();
double FindLOS2();
double NextAOS();
