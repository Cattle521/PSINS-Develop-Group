
/* PSINS(Precise Strapdown Inertial Navigation System) C++ alogrithm hearder file PSINS.h */
/*
	Copyright(c) 2015-2018, by Gongmin Yan, All rights reserved.
	Northwestern Polytechnical University, Xi'an, P.R.China
	Date: 17/02/2015, 19/07/2017, 08/10/2018
*/

#ifndef _PSINS_H
#define _PSINS_H

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>

/************** compile control !!! ***************/
//#define MAT_COUNT_STATISTIC
#define PSINS_IO_FILE
//#define PSINS_RMEMORY
#define PSINS_AHRS_MEMS
//#define PSINS_LOW_GRADE_MEMS
//#define CSINSGPSOD_INIT_BY_USER
//#define PSINS_psinsassert
//#define PSINS_UART_PUSH_POP

#define disp(k, nSec, frq)  if(k%(nSec*frq)==0) printf("%d\n", k/frq)

// type re-define
#ifndef BOOL
typedef int		BOOL;
#endif

#ifndef BYTE
typedef unsigned char BYTE;
#endif

// constant define
#ifndef TRUE
#define TRUE	1
#define FALSE	0
#endif

#ifndef NULL
#define NULL	((void *)0)
#endif

#ifndef PI
#define PI		3.14159265358979
#endif
#define PI_2	(PI/2.0)
#define PI_4	(PI/4.0)
#define _2PI	(2.0*PI)

#ifndef DEG
#define DEG		(PI/180.0)
#endif
#define DPS		(DEG/1.0)		// deg/s
#define SEC		(DEG/3600.0)	// arcsec
#define DPH		(DEG/3600.0)	// deg/h
#define G0		9.7803267714
#define UG		(G0/1.0e6)		// ug
#define RE		6378137.0

#ifndef EPS
#define EPS		2.220446049e-16F
#endif
#ifndef INF
#define INF		3.402823466e+30F
#endif

#ifdef PSINS_psinsassert
	BOOL	psinspsinsassert(BOOL b);
#else
	#define psinsassert(b)  {};
#endif

#ifndef max
#define max(x,y)        ( (x)>=(y)?(x):(y) )
#endif
#ifndef min
#define min(x,y)        ( (x)<=(y)?(x):(y) )
#endif

#define CC180toC360(yaw)  ( (yaw)>0.0 ? (_2PI-(yaw)) : -(yaw) )   // counter-clockwise +-180deg -> clockwise 0~360deg for yaw
#define C360toCC180(yaw)  ( (yaw)>=PI ? (_2PI-(yaw)) : -(yaw) )   // clockwise 0~360deg -> counter-clockwise +-180deg for yaw

// class define
class CGLV;
class CVect3;	class CMat3;	class CQuat;
class CEarth;	class CIMU;		class CSINS;	class CAligni0;
class CVect;	class CMat;		class CKalman;	class CSINSKF;		class CSINSTDKF;  class CSINSGPS;
class CIIR;		class CIIRV3;	class CRMemory; class CFileRdWt;	class CUartPP;

// function define
double	r2dm(double r);
double	dm2r(double dm);
int		sign(double val, double eps=EPS);
double	range(double val, double minVal, double maxVal);
double	atan2Ex(double y, double x);
double  diffYaw(double yaw, double yaw0);
double	MKQt(double sR, double tau);
#define pow2(x)			((x)*(x))
#define asinEx(x)		asin(range(x, -1.0, 1.0))
#define acosEx(x)		acos(range(x, -1.0, 1.0))
CVect3	q2att(const CQuat &qnb);
CMat3	diag(const CVect3 &v);
void	IMURFU(CVect3 *pwm, int nSamples, const char *str);
void	IMURFU(CVect3 *pwm, CVect3 *pvm, int nSamples, const char *str);

// Max Matrix Dimension define
#define MMD		20
#define MMD2	(MMD*MMD)

// global variables and functions
extern const CVect3	I31, O31, Ipos;
extern const CQuat	qI;
extern const CMat3	I33, O33;
extern const CVect  On1, O1n;
extern const CGLV	glv;

#define LLH(latitude,longitude,height)		CVect3((latitude)*PI/180,(longitude)*PI/180,height)
#define PRY(pitch,roll,yaw)					CVect3((pitch)*PI/180,(roll)*PI/180,(yaw)*PI/180)

class CGLV
{
public:
	double Re, f, g0, wie;										// the Earth's parameters
	double e, e2, ep, ep2, Rp;
	double mg, ug, deg, min, sec, hur, ppm, ppmpsh;					// commonly used units
	double dps, dph, dpsh, dphpsh, ugpsh, ugpsHz, mpsh, mpspsh, secpsh;

	CGLV(double Re=6378137.0, double f=(1.0/298.257), double wie0=7.2921151467e-5, double g0=9.7803267714);
};

class CVect3 
{
public:
	double i, j, k;

	CVect3(void);
	CVect3(double xyz);
	CVect3(double xx, double yy, double zz);
	CVect3(const double *pdata);

	CVect3& operator=(double f);							// every element equal to a same double
	friend BOOL IsZero(const CVect3 &v, double eps=EPS);	// psinsassert if all elements are zeros
	friend BOOL IsZeroXY(const CVect3 &v, double eps=EPS);	// psinsassert if x&&y-elements are zeros
	friend BOOL IsNaN(const CVect3 &v);						// psinsassert if any element is NaN
	CVect3 operator+(const CVect3 &v) const;				// vector addition
	CVect3 operator-(const CVect3 &v) const;				// vector subtraction
	CVect3 operator*(const CVect3 &v) const;				// vector cross multiplication
	CVect3 operator*(const CMat3 &m) const;					// row-vector multiply matrix
	CVect3 operator*(double f) const;						// vector multiply scale
	CVect3 operator/(double f) const;						// vector divide scale
	CVect3 operator/(const CVect3 &v) const;				// vector divide vect3 element by element
	CVect3& operator=(const double *pf);					// vector equal to a array
	CVect3& operator+=(const CVect3 &v);					// vector addition
	CVect3& operator-=(const CVect3 &v);					// vector subtraction
	CVect3& operator*=(double f);							// vector multiply scale
	CVect3& operator/=(double f);							// vector divide scale
	CVect3& operator/=(const CVect3 &v);					// vector divide vect3 element by element
	friend CVect3 operator*(double f, const CVect3 &v);		// scale multiply vector
	friend CVect3 operator-(const CVect3 &v);				// minus
	friend CVect3 abs(const CVect3 &v);						// abs
	friend double norm(const CVect3 &v);					// vector norm
	friend double normXY(const CVect3 &v);					// vector norm of X & Y components
	friend CVect3 sqrt(const CVect3 &v);					// sqrt
	friend CVect3 pow(const CVect3 &v, int k=2);			// power
	friend double dot(const CVect3 &v1, const CVect3 &v2);	// vector dot multiplication
	friend CMat3 a2mat(const CVect3 &att);					// Euler angles to DCM 
	friend CVect3 m2att(const CMat3 &Cnb);					// DCM to Euler angles 
	friend CQuat a2qua(double pitch, double roll, double yaw);	// Euler angles to quaternion
	friend CQuat a2qua(const CVect3 &att);					// Euler angles to quaternion
	friend CVect3 q2att(const CQuat &qnb);					// quaternion to Euler angles 
	friend CQuat rv2q(const CVect3 &rv);					// rotation vector to quaternion
	friend CVect3 q2rv(const CQuat &q);						// quaternion to rotation vector
	friend CMat3 askew(const CVect3 &v);					// askew matrix;
	friend CMat3 pos2Cen(const CVect3 &pos);				// to geographical position matrix
	friend CVect3 pp2vn(const CVect3 &pos1, const CVect3 &pos0, double ts=1.0, CEarth *pEth=(CEarth*)NULL);  // position difference to velocity
	friend CVect3 MKQt(const CVect3 &sR, const CVect3 &tau);// first order Markov white-noise variance calculation
	friend CMat3 dv2att(CVect3 &va1, const CVect3 &va2, CVect3 &vb1, const CVect3 &vb2);  // attitude determination using double-vector
	friend CVect3 Alignsb(CVect3 wmm, CVect3 vmm, double latitude);  // align in static-base
	friend double MagYaw(const CVect3 &mag, const CVect3 &att, double declination=0);
	friend CVect3 xyz2blh(const CVect3 &xyz);				// ECEF X/Y/Z to latitude/longitude/height
	friend CVect3 blh2xyz(const CVect3 &blh);				// latitude/longitude/height to ECEF X/Y/Z 
	friend CVect3 Vxyz2enu(const CVect3 &Vxyz, const CVect3 &pos);  // ECEF Vx/Vy/Vz to Ve/Vn/Vu
};

class CQuat
{
public:
	double q0, q1, q2, q3;

	CQuat(void);
	CQuat(double qq0, double qq1=0.0, double qq2=0.0, double qq3=0.0);
	CQuat(const double *pdata);

	CQuat operator+(const CVect3 &phi) const;	// true quaternion add misalign angles
	CQuat operator-(const CVect3 &phi) const;	// calculated quaternion delete misalign angles
	CVect3 operator-(CQuat &quat) const;		// get misalign angles from calculated quaternion & true quaternion
	CQuat operator*(const CQuat &q) const;		// quaternion multiplication
	CVect3 operator*(const CVect3 &v) const;	// quaternion multiply vector
	CQuat& operator*=(const CQuat &q);			// quaternion multiplication
	CQuat& operator-=(const CVect3 &phi);		// calculated quaternion delete misalign angles
	void SetYaw(double yaw=0.0);				// set Euler angles to designated yaw
	friend void normlize(CQuat *q);				// quaternion norm
	friend CQuat operator~(const CQuat &q);		// quaternion conjugate
	friend CQuat UpDown(const CQuat &q);		// Up-Down the quaternion
};

class CMat3 
{
public:
	double e00, e01, e02, e10, e11, e12, e20, e21, e22;

	CMat3(void);
	CMat3(double xx, double xy, double xz,
		  double yx, double yy, double yz,
		  double zx, double zy, double zz );
	CMat3(const CVect3 &v0, const CVect3 &v1, const CVect3 &v2);  // M = [v0; v1; v2]

	CMat3 operator+(const CMat3 &m) const;					// matirx addition
	CMat3 operator-(const CMat3 &m) const;					// matirx subtraction
	CMat3 operator*(const CMat3 &m) const;					// matirx multiplication
	CMat3 operator*(double f) const;						// matirx multiply scale
	CVect3 operator*(const CVect3 &v) const;				// matirx multiply vector
	friend CMat3 operator-(const CMat3 &m);					// minus
	friend CMat3 operator~(const CMat3 &m);					// matirx transposition
	friend CMat3 operator*(double f, const CMat3 &m);		// scale multiply matirx
	friend double det(const CMat3 &m);						// matirx determination
	friend CMat3 inv(const CMat3 &m);						// matirx inverse
	friend CVect3 diag(const CMat3 &m);						// diagonal of a matrix
	friend CMat3 diag(const CVect3 &v);						// diagonal matrix
	friend CQuat  m2qua(const CMat3 &Cnb);					// DCM to quaternion
	friend CMat3  q2mat(const CQuat &qnb);					// quaternion to DCM
};

class CVect
{
public:
	int row, clm, rc;
	double dd[MMD];

	CVect(void);
	CVect(int row0, int clm0=1);
	CVect(int row0, double f);
	CVect(int row0, double f, double f1, ...);
	CVect(int row0, const double *pf);
	CVect(const CVect3 &v);
	CVect(const CVect3 &v1, const CVect3 v2);

	void Set(double f, ...);
	void Set2(double f, ...);
	CVect operator+(const CVect &v) const;		// vector addition
	CVect operator-(const CVect &v) const;		// vector subtraction
	CVect operator*(double f) const;			// vector multiply scale
	CVect& operator=(double f);					// every element equal to a same double
	CVect& operator=(const double *pf);			// vector equal to a array
	CVect& operator+=(const CVect &v);			// vector addition
	CVect& operator-=(const CVect &v);			// vector subtraction
	CVect& operator*=(double f);				// vector multiply scale
	CVect operator*(const CMat &m) const;		// row-vector multiply matrix
	CMat operator*(const CVect &v) const;		// 1xn vector multiply nx1 vector, or nx1 vector multiply 1xn vector
	double& operator()(int r);					// vector element
	friend CVect operator~(const CVect &v);		// vector transposition
	friend CVect abs(const CVect &v);			// vector abs
	friend double norm(const CVect &v);			// vector norm
	friend CVect pow(const CVect &v, int k=2);	// power
};

class CMat
{
public:
	int row, clm, rc;
	double dd[MMD2];

	CMat(void);
	CMat(int row0, int clm0);
	CMat(int row0, int clm0, double f);
	CMat(int row0, int clm0, const double *pf);

	void SetDiag(double f, ...);
	void SetDiag2(double f, ...);
	CMat operator+(const CMat &m) const;				// matirx addition
	CMat operator-(const CMat &m) const;				// matirx subtraction
	CMat operator*(double f) const;						// matirx multiply scale
	CVect operator*(const CVect &v) const;				// matirx multiply vector
	CMat operator*(const CMat &m) const;				// matirx multiplication
	CMat& operator+=(const CMat &m0);					// matirx addition
	CMat& operator+=(const CVect &v);					// matirx + diag(vector)
	CMat& operator-=(const CMat &m0);					// matirx subtraction
	CMat& operator*=(double f);							// matirx multiply scale
	CMat& operator++();									// 1.0 + diagonal
	double& operator()(int r, int c);					// get element m(r,c)
	void ZeroRow(int i);								// set i-row to 0
	void ZeroClm(int j);								// set j-column to 0
	void SetRow(int i, double f, ...);
	void SetRow(int i, const CVect &v);					// set i-row from vector
	void SetClm(int j, double f, ...);
	void SetClm(int j, const CVect &v);					// set j-column from vector
	CVect GetRow(int i) const;							// get i-row from matrix
	CVect GetClm(int j) const;							// get j-column from matrix
	void SetRowVect3(int i, int j, const CVect3 &v);	// set i-row&j...(j+2)-column from CVect3
	void SetClmVect3(int i, int j, const CVect3 &v);	// set i...(i+2)-row&j-column from CVect3
	void SetDiagVect3(int i, int j, const CVect3 &v);	// m(i,j)=v.i, m(i+1,j+1)=v.j, m(i+2,j+2)=v.k;
	void SetMat3(int i, int j, const CMat3 &m);			// set i...(i+2)-row&j...(j+2)-comumn from CMat3
	CMat3 GetMat3(int i, int j) const;					// get CMat3 from i...(i+2)-row&j...(j+2)-comumn
	void SubAddMat3(int i, int j, const CMat3 &m);		// add i...(i+2)-row&j...(j+2)-comumn with CMat3 m
	friend CMat operator~(const CMat &m);				// matirx transposition
	friend void symmetry(CMat &m);						// matirx symmetrization
	friend double norm1(CMat &m);						// 1-norm
	friend CVect diag(const CMat &m);					// diagonal of a matrix
	friend CMat diag(const CVect &v);					// diagonal matrix
	friend void RowMul(CMat &m, const CMat &m0, const CMat &m1, int r); // m(r,:)=m0(r,:)*m1
	friend void RowMulT(CMat &m, const CMat &m0, const CMat &m1, int r); // m(r,:)=m0(r,:)*m1'
	friend void DVMDVafa(const CVect &V, CMat &M, double afa=1.0);	// M = diag(V)*M*diag(V)*afa
#ifdef MAT_COUNT_STATISTIC
	static int iCount, iMax;
	~CMat(void);
#endif
};

class CRAvar
{
public:
	#define RAMAX MMD
	int nR0, Rmaxcount[RAMAX], Rmaxflag[RAMAX];
	double ts, R0[RAMAX], Rmax[RAMAX], Rmin[RAMAX], tau[RAMAX], r0[RAMAX];

	CRAvar(void);
	CRAvar(int nR0, int maxCount0=2);
	void set(double r0, double tau, double rmax=0.0, double rmin=0.0, int i=0);
	void set(const CVect3 &r0, const CVect3 &tau, const CVect3 &rmax=O31, const CVect3 &rmin=O31);
	void set(const CVect &r0, const CVect &tau, const CVect &rmax=On1, const CVect &rmin=On1);
	void Update(double r, double ts, int i=0);
	void Update(const CVect3 &r, double ts);
	void Update(const CVect &r, double ts);
	double operator()(int k);			// get element sqrt(R0(k))
};

class CVAR
{
public:
	#define VARMAX 50
	int ipush, imax;
	float array[VARMAX], mean, var;

	CVAR(int imax0=10, float data0=0.0);
	float Update(float data);
};

class CVARn {
public:
	int row, clm, idxpush, rowcnt;
	double **pData, *pd, *Sx, *Sx2, *mx, *stdx;  // sum(x), sum(x^2), mean(x), std(x)
	double stdsf;  // std scalefactor
	CVARn(void);
	CVARn(int row0, int clm0);
	~CVARn(void);
	void Reset(void);
	BOOL Update(const double *pd);
	BOOL Update(double f, ...);
};

class CEarth
{
public:
	double a, b;
	double f, e, e2;
	double wie;

	double sl, sl2, sl4, cl, tl, RMh, RNh, clRNh, f_RMh, f_RNh, f_clRNh;
	CVect3 pos, vn, wnie, wnen, wnin, gn, gcc, *pgn;

	CEarth(double a0=glv.Re, double f0=glv.f, double g0=glv.g0);
	void Update(const CVect3 &pos, const CVect3 &vn=O31);
	CVect3 vn2dpos(const CVect3 &vn, double ts=1.0) const;
};

class CEGM  // Earth Gravitational Model 
{
public:
	CEGM(void);
	int Init(const char *fileModel);
	int Init(double GM0, double Re0, double ff0, double wie0, double *pC0, double *pS0, int N0);
	~CEGM(void);
	void Update(double *blh, double *gn, int degree=-1);
};

class CIMU
{
public:
	int nSamples, prefirst;
	CVect3 phim, dvbm, wm_1, vm_1;

	CIMU(void);
	void Update(const CVect3 *pwm, const CVect3 *pvm, int nSamples);
	friend void IMURFU(CVect3 *pwm, int nSamples, const char *str);
	friend void IMURFU(CVect3 *pwm, CVect3 *pvm, int nSamples, const char *str);
};

class CSINS
{
public:
	double ts, nts, tk, velMax, hgtMin, hgtMax;
	CEarth eth;
	CIMU imu;
	CQuat qnb;
	CMat3 Cnb, Cnb0, Cbn, Kg, Ka;
	CVect3 wib, fb, fn, an, web, wnb, att, vn, vb, pos, eb, db, tauGyro, tauAcc, _betaGyro, _betaAcc;
	CMat3 Maa, Mav, Map, Mva, Mvv, Mvp, Mpv, Mpp;	// for etm
	CVect3 lvr, vnL, posL; CMat3 CW, MpvCnb;		// for lever arm
	CQuat qnbE; CVect3 attE, vnE, posE;				// for extrapolation
	CRAvar Rwfa;

	CSINS(const CQuat &qnb0=qI, const CVect3 &vn0=O31, const CVect3 &pos0=O31, double tk0=0.0);    // initialization using quat attitude, velocity & position
	void SetTauGA(const CVect3 &tauG, const CVect3 &tauA);
	void Update(const CVect3 *pwm, const CVect3 *pvm, int nSamples, double ts);		// SINS update using Gyro&Acc samples
	void Extrap(const CVect3 &wm=O31, const CVect3 &vm=O31, double ts=0.0);			// SINS fast extrapolation using 1 Gyro&Acc sample
	void lever(const CVect3 &dL=O31);		// lever arm
	void etm(void);							// SINS error transform matrix coefficients
};

class CIIR
{
public:
	#define IIRnMax 6
	int n;
	double b[IIRnMax], a[IIRnMax], x[IIRnMax], y[IIRnMax];

	CIIR(void);
	CIIR(double *b0, double *a0, int n0);
	double Update(double x0);
};

class CIIRV3
{
public:
	CIIR iir0, iir1, iir2;
	CVect3 y;

	CIIRV3(void);
	CIIRV3(double *b0, double *a0, int n0,
		   double *b1=(double*)NULL, double *a1=(double*)NULL, int n1=0, 
		   double *b2=(double*)NULL, double *a2=(double*)NULL, int n2=0);
	CVect3 Update(CVect3 &x);
};

class CAligni0
{
public:
	int velAid, t0, t1, t2;
	CVect3 vel0, wmm, vmm, vib0, vi0, Pib01, Pib02, Pi01, Pi02, tmpPib0, tmpPi0;
	CQuat qib0b;
	CEarth eth;
	CIMU imu;
	double tk;
	CQuat qnb0, qnb;

	CAligni0(const CVect3 &pos0=O31, const CVect3 &vel0=O31, int velAid=0);
	CQuat Update(const CVect3 *wm, const CVect3 *vm, int nSamples, double ts, const CVect3 &vel=O31);
};

class CKalman
{
public:
	double kftk, zfdafa;
	int nq, nr, measflag;
	CMat Ft, Pk, Hk, Fading;
	CVect Xk, Zk, Qt, Rt, rts, Xmax, Pmax, Pmin, Zfd, Zfd0,
		Rmax, Rmin, Rbeta, Rb,				// measurement noise R adaptive
		FBTau, FBMax, FBXk, FBTotal;		// feedback control
	int Rmaxcount[MMD], Rmaxcount0[MMD];
//	CRAvar Ravar;

	CKalman(int nq0, int nr0);
	virtual void Init(void) = 0;				// initialize Qk,Rk,P0...
	virtual void SetFt(void) = 0;				// process matrix setting
	virtual void SetHk(void) = 0;				// measurement matrix setting
	virtual void SetMeas(void) = 0;				// set measurement
	virtual void Feedback(double fbts);			// feed back
	void TimeUpdate(double kfts, int fback=1);	// time update
	int MeasUpdate(double fading=1.0);			// measurement update
	int RAdaptive(int i, double r, double Pr); // Rt adaptive
	void RPkFading(int i);						// multiple fading
	void SetMeasFlag(int flag);					// measurement flag setting
	void XPConstrain(void);						// Xk & Pk constrain: -Xmax<Xk<Xmax, Pmin<diag(Pk)<Pmax
	virtual void Miscellanous(void) {};
};

class CSINSKF:public CKalman
{
public:
	CSINS sins;

	CSINSKF(int nq0, int nr0);
	virtual void Init(void) {} ;
	virtual void Init(const CSINS &sins0, int grade=0);
	virtual void SetFt(void);
	virtual void SetHk(void);
	virtual void Feedback(double fbts);
	int Update(const CVect3 *pwm, const CVect3 *pvm, int nSamples, double ts);	// KF Time&Meas Update
	void QtMarkovGA(const CVect3 &tauG, const CVect3 &sRG, const CVect3 &tauA, const CVect3 &sRA); 
};

class CSINSTDKF:public CSINSKF
{
public:
	double tdts, Pz0, innovation;
	int iter, ifn, adptOKi, measRes, tdStep, maxStep;
	CMat Fk, Pk1; 
	CVect Pxz, Qk, Kk, Hi, tmeas;
	CVect3 meanfn;

	CSINSTDKF(int nq0, int nr0);
	void TDReset(void);
	int TDUpdate(const CVect3 *pwm, const CVect3 *pvm, int nSamples, double ts, int nStep=1);  // Time-Distributed Update
};

class CSINSGPSOD:public CSINSTDKF
{
public:
	CVect3 vbINS, vbOD;
	CMat3 Cbo;			// from body-frame to OD-frame
	double tODInt;
	BOOL measGPSvnValid, measGPSposValid, measODvbValid, measMAGyawValid;

	CSINSGPSOD(void);
	virtual void Init(const CSINS &sins0, int grade=0);
	virtual void SetMeas(void);
	void SetMeasGPS(const CVect3 &pgps=O31, const CVect3 &vgps=O31);
	void SetMeasOD(double dSod, double ts);
	void SetMeasYaw(double ymag);
};

#ifdef PSINS_AHRS_MEMS

class CMahony
{
public:
	double tk, Kp, Ki;
	CQuat qnb;
	CMat3 Cnb;
	CVect3 exyzInt;

	CMahony(double tau=4.0, const CQuat &qnb0=qI);
	void SetTau(double tau=4.0);
	void Update(const CVect3 &wm, const CVect3 &vm, double ts, const CVect3 &mag=O31);
	void Update(const CVect3 &gyro, const CVect3 &acc, const CVect3 &mag, double ts);
};

class CQEAHRS:public CKalman
{
public:
	CMat3 Cnb;

	CQEAHRS(double ts);
	void Update(const CVect3 &gyro, const CVect3 &acc, const CVect3 &mag, double ts);
};

#endif // PSINS_AHRS_MEMS

#ifdef PSINS_IO_FILE

class CFileRdWt
{
	static char dirIn[256], dirOut[256];
public:
	FILE *f;
	char fname[256], line[512], sstr[64*4];
	double buff[64];
	int columns, linelen;

	static void Dir(const char *dirI, const char *dirO=(const char*)NULL);
	CFileRdWt(void);
	CFileRdWt(const char *fname0, int columns0=0);
	void Init(const char *fname0, int columns0=0);
	int load(int lines=1, BOOL txtDelComma=1);
	int loadf32(int lines=1);
	int getl(void);  // get a line
	CFileRdWt& operator<<(double d);
	CFileRdWt& operator<<(const CVect3 &v);
	CFileRdWt& operator<<(const CVect &v);
	CFileRdWt& operator<<(const CMat &m);
	CFileRdWt& operator<<(const CRAvar &R);
	CFileRdWt& operator<<(const CAligni0 &aln);
	CFileRdWt& operator<<(const CSINS &sins);
#ifdef PSINS_AHRS_MEMS
	CFileRdWt& operator<<(const CMahony &ahrs);
	CFileRdWt& operator<<(const CQEAHRS &ahrs);
#endif
#ifdef PSINS_UART_PUSH_POP
	CFileRdWt& operator<<(const CUartPP &uart);
#endif
	CFileRdWt& operator<<(const CKalman &kf);
	CFileRdWt& operator>>(double &d);
	CFileRdWt& operator>>(CVect3 &v);
	CFileRdWt& operator>>(CVect &v);
	CFileRdWt& operator>>(CMat &m);
	~CFileRdWt();
};

char* time2fname(void);

#endif // PSINS_IO_FILE

#ifdef PSINS_RMEMORY

#define MAX_RECORD_BYTES 128

class CRMemory
{
	BYTE *pMemStart, *pMemEnd;
	BYTE pushLen, popLen, recordLen;
	long memLen, dataLen;
public:
	BYTE *pMemPush, *pMemPop, pushBuf[MAX_RECORD_BYTES], popBuf[MAX_RECORD_BYTES];

	CRMemory(BYTE *pMem, long memLen0, BYTE recordLen0=0);
	BOOL push(const BYTE *p=(const BYTE*)NULL);
	BYTE pop(BYTE *p=(BYTE*)NULL);
};

#endif // PSINS_RMEMORY

#ifdef PSINS_UART_PUSH_POP

class CUartPP
{
public:
#define UARTFRMLEN  (50*4)
#define UARTBUFLEN  (UARTFRMLEN*20)
	unsigned char head[2], popbuf[UARTFRMLEN], buf[UARTBUFLEN], chksum;
	int pushIdx, popIdx, frameLen, overflow, getframe;
	int csflag, cs0, cs1, css;   // 0: no checksum, 1: unsigned char sum, 2: crc8; popbuf[css]==popbuf[cs0]+...+popbuf[cs1] 
//unsigned short chksum;

	BOOL checksum(const unsigned char *pc);
	CUartPP(int frameLen0, unsigned short head0=0x55aa);  // NOTE: char {0xaa 0x55}
	int nextIdx(int idx);
	int push(const unsigned char *buf0, int len);
	int pop(unsigned char *buf0=(unsigned char*)NULL);
};

#endif // PSINS_UART_PUSH_POP

class CPTimer	// PSINS Timer
{
	double tCurrent, tMax;
	int isStart, autoReset;

public:
	BOOL overflow;
	CPTimer(double tMax0=1.0, BOOL isStart0=0, BOOL autoReset0=0);
	void Start(double tMax0=0.0);
	void Stop(void);
	BOOL Update(double tStep);
};

typedef struct {
	unsigned short head, chksum; 
	float t, gx, gy, gz, ax, ay, az, magx, magy, magz, bar, pch, rll, yaw, ve, vn, vu,
		lon0, lon1, lat0, lat1, hgt, gpsve, gpsvn, gpsvu, gpslon0, gpslon1, gpslat0, gpslat1, gpshgt,
		gpsstat, gpsdly, tmp, rsv;
} PSINSBoard;

#endif
