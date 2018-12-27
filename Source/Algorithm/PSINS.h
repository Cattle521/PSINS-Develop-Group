/* Copyright(c) 2015-2018, by Gongmin Yan, All rights reserved. */
/**
* @file PSINS.h
* @brief PSINS(Precise Strapdown Inertial Navigation System) C++ alogrithm hearder file
* @author Gongmin Yan(Northwestern Polytechnical University, Xi'an, P.R.China)
* @version 1.0
* @date 17/02/2015, 19/07/2017, 11/12/2018
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
#include <glog/logging.h>

#include "MathBase.h"

/************** compile control !!! ***************/
//#define MAT_COUNT_STATISTIC
#define PSINS_IO_FILE
//#define PSINS_RMEMORY
#define PSINS_AHRS_MEMS
#define PSINS_LOW_GRADE_MEMS
//#define CSINSGPSOD_INIT_BY_USER
//#define PSINS_psinsassert
//#define PSINS_UART_PUSH_POP

#define disp(k, nSec, frq)  if(k%(nSec*frq)==0) printf("%d\n", k/frq)

// type re-define
#ifndef BOOL
typedef int     BOOL;
#endif

#ifndef BYTE
typedef unsigned char BYTE;
#endif

// constant define
#ifndef TRUE
#define TRUE    1
#define FALSE   0
#endif

#ifndef NULL
#define NULL    ((void *)0)
#endif

#ifndef PI
#define PI      3.14159265358979
#endif
#define PI_2    (PI/2.0)
#define PI_4    (PI/4.0)
#define _2PI    (2.0*PI)

#ifdef PSINS_psinsassert
    BOOL    psinspsinsassert(BOOL b);
#else
    #define psinsassert(b)  {};
#endif

#ifndef max
#define max(x,y)        ( (x)>=(y)?(x):(y) )
#endif
#ifndef min
#define min(x,y)        ( (x)<=(y)?(x):(y) )
#endif

// class define
class CGLV;
class CEarth;   class CIMU;     class CSINS;    class CAligni0;
class CKalman;  class CSINSKF;  class CSINSTDKF;  class CSINSGPS;
class CIIR;     class CIIRV3;   class CRMemory; class CFileRdWt;    class CUartPP;

// function define
double  diffYaw(double yaw, double yaw0);
double  MKQt(double sR, double tau);
void    IMURFU(CVect3 *pwm, int nSamples, const char *str);
void    IMURFU(CVect3 *pwm, CVect3 *pvm, int nSamples, const char *str);

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
    double operator()(int k);           // get element sqrt(R0(k))
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
    CMat3 Maa, Mav, Map, Mva, Mvv, Mvp, Mpv, Mpp;   // for etm
    CVect3 lvr, vnL, posL; CMat3 CW, MpvCnb;        // for lever arm
    CQuat qnbE; CVect3 attE, vnE, posE;             // for extrapolation
    CRAvar Rwfa;

    CSINS(const CQuat &qnb0=qI, const CVect3 &vn0=O31, const CVect3 &pos0=O31, double tk0=0.0);    // initialization using quat attitude, velocity & position
    void SetTauGA(const CVect3 &tauG, const CVect3 &tauA);
    void Update(const CVect3 *pwm, const CVect3 *pvm, int nSamples, double ts);     // SINS update using Gyro&Acc samples
    void Extrap(const CVect3 &wm=O31, const CVect3 &vm=O31, double ts=0.0);         // SINS fast extrapolation using 1 Gyro&Acc sample
    void lever(const CVect3 &dL=O31);       // lever arm
    void etm(void);                         // SINS error transform matrix coefficients
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
    double kftk;    ///< Kalman update time interval
    double zfdafa;  ///<
    int nq;         ///< Number of states
    int nr;         ///< Number of observations
    int measflag;   ///< 
    CMat Ft;        ///< System matrix in continuous system: x_dot = Ft * x
    CMat Pk;        ///< State's variance-covariance matrix at moment k
    CMat Hk;        ///< Measurement matrix
    CMat Fading;    ///< 
    CVect Xk;       ///< State vector
    CVect Zk;       ///< Observation vector
    CVect Qt;       ///< Main diagonal vector of system noise matrix
    CVect Rt;       ///< Observation noise vector
    CVect rts;      ///< Realtime time span

    CVect Xmax;
    CVect Pmax;
    CVect Pmin;
    CVect Zfd; 
    CVect Zfd0;

    CVect Rmax;     ///< 
    CVect Rmin;     ///< 
    CVect Rbeta;    ///< 
    CVect Rb;              // measurement noise R adaptive

    CVect FBTau;
    CVect FBMax;
    CVect FBXk;
    CVect FBTotal;        // feedback control
    int Rmaxcount[MMD], Rmaxcount0[MMD];
//  CRAvar Ravar;

    CKalman(int nq0, int nr0);
    virtual void Init(void) = 0;                // initialize Qk,Rk,P0...
    virtual void SetFt(void) = 0;               // process matrix setting
    virtual void SetHk(void) = 0;               // measurement matrix setting
    virtual void SetMeas(void) = 0;             // set measurement
    virtual void Feedback(double fbts);         // feed back
    void TimeUpdate(double kfts, int fback=1);  // time update
    int MeasUpdate(double fading=1.0);          // measurement update
    int RAdaptive(int i, double r, double Pr);  // Rt adaptive
    void RPkFading(int i);                      // multiple fading
    void SetMeasFlag(int flag);                 // measurement flag setting
    void XPConstrain(void);                     // Xk & Pk constrain: -Xmax<Xk<Xmax, Pmin<diag(Pk)<Pmax
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
    int Update(const CVect3 *pwm, const CVect3 *pvm, int nSamples, double ts);  // KF Time&Meas Update
    void QtMarkovGA(const CVect3 &tauG, const CVect3 &sRG, const CVect3 &tauA, const CVect3 &sRA); 
    virtual void Miscellanous(void) {};
    virtual void Secret(void);
};

class CSINSTDKF:public CSINSKF
{
public:
    double tdts;
    double Pz0;
    double innovation;
    int iter;
    int ifn;
    int adptOKi;
    int measRes;
    int tdStep;
    int maxStep;
    CMat Fk;
    CMat Pk1; 
    CVect Pxz;
    CVect Qk; 
    CVect Kk;
    CVect Hi;
    CVect tmeas;
    CVect3 meanfn;

    CSINSTDKF(int nq0, int nr0);
    void TDReset(void);
    int TDUpdate(const CVect3 *pwm, const CVect3 *pvm, int nSamples, double ts, int nStep=1);  // Time-Distributed Update
};

class CSINSGPSOD:public CSINSTDKF
{
public:
    CVect3 vbINS, vbOD;
    CMat3 Cbo;          // from body-frame to OD-frame
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
    CVect3 exyzInt, ebMax;

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

/**
 * @brief File read and write operations
 */
class CFileRdWt
{
    static char dirIn[256];     ///< Input file direcotry
    static char dirOut[256];    ///< Output file direcotry
public:
    FILE *f;            ///< file pointer
    char fname[256];    ///< file name
    char line[512];     ///< a line content of current file
    char sstr[64*4];    ///< string format
    double buff[64];    ///< data in one line(double)
    float buff32[64];   ///< data in one line(float)
    int columns;        ///< number of data columns, when read binary, columns<0
    int linelen;        ///< lenth of line

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

class CPTimer   // PSINS Timer
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
