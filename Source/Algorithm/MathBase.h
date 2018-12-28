/**
 * @file MathBase.h
 * @brief 1. Vector and matrix operations 2. Earth and Gravity definition 3. Constant 
 *      definition
 * @author Yan GongMing, yinflying
 * @version 1.0
 * @date 2018-12-27
 */

#ifndef _MATHBATH_H
#define _MATHBATH_H

#include <math.h>
#include <stdarg.h>
#include <string.h>

#undef NDEBUG       // Open the assert
#include <assert.h>

// Matrix Max Dimension define
#define MMD     20
#define MMD2    (MMD*MMD)

#ifndef BOOL
typedef int     BOOL;
#endif

#ifndef EPS
#define EPS     2.220446049e-16F
#endif
#define DPS     (DEG/1.0)       // deg/s
#define SEC     (DEG/3600.0)    // arcsec
#define DPH     (DEG/3600.0)    // deg/h
#define G0      9.7803267714
#define UG      (G0/1.0e6)      // ug
#define RE      6378137.0

#define PI      3.14159265358979
#define PI_2    (PI/2.0)
#define PI_4    (PI/4.0)
#define _2PI    (2.0*PI)
#define INF     3.402823466e+30F

#ifndef max
#define max(x,y)        ( (x)>=(y)?(x):(y) )
#endif
#ifndef min
#define min(x,y)        ( (x)<=(y)?(x):(y) )
#endif

#ifndef DEG
#define DEG     (PI/180.0)
#endif

#define pow2(x)         ((x)*(x))
#define asinEx(x)       asin(range(x, -1.0, 1.0))
#define acosEx(x)       acos(range(x, -1.0, 1.0))

// counter-clockwise +-180deg -> clockwise 0~360deg for yaw
#define CC180toC360(yaw)  ( (yaw)>0.0 ? (_2PI-(yaw)) : -(yaw) )
// clockwise 0~360deg -> counter-clockwise +-180deg for yaw
#define C360toCC180(yaw)  ( (yaw)>=PI ? (_2PI-(yaw)) : -(yaw) )

#define LLH(latitude,longitude,height) CVect3((latitude)*PI/180,(longitude)*PI/180,height)
#define PRY(pitch,roll,yaw)            CVect3((pitch)*PI/180,(roll)*PI/180,(yaw)*PI/180)

class CVect3; class CMat3;  class CQuat;  class CVect;  class CMat;  class CGLV;
class CEarth;

// Friend function declination with default parameter
BOOL    IsZero(const CVect3 &v, double eps=EPS);
BOOL    IsZeroXY(const CVect3 &v, double eps=EPS);
CVect3  pow(const CVect3 &v, int k=2);
CVect3  pp2vn(const CVect3 &pos1, const CVect3 &pos0, double ts = 1,
            CEarth *pEth=(CEarth*)NULL);
double  MagYaw(const CVect3 &mag, const CVect3 &att, double declination=0);
CVect   pow(const CVect &v, int k=2);  // power
void    DVMDVafa(const CVect &V, CMat &M, double afa=1.0);  // M = diag(V)*M*diag(V)*afa

double  r2dm(double r);
double  dm2r(double dm);
int     sign(double val, double eps=EPS);
double  range(double val, double minVal, double maxVal);
double  atan2Ex(double y, double x);

extern const CVect3 I31, O31, Ipos;
extern const CQuat  qI;
extern const CMat3  I33, O33;
extern const CVect  On1, O1n;
extern CGLV glv;

// CVect3 friend function
CVect3 abs(const CVect3 &v);
double norm(const CVect3 &v);
double norm1(const CVect3 &v);
double norminf(const CVect3 &v);
double normXY(const CVect3 &v);
CVect3 sqrt(const CVect3 &v);
CVect3 pow(const CVect3 &v, int k);
double dot(const CVect3 &v1, const CVect3 &v2);
CMat3  a2mat(const CVect3 &att);
CVect3 m2att(const CMat3 &Cnb);
CQuat  a2qua(double pitch, double roll, double yaw);
CQuat  a2qua(const CVect3 &att);
CVect3 q2att(const CQuat &qnb);
CQuat  rv2q(const CVect3 &rv);
CVect3 q2rv(const CQuat &q);
CMat3  askew(const CVect3 &v);
CMat3  pos2Cen(const CVect3 &pos);
CVect3 pp2vn(const CVect3 &pos1, const CVect3 &pos0, double ts,CEarth *pEth);
CVect3 MKQt(const CVect3 &sR, const CVect3 &tau);
CMat3  dv2att(const CVect3 &va1,const CVect3 &va2,CVect3 &vb1,const CVect3 &vb2);
CVect3 Alignsb(CVect3 wmm, CVect3 vmm, double latitude);  
double MagYaw(const CVect3 &mag, const CVect3 &att, double declination);
CVect3 xyz2blh(const CVect3 &xyz);
CVect3 blh2xyz(const CVect3 &blh);
CVect3 Vxyz2enu(const CVect3 &Vxyz, const CVect3 &pos);  
CVect3 operator*(double f, const CVect3 &v);
CVect3 operator-(const CVect3 &v);

/**
 * @brief Three Dimension Vector representation and operations
 * @details 3D vecotr could denote various type physical quantities, such as: 
 *  position, velocity, euler angle, this class contain some functions 
 *  belongs to specific physical quantities.
 */
class CVect3 
{
public:
    double i;   ///< First elements of vector
    double j;   ///< Second element of vector
    double k;   ///< Third element of vector

    CVect3(void);
    CVect3(double xyz);
    CVect3(double xx, double yy, double zz);
    CVect3(const double *pdata);
    CVect3(const float *pdata);

    CVect3& operator=(double f);                            // every element equal to a same double
    CVect3& operator=(const double *pf);                    // vector equal to a array
    friend BOOL IsZero(const CVect3 &v, double eps);        // psinsassert if all elements are zeros
    friend BOOL IsZeroXY(const CVect3 &v, double eps);      // psinsassert if x&&y-elements are zeros
    friend BOOL IsNaN(const CVect3 &v);                     // psinsassert if any element is NaN
    CVect3 operator+(const CVect3 &v) const;                // vector addition
    CVect3 operator-(const CVect3 &v) const;                // vector subtraction
    CVect3 operator*(const CVect3 &v) const;                // vector cross multiplication
    CVect3 operator*(const CMat3 &m) const;                 // row-vector multiply matrix
    CVect3 operator*(double f) const;                       // vector multiply scale
    CVect3 operator/(double f) const;                       // vector divide scale
    CVect3 operator/(const CVect3 &v) const;                // vector divide vect3 element by element
    CVect3& operator+=(const CVect3 &v);                    // vector addition
    CVect3& operator-=(const CVect3 &v);                    // vector subtraction
    CVect3& operator*=(double f);                           // vector multiply scale
    CVect3& operator/=(double f);                           // vector divide scale
    CVect3& operator/=(const CVect3 &v);                    // vector divide vect3 element by element
    friend CVect3 operator*(double f, const CVect3 &v);     // scale multiply vector
    friend CVect3 operator-(const CVect3 &v);               // minus
    friend CVect3 abs(const CVect3 &v);                     // abs
    friend double norm(const CVect3 &v);                    // vector norm
    friend double norm1(const CVect3 &v);                   // vector 1-norm
    friend double norminf(const CVect3 &v);
    friend double normXY(const CVect3 &v);                  // vector norm of X & Y components
    friend CVect3 sqrt(const CVect3 &v);                    // sqrt
    friend CVect3 pow(const CVect3 &v, int k);              // power
    friend double dot(const CVect3 &v1, const CVect3 &v2);  // vector dot multiplication
    friend CMat3  a2mat(const CVect3 &att);                  // Euler angles to DCM 
    friend CVect3 m2att(const CMat3 &Cnb);                  // DCM to Euler angles
    friend CQuat  a2qua(double pitch, double roll, double yaw);  // Euler angles to quaternion
    friend CQuat  a2qua(const CVect3 &att);                  // Euler angles to quaternion
    friend CVect3 q2att(const CQuat &qnb);                  // quaternion to Euler angles 
    friend CQuat  rv2q(const CVect3 &rv);                    // rotation vector to quaternion
    friend CVect3 q2rv(const CQuat &q);                     // quaternion to rotation vector
    friend CMat3  askew(const CVect3 &v);                    // askew matrix;
    friend CMat3  pos2Cen(const CVect3 &pos);                // to geographical position matrix
    friend CVect3 pp2vn(const CVect3 &pos1, const CVect3 &pos0, double ts,
            CEarth *pEth);                                  // position difference to velocity
    friend CVect3 MKQt(const CVect3 &sR, const CVect3 &tau);// first order Markov white-noise variance calculation
    friend CMat3  dv2att(const CVect3 &va1, const CVect3 &va2, CVect3 &vb1,
            const CVect3 &vb2);                              // attitude determination using double-vector
    friend CVect3 Alignsb(CVect3 wmm, CVect3 vmm, double latitude);  // align in static-base
    friend double MagYaw(const CVect3 &mag, const CVect3 &att,
            double declination);
    friend CVect3 xyz2blh(const CVect3 &xyz);               // ECEF X/Y/Z to latitude/longitude/height
    friend CVect3 blh2xyz(const CVect3 &blh);               // latitude/longitude/height to ECEF X/Y/Z 
    friend CVect3 Vxyz2enu(const CVect3 &Vxyz, const CVect3 &pos);  // ECEF Vx/Vy/Vz to Ve/Vn/Vu
};

CVect operator~(const CVect &v);      // vector transposition
CVect  abs(const CVect &v);           // vector abs
double norm(const CVect &v);          // vector norm
CVect  pow(const CVect &v, int k);    // power
/**
 * @brief Multidimensional vector representation and operations
 */
class CVect
{
public:
    int row;        ///< Number of vector row(one of row or column must be 1)
    int clm;        ///< Number of vector column
    int rc;         ///< Size of vector
    double dd[MMD]; ///< Vector data

    CVect(void);
    CVect(int row0, int clm0=1);
    CVect(int row0, double f);
    CVect(int row0, double f, double f1, ...);
    CVect(int row0, const double *pf);
    CVect(const CVect3 &v);
    CVect(const CVect3 &v1, const CVect3 v2);

    void Set(double f, ...);
    void Set2(double f, ...);
    CVect operator+(const CVect &v) const;      // vector addition
    CVect operator-(const CVect &v) const;      // vector subtraction
    CVect operator*(double f) const;            // vector multiply scale
    CVect& operator=(double f);                 // every element equal to a same double
    CVect& operator=(const double *pf);         // vector equal to a array
    CVect& operator+=(const CVect &v);          // vector addition
    CVect& operator-=(const CVect &v);          // vector subtraction
    CVect& operator*=(double f);                // vector multiply scale
    CVect operator*(const CMat &m) const;       // row-vector multiply matrix
    CMat operator*(const CVect &v) const;       // 1xn vector multiply nx1 vector, or nx1 vector multiply 1xn vector
    double& operator()(int r);                  // vector element
    friend CVect operator~(const CVect &v);     // vector transposition
    friend CVect abs(const CVect &v);           // vector abs
    friend double norm(const CVect &v);         // vector norm
    friend CVect pow(const CVect &v, int k);  // power
};

void normlize(CQuat *q);             // quaternion norm
CQuat operator~(const CQuat &q);     // quaternion conjugate
CQuat UpDown(const CQuat &q);        // Up-Down the quaternion
/**
 * @brief Quaternion representation and operations
 */
class CQuat
{
public:
    double q0, q1, q2, q3;

    CQuat(void);
    CQuat(double qq0, double qq1=0.0, double qq2=0.0, double qq3=0.0);
    CQuat(const double *pdata);

    CQuat operator+(const CVect3 &phi) const;   // true quaternion add misalign angles
    CQuat operator-(const CVect3 &phi) const;   // calculated quaternion delete misalign angles
    CVect3 operator-(CQuat &quat) const;        // get misalign angles from calculated quaternion & true quaternion
    CQuat operator*(const CQuat &q) const;      // quaternion multiplication
    CVect3 operator*(const CVect3 &v) const;    // quaternion multiply vector
    CQuat& operator*=(const CQuat &q);          // quaternion multiplication
    CQuat& operator-=(const CVect3 &phi);       // calculated quaternion delete misalign angles
    void SetYaw(double yaw=0.0);                // set Euler angles to designated yaw
    friend void normlize(CQuat *q);             // quaternion norm
    friend CQuat operator~(const CQuat &q);     // quaternion conjugate
    friend CQuat UpDown(const CQuat &q);        // Up-Down the quaternion
};

CMat3 operator-(const CMat3 &m);                 // minus
CMat3 operator~(const CMat3 &m);                 // matrix transposition
CMat3 operator*(double f, const CMat3 &m);       // scale multiply matrix
double  det(const CMat3 &m);                     // matrix determination
CMat3   inv(const CMat3 &m);                     // matrix inverse
CVect3  diag(const CMat3 &m);                    // diagonal of a matrix
CMat3   diag(const CVect3 &v);                   // diagonal matrix
CQuat   m2qua(const CMat3 &Cnb);                 // DCM to quaternion
CMat3   q2mat(const CQuat &qnb);                 // quaternion to DCM

/*
 * @brief Three Dimension matrix representation and operations
 * @details 3D matrix could denote directional cosine matrix(DCM), os the class
 *  contains some function about attitude transform.
 */
class CMat3 
{
public:
    double e00, e01, e02, e10, e11, e12, e20, e21, e22;

    CMat3(void);
    CMat3(double xx, double xy, double xz,
          double yx, double yy, double yz,
          double zx, double zy, double zz );
    CMat3(const CVect3 &v0, const CVect3 &v1, const CVect3 &v2);  // M = [v0; v1; v2]

    CMat3 operator+(const CMat3 &m) const;                  // matrix addition
    CMat3 operator-(const CMat3 &m) const;                  // matrix subtraction
    CMat3 operator*(const CMat3 &m) const;                  // matrix multiplication
    CMat3 operator*(double f) const;                        // matrix multiply scale
    CVect3 operator*(const CVect3 &v) const;                // matrix multiply vector
    friend CMat3 operator-(const CMat3 &m);                 // minus
    friend CMat3 operator~(const CMat3 &m);                 // matrix transposition
    friend CMat3 operator*(double f, const CMat3 &m);       // scale multiply matrix
    friend double det(const CMat3 &m);                      // matrix determination
    friend CMat3 inv(const CMat3 &m);                       // matrix inverse
    friend CVect3 diag(const CMat3 &m);                     // diagonal of a matrix
    friend CMat3 diag(const CVect3 &v);                     // diagonal matrix
    friend CQuat  m2qua(const CMat3 &Cnb);                  // DCM to quaternion
    friend CMat3  q2mat(const CQuat &qnb);                  // quaternion to DCM
};

CMat operator~(const CMat &m);               // matrix transposition
void symmetry(CMat &m);                      // matrix symmetrization
double norm1(const CMat &m);                 // 1-norm
CVect diag(const CMat &m);                   // diagonal of a matrix
CMat diag(const CVect &v);                   // diagonal matrix
void RowMul(CMat &m, const CMat &m0, const CMat &m1, int r); // m(r,:)=m0(r,:)*m1
void RowMulT(CMat &m, const CMat &m0, const CMat &m1, int r); // m(r,:)=m0(r,:)*m1'
void DVMDVafa(const CVect &V, CMat &M, double afa);  // M = diag(V)*M*diag(V)*afa
/**
* @brief Multidimensional matrix representation and operations
*/
class CMat
{
public:
    int row;        ///< Number of matrix row
    int clm;        ///< Number of matrix column
    int rc;         ///< Number of matrix size(row times column)
    double dd[MMD2];///< Matrix elements(Max size depend on Macro MMD)

    CMat(void);
    CMat(int row0, int clm0);
    CMat(int row0, int clm0, double f);
    CMat(int row0, int clm0, const double *pf);

    void SetDiag(double f, ...);
    void SetDiag2(double f, ...);
    CMat operator+(const CMat &m) const;                // matrix addition
    CMat operator-(const CMat &m) const;                // matrix subtraction
    CMat operator*(double f) const;                     // matrix multiply scale
    CVect operator*(const CVect &v) const;              // matrix multiply vector
    CMat operator*(const CMat &m) const;                // matrix multiplication
    CMat& operator+=(const CMat &m0);                   // matrix addition
    CMat& operator+=(const CVect &v);                   // matrix + diag(vector)
    CMat& operator-=(const CMat &m0);                   // matrix subtraction
    CMat& operator*=(double f);                         // matrix multiply scale
    CMat& operator++();                                 // 1.0 + diagonal
    double& operator()(int r, int c);                   // get element m(r,c)
    void ZeroRow(int i);                                // set i-row to 0
    void ZeroClm(int j);                                // set j-column to 0
    void SetRow(int i, double f, ...);                  // set i-row from n-double
    void SetRow(int i, const CVect &v);                 // set i-row from vector
    void SetClm(int j, double f, ...);                  // set j-column from n-double
    void SetClm(int j, const CVect &v);                 // set j-column from vector
    CVect GetRow(int i) const;                          // get i-row from matrix
    CVect GetClm(int j) const;                          // get j-column from matrix
    void SetRowVect3(int i, int j, const CVect3 &v);    // set i-row&j...(j+2)-column from CVect3
    void SetClmVect3(int i, int j, const CVect3 &v);    // set i...(i+2)-row&j-column from CVect3
    void SetDiagVect3(int i, int j, const CVect3 &v);   // m(i,j)=v.i, m(i+1,j+1)=v.j, m(i+2,j+2)=v.k;
    void SetMat3(int i, int j, const CMat3 &m);         // set i...(i+2)-row&j...(j+2)-comumn from CMat3
    CMat3 GetMat3(int i, int j) const;                  // get CMat3 from i...(i+2)-row&j...(j+2)-comumn
    void SubAddMat3(int i, int j, const CMat3 &m);      // add i...(i+2)-row&j...(j+2)-comumn with CMat3 m
    friend CMat operator~(const CMat &m);               // matrix transposition
    friend void symmetry(CMat &m);                      // matrix symmetrization
    friend double norm1(const CMat &m);                 // 1-norm
    friend CVect diag(const CMat &m);                   // diagonal of a matrix
    friend CMat diag(const CVect &v);                   // diagonal matrix
    friend void RowMul(CMat &m, const CMat &m0, const CMat &m1, int r); // m(r,:)=m0(r,:)*m1
    friend void RowMulT(CMat &m, const CMat &m0, const CMat &m1, int r); // m(r,:)=m0(r,:)*m1'
    friend void DVMDVafa(const CVect &V, CMat &M, double afa);  // M = diag(V)*M*diag(V)*afa

#ifdef MAT_COUNT_STATISTIC
    static int iCount, iMax;
    ~CMat(void);
#endif
};

/**
 * @brief CGLV(what does it mean?). 
 */
class CGLV
{
public:
    double Re;      ///< Earth equatorial radius(m)
    double Rp;      ///< Earth polar radius(m)
    double f;       ///< Earth oblateness
    double g0;      ///< Earth gravity acceleration(m^2/s) 
    double wie;     ///< Earth average rate of rotation angle(rad/s)
    double e;       ///< Earth first eccentricity
    double e2;      ///< Square of earth first eccentricity
    double ep;      ///< Earth second eccentricity
    double ep2;     ///< Square of earth second eccentricity
    double mg;      ///< Thousandth g(gravitional acceleration unit)(10^-3)
    double ug;      ///< Millionth g(10^-6)
    double deg;     ///< Angle of degree:(PI/180)
    double min;     ///< Angle/time of minute(60)
    double sec;     ///< Angle/time of second(1/3600)
    double hur;     ///< Time of hour(3600)
    double ppm;     ///< Parts per million(10^-6)
    double ppmpsh;  ///< (10^-6/sqrt(3600))
    double dps;     ///< Degree per second(1)
    double dph;     ///< Degree per hour(1/3600)
    double dpsh;    ///< Degree per squared hour(1/sqrt(3600))
    double dphpsh;  ///< Degree per hour per squared hour(1/3600/sqrt(3600))
    double ugpsh;   ///< ug per squared hour(1/sqrt(3600))
    double ugpsHz;  ///< ug per second(1)
    double mpsh;    ///< meter per squared hour(1/sqrt(3600))
    double mpspsh;  ///< meter per second per squared hour(1/1/sqrt(3600))
    double secpsh;  ///< second per squared hour(1/3600/sqrt(3600))

    CGLV(double Re=6378137.0, double f=(1.0/298.257), double wie0=7.2921151467e-5, double g0=9.7803267714);
};

/**
 * @brief 
 */
class CEarth
{
public:
    double a;       ///< Equatorial radius
    double b;       ///< Polar radius
    double f;       ///< Flattening 
    double e;       ///< Eccentricity
    double e2;      ///< Square of eccentricity
    double wie;     ///< Rotational angular velocity
    double g0;      ///< gravity force

    double sl;      ///< sin(lat)
    double sl2;     ///< (sin(lat))^2
    double sl4;     ///< (sin(lat))^4
    double cl;      ///< cos(lat)
    double tl;      ///< tan(lat)
    double RMh;     ///< Meridian radius of curvature
    double RNh;     ///< Prime vertial radius of curvature
    double clRNh;   ///< cos(lat) * RNh
    double f_RMh;   ///< 1/RMh
    double f_RNh;   ///< 1/RNh
    double f_clRNh; ///< 1/(cos(lat)*RNh)

    CVect3 pos;     ///< Position in geodetic coordinate
    CVect3 vn;      ///< Velocity (n-frame to e-frame project in n-frame)
    CVect3 wnie;    ///< w_ie^n, rotation e-frame to i-frame project in n-frame
    CVect3 wnen;    ///< w_en^n, rotation n-frame to e-frame project in n-frame
    CVect3 wnin;    ///< w_in^n, rotation n-frame to i-frame project in n-frame
    CVect3 gn;      ///< g^n, gravity on the surface of earth(Simple gravity field model)
    CVect3 gcc;     ///< Unwanted acceleration under n-frame(a_en^n = f_ib^n - gcc)
    CVect3 *pgn;    ///< g^n, if setting value, gcc vill use it.

    CEarth(double a0=glv.Re, double f0=glv.f, double wie=glv.wie,
            double g0=glv.g0);
    void Update(const CVect3 &pos, const CVect3 &vn=O31);
    CVect3 vn2dpos(const CVect3 &vn, double ts=1.0) const;
};

/**
 * @brief Earth Gravitational Model(Unfinished)
 */
class CEGM
{
public:
    CEGM(void);
    int Init(const char *fileModel);
    int Init(double GM0, double Re0, double ff0, double wie0, double *pC0, double *pS0, int N0);
    ~CEGM(void);
    void Update(double *blh, double *gn, int degree=-1);
};

#endif /* ifndef _MATHBATH_H */
