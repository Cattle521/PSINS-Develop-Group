/**
 * @file MathBase.cpp
 * @brief 1. Vector and matrix operations 2. Earth and Gravity definition 3. Constant 
 *      definition
 * @author Yan GongMing, yinflying
 * @version 1.0
 * @date 2018-12-27
 */
#include "PSINSBase.h"

const CVect3 I31(1.0), O31(0.0), Ipos(1.0/RE,1.0/RE,1.0);
const CQuat  qI(1.0,0.0,0.0,0.0);
const CMat3  I33(1,0,0, 0,1,0, 0,0,1), O33(0,0,0, 0,0,0, 0,0,0);
const CVect  On1(MMD,0.0), O1n=~On1;
CGLV glv;

/**
 * @brief Initialize the global constant
 * @param[in] Re    Earth equatorial radius(m), default: 6378137.0
 * @param[in] f     Earth oblateness, default: 1.0/298.257
 * @param[in] wie0  Earth average rate of rotation angle(rad/s), 
 *                  default:7.2921151467e-5
 * @param[in] g0    Earth gravity acceleration(m^2/s), default: 9.7803267714
 */
CGLV::CGLV(double Re, double f, double wie0, double g0)
{
    this->Re = Re; this->f = f; this->wie = wie0; this->g0 = g0;
    Rp = (1-f)*Re;
    e = sqrt(2*f-f*f); e2 = e*e;
    ep = sqrt(Re*Re-Rp*Rp)/Rp; ep2 = ep*ep;
    mg     = g0/1000.0;
    ug     = mg/1000.0;
    deg    = PI/180.0;
    min    = deg/60.0;
    sec    = min/60.0;
    ppm    = 1.0e-6;
    hur    = 3600.0;
    dps    = deg/1.0;
    dph    = deg/hur;
    dpsh   = deg/sqrt(hur);
    dphpsh = dph/sqrt(hur);
    ugpsHz = ug/sqrt(1.0);
    ugpsh  = ug/sqrt(hur);
    mpsh   = 1/sqrt(hur);
    mpspsh = 1/1/sqrt(hur);
    ppmpsh = ppm/sqrt(hur);
    secpsh = sec/sqrt(hur);
}

/***************************  class CVect3  *********************************/
CVect3::CVect3(void)
{
}
/**
 * @brief Create 3D-vector three elements are equal: [xyz xyz xyz]
 * @param[in] xyz Input elements value
 */
CVect3::CVect3(double xyz)
{
    i=j=k=xyz;
}

/**
 * @brief Create 3D-vector: [xx yy zz]
 * @param[in] xx Fisrt input element
 * @param[in] yy Second input element
 * @param[in] zz Third input element
 */
CVect3::CVect3(double xx, double yy, double zz)
{
    i=xx, j=yy, k=zz;
}

/**
 * @brief Create 3D-vector from a point to an array
 * @param[in] pdata input array
 */
CVect3::CVect3(const double *pdata)
{
    i = *pdata++, j = *pdata++, k = *pdata;
}

/**
 * @see CVect3::CVect3
 */
CVect3::CVect3(const float *pdata)
{
    i = *pdata++, j = *pdata++, k = *pdata;
}


/**
 * @brief Determine if the vector is zero vector
 * @param[in] v Input vector
 * @param[in] eps Zero threshold, e.g. 1e-20
 * @return True if v is zero vector, False if not
 */
BOOL IsZero(const CVect3 &v, double eps)
{
    return (v.i<eps&&v.i>-eps && v.j<eps&&v.j>-eps && v.k<eps&&v.k>-eps);
}

/**
 * @brief Determine if first two elements of vector is zero
 * @param v     Input vector
 * @param eps   Zero threshold, e.g. 1e-20
 * @return True if v meet the condition, Flase if not 
 */
BOOL IsZeroXY(const CVect3 &v, double eps)
{
    return (v.i<eps&&v.i>-eps && v.j<eps&&v.j>-eps);
}

/**
 * @brief Determine if vector is NaN
 * @param[in] v Input vector
 * @return True if v is NaN, Flase if not 
 * @warning  This function do NOT work at present
 */
BOOL IsNaN(const CVect3 &v)
{
    return 0; //(_isnan(i) || _isnan(j) || _isnan(k));
}

/**
 * @brief Overload operator +, [x1+x2 y1+y2 z1+z2] = [x1 y1 z1] + [x2 y2 z2]
 * @param[in] v Input vector
 * @return Plus result
 */
CVect3 CVect3::operator+(const CVect3 &v) const
{
    return CVect3(this->i+v.i, this->j+v.j, this->k+v.k);
}

/**
 * @brief Overload operator -, [x1-x2 y1-y2 z1-z2] = [x1 y1 z1] - [x2 y2 z2]
 * @param[in] v Input vector
 * @return Subtraction result
 */
CVect3 CVect3::operator-(const CVect3 &v) const
{
    return CVect3(this->i-v.i, this->j-v.j, this->k-v.k);
}

/**
 * @brief Overload operator * as cross product of vector
 * @param[in] v Input vector
 * @return Cross product result
 */
CVect3 CVect3::operator*(const CVect3 &v) const
{
    return CVect3(this->j*v.k - this->k*v.j, this->k*v.i - this->i*v.k, 
            this->i*v.j - this->j*v.i);
}
    
/**
 * @brief Overload operator * as dot product
 * @param[in] f Input number
 * @return Dot product result
 */
CVect3 CVect3::operator*(double f) const
{
    return CVect3(i*f, j*f, k*f);
}

/**
 * @brief Overload operator * as Matrix product(row vector * Mastrix)
 * @param[in] m Input 3D-matrix
 * @return Matrix product result(row vector)
 */
CVect3 CVect3::operator*(const CMat3 &m) const
{
    return CVect3(i*m.e00 + j*m.e10 + k*m.e20,i*m.e01 + j*m.e11 + k*m.e21,
            i*m.e02 + j*m.e12 + k*m.e22);
}
    
/**
 * @brief Overload operator / as right-array division, [x/f y/f z/f] = [x y z]/f
 * @param[in] f Input number
 * @return Right-array division result
 */
CVect3 CVect3::operator/(double f) const
{
    return CVect3(i/f, j/f, k/f);
}

/**
 * @brief Overload operator / as right-arrary division, [x1/x2 y1/y2 z1/z2] = 
 * [x1 y1 z1] - [x2 y2 z2]
 * @param[in] v Input vector
 * @return Right-array division result
 */
CVect3 CVect3::operator/(const CVect3 &v) const
{
    return CVect3(i/v.i, j/v.j, k/v.k);
}

/**
 * @brief Overload operator =, assign equal value [f f f] for vector
 * @param[in] f Input number
 * @return Assigned vector
 */
CVect3& CVect3::operator=(double f)
{
    i = j = k = f;
    return *this;
}

/**
 * @brief Overload operator =, assign value from an array(size 3)
 * @param[in] pf Input arrary
 * @return Assigned vector
 */
CVect3& CVect3::operator=(const double *pf)
{
    i = *pf++, j = *pf++, k = *pf;
    return *this;
}

/**
 * @brief Overload operator +=, Add and assign for two vector
 * @param[in] v  Input vector
 * @return Added vector
 */
CVect3& CVect3::operator+=(const CVect3 &v)
{ 
    i += v.i, j += v.j, k += v.k;
    return *this;
}

/**
 * @brief Overload ooperator -=
 * @param[in] v Input vector
 * @return Substracted vector
 */
CVect3& CVect3::operator-=(const CVect3 &v)
{ 
    i -= v.i, j -= v.j, k -= v.k;
    return *this;
}

/**
 * @brief  Overload operator *= for dot product
 * @param[in] f Input number
 * @return Dot producted vector
 */
CVect3& CVect3::operator*=(double f)
{ 
    i *= f, j *= f, k *= f;
    return *this;
}

/**
 * @brief  Overload operator /= for right division
 * @param[in] f Input number
 * @return right divided vector
 */
CVect3& CVect3::operator/=(double f)
{ 
    i /= f, j /= f, k /= f;
    return *this;
}

/**
 * @brief  Overload operator /= for right division
 * @param[in] v Input vector
 * @return right divided vector
 */
CVect3& CVect3::operator/=(const CVect3 &v)
{ 
    i /= v.i, j /= v.j, k /= v.k;
    return *this;
}

/**
 * @brief Overload operator * for dot production(number x vector)
 * @param[in] f Input number
 * @param[in] v Input vector
 * @return Overload vector
 */
CVect3 operator*(double f, const CVect3 &v)
{
    return CVect3(v.i*f, v.j*f, v.k*f);
}
    
/**
 * @brief Take a negative on vector
 * @param[in] v Input vector
 * @return Negative vector
 */
CVect3 operator-(const CVect3 &v)
{
    return CVect3(-v.i, -v.j, -v.k);
}

/**
 * @brief Square root for each element in the vector
 * @param[in] v Input vector
 * @return Square root vector 
 */
CVect3 sqrt(const CVect3 &v)
{
    return CVect3(sqrt(v.i), sqrt(v.j), sqrt(v.k));
}

/** 
 * @brief Power for each element in the vector
 * @param[in] v Input vector
 * @param[in] k order of the power(optional, default 2)
 * @return  Power of vector
 */
CVect3 pow(const CVect3 &v, int k)
{
    CVect3 pp = v;
    for(int i=1; i<k; i++)
    {
        pp.i *= v.i, pp.j *= v.j, pp.k *= v.k;
    }
    return pp;
}

/**
 * @brief Absolute value for each element in the vector
 * @param[in] v Input vector
 * @return Absolute value of vector
 */
CVect3 abs(const CVect3 &v)
{
    CVect3 res;
    res.i = v.i>0.0 ? v.i : -v.i;
    res.j = v.j>0.0 ? v.j : -v.j;
    res.k = v.k>0.0 ? v.k : -v.k;
    return res;
}

/**
 * @brief Euclidean distance of vector
 * @param[in] v Input vector
 * @return Euclidean distance of vector 
 */
double norm(const CVect3 &v)
{
    return sqrt(v.i*v.i + v.j*v.j + v.k*v.k);
}

/**
 * @brief Infinity norm: Maximum value of the absolute values of all elements
 * @param[in] v Input vector
 * @return Infinity norm result
 */
double norminf(const CVect3 &v)
{
    double nm=0;
    if(v.i>0)  nm = max(nm, v.i);   else   nm = max(nm, -v.i);
    if(v.j>0)  nm = max(nm, v.j);   else   nm = max(nm, -v.j);
    if(v.k>0)  nm = max(nm, v.k);   else   nm = max(nm, -v.k);
    return nm;
}
/**
 * @brief 1_norm, sum of absolute values of all vector elements
 * @param[in] v Input vector 
 * @return 1_norm result
 */
double norm1(const CVect3 &v)
{
    double sum = 0;
    sum += (v.i>0.0?v.i:-v.i);
    sum += (v.j>0.0?v.j:-v.j);
    sum += (v.k>0.0?v.k:-v.k);
    return sum;
}
/**
 * @brief Vector norm of X & Y components
 * @param[in] v Input vector
 * @return Vector norm of X & Y components
 */
double normXY(const CVect3 &v)
{
    return sqrt(v.i*v.i + v.j*v.j);
}

/**
 * @brief Dot product of two 3-D vector
 * @param[in] v1 First vector
 * @param[in] v2 Second vector
 * @return Dot product to two 3-D vector
 */
double dot(const CVect3 &v1, const CVect3 &v2)
{
    return (v1.i*v2.i + v1.j*v2.j + v1.k*v2.k);
}

/**
 * @brief Convert rotation vector to quaternion
 * @details if a small angle is rotated(<1 deg), use Taylor series expansion
 *      of sine or cosine, otherwise, use sine or cosine
 * @param[in] rv Rotation vector
 * @return Quaternion
 */
CQuat rv2q(const CVect3 &rv)
{
    const double F1 = 2 * 1;    // define Fk = 2^k * k!
    const double F2 = F1*2 * 2;
    const double F3 = F2*2 * 3;
    const double F4 = F3*2 * 4;
    const double F5 = F4*2 * 5;
    double n2 = rv.i*rv.i+rv.j*rv.j+rv.k*rv.k, c, f;
    if(n2<(PI/180.0*PI/180.0))  // 0.017^2 
    {
        double n4=n2*n2;
        c = 1.0 - n2*(1.0/F2) + n4*(1.0/F4);
        f = 0.5 - n2*(1.0/F3) + n4*(1.0/F5);
    }
    else
    {
        double n_2 = sqrt(n2)/2.0;
        c = cos(n_2);
        f = sin(n_2)/n_2*0.5;
    }
    return CQuat(c, f*rv.i, f*rv.j, f*rv.k);
}

/**
 * @brief Get skew-symmetric matrix of input 3D-vector
 * @param[in] v Input vector
 * @return Skew-sysmmetric matrix
 */
CMat3 askew(const CVect3 &v)
{
    return CMat3(0,  -v.k, v.j, 
                 v.k, 0.0,  -v.i,
                -v.j, v.i, 0);
}

/**
 * @brief Get DCM of e-frame to n-frame from geodetic coordinate postion
 * @param[in] pos Geodetic coordinate postion.
 * @return Cbn
 */
CMat3 pos2Cen(const CVect3 &pos)
{
    double si = sin(pos.i), ci = cos(pos.i), sj = sin(pos.j), cj = cos(pos.j);
    return CMat3(   -sj, -si*cj,  ci*cj,  
                     cj, -si*sj,  ci*sj,  
                     0,   ci,     si      );    //Cen
}

/**
 * @brief Calculate average velocity via position difference
 * @param[in] pos1 End postion, under geodetic coordinate
 * @param[in] pos0 Start postion, under geodetic coordinate
 * @param[in] ts Time span between two position(s)
 * @param[in] pEth Earth parameter
 * @return Velocity under n-frame(m/s)
 */
CVect3 pp2vn(const CVect3 &pos1, const CVect3 &pos0, double ts, CEarth *pEth)
{
    double sl, cl, sl2, sq, sq2, RMh, RNh, clRNh;
    if(pEth)
    {
        RMh = pEth->RMh; clRNh = pEth->clRNh;
    }
    else
    {
        sl=sin(pos0.i); cl=cos(pos0.i); sl2=sl*sl;
        sq = 1-glv.e2*sl2; sq2 = sqrt(sq);
        RMh = glv.Re*(1-glv.e2)/sq/sq2+pos0.k;
        RNh = glv.Re/sq2+pos0.k;    clRNh = cl*RNh;
    }
    CVect3 vn = pos1 - pos0;
    return CVect3(vn.j*clRNh/ts, vn.i*RMh/ts, vn.k/ts);
}

/**
 * @brief Calculate yaw angle using geomagnetic field
 * @param[in] mag
 * @param[in] att Euler attitude
 * @param[in] declination
 * @return 
 */
double MagYaw(const CVect3 &mag, const CVect3 &att, double declination)
{
    CVect3 attH(att.i, att.j, 0.0);
    CVect3 magH = a2mat(attH)*mag;
    double yaw = 0.0;
    if(attH.i<(80.0*DEG)&&attH.i>-(80.0*DEG))
    {
        yaw = atan2Ex(magH.i, magH.j) + declination;
        if(yaw>PI)       yaw -= _2PI;
        else if(yaw<-PI) yaw += _2PI;
    }
    return yaw;
}

/**
 * @brief Convert ECEF Cartesion coordinate to geodetic coordinate
 * @param[in] xyz  ECEF Castesion coordinate(m)
 * @return Geodetic coordinate
 */
CVect3 xyz2blh(const CVect3 &xyz)
{
    double s = normXY(xyz), theta = atan2(xyz.k*glv.Re, s*glv.Rp),
        s3 = sin(theta), c3 = cos(theta); s3 = s3*s3*s3, c3 = c3*c3*c3;
    if(s<(6378137.0*1.0*DEG))  return O31;
    double L = atan2(xyz.j, xyz.i);
    double B = atan2(xyz.k+glv.ep2*glv.Rp*s3, s-glv.e2*glv.Re*c3);
    double sB = sin(B), cB = cos(B), N = glv.Re/sqrt(1-glv.e2*sB*sB);
    return CVect3(B, L, s/cB-N);
}

/**
 * @brief Convert geodetic coordinate to ECEF cartesion coordinate
 * @param[in] blh geodetic coordinate
 * @return Cartesion coordinate
 */
CVect3 blh2xyz(const CVect3 &blh)
{
    double sB = sin(blh.i), cB = cos(blh.i), sL = sin(blh.j), cL = cos(blh.j),
        N = glv.Re/sqrt(1-glv.e2*sB*sB);
    return CVect3((N+blh.k)*cB*cL, (N+blh.k)*cB*sL, (N*(1-glv.e2)+blh.k)*sB);
}

/**
 * @brief Project ECEF vector to ENU coordinate system
 * @param[in] Vxyz Vector under ECEF 
 * @param[in] pos Geodetic coordinate pos
 * @return vector under ENU-coordinate system
 */
CVect3 Vxyz2enu(const CVect3 &Vxyz, const CVect3 &pos)
{
    return Vxyz*pos2Cen(pos);
}

/***************************  class CQuat  *********************************/
CQuat::CQuat(void){}

/**
 * @brief Quaternion initialization
 */
CQuat::CQuat(double q0, double q1, double q2, double q3):
    q0(q0),q1(q1),q2(q2),q3(q3){}

/**
 * @brief Quaternion initialization by array
 * @param[in] pdata array containing four elements
 */
CQuat::CQuat(const double *pdata)
{
    q0=*pdata++, q1=*pdata++, q2=*pdata++, q3=*pdata;
}

/**
 * @brief Quaternion add rotation vector(Qac = Qab + phi_bc)
 * @param[in] phi Rotation vector
 * @return Quaternion
 */
CQuat CQuat::operator+(const CVect3 &phi) const
{
    CQuat qtmp = rv2q(-phi);
    return qtmp*(*this);
}

/**
 * @brief  Quaternion minus rotation vector(Qac = Qab - phi_cb)
 * @param[in] phi Rotation vector
 * @return Quaternion
 */
CQuat CQuat::operator-(const CVect3 &phi) const
{
    CQuat qtmp = rv2q(phi);
    return qtmp*(*this);
}

/**
 * @brief Get the angle difference between two quaternions(phi_ac = Qab - Qbc)
 * @param[in] quat
 * @return Rotation angle
 */
CVect3 CQuat::operator-(CQuat &quat) const
{
    CQuat dq;
    
    dq = quat*(~(*this));
    if(dq.q0<0)
    {
        dq.q0=-dq.q0, dq.q1=-dq.q1, dq.q2=-dq.q2, dq.q3=-dq.q3;
    }
    double n2 = acos(dq.q0), f;
    if(sign(n2)!=0 )
    {
        f = 2.0/(sin(n2)/n2);
    }
    else
    {
        f = 2.0;
    }
    return CVect3(dq.q1,dq.q2,dq.q3)*f;
}

/**
 * @brief Quaternion multiplication
 * @param[in] quat multiplied Quaternion
 * @return Multiplied result
 */
CQuat CQuat::operator*(const CQuat &quat) const
{
    CQuat qtmp;
    qtmp.q0 = q0*quat.q0 - q1*quat.q1 - q2*quat.q2 - q3*quat.q3;
    qtmp.q1 = q0*quat.q1 + q1*quat.q0 + q2*quat.q3 - q3*quat.q2;
    qtmp.q2 = q0*quat.q2 + q2*quat.q0 + q3*quat.q1 - q1*quat.q3;
    qtmp.q3 = q0*quat.q3 + q3*quat.q0 + q1*quat.q2 - q2*quat.q1;
    return qtmp;
}

/**
 * @brief Quaternion multiply and assign itself
 */
CQuat& CQuat::operator*=(const CQuat &quat)
{
    return (*this=*this*quat);
}

/**
 * @brief Quaternion minus the rotation angle and assign itself
 */
CQuat& CQuat::operator-=(const CVect3 &phi)
{
    CQuat qtmp = rv2q(phi);
    return (*this=qtmp*(*this));
}

/**
 * @brief Quaternion inverse(or conjugate) operator
 */
CQuat operator~(const CQuat &q)
{
    return CQuat(q.q0,-q.q1,-q.q2,-q.q3);
}

/**
 * @brief Quaternion multiply rotation vector
 * @param[in] v Rotation vector
 * @return rotation vector
 */
CVect3 CQuat::operator*(const CVect3 &v) const
{
    CQuat qtmp;
    CVect3 vtmp;
    qtmp.q0 =         - q1*v.i - q2*v.j - q3*v.k;
    qtmp.q1 = q0*v.i           + q2*v.k - q3*v.j;
    qtmp.q2 = q0*v.j           + q3*v.i - q1*v.k;
    qtmp.q3 = q0*v.k           + q1*v.j - q2*v.i;
    vtmp.i = -qtmp.q0*q1 + qtmp.q1*q0 - qtmp.q2*q3 + qtmp.q3*q2;
    vtmp.j = -qtmp.q0*q2 + qtmp.q2*q0 - qtmp.q3*q1 + qtmp.q1*q3;
    vtmp.k = -qtmp.q0*q3 + qtmp.q3*q0 - qtmp.q1*q2 + qtmp.q2*q1;
    return vtmp;
}

/**
 * @brief Set yaw angle
 * @param[in] yaw input yaw angle
 */
void CQuat::SetYaw(double yaw)
{
    CVect3 att = q2att(*this);
    att.k = yaw;
    *this = a2qua(att);
}

/**
 * @brief Quaternion normalization
 * @param[in,out] q Input and output quaternion
 */
void normlize(CQuat *q)
{
    double nq=sqrt(q->q0*q->q0 + q->q1*q->q1 + q->q2*q->q2 + q->q3*q->q3);
    q->q0 /= nq, q->q1 /= nq, q->q2 /= nq, q->q3 /= nq;
}

/**
 * @brief Convert quaternion to rotation angel
 * @param[in] q Input quaternion
 * @return Rotation angel
 */
CVect3 q2rv(const CQuat &q)
{
    CQuat dq;
    dq = q;
    if(dq.q0<0)  { dq.q0=-dq.q0, dq.q1=-dq.q1, dq.q2=-dq.q2, dq.q3=-dq.q3; }
    if(dq.q0>1.0) dq.q0=1.0;
    double n2 = acos(dq.q0), f;
    if(n2>1.0e-20)
    {
        f = 2.0/(sin(n2)/n2);
    }
    else
    {
        f = 2.0;
    }
    return CVect3(dq.q1,dq.q2,dq.q3)*f;
}

/**
 * @brief Turn upside down the attitude, keep yaw unchanged
 * @param[in] q Input quaternion
 * @return Turned quaternion
 */
CQuat UpDown(const CQuat &q)
{
    CVect3 att = q2att(q);
    att.i = -att.i; att.j += PI;
    return a2qua(att);
}

/***************************  class CMat3  *********************************/
CMat3::CMat3(void)
{
}

/**
 * @brief Initalization 3D matrix by 9 number
 */
CMat3::CMat3(double xx, double xy, double xz, 
          double yx, double yy, double yz,
          double zx, double zy, double zz )
{
    e00=xx,e01=xy,e02=xz; e10=yx,e11=yy,e12=yz; e20=zx,e21=zy,e22=zz;
}

/**
 * @brief  Initalization 3D matrix by three row vector
 */
CMat3::CMat3(const CVect3 &v0, const CVect3 &v1, const CVect3 &v2)
{
    e00 = v0.i, e01 = v0.j, e02 = v0.k;
    e10 = v1.i, e11 = v1.j, e12 = v1.k;
    e20 = v2.i, e21 = v2.j, e22 = v2.k;
}

/**
 * @brief Determine attitude by double vector
 * @param[in] va1 Projection of vector "1" in the coordinate "a"
 * @param[in] va2 Projection of vector "2" in the coordinate "b"
 * @param[in] vb1 Projection of vector "1" in the coordinate "a"
 * @param[in] vb2 Projection of vector "2" in the coordinate "b"
 * @return Attitude a to b in DCM(Cab)
 */
CMat3 dv2att(const CVect3 &va1, const CVect3 &va2, CVect3 &vb1, const CVect3 &vb2)
{
    CVect3 a=va1*va2, b=vb1*vb2, aa=a*va1, bb=b*vb1;
    CMat3 Ma(va1/norm(va1),a/norm(a),aa/norm(aa));
    CMat3 Mb(vb1/norm(vb1),b/norm(b),bb/norm(bb));
    return (~Ma)*(Mb);
}

/**
 * @brief Negative full matrix
 */
CMat3 operator-(const CMat3 &m)
{
    return CMat3(-m.e00,-m.e01,-m.e02,-m.e10,-m.e11,-m.e12,-m.e20,-m.e21,-m.e22);
}

/**
 * @brief Matrix transpose
 */
CMat3 operator~(const CMat3 &m)
{
    return CMat3(m.e00,m.e10,m.e20, m.e01,m.e11,m.e21, m.e02,m.e12,m.e22);
}

/**
 * @brief Matrix multiplication
 */
CMat3 CMat3::operator*(const CMat3 &mat) const
{
    CMat3 mtmp;
    mtmp.e00 = e00*mat.e00 + e01*mat.e10 + e02*mat.e20;
    mtmp.e01 = e00*mat.e01 + e01*mat.e11 + e02*mat.e21;
    mtmp.e02 = e00*mat.e02 + e01*mat.e12 + e02*mat.e22;
    mtmp.e10 = e10*mat.e00 + e11*mat.e10 + e12*mat.e20;
    mtmp.e11 = e10*mat.e01 + e11*mat.e11 + e12*mat.e21;
    mtmp.e12 = e10*mat.e02 + e11*mat.e12 + e12*mat.e22;
    mtmp.e20 = e20*mat.e00 + e21*mat.e10 + e22*mat.e20;
    mtmp.e21 = e20*mat.e01 + e21*mat.e11 + e22*mat.e21;
    mtmp.e22 = e20*mat.e02 + e21*mat.e12 + e22*mat.e22;
    return mtmp;
}

/**
 * @brief Matrix addition
 */
CMat3 CMat3::operator+(const CMat3 &mat) const
{
    CMat3 mtmp;
    mtmp.e00 = e00 + mat.e00;  mtmp.e01 = e01 + mat.e01;  mtmp.e02 = e02 + mat.e02;  
    mtmp.e10 = e10 + mat.e10;  mtmp.e11 = e11 + mat.e11;  mtmp.e12 = e12 + mat.e12;  
    mtmp.e20 = e20 + mat.e20;  mtmp.e21 = e21 + mat.e21;  mtmp.e22 = e22 + mat.e22;  
    return mtmp;
}

/**
 * @brief Matrix subtraction
 */
CMat3 CMat3::operator-(const CMat3 &mat) const
{
    CMat3 mtmp;
    mtmp.e00 = e00 - mat.e00; mtmp.e01 = e01 - mat.e01; mtmp.e02 = e02 - mat.e02;
    mtmp.e10 = e10 - mat.e10; mtmp.e11 = e11 - mat.e11; mtmp.e12 = e12 - mat.e12;
    mtmp.e20 = e20 - mat.e20; mtmp.e21 = e21 - mat.e21; mtmp.e22 = e22 - mat.e22;
    return mtmp;
}

/**
 * @brief Matrix dot product(matrix times number)
 */
CMat3 CMat3::operator*(double f) const
{
    return CMat3(e00*f,e01*f,e02*f, e10*f,e11*f,e12*f, e21*f,e20*f,e22*f);
}

/**
 * @brief Matrix dot product(number times matrix)
 */
CMat3 operator*(double f, const CMat3 &m)
{
    return CMat3(m.e00*f,m.e01*f,m.e02*f, m.e10*f,m.e11*f,m.e12*f, 
            m.e20*f,m.e21*f,m.e22*f);
}

/**
 * @brief Matrix multiplied by column matrix
 * @param[in] v 3D column vector 
 * @return 3D column vector
 */
CVect3 CMat3::operator*(const CVect3 &v) const
{
    return CVect3(e00*v.i + e01*v.j + e02*v.k, e10*v.i + e11*v.j + e12*v.k,
            e20*v.i + e21*v.j + e22*v.k);
}

/**
 * @brief Determinant of matrix
 * @param[in] m Input matrix
 * @return Determinant(|m|)
 */
double det(const CMat3 &m)
{
    return m.e00*(m.e11*m.e22-m.e12*m.e21) - m.e01*(m.e10*m.e22-m.e12*m.e20) +
        m.e02*(m.e10*m.e21-m.e11*m.e20);
}

/**
 * @brief Convert Euler angle to quaternion
 * @param[in] pitch 
 * @param[in] roll
 * @param[in] yaw
 * @return Quaternion
 */
CQuat a2qua(double pitch, double roll, double yaw)
{
    pitch /= 2.0, roll /= 2.0, yaw /= 2.0;
    double  sp = sin(pitch), sr = sin(roll), sy = sin(yaw), 
            cp = cos(pitch), cr = cos(roll), cy = cos(yaw);
    CQuat qnb;
    qnb.q0 = cp*cr*cy - sp*sr*sy;
    qnb.q1 = sp*cr*cy - cp*sr*sy;
    qnb.q2 = cp*sr*cy + sp*cr*sy;
    qnb.q3 = cp*cr*sy + sp*sr*cy;
    return qnb;
}

/**
 * @brief Convert Euler angle to quaternion
 * @param[in] att Euler angle in 3D vector(pitch-roll-yaw)
 * @return Quaternion
 */
CQuat a2qua(const CVect3 &att)
{
    return a2qua(att.i, att.j, att.k);
}

/**
 * @brief Convert Euler angle to DCM
 * @param[in] att Euler angle in 3D vector(pitch-roll-yaw)
 * @return DCM(Cnb)
 */
CMat3 a2mat(const CVect3 &att)
{
    double  si = sin(att.i), ci = cos(att.i),
            sj = sin(att.j), cj = cos(att.j),
            sk = sin(att.k), ck = cos(att.k);
    CMat3 Cnb;
    Cnb.e00 =  cj*ck - si*sj*sk;    Cnb.e01 =  -ci*sk;  Cnb.e02 = sj*ck + si*cj*sk;
    Cnb.e10 =  cj*sk + si*sj*ck;    Cnb.e11 =  ci*ck;   Cnb.e12 = sj*sk - si*cj*ck;
    Cnb.e20 = -ci*sj;               Cnb.e21 =  si;      Cnb.e22 = ci*cj;
    return Cnb;
}

/**
 * @brief Convert DCM to Euler angle
 * @param[in] Cnb Input DCM
 * @return Euler angle in 3D vector
 */
CVect3 m2att(const CMat3 &Cnb)
{
    CVect3 att;
    att.i = asinEx(Cnb.e21);
    att.j = atan2Ex(-Cnb.e20, Cnb.e22);
    att.k = atan2Ex(-Cnb.e01, Cnb.e11);
    return att;
}

/**
 * @brief Convert DCM to quaternion
 * @param[in] Cnb Input DCM
 * @return Quaternion
 */
CQuat m2qua(const CMat3 &Cnb)
{
    double q0, q1, q2, q3, qq4;
    if(Cnb.e00>=Cnb.e11+Cnb.e22)
    {
        q1 = 0.5*sqrt(1+Cnb.e00-Cnb.e11-Cnb.e22);  qq4 = 4*q1;
        q0 = (Cnb.e21-Cnb.e12)/qq4; q2 = (Cnb.e01+Cnb.e10)/qq4; q3 = (Cnb.e02+Cnb.e20)/qq4;
    }
    else if(Cnb.e11>=Cnb.e00+Cnb.e22)
    {
        q2 = 0.5*sqrt(1-Cnb.e00+Cnb.e11-Cnb.e22);  qq4 = 4*q2;
        q0 = (Cnb.e02-Cnb.e20)/qq4; q1 = (Cnb.e01+Cnb.e10)/qq4; q3 = (Cnb.e12+Cnb.e21)/qq4;
    }
    else if(Cnb.e22>=Cnb.e00+Cnb.e11)
    {
        q3 = 0.5*sqrt(1-Cnb.e00-Cnb.e11+Cnb.e22);  qq4 = 4*q3;
        q0 = (Cnb.e10-Cnb.e01)/qq4; q1 = (Cnb.e02+Cnb.e20)/qq4; q2 = (Cnb.e12+Cnb.e21)/qq4;
    }
    else
    {
        q0 = 0.5*sqrt(1+Cnb.e00+Cnb.e11+Cnb.e22);  qq4 = 4*q0;
        q1 = (Cnb.e21-Cnb.e12)/qq4; q2 = (Cnb.e02-Cnb.e20)/qq4; q3 = (Cnb.e10-Cnb.e01)/qq4;
    }
    double nq = sqrt(q0*q0+q1*q1+q2*q2+q3*q3);
    q0 /= nq; q1 /= nq; q2 /= nq; q3 /= nq;
    return CQuat(q0, q1, q2, q3);
}

/**
 * @brief Convert quaternion to Euler angle
 * @param[in] qnb Input quaternion
 * @return Euler angle in 3D vector
 */
CVect3 q2att(const CQuat &qnb)
{
    double  q11 = qnb.q0*qnb.q0, q12 = qnb.q0*qnb.q1, q13 = qnb.q0*qnb.q2, q14 = qnb.q0*qnb.q3, 
            q22 = qnb.q1*qnb.q1, q23 = qnb.q1*qnb.q2, q24 = qnb.q1*qnb.q3,     
            q33 = qnb.q2*qnb.q2, q34 = qnb.q2*qnb.q3,  
            q44 = qnb.q3*qnb.q3;
    CVect3 att;
    att.i = asinEx(2*(q34+q12));
    att.j = atan2Ex(-2*(q24-q13), q11-q22-q33+q44);
    att.k = atan2Ex(-2*(q23-q14), q11-q22+q33-q44);
    return att;
}

/**
 * @brief Convert quaternion to DCM
 * @param[in] qnb Input quaternion
 * @return DCM
 */
CMat3 q2mat(const CQuat &qnb)
{
    double  q11 = qnb.q0*qnb.q0, q12 = qnb.q0*qnb.q1, q13 = qnb.q0*qnb.q2, q14 = qnb.q0*qnb.q3, 
            q22 = qnb.q1*qnb.q1, q23 = qnb.q1*qnb.q2, q24 = qnb.q1*qnb.q3,     
            q33 = qnb.q2*qnb.q2, q34 = qnb.q2*qnb.q3,  
            q44 = qnb.q3*qnb.q3;
    CMat3 Cnb;
    Cnb.e00 = q11+q22-q33-q44,  Cnb.e01 = 2*(q23-q14),     Cnb.e02 = 2*(q24+q13),
    Cnb.e10 = 2*(q23+q14),      Cnb.e11 = q11-q22+q33-q44, Cnb.e12 = 2*(q34-q12),
    Cnb.e20 = 2*(q24-q13),      Cnb.e21 = 2*(q34+q12),     Cnb.e22 = q11-q22-q33+q44;
    return Cnb;
}

/**
 * @brief Matrix inverse
 * @param[in] m Input matrix
 * @return Inverse of matrix
 */
CMat3 inv(const CMat3 &m)
{
    double nm;
    nm = m.e00*(m.e11*m.e22-m.e12*m.e21) - m.e01*(m.e10*m.e22-m.e12*m.e20) +
        m.e02*(m.e10*m.e21-m.e11*m.e20);
    CMat3 mtmp;
    // if(!sign(nm)) printf("CMat3::inv Determinant of Marix near zero\n");
    mtmp.e00 =  (m.e11*m.e22-m.e12*m.e21)/nm;
    mtmp.e10 = -(m.e10*m.e22-m.e12*m.e20)/nm;
    mtmp.e20 =  (m.e10*m.e21-m.e11*m.e20)/nm;
    mtmp.e01 = -(m.e01*m.e22-m.e02*m.e21)/nm;
    mtmp.e11 =  (m.e00*m.e22-m.e02*m.e20)/nm;
    mtmp.e21 = -(m.e00*m.e21-m.e01*m.e20)/nm;
    mtmp.e02 =  (m.e01*m.e12-m.e02*m.e11)/nm;
    mtmp.e12 = -(m.e00*m.e12-m.e02*m.e10)/nm;
    mtmp.e22 =  (m.e00*m.e11-m.e01*m.e10)/nm;
    return mtmp;
}

/**
 * @brief Get the main diagonal elements of matrix
 * @param[in] m Input matrix
 * @return Main diagonal elements vector
 */
CVect3 diag(const CMat3 &m)
{
    return CVect3(m.e00, m.e11, m.e22);
}

/**
 * @brief Create a main diagonal matrix by vector, 
 * @param[in] v Input vector 
 * @return Main diagonal matrix(other elements are zero)
 */
CMat3 diag(const CVect3 &v)
{
    return CMat3(v.i,0,0, 0,v.j,0, 0,0,v.k);
}

/***************************  class CMat  *********************************/
CMat::CMat(void)
{
#ifdef MAT_COUNT_STATISTIC
    #pragma message("  MAT_COUNT_STATISTIC")
    if(iMax<++iCount) iMax = iCount;
#endif
}

/**
 * @brief Initialize the size the of CMat
 * @param[in] row0 Number of row
 * @param[in] clm0 Number of column
 */
CMat::CMat(int row0, int clm0)
{
#ifdef MAT_COUNT_STATISTIC
    if(iMax<++iCount) iMax = iCount;
#endif
    row=row0; clm=clm0; rc=row*clm;
}

/*
 * @brief Initialize CMat's size and elements(same)
 * @param[in] row0 Number of row
 * @param[in] clm0 Number of column
 * @param[in] f Input number
 */
CMat::CMat(int row0, int clm0, double f)
{
#ifdef MAT_COUNT_STATISTIC
    if(iMax<++iCount) iMax = iCount;
#endif
    row=row0; clm=clm0; rc=row*clm;
    for(double *pd=dd, *pEnd=&dd[rc]; pd<pEnd; pd++)  *pd = f;
}

/**
 * @brief Initialize CMat's size and elements by array
 * @param[in] row0 Number of row
 * @param[in] clm0 Number of column
 * @param[in] pf Input array
 */
CMat::CMat(int row0, int clm0, const double *pf)
{
#ifdef MAT_COUNT_STATISTIC
    if(iMax<++iCount) iMax = iCount;
#endif
    row=row0; clm=clm0; rc=row*clm;
    memcpy(dd, pf, rc*sizeof(double));
}

#ifdef MAT_COUNT_STATISTIC
int CMat::iCount=0, CMat::iMax=0;
CMat::~CMat(void)
{
    iCount--;
}
#endif

/**
 * @brief Matrix multiplication
 */
CMat CMat::operator*(const CMat &m0) const
{
#ifdef MAT_COUNT_STATISTIC
    ++iCount;
#endif
    assert(this->clm==m0.row);
    CMat mtmp(this->row,m0.clm);
    int m=this->row, k=this->clm, n=m0.clm;
    double *p=mtmp.dd; const double *p1i=this->dd, *p2=m0.dd;
    for(int i=0; i<m; i++,p1i+=k)
    {
        for(int j=0; j<n; j++)
        {
            double f=0.0; const double *p1is=p1i, *p1isEnd=&p1i[k], *p2sj=&p2[j];
            for(; p1is<p1isEnd; p1is++,p2sj+=n)
                f += (*p1is) * (*p2sj);
            *p++ = f;
        }
    }
    return mtmp;
}

/**
 * @brief Matrix multiply by column vector
 */
CVect CMat::operator*(const CVect &v) const
{
    assert(this->clm==v.row);
    CVect vtmp(this->row);
    double *p=vtmp.dd, *pEnd=&vtmp.dd[vtmp.row]; const double *p1ij=this->dd, *p2End=&v.dd[v.row];
    for(; p<pEnd; p++)
    {
        double f=0.0; const double *p2j=v.dd;
        for(; p2j<p2End; p1ij++,p2j++)  f += (*p1ij) * (*p2j);
        *p = f;
    }
    return vtmp;
}

/**
 * @brief Matrix addition
 */
CMat CMat::operator+(const CMat &m0) const
{
#ifdef MAT_COUNT_STATISTIC
    ++iCount;
#endif
    assert(row==m0.row&&clm==m0.clm);
    CMat mtmp(row,clm);
    double *p=mtmp.dd, *pEnd=&mtmp.dd[rc]; const double *p1=this->dd, *p2=m0.dd;
    while(p<pEnd)
    { *p++ = (*p1++) + (*p2++); } 
    return mtmp;
}

/**
 * @brief Matrix add main diagonal vector and assign
 */
CMat& CMat::operator+=(const CVect &v)
{
    assert(row==v.row||clm==v.clm);
    int row1 = row+1;
    double *p=dd, *pEnd=&dd[rc];
    for(const double *p1=v.dd; p<pEnd; p+=row1, p1++)   *p += *p1;
    return *this;
}

/**
 * @brief Matrix Subtraction
 */
CMat CMat::operator-(const CMat &m0) const
{
#ifdef MAT_COUNT_STATISTIC
    ++iCount;
#endif
    assert(row==m0.row&&clm==m0.clm);
    CMat mtmp(row,clm);
    double *p=mtmp.dd, *pEnd=&mtmp.dd[rc]; const double *p1=this->dd, *p2=m0.dd;
    while(p<pEnd)
    { *p++ = (*p1++) - (*p2++); } 
    return mtmp;
}

/**
 * @brief Matrix dot product
 */
CMat CMat::operator*(double f) const
{
#ifdef MAT_COUNT_STATISTIC
    ++iCount;
#endif
    CMat mtmp(row,clm);
    double *p=mtmp.dd, *pEnd=&mtmp.dd[rc]; const double *p1=this->dd;
    while(p<pEnd)
    { *p++ = (*p1++) * f; } 
    return mtmp;
}

/**
 * @brief Matrix add and assign
 */
CMat& CMat::operator+=(const CMat &m0)
{
    assert(row==m0.row&&clm==m0.clm);
    double *p=dd, *pEnd=&dd[rc]; const double *p1=m0.dd;
    while(p<pEnd)
    { *p++ += *p1++; } 
    return *this;
}

/**
 * @brief Matrix minus and assign
 */
CMat& CMat::operator-=(const CMat &m0)
{
    assert(row==m0.row&&clm==m0.clm);
    double *p=dd, *pEnd=&dd[rc]; const double *p1=m0.dd;
    while(p<pEnd)
    { *p++ -= *p1++; } 
    return *this;
}

/**
 * @brief  Matrix dot product and assign
 */
CMat& CMat::operator*=(double f)
{
    double *p=dd, *pEnd=&dd[rc];
    while(p<pEnd)
    { *p++ *= f; } 
    return *this;
}

/**
 * @brief Add 1.0 to matrix diagonal elements
 */
CMat& CMat::operator++()
{
    int row1=row+1;
    for(double *p=dd, *pEnd=&dd[rc]; p<pEnd; p+=row1)   *p += 1.0;
    return *this;
}

/**
 * @brief Matrix transpose
 */
CMat operator~(const CMat &m0)
{
#ifdef MAT_COUNT_STATISTIC
    ++CMat::iCount;
#endif
    CMat mtmp(m0.clm,m0.row);
    const double *pm=m0.dd;
    for(int i=0; i<m0.row; i++)
        for(int j=i; j<m0.rc; j+=m0.row) 
            mtmp.dd[j] = *pm++;
    return mtmp;
}

/**
 * @brief Matrix symmetrization
 */
void symmetry(CMat &m)
{
    assert(m.row==m.clm);
    for(int i=0; i<m.clm; i++)
    {
        double *prow    = &m.dd[i*m.clm+i+1];
        double *prowEnd = &m.dd[i*m.clm+m.clm];
        double *pclm    = &m.dd[i*m.clm+i+m.clm];
        for(; prow<prowEnd; prow++,pclm+=m.clm)  *prow=*pclm=(*prow+*pclm)*0.5;
    }
}

/**
 * @brief Get Matrix element
 * @param[in] r row
 * @param[in] c column
 * @return the matrix element locate at (r,c)
 */
double& CMat::operator()(int r, int c)
{
    return this->dd[r*this->clm+c];
}

/**
 * @brief Set matrix specified row elements by a set of numbers
 * @param[in] i Specified row number
 * @param[in] f First number of the input number set
 * @param[in] ... Rest of numbers of the input number set
 */
void CMat::SetRow(int i, double f, ...)
{
    va_list vl;
    va_start(vl, f);
    for(double *p=&dd[i*clm], *pEnd=p+clm; p<pEnd; p++)
    { *p = f;  f = va_arg(vl, double);  }
    va_end(vl);
    return;
}

/**
 * @brief Set matrix specified row elements by a row vector
 * @param[in] i Specified row number
 * @param[in] v Input vector
 */
void CMat::SetRow(int i, const CVect &v)
{
    assert(clm==v.clm);
    const double *p=v.dd;
    for(double *p1=&dd[i*clm],*pEnd=p1+clm; p1<pEnd; p++,p1++) *p1 = *p;
    return;
}

/**
 * @brief Set matrix specified column elements by a set of numbers
 * @param[in] j Specified column number
 * @param[in] f First number of the input number set
 * @param[in] ... Rest of number of the input number set
 */
void CMat::SetClm(int j, double f, ...)
{
    va_list vl;
    va_start(vl, f);
    for(double *p=&dd[j], *pEnd=&p[rc]; p<pEnd; p+=clm)
    { *p = f;  f = va_arg(vl, double);  }
    va_end(vl);
    return;
}

/**
 * @brief Set matrix specified column elements by a column vector
 * @param[in] j Specified column number
 * @param[in] v Input vector
 */
void CMat::SetClm(int j, const CVect &v)
{
    assert(row==v.row);
    const double *p=v.dd;
    for(double *p1=&dd[j],*pEnd=&dd[rc]; p1<pEnd; p++,p1+=clm) *p1 = *p;
    return;
}

/**
 * @brief Set part of matrix specified column by a 3D vector
 * @param[in] i Specified row number
 * @param[in] j Specified column number
 * @param[in] v Input vector
 */
void CMat::SetClmVect3(int i, int j, const CVect3 &v)
{
    double *p=&dd[i*clm+j];
    *p = v.i; p += clm;
    *p = v.j; p += clm;
    *p = v.k;
}

/**
 * @brief Set part of matrix specified row by a 3D vector
 * @param[in] i Specified row number
 * @param[in] j Specified column number
 * @param[in] v Input vector
 */
void CMat::SetRowVect3(int i, int j, const CVect3 &v)
{
    *(CVect3*)&dd[i*clm+j] = v;
}

/**
 * @brief Set part of matrix in diagonal by a 3D vector
 * @param[in] i Specified row number
 * @param[in] j Specified column number
 * @param[in] v Input vector
 */
void CMat::SetDiagVect3(int i, int j, const CVect3 &v)
{
    double *p=&dd[i*clm+j];
    *p = v.i;  p += clm+1;
    *p = v.j;  p += clm+1;
    *p = v.k;
}

/**
 * @brief Set part of matrix by a 3x3 matrix
 * @param[in] i Specified row number
 * @param[in] j Specified column number
 * @param[in] m Input 3x3 matrix
 */
void CMat::SetMat3(int i, int j, const CMat3 &m)
{
    double *p=&dd[i*clm+j];
    *(CVect3*)p = *(CVect3*)&m.e00;  p += clm;
    *(CVect3*)p = *(CVect3*)&m.e10;  p += clm;
    *(CVect3*)p = *(CVect3*)&m.e20;
}

/**
 * @brief Get a 3x3 matrix from whole matrix
 * @param[in] i Specified row number
 * @param[in] j Specified column number
 * @return 3x3 matrix
 */
CMat3 CMat::GetMat3(int i, int j) const
{
    CMat3 m;
    const double *p=&dd[i*clm+j];
    *(CVect3*)&m.e00 = *(CVect3*)p;  p += clm;
    *(CVect3*)&m.e10 = *(CVect3*)p;  p += clm;
    *(CVect3*)&m.e20 = *(CVect3*)p;
    return m;
}

/**
 * @brief Part of matrix add 3x3 matrix
 * @param[in] i Specified row number
 * @param[in] j Specified column number
 * @param[in] m Input 3x3 matrix
 */
void CMat::SubAddMat3(int i, int j, const CMat3 &m)
{
    double *p=&dd[i*clm+j];
    *(CVect3*)p += *(CVect3*)&m.e00;  p += clm;
    *(CVect3*)p += *(CVect3*)&m.e10;  p += clm;
    *(CVect3*)p += *(CVect3*)&m.e20;
}

/**
 * @brief Get row vector from matrix
 * @param[in] i Specified row number
 * @return Row vector
 */
CVect CMat::GetRow(int i) const
{
    CVect v(1, clm);
    const double *p1=&dd[i*clm], *pEnd=p1+clm;
    for(double *p=v.dd; p1<pEnd; p++,p1++) *p = *p1;
    return v;
}

/**
 * @brief Get column vector from matrix
 * @param[in] Specified column number
 * @return Column matrix 
 */
CVect CMat::GetClm(int j) const
{
    CVect v(row, 1);
    const double *p1=&dd[j], *pEnd=&dd[rc];
    for(double *p=v.dd; p1<pEnd; p++,p1+=clm) *p = *p1;
    return v;
}

/**
 * @brief Set matrix row to 0
 * @param[in] i Specified column number
 */
void CMat::ZeroRow(int i)
{
    for(double *p=&dd[i*clm],*pEnd=p+clm; p<pEnd; p++) *p = 0.0;
    return;
}

/**
 * @brief Set matrix column to 0
 * @param[in] j Specified row number
 */
void CMat::ZeroClm(int j)
{
    for(double *p=&dd[j],*pEnd=&dd[rc]; p<pEnd; p+=clm) *p = 0.0;
    return;
}

/**
 * @brief Set matrix main diagonal elements by a set of number
 * @param[in] f  First number
 * @param[in] ... The rest of number
 */
void CMat::SetDiag(double f, ...)
{
    *this = CMat(this->row, this->clm, 0.0);
    va_list vl;
    va_start(vl, f);
    double *p=dd, *pEnd=&dd[rc];
    for(int row1=row+1; p<pEnd; p+=row1)
    { 
        *p = f;
        f = va_arg(vl, double);
    }
    va_end(vl);
}

/**
 * @brief Set Matrix main diagonal elements by Squared numbers(f^2)
 * @param[in] f First number
 * @param[in] ... The rest of number
 */
void CMat::SetDiag2(double f, ...)
{
    *this = CMat(this->row, this->clm, 0.0);
    va_list vl;
    va_start(vl, f);
    double *p=dd, *pEnd=&dd[rc];
    for(int row1=row+1; p<pEnd; p+=row1)
    { 
        *p = f*f; 
        f = va_arg(vl, double);
    }
    va_end(vl);
}

/**
 * @brief 1-norm of matrix, Calculate the sum of absolute number of every matrix
 *      element
 * @param[in] m Input matrix
 * @return 1-norm of the matrix
 */
double norm1(const CMat &m)
{
    double n1=0.0;
    for(const double *p=m.dd,*pEnd=&m.dd[m.rc]; p<pEnd; p++)
    {
        if(*p>0.0)   n1 += *p;
        else  n1 -= *p;
    }
    return n1;
}

/**
 * @brief Get diagonal elements of matrix to a column vector
 * @param[in] m Input Matrix
 * @return Row vector consist of matrix main diagonal elements
 */
CVect diag(const CMat &m)
{
    int row1 = m.row+1;
    CVect vtmp(m.row,1);
    double *p=vtmp.dd, *pEnd=&vtmp.dd[vtmp.row];
    for(const double *p1=m.dd; p<pEnd; p++, p1+=row1)   *p = *p1;
    return vtmp;
}

/**
 * @brief  Calculate the vector(r row of matrix m0) multiply by matrix m1:
 *      m(r,:)=m0(r,:)*m1
 * @param[out] m Matrix result
 * @param[in] m0 Input matrix providing row vector
 * @param[in] m1 By matrix
 * @param[in] r  The row of m0
 */
void RowMul(CMat &m, const CMat &m0, const CMat &m1, int r)
{
    assert(m0.clm==m1.row);
    int rc0=r*m0.clm;
    double *p=&m.dd[rc0], *pEnd=p+m0.clm;
    const double *p0=&m0.dd[rc0], *p0End=p0+m0.clm, *p1j=m1.dd;
    for(; p<pEnd; p++)
    {
        double f=0.0; const double *p0j=p0, *p1jk=p1j++;
        for(; p0j<p0End; p0j++,p1jk+=m1.clm)     f += (*p0j) * (*p1jk);
        *p = f;
    }
}

/**
 * @brief  Calculate the vector(r row of matrix m0) multiply by transpose of 
 *      matrix m1: m(r,:)=m0(r,:)*m1'
 * @param[out] m Matrix result
 * @param[in] m0 Input matrix providing row vector
 * @param[in] m1 By matrix
 * @param[in] r  The row of m0
 */
void RowMulT(CMat &m, const CMat &m0, const CMat &m1, int r)
{
    assert(m0.clm==m1.clm);
    int rc0=r*m0.clm;
    double *p=&m.dd[rc0], *pEnd=p+m0.clm;
    const double *p0=&m0.dd[rc0], *p0End=p0+m0.clm, *p1jk=m1.dd;
    for(; p<pEnd; p++)
    {
        double f=0.0; const double *p0j=p0;
        for(; p0j<p0End; p0j++,p1jk++)   f += (*p0j) * (*p1jk);
        *p = f;
    }
}

/**
 * @brief Create main diagonal matrix by a vector
 * @param[in] v Input vector
 * @return Mian diagonal matrix
 */
CMat diag(const CVect &v)
{
#ifdef MAT_COUNT_STATISTIC
    ++CMat::iCount;
#endif
    int rc = v.row>v.clm ? v.row : v.clm, rc1=rc+1;
    CMat mtmp(rc,rc,0.0);
    double *p=mtmp.dd;
    for(const double *p1=v.dd, *p1End=&v.dd[rc]; p1<p1End; p+=rc1, p1++)
        *p = *p1;
    return mtmp;
}
/**
* @brief M = diag(V)*M*diag(V)*afa
* @param[in] V Input Vector
* @param[in] M Input Matrix
* @param[in] afa Input factor
*/
void DVMDVafa(const CVect &V, CMat &M, double afa)
{
    assert(V.rc==M.row&&M.row==M.clm);
    int i = 0;
    const double *pv = V.dd;
    for(double vi=*pv, viafa=vi*afa; i<M.clm; i++,pv++,vi=*pv,viafa=vi*afa)
    {
        double *prow=&M.dd[i*M.clm],*prowEnd=prow+M.clm,*pclm=&M.dd[i];
        for(;prow<prowEnd; prow++,pclm+=M.row)
        {
            *prow *= vi;
            *pclm *= viafa;
        }
    }
}

/***************************  class CVect  *********************************/
CVect::CVect(void)
{
}

/**
 * @brief Specified the vector Size, either row0 or clm0 is 1
 * @param[in] row0 Number of row
 * @param[in] clm0 Number of column
 */
CVect::CVect(int row0, int clm0)
{
    if(clm0==1) { row=row0; clm=1;   }
    else        { row=1;    clm=clm0;}
    rc = row*clm;
 }

/**
 * @brief Initialize column vector by same number
 * @param[in] row0 Number of row
 * @param[in] f Specified number
 */
CVect::CVect(int row0, double f)
{
    row=row0; clm=1; rc=row*clm;
    for(int i=0;i<row;i++) dd[i]=f;
}

/**
 * @brief Initialize column vector by an array
 * @param[in] row0  Number of row
 * @param[in] pf    Input arrary
 */
CVect::CVect(int row0, const double *pf)
{
    row=row0; clm=1; rc=row*clm;
    memcpy(dd, pf, row*sizeof(double));
}

/**
 * @brief Initialize column vector by a set of number(at least two)
 * @param[in] row0 Number of row
 * @param[in] f  First number
 * @param[in] ... The rest of number
 */
CVect::CVect(int row0, double f, double f1, ...)
{
    row=row0; clm=1; rc=row*clm;
    assert(row<=MMD&&clm<=MMD);
    va_list vl;
    va_start(vl, f);
    for(int i=0, rc=row>clm?row:clm; i<rc; i++)
    { dd[i] = f;  f = va_arg(vl, double);   }
    va_end(vl);
}

/**
 * @brief Initialize column vector by a 3D vector
 * @param[in] v Input 3D vector
 */
CVect::CVect(const CVect3 &v)
{
    row=3; clm=1; rc=row*clm;
    dd[0]=v.i; dd[1]=v.j; dd[2]=v.k;
}

/**
 * @brief Initialize 6x1 vector by two 3D vectors
 * @param[in] v1 First 3D vector
 * @param[in] v2 Second 3D vector
 */
CVect::CVect(const CVect3 &v1, const CVect3 v2)
{
    row=6; clm=1; rc=row*clm;
    dd[0]=v1.i; dd[1]=v1.j; dd[2]=v1.k;
    dd[3]=v2.i; dd[4]=v2.j; dd[5]=v2.k;
}

/**
 * @brief Transpose of vector
 */
CVect operator~(const CVect &v)
{
    CVect vtmp=v;
    vtmp.row=v.clm; vtmp.clm=v.row;
    return vtmp;
}

/**
 * @brief Row vector multiply by matrix 
 */
CVect CVect::operator*(const CMat &m) const
{
    assert(clm==m.row);
    CVect vtmp(row,clm);
    double *p=vtmp.dd; const double *p1End=&dd[clm];
    for(int j=0; j<clm; p++,j++)
    {
        double f=0.0; const double *p1j=dd, *p2jk=&m.dd[j];
        for(; p1j<p1End; p1j++,p2jk+=m.clm)  f += (*p1j) * (*p2jk);
        *p = f;
    }
    return vtmp;
}

/**
 * @brief Two vector production: (1xn)*(nx1) or (nx1)*(1*n)
 * @param[in] v Multiplied matrix
 * @return (n*n) matrix or (1x1) matrix
 */
CMat CVect::operator*(const CVect &v) const
{
#ifdef MAT_STATISTIC
    ++CMat::iCount;
#endif
    assert(clm==v.row);
    CMat mtmp(row,v.clm);
    if(row==1 && v.clm==1)  // (1x1) = (1xn)*(nx1)
    {
        double f = 0.0;
        for(int i=0; i<clm; i++)  f += dd[i]*v.dd[i];
        mtmp.dd[0] = f;
    }
    else    // (nxn) = (nx1)*(1xn)
    {
        double *p = mtmp.dd;
        for(int i=0; i<row; i++)
        {
            for(int j=0; j<v.clm; j++)  *p++ = dd[i]*v.dd[j];
        }
    }
    return mtmp;
}

/**
 * @brief Vector addition, adding the cooresponding elements.
 */
CVect CVect::operator+(const CVect &v) const
{
    assert(row==v.row&&clm==v.clm);
    const double *p2=v.dd, *p1=dd, *p1End=&dd[rc];
    CVect vtmp(row,clm);
    for(double *p=vtmp.dd; p1<p1End; p++,p1++,p2++)  { *p=*p1+*p2; }
    return vtmp;
}

/**
 * @brief Vector Substraction, Subtraction of cooresponding elements
 */
CVect CVect::operator-(const CVect &v) const
{
    assert(row==v.row&&clm==v.clm);
    const double *p2=v.dd, *p1=dd, *p1End=&dd[rc];
    CVect vtmp(row,clm);
    for(double *p=vtmp.dd; p1<p1End; p++,p1++,p2++)  { *p=*p1-*p2; }
    return vtmp;
}

/**
 * @brief Vector dot product
 */
CVect CVect::operator*(double f) const
{
    CVect vtmp(row,clm);
    const double *p1=dd,*p1End=&dd[rc];
    for(double *p=vtmp.dd; p1<p1End; p++,p1++)  { *p=*p1*f; }
    return vtmp;
}

/**
 * @brief Let all elements in the vector be the equal to the number f
 * @param[in] f Input number
 * @return Vector
 */
CVect& CVect::operator=(double f)
{
    for(double *p=dd, *pEnd=&dd[rc]; p<pEnd; p++)  { *p = f; }
    return *this;
}

/**
 * @brief Assign an array to a vector
 * @param[in] pf Input array(lenth equal to vector)
 * @return Vector
 */
CVect& CVect::operator=(const double *pf)
{
    for(double *p=dd, *pEnd=&dd[rc]; p<pEnd; p++,pf++)  { *p = *pf; }
    return *this;
}

/**
 * @brief Add and assign a vector
 */
CVect& CVect::operator+=(const CVect &v)
{
    assert(row==v.row&&clm==v.clm);
    const double *p1 = v.dd;
    for(double *p=dd, *pEnd=&dd[rc]; p<pEnd; p++,p1++)  { *p += *p1; }
    return *this;
}

/**
 * @brief Substract and assign a vector
 */
CVect& CVect::operator-=(const CVect &v)
{
    assert(row==v.row&&clm==v.clm);
    const double *p1 = v.dd;
    for(double *p=dd, *pEnd=&dd[rc]; p<pEnd; p++,p1++)  { *p -= *p1; }
    return *this;
}

/**
 * @brief Vector dot product and assign
 */
CVect& CVect::operator*=(double f)
{
    for(double *p=dd, *pEnd=&dd[rc]; p<pEnd; p++)  { *p *= f; }
    return *this;
}

/**
 * @brief The power of each element in the vector
 * @param[in] v Input vector
 * @param[in] k Order of power
 * @return Vector
 */
CVect pow(const CVect &v, int k)
{
    CVect pp = v;
    double *p, *pEnd=&pp.dd[pp.rc];
    for(int i=1; i<k; i++)
    {
        p=pp.dd;
        for(const double *p1=v.dd; p<pEnd; p++,p1++)
            *p *= *p1;
    }
    return pp;
}

/**
 * @brief The absolute value of each element in the vector
 */
CVect abs(const CVect &v)
{
    CVect res(v.row,v.clm);
    const double *p=v.dd, *pEnd=&v.dd[v.rc];
    for(double *p1=res.dd; p<pEnd; p++,p1++)  { *p1 = *p>0 ? *p : -*p; }
    return res;
}

/**
 * @brief 2-norm of vector(Euclidean Distance)
 */
double norm(const CVect &v)
{
    const double *p=v.dd, *pEnd=&v.dd[v.rc];
    double f=0.0;
    for(; p<pEnd; p++)  { f += (*p)*(*p); }
    return sqrt(f);
}

/**
 * @brief  Get vector element value
 * @param[in] r Position
 * @return Element r
 */
double& CVect::operator()(int r)
{
    return this->dd[r];
}

/**
 * @brief Set vector by a set of number
 */
void CVect::Set(double f, ...)
{
    assert(rc<=MMD);
    va_list vl;
    va_start(vl, f);
    for(int i=0; i<rc; i++)
    { dd[i] = f;  f = va_arg(vl, double);   }
    va_end(vl);
}

/**
 * @brief Set vector elements by a set of squared number
 */
void CVect::Set2(double f, ...)
{
    assert(rc<=MMD);
    va_list vl;
    va_start(vl, f);
    for(int i=0; i<rc; i++)
    { dd[i] = f*f;  f = va_arg(vl, double); }
    va_end(vl);
}

/***************************  class CEarth  *********************************/
/**
 * @brief Initialize CEarth class
 * @param[in] a0 Earth equatorial radius
 * @param[in] f0 Earth flattening
 * @param[in] wie Earth rotational angular velocity 
 * @param[in] g0 Earth gravity magnitude
 */
CEarth::CEarth(double a0, double f0, double wie, double g0)
{
    a = a0; f = f0; 
    this->wie = wie; 
    this->g0 = g0;
    b = (1-f)*a;
    e = sqrt(a*a-b*b)/a;    e2 = e*e;
    gn = O31;  pgn = 0;
    Update(O31);
}

/**
 * @brief Update mamber variables
 * @param[in] pos Geodetic coordinate postion
 * @param[in] vn  Velocity: \f$ v_{en}^n \f$
 */
void CEarth::Update(const CVect3 &pos, const CVect3 &vn)
{
#ifdef PSINS_LOW_GRADE_MEMS
    this->pos = pos;  this->vn = vn;
    sl = sin(pos.i), cl = cos(pos.i), tl = sl/cl;
    double sq = 1-e2*sl*sl, sq2 = sqrt(sq);
    RMh = a*(1-e2)/sq/sq2+pos.k;    f_RMh = 1.0/RMh;
    RNh = a/sq2+pos.k;    clRNh = cl*RNh;  f_RNh = 1.0/RNh; f_clRNh = 1.0/clRNh;
//  wnie.i = 0.0,           wnie.j = wie*cl,        wnie.k = wie*sl;
//  wnen.i = -vn.j*f_RMh,   wnen.j = vn.i*f_RNh,    wnen.k = wnen.j*tl;
    wnin = wnie = wnen = O31;
    sl2 = sl*sl;
    gn.k = -( g0*(1+5.27094e-3*sl2)-3.086e-6*pos.k );
    gcc = pgn ? *pgn : gn;
#else
    this->pos = pos;  this->vn = vn;
    sl = sin(pos.i), cl = cos(pos.i), tl = sl/cl;
    double sq = 1-e2*sl*sl, sq2 = sqrt(sq);
    RMh = a*(1-e2)/sq/sq2+pos.k;    f_RMh = 1.0/RMh;
    RNh = a/sq2+pos.k;    clRNh = cl*RNh;  f_RNh = 1.0/RNh; f_clRNh = 1.0/clRNh;
    // ref Yan2016(P70,4.1-3)
    wnie.i = 0.0,           wnie.j = wie*cl,        wnie.k = wie*sl;
    // ref Yan2016(P70,4.1-4)
    wnen.i = -vn.j*f_RMh,   wnen.j = vn.i*f_RNh,    wnen.k = wnen.j*tl;
    wnin = wnie + wnen;
    sl2 = sl*sl, sl4 = sl2*sl2;
    // ref Yan2016(P49,3.2-23)
    gn.k = -( g0*(1+5.27094e-3*sl2+2.32718e-5*sl4)-3.086e-6*pos.k );
    gcc = pgn ? *pgn : gn;
    // ref Yan2016(P71,4.1-20)
    gcc -= (wnie+wnin)*vn;
#endif
}

/**
 * @brief Calculate position change in a short time
 * @param[in] vn Velocity under n-frame
 * @param[in] ts Time span(sec)
 * @return Position Change under n-frame
 */
CVect3 CEarth::vn2dpos(const CVect3 &vn, double ts) const
{
    return CVect3(vn.j*f_RMh, vn.i*f_clRNh, vn.k)*ts;
}

/**
 * @brief The same as atan2, ref: https://en.wikipedia.org/wiki/Atan2
 * @param[in] y Input y
 * @param[in] x Input x
 * @return arctan(y/x)
 */
double atan2Ex(double y, double x)
{
    if((sign(y)==0) && (sign(x)==0))
        return  0.0;
    else
        return atan2(y, x);
}

/**
* @brief Determine the sign of 'val' with the sensitivity of 'eps'
* @param[in] val Input value
* @param[in] eps Zero threshold
* @return 
*   @retval -1 val is a negative number
*   @retval 0  val equalval to zero
*   @retval 1 val is a positive number
*/
int sign(double val, double eps)
{
    if(val<-eps) 
        return -1;
    else if(val>eps) 
        return  1;
    else 
        return 0;
}

/**
 * @brief Convert radian to deg_min notation(a unique method of unit, e.g.
 *      1234.56 deg_min = 12deg + 34.56min)
 * @param[in] r Input angle in radian
 * @return Angle in deg_min
 */
double r2dm(double r)
{
    int sgn=1;
    if(r<0.0) { r=-r; sgn=0; }
    double deg = r/DEG;
    int ideg = (int)deg;
    double dm = ideg*100 + (deg-ideg)*60.0;
    return sgn ? dm : -dm;
}

/**
 * @brief Convert deg_min notation(a unique method of unit, e.g.
 *      1234.56 deg_min = 12deg + 34.56min) to radian
 * @param[in] dm Input angle in deg_min
 * @return Angle in radian
 */
double dm2r(double dm)
{
    int sgn=1;
    if(dm<0.0) { dm=-dm; sgn=0; }
    int ideg = (int)(dm/100);
    double r = ideg*DEG + (dm-ideg*100)*(DEG/60);
    return sgn ? r : -r;
}

/**
 * @brief Return value 'val' between range 'minVal' and 'maxVal'
 * @param[in] val       Input value
 * @param[in] minVal    Upper bound
 * @param[in] maxVal    Lower bound
 * @return
 *      @retval val     minVal < val < maxVal
 *      @retval minVal  val < minVal
 *      @retval maxVal  val > maxVal
 */
double range(double val, double minVal, double maxVal)
{
    if(val<maxVal)
        return val>minVal?val:minVal;
    else
        return maxVal;
}

