/**
 * @file PSINS.cpp
 * @brief PSINS Main Definetaion File
 * @author Yan GongMing
 * @version 1.0
 * @date 2018-12-13
 */
#include "PSINS.h"

/***************************  class CIIR  *********************************/
CIIR::CIIR(void){}

CIIR::CIIR(double *b0, double *a0, int n0)
{
    psinsassert(n0<IIRnMax);
    for(int i=0; i<n0; i++)  { b[i]=b0[i]/a0[0]; a[i]=a0[i]; x[i]=y[i]=0.0; }
    n = n0;
}

/**
 * @brief 
 * @param[in] x0
 * @return 
 */
double CIIR::Update(double x0)
{
//  a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
//                        - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
    double y0 = 0.0;
    for(int i=n-1; i>0; i--)
    {
        x[i] = x[i-1]; y[i] = y[i-1];
        y0 += b[i]*x[i] - a[i]*y[i];
    }
    x[0] = x0;
    y0 += b[0]*x0;
    y[0] = y0;
    return y0;
}

/***************************  class CV3IIR  *********************************/
CIIRV3::CIIRV3(void)
{
}

CIIRV3::CIIRV3(double *b0, double *a0, int n0, double *b1, double *a1, int n1, 
               double *b2, double *a2, int n2)
{
    iir0 = CIIR(b0, a0, n0);
    if(n1==0)   iir1 = iir0;    // iir1 the same as iir0
    else        iir1 = CIIR(b1, a1, n1);
    if(n2==0)   iir2 = iir0;    // iir2 the same as iir0
    else        iir2 = CIIR(b2, a2, n2);
}

CVect3 CIIRV3::Update(CVect3 &x)
{
    return y = CVect3(iir0.Update(x.i), iir1.Update(x.j), iir2.Update(x.k));
}

/***************************  class CRAvar  *********************************/
CRAvar::CRAvar()
{
}

CRAvar::CRAvar(int nR0, int maxCount0)
{
    psinsassert(nR0<RAMAX);
    this->nR0 = nR0;
    for(int i=0; i<RAMAX; i++)  { Rmaxcount[i]=maxCount0, tau[i]=INF; }
}

void CRAvar::set(double r0, double tau, double rmax, double rmin, int i)
{
    this->R0[i] = r0*r0;
    this->tau[i] = tau;
    this->r0[i] = 0.0;  Rmaxflag[i] = Rmaxcount[i];
    this->Rmax[i] = rmax==0.0 ? 100.0*this->R0[i] : rmax*rmax;
    this->Rmin[i] = rmin==0.0 ?  0.01*this->R0[i] : rmin*rmin;
}

void CRAvar::set(const CVect3 &r0, const CVect3 &tau, const CVect3 &rmax, const CVect3 &rmin)
{
    const double *pr0=&r0.i, *ptau=&tau.i, *prmax=&rmax.i, *prmin=&rmin.i;
    for(int i=0; i<3; i++,pr0++,ptau++,prmax++,prmin++)
        set(*pr0, *ptau, *prmax, *prmin, i);
}

void CRAvar::set(const CVect &r0, const CVect &tau, const CVect &rmax, const CVect &rmin)
{
    const double *pr0=r0.dd, *ptau=tau.dd, *prmax=rmax.dd, *prmin=rmin.dd;
    for(int i=0; i<nR0; i++,pr0++,ptau++,prmax++,prmin++)
        set(*pr0, *ptau, *prmax, *prmin, i);
}

void CRAvar::Update(double r, double ts, int i)
{
    if(tau[i]>INF/2) return;
    double tstau = ts>tau[i] ? 1.0 : ts/tau[i];
    double dr2=r-r0[i]; dr2=dr2*dr2; r0[i]=r;
    if(dr2>R0[i]) R0[i]=dr2; else R0[i]=(1.0-tstau)*R0[i]+tstau*dr2;
    if(R0[i]<Rmin[i]) R0[i]=Rmin[i];
    if(R0[i]>Rmax[i]) {R0[i]=Rmax[i];Rmaxflag[i]=Rmaxcount[i];} else {Rmaxflag[i]-=Rmaxflag[i]>0;}
}

void CRAvar::Update(const CVect3 &r, double ts)
{
    const double *pr=&r.i;
    for(int i=0; i<3; i++,pr++)
        Update(*pr, ts, i);
}

void CRAvar::Update(const CVect &r, double ts)
{
    const double *pr=r.dd;
    for(int i=0; i<nR0; i++,pr++)
        Update(*pr, ts, i);
}

double CRAvar::operator()(int k)
{
    return Rmaxflag[k] ? INF : sqrt(R0[k]);
}

/***************************  class CVAR  ***********************************/
CVAR::CVAR(int imax0, float data0)
{
    imax = min(imax0, VARMAX);
    for(ipush=0; ipush<imax; ipush++)   array[ipush] = data0;
    ipush = 0;
    mean = data0; var = 0.0;
}

float CVAR::Update(float data)
{
    array[ipush] = data;
    if(++ipush==imax) ipush=0;
    float *p0, *p1;
    for(mean=0.0,p0=&array[0],p1=&array[imax]; p0<p1; p0++)  mean+=*p0;
    mean /= imax;
    for(var=0.0,p0=&array[0],p1=&array[imax]; p0<p1; p0++)  { float vi=*p0-mean; var+=vi*vi; }
    var /= imax-1;
    return var;
}

/***************************  class CVARn  ***********************************/
CVARn::CVARn(void)
{
    pData = NULL;
}

CVARn::CVARn(int row0, int clm0)
{
    row = row0, clm = clm0;
    pData = new double*[clm];
    pData[0] = new double[row*clm+5*clm];
    for (int i = 1; i < clm; i++) {
        pData[i] = pData[i - 1] + row;
    }
    pd = pData[clm-1]+row;  Sx = pd + clm;  Sx2 = Sx + clm;  mx = Sx2 + clm;  stdx = mx + clm;
    stdsf = sqrt((row - 1) / (row - 2.0));
    Reset();
}

CVARn::~CVARn(void)
{
    if (pData) {
//      delete pData[0];  ?
        delete pData; pData = NULL;
    }
}

void CVARn::Reset(void)
{
    idxpush = rowcnt = 0;
    memset(pData[0], 0, row*clm*sizeof(double));
    memset(Sx, 0, 5*clm*sizeof(double));
}

BOOL CVARn::Update(const double *pf)
{
    if (!pData[0])  return FALSE;
    if (++rowcnt > row) rowcnt = row;
    int idxnext = (idxpush >= row - 1) ? 0 : idxpush + 1;
    for (int i = 0; i < clm; i++)
    {
        double f=*pf++;  if(f>1e5) f=1e5; else if(f<-1e5) f=-1e5;
        pData[i][idxpush] = f;
        Sx[i] += f - pData[i][idxnext];
        Sx2[i] += f*f - pData[i][idxnext] * pData[i][idxnext];
        mx[i] = Sx[i] / rowcnt;
        stdx[i] = sqrt(Sx2[i] / rowcnt - mx[i] * mx[i]) * stdsf;   // Dx = E(x^2) - (Ex)^2
    }
    if (++idxpush == row) {
        idxpush = 0;
    }
    return idxpush == 0;
}

BOOL CVARn::Update(double f, ...)
{
    va_list vl;
    va_start(vl, f);
    for (int i = 0; i < clm; i++)
    {
        pd[i] = f;
        f = va_arg(vl, double);
    }
    va_end(vl);
    return Update(pd);
}

/***************************  class CKalman  *********************************/
/**
 * @brief 
 * @param nq0
 * @param nr0
 */
CKalman::CKalman(int nq0, int nr0)
{
    psinsassert(nq0<=MMD&&nr0<=MMD);
    kftk = 0.0;
    nq = nq0; nr = nr0;
    Ft = Pk = CMat(nq,nq,0.0);
    Hk = CMat(nr,nq,0.0);  Fading = CMat(nr,nq,1.0); zfdafa = 0.1;
    Qt = Pmin = Xk = CVect(nq,0.0);  Xmax = Pmax = CVect(nq,INF);
    Zk = CVect(nr,0.0);  Rt = CVect(nr,INF); rts = CVect(nr,1.0);  Zfd = CVect(nr,0.0); Zfd0 = CVect(nr,INF);
    Rmax = CVect(nr,INF); Rmin = Rb = CVect(nr,0.0); Rbeta = CVect(nr,1.0);
    for(int i=0; i<nr; i++) { Rmaxcount[i]=0, Rmaxcount0[i]=5; }
    FBTau = FBMax = CVect(nq,INF); FBXk = FBTotal = CVect(nq,0.0);
    measflag = 0;
}

void CKalman::TimeUpdate(double kfts, int fback)
{
/*  CMat Fk, FtPk;
    kftk += kfts;
    SetFt();
    Fk=++(Ft*kfts);  // Fk = I+Ft*ts
    Xk = Fk * Xk;
    FtPk = Ft*Pk*kfts;
    Pk = Pk + FtPk + (~FtPk);// + FtPk*(~Ft)*kfts;
    Pk += Qt*kfts;
    if(fback)  Feedback(kfts);
*/
    CMat Fk;
    kftk += kfts;
    SetFt();
    Fk=++(Ft*kfts);  // Fk = I+Ft*ts
    Xk = Fk * Xk;
    Pk = Fk*Pk*(~Fk);  Pk += Qt*kfts;
    if(fback)  Feedback(kfts);
}

void CKalman::SetMeasFlag(int flag)
{
    measflag = (flag==0) ? 0 : (measflag|flag);
}

int CKalman::MeasUpdate(double fading)
{
    CVect Pxz, Kk, Hi;
    SetMeas();
    for(int i=0; i<nr; i++)
    {
        if(measflag&(0x01<<i))
        {
            Hi = Hk.GetRow(i);
            Pxz = Pk*(~Hi);
            double Pz0 = (Hi*Pxz)(0,0), r=Zk(i)-(Hi*Xk)(0,0);
            if(Rb.dd[i]>EPS)
                RAdaptive(i, r, Pz0);
            if(Zfd.dd[i]<INF/2)
                RPkFading(i);
            double Pzz = Pz0+Rt.dd[i]/rts.dd[i];
            Kk = Pxz*(1.0/Pzz);
            Xk += Kk*r;
            Pk -= Kk*(~Pxz);
        }
    }
    if(fading>1.0) Pk *= fading;
    XPConstrain();
    symmetry(Pk);
    SetMeasFlag(0);
    return measflag;
}

/**
 * @brief Measurement noise adaptive(Sage-Husa adaptive filter)
 * @param[in] i 
 * @param[in] r 
 * @param[in] Pr
 * @return 
 *      @retval 0
 *      @retval 1
 */
int CKalman::RAdaptive(int i, double r, double Pr)
{
    double rr=r*r-Pr;
    if(rr<Rmin.dd[i])   rr = Rmin.dd[i];
//  if(rr>Rmax.dd[i])   Rt.dd[i] = Rmax.dd[i];  
//  else                Rt.dd[i] = (1.0-Rbeta.dd[i])*Rt.dd[i]+Rbeta.dd[i]*rr;
    if(rr>Rmax.dd[i])   
    { 
        Rt.dd[i] = Rmax.dd[i];
        Rmaxcount[i]++; 
    }
    else
    {
        Rt.dd[i] = (1.0-Rbeta.dd[i])*Rt.dd[i]+Rbeta.dd[i]*rr; 
        Rmaxcount[i] = 0; 
    }
    Rbeta.dd[i] = Rbeta.dd[i]/(Rbeta.dd[i]+Rb.dd[i]);  // beta = beta / (beta+b)
    int adptOK = (Rmaxcount[i]==0||Rmaxcount[i]>Rmaxcount0[i]) ? 1: 0;
    return adptOK;
}

void CKalman::RPkFading(int i)
{
    Zfd.dd[i] = Zfd.dd[i]*(1-zfdafa) + Zk.dd[i]*zfdafa;
    if(Zfd.dd[i]>Zfd0.dd[i] || Zfd.dd[i]<-Zfd0.dd[i])
        DVMDVafa(Fading.GetRow(i), Pk);
}

void CKalman::XPConstrain(void)
{
    int i=0, nq1=nq+1;
    for(double *px=Xk.dd,*pxmax=Xmax.dd,*p=Pk.dd,*pmin=Pmin.dd,*pminEnd=&Pmin.dd[nq],*pmax=Pmax.dd;
        pmin<pminEnd; px++,pxmax++,p+=nq1,pmin++,pmax++)
    {
        if(*px>*pxmax)          // Xk constrain
        {
            *px = *pxmax;
        }
        else if(*px<-*pxmax)
        {
            *px = -*pxmax;
        }
        if(*p<*pmin)    // Pk constrain
        {
            *p = *pmin;
        }
        else if(*p>*pmax)
        {
            double sqf=sqrt(*pmax/(*p))*0.95;
            for(double *prow=&Pk.dd[i*Pk.clm],*prowEnd=prow+nq,*pclm=&Pk.dd[i]; prow<prowEnd; prow++,pclm+=nq)
            {
                *prow *= sqf;
                *pclm *= sqf;
            }
        }
        i++;
    }
}

void CKalman::Feedback(double fbts)
{
    double *pTau=FBTau.dd, *pTotal=FBTotal.dd, *pMax=FBMax.dd, *pXk=FBXk.dd, *p=Xk.dd;
    for(int i=0; i<nq; i++, pTau++,pTotal++,pMax++,pXk++,p++)
    {
        if(*pTau<INF/2)
        {
            double afa = fbts<*pTau ? fbts/(*pTau) : 1.0;
            *pXk = *p*afa;
            if(*pMax<INF/2)
            {
                if(*pTotal+*pXk>*pMax)          *pXk = *pMax-*pTotal;
                else if(*pTotal+*pXk<-*pMax)    *pXk = -*pMax-*pTotal;
            }
            *p -= *pXk;
            *pTotal += *pXk;
        }
        else
        {
            *pXk = 0.0;
        }
    }
}

/***************************  class CSINSKF  *********************************/
CSINSKF::CSINSKF(int nq0, int nr0):CKalman(nq0,nr0)
{
    sins = CSINS(qI, O31, O31);
    SetFt();
    SetHk(); 
}

void CSINSKF::Init(const CSINS &sins0, int grade)
{
    sins = sins0;
    sins.Rwfa.set(
        CVect(9, 100*glv.dps,100*glv.dps,100*glv.dps, 1*glv.g0,1*glv.g0,1*glv.g0, 1*glv.g0,1*glv.g0,1*glv.g0),
        CVect(9, 1.0,1.0,1.0, 1.0,1.0,1.0, 1.0,1.0,1.0),
        CVect(9, 100*glv.dps,100*glv.dps,100*glv.dps, 1*glv.g0,1*glv.g0,1*glv.g0, 1*glv.g0,1*glv.g0,1*glv.g0),
        CVect(9, 0.01*glv.dps,0.01*glv.dps,0.01*glv.dps, 0.1*glv.mg,0.1*glv.mg,0.1*glv.mg, 0.1*glv.mg,0.1*glv.mg,0.1*glv.mg)
        );
    // a example for 15-state(phi,dvn,dpos,eb,db) KF setting
    if(grade==0) // inertial-grade
    {
    Pmax.Set2(10.0*glv.deg,10.0*glv.deg,30.0*glv.deg,    50.0,50.0,50.0,    1.0e4/glv.Re,1.0e4/glv.Re,1.0e4, 
        10.0*glv.dph,10.0*glv.dph,10.0*glv.dph,    10.0*glv.mg,10.0*glv.mg,10.0*glv.mg);
    Pmin.Set2(0.01*glv.min,0.01*glv.min,0.1*glv.min,    0.01,0.01,0.1,    1.0/glv.Re,1.0/glv.Re,0.1, 
        0.001*glv.dph,0.001*glv.dph,0.001*glv.dph,    1.0*glv.ug,1.0*glv.ug,5.0*glv.ug);
    Pk.SetDiag2(1.0*glv.deg,1.0*glv.deg,10.0*glv.deg,    1.0,1.0,1.0,     100.0/glv.Re,100.0/glv.Re,100.0, 
        0.1*glv.dph,0.1*glv.dph,0.1*glv.dph,    1.0*glv.mg,1.0*glv.mg,1.0*glv.mg);
    Qt.Set2(0.001*glv.dpsh,0.001*glv.dpsh,0.001*glv.dpsh,    1.0*glv.ugpsHz,1.0*glv.ugpsHz,1.0*glv.ugpsHz,    0.0,0.0,0.0,
        0.0*glv.dphpsh,0.0*glv.dphpsh,0.0*glv.dphpsh,    0.0*glv.ugpsh,0.0*glv.ugpsh,0.0*glv.ugpsh);
    Xmax.Set(INF,INF,INF, INF,INF,INF, INF,INF,INF,  0.1*glv.dps,0.1*glv.dps,0.1*glv.dps,  200.0*glv.ug,200.0*glv.ug,200.0*glv.ug);
    Rt.Set2(0.2,0.2,0.6,   10.0/glv.Re,10.0/glv.Re, 30.0);
    Rmax = Rt*10000;  Rmin = Rt*0.01;  Rb = 0.9;
    FBTau.Set(1.0,1.0,10.0,     1.0,1.0,1.0,     1.0,1.0,1.0,    10.0,10.0,10.0,    10.0,10.0,10.0);
    }
    else if(grade==1) // MEMS-grade
    {
    Pmax.Set2(50.0*glv.deg,50.0*glv.deg,100.0*glv.deg,    500.0,500.0,500.0,    1.0e6/glv.Re,1.0e6/glv.Re,1.0e6, 
        5000.0*glv.dph,5000.0*glv.dph,5000.0*glv.dph,    50.0*glv.mg,50.0*glv.mg,50.0*glv.mg);
    Pmin.Set2(0.5*glv.min,0.5*glv.min,3.0*glv.min,    0.01,0.01,0.1,    1.0/glv.Re,1.0/glv.Re,0.1, 
        1.0*glv.dph,1.0*glv.dph,1.0*glv.dph,    50.0*glv.ug,50.0*glv.ug,50.0*glv.ug);
    Pk.SetDiag2(10.0*glv.deg,10.0*glv.deg,100.0*glv.deg,    10.0,10.0,10.0,     100.0/glv.Re,100.0/glv.Re,100.0, 
        100.0*glv.dph,101.0*glv.dph,102.0*glv.dph,    1.0*glv.mg,1.01*glv.mg,10.0*glv.mg);
    Qt.Set2(1.0*glv.dpsh,1.0*glv.dpsh,1.0*glv.dpsh,    10.0*glv.ugpsHz,10.0*glv.ugpsHz,10.0*glv.ugpsHz,    0.0,0.0,0.0,
        0.0*glv.dphpsh,0.0*glv.dphpsh,0.0*glv.dphpsh,    0.0*glv.ugpsh,0.0*glv.ugpsh,0.0*glv.ugpsh);
    Xmax.Set(INF,INF,INF, INF,INF,INF, INF,INF,INF,  1.0*glv.dps,1.0*glv.dps,1.0*glv.dps,  50.0*glv.mg,50.0*glv.mg,50.0*glv.mg);
    Rt.Set2(0.2,0.2,0.6,   10.0/glv.Re,10.0/glv.Re,30.0);
    Rmax = Rt*10000;  Rmin = Rt*0.01;  Rb = 0.9;
    FBTau.Set(1.0,1.0,1.0,     1.0,1.0,1.0,     1.0,1.0,1.0,    10.0,10.0,10.0,    10.0,10.0,10.0);
    }
    Zfd0.Set(1.0,15.0,1.0, 100.0/RE,100.0/RE,100.0, INF);   zfdafa = 0.1;
    Fading(0,1)=1.01;  Fading(0,3)=1.01;    Fading(0,7)=1.01;
    Fading(1,0)=1.01;  Fading(1,4)=1.01;    Fading(1,6)=1.01;
    Fading(2,5)=1.01;  Fading(2,8)=1.01;    Fading(2,14)=1.01;
    Fading(3,0)=1.01;  Fading(3,4)=1.01;    Fading(3,6)=1.01;
    Fading(4,1)=1.01;  Fading(4,3)=1.01;    Fading(4,7)=1.01;
    Fading(5,5)=1.01;  Fading(5,8)=1.01;    Fading(5,14)=1.01;
    Zfd0 = INF;
}

void CSINSKF::SetFt(void)
{
    sins.etm();
//  Ft = [ Maa    Mav    Map    -Cnb    O33 
//         Mva    Mvv    Mvp     O33    Cnb 
//         O33    Mpv    Mpp     O33    O33
//         zeros(6,9)  diag(-1./[ins.tauG;ins.tauA]) ];
    Ft.SetMat3(0,0,sins.Maa), Ft.SetMat3(0,3,sins.Mav), Ft.SetMat3(0,6,sins.Map), Ft.SetMat3(0,9,-sins.Cnb); 
    Ft.SetMat3(3,0,sins.Mva), Ft.SetMat3(3,3,sins.Mvv), Ft.SetMat3(3,6,sins.Mvp), Ft.SetMat3(3,12,sins.Cnb); 
                              Ft.SetMat3(6,3,sins.Mpv), Ft.SetMat3(6,6,sins.Mpp);
    Ft.SetDiagVect3( 9, 9, sins._betaGyro);
    Ft.SetDiagVect3(12,12, sins._betaAcc);
}

void CSINSKF::SetHk(void)
{
    // an example for SINS/GPS vn&pos measurement
    Hk(0,3) = Hk(1,4) = Hk(2,5) = 1.0;
    Hk(3,6) = Hk(4,7) = Hk(5,8) = 1.0;
}

int CSINSKF::Update(const CVect3 *pwm, const CVect3 *pvm, int nSamples, double ts)
{
    sins.Update(pwm, pvm, nSamples, ts);
    TimeUpdate(sins.nts);  kftk = sins.tk;
    return MeasUpdate();
}

void CSINSKF::Feedback(double fbts)
{
    CKalman::Feedback(fbts);
    sins.qnb -= *(CVect3*)&FBXk.dd[0];  sins.vn -= *(CVect3*)&FBXk.dd[ 3];  sins.pos -= *(CVect3*)&FBXk.dd[6];
    sins.eb  += *(CVect3*)&FBXk.dd[9];  sins.db += *(CVect3*)&FBXk.dd[12]; 
}

void CSINSKF::QtMarkovGA(const CVect3 &tauG, const CVect3 &sRG, const CVect3 &tauA=I31, const CVect3 &sRA=O31)
{
    sins.SetTauGA(tauG, tauA);
    *(CVect3*)&Qt.dd[ 9] = MKQt(sRG, sins.tauGyro);
    *(CVect3*)&Qt.dd[12] = MKQt(sRA, sins.tauAcc);
}

void CSINSKF::Secret(void)
{
#define SC_PCH (12*DEG)
#define SC_RLL (34*DEG)
#define SC_YAW (56*DEG)
#define SC_THR (0.1*DEG)
    if( sins.att.i<SC_PCH+SC_THR && sins.att.i>SC_PCH-SC_THR &&
        sins.att.j<SC_RLL+SC_THR && sins.att.j>SC_RLL-SC_THR &&
        sins.att.k<SC_YAW+SC_THR && sins.att.k>SC_YAW-SC_THR )
    {
        sins.att.k = SC_YAW;
        sins.qnb.SetYaw(sins.att.k); 
    }
#undef SC_PCH
#undef SC_RLL
#undef SC_YAW
#undef SC_THR
}

/***************************  class CSINSTDKF  *********************************/
CSINSTDKF::CSINSTDKF(int nq0, int nr0):CSINSKF(nq0,nr0)
{
    Fk = Pk1 = CMat(nq,nq, 0.0);
    Pxz = Qk = Kk = tmeas = CVect(nr, 0.0);
    tdts = 0.0;
    maxStep = 2*(nq+nr)+3;
    TDReset();
}

void CSINSTDKF::TDReset(void)
{
    iter = -2;
    ifn = 0;    meanfn = O31;
    SetMeasFlag(0);
}

/**
 * @brief 
 * @param[in] pwm       Angle increment
 * @param[in] pvm       Velocity increment
 * @param[in] nSamples  Number of samples
 * @param[in] ts        The Sampling period
 * @param[in] nStep     
 * @return 
 */
int CSINSTDKF::TDUpdate(const CVect3 *pwm, const CVect3 *pvm, 
        int nSamples, double ts, int nStep)
{
//  if(sins.tk>99)
//      int debugi = 1;

    sins.Update(pwm, pvm, nSamples, ts);
    Feedback(sins.nts);

    measRes = 0;

    if(nStep<=0||nStep>=maxStep) { nStep=maxStep; }
    tdStep = nStep;

    tdts += sins.nts; kftk = sins.tk;
    meanfn = meanfn+sins.fn; ifn++;
    for(int i=0; i<nStep; i++)
    {
        if(iter==-2)            // -2: set measurements
        {
            if(ifn==0)  break;
            CVect3 vtmp=meanfn*(1.0/ifn); meanfn = O31; ifn = 0;
            sins.fn=vtmp; SetFt(); sins.fn = vtmp;          
            SetMeas();
        }
        else if(iter==-1)           // -1: discrete
        {
            Fk = ++(Ft*tdts); // Fk = I+Ft*ts
            Qk = Qt*tdts;
            Xk = Fk*Xk;
            tdts = 0.0;
        }
        else if(iter<nq)        // 0 -> (nq-1): Fk*Pk
        {
            int row=iter;
            RowMul(Pk1, Fk, Pk, row);
        }
        else if(iter<2*nq)      // nq -> (2*nq-1): Fk*Pk*Fk+Qk
        {
            int row=iter-nq;
            RowMulT(Pk, Pk1, Fk, row);
            Pk.dd[nq*row+row] += Qk.dd[row];
//          if(row==nq-1) { Pk += Qk; }
        }
        else if(iter<2*(nq+nr)) // (2*nq) -> (2*(nq+nr)-1): sequential measurement updating
        {
            int row=(iter-2*Ft.row)/2;
            int flag = measflag&(0x01<<row);
            if(flag)
            {
                if((iter-2*Ft.row)%2==0)
                {
                    Hi = Hk.GetRow(row);
                    Pxz = Pk*(~Hi);
                    Pz0 = (Hi*Pxz)(0,0);
                    innovation = Zk(row)-(Hi*Xk)(0,0);
                    adptOKi = 1;
                    if(Rb.dd[row]>EPS)
                        adptOKi = RAdaptive(row, innovation, Pz0);
                    double Pzz = Pz0 + Rt(row)/rts(row);
                    Kk = Pxz*(1.0/Pzz);
                }
                else
                {
                    measflag ^= flag;
                    if(adptOKi)
                    {
                        measRes |= flag;
                        Xk += Kk*innovation;
                        Pk -= Kk*(~Pxz);
                    }
                    if(Zfd0.dd[row]<INF/2)
                    {
                        RPkFading(row);
                    }
                }
            }
            else
            {
                nStep++;
            }
        }
        else if(iter==2*(nq+nr))    // 2*(nq+nr): Xk,Pk constrain & symmetry
        {
            XPConstrain();
            symmetry(Pk);
        }
        else if(iter>=2*(nq+nr)+1)  // 2*(nq+nr)+1: Miscellanous
        {
            Miscellanous();
            iter = -3;
        }
        iter++;
    }
    Secret();

    return measRes;
}

/***************************  class CSINSGPS  *********************************/
CSINSGPSOD::CSINSGPSOD(void):CSINSTDKF(15, 10)
{
    vbINS = vbOD = O31;  tODInt = 0.0;
    Cbo = I33;
    Hk(9,2) = 1.0;
    measGPSvnValid = measGPSposValid = measODvbValid = measMAGyawValid = 0;
}

#ifndef CSINSGPSOD_INIT_BY_USER
// may be copied & implemented by user
void CSINSGPSOD::Init(const CSINS &sins0, int grade)
{
    CSINSKF::Init(sins0, grade);
    Qt.Set2(1.0*glv.dpsh,1.0*glv.dpsh,1.0*glv.dpsh,    10.0*glv.ugpsHz,10.0*glv.ugpsHz,10.0*glv.ugpsHz,    0.0,0.0,0.0,
        0.0*glv.dphpsh,0.0*glv.dphpsh,0.0*glv.dphpsh,    0.0*glv.ugpsh,0.0*glv.ugpsh,0.0*glv.ugpsh);
    Rt.Set2(0.2,0.2,0.6,   10.0/glv.Re,10.0/glv.Re,30.0, 0.1,1000.0,0.1, 1.0*glv.deg);
    Rmax = Rt*10000;  Rmin = Rt*0.01;  Rb = 0.9;
}
#else
    #pragma message("  CSINSGPSOD_INIT_BY_USER")
#endif

void CSINSGPSOD::SetMeas(void)
{
    if(measGPSposValid)
    {
        SetMeasFlag(0070);
        measGPSposValid = 0;
    }
    if(measGPSvnValid)
    {
        SetMeasFlag(0007);
        measGPSvnValid = 0;
    }
    if(measODvbValid)
    {
        SetMeasFlag(0700);
        measODvbValid = 0;
    }
    if(measMAGyawValid)
    {
        SetMeasFlag(01000);
        measMAGyawValid = 0;
    }
}

void CSINSGPSOD::SetMeasGPS(const CVect3 &pgps, const CVect3 &vgps)
{
    if(!IsZero(pgps))
    {
        *(CVect3*)&Zk.dd[3] = sins.pos - pgps;
        measGPSposValid = 1;
    }
    if(!IsZero(vgps))
    {
        *(CVect3*)&Zk.dd[0] = sins.vn - vgps;
        measGPSvnValid = 1;
    }
}

void CSINSGPSOD::SetMeasOD(double dSod, double ts)
{
    if(sins.wnb.k<5.0*glv.dps)
    {
        tODInt += ts;
        vbINS += sins.vb;
        vbOD += Cbo*CVect3(0,dSod/ts,0);
        if(tODInt>1.0)
        {
            *(CVect3*)&Zk.dd[6] = vbINS - vbOD;
            Hk.SetMat3(6, 3, sins.Cbn);
            vbINS = vbOD = O31; tODInt = 0.0;
            measODvbValid = 1;
        }
    }
    else
    {
        vbINS = vbOD = O31; tODInt = 0.0;
    }
}

void CSINSGPSOD::SetMeasYaw(double ymag)
{
    if(ymag>EPS || ymag<-EPS)  // !IsZero(yawGPS)
        if(sins.att.i<30*DEG && sins.att.i>-30*DEG)
        {
            *(CVect3*)&Zk.dd[9] = -diffYaw(sins.att.k, ymag);
            measMAGyawValid = 1;
        }
}

/***************************  class CIMU  *********************************/
CIMU::CIMU(void)
{
    nSamples = prefirst = 1;
    phim = dvbm = wm_1 = vm_1 = O31;
}

void CIMU::Update(const CVect3 *pwm, const CVect3 *pvm, int nSamples)
{
    static double conefactors[5][4] = {             // coning coefficients
        {2./3},                                     // 2
        {9./20, 27./20},                            // 3
        {54./105, 92./105, 214./105},               // 4
        {250./504, 525./504, 650./504, 1375./504}   // 5
        };
    int i;
    double *pcf = conefactors[nSamples-2];
    CVect3 cm(0.0), sm(0.0), wmm(0.0), vmm(0.0);

    this->nSamples = nSamples;
    if(nSamples==1)  // one-plus-previous sample
    {
        if(prefirst==1) {wm_1=pwm[0]; vm_1=pvm[0]; prefirst=0;}
        cm = 1.0/12*wm_1; wm_1=pwm[0]; 
        sm = 1.0/12*vm_1; vm_1=pvm[0];
    }
    if(nSamples>1) prefirst=1;
    for(i=0; i<nSamples-1; i++)
    {
        cm += pcf[i]*pwm[i];
        sm += pcf[i]*pvm[i];
        wmm += pwm[i];
        vmm += pvm[i];
    }
    wmm += pwm[i];
    vmm += pvm[i];
    phim = wmm + cm*pwm[i];
    dvbm = vmm + 1.0/2*wmm*vmm + (cm*pvm[i]+sm*pwm[i]);
}

void IMURFU(CVect3 *pwm, int nSamples, const char *str)
{
    for(int n=0; n<nSamples; n++)
    {
        CVect3 tmpwm;
        double *pw=(double*)&pwm[n].i;
        for(int i=0; i<3; i++,pw++)
        {
            switch(str[i])
            {
            case 'R':  tmpwm.i= *pw;  break;
            case 'L':  tmpwm.i=-*pw;  break;
            case 'F':  tmpwm.j= *pw;  break;
            case 'B':  tmpwm.j=-*pw;  break;
            case 'U':  tmpwm.k= *pw;  break;
            case 'D':  tmpwm.k=-*pw;  break;
            }
        }
        pwm[n] = tmpwm;
    }
}

void IMURFU(CVect3 *pwm, CVect3 *pvm, int nSamples, const char *str)
{
    IMURFU(pwm, nSamples, str);
    IMURFU(pvm, nSamples, str);
}

/***************************  class CSINS  *********************************/
CSINS::CSINS(const CQuat &qnb0, const CVect3 &vn0, const CVect3 &pos0, double tk0)
{
    tk = tk0;  ts = nts = 1.0;
    velMax = 400.0; hgtMin = -RE*0.01, hgtMax = -hgtMin;
    qnb = qnb0; vn = vn0, pos = pos0;
    Kg = Ka = I33; eb = db = O31;
    Maa = Mav = Map = Mva = Mvv = Mvp = Mpv = Mpp = O33;
    SetTauGA(CVect3(INF),CVect3(INF));
    CVect3 wib(0.0), fb=(~qnb)*CVect3(0,0,glv.g0);
    lvr = an = O31;
    Rwfa = CRAvar(9);
    Update(&wib, &fb, 1, 1.0); tk = tk0;  ts = nts = 1.0; qnb = qnb0;   vn = vn0, pos = pos0;
    etm(); lever(); Extrap();
}

void CSINS::SetTauGA(const CVect3 &tauG, const CVect3 &tauA)
{
    tauGyro = tauG, tauAcc = tauA;
    _betaGyro.i = tauG.i>INF/2 ? 0.0 : -1.0/tauG.i;   // Gyro&Acc inverse correlation time for AR(1) model
    _betaGyro.j = tauG.j>INF/2 ? 0.0 : -1.0/tauG.j;
    _betaGyro.k = tauG.k>INF/2 ? 0.0 : -1.0/tauG.k;
    _betaAcc.i  = tauA.i>INF/2 ? 0.0 : -1.0/tauA.i;
    _betaAcc.j  = tauA.j>INF/2 ? 0.0 : -1.0/tauA.j;
    _betaAcc.k  = tauA.k>INF/2 ? 0.0 : -1.0/tauA.k;
}

void CSINS::Update(const CVect3 *pwm, const CVect3 *pvm, int nSamples, double ts)
{
#ifdef PSINS_LOW_GRADE_MEMS
#pragma message("  PSINS_LOW_GRADE_MEMS")
    this->ts = ts;  nts = nSamples*ts;  tk += nts;
    double nts2 = nts/2;
    imu.Update(pwm, pvm, nSamples);
    imu.phim = Kg*imu.phim - eb*nts; imu.dvbm = Ka*imu.dvbm - db*nts;  // IMU calibration
//  CVect3 vn01 = vn+an*nts2, pos01 = pos+eth.vn2dpos(vn01,nts2);
    eth.Update(pos);
    wib = imu.phim/nts; fb = imu.dvbm/nts;
    web = wib;
    wnb = wib;
    fn = qnb*fb;
    an = fn+eth.gcc;
    CVect3 vn1 = vn + an*nts;
    pos = pos + eth.vn2dpos(vn+vn1, nts2);  vn = vn1;
    Cnb0 = Cnb;
    qnb = qnb*rv2q(imu.phim);
    Cnb = q2mat(qnb); att = m2att(Cnb); Cbn = ~Cnb; vb = Cbn*vn;
#else
    this->ts = ts;  nts = nSamples*ts;  tk += nts;
    double nts2 = nts/2;
    imu.Update(pwm, pvm, nSamples);
    imu.phim = Kg*imu.phim - eb*nts; imu.dvbm = Ka*imu.dvbm - db*nts;  // IMU calibration
    CVect3 vn01 = vn+an*nts2, pos01 = pos+eth.vn2dpos(vn01,nts2);
    eth.Update(pos01, vn01);
    wib = imu.phim/nts; fb = imu.dvbm/nts;
    web = wib - Cbn*eth.wnie;
    wnb = wib - (qnb*rv2q(imu.phim/2))*eth.wnin;
    fn = qnb*fb;
    an = rv2q(-eth.wnin*nts2)*fn+eth.gcc;
    CVect3 vn1 = vn + an*nts;
    pos = pos + eth.vn2dpos(vn+vn1, nts2);  vn = vn1;
    Cnb0 = Cnb;
    qnb = rv2q(-eth.wnin*nts)*qnb*rv2q(imu.phim);
    Cnb = q2mat(qnb); att = m2att(Cnb); Cbn = ~Cnb; vb = Cbn*vn;
#endif
    psinsassert(pos.i<85.0*PI/180 && pos.i>-85.0*PI/180);
    if(vn.i>velMax) vn.i=velMax; else if(vn.i<-velMax) vn.i=-velMax;
    if(vn.j>velMax) vn.j=velMax; else if(vn.j<-velMax) vn.j=-velMax;
    if(vn.k>velMax) vn.k=velMax; else if(vn.k<-velMax) vn.k=-velMax;
    if(pos.j>PI) pos.j-=_2PI; else if(pos.j<-PI) pos.j+=_2PI;
    if(pos.k>hgtMax) pos.k=hgtMax; else if(pos.k<hgtMin) pos.k=hgtMin;
    CVect wfa(9);
    *(CVect3*)&wfa.dd[0]=wib, *(CVect3*)&wfa.dd[3]=fb, *(CVect3*)&wfa.dd[6]=an;
    Rwfa.Update(wfa, nts);
}

void CSINS::Extrap(const CVect3 &wm, const CVect3 &vm, double ts)
{
    if(ts<1.0e-6)  // reset
    {
        qnbE = qnb, vnE = vn, posE = pos, attE = att;
    }
    else
    {
        vnE = vnE + qnbE*vm + eth.gcc*ts;
        posE = posE + eth.vn2dpos(vnE,ts);
        qnbE = qnbE*rv2q(wm); attE = q2att(qnbE);
    }
}

void CSINS::lever(const CVect3 &dL)
{
    if(!IsZero(dL)) lvr = dL;
    Mpv = CMat3(0,eth.f_RMh,0, eth.f_clRNh,0,0, 0,0,1);
    CW = Cnb*askew(web), MpvCnb = Mpv*Cnb;
    vnL = vn + CW*lvr; posL = pos + MpvCnb*lvr;
}

void CSINS::etm(void)
{
#ifdef PSINS_LOW_GRADE_MEMS
    Mva = askew(fn);
    Mpv = CMat3(0,eth.f_RMh,0, eth.f_clRNh,0,0, 0,0,1);
#else
    double tl=eth.tl, secl=1.0/eth.cl, secl2=secl*secl, 
        wN=eth.wnie.j, wU=eth.wnie.k, vE=vn.i, vN=vn.j;
    double f_RMh=eth.f_RMh, f_RNh=eth.f_RNh, f_clRNh=eth.f_clRNh, 
        f_RMh2=f_RMh*f_RMh, f_RNh2=f_RNh*f_RNh;
    CMat3 Avn=askew(vn),
        Mp1(0,0,0, -wU,0,0, wN,0,0),
        Mp2(0,0,vN*f_RMh2, 0,0,-vE*f_RNh2, vE*secl2*f_RNh,0,-vE*tl*f_RNh2);
    Maa = askew(-eth.wnin);
    Mav = CMat3(0,-f_RMh,0, f_RNh,0,0, tl*f_RNh,0,0);
    Map = Mp1+Mp2;
    Mva = askew(fn);
    Mvv = Avn*Mav - askew(eth.wnie+eth.wnin);
    Mvp = Avn*(Mp1+Map);
    double scl = eth.sl*eth.cl;
    Mvp.e20 = Mvp.e20-glv.g0*(5.27094e-3*2*scl+2.32718e-5*4*eth.sl2*scl); Mvp.e22 = Mvp.e22+3.086e-6;
    Mpv = CMat3(0,f_RMh,0, f_clRNh,0,0, 0,0,1);
    Mpp = CMat3(0,0,-vN*f_RMh2, vE*tl*f_clRNh,0,-vE*secl*f_RNh2, 0,0,0);
#endif
}

#ifdef PSINS_AHRS_MEMS
#pragma message("  PSINS_AHRS_MEMS")

/*********************  class CMahony AHRS  ************************/
CMahony::CMahony(double tau, const CQuat &qnb0)
{
    SetTau(tau);
    qnb = qnb0;
    Cnb = q2mat(qnb);
    exyzInt = O31;  ebMax = I31*glv.dps;
    tk = 0.0;
}

void CMahony::SetTau(double tau)
{
    double beta = 2.146/tau;
    Kp = 2.0f*beta, Ki = beta*beta;
}

void CMahony::Update(const CVect3 &wm, const CVect3 &vm, double ts, const CVect3 &mag)
{
    double nm;
    CVect3 acc0, mag0, exyz, bxyz, wxyz;

    nm = norm(vm)/ts;  // f (in m/s^2)
    acc0 = nm>0.1 ? vm/(nm*ts) : O31;
    nm = norm(mag);    // mag (in Gauss)
    if(nm>0.1)
    {
        mag0 = mag/nm;
        bxyz = Cnb*mag0;
        bxyz.j = normXY(bxyz); bxyz.i = 0.0;
        wxyz = (~Cnb)*bxyz;
    }
    else
    {
        mag0 = wxyz = O31;
    }
    exyz = *((CVect3*)&Cnb.e20)*acc0 + wxyz*mag0;
    exyzInt += exyz*(Ki*ts);
    qnb *= rv2q(wm-(Kp*exyz+exyzInt)*ts);
    Cnb = q2mat(qnb);
    tk += ts;
    if(exyzInt.i>ebMax.i)  exyzInt.i=ebMax.i;  else if(exyzInt.i<-ebMax.i)  exyzInt.i=-ebMax.i;
    if(exyzInt.j>ebMax.j)  exyzInt.j=ebMax.j;  else if(exyzInt.j<-ebMax.j)  exyzInt.j=-ebMax.j;
    if(exyzInt.k>ebMax.k)  exyzInt.k=ebMax.k;  else if(exyzInt.k<-ebMax.k)  exyzInt.k=-ebMax.k;
}

void CMahony::Update(const CVect3 &gyro, const CVect3 &acc, const CVect3 &mag, double ts)
{
    double nm;
    CVect3 acc0, mag0, exyz, bxyz, wxyz;

    nm = norm(acc);
    acc0 = nm>0.1 ? acc/nm : O31;
    nm = norm(mag);
    mag0 = nm>0.1 ? mag/nm : O31;
    bxyz = Cnb*mag0;
    bxyz.j = normXY(bxyz); bxyz.i = 0.0;
    wxyz = (~Cnb)*bxyz;
    exyz = *((CVect3*)&Cnb.e20)*acc0 + wxyz*mag0;
    exyzInt += exyz*(Ki*ts);
    qnb *= rv2q((gyro*glv.dps-Kp*exyz-exyzInt)*ts);
    Cnb = q2mat(qnb);
    tk += ts;
}

/*********************  class Quat&EKF based AHRS  ************************/
CQEAHRS::CQEAHRS(double ts):CKalman(7,3)
{
    double sts = sqrt(ts);
    Pmax.Set2(2.0,2.0,2.0,2.0, 1000*glv.dph,1000.0*glv.dph,1000.0*glv.dph);
    Pmin.Set2(0.001,0.001,0.001,0.001, 10.0*glv.dph,10.0*glv.dph,10.0*glv.dph);
    Pk.SetDiag2(1.0,1.0,1.0,1.0, 1000.0*glv.dph,1000.0*glv.dph,1000.0*glv.dph);
    Qt.Set2(10.0*glv.dpsh,10.0*glv.dpsh,10.0*glv.dpsh, 10.0*glv.dphpsh,10.0*glv.dphpsh,10.0*glv.dphpsh);
    Rt.Set2(100.0*glv.mg/sts,100.0*glv.mg/sts, 1.0*glv.deg/sts);
    Xk(0) = 1.0;
    Cnb = q2mat(*(CQuat*)&Xk.dd[0]);
}

void CQEAHRS::Update(const CVect3 &gyro, const CVect3 &acc, const CVect3 &mag, double ts)
{
    double q0, q1, q2, q3, wx, wy, wz, fx, fy, fz, mx, my, mz, h11, h12, h21, h22; 
    q0 = Xk.dd[0],       q1 = Xk.dd[1],       q2 = Xk.dd[2],        q3 = Xk.dd[3];
    wx = gyro.i*glv.dps, wy = gyro.j*glv.dps, wz = gyro.k*glv.dps; 
    fx = acc.i,          fy = acc.j,          fz = acc.k; 
    mx = mag.i,          my = mag.j,          mz = mag.k; 
    // Ft
                    0, Ft.dd[ 1] = -wx/2, Ft.dd[ 2] = -wy/2, Ft.dd[ 3] = -wz/2,  Ft.dd[ 4] =  q1/2, Ft.dd[ 5] =  q2/2, Ft.dd[ 6] =  q3/2; 
    Ft.dd[ 7] =  wx/2,                 0, Ft.dd[ 9] =  wz/2, Ft.dd[10] = -wy/2,  Ft.dd[11] = -q0/2, Ft.dd[12] =  q3/2, Ft.dd[13] = -q2/2; 
    Ft.dd[14] =  wy/2, Ft.dd[15] = -wz/2,                 0, Ft.dd[17] =  wx/2,  Ft.dd[18] = -q3/2, Ft.dd[18] = -q0/2, Ft.dd[20] =  q1/2; 
    Ft.dd[21] =  wz/2, Ft.dd[22] =  wy/2, Ft.dd[23] = -wx/2,                 0,  Ft.dd[25] =  q2/2, Ft.dd[26] = -q1/2, Ft.dd[27] = -q0/2; 
    // Hk
    h11 = fx*q0-fy*q3+fz*q2;  h12 = fx*q1+fy*q2+fz*q3;
    h21 = fx*q3+fy*q0-fz*q1;  h22 = fx*q2-fy*q1-fz*q0;
    Hk.dd[ 0] = h11*2,  Hk.dd[ 1] = h12*2,  Hk.dd[ 2] = -h22*2,  Hk.dd[ 3] = -h21*2;
    Hk.dd[ 7] = h21*2,  Hk.dd[ 8] = h22*2,  Hk.dd[ 9] =  h12*2,  Hk.dd[10] =  h11*2;
/*  CVect3 magH = Cnb*mag;
    double C11=Cnb.e11, C01=Cnb.e01, CC=C11*C11+C01*C01;
    if(normXY(magH)>0.01 && CC>0.25)  // CC>0.25 <=> pitch<60deg
    {
        double f2=2.0/CC;
        Hk.dd[14] = (q3*C11+q0*C01)*f2,  Hk.dd[15] = (-q2*C11-q1*C01)*f2,  Hk.dd[16] = (-q1*C11+q2*C01)*f2,  Hk.dd[17] = (q0*C11-q3*C01)*f2;
        Zk.dd[2] = atan2(magH.i, magH.j);
    }
    else
    {
        Hk.dd[14] = Hk.dd[15] = Hk.dd[16] = Hk.dd[17] = 0.0;
        Zk.dd[2] = 0.0;
    }*/

    SetMeasFlag(0x03);
    TimeUpdate(ts);
    MeasUpdate();
    XPConstrain();
    normlize((CQuat*)&Xk.dd[0]);
    Cnb = q2mat(*(CQuat*)&Xk.dd[0]);
}

#endif  // PSINS_AHRS_MEMS

/******************************  File Read or Write *********************************/
#ifdef PSINS_IO_FILE
#pragma message("  PSINS_IO_FILE")

//#include "io.h"
#include <stdio.h>
char* time2fname(void)
{
    static char PSINSfname[32];
    time_t tt;  time(&tt);
    tm *Time = localtime(&tt);
    strftime(PSINSfname, 32, "PSINS%Y%m%d_%H%M%S.bin", Time);
    return PSINSfname;
}

char CFileRdWt::dirIn[256] = {0}, CFileRdWt::dirOut[256] = {0};

void CFileRdWt::Dir(const char *dirI, const char *dirO)  // set dir
{
    int len = strlen(dirI);
    memcpy(dirIn, dirI, len);
    if(dirIn[len-1]!='/')
    { 
        dirIn[len]='/';
        dirIn[len+1]='\0';
    }
    if(dirO)
    {
        len = strlen(dirO);
        memcpy(dirOut, dirO, len);
        if(dirOut[len-1]!='/') { dirOut[len]='/'; dirOut[len+1]='\0'; }
    }
    else
        memcpy(dirOut, dirIn, strlen(dirIn));
}

CFileRdWt::CFileRdWt()
{
    f = 0;
}

CFileRdWt::CFileRdWt(const char *fname0, int columns0)
{
    Init(fname0, columns0);
    memset(buff, 0, sizeof(buff));
}

void CFileRdWt::Init(const char *fname0, int columns0)
{
    fname[0]='\0';
    int findc=0, len0=strlen(fname0);
    for(int i=0; i<len0; i++)   { if(fname0[i]=='\\') { findc=1; break; } }
    columns = columns0;
    if(columns==0)      // file write
    {   if(dirOut[0]!=0&&findc==0)  { strcat(fname, dirOut); } }
    else                // file read
    {   if(dirIn[0]!=0&&findc==0)   { strcat(fname, dirIn); } }
    strcat(fname, fname0);
    if(columns==0)              // bin file write
    {
        f = fopen(fname, "wb");
    }
    else if(columns<0)          // bin file read
    {
        f = fopen(fname, "rb");
    }
    else if(columns>0)          // txt file read
    {
        f = fopen(fname, "rt");
        if(!f){
            LOG(WARNING)<<"file "<<fname<<" do NOT exist";
            return;
        }
        fpos_t pos;
        while(1)  // skip txt-file comments
        {
            // pos = ftell(f);
            fgetpos(f,&pos);
            fgets(line, sizeof(line), f);
            if(feof(f)) break;
            int allSpace=1, allDigital=1;
            for(int i=0; i<sizeof(line); i++)
            {
                char c = line[i];
                if(c=='\n') break;
                if(c!=' ') allSpace = 0;
                if( !(c==' '||c==','||c==';'||c==':'||c=='+'||c=='-'||c=='.'||c=='\t'
                    ||c=='e'||c=='E'||c=='d'||c=='D'||('9'>=c&&c>='0')) )
                    allDigital = 0;
            }
            if(!allSpace && allDigital) break;
        }
        fsetpos(f, &pos);
        // this->columns = columns;
        for(int i=0; i<columns; i++)
        { 
            sstr[4*i+0]='%'; 
            sstr[4*i+1]='l'; 
            sstr[4*i+2]='f';
            sstr[4*i+3]=' ';
            sstr[4*i+4]='\0';
        } 
    }
    else
    {
        f = 0;
    }
    linelen = 0;
}

/**
 * @brief Parse  
 * @param[in] lines 
 * @param[in] txtDelComma
 * @return 
 *      @retval 0 
 *      @retval 1
 */
int CFileRdWt::load(int lines, BOOL txtDelComma)
{
    if(columns<0)           // bin file read
    {
        if(lines>1)
            fseek(f, (lines-1)*(-columns)*sizeof(double), SEEK_CUR);
        fread(buff, -columns, sizeof(double), f);
    }
    else                    // txt file read
    {
        for(int i=0; i<lines; i++)  fgets(line, sizeof(line), f);
        // replace other Separator with ' '
        if(txtDelComma)
        {
            for(char *pc=line, *pend=line+sizeof(line); pc<pend; pc++)
            {
                if(*pc==','||*pc==';'||*pc==':'||*pc=='\t') *pc=' ';
                else if(*pc=='\0') break;
            }
        }
        if(columns<10)
            sscanf(line, sstr,
                &buff[ 0], &buff[ 1], &buff[ 2], &buff[ 3], &buff[ 4],
                &buff[ 5], &buff[ 6], &buff[ 7], &buff[ 8], &buff[ 9] ); 
        else if(columns<20)
            sscanf(line, sstr,
                &buff[ 0], &buff[ 1], &buff[ 2], &buff[ 3], &buff[ 4],
                &buff[ 5], &buff[ 6], &buff[ 7], &buff[ 8], &buff[ 9],
                &buff[10], &buff[11], &buff[12], &buff[13], &buff[14],
                &buff[15], &buff[16], &buff[17], &buff[18], &buff[19] ); 
        else if(columns<40)
            sscanf(line, sstr,
                &buff[ 0], &buff[ 1], &buff[ 2], &buff[ 3], &buff[ 4],
                &buff[ 5], &buff[ 6], &buff[ 7], &buff[ 8], &buff[ 9],
                &buff[10], &buff[11], &buff[12], &buff[13], &buff[14],
                &buff[15], &buff[16], &buff[17], &buff[18], &buff[19],
                &buff[20], &buff[21], &buff[22], &buff[23], &buff[24],
                &buff[25], &buff[26], &buff[27], &buff[28], &buff[29],
                &buff[30], &buff[31], &buff[32], &buff[33], &buff[34],
                &buff[35], &buff[36], &buff[37], &buff[38], &buff[39] ); 
    }
    linelen += lines;
    if(feof(f))  return 0;
    else return 1;
}

int CFileRdWt::loadf32(int lines)   // float32 bin file read
{
    if(lines>1)
        fseek(f, (lines-1)*(-columns)*sizeof(float), SEEK_CUR);
    fread(buff32, -columns, sizeof(float), f);
    for(int i=0; i<-columns; i++) buff[i]=buff32[i];
    linelen += lines;
    if(feof(f))  return 0;
    else return 1;
}

int CFileRdWt::getl(void)   // txt file get a line
{
    fgets(line, sizeof(line), f);
    return strlen(line);
}

CFileRdWt& CFileRdWt::operator<<(double d)
{
    fwrite(&d, 1, sizeof(double), f);
    return *this;
}

CFileRdWt& CFileRdWt::operator<<(const CVect3 &v)
{
    fwrite(&v, 1, sizeof(v), f);
    return *this;
}

CFileRdWt& CFileRdWt::operator<<(const CVect &v)
{
    fwrite(v.dd, v.clm*v.row, sizeof(double), f);
    return *this;
}

CFileRdWt& CFileRdWt::operator<<(const CMat &m)
{
    fwrite(m.dd, m.clm*m.row, sizeof(double), f);
    return *this;
}

CFileRdWt& CFileRdWt::operator<<(const CRAvar &R)
{
    fwrite(R.R0, R.nR0, sizeof(double), f);
    return *this;
}

CFileRdWt& CFileRdWt::operator<<(const CAligni0 &aln)
{
    return *this<<q2att(aln.qnb)<<aln.vib0<<aln.Pi02<<aln.tk;
}

CFileRdWt& CFileRdWt::operator<<(const CSINS &sins)
{
    return *this<<sins.att<<sins.vn<<sins.pos<<sins.eb<<sins.db<<sins.tk;
}

#ifdef PSINS_AHRS_MEMS
CFileRdWt& CFileRdWt::operator<<(const CMahony &ahrs)
{
    return *this<<m2att(ahrs.Cnb)<<ahrs.exyzInt<<ahrs.tk;
}

CFileRdWt& CFileRdWt::operator<<(const CQEAHRS &ahrs)
{
    return *this<<m2att(ahrs.Cnb)<<*(CVect3*)&ahrs.Xk.dd[4]<<diag(ahrs.Pk)<<ahrs.kftk;
}
#endif

#ifdef PSINS_UART_PUSH_POP
CFileRdWt& CFileRdWt::operator<<(const CUartPP &uart)
{
    fwrite(uart.popbuf, uart.frameLen, sizeof(char), f);
    return *this;
}
#endif

CFileRdWt& CFileRdWt::operator<<(const CKalman &kf)
{
    return *this<<kf.Xk<<diag(kf.Pk)<<kf.Rt<<kf.kftk;
}

CFileRdWt::~CFileRdWt()
{
    if(f) { fclose(f); f=(FILE*)NULL; } 
}

CFileRdWt& CFileRdWt::operator>>(double &d)
{
    fread(&d, 1, sizeof(double), f);
    return *this;
}

CFileRdWt& CFileRdWt::operator>>(CVect3 &v)
{
    fread(&v, 1, sizeof(v), f);
    return *this;
}

CFileRdWt& CFileRdWt::operator>>(CVect &v)
{
    fread(v.dd, v.clm*v.row, sizeof(double), f);
    return *this;
}

CFileRdWt& CFileRdWt::operator>>(CMat &m)
{
    fread(m.dd, m.clm*m.row, sizeof(double), f);
    return *this;
}

#endif // PSINS_IO_FILE

#ifdef PSINS_RMEMORY
#pragma message("  PSINS_RMEMORY")

CRMemory::CRMemory(BYTE *pMem, long memLen0, BYTE recordLen0)
{
    psinsassert(recordLen0<=MAX_RECORD_BYTES);
    pMemStart = pMemPush = pMemPop = pMem;
    pMemEnd = pMemStart + memLen0;
    pushLen = popLen = recordLen = recordLen0;
    memLen = memLen0;
    dataLen = 0;
}

BYTE CRMemory::pop(BYTE *p)
{
    if(dataLen==0) return 0;
    popLen = recordLen==0 ? *pMemPop : recordLen;
    if(p==(BYTE*)NULL) p = popBuf;
    BYTE i;
    for(i=0; i<popLen; i++,dataLen--)
    {
        *p++ = *pMemPop++;
        if(pMemPop>=pMemEnd)  pMemPop = pMemStart;
    }
    return i;
}

BOOL CRMemory::push(const BYTE *p)
{
    BOOL res = 1;
    if(p==(BYTE*)NULL) p = pushBuf;
    pushLen = recordLen==0 ? *p : recordLen;
    psinsassert(pushLen<=MAX_RECORD_BYTES);
    for(BYTE i=0; i<pushLen; i++,dataLen++)
    {
        *pMemPush++ = *p++;
        if(pMemPush>=pMemEnd)  pMemPush = pMemStart;
        if(pMemPush==pMemPop) { res=0; pop(); }
    }
    return res;
}

#endif // PSINS_RMEMORY

/***************************  class CPTimer  *********************************/
CPTimer::CPTimer(double tMax0, int isStart0, int autoReset0)
{
    Start(tMax0);
    isStart = isStart0, autoReset = autoReset0;
}

void CPTimer::Start(double tMax0)
{
    overflow = 0;
    tCurrent = 0.0;
    isStart = 1;
    if(tMax0>1.0e-6) tMax=tMax0;
}

void CPTimer::Stop(void)
{
    isStart = 0;
}

BOOL CPTimer::Update(double tStep)
{
    overflow = 0;
    if(isStart)
    {
        tCurrent += tStep;
        if(tCurrent>=tMax)
        {
            overflow = 1;
            if(autoReset)
                tCurrent = 0.0;
        }
    }
    return overflow;
}

/***************************  function AlignCoarse  *********************************/
CVect3 Alignsb(CVect3 wmm, CVect3 vmm, double latitude)
{
    double T11, T12, T13, T21, T22, T23, T31, T32, T33;
    double cl = cos(latitude), tl = tan(latitude), nn;
    CVect3 wbib = wmm / norm(wmm),  fb = vmm / norm(vmm);
    T31 = fb.i,              T32 = fb.j,                T33 = fb.k;
    T21 = wbib.i/cl-T31*tl,  T22 = wbib.j/cl-T32*tl,    T23 = wbib.k/cl-T33*tl;     nn = sqrt(T21*T21+T22*T22+T23*T23);  T21 /= nn, T22 /= nn, T23 /= nn;
    T11 = T22*T33-T23*T32,   T12 = T23*T31-T21*T33,     T13 = T21*T32-T22*T31;      nn = sqrt(T11*T11+T12*T12+T13*T13);  T11 /= nn, T12 /= nn, T13 /= nn;
    CMat3 Cnb(T11, T12, T13, T21, T22, T23, T31, T32, T33);
    return m2att(Cnb);
}

CAligni0::CAligni0(const CVect3 &pos0, const CVect3 &vel0, int velAid0)
{
    eth.Update(pos0);
    this->vel0 = vel0, velAid = velAid0;
    tk = 0;
    t0 = t1 = 10, t2 = 0; 
    wmm = vmm = vib0 = vi0 = Pib01 = Pib02 = Pi01 = Pi02 = O31;
    qib0b = CQuat(1.0);
}

CQuat CAligni0::Update(const CVect3 *pwm, const CVect3 *pvm, int nSamples, double ts, const CVect3 &vel)
{
    double nts = nSamples*ts;
    imu.Update(pwm, pvm, nSamples);
    wmm = wmm + imu.phim;  vmm = vmm + imu.dvbm;
    // vtmp = qib0b * (vm + 1/2 * wm X vm)
    CVect3 vtmp = qib0b*imu.dvbm;
    // vtmp1 = qni0' * [dvn+(wnin+wnie)Xvn-gn] * ts;
    tk += nts;
    CMat3 Ci0n = pos2Cen(CVect3(eth.pos.i,eth.wie*tk,0.0));
    CVect3 vtmp1 = Ci0n*(-eth.gn*nts);
    if(velAid>0)
    {
        CVect3 dv=vel-vel0;  vel0 = vel;
        if(velAid==1)       vtmp1 += Ci0n*dv;               // for GPS vn-aided
        else if(velAid==2)  vtmp -= qib0b*dv+imu.phim*vel;  // for OD vb-aided
    }
    // Pib02 = Pib02 + vib0*ts, Pi02 = Pi02 + vi0*ts
    vib0 = vib0 + vtmp,      vi0 = vi0 + vtmp1;
    Pib02 = Pib02 + vib0*nts, Pi02 = Pi02 + vi0*nts;
    //
    if(++t2>3*t0)
    {
        t0 = t1, Pib01 = tmpPib0, Pi01 = tmpPi0;
    }
    else if(t2>2*t0 && t1==t0)
    {
        t1 = t2, tmpPib0 = Pib02, tmpPi0 = Pi02;
    }
    //
    qib0b = qib0b*rv2q(imu.phim);
    // qnb=qni0*qiib0*qib0b
    if(t2<100)
    {
        qnb0 = qnb = CQuat(1.0);
    }
    else if(t2<1000)
    {
        qnb0 = qnb = a2qua(Alignsb(wmm, vmm, eth.pos.i));
    }
    else
    {
        CQuat qi0ib0 = m2qua(dv2att(Pi01, Pi02, Pib01, Pib02));
        qnb0 = (~m2qua(pos2Cen(CVect3(eth.pos.i,0.0,0.0))))*qi0ib0;
        qnb = (~m2qua(Ci0n))*qi0ib0*qib0b;
    }
    return qnb;
}


/***************************  class CUartPP  *********************************/
#ifdef PSINS_UART_PUSH_POP

#pragma message("  PSINS_UART_PUSH_POP")

CUartPP::CUartPP(int frameLen0, unsigned short head0)
{
    popIdx = 0;
    pushIdx = 1;
    unsigned char *p = (unsigned char*)&head0;
    head[0] = p[1], head[1] = p[0];
//  *(unsigned short*)head = head0;
    frameLen = frameLen0;
    csflag = 1, cs0 = 4, cs1 = frameLen-1, css = 2;
    overflow = getframe = 0;
}

BOOL CUartPP::checksum(const unsigned char *pc)
{
    chksum = 0;
    if(csflag==0) return 1;
    for(int i=cs0; i<=cs1; i++)
    {
        if(csflag==1)   chksum += pc[i];
        else if(csflag==2)  chksum ^= pc[i];
    }
    return chksum==pc[css];
}

int CUartPP::nextIdx(int idx)
{
    return idx >= UARTBUFLEN - 1 ? 0 : idx + 1;
}

int CUartPP::push(const unsigned char *buf0, int len)
{
    int overflowtmp = 0;
    for (int i=0; i<len; i++)
    {
        buf[pushIdx] = buf0[i];
        pushIdx = nextIdx(pushIdx);
        if (pushIdx==popIdx) {
            overflow++;
            overflowtmp++;
            popIdx = nextIdx(popIdx);
        }
    }
    return overflowtmp;
}

int CUartPP::pop(unsigned char *buf0)
{
    int getframetmp = 0;
    while(1)
    {
        if((pushIdx>popIdx&&popIdx+frameLen<pushIdx) || (pushIdx<popIdx&&popIdx+frameLen<pushIdx+UARTBUFLEN))
        {
            int popIdx1 = nextIdx(popIdx);
            if(buf[popIdx]==head[0] && buf[popIdx1]==head[1])
            {
                if(!buf0) buf0=&popbuf[0];
                for (int i=0; i<frameLen; i++)
                {
                    buf0[i] = buf[popIdx];
                    popIdx = nextIdx(popIdx);
                }
                if(checksum(buf0))  // checksum
                {
                    getframe++;
                    getframetmp = 1;
                    break;
                }
                else
                    popIdx = popIdx1; 
            }
            else
            {
                popIdx = popIdx1; 
            }
        }
        else
            break;
    }
    return getframetmp;
}
#endif

#ifdef PSINS_psinsassert

#pragma message("  psinsassert();")

BOOL psinsassert(BOOL b)
{
    int res;

    if(b)
    {
        res = 1;
    }
    else
    {
        res = 0;
    }
    return res;
}

#endif

double diffYaw(double yaw, double yaw0)
{
    double dyaw = yaw-yaw0;
    if(dyaw>=PI) dyaw-=_2PI;
    else if(dyaw<=-PI) dyaw+=_2PI;
    return dyaw;
}

double MKQt(double sR, double tau)
{
    return sR*sR*2.0/tau;
}

CVect3 MKQt(const CVect3 &sR, const CVect3 &tau)
{
    return CVect3(sR.i*sR.i*2.0/tau.i, sR.j*sR.j*2.0/tau.j, sR.k*sR.k*2.0/tau.k);
}
