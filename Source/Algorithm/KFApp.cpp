#include "KFApp.h"

CKFApp::CKFApp(void):CSINSTDKF(15, 13)
{
	vbINS = vbOD = O31;  tODInt = 0.0;
	Cbo = I33;
	Hk(9,3) = Hk(10,4) = Hk(11,5) = 1.0;
	Hk(12,2) = 1.0;
	measGPSvnValid = measGPSposValid = measODvbValid = measVGValid = measMAGyawValid = 0;
}

void CKFApp::Init(const CSINS &sins0)
{
	CSINSKF::Init(sins0, 1);
	Pmax.Set2(100.0*glv.deg,100.0*glv.deg,100.0*glv.deg,    50.0,50.0,50.0,    1.0e6/glv.Re,1.0e6/glv.Re,1.0e6, 
		5000.0*glv.dph,5000.0*glv.dph,5000.0*glv.dph,    50.0*glv.mg,50.0*glv.mg,50.0*glv.mg);
	Pmin.Set2(3.*glv.min,3.*glv.min,3.0*glv.min,    0.01,0.01,0.1,    1.0/glv.Re,1.0/glv.Re,0.1, 
		1.0*glv.dph,1.0*glv.dph,1.0*glv.dph,    150.0*glv.ug,150.0*glv.ug,150.0*glv.ug);
	Pk.SetDiag2(100.0*glv.deg,100.0*glv.deg,100.0*glv.deg,    10.0,10.0,10.0,     100.0/glv.Re,100.0/glv.Re,100.0, 
		100.0*glv.dph,101.0*glv.dph,102.0*glv.dph,    1.0*glv.mg,1.01*glv.mg,10.0*glv.mg);
	Qt.Set2(1.0*glv.dpsh,1.0*glv.dpsh,1.0*glv.dpsh,    10.0*glv.ugpsHz,10.0*glv.ugpsHz,10.0*glv.ugpsHz,    0.0,0.0,0.0,
		0.0*glv.dphpsh,0.0*glv.dphpsh,0.0*glv.dphpsh,    0.0*glv.ugpsh,0.0*glv.ugpsh,0.0*glv.ugpsh);
	Xmax.Set(INF,INF,INF, INF,INF,INF, INF,INF,INF,  3600.0*glv.dps,3600.0*glv.dps,3600.0*glv.dps,  10.0*glv.mg,10.0*glv.mg,50.0*glv.mg);
	FBTau.Set(.20,.20,1.0,     1.0,1.0,1.0,     1.0,1.0,1.0,    10.0,10.0,10.0,    10.0,10.0,10.0);
	Rt.Set2(0.1,0.1,0.3,   10.0/glv.Re,10.0/glv.Re,30.0, 0.1,10.0,0.1, .10,.10,.10, 1.0*glv.deg);
	Rmax = Rt*10.0;  Rmin = Rt*0.01;  
	Rb = 0.9;
}

void CKFApp::SetMeas(void)
{
	if(measGPSvnValid)	{ SetMeasFlag(000007);	}
	if(measGPSposValid)	{ SetMeasFlag(000070);	}
	if(measODvbValid)	{ SetMeasFlag(000700);	}
	if(measVGValid)		{ SetMeasFlag(007000);	}
	if(measMAGyawValid)	{ SetMeasFlag(010000);	}
	measGPSvnValid = measGPSposValid = measODvbValid = measVGValid = measMAGyawValid = 0;
}

void CKFApp::SetMeasGPS(const CVect3 &pgps, const CVect3 &vgps)
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

void CKFApp::SetMeasOD(double dSod, double ts)
{
	if(fabs(sins.wnb.k)<5.0*glv.dps)
	{
		tODInt += ts;
		vbINS += sins.vb*ts;
		vbOD += Cbo*CVect3(0,dSod,0);
		if(tODInt>1.0)
		{
			*(CVect3*)&Zk.dd[6] = (vbINS - vbOD)/tODInt;
			*(CVect3*)&rts.dd[6] = tODInt;
			Hk.SetMat3(6, 3, sins.Cbn);
			vbINS = vbOD = O31;	tODInt = 0.0;
			measODvbValid = 1;
		}
	}
	else
	{
		vbINS = vbOD = O31;	tODInt = 0.0;
	}
}

void CKFApp::SetMeasVG(void)
{
	*(CVect3*)&Zk.dd[9] = sins.vn;
	*(CVect3*)&rts.dd[9] = 0.01;
	measVGValid = 1;
}

void CKFApp::SetMeasYaw(double ymag)
{
	if(ymag>EPS || ymag<-EPS)  // !IsZero(yawGPS)
		if(sins.att.i<30*DEG && sins.att.i>-30*DEG)
		{
			*(CVect3*)&Zk.dd[12] = -diffYaw(sins.att.k, ymag);
			measMAGyawValid = 1;
		}
}

int CKFApp::Update(const CVect3 *pwm, const CVect3 *pvm, int nn, double ts)
{
	int res=TDUpdate(pwm, pvm, nn, ts, 5);  
	return res;
}
