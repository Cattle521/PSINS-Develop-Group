#include "KFApp.h"

#define FRQ 100
#define Ts (1.0/FRQ)

struct SSensor { double t; CVect3 Acc, Gyro, MtiAtt; } *ps;
void convert(void)
{
	ps->t *= 0.01; ps->Gyro*=glv.dps*Ts; ps->Acc*=Ts;
	ps->MtiAtt*=DEG; ps->MtiAtt=CVect3(ps->MtiAtt.i, ps->MtiAtt.j, C360toCC180(ps->MtiAtt.k));
}

int main()
{
	CKFApp kf;

	CFileRdWt::Dir("./");
	CFileRdWt fin("MTI_gaotie.txt", 10);  ps=(SSensor*)fin.buff;
	CFileRdWt fins("ins.bin");
    CFileRdWt fkf("kf.bin"); 

	fin.load(1); convert();
	kf.Init(CSINS(a2qua(0,0,ps->MtiAtt.k*DEG), O31, O31)); 
	
	for(int i=0; i<10000*FRQ; i+=1)
	{
		if(!fin.load()) break;  convert();
		kf.SetMeasVG();
		kf.Update(&ps->Gyro, &ps->Acc, 1, Ts);
		if(i%5==0)
		{
			fins<<kf.sins<<ps->MtiAtt;
			fkf<<kf;
		}
		disp(i, 100, FRQ);
	}
}
