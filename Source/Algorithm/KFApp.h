/**
 * @file KFApp.h
 * @brief KFApp c++ hearder file
 * @author Yan Gongmin,College of Automation,
 *      Northwestern Polytechnical University,Xi'an 710072, China
 * @version 1.0
 * @date 2017-04-29
 */
#ifndef _KFAPP_H
#define _KFAPP_H

#include "PSINS.h"

class CKFApp:public CSINSTDKF
{
public:
    CVect3 vbINS, vbOD;
    CMat3 Cbo;              ///< from body-frame to OD-frame
    double tODInt;          ///<
    BOOL measGPSvnValid;    ///<
    BOOL measGPSposValid;   ///<
    BOOL measODvbValid;     ///<
    BOOL measVGValid;       ///<
    BOOL measMAGyawValid;   ///<
    CRAvar RGV;             ///<
    CRAvar RGP;             ///<
    CRAvar RVG;             ///<

    CKFApp(void);
    virtual void Init(const CSINS &sins0);
    virtual void SetMeas(void);
    void SetMeasGPS(const CVect3 &pgps=O31, const CVect3 &vgps=O31);
    void SetMeasOD(double dSod, double ts);
    void SetMeasVG(void);
    void SetMeasYaw(double ymag);
    int Update(const CVect3 *pwm, const CVect3 *pvm, int nn, double ts);
};

#endif

