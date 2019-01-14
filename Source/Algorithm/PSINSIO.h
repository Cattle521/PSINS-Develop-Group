#ifndef _PSINSIO_H
#define _PSINSIO_H

#include "PSINS.h"
#include <stdio.h>

class CFileRdWt;

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

#endif /* ifndef _PSINSIO_H */
