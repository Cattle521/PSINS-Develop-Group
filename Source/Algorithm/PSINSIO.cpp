/**
 * @file PSINSIO.cpp
 * @brief Input(Read File) and Output(Write File) related classes and functions
 * @author Yan Gongmin
 * @version 1.0
 * @date 2019-01-04
 */

#include "PSINSIO.h"

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
