#include <iostream>
#include <complex>
#include <math.h>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <iterator>

using namespace std;
using Complex = complex<double>;

#define MAX 200
#define PI  3.14159265358979323846264338327950288419716939937510    //Pi, 50 decimal places
#define NFFT 256.0
#define NSC 72.0
#define Nbitpersymbol 72
#define Nsym 77
#define NCP 30.0
#define N_itr 500
#define N_frame 1
#define EbN0dB 10.0

void transpose(Complex M_t[],Complex *M);
void print(const char * prompt, Complex A[], int N);
void FFT(Complex f[], Complex ftilde[], int N, int step = 1);
void DFT(Complex f[], Complex ftilde[], int N, int step = 1);
void iFFT(Complex ftilde[], Complex f[], int N);




int main()
{

    double EsN0dB = EbN0dB + 10*log10(NSC/NFFT) + 10*log10(NFFT/(NFFT+NCP)); // converting to symbol to noise ratio;
    double snr= EsN0dB- 10*log10(NFFT/(NFFT+NCP));
    double EST_SYNCP=0;
    double EST_SYNCM=0;
    int N;
    Complex ip_BPSK[Nsym*Nbitpersymbol]  ;
    Complex ip_BPSK1[Nsym][Nbitpersymbol];
    std::complex<double> Phi1[15];
    std::complex<double> SYNCP_f[256]={{0,0}};
    double cosinus,sinus;
    Complex Preamble[2432];//Preamble Contains 8 Symbol SYNCP and 1,5 Symbols SYNCM
    /*-----------------------Preamble Generation----------------------------------------------------------*/
    /*------------SYNCP-----------------------------------------------------------------------------------*/

    for(int k=1; k< 16;k++){
        cosinus = cos(k*PI/8);
        sinus = sin(k*PI/8);
        if (std::abs(cosinus) < 0.0001) cosinus= 0;
        if (std::abs(sinus) < 0.0001) sinus = {0};
        Phi1[k]={{cosinus,sinus}};
        //std::cout<<"real "<<setprecision(4)<<real(Phi1[k]) <<" image "<<imag(Phi1[k])<<endl;
    }
    std::complex<double> Phi1_f[72]={Phi1[2],Phi1[1],Phi1[1],{1,0},{1,0},Phi1[15],Phi1[14],Phi1[12],Phi1[11],Phi1[9],Phi1[7],Phi1[4],Phi1[1],Phi1[15],Phi1[12],Phi1[9],Phi1[5],Phi1[1],Phi1[14],Phi1[10],Phi1[5],{1,0},Phi1[12],Phi1[6],Phi1[1],Phi1[12],Phi1[6],{1,0},Phi1[10],Phi1[3],Phi1[13],Phi1[6],Phi1[15],Phi1[7],{1,0},Phi1[8],{1,0},Phi1[8],Phi1[15],Phi1[6],Phi1[14],Phi1[4],Phi1[11],Phi1[2],Phi1[8],Phi1[14],Phi1[3],Phi1[9],Phi1[15],Phi1[3],Phi1[8],Phi1[13],Phi1[1],Phi1[5],Phi1[9],Phi1[13],Phi1[1],Phi1[4],Phi1[7],Phi1[10],Phi1[13],Phi1[15],Phi1[1],Phi1[3],Phi1[4],Phi1[5],Phi1[7],Phi1[7],Phi1[8],Phi1[9],Phi1[10],Phi1[10]};
    for(int i=33;i<105;i++){
        SYNCP_f[i]=Phi1_f[i-33];
    }
    /*for(int i=0;i<256;i++){
        std::cout<<setprecision(4)<<SYNCP_f[i]<<endl;
    }*/
    N= sizeof(SYNCP_f)/sizeof(SYNCP_f[0]);
    Complex * SYNCP = new Complex[N];
    iFFT(SYNCP_f,SYNCP, N);//

    for(int i=0;i<N;i++){
        EST_SYNCP+= pow(sqrt(norm(SYNCP[i])),2);
    }
    EST_SYNCP/=N;
    for(int i=0;i<N;i++){
           SYNCP[i]/=sqrt(EST_SYNCP);
    }


    //print("SYNCP",SYNCP,N);

    /*------------SYNCM-----------------------------------------------------------------------------------
    ------------------------------------------------------------------------------------------------------*/

    std::complex<double> Phi2[15];
    std::complex<double> SYNCM_f[256]={{0,0}};
    for(int k=1; k< 16;k++){
        cosinus = cos(k*(PI+PI/8));
        sinus = sin(k*(PI+PI/8));
        if (std::abs(cosinus) < 0.0001) cosinus= 0;
        if (std::abs(sinus) < 0.0001) sinus = {0};
        Phi2[k]={{cosinus,sinus}};
        //std::cout<<"real "<<setprecision(4)<<real(Phi2[k]) <<" image "<<imag(Phi2[k])<<endl;
    }
    std::complex<double> Phi2_f[72]={Phi2[2],Phi2[1],Phi2[1],{1,0},{1,0},Phi2[15],Phi2[14],Phi2[12],Phi2[11],Phi2[9],Phi2[7],Phi2[4],Phi2[1],Phi2[15],Phi2[12],Phi2[9],Phi2[5],Phi2[1],Phi2[14],Phi2[10],Phi2[5],{1,0},Phi2[12],Phi2[6],Phi2[1],Phi2[12],Phi2[6],{1,0},Phi2[10],Phi2[3],Phi2[13],Phi2[6],Phi2[15],Phi2[7],{1,0},Phi2[8],{1,0},Phi2[8],Phi2[15],Phi2[6],Phi2[14],Phi2[4],Phi2[11],Phi2[2],Phi2[8],Phi2[14],Phi2[3],Phi2[9],Phi2[15],Phi2[3],Phi2[8],Phi2[13],Phi2[1],Phi2[5],Phi2[9],Phi2[13],Phi2[1],Phi2[4],Phi2[7],Phi2[10],Phi2[13],Phi2[15],Phi2[1],Phi2[3],Phi2[4],Phi2[5],Phi2[7],Phi2[7],Phi2[8],Phi2[9],Phi2[10],Phi2[10]};
    for(int i=33;i<105;i++){
        SYNCM_f[i]=Phi2_f[i-33];
    }
    N= sizeof(SYNCM_f)/sizeof(SYNCM_f[0]);
    Complex * SYNCM = new Complex[N];
    iFFT(SYNCM_f,SYNCM, N);//
    //print("SYNCM",SYNCM,N);
    for(int i=0;i<N;i++){
        EST_SYNCM+= pow(sqrt(norm(SYNCM[i])),2);
    }
    EST_SYNCM/=N;
    for(int i=0;i<N;i++){
           SYNCM[i]/=sqrt(EST_SYNCM);
    }

 /*--------Preamble Generation-----------------------------------------
 ----------------------------------------------------------------------*/
    for(int i=0;i<(int)N;i++){
        Preamble[i]=Preamble[i+N]=Preamble[i+2*N]=Preamble[i+3*N]=Preamble[i+4*N]=Preamble[i+5*N]=Preamble[i +6*N]=Preamble[i+7*N]=SYNCP[i];
        Preamble[i+8*N]=SYNCM[i];
        if(i<128)Preamble[i+9*N]=SYNCM[i];
    }

    /*--------*/

    for(int k=0; k<(sizeof(Preamble)/sizeof(Preamble[0]));k++){
           /* cout<<k<<endl;
            cout<<Preamble[k]<<endl; */

    }

    for(int i=0;i< Nsym*Nbitpersymbol ;i++){
       ip_BPSK[i] = {{2*(rand() %2 )-1,0}};

   }
   //cout<<sizeof(ip_BPSK)<<endl ;
    for(int i=0;i<Nsym ;i++){
        for(int j=0;j<Nbitpersymbol  ;j++){
            ip_BPSK1[i][j]=ip_BPSK[j+Nbitpersymbol*i];
            cout<<ip_BPSK1[i][j]<<endl ;
        }
    }
    cout<<sizeof(ip_BPSK1)/sizeof(ip_BPSK1[0])<<endl ;
    cout<<sizeof(ip_BPSK1[0])/sizeof(ip_BPSK1[0][0])<<endl ;
    Complex M_t[Nbitpersymbol][Nsym];
    //Complex *T=ip_BPSK1;
    //transpose(M_t,T);
    return 0;
}


 //======================================================================

void print(const char * prompt, Complex A[], int N)
{
	cout << prompt << '\n' << fixed;
	for (int i = 0; i < N; i++) cout << A[i] << '\n';
}

//===========FFT Function===========================================================

void FFT(Complex f[], Complex ftilde[], int N, int step)           // Fast Fourier Transform
{
	if (N > 2 && N % 2 == 0)                                        // Recursion if a multiple of 2
	{
		int N2 = N >> 1;
		Complex* g = new Complex[N];

		int step2 = step << 1;
		FFT(f, g, N2, step2);
		FFT(f + step, g + N2, N2, step2);

		Complex w = polar(1.0, -2.0 * PI / N);
		Complex wn = 1.0;
		for (int n = 0; n < N2; n++)
		{
			if (n > 0) wn *= w;
			ftilde[n] = g[n] + wn * g[n + N2];
			ftilde[n + N2] = g[n] - wn * g[n + N2];
		}

		delete[] g;
	}
	else                                                              // Otherwise just do a DFT
	{
		DFT(f, ftilde, N, step);
	}
}

//==================DFT Function====================================================

void DFT(Complex f[], Complex ftilde[], int N, int step)           // Basic Discrete Fourier Transform
{
	Complex w = polar(1.0, -2.0 * PI / N);
	Complex wn = 1.0;
	for (int n = 0; n < N; n++)
	{
		if (n > 0) wn *= w;
		Complex wm = 1.0;
		ftilde[n] = 0.0;
		for (int m = 0, mpos = 0; m < N; m++, mpos += step)
		{
			if (m > 0) wm *= wn;
			ftilde[n] += f[mpos] * wm;
		}
	}
}

//=================iFFT Function=====================================================

void iFFT(Complex ftilde[], Complex f[], int N)                    // Inverse Fast Fourier Transform
{
	Complex* ftildeConjugate = new Complex[N];
	for (int m = 0; m < N; m++) ftildeConjugate[m] = conj(ftilde[m]);

	FFT(ftildeConjugate, f, N);

	double factor = 1.0 / N;
	for (int m = 0; m < N; m++) f[m] = conj(f[m]) * factor;

	delete[] ftildeConjugate;
}

//=====================Transpose Function====================================

void transpose(Complex M_t[],Complex M[]){
        int rows=sizeof(M)/sizeof(M[0]);
        int coloms=sizeof(M_t)/sizeof(M_t[0]);

        for(int i=0;i<coloms;i++){
                for(int j=0;j<rows;j++){
                    M_t[i][j]=M[j][i];
                }

        }

}

