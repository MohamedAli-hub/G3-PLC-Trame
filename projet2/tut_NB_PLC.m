clear;
channel_flag = 1;   % 0 for no channel
AWGN = 0;   % 1 for AWGN
warning off
clc
NFFT=256;
NSC=72;
Nbitpersymbol=72;
Nsym=77;
NCP=30;
N_itr=500;
N_frame=1;
j=sqrt(-1);
EbN0dB=10;%-5:10;
EsN0dB = EbN0dB + 10*log10(NSC/NFFT) + 10*log10(NFFT/(NFFT+NCP)); % converting to symbol to noise ratio
snr=EsN0dB- 10*log10(NFFT/(NFFT+NCP));
%% PREMABLE GENERATION
%% SYNCP
phic=pi/8;
for  k=1:15
    phi(k)=exp(j*k*phic);
end
SYNCP_f=[zeros(1,33) phi(2) phi(1) phi(1) 1 1 phi(15) phi(14) phi(12) phi(11) phi(9) phi(7) phi(4) phi(1) phi(15) phi(12) phi(9) phi(5) phi(1) phi(14) phi(10) phi(5) 1 phi(12) phi(6) phi(1) phi(12) phi(6) 1 phi(10) phi(3) phi(13) phi(6) phi(15) phi(7) 1 phi(8) 1 phi(8) phi(15) phi(6) phi(14) phi(4) phi(11) phi(2) phi(8) phi(14) phi(3) phi(9) phi(15) phi(3) phi(8) phi(13) phi(1) phi(5) phi(9) phi(13) phi(1) phi(4) phi(7) phi(10) phi(13) phi(15) phi(1) phi(3) phi(4) phi(5) phi(7) phi(7) phi(8) phi(9) phi(10) phi(10) zeros(1,151)];
SYNCP=ifft(SYNCP_f);
EST_SYNCP= sum(abs(SYNCP).^2)/length(SYNCP);
SYNCP=SYNCP./sqrt(EST_SYNCP);% Power Normalization

%% SYNCM
phic1=(pi/8)+pi;
for jj=1:15
    phi1(jj)=exp(j*jj*phic1);
end
SYNCM_f=[zeros(1,33) phi1(2) phi1(1) phi1(1) 1 1 phi1(15) phi1(14) phi1(12) phi1(11) phi1(9) phi1(7) phi1(4) phi1(1) phi1(15) phi1(12) phi1(9) phi1(5) phi1(1) phi1(14) phi1(10) phi1(5) 1 phi1(12) phi1(6) phi1(1) phi1(12) phi1(6) 1 phi1(10) phi1(3) phi1(13) phi1(6) phi1(15) phi1(7) 1 phi1(8) 1 phi1(8) phi1(15) phi1(6) phi1(14) phi1(4) phi1(11) phi1(2) phi1(8) phi1(14) phi1(3) phi1(9) phi1(15) phi1(3) phi1(8) phi1(13) phi1(1) phi1(5) phi1(9) phi1(13) phi1(1) phi1(4) phi1(7) phi1(10) phi1(13) phi1(15) phi1(1) phi1(3) phi1(4) phi1(5) phi1(7) phi1(7) phi1(8) phi1(9) phi1(10) phi1(10) zeros(1,151)];
SYNCM=ifft(SYNCM_f);
EST_SYNCM= sum(abs(SYNCM).^2)/length(SYNCM);
SYNCM=SYNCM./sqrt(EST_SYNCM);% Power Normalization
%% Preamble generation
preamble=[repmat(SYNCP,1,8) SYNCM SYNCM(1:128) ];
for ii=1:length(EbN0dB)
    ChMSE_LMMSE=0;
    No_of_errors_plc = 0;
    for itr=1:N_itr
        for frame=1:N_frame
            %% DATA GENERATION
            %% BPSK
            ip = rand(1,Nbitpersymbol*Nsym) > 0.5; % random 1's and 0's
            ip_BPSK = 2*ip-1; % BPSK modulation 0 --> -1, 1 --> +1
            ip_BPSK = reshape(ip_BPSK,Nbitpersymbol,Nsym).'; % grouping into multiple symbols,bit mapper
            X_Freq = [zeros(Nsym,33) ip_BPSK zeros(Nsym,151)] ;
            
            X_temps =sqrt(NFFT)*ifft((X_Freq.')).';
            
            % Appending cyclic prefix
            X_CP = [X_temps(:,[NFFT-NCP+1:NFFT]) X_temps];
            transmit_data=reshape(X_CP.',Nsym *286,1).';
            %x_data1=sqrt((nFFT+nG)/nFFT)*x_data;
            %% Packet generation
            %transmit_data(1:277)
            signal_oneframe=[transmit_data(1:256) preamble transmit_data];
            
            
            
            
            
            %% CHANNEL
            N=10; %% path number
            k=1;
            a0=1e-3;
            a1=2.5e-9;
            g=randn(1,N);
            d=1000+400*randn(1,N);
            ff=200;
            vp=3e8/4;
            H=zeros(1,length(ff));
            for m=1:N
                H(m)=g(m).*exp(-(a0+a1.*((ff.*1e3).^k)).*d(m)).*exp(-2i.*pi.*(ff.*1e3).*(d(m)./vp)) ;
                
            end
            % H0(ff) = sum((H(ff,1:20))');
            h=ifft(H);
            h=h./norm(h);
            H_plc=fft(h,NFFT);
            %%
            
            if channel_flag == 0
                x_received=signal_oneframe;%conv(h,signal_oneframe); %%%%%% Received signal
            else
                x_received=conv(h,signal_oneframe); %%%%%% Received signal
            end
            %% noise
            
            
            
            % AWGN
            yt_PLC=awgn(x_received,EbN0dB(ii),'measured'); %%% With AWGN noise
            
            
            
%             
%             %%  Cyclostationary noise generation
%             load ('LV14Params.mat');
%             % Mains period, sampling rate, Number of samples
%             Tac=1/(2*50); fs=1.25E6;
%             SampleNum=fs*Tac; % Number of sample for a Period
%             % Filter for period 1
%             filter1.InputNum= (10/8.3)*(timeRange1(2)-timeRange1(1))*fs*1E-3;
%             filter1.n=[0.25796067 0.515921 0.071 0.2579];
%             filter1.d=[1.000 -0.757914 0.275010];
%             % Filter for period 2
%             filter2.InputNum= (10/8.3)*(timeRange2(2)-timeRange2(1))*fs*1E-3;
%             filter2.n=[0.5 0.5];
%             filter2.d=[1,-0.89968];
%             % Filter for period 3
%             filter3.InputNum= (10/8.3)*(timeRange3(2)-timeRange3(1))*fs*1E-3;
%             filter3.n=[1 -0.092376 -0.63353*sqrt(-1)];
%             filter3.d=[1,-0.575198];
%             %------------------
%             Noise_sim=[];nIter=2;
%             for i =1:nIter
%                 noise =(randn(1,SampleNum)+1i*randn(1,SampleNum))/sqrt(2);
%                 NoiseT1=noise(1:filter1.InputNum)*sqrt(1.1E-4);
%                 FilteredNoise1=filter(filter1.n,filter1.d,NoiseT1); %-For 1st time range: thefiltered noise
%                 NoiseT2=noise(filter1.InputNum +1 : filter1.InputNum + filter2.InputNum)*sqrt(6.6E-4);
%                 FilteredNoise2=filter(filter2.n,filter2.d,NoiseT2); %-For 2nd time range: thefiltered noise
%                 NoiseT3=noise(end-filter3.InputNum+1: end)*sqrt(6E-3);
%                 FilteredNoise3=filter(filter3.n,filter3.d,NoiseT3); %-For 3rd time range: thefiltered noise
%                 Noise_1_Period=[FilteredNoise1,FilteredNoise2,FilteredNoise3];
%                 % Noise_sim=ifft(Noise_1_Period,SampleNum);
%                 Noise_sim=[Noise_sim,Noise_1_Period];
%             end
%             %% Normalization
%   
%             noise_power = var(Noise_sim(1:length(x_received)));
%             SNR=10.^(EbN0dB(ii)/10);
%             %%
%             if AWGN == 1
%                 %% test AWGN
%                 noise =(randn(1,length(x_received))+1i*randn(1,length(x_received)))/sqrt(2);
%                 yt_plc_cyc=x_received+noise/sqrt(SNR); %% TO VALIDATE THE WORK WITH THEORETICAL CASE
%             else
%                 
%                 yt_plc_cyc=x_received+(((Noise_sim(1:length(x_received))))/sqrt(noise_power))/sqrt(SNR);
%             end
% 
%             
%             %% Channel estimation
%             H_orig=H_plc(:,[34:105])';
%             beta=1;
%             % LS
%             prechannel=0;
%             for n=2:9
%                 prechannel=prechannel+yt_plc_cyc((n-1)*NFFT+1:n*NFFT);
%             end
%             prechannel=prechannel/8;
%             prechannel=prechannel*sqrt(EST_SYNCP);
%             averageinfreq=fft(prechannel,NFFT);
%             channel_est=averageinfreq(:,[34:105])./SYNCP_f(:,[34:105]);
%             % MMSE
%             SNR=10.^(EbN0dB/10);
%             Rhh = H_orig*H_orig';
%             W = Rhh*(inv(Rhh+(beta/SNR(ii))*eye(size(Rhh))));
%             HhatLMMSE = W*channel_est';
%             ChMSE_LMMSE = ChMSE_LMMSE + ((H_orig -HhatLMMSE)'*(H_orig-HhatLMMSE))/72;
%             
%             
%             
%             
%             %% demodulation
%             
%             
%             if channel_flag == 0
%                 yt_Data=yt_plc_cyc(2689:length(yt_plc_cyc));%-(length(h)-1));
%             else
%                 yt_Data=yt_plc_cyc(2689:length(yt_plc_cyc)-(length(h)-1));
%             end
%             
%             % %yt_data=yt(pp:length(yt)-(nTap-1));
%             yt_PLC = reshape(yt_Data.',286,Nsym).'; % formatting the received vector into symbols
            yt_rem_cp_plc= yt_PLC(:,[31:286]); % removing cyclic prefix
            yt_cp=yt_PLC(:,[1:30]);
            sigma=mean(abs(yt_cp(1,:)).^2);
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %                       FFT BLOCK
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            yF_plc = (1/sqrt(NFFT))*(fft(yt_rem_cp_plc.')).';
            
            
            %%
            %%%%%%%%%%%%%%%%%%% Trials for Noise PSD estimation
            
            %% taking symbols in frequency domain
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dta=reshape(yF_plc.',NFFT*Nsym,1).';
            %%%%%%%%  R1  %%%%%%%%%%%%%%
            region1_bound=10030;
            index1=((2688)+1):region1_bound;
            y1=dta(index1);
            %%%%%%%%  R2  %%%%%%%%%%%%%%%
            region2_bound=region1_bound+2078;
            index2=((region1_bound)+1):region2_bound;
            y2=dta(index2);
            %%%%%%%   R3  %%%%%%%%%%%%%%%%%%
            region3_bound=region2_bound+391;
            index3=((region2_bound)+1):region3_bound;
            y3=dta(index3);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PSD estimation
            Mtx=zeros(25,256);
            %  Mtx=[];
            for t=1:25
                psd1=(mean((abs(y1((t-1)*NFFT+1:t*NFFT))).^2))-(((abs(H_plc)).^2));
                % Mtx(t,:)=psd1(t);
            end
            % Mtx;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            PSd1 = mean(abs((y1).^2))-((abs((H_plc).^2)));
            B=(1/10030)*ones(10030,1);
            out=filter(B,1,PSd1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            PSd2 = mean(abs((y2).^2))-((abs((H_plc).^2)));
            C=(1/2078)*ones(2078,1);
            out1=filter(C,1,PSd2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            PSd3 = mean(abs((y3).^2))-((abs((H_plc).^2)));
            D=(1/391)*ones(391,1);
            out2=filter(D,1,PSd2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Preamble_freq=fft(yt_plc_cyc(257:2688)*sqrt(EST_SYNCP));
            
            Y=[Preamble_freq dta];
            
            %% Region 1
            for l=1:34
                psd_plc_est_r1(l)=(mean((abs(Y((l-1)*NFFT+1:l*NFFT))).^2))-(mean((abs(H_plc)).^2));
            end
            var_r1=mean(psd_plc_est_r1);
            %% Region 2
            for k=35:41
                psd_plc_est_r2(k)=(mean((abs(Y((k-1)*NFFT+1:k*NFFT))).^2))-(mean((abs(H_plc)).^2));
            end
            var_r2=mean(psd_plc_est_r2);
            %% Region 3
            for m=42:43
                psd_plc_est_r3(m)=(mean((abs(Y((m-1)*NFFT+1:m*NFFT))).^2))-(mean((abs(H_plc)).^2));
            end
            var_r3=mean(psd_plc_est_r3);
            
            %%
            %%%%%%%%%%%%%%%%%% End of trials
            
            
            %%
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %                   Selecting 72 Data Bits from array
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            yF_plc1=yF_plc(:,[34:105]);
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %    Equalization by the estimated channel frequency response in data
            %    subcarriers
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            for t=1:Nsym
                if channel_flag == 0
                    yF_plc1(t,:) = yF_plc1(t,:);%./channel_est(1,:);
                else
                    yF_plc1(t,:) = yF_plc1(t,:)./channel_est(1,:);
                end
            end
            
            yMod_plc = yF_plc1;
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %                       BPSK Demodulation
            
            % BPSK demodulation for Rayleigh
            % +ve value --> 1, -ve value --> -1
            
            ipModHat_plc = 2*floor(real(yMod_plc/2)) + 1;
            ipModHat_plc(find(ipModHat_plc>1)) = +1;
            ipModHat_plc(find(ipModHat_plc<-1)) = -1;
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %Converting the modulated values into bits and rehaping to count errors
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % converting modulated values into bits
            ipBitHat_plc = (ipModHat_plc+1)/2;
            ipBitHat_plc = reshape(ipBitHat_plc.',Nbitpersymbol*Nsym,1).';
            
            % counting the errors
            nErr_plc = size(find(ipBitHat_plc - ip),2);
            
            No_of_errors_plc = No_of_errors_plc + nErr_plc;
        end
    end
    T_Errors_plc(ii) = No_of_errors_plc;
end
simBer_plc= T_Errors_plc./(Nsym*Nbitpersymbol*N_itr);


%% Figures

semilogy(EbN0dB,simBer_plc,'r*-','LineWidth',2);

xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
%axis([0 25 10^-6 1])
grid on

