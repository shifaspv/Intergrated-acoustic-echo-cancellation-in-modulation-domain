close all;
clear all;
clc;

delete 'clean.wav'; delete 'enhanced.wav'
%% impuls response generation
fs=8000;
mic=[2 1.5 .5];
n=8;
r1=0.6;
rm=[5 4 3];
src=[5 2 1];
h=rir(fs, mic, n, r1, rm, src);
L=length(h)
%input
[x,fs] = audioread('S_01_01.wav'); %Far-end signal
[v,fs] = audioread('S_01_02.wav');
%  x=[zeros(2000,1);x];
v=[zeros(8000,1);v(8000:end)];
EchoSignal=filter(h,1,x); %Echo signal
P_near=bandpower(v);
P_echo=bandpower(EchoSignal);
while (P_echo>=0.448*P_near)
    x=0.99*x;
    EchoSignal=filter(h,1,x);
    P_echo=bandpower(EchoSignal);
end

% EchoSignal=[zeros(1000,1);EchoSignal];
M=length(v);
N=length(EchoSignal);
v=[zeros(N-M,1);v];
load babble
[v1,noise]=addnoise(v,babble(1:length(v)),0.001);
% disp('echo in the microphone')
% sound(EchoSignal)
% pause
% disp('near end speeech')
% sound(v)
% pause
% disp('combined input')

y=v1+EchoSignal; %(Microphone Signal)
% sound(y)
% DELAY ESTIMATION
delay=finddelay(x,y);
x1=[zeros(ceil(delay),1);x];
DTbegin=1;
alpha=0.01;
T1=0.75;  %for DTD
T2=0.65;
%% AES based on soft decision
% default algorithm constants for MMSE estimator
qq.of=8;        % overlap factor = (fft length)/(frame increment)
qq.ti=4e-3;    % desired frame increment (16 ms)
% qq.ri=0;        % round ni to the nearest power of 2
qq.ta=0.396;    % Time const for smoothing SNR estimate = -tinc/log(0.98) from [1]
qq.gx=100;     % maximum posterior SNR = 30dB
qq.gn=0.1;        % min posterior SNR as a power ratio when estimating prior SNR [1]
qq.gz=0.001;    % min posterior SNR as a power ratio [0.001 = -30dB]
qq.xn=0;        % minimum prior SNR = -Inf dB
% qq.xb=1;        % bias compensation factor for prior SNR [1]
% qq.lg=1;        % use log-domain estimator by default
qq.ne=0;        % noise estimation: 0=min statistics, 1=MMSE [0]
% qq.bt=-1;       % suppress binary masking
% qq.mx=0;        % no input mixing
% qq.tf='g';      % output the gain time-frequency plane by default
% qq.rf=0;
qq.tn=0.9;     % smoothing constant for noise estimation [500 ms]
qq.le=0.15;    % VAD threshold; use -Inf to prevent updating
qq.tx=0.06;    % initial noise interval [60 ms]
% parameters for DTD
c=0.998; % parameters for coloration filter(echo estimation)
r=0.998;
alpha_q=0.89; % smoothing parameter for q
beta=1.35; % over estimation facor of echo
% eta_e=0; % echo smoothing parameter for echo PSD eq(8)
% eta=0; %  smoothing parameter combined PSD eq(11)
q=0; % initial value of q

ni=round(qq.ti*fs);    % frame increment in samples
tinc=ni/fs;         % true frame increment time
a=exp(-tinc/qq.ta);
a=0.92
% a=0.92;              % SNR smoothing coefficient
gx=qq.gx;           % max posterior SNR as a power ratio
gz=qq.gz;           % min posterior SNR as a power ratio
xn=qq.xn;           % floor for prior SNR, xi
gn1=max(qq.gn-1,0); % floor for posterior SNR when estimating prior SNR
le=qq.le; % VAD threshold
% xb=qq.xb;
% tf=qq.tf;
% rf=qq.rf;
nd=max(1,round(qq.tx/tinc)); % number of frames for initial noise estimate
an=exp(-tinc/qq.tn); % Noise spectrum smoothing coefficient

% calculate power spectrum in frames
no=round(qq.of);                  	% integer overlap factor
nf=ni*no;     % fft length
tic
w=sqrt(hamming(nf+1))'; w(end)=[];  % for now always use sqrt hamming window
w=w/sqrt(sum(w(1:ni:nf).^2));       % normalize to give overall gain of 1
X1=fft(enframe(x1,w,ni),nf,2);
Y1=fft(enframe(y,w,ni),nf,2);
V1=fft(enframe(v,w,ni),nf,2);
E1=fft(enframe(noise,w,ni),nf,2);
[a1,b]=size(Y1);
% alf=0.8;
% alpha_q=0.7; % smoothing parameter for q
xphase=angle(X1);
xspec=abs(X1);
yphase=angle(Y1);
yspec=abs(Y1);
vphase=angle(V1);
vspec=abs(V1);
ephase=angle(E1);
espec=abs(E1);
xspec1=xspec';
yspec1=yspec';
vspec1=vspec';
espec1=espec';
[u,r1]=size(xspec1);
ni2=2;
nf1=1;
w1=sqrt(hamming(nf1+1))'; w1(end)=[];  % for now always use sqrt hamming window
w1=w1/sqrt(sum(w1(1:ni:nf1).^2));       % normalize to give overall gain of 1
nf2=256;

for k=1:nf
    xms=enframe(xspec1(k,:),w1,ni2);
    xmsf=(fft(xms,nf2,2));
    yms=enframe(yspec1(k,:),w1,ni2);
    ymsf=(fft(yms,nf2,2));
    vms=enframe(vspec1(k,:),w1,ni2);
    vmsf=(fft(vms,nf2,2));
    ems=enframe(espec1(k,:),w1,ni2);
    emsf=(fft(ems,nf2,2));
    xyp=xmsf.*conj(xmsf);
    yyp=ymsf.*conj(ymsf);
    [c1,d]=size(yyp);
    dpi=0;   % noise estimate
    ndp=0;              % noise estimate based on ndp frames
    
    
    xu=1;                           % dummy unsmoothed SNR from previous frame
    %     nd=8;
    %     xn=0;
    %     gn1=0;
    %     a=0.98;
    %     gx=100;
    %     gz=0.001;
    %     le=0.15;
    %     an=0.9912;
    %     kk=sqrt(2*pi);
    if ndp<nd
        ndx=min(c1,nd-ndp);         % number of frames to use
        %         sq=squeeze(yp(k,:,:));
        dpi=ndp/(ndp+ndx)*dpi+sum(yyp(1:ndx,:),1)/(ndp+ndx);
        ndp=ndp+ndx;
    end
    g=zeros(c1,d);                % create space for gain matrix
    x=zeros(c1,d);                % create space for prior SNR
    dp=zeros(c1,d);               % create space for noise power spectrum estimate
    %     switch qq.lg
    %         case 0                      % use amplitude domain estimator from [1]
    clear C;clear R;clear lamda1;clear lamda_e
    flag=0;
    q=0;
    lamda=0;
    eta=0.3;
    eta_e=0;
    px=10^-30;
    pzd=10^-30;
    pd=10^-30;
    pz=10^-30;
    pze=10^-30;
    pe=10^-30;
    ps=10^-30;
    for z=1:c1
        ypi=yyp(z,:);
        if flag==0
            if z==1
                C=(1-c)*abs(ymsf(z,:).*conj(xmsf(z,:)))+1e-10;
                R=(1-r)*abs(xmsf(z,:).*conj(xmsf(z,:)))+1e-10;
            else
                
                C=c*C+(1-c)*abs(ymsf(z,:).*conj(xmsf(z,:)));
                R=r*R+(1-r)*abs(xmsf(z,:).*conj(xmsf(z,:)));
                
            end
            H1=C./R;
        end
        E_hat=(H1).*abs(xmsf(z,:));
        E_hat=beta*E_hat;
        %         if flag==0
        if z==1
            lamda_e=((E_hat).^2);
        else
            lamda_e=eta_e*lamda_e+(1-eta_e)*((E_hat).^2);
        end
        %         end
        %          echo(z,:)=lamda_e;
        lamda1=lamda_e+dpi;
        if z==1
            lamda=((lamda1));
        else
            lamda=eta*lamda+(1-eta)*((lamda1));
        end
        %          end
        %         % mmse gain
        gami=max(min(ypi./(lamda),gx),gz);     % gamma = posterior SNR
        xi=max(a*xu+(1-a)*max(gami-1,gn1),xn);  % prior SNR
        xir=xi./(1+xi);
        if sum(gami.*xir-log(1+xi))<le*nf2 % noise frame
            if sum(E_hat)==0
                dpi=dpi*an+(1-an)*ypi;
            end
        end
        
         % near end speech absent probability
        f=1./(1+xi).*exp(xir.*(ypi./(lamda)));
        I=gami>3;
        q=alpha_q*q+(1-alpha_q)*I;
        Pr=1./(1+q.*f);
        gi=xir.*(1-Pr);
        xu=gami.*gi.^2;         % unsmoothed prior SNR
        
        %         gi=ones(1,d)
        %         gi=((max((abs(ymsf(z,:)).^2-Y_hat.^2),0)./(abs(ymsf(z,:)).^2)).^0.5);
        S=ymsf(z,:).*gi;
        
        pzd=(1-alpha)*pzd+alpha*abs(ymsf(z,:)*(E_hat'));
        pd=(1-alpha)*pd+alpha*(abs(E_hat*E_hat'));
        pz=(1-alpha)*pz+alpha*(abs(ymsf(z,:)*ymsf(z,:)'));
        pze=(1-alpha)*pze+alpha*abs(ymsf(z,:)*(S)');
        pe=(1-alpha)*pe+alpha*(abs(S*S'));
        l=z;
        rzd=pzd/(sqrt(pz*pd));
        rze=pze/(sqrt(pz*pe));
        m1(l)=(rzd);
        m2(l)=(rze);
        
        if (l>DTbegin)
            if abs(m1(l))<T1
                if abs(m2(l))>T2
                    
                    flag=1;
                    f1(l)=1;
                    
                else
                    
                    flag=0;
                    f1(l)=0;
                end
            else
                flag=0;
                f1(l)=0;
            end
        end
        
        g(z,:)=gi;              % save gain for later
        %                 x(z,:)=xi;              % save prior SNR
        %                 xu=gami.*gi.^2;         % unsmoothed prior SNR
        
    end
   
    se=(ifft((ymsf.*g),nf2,2));     % inverse dft and apply output window
    sv=(ifft((vmsf.*g),nf2,2));
    esv=(ifft((emsf.*g),nf2,2));
    [y1,p]=overlapadd(se(:,1),w1,ni2);
    [v2,p]=overlapadd(sv(:,1),w1,ni2);
    echo_est=overlapadd(esv(:,1),w1,ni2);
    ss(k,:)=y1;
    svv(k,:)=v2;
    echos(k,:)=echo_est;
end
%  Y=ss;
%  ni=ni;
% % Yfi=vec2mat(Y,nf);
%    [c,np]=size(spec);
p_rec=ss.';
[n,m]=size(yphase);
se2=(ifft(p_rec(1:n,:).*exp(1i*yphase),nf,2));     % inverse dft and apply output window
[y2,p]=overlapadd(se2,w,ni);
toc
ss1=2.5*y2;
% pause(3)
disp('enhanced')
 sound(real(ss1))

for k=1:length(ss1)
    powerY(k) = abs(y(k))^2; %Power of Microphone
    powerE(k)=abs((ss1(k)))^2; %power of Error
end
L=1000;
for k=1:N-2*L
    %Echo Return Loss Enhancement
    ERLE(k)=10*log10(mean(powerY(k:k+L))/mean(powerE(k:k+L)));
end

 audiowrite('microphone_input.wav',y,fs)
 audiowrite('echo_noise_cancelled.wav',(ss1),fs)
figure
subplot(311)
spectrogram(v,hamming(nf),ni,nf,fs,'yaxis')
title('Clean Near-end')
subplot(312)
spectrogram(y,hamming(nf),ni,nf,fs,'yaxis')
title('Mixture at microphone')
subplot(313)
spectrogram(real(ss1),hamming(nf),ni,nf,fs,'yaxis')
title('Echo_noise cancelled output')
