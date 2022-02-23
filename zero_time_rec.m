clear all, clc
[x,y,Par]=eprload('Q_6hSACN_50K_2PESEEM_tau200');
tau=0.200;
x=(x/2);
x = x/10^3;
y=abs(y);
corry = bcp(y,'c',3);
figure
plot(x,y)
figure; hold on
plot(x,corry)
% corry_ = smooth(corry,5);
% plot(x, corry_)
%%
figure(1)
plot(x,y)
title('select the points')
[n,p] = ginput;
%% a
fo=fitoptions ('Method','NonlinearLeastSquares',......
               'Lower',[max(y),0.01,0],.....
               'Startpoint',[max(y) 1 0]);
ft=fittype('a*exp(-x/t)+y0','dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a', 't','y0'});
[cf]=fit(n,p,ft,fo);
        error=confint(cf,0.68);
        ae=error(2,1)-error(1,1);
        te=error(2,2)-error(1,2);
        y0e=error(2,3)-error(1,3);
        yfit=cf.a*exp(-x/cf.t)+cf.y0;
figure
plot(n,p)
xlim([0, 3])
hold on
plot(x,yfit)
hold off
%% Ist reconstruction
xn=x+tau;
ki = y(1);
k0 = yfit(1);
tauD = 0.2; %us
t = linspace(0,(tau-0.008),(tau/0.008));
k = k0-(k0-ki)*sin((pi*t)/(2*tauD));
y=y';
xn=xn';
rec1x=[t,xn];
rec1=[k,y];
figure
plot(rec1x,rec1)
%%
trasformata = 'fft';
switch trasformata
    case 'fft' 
        s1=bcp(rec1,'r',3);
        spec=ham(s1);
        s1=zf(spec,600,'r');
        spec=fft(s1);
        s1=abs(spec);
        trec1=s1(1:300); 
    case 'dct' 
        s1=bcp(rec1,'r',3);
        spec=ham(s1);
        s1=zf(spec,300,'r');
        spec=dct(s1);
        trec1=abs(spec); 
end
w = 1 / (2 * (rec1x(1, 2))); %define max frequency
v = 0:w/(numel(trec1)-1):w; %frequency axes
figure
plot(v, trec1)
title(trasformata)
%%
plot(v,trec1)
xlim([0 20])
title('select the frequency')
[v0,~]=ginput;
eta=0.86; %chosen randomly
dv=(eta*0.35)/tauD; % ?
gw = exp(-((v-v0)/dv).^4);
wind = trec1.*gw;
figure
plot(v,wind)
xlim([0 20]) 
%%
inversa = 'ifft';
switch inversa
    case 'ifft'
        back=ifft(wind,600);
%         back=abs(back);
    case 'idct'
        back=idct(wind);
        back=abs(back);
end
figure; hold on
plot(real(back))
plot(imag(back))
add = real(back(1:25));

rec2_ = [add,corry'];
f = rec2_(26)/rec2_(25);
rec2 = [add*f, corry'];

figure
rec2x=rec1x;
plot(rec2)
%%
ta = fdaxis(rec1x(2)-rec1x(1), 650);
figure; hold on
plot(ta, real(fftshift(fft(rec2,650))))
plot(ta, smooth(real(fftshift(fft(rec2,650)))))
