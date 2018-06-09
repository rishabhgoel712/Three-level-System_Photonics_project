clear;
clc;
close all;
tic
n=1; %current iteration
n1=2500000;%total number of iterations
p=10^7; %pumping rate
LT=10e-9;  %initialised here, follows exponential distribution 
T=exprnd(LT,1,n1);
gnr1=145; %rate of going from excited to trap state
gnr2=2.84;  %rate of going from trap to ground state
bin_size=10^-2;  %size of each bin

t(1)=0; % stores time at which intensity is measured
intensity(1)=0; %intensity at a particular time
i=2; %index for time array
j=2; %index for intensity array
k=1; %index for choose array
q=1; %index for gsp1 array(array of gsp) 
c=1; %index for iterating over time while binning
count=0; %number of time it goes in trap state
while n<=n1
    gsp1=1./T(n); %random number generated for gsp from exponential distribution
   %{
       while gsp<10^7  % put limit to value of gsp
        gsp=exprnd(10^8);
       end
        %}
    S3=sqrt(p.^2-2*p*gnr1-2*p*gnr2+2*p*gsp1+gnr1.^2-2*gnr1*gnr2+2*gnr1*gsp1+gnr2.^2-2*gnr2*gsp1+gsp1.^2)/2;
    S2=p/2 +S3 +gsp1/2 +gnr1/2 +gnr2/2;
    S1=p/2 -S3 +gsp1/2 +gnr1/2 +gnr2/2;
    A1=(p*gnr1)/(p*gnr1+p*gnr2+gnr1*gnr2+gnr2*gsp1);
    A2=(S1/(2*S3))*A1;
    A3=-(S2/(2*S3))*A1;
    gsp2(q)=gsp1; %store gsp in gsp1
    choose(k)=rand();  %random number from 0 to 1 for choosing which state it goes
    if choose(k)<(gnr1/(gsp1+gnr1)) % case when carrier goes into trap state
        intensity(j)=0; 
        j=j+1;
        t(i)=t(i-1)+ 1/gnr1 +1/gnr2; %rise in time
        i=i+1;
        count=count+1;
    else  %case when it goes to ground state directly
        t(i)=t(i-1)+1/gsp1; %rise in time
        intensity(j)=(gsp1^2/(gsp1+gnr1))*((gnr2*A1/gnr1)+((gnr2-S2)/gnr1)*A2*exp(-S2/gsp1) + ((gnr2-S1)/gnr1)*A3*exp(-S1/gsp1));  % function for intensity
        i=i+1; %increment in time index
        j=j+1;  %increment in intensity index
    end
        
    k=k+1; 
    n=n+1;
    q=q+1;
end
plot(t, intensity,'o')
figure()
n=fix(max(t)/bin_size); %number of bins formed
i1=zeros(n,1); %array initialized with zero to store intensity after binning
for m=1:n %iterating over each bin
 cprev=c;
   while t(c)<m*bin_size
      i1(m)=i1(m)+intensity(c);
      c=c+1;
   end
   diff=c-cprev;
 i1(m)=i1(m)/diff; %average of intensity is stored after binning
end
i1(isnan(i1)) =0;%replace all NaN to 0 in i1 array
t1 = 0:bin_size:(n-1)*bin_size;

c = [[0,0,0];[0.7,0.7,0.7];[1,0,0];[0,1,0];[0,0,1]];
line_style = {'-','--','-.','-'};
plot(t1,i1,'LineWidth',2,'Color', [c(3,:)],'LineStyle',line_style{1});
xlim([0, max(t)]);

h_legend = legend('3 level','Location','NorthEast');
title(['n=',num2str(n1),', gnr1=',num2str(gnr1),', gnr2=',num2str(gnr2), ', Noff =',num2str(count)],'fontsize',16)
set(h_legend,'fontsize',16, 'box', 'off');
xlabel('time (s)','fontsize',24);
ylabel('Intensity (a.u.)','fontsize',24);
set(gca, 'Fontsize',24);

figname_png = ['Intensity and time for 3 level system9.png'];
figname = ['Intensity and time for 3 level system9'];
width = 15;
height = 4;
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPosition', [0 0 width height]);
set(gca,'position',[0.1 0.19 .85 .7]);% specify these as the fraction of the total.. between 0 and 1

print('-dpng','-r125',figname_png);
toc