function [T1_new,T2_new,T3_new,ssd_initial]=PJDCV_fixedBC_3D(f0,g01,g02,g03,N,tstep)
% construct a uniform square grid on [0,1]*[0,1]*[0,1]
% Physican Domain [0,1]*[0,1]*[0,1]

% N=33;% actual grid
Npts=N-2; % interior grid
% number of points each direction
% alpha=1;
h=1;

%initial value for T=x
[T1,T2,T3]=ndgrid(1:N,1:N,1:N);

% read jacobian file
f0=f0(2:Npts+1,2:Npts+1,2:Npts+1);
g01=g01(2:Npts+1,2:Npts+1,2:Npts+1);
g02=g02(2:Npts+1,2:Npts+1,2:Npts+1);
g03=g03(2:Npts+1,2:Npts+1,2:Npts+1);

%initial value for control function
f1=zeros(Npts,Npts,Npts);
f2=zeros(Npts,Npts,Npts);
f3=zeros(Npts,Npts,Npts);
u1=zeros(Npts,Npts,Npts);
u2=zeros(Npts,Npts,Npts);
u3=zeros(Npts,Npts,Npts);
%initial value for ssd
[T1y,T1x,T1z]=gradient(T1,h);
[T2y,T2x,T2z]=gradient(T2,h);
[T3y,T3x,T3z]=gradient(T3,h);


JT=T1x.*(T2y.*T3z-T3y.*T2z)-T2x.*(T1y.*T3z-T3y.*T1z)+T3x.*(T1y.*T2z-T2y.*T1z);
curlT1=T3y-T2z;
curlT2=T1z-T3x;
curlT3=T2x-T1y;
JT_interior=JT(2:Npts+1,2:Npts+1,2:Npts+1);
T1y=T1y(2:Npts+1,2:Npts+1,2:Npts+1);
T1x=T1x(2:Npts+1,2:Npts+1,2:Npts+1);
T1z=T1z(2:Npts+1,2:Npts+1,2:Npts+1);
T2y=T2y(2:Npts+1,2:Npts+1,2:Npts+1);
T2x=T2x(2:Npts+1,2:Npts+1,2:Npts+1);
T2z=T2z(2:Npts+1,2:Npts+1,2:Npts+1);
T3y=T3y(2:Npts+1,2:Npts+1,2:Npts+1);
T3x=T3x(2:Npts+1,2:Npts+1,2:Npts+1);
T3z=T3z(2:Npts+1,2:Npts+1,2:Npts+1);

%P=JT-f0
%only compute the interior part
P=JT_interior-f0;
Q1=curlT1(2:Npts+1,2:Npts+1,2:Npts+1)-g01;
Q2=curlT2(2:Npts+1,2:Npts+1,2:Npts+1)-g02;
Q3=curlT3(2:Npts+1,2:Npts+1,2:Npts+1)-g03;
ssd_old=sum(sum(sum(P.^2+Q1.^2+Q2.^2+Q3.^2)));
ssd_initial=ssd_old;
ratio=1;
%iteration parameters
imax=1.0e4;
% tstep=1;
tstepratio=1e-7;
better=1;
iter=0;
ratiotol=1e-5;
% ssd_tol=1e-5;
while iter<imax && tstep>tstepratio && ratio>ratiotol
    iter=iter+1;
    if better
   
        %ssd decreases
        b1=P.*(T2y.*T3z-T3y.*T2z);
        b2=(P.*(T3x.*T2z-T2x.*T3z)-Q3);
        b3=(P.*(T2x.*T3y-T2y.*T3x)+Q2);
        [~,b1x,~]=gradient(b1,h);
        [b2y,~,~]=gradient(b2,h);
        [~,~,b3z]=gradient(b3,h);
        a1=b1x+b2y+b3z;
        
        g1=pois3fft(-a1);

        b1=(P.*(T3y.*T1z-T1y.*T3z)+Q3);
        b2=P.*(T1x.*T3z-T3x.*T1z);
        b3=(P.*(T3x.*T1y-T3y.*T1x)-Q1);
        [~,b1x,~]=gradient(b1,h);
        [b2y,~,~]=gradient(b2,h);
        [~,~,b3z]=gradient(b3,h);
        a2=b1x+b2y+b3z;
        
        g2=pois3fft(-a2);
        
        b1=(P.*(T1y.*T2z-T2y.*T1z)-Q2);
        b2=(P.*(T2x.*T1z-T1x.*T2z)+Q1);
        b3=P.*(T1x.*T2y-T1y.*T2x);
        [~,b1x,~]=gradient(b1,h);
        [b2y,~,~]=gradient(b2,h);
        [~,~,b3z]=gradient(b3,h);
        a3=b1x+b2y+b3z;
        
        g3=pois3fft(-a3);
        
    end
    %update f1n,f2n to calculate u1,u2
    f1n=f1-tstep*g1;
    f2n=f2-tstep*g2;
    f3n=f3-tstep*g3;
    %update u1,u2
    u1=pois3fft(f1n);
    u2=pois3fft(f2n);
    u3=pois3fft(f3n);
    %update T1_new,T2_new to compute ssd
    T1_new=T1;
    T2_new=T2;
    T3_new=T3;
    T1_new(2:Npts+1,2:Npts+1,2:Npts+1)=T1_new(2:Npts+1,2:Npts+1,2:Npts+1)+u1;
    T2_new(2:Npts+1,2:Npts+1,2:Npts+1)=T2_new(2:Npts+1,2:Npts+1,2:Npts+1)+u2;
    T3_new(2:Npts+1,2:Npts+1,2:Npts+1)=T3_new(2:Npts+1,2:Npts+1,2:Npts+1)+u3;
    %update JT
    [T1y,T1x,T1z]=gradient(T1_new,h);
    [T2y,T2x,T2z]=gradient(T2_new,h);
    [T3y,T3x,T3z]=gradient(T3_new,h);
    
    JT=T1x.*(T2y.*T3z-T3y.*T2z)-T2x.*(T1y.*T3z-T3y.*T1z)+T3x.*(T1y.*T2z-T2y.*T1z);
    curlT1=T3y-T2z;
    curlT2=T1z-T3x;
    curlT3=T2x-T1y;
    JT_interior=JT(2:Npts+1,2:Npts+1,2:Npts+1);
    T1y=T1y(2:Npts+1,2:Npts+1,2:Npts+1);
    T1x=T1x(2:Npts+1,2:Npts+1,2:Npts+1);
    T1z=T1z(2:Npts+1,2:Npts+1,2:Npts+1);
    T2y=T2y(2:Npts+1,2:Npts+1,2:Npts+1);
    T2x=T2x(2:Npts+1,2:Npts+1,2:Npts+1);
    T2z=T2z(2:Npts+1,2:Npts+1,2:Npts+1);
    T3y=T3y(2:Npts+1,2:Npts+1,2:Npts+1);
    T3x=T3x(2:Npts+1,2:Npts+1,2:Npts+1);
    T3z=T3z(2:Npts+1,2:Npts+1,2:Npts+1);
    P=JT_interior-f0;
    Q1=curlT1(2:Npts+1,2:Npts+1,2:Npts+1)-g01;
    Q2=curlT2(2:Npts+1,2:Npts+1,2:Npts+1)-g02;
    Q3=curlT3(2:Npts+1,2:Npts+1,2:Npts+1)-g03;
    ssd_new=sum(sum(sum(P.^2+Q1.^2+Q2.^2+Q3.^2)));
    % if ssd decrease, then update f,T, increase tstep
    % if not, decrese tstep
    if (ssd_new<ssd_old)
        tstep=tstep*1.1;
        f1=f1n;
        f2=f2n;
        f3=f3n;
        ssd_old=ssd_new;
        ratio=ssd_new/ssd_initial;
        better=1;
    else
        better=0;
        tstep=tstep*0.9;
    end 
%     display([' tstep:',num2str(tstep),' ssd old:',num2str(ssd_initial),' ssd:',num2str(ssd_old),' ratio:',num2str(ratio),' iter:',num2str(iter)]);
end
display([' tstep:',num2str(tstep),' ssd old:',num2str(ssd_initial),' ssd:',num2str(ssd_old),' ratio:',num2str(ratio),' iter:',num2str(iter)]);
T1_new=T1;
T2_new=T2;
T3_new=T3;
T1_new(2:Npts+1,2:Npts+1,2:Npts+1)=T1_new(2:Npts+1,2:Npts+1,2:Npts+1)+u1;
T2_new(2:Npts+1,2:Npts+1,2:Npts+1)=T2_new(2:Npts+1,2:Npts+1,2:Npts+1)+u2;
T3_new(2:Npts+1,2:Npts+1,2:Npts+1)=T3_new(2:Npts+1,2:Npts+1,2:Npts+1)+u3;
