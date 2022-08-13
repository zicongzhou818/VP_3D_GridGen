function [phi1_comp,phi2_comp,phi3_comp,U1,U2,U3,phi_m1, phi_m2,phi_m3,U1_m,U2_m,U3_m,ratio]=PJCD_on_given_mesh3D2(JD_presb,CV_presb1,CV_presb2,CV_presb3,N,given_x,given_y,given_z)
%% initialization
better = 1;
ei = 0;
ratio = 1;
ti = 0;
tstep = 1e-4;
ts_up = 1.1;
ts_dn = 0.9;
Npts = N-2;
imax = 1.0e4;
tstepratio = 1e-7;
ratiotol = 1e-5;
%% preset control functions corresponsive to the identity map
f1 = zeros(N-2,N-2,N-2);
f2 = zeros(N-2,N-2,N-2);
f3 = zeros(N-2,N-2,N-2);
%% preset new deformation based on control functions that are corresponsive to the identity map
[X,Y,Z] = ndgrid(1:N,1:N,1:N);
phi_m1 = X;
phi_m2 = Y;
phi_m3 = Z;
%% reading prescriptions
fo=JD_presb(2:Npts+1,2:Npts+1,2:Npts+1);
go1=CV_presb1(2:Npts+1,2:Npts+1,2:Npts+1);
go2=CV_presb2(2:Npts+1,2:Npts+1,2:Npts+1);
go3=CV_presb3(2:Npts+1,2:Npts+1,2:Npts+1);
%% reading given mesh 
h=1;
phi_o1=given_x;
phi_o2=given_y;
phi_o3=given_z;
[phi_o1_y,phi_o1_x,phi_o1_z]=gradient(phi_o1,h);
[phi_o2_y,phi_o2_x,phi_o2_z]=gradient(phi_o2,h);
[phi_o3_y,phi_o3_x,phi_o3_z]=gradient(phi_o3,h);
%% (ps: the final phi initially is the same as given phi here because the current phi is id)
phi1_comp=given_x;
phi2_comp=given_y;
phi3_comp=given_z;
%% compute ssd
[phi_o_jd, phi_o_cv1,phi_o_cv2,phi_o_cv3]=compute_JD_and_Curl3D(phi_o1,phi_o2,phi_o3,h);
Pn=(phi_o_jd(2:Npts+1,2:Npts+1,2:Npts+1)-fo);
Qn1=(phi_o_cv1(2:Npts+1,2:Npts+1,2:Npts+1)-go1);
Qn2=(phi_o_cv2(2:Npts+1,2:Npts+1,2:Npts+1)-go2);
Qn3=(phi_o_cv3(2:Npts+1,2:Npts+1,2:Npts+1)-go3);
ssd_initial=sum(sum(sum(Pn.^2+Qn1.^2+Qn2.^2+Qn3.^2)));
ssd_old=ssd_initial;

%% main loop
while ei<imax && tstep>tstepratio && ratio>ratiotol
    ti=ti+1;
    if better % form gradients wrt F when ssd decreases
        ei=ei+1;
        %% form B
        [B1,B2,B3]=computingB(phi1_comp, phi2_comp,phi3_comp,phi_m1,phi_m2,phi_m3,JD_presb,CV_presb1,CV_presb2,CV_presb3,phi_o_jd,phi_o1_y,phi_o1_x, phi_o1_z,phi_o2_y,phi_o2_x, phi_o2_z,phi_o3_y,phi_o3_x, phi_o3_z);
    end
    %update f1n,f2n to calculate u1,u2
    f1n=f1-tstep*B1(2:Npts+1,2:Npts+1,2:Npts+1);
    f2n=f2-tstep*B2(2:Npts+1,2:Npts+1,2:Npts+1);
    f3n=f3-tstep*B3(2:Npts+1,2:Npts+1,2:Npts+1);
    %update u1,u2
    u1=pois3fft(f1n);
    u2=pois3fft(f2n);
    u3=pois3fft(f3n);
    U1_m = matrixpad3D(u1,0);
    U2_m = matrixpad3D(u2,0);
    U3_m = matrixpad3D(u3,0);
    %update T_new to compute ssd
%     k1=25;k2=17;h2=2;
%     figure(11)
%     gridplot3D_flexibleR(phi1_comp,phi2_comp,phi3_comp,h2,h2,h2),xlabel('x'),ylabel('y'),zlabel('z');
%     figure(12)
%     gridplot3D_framed_axis(phi1_comp,phi2_comp,phi3_comp,1,k1,1),
%     gridplot3D_framed_axis(phi1_comp,phi2_comp,phi3_comp,1,k2,3), grid on
%     xlabel('x'),ylabel('y'),zlabel('z');
    phi_m1=X+U1_m;
    phi_m2=Y+U2_m;
    phi_m3=Z+U3_m;
    
%     figure(21)
%     gridplot3D_flexibleR(phi_m1,phi_m2,phi_m3,h2,h2,h2),xlabel('x'),ylabel('y'),zlabel('z');
%     figure(22)
%     gridplot3D_framed_axisR(phi_m1,phi_m2,phi_m3,1,k1,1),
%     gridplot3D_framed_axisR(phi_m1,phi_m2,phi_m3,1,k2,3), grid on
%     xlabel('x'),ylabel('y'),zlabel('z');
%     %%
    phi1_comp = interpn(X, Y, Z, phi_o1, phi_m1, phi_m2, phi_m3, 'makima'); 
    phi2_comp = interpn(X, Y, Z, phi_o2, phi_m1, phi_m2, phi_m3, 'makima'); 
    phi3_comp = interpn(X, Y, Z, phi_o3, phi_m1, phi_m2, phi_m3, 'makima'); 
    %%
%     figure(31)
%     gridplot3D_flexibleR(phi1_comp,phi2_comp,phi3_comp,h2,h2,h2),xlabel('x'),ylabel('y'),zlabel('z');
%     figure(12)
%     gridplot3D_framed_axisR(phi1_comp,phi2_comp,phi3_comp,1,k1,1),
%     gridplot3D_framed_axisR(phi1_comp,phi2_comp,phi3_comp,1,k2,3), grid on
%     xlabel('x'),ylabel('y'),zlabel('z');
    %update JT,curlT
    [phi_jd, phi_cv1,phi_cv2,phi_cv3]=compute_JD_and_Curl3D(phi1_comp,phi2_comp,phi3_comp,h);
    Pn=(phi_jd(2:Npts+1,2:Npts+1,2:Npts+1)-fo);
    Qn1=(phi_cv1(2:Npts+1,2:Npts+1,2:Npts+1)-go1);
    Qn2=(phi_cv2(2:Npts+1,2:Npts+1,2:Npts+1)-go2);
    Qn3=(phi_cv3(2:Npts+1,2:Npts+1,2:Npts+1)-go3);
%     Pn=(phi_jd-JD_prescribed);
%     Qn=(phi_cv-CV_prescribed);
    ssd_new=sum(sum(sum(Pn.^2+Qn1.^2+Qn2.^2+Qn3.^2)));
    ratio=ssd_new/ssd_initial;
    % if ssd decrease, then update f,T, increase tstep
    % if not, decrese tstep
%     display([' ts: ',num2str(tstep),' r: ',num2str(ratio), ' ei: ',num2str(ei), ' ti: ',num2str(ti)]);
    if (ssd_new<ssd_old)
        tstep=tstep*ts_up;
        f1=f1n;
        f2=f2n;
        f3=f3n;
        ssd_old=ssd_new;
        better=1;
    else
        better=0;
        tstep=tstep*ts_dn;
    end 
end
display([' ts: ',num2str(tstep),' r: ',num2str(ratio), ' ei: ',num2str(ei), ' ti: ',num2str(ti)]);
U1 = phi1_comp-X;
U2 = phi2_comp-Y;
U3 = phi3_comp-Z;


end
