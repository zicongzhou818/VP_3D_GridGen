close all
clear
clc


N=51;k1=37;k2=25;
Npts=N-2; 

h=1;
[x1,x2,x3]=ndgrid(1:N,1:N,1:N);
phi1=x1;
phi2=x2;
phi3=x3;

Phi1=x1;
Phi2=x2;
Phi3=x3;
for i=1:N
    for j=1:N
        for k=1:N
            x=phi1(i,j,k);
            y=phi2(i,j,k);
            z=phi3(i,j,k);
            [x_new, y_new, z_new]=cut_off_small_3D(x,y,z,N,-pi/8);
            phi1(i,j,k)=x_new;
            phi2(i,j,k)=y_new;
            phi3(i,j,k)=z_new;
        end
    end
end
for i=1:N
    for j=1:N
        for k=1:N
            x=phi1(i,j,k);
            y=phi2(i,j,k);
            z=phi3(i,j,k);
            [x_new, y_new, z_new]=cut_off_3D(x,y,z,N,pi/3);
            [x_new, y_new, z_new]=shift_var(x_new, y_new, z_new,N,0.5,1);
            [x_new, y_new, z_new]=shift_var(x_new, y_new, z_new,N,0.5,2);
            [x_new, y_new, z_new]=shift_var(x_new, y_new, z_new,N,0.5,3);
            Phi1(i,j,k)=x_new;
            Phi2(i,j,k)=y_new;
            Phi3(i,j,k)=z_new;
        end
    end
end

pPhi1=Phi1;
pPhi2=Phi2;
pPhi3=Phi3;
Up1 = pPhi1-x1;
Up2 = pPhi2-x2;
Up3 = pPhi3-x3;
[JD_Phi, CV_Phi1, CV_Phi2, CV_Phi3]=compute_JD_and_Curl3D(pPhi1,pPhi2,pPhi3,h);
min(min(min(JD_Phi)))
h2=1;
% figure(1)
% gridplot3D_flexible(pPhi1,pPhi2,pPhi3,h2,h2,h2),xlabel('x'),ylabel('y'),zlabel('z');
figure(2)
gridplot3D_framed_axis(pPhi1,pPhi2,pPhi3,1,k1,1),
gridplot3D_framed_axis(pPhi1,pPhi2,pPhi3,1,k2,3), grid on
axis([0, N+1, 0, N+1, 0, N+1]);
xlabel('x'),ylabel('y'),zlabel('z');
figure(3)
quiver3(x1(:, :, k2), x2(:, :, k2), x3(:, :, k2), Up1(:, :, k2),Up2(:, :, k2),Up3(:, :, k2), 0,'g'), hold on
quiver3(x1(k1, :, :), x2(k1, :, :), x3(k1, :, :), Up1(k1, :, :),Up2(k1, :, :),Up3(k1, :, :), 0,'g'), grid on;
axis([0, N+1, 0, N+1, 0, N+1]);
xlabel('x'),ylabel('y'),zlabel('z');


%% iterative method old
tstep=1e-4;
tic
[T1,T2,T3]=PJDCV_fixedBC_3D(JD_Phi,CV_Phi1, CV_Phi2, CV_Phi3,N,tstep);
toc
[JD_T, CV_T1, CV_T2, CV_T3]=compute_JD_and_Curl3D(T1,T2,T3,h);
U1 = T1-x1;
U2 = T2-x2;
U3 = T3-x3;

figure(4)
gridplot3D_framed_axisR(T1,T2,T3,1,k1,1),
gridplot3D_framed_axisR(T1,T2,T3,1,k2,3), grid on
axis([0, N+1, 0, N+1, 0, N+1]);
xlabel('x'),ylabel('y'),zlabel('z');

figure(5)
quiver3(x1(:, :, k2), x2(:, :, k2), x3(:, :, k2), U1(:, :, k2),U2(:, :, k2),U3(:, :, k2), 0,'b'), hold on
quiver3(x1(k1, :, :), x2(k1, :, :), x3(k1, :, :), U1(k1, :, :),U2(k1, :, :),U3(k1, :, :), 0,'b'), grid on;
axis([0, N+1, 0, N+1, 0, N+1]);
xlabel('x'),ylabel('y'),zlabel('z');

figure(6)
gridplot3D_framed_axis(pPhi1,pPhi2,pPhi3,1,k1,1)
gridplot3D_framed_axis(pPhi1,pPhi2,pPhi3,1,k2,3)
gridplot3D_framed_axisR(T1,T2,T3,1,k1,1),
gridplot3D_framed_axisR(T1,T2,T3,1,k2,3), grid on
axis([0, N+1, 0, N+1, 0, N+1]);
xlabel('x'),ylabel('y'),zlabel('z');

figure(7)
quiver3(x1(:, :, k2), x2(:, :, k2), x3(:, :, k2), Up1(:, :, k2),Up2(:, :, k2),Up3(:, :, k2), 0,'g'), hold on
quiver3(x1(k1, :, :), x2(k1, :, :), x3(k1, :, :), Up1(k1, :, :),Up2(k1, :, :),Up3(k1, :, :), 0,'g'), hold on;
quiver3(x1(:, :, k2), x2(:, :, k2), x3(:, :, k2), U1(:, :, k2),U2(:, :, k2),U3(:, :, k2), 0,'b'), hold on
quiver3(x1(k1, :, :), x2(k1, :, :), x3(k1, :, :), U1(k1, :, :),U2(k1, :, :),U3(k1, :, :), 0,'b'), grid on;
axis([0, N+1, 0, N+1, 0, N+1]);
xlabel('x'),ylabel('y'),zlabel('z');
view(0.0, 0.0);
diff_x1T=pPhi1-T1;
diff_x2T=pPhi2-T2;
diff_x3T=pPhi3-T3;

diff_JDT=abs(JD_T-JD_Phi);
diff_CV1T=CV_T1-CV_Phi1;
diff_CV2T=CV_T2-CV_Phi2;
diff_CV3T=CV_T3-CV_Phi3;

max_JDT=max(max(max(diff_JDT)))
mag_CVT=(diff_CV1T.^2+diff_CV2T.^2+diff_CV3T.^2).^(.5);
max_CVT=max(max(max(mag_CVT)))

mag_distT=(diff_x1T.^2+diff_x2T.^2+diff_x3T.^2).^(.5);
max_mag_distT=max(max(max(mag_distT)))
max_rel_mag2T=max(max(max(mag_distT./((N-1)^3))))

%  tstep:9.7956e-08 ssd old:65544.7911 ssd:34.2581 ratio:0.00052267 iter:3016
% Elapsed time is 1765.896838 seconds.
% 
% max_JDT =
% 
%     0.7279
% 
% 
% max_CVT =
% 
%     0.0866
% 
% 
% max_mag_distT =
% 
%     1.9454
% 
% 
% max_rel_mag2T =
% 
%    1.5563e-05

%% iterative method new
tstep=1e-3;
tic
[T1_new,T2_new,T3_new,U1_new,U2_new,U3_new,phi_m1, phi_m2,phi_m3,U1_m,U2_m,U3_m,ratio111]=PJCD_on_given_mesh3D2(JD_Phi,CV_Phi1, CV_Phi2, CV_Phi3,N,x1,x2,x3);
toc
[JD_T, CV_T1, CV_T2, CV_T3]=compute_JD_and_Curl3D(T1_new,T2_new,T3_new,h);
diff_x1=pPhi1-T1_new;
diff_x2=pPhi2-T2_new;
diff_x3=pPhi3-T3_new;

diff_JD=abs(JD_T-JD_Phi);
diff_CV1=CV_T1-CV_Phi1;
diff_CV2=CV_T2-CV_Phi2;
diff_CV3=CV_T3-CV_Phi3;


figure(8)
gridplot3D_framed_axisR(T1_new,T2_new,T3_new,1,k1,1),
gridplot3D_framed_axisR(T1_new,T2_new,T3_new,1,k2,3), grid on
axis([0, N+1, 0, N+1, 0, N+1]);
xlabel('x'),ylabel('y'),zlabel('z');

figure(9)
quiver3(x1(:, :, k2), x2(:, :, k2), x3(:, :, k2), U1_new(:, :, k2),U2_new(:, :, k2),U3_new(:, :, k2), 0,'b'), hold on
quiver3(x1(k1, :, :), x2(k1, :, :), x3(k1, :, :), U1_new(k1, :, :),U2_new(k1, :, :),U3_new(k1, :, :), 0,'b'), grid on;
axis([0, N+1, 0, N+1, 0, N+1]);
xlabel('x'),ylabel('y'),zlabel('z');

figure(10)
gridplot3D_framed_axis(pPhi1,pPhi2,pPhi3,1,k1,1)
gridplot3D_framed_axis(pPhi1,pPhi2,pPhi3,1,k2,3)
gridplot3D_framed_axisR(T1_new,T2_new,T3_new,1,k1,1),
gridplot3D_framed_axisR(T1_new,T2_new,T3_new,1,k2,3), grid on
axis([0, N+1, 0, N+1, 0, N+1]);
xlabel('x'),ylabel('y'),zlabel('z');

figure(11)
quiver3(x1(:, :, k2), x2(:, :, k2), x3(:, :, k2), Up1(:, :, k2),Up2(:, :, k2),Up3(:, :, k2), 0,'g'), hold on
quiver3(x1(k1, :, :), x2(k1, :, :), x3(k1, :, :), Up1(k1, :, :),Up2(k1, :, :),Up3(k1, :, :), 0,'g'), hold on;
quiver3(x1(:, :, k2), x2(:, :, k2), x3(:, :, k2), U1_new(:, :, k2),U2_new(:, :, k2),U3_new(:, :, k2), 0,'b'), hold on
quiver3(x1(k1, :, :), x2(k1, :, :), x3(k1, :, :), U1_new(k1, :, :),U2_new(k1, :, :),U3_new(k1, :, :), 0,'b'), grid on;
axis([0, N+1, 0, N+1, 0, N+1]);
xlabel('x'),ylabel('y'),zlabel('z');

max_JD=max(max(max(diff_JD)))

mag_CV=(diff_CV1.^2+diff_CV2.^2+diff_CV3.^2).^(.5);
max_CV=max(max(max(mag_CV)))

mag_dist=(diff_x1.^2+diff_x2.^2+diff_x3.^2).^(.5);
max_mag_dist=max(max(max(mag_dist)))
max_rel_mag2=max(max(max(mag_dist./((N-1)^3))))

% ts: 0.021685 r: 9.9984e-06 ei: 576 ti: 1046
% Elapsed time is 788.216963 seconds.
% 
% max_JD =
% 
%     0.0104
% 
% 
% max_CV =
% 
%     0.0136
% 
% 
% max_mag_dist =
% 
%     0.0403
% 
% 
% max_rel_mag2 =
% 
%    3.2218e-07


