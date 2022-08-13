function [B1,B2,B3]=computingB(phi1_comp, phi2_comp,phi3_comp,phi_m1,phi_m2,phi_m3,JD_presb,CV_prescb1,CV_prescb2,CV_prescb3,phi_o_jd,phi_o1_y,phi_o1_x, phi_o1_z,phi_o2_y,phi_o2_x, phi_o2_z,phi_o3_y,phi_o3_x, phi_o3_z)
h=1;
%% find P and Q on the current solution phi_comp
[phi_jd, phi_cv1, phi_cv2, phi_cv3]=compute_JD_and_Curl3D(phi1_comp, phi2_comp,phi3_comp,h);
[phi_m1_y,phi_m1_x,phi_m1_z]=gradient(phi_m1,h);
[phi_m2_y,phi_m2_x,phi_m2_z]=gradient(phi_m2,h);
[phi_m3_y,phi_m3_x,phi_m3_z]=gradient(phi_m3,h);

PdetG=(phi_jd - JD_presb).*phi_o_jd;
Q1=phi_cv1 - CV_prescb1;
Q2=phi_cv2 - CV_prescb2;
Q3=phi_cv3 - CV_prescb3;

%% form A
P_11=PdetG.*(phi_m2_y.*phi_m3_z-phi_m2_z.*phi_m3_y);
Q_11=Q2.*phi_o1_z-Q3.*phi_o1_y;
A11=P_11+Q_11;

P_12=PdetG.*(phi_m3_x.*phi_m2_z-phi_m2_x.*phi_m3_z);
Q_12=Q2.*phi_o2_z-Q3.*phi_o2_y;
A12=P_12+Q_12;

P_13=PdetG.*(phi_m2_x.*phi_m3_y-phi_m2_y.*phi_m3_x);
Q_13=Q2.*phi_o3_z-Q3.*phi_o3_y;
A13=P_13+Q_13;

P_21=PdetG.*(phi_m3_y.*phi_m1_z-phi_m1_y.*phi_m3_z);
Q_21=-Q1.*phi_o1_z+Q3.*phi_o1_x;
A21=P_21+Q_21;

P_22=PdetG.*(phi_m1_x.*phi_m3_z-phi_m1_z.*phi_m3_x);
Q_22=-Q1.*phi_o2_z+Q3.*phi_o2_x;
A22=P_22+Q_22;

P_23=PdetG.*(phi_m3_x.*phi_m1_y-phi_m1_x.*phi_m3_y);
Q_23=-Q1.*phi_o3_z+Q3.*phi_o3_x;
A23=P_23+Q_23;

P_31=PdetG.*(phi_m3_y.*phi_m2_z-phi_m2_y.*phi_m1_z);
Q_31=Q1.*phi_o1_y-Q2.*phi_o1_x;
A31=P_31+Q_31;

P_32=PdetG.*(phi_m2_x.*phi_m1_z-phi_m1_x.*phi_m2_z);
Q_32=Q1.*phi_o2_y-Q2.*phi_o2_x;
A32=P_32+Q_32;

P_33=PdetG.*(phi_m1_x.*phi_m2_y-phi_m2_x.*phi_m1_y);
Q_33=Q1.*phi_o3_y-Q2.*phi_o3_x;
A33=P_33+Q_33;
%% form B
[~,A11_x,~]=gradient(A11,h);
[A12_y,~,~]=gradient(A12,h);
[~,~,A13_z]=gradient(A13,h);
lap_B1=A11_x+A12_y+A13_z;

[~,A21_x,~]=gradient(A21,h);
[A22_y,~,~]=gradient(A22,h);
[~,~,A23_z]=gradient(A23,h);
lap_B2=A21_x+A22_y+A23_z;

[~,A31_x,~]=gradient(A31,h);
[A32_y,~,~]=gradient(A32,h);
[~,~,A33_z]=gradient(A33,h);
lap_B3=A31_x+A32_y+A33_z;

B1=pois3fft(-lap_B1);
B2=pois3fft(-lap_B2);
B3=pois3fft(-lap_B3);
end
