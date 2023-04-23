function [R3_w,R1_i,R3_W,Q_pX,r,v]= sv_from_coe(coe,mu)

h = coe(1);
e = coe(2);
RA = coe(3);
incl = coe(4);
w = coe(5);
TA = coe(6);
rp = (h^2/mu) * (1/(1 + e*cos(TA))) * (cos(TA)*[1;0;0] + sin(TA)*[0;1;0]);
vp = (mu/h) * (-sin(TA)*[1;0;0] + (e + cos(TA))*[0;1;0]);

R3_W = [ cos(RA) sin(RA) 0
-sin(RA) cos(RA) 0
0 0 1];
%...Equation 4.32:
R1_i = [1 0 0
0 cos(incl) sin(incl)
0 -sin(incl) cos(incl)];
%...Equation 4.34:
R3_w = [ cos(w) sin(w) 0
-sin(w) cos(w) 0
0 0 1];
Q_pX = (R3_w*R1_i*R3_W)';
r = Q_pX*rp;
v = Q_pX*vp;
end