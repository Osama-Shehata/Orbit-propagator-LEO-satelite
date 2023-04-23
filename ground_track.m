function ground_track
%% { First, the program asks the used to input the six orbital elements of the satellite orbit. }
clear all; close all; clc
global ra dec n_curves RA Dec
%...Constants
    
deg = pi/180;
mu = 398600;
J2 = 0.00108263;
Re = 6378;
earths_angular_velocity= (2*pi + 2*pi/365.26)/(24*3600);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v=input('Please enter 1 or 2 or 3 : \n');
if (v==1)
perigee_of_orbit = 6700;
apogee_of_orbit = 10000;
TAo = 230*deg;
Wo = 270*deg;
incl = 60*deg;
e = (apogee_of_orbit - perigee_of_orbit)/(apogee_of_orbit + perigee_of_orbit);
a = (apogee_of_orbit + perigee_of_orbit)/2;
elseif (v==2)
  % TLE file name           
            % Open the TLE file and read TLE elements
            fid = fopen('TLE.txt', 'rb');
            while 1
            % read second line
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            incl = str2double(tline(10:15));                               % Orbit Inclination (degrees)
            Wo= str2double(tline(17:22));                                  % Right Ascension of Ascending Node (degrees)
            e = str2double(tline(24:28));                                  % Eccentricity
            wp = str2double(tline(30:37));                                 % Argument of Perigee (degrees)
            TAo = str2double(tline(40:44));                                % Mean Anomaly (degrees) 
            a = (str2double(tline(47:52))) ;                               % semi major axis 
            perigee_of_orbit = 6700;
            apogee_of_orbit = 10000;
            end
            fclose(fid);
elseif (v==3)
    perigee_of_orbit=input('Please enter a value for the perigee radius: \n');
    apogee_of_orbit=input('Please enter a value for the apogee radius: \n');
    TAo=input('Please enter a value for the true anamoly: \n');
    Wo = input('Please enter a value for Right ascension of the ascending node: \n');
    incl = input('Please enter a value for Right ascension of the inclination: \n'); 
    wp = input('Please enter a value for the Argument of perigeee : \n');
    n_Periods = input('Please enter a value for number of periods of the ground track: \n');
e = (apogee_of_orbit - perigee_of_orbit)/(apogee_of_orbit + perigee_of_orbit);
a = (apogee_of_orbit + perigee_of_orbit)/2;
end
%...End data declaration
%...Compute the initial time (since perigee) and
% the rates of node regression and perigee advance
%% { Third, the program does some calculation to find the time and the true anomaly of the satellite }
n_periods = 3.25;
T = 2*pi/sqrt(mu)*a^(3/2);

h = sqrt(mu*a*(1 - e^2));
Eo = 2*atan(tan(TAo/2)*sqrt((1-e)/(1+e)));
Mo = Eo - e*sin(Eo);
to = Mo*(T/2/pi);
tf = to + n_periods*T;
to = Mo*(T/2/pi);
wpo = 45*deg;



fac = -3/2*sqrt(mu)*J2*Re^2/(1-e^2)^2/a^(7/2);
Wdot = fac*cos(incl);
wpdot = fac*(5/2*sin(incl)^2 - 2);
find_ra_and_dec
form_separate_curves
plot_ground_track
print_data
return


function find_ra_and_dec               %This function calculates the right ascension and the
                                       %declination from the geocentric equatorial position vector.

% Propagates the orbit over the specified time interval, transforming
% the position vector into the earth-fixed frame and, from that,
% computing the right ascension and declination histories.

%% Fifth, for loop is used to collect the satellite data from its initial position to its final position.
times = linspace(to,tf,1000);
ra =[];
dec = [];
theta = 0;

for i = 1:length(times)
t = times(i);
%% Forth, the satellite is orbiting in an elliptic orbit, so the equation which would be used,is elliptic kepler's equation.
M = 2*pi/T*t;
E = kepler_E(e, M);
TA = 2*atan(tan(E/2)*sqrt((1+e)/(1-e)));
r = h^2/mu/(1 + e*cos(TA))*[cos(TA) sin(TA) 0]';
W = Wo + Wdot*t;
wp = wpo + wpdot*t;
%% Sixth, doing transformation from perifocal frame to geocentric frame. Then, transforming from geocentric frame to the rotating frame.
R1 = [ cos(W) sin(W) 0
-sin(W) cos(W) 0
0 0 1];
R2 = [1 0 0
0 cos(incl) sin(incl)
0 -sin(incl) cos(incl)];
R3 = [ cos(wp) sin(wp) 0
-sin(wp) cos(wp) 0
0 0 1];
QxX = (R3*R2*R1)';
R = QxX*r;
theta = earths_angular_velocity*(t - to);
Q = [ cos(theta) sin(theta) 0
-sin(theta) cos(theta) 0
0 0 1];
r_rel = Q*R;
[alpha delta] = ra_and_dec_from_r(r_rel);
ra = [ra; alpha];
dec = [dec; delta];
end
end %find_ra_and_dec

%% Breaks the ground track up into separate curves which start and terminate at right ascensions in the range [0,360 deg].
function form_separate_curves

tol = 100;
curve_no = 1;
n_curves = 1;
k = 0;
ra_prev = ra(1);
for i = 1:length(ra)
if abs(ra(i) - ra_prev) > tol
curve_no = curve_no + 1;
n_curves = n_curves + 1;
k = 0;
end
k = k + 1;
RA{curve_no}(k) = ra(i);
Dec{curve_no}(k) = dec(i);
ra_prev = ra(i);
end
end %form_separate_curves
% wwwwwwwwwwwwwwwwwwwwwwww
%% plot_ground_track
function plot_ground_track
% wwwwwwwwwwwwwwwwwwwwwwww
hold on
xlabel('East longitude (degrees)')
ylabel('Latitude (degrees)')
axis equal
grid on
for i = 1:n_curves
plot(RA{i}, Dec{i})
end

axis ([0 360 -90 90])
text( ra(1), dec(1), 'o Start','color','blue')
text(ra(end), dec(end), 'o Finish','color','red')
line([min(ra) max(ra)],[0 0], 'Color','k') %the equator
 I=imread ('earth.jpg');
 m=image(xlim,-ylim,I);
  uistack(m,'bottom')

end
%% print the data
    function print_data
coe = [h e Wo incl wpo TAo];
[ro, vo] = sv_from_coe(coe, mu);
disp('\n Angular momentum = km^2/s');
disp(h);
disp('\n Eccentricity = %g');
disp(e);
disp('\n Semimajor axis = %g km' );
disp(a);
disp('\n Perigee radius = %g km');
disp(perigee_of_orbit);
disp('\n Apogee radius = %g km');
disp(apogee_of_orbit);
disp('\n Period = %g hours');
disp(T/3600);
disp('\n Inclination = %g deg');
disp(incl/deg);
disp('\n Initial true anomaly = %g deg');
disp(TAo/deg);
disp('\n Time since perigee = %g hours');
disp(to/3600);
disp('\n Initial apogee_of_orbit = %g deg');
disp( Wo/deg);
disp('\n RA_dot = %g deg/period');
disp( Wdot/deg*T);
disp('\n Initial wp = %g deg');
disp(wpo/deg);
disp('\n wp_dot = %g deg/period');
disp(wpdot/deg*T);


fprintf('\n r0 = [%12g, %12g, %12g] (km)', ro(1), ro(2), ro(3))
fprintf('\n magnitude = %g km\n', norm(ro))
fprintf('\n v0 = [%12g, %12g, %12g] (km)', vo(1), vo(2), vo(3))
fprintf('\n magnitude = %g km\n', norm(vo))

    end
end %ground_track