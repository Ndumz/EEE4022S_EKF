%This is an implementation of the kalman filter for postion and velocity
%navigation

%--------------- variables and Matrices ---------------%
e= 0.08181919084;
a = 6378137;
%Sample time
%T= 0.06968;
T = 1;%0.069;
data = readtable('roadtrip.csv');
truedata = readtable('roadtrip_true.csv');
N = height(data);

%intial states
vx=0;
vy=0;
vz=0;
 px=5.023801748007990e+06;
 py=1.676872550947940e+06;
 pz=-3.5422e+06;
%px  = 5.023370826601166e+06;
%py = 1.677145974396276e+06;
%pz = -3.542608037303097e+06;


% %a prior estimate
% x_hat = [0;0;0;5.0239e+06;1.6769e+06;-3.5422e+06];
 x_hat = [0;0;0;px;py;pz];
 x_rms = [0;0;0;px;py;pz];
%posterior estimate
% x= [0;0;0;5.0239e+06;1.6769e+06;-3.5422e+06];
x= [0;0;0;px;py;pz];

%to hold NEU velocity and position
X_neu= [0;0;0;0;0;0];
%output vector: n_y X 1
y = zeros(1,N);
z= zeros(3,N);
z2= zeros(3,114);
z(:,1) = x(4:6,1);
z2(:,1) = x(4:6,1);
%input vector:  n_u X 1




%Processed_attitude
%acceleration noise
a_noise = 0.25;

%process noise vector: n_x X 1
w = [T*a_noise;T*a_noise;T*a_noise;T^2*a_noise;T^2*a_noise;T^2*a_noise];
%measurement noise vector:  ny X 1
v = [0.01;0.01;0.01;0.01;0.01;0.01];
%Non-liner system state matrix : n_x X n_x
A = [1 0 0 0 0 0;
    0 1 0 0 0 0; 
    0 0 1 0 0 0; 
    0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1];

%System input matrix: n_x X n_u
G = [T      0   0
     0      T   0
     0      0   T
     T^2    0   0
     0      T^2 0
     0      0   T^2];
%observation matrix n_y X n_u
c = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];
%Linearized Observation matrix: n_y X n_u
H = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];


%position and noise 
pv_noise = [T*a_noise;T*a_noise;T*a_noise;T^2*a_noise;T^2*a_noise;T^2*a_noise];
%Process covariance matrix
Q =[(T*a_noise)^2 0 0 0 0 0;
   0 (T*a_noise)^2 0 0 0 0;
   0 0 (T*a_noise)^2 0 0 0;
   0 0 0 (T^2*a_noise)^2 0 0;
   0 0 0 0 (T^2*a_noise)^2 0;
   0 0 0 0 0 (T^2*a_noise)^2];

P=eye(6);

%measurement noise
meas_noise = 0.01; %this is metres
%Measurement covariance matrix

R = [meas_noise^2 0 0;
    0 meas_noise^2  0;
    0 0 meas_noise^2];



%%%other usefule variables%%%
I = eye(6);

%%%%%%%%%%%%%% Coord frames %%%%%%%%%%%%%%%
%Yaw,pitch and roll are resolved in body frame


%-----------Euler rotations from body to inertial-----------%
syms roll pitch yaw;
%yaw rotation
yaw_r = [cos(yaw) sin(yaw) 0;
        -sin(yaw) cos(yaw) 0;
         0          0      1 ];

 %pitch rotation
 pitch_r = [cos(pitch)  0   -sin(pitch);
            0           1   0
            sin(pitch)  0   cos(pitch)];

 %roll rotation
 roll_r = [1    0           0;
           0    cos(roll)   sin(roll);
           1    -sin(roll)  cos(roll)];

%transforamtion matrix from body to inertial fram
C_ib = roll_r*pitch_r*yaw_r;
C_bi = C_ib;

%----------- end Euler rotations from body to inertial-----------%

%-------------Converting from WGS84 ECEF to NEU Frame------------%
syms lat long alt

%raotions matrix from ECEF to navigatino frame NEU
% C_en = [-sin(lat)           cos(lat)            0;
%         -cos(lat)*sin(long) -sin(lat)*sin(long) cos(long);
%         cos(lat)*cos(long)  sin(lat)*cos(long)  sin(long);];
%     
%     
    
C_en = [-sin(long)           cos(long)            0;
        -cos(long)*sin(lat) -sin(lat)*sin(long) cos(lat);
        cos(lat)*cos(long)  sin(long)*cos(lat)  sin(lat);];

nav_vel = [0;0;0];
nav_pos = [0;0;0];
nav_gps = [0;0;0];
nav_true = [0;0;0];

C_en = subs(C_en,{lat,long},{-33.952548*pi/180,18.458266*pi/180});
nav_vel_init = C_en*x(1:3,1);
nav_pos_init = C_en*x(4:6,1);



%------------end of conversion from WGS84 ECEF to NEU-------------%





%%%%%%%%%%%%%% navigation equations %%%%%%%%%%%%%%%

lat_news = zeros(1,N);
long_news = zeros(1,N);
hs = zeros(1,N);
GPS = zeros(3,N);
GPS(:,1) = [-33.952548;18.458266;197.3];
k=1;

GPS_count =0;
GPS_count2 =0;
%intialise euler angles
    %yaw = data.Yaw(1);
    %pitch = data.Pitch(1);
    %roll = data.Roll(1);
% subs(C_bi ,{roll, pitch, yaw},{data.Roll,data.Pitch,data.Yaw(1)});
for k = 1:N
    
    %------------ get attitude ----------------------%
     %yaw = data.Yaw(k);
     %pitch = data.Pitch(k);
     %roll = data.Roll(k);
     
     C_bi = subs(C_bi ,{roll, pitch, yaw},{data.Roll(k),data.Pitch(k),data.Yaw(k)});
    
    %-----------interial frame acceleration-----------%
     ax = (data.XAccel(k)-0.010015263)*9.8;
     ay = (data.YAccel(k)+0.029243295)*9.8;
     az = (data.ZAccel(k)+0.112089019)*9.8;
     
    % subs(C_bi,{roll,pitch, yaw,}
    
     %u = double(C_bi)*[ax;ay;az];
     u = double(subs(C_bi ,{roll, pitch, yaw},{data.Roll(k),data.Pitch(k),data.Yaw(k)}))*[ax;ay;az];
     %compensate for gravity
      u = u +[0;0;9.8];
     %-----------end of interial frame acceleration----------%
   
     
     %----------------Prediciton stage---------------------% 
     x_hat(:,k+1) =   A* x(:,k) + G*(u);
     x_rms(:,k+1) =   A*x_rms(:,k)+ G*u;
     %error covariance prediction
     P = A*P*A'+ Q;
     %--------------end of Prediction stage----------------%
   
     %Kalman gain
     K = P*H'/(H*P*H'+R);
     %measurement
     %z = [data.X_meas(k);data.Y_meas(k);0];
     
     %-------- Transforming GPS(lat,long,h)->ECEF----------%   
     X = sqrt(1-(e^2)*sin(lat)^2);
     GPS_x = (a/X +alt)*cos(lat)*cos(long);
     GPS_y = (a/X + alt)*cos(lat)*sin(long);
     GPS_z = (a*(1-e^2)/X +alt)*sin(lat);
     GPS_x = subs(GPS_x,{lat,long,alt },{data.lat(k)*pi/180,data.long(k)*pi/180,data.height(k)});
     GPS_y = subs(GPS_y,{lat,long,alt },{data.lat(k)*pi/180,data.long(k)*pi/180,data.height(k)});
     GPS_z = subs(GPS_z,{lat,alt },{data.lat(k)*pi/180,data.height(k)});
     
     if(k<114)

         GPS_x2 = (a/X +alt)*cos(lat)*cos(long);
         GPS_y2 = (a/X + alt)*cos(lat)*sin(long);
         GPS_z2 = (a*(1-e^2)/X +alt)*sin(lat);
         GPS_x2 = subs(GPS_x2,{lat,long,alt },{truedata.lat(k)*pi/180,truedata.long(k)*pi/180,truedata.height(k)});
         GPS_y2 = subs(GPS_y2,{lat,long,alt },{truedata.lat(k)*pi/180,truedata.long(k)*pi/180,truedata.height(k)});
         GPS_z2 = subs(GPS_z2,{lat,alt },{truedata.lat(k)*pi/180,truedata.height(k)});
         z2(:,k+1) =[double(GPS_x2);double(GPS_y2);double(GPS_z2)];
     end      
        Z = [double(GPS_x);double(GPS_y);double(GPS_z)];
     if(GPS_count2 ==10)
        z(:,k+1)= [double(GPS_x);double(GPS_y);double(GPS_z)];
        GPS(:,k+1) = [data.lat(k),data.long(k),data.height(k)];
        GPS_count2 =0;
     else
         z(:,k+1) = z(:,k);
         GPS(:,k+1) = GPS(:,k);
         GPS_count2 = GPS_count2 +1;
     end
 
     %y(k) = double(GPS_x);
     %-------------------Innovation step--------------------%
     %assume 1Hz GPS update rate
    
     if(GPS_count == 1)
        x(:,k+1) = x_hat(:,k+1) +K*((Z) - c*x_hat(:,k+1));
        GPS_count =0;
        
     else
         x(:,k+1) = x_hat(:,k+1)+K*(x(4:6,k) - c*x_hat(:,k+1));
         GPS_count = GPS_count+1;
    
     end
     
            
     %take velocity and postion to navigation frame%
        
        %covert current EKF estimate from ECEF to lat long
    
        [lat_new ,long_new, h ]= ECEF_LLh(x(4,k+1),x(5,k+1),x(6,k+1));
        [lat_news(k) ,long_news(k), hs(k) ]= ECEF_LLh(x(4,k+1),x(5,k+1),x(6,k+1));
        %update transformation matrix 
        C_en = subs(C_en,{lat,long},{lat_new*pi/180,long_new*pi/180});
        %nav_vel(:,k+1) = C_en*x(1:3,k+1);
        nav_pos(:,k+1) = C_en*(x(4:6,k)- x(4:6,1));
        C_en = subs(C_en,{lat,long},{data.lat(k)*pi/180,data.long(k)*pi/180});
        nav_gps(:,k+1) = C_en*(z(:,k)-z(:,1));
        if(k<114)
          nav_true(:,k+1) = C_en*(z2(:,k)-z2(:,1));  
        end
     %update the state covariance
     P = (I-K*H)*P;
    
    
end 

%----------take the all ECEF data to navigation frame---------%




%Graphs
% predic_Px = x_hat(4,:);
% EKF_Px = x(4,:);
% %start(1:1,1:N) = 5.0235e+06;
% 
% figure;
% plot(y);
% title("X postion");
% 
% hold on
% 
% plot(EKF_Px);
% 
% legend({'predicted','EKF'},'location','southwest');
% hold off

%plotting north ,east up 


figure

plot(nav_pos(1,:),nav_pos(2,:));

hold on

plot(nav_gps(1,:),nav_gps(2,:));
plot(nav_true(1,:),nav_true(2,:));
 title("North EAST Map")
 ylabel("Distance East(m)")
 xlabel("Distance North(m)")
 legend({'EKF','GPS','true'},'location','southwest');
 hold off



figure

plot(nav_pos(2,:));
title("EAST")
xlabel("Sample Point")
ylabel("Distance(m)")
hold on 
plot(nav_gps(2,:));
legend({'EKF','GPS'},'location','southwest');
hold off

figure

plot(nav_pos(3,:));
title("UP")
xlabel("Sample Point")
ylabel("Distance(m)")
hold on
plot(nav_gps(3,:));
legend({'EKF','GPS'},'location','southwest');
hold off

%-------------ditance travelled----------------------%

%x direction in metres
distance = x(4,N)-x(4,1)



%------------- some useful functions ----------------%
function [LAT, LONG, height] = ECEF_LLh(x,y,z)
%this functions converts coordiantes in ECEF to geodetic lat long and
%height

%GPS constants: WGS84
%semi major axis
a = 6378137;
e= 0.08181919084;
b = 6356752.314;

w = sqrt(x^2+y^2);
l = (e^2)/2;
m = (w/a)^2;
n= [(1-e^2)*z/b]^2;
i = -(2*l^2 + m+n)/2;
k = (l^2)*(l^2 - m -n);
q = ((m+n-4*l^2)^3)/216 + m*n*l^2;
D = sqrt((2*q-m*n*l^2)*m*n*l^2);
beta = i/3 - nthroot(q+D,3)-nthroot(q-D,3);
t= sqrt(sqrt(beta^2-k)-(beta+i)/2) - sign(m-n)*sqrt((beta-i)/2);
w1 = w/(t+l);
z1 = (1-e^2)*z/(t-l);
LAT= atan(z1/((1-e^2)*w1));
LAT = rad2deg(LAT);
LONG = 2*atan((w-x)/y);
LONG = rad2deg(LONG);
height = sign(t-1+l)*sqrt((w-w1)^2+(z-z1)^2);
end 
