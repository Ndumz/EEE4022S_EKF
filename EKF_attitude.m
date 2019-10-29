%Extended kalman filter : For attitude estimation

%Input data containing IMU data
data = readtable('pitch_2_30_45.csv');
%number of samles
N = height(data);

%sampling interval
T = 0.0624422;
%T = 1;
%Extended kalamn filter states
syms roll_s pitch_s yaw_s;

%intial states: assuming the IMU is level and pointing true north
roll = 0;
pitch = 0;%-31.6896283091164;
yaw = 0;

%state matrix apior
x_hat=[roll; pitch; yaw];
%x_hat = zeros(3,N);
%state matrix posterior
x =   [roll;pitch;yaw];
%x = zeros(3,N);
%body to interial frame transformation for angular rates
C_bi = [1 sin(roll_s)*tan(pitch_s) cos(roll_s)*tan(pitch_s);
        0 cos(roll_s)            -sin(roll_s);
        0 sin(roll_s)*sec(pitch_s) cos(roll_s)*sec(pitch_s)];
    
%state transition matrix
A = eye(3);

%input matrix
B = C_bi;


%Measurement matrix
H = eye(3);
y  = [0;0;0];
y_triad =[0;0;0];


%Process noise Covariance: Gyro has 0.06deg/s
%noise covariance
noise_q = (0.06);
Q = noise_q*eye(3);

%Measurement covariance: this is an estimate
R = 2*eye(3);

%Error covariance matrix
P =5*eye(3);
%non linear state transition

for k = 1:N-1
    %Reading in inputs from gyro;
    P_r = data.XGyro(k+1);%-0.018;
    Q_r = data.YGyro(k+1);%-0.013;
    R_r = data.ZGyro(k+1);%-0.008;
    u = [P_r;Q_r;R_r];
    
    %-------------------- Prediction stage --------------------%
     B = subs(B,{roll_s,pitch_s,yaw_s},{x(1,k), x(2,k),x(3,k)});
    x_hat(:,k+1) = A*x(:,k) + double(subs(B,{roll_s,pitch_s,yaw_s},{x(1,k), x(2,k),x(3,k)}))*u*T;
    G = T*B*u;
    f =[G(1,:),G(2,:),G(3,:)];
    v = [yaw_s pitch_s roll_s];
    %linearize nonllinear state equations
    F = jacobian(f,v);
    
    F= subs(F,{roll_s, pitch_s ,yaw_s},{x(1,k), x(2,k) ,x(3,k)});
    S= 300*eye(3);
    
    %predict covariance matrix
    P = (double(F)*P*double(F'))+Q;
    
    %------------------- end Prediction stage -------------------%
    
    

     
     %------------------- traid algorithm ------------------------%
     
     %body frame vectors
     a_b = [data.XAccel(k)-0.01,data.YAccel(k)+0.0292,data.ZAccel(k)+0.11208]';
     a_b= (a_b)/norm(a_b);
     m_b =[data.XMag(k),data.YMag(k),data.ZMag(k)]'; 
     m_b = (m_b)/norm(m_b);
     
     %inertail frame vectors
     a_i = [0,0,1]';
     a_i= (a_i)/norm(a_i);
     m_i = [9486*10^-9,0,-230764.4*10^-6]';
     m_i = m_i/norm(m_i);
     
     t_1b = a_b;
     t_2b = (cross(a_b,m_b))/(norm(cross(a_b,m_b)));
     t_3b = cross(t_1b,t_2b);
     
     t_1i = a_i;
     t_2i = (cross(a_i,m_i))/(norm(cross(a_i,m_i)));
     t_3i = cross(t_1i,t_2i);
     
     t_b = [t_1b t_2b t_3b];
     t_i = [t_1i t_2i t_3i];
     
     R_bi = t_b*t_i';
     %roll , pitch , yaw
     r =180* atan(R_bi(2,3)/R_bi(3,3))/pi;
     p = 180*asin(R_bi(1,3))/pi;
     q = 180*atan(R_bi(1,2)/R_bi(1,1))/pi;
     
     %absolute measurement
     y(:,k) = [r;p;q];
     
     %---------------- end traid algorithm ------------------------%
     
     %kalman Gain
     K = P*H'/(H*P*H'+R);
     
     %--------------- innovation step -----------------------%
     x(:,k+1)= x_hat(:,k+1 )+K*(y(:,k) - H*x_hat(:,k+1));
     
     %update covariance matrix
     P = (eye(3)-K*H)*P';
    
end 

% predic_roll  = x_hat(1,:);
% predic_pitch = x_hat(2,:);
% predic_yaw = x_hat(3,:);
% actual_pitch = data.actual;

meas_roll = y(1,:);
meas_pitch = y(2,:);
meas_yaw = y(3,:);

EKF_roll = x(1,:);
EKF_pitch = x_hat(2,:);
EKF_yaw = x(3,:);

t = linspace(0,1,808);

figure(1);
plot(EKF_roll)
title("Roll");

hold on

plot(meas_roll)
xlabel("Sample point");
ylabel("Angle(degrees)");

% plot(EKF_roll)

legend({'EKF','TRIAD'},'Location','southwest')

hold off
h1 = figure(1);
print(h1,'roll_accel','-depsc') 

figure(2);
plot(EKF_pitch)
title("Pitch");

hold on

plot(meas_pitch)

%plot(data.actual)

xlabel("Sample point");
ylabel("Angle(degrees)");
%plot(EKF_pitch)

% legend({'Actual','measured','EKF'},'Location','southwest')
legend({'EKF','TRIAD'},'Location','southwest')

hold off
% h = figure(2);
 %saveas(h,'improvedPitchD.pdf');
 h2=figure(2);
%set(h,'PaperSize',[10 2]); %set the paper size to what you want  
print(h2,'pitch_accel','-depsc') % then print it

 
 
 figure;
plot(EKF_yaw)
title("Yaw");

hold on

plot(meas_yaw)

% plot(EKF_yaw)

% legend({'predicted','measured','EKF'},'Location','southwest')
legend({'EKF','TRIAD'},'Location','northhwest')

hold off

%plot graphs




%%%old code
%%code to get roll , pitch , yaw
%     m_roll = -180 * atan2(data.YAccel(k), sqrt(data.XAccel(k)*data.XAccel(k) + data.ZAccel(k)*data.ZAccel(k)))/pi;
%     m_pitch = -180 * atan2(data.YAccel(k), sqrt(data.YAccel(k)*data.YAccel(k) + data.ZAccel(k)*data.ZAccel(k)))/pi;
%     
%     mag_x = data.XMag(k)*cos(m_pitch) + data.YMag(k)*sin(m_roll)*sin(m_pitch) + data.ZMag(k)*cos(m_roll)*sin(m_pitch);
%     mag_y = data.YMag(k) * cos(m_roll) - data.ZMag(k) * sin(m_roll);
%     m_yaw = 180 * atan2(-mag_y,mag_x)/pi;
%     
%      y(:,k) = [m_roll;m_pitch;m_yaw];