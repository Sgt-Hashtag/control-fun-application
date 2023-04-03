clear; close all;
load("measurement_data.mat")
tire_radius=0.145;
L=tire_radius*n_L.Data;
R=tire_radius*n_R.Data;
timeStep=IMU_GZ.Time(2)-IMU_GZ.Time(1);
Yaw=cumtrapz(IMU_GZ.Time,IMU_GZ.Data);
P=[0 0]';

Vel_list=[];
for b =1:size(n_L.Time,1)
    Vel=(L(b)+R(b))/2;
    Vel_comp=Vel*[cos(Yaw(b)) sin(Yaw(b))]';
    Vel_list = [Vel_list;Vel];
    pos=P(:,end) + Vel_comp*timeStep;
    P = [P,pos];
end


subplot(2,2,1);
plot(P(1,:),P(2,:));
xlabel("x(m)"),ylabel("y(m)");
title("Position")

subplot(2,2,2);
plot(n_L.Time,Vel_list);
xlabel("Time(sec)"),ylabel("Velocity(m/s)")
title("Velocity Profile")

subplot(2,2,3)
plot(n_L.Time,n_L.Data)
hold on
plot(n_R.Time,n_R.Data)
hold off
xlabel("Time(sec)"),ylabel("Angular Velocity(m/s)")
title("Angular Velocity")

subplot(2,2,4)
plot(IMU_GZ.Time,Yaw)
xlabel("Time(sec)"),ylabel("Yaw")
title("Cumulative Yaw")

% function [pf,theta2] = orient(P,b,x,IMU_GZ,n_L,n_R,A)
%     theta1=x;
%     theta2=theta1+getdatasamples(IMU_GZ,b)*IMU_GZ.Time(b);
%     if getdatasamples(n_L,b)>getdatasamples(n_R,b)
%         c=P-getdatasamples(A,b)*[-sin(theta1) cos(theta1)]';
%         pf=c+getdatasamples(A,b)*[-sin(theta2) cos(theta2)]';
%     elseif getdatasamples(n_L,b)<getdatasamples(n_R,b)
%         c=P-getdatasamples(A,b)*[sin(theta1) -cos(theta1)]';
%         pf=c+getdatasamples(A,b)*[sin(theta2) -cos(theta2)]';
%     elseif getdatasamples(n_L,b)==getdatasamples(n_R,b)
%         pf=P+getdatasamples(X,b)*[sin(theta2) cos(theta2)]';
%     end
% end
