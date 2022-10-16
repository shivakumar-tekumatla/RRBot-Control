clc; clear

% ROS Setup
rosinit;
j1_effort = rospublisher('/rrbot/joint1_effort_controller/command');
j2_effort = rospublisher('/rrbot/joint2_effort_controller/command');
JointStates = rossubscriber('/rrbot/joint_states');
tau1 = rosmessage(j1_effort);
tau2 = rosmessage(j2_effort);
tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);
client = rossvcclient('/gazebo/set_model_configuration');
req = rosmessage(client);
req.ModelName = 'rrbot';
req.UrdfParamName = 'robot_description';
req.JointNames = {'joint1','joint2'};
req.JointPositions = [deg2rad(200), deg2rad(125)];
resp = call(client,req,'Timeout',3);

t = 0;

m1 = 1; m2 = 1; %kg
m1_hat = .75; m2_hat = .75; %kg

l1 = 1; l2 = 1; r1 = .45; r2 = .45; %m

I1 = .084; I2 = .084; %kq*m^2
I1_hat = .063; I2_hat = .063; %kg*m^2

g = 9.81; %m/s^2

old_t= 0;
old_dq1 = 0;
old_dq2 = 0;

alpha = [m2*l1^2 + m1*r1^2 + m2*r2^2 + I1 + I2;
         m2*l1*r2
         m2*r2^2 + I2
         m1*r1 + m2*l1
         m2*r2];

alpha = .75*alpha;

K = [30  0 11  0;
      0 30  0 11];

a1 = [ 3.1416;
            0;
      -0.0942;
       0.0063];

a2 = [ 1.5708;
            0;
      -0.0471;
       0.0031];

B = [0 0;
     0 0;
     1 0;
     0 1];

P =[1.5924         0    0.0167         0;
         0    1.5924         0    0.0167;
    0.0167         0    0.0470         0;
         0    0.0167         0    0.0470];  

Gamma = eye(5).*10;

Timey = [];
DesiredTraj_pos = [];
DesiredTraj_vel = [];
MeasuredTraj = [];
Torques = [];
Error = [];
Uncertainties = [];

tic
while(t < 10)

    t = toc;
    % read the joint states
    jointData = receive(JointStates);
    % inspect the "jointData" variable in MATLAB to get familiar with its structure
    % design your state feedback controller in the following
    q1 = jointData.Position(1);
    dq1 = jointData.Velocity(1);
    q2 = jointData.Position(2);
    dq2 = jointData.Velocity(2);

    X = [q1;q2;dq1;dq2];

    q1d = a1(1) + a1(2)*t + a1(3)*t^2 + a1(4)*t^3;
    dq1d = a1(2) + a1(3)*2*t + a1(4)*3*t^2;
    ddq1d = a1(3)*2 + a1(4)*6*t;

    q2d = a2(1) + a2(2)*t + a2(3)*t^2 + a2(4)*t^3;
    dq2d = a2(2) + a2(3)*2*t + a2(4)*3*t^2;
    ddq2d = a2(3)*2 + a2(4)*6*t;

    qd = [q1d; 
          q2d];
    dqd = [dq1d; 
           dq2d];
    ddqd = [ddq1d;
            ddq2d];

    Xref = [qd;dqd];

    e = X - Xref;
    e(1:2,:) = wrapToPi(e(1:2,:));

    v = ddqd - K*e;

    % uncertain values
    Y = [v(1), ...
         cos(q2)*(2*v(1)+v(2))-2*sin(q2)*dq1*dq2-sin(q2)*dq2^2, ...
         v(2), ...
         -sin(q1)*g, ...
         -sin(q1+q2)*g; ...
         0, ...
         sin(q2)*dq1^2+cos(q2)*v(1), ...
         v(1)+v(2), ...
         0, ...
         -sin(q1+q2)*g];

    T = Y*alpha;

    %simulate with measured values
    deltaT = t-old_t;
    deltaV = [dq1-old_dq1;
              dq2-old_dq2];
    ddq = deltaV/deltaT;

    Y = [ddq(1), ...
         cos(q2)*(2*ddq(1)+ddq(2))-2*sin(q2)*dq1*dq2-sin(q2)*dq2^2, ...
         ddq(2), ...
         -sin(q1)*g, ...
         -sin(q1+q2)*g; ...
         0, ...
         sin(q2)*dq1^2+cos(q2)*ddq(1), ...
         ddq(1)+ddq(2), ...
         0, ...
         -sin(q1+q2)*g];

    a = alpha(1);
    b = alpha(2);
    d = alpha(3);
    M_hat = [a+2*b*cos(q2), d+b*cos(q2); d+b*cos(q2), d];

    Phi = M_hat\Y;
    dalpha = -Gamma\(Phi'*B'*P*e);

    alpha = alpha + dalpha*deltaT;

    tau1.Data = T(1);
    tau2.Data = T(2);
    send(j1_effort,tau1);
    send(j2_effort,tau2);
    
    % you can sample data here to be plotted at the end
    Timey = [t Timey];
    DesiredTraj_pos = [qd DesiredTraj_pos];
    DesiredTraj_vel = [dqd DesiredTraj_vel];
    MeasuredTraj = [X MeasuredTraj];
    Torques = [T Torques];
    Error = [e Error];
    Uncertainties = [alpha Uncertainties];

    old_t= t;
    old_dq1 = dq1;
    old_dq2 = dq2;

end

tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);

% disconnect from roscore
rosshutdown;

figure,
C = tiledlayout(2,1);
nexttile
hold on
plot(Timey, MeasuredTraj(1,:), "red")
plot(Timey, DesiredTraj_pos(1,:), 'blue')
title('theta1 trajectory')
xlabel('time (s)')
ylabel('theta (rad)')
nexttile
hold on
plot(Timey, MeasuredTraj(2,:), "red")
plot(Timey, DesiredTraj_pos(2,:), 'blue')
title('theta2 trajectory')
xlabel('time (s)')
ylabel('theta (rad)')

figure,
B = tiledlayout(2,1);
nexttile
hold on
plot(Timey, MeasuredTraj(3,:), "red")
plot(Timey, DesiredTraj_vel(1,:), 'blue')
title('theta1d')
xlabel('time (s)')
ylabel('thetad (rad/s)')
nexttile
hold on
plot(Timey, MeasuredTraj(4,:), "red")
plot(Timey, DesiredTraj_vel(2,:), 'blue')
title('theta2d')
xlabel('time (s)')
ylabel('thetad (rad/s)')

figure,
B = tiledlayout(2,1);
nexttile
hold on
plot(Timey, Error(1,:), "red")
plot(Timey, Error(2,:), 'blue')
title('Position Error')
xlabel('time (s)')
ylabel('Angle (rad)')
nexttile
hold on
plot(Timey, Error(3,:), "red")
plot(Timey, Error(4,:), 'blue')
title('Velocity Error')
xlabel('time (s)')
ylabel('Angular Velocity (rad/s)')

figure,
C = tiledlayout(2,1);
nexttile
plot(Timey, Torques(1,:), "red")
title('Torque Joint 1')
xlabel('time (s)')
ylabel('torque (N*m)')
nexttile
plot(Timey, Torques(2,:), "blue") 
title('Torque Joint 2')
xlabel('time (s)')
ylabel('torque (N*m)')
