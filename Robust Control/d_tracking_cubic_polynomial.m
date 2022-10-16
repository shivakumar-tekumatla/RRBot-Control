clc;clear;close all;
%% Generalized coordinates 
% z =[theta1;theta2]
n = 2; %Number of generaized coordinates;
%%
c_robust_control_law;
 
%  Uncertainity = inv(Mmat)*Mmat_hat - eye(size(M))
%% ODE Implementation 
phi = 0;%0.01;
v=[0;0]; % Intial torques 
rho_0 = 0.5; % this will be added to the uncertainity to compute rho 
[time,y] = ode45(@(t,x) robust(t,x,Kp,Kd,P,B,phi,v,rho_0), [0,10],[deg2rad(180);deg2rad(90);0;0]);
 
%% Plotting  
 theta1_desired=zeros(size(time));
 theta2_desired=zeros(size(time));
 theta1_dot_desired=zeros(size(time));
 theta2_dot_desired=zeros(size(time));
 u1 = zeros(size(time));
 u2 = zeros(size(time));
 e1 = zeros(size(time));
 e2 = zeros(size(time));
 e1_dot = zeros(size(time));
 e2_dot = zeros(size(time));
 for i = 1:length(time)
     t = time(i);
     
     theta1_desired(i) =(pi*t^3)/500 - (3*pi*t^2)/100 + pi;
     theta2_desired(i) =  (pi*t^3)/1000 - (3*pi*t^2)/200 + pi/2;
     theta1_dot_desired(i) = (3*pi*t^2)/500 - (3*pi*t)/50;
     theta2_dot_desired(i) =    (3*pi*t^2)/1000 - (3*pi*t)/100;
     xref = [theta1_desired(i);theta2_desired(i);theta1_dot_desired(i);theta2_dot_desired(i)];
     x = [y(i,1);y(i,2);y(i,3);y(i,4)];
     q1 = x(1);
     q2 = x(2);
     dq1 = x(3);
     dq2 = x(4);
     e = x-xref;
     e1(i) = e(1);
     e2(i) = e(2);
     e1_dot(i) = e(3);
     e2_dot(i) = e(4);
     %Actual physical parameters 
    Mmat  = [(9*cos(q2))/10 + 1573/1000  (9*cos(q2))/20 + 573/2000;
             (9*cos(q2))/20 + 573/2000                   573/2000];
    Cmat =[-(9*dq2*sin(q2))/20, -(9*sin(q2)*(dq1 + dq2))/20;
        (9*dq1*sin(q2))/20,                           0];
    Gmat =[- (8829*sin(q1 + q2))/2000 - (28449*sin(q1))/2000;
                        -(8829*sin(q1 + q2))/2000];
    %Nominal parameters
    
    Mmat_hat =[(27*cos(q2))/40 + 4719/4000, (27*cos(q2))/80 + 1719/8000;
        (27*cos(q2))/80 + 1719/8000,                   1719/8000];
    Cmat_hat =[-(27*dq2*sin(q2))/80, -(27*sin(q2)*(dq1 + dq2))/80;
            (27*dq1*sin(q2))/80,                            0];
    Gmat_hat =[- (26487*sin(q1 + q2))/8000 - (85347*sin(q1))/8000;
                        -(26487*sin(q1 + q2))/8000];
     ddqd = [(3*pi*t)/250 - (3*pi)/50;
            (3*pi*t)/500 - (3*pi)/100];
    Cmat_tilde = Cmat_hat-Cmat;
    Gmat_tilde = Gmat_hat-Gmat;

    uncertainity = (inv(Mmat)*Mmat_hat-eye(size(Mmat)))*v+inv(Mmat)*(Cmat_tilde*[dq1;dq2]+Gmat_tilde);
%     disp(uncertainity);
    rho = norm(uncertainity)+rho_0; 
    if phi>0
        if norm(e'*P*B)>phi
            vr = -rho*(e'*P*B)/(norm(e'*P*B));
        else
            vr = -rho*(e'*P*B)/phi;
        end
    else
        if norm(e'*P*B)~=0
            vr = -rho*(e'*P*B)/(norm(e'*P*B));
        else
            vr=0;
        end
    end
     
     v = ddqd- K*e +vr';
%     disp(v)
    tau = Mmat_hat*v + Cmat_hat*[dq1;dq2]+Gmat_hat;
    u1(i) =tau(1);
    u2(i) =tau(2);
 end
 
 
 
figure(1);
plot(time,y,'linewidth',2);
title('Time vs Theta1, Theta2, Theta1-dot, Theta2-dot')
lgd = legend('Theta1', 'Theta2', 'Theta1-dot', 'Theta2-dot');
lgd.FontSize = 14;

figure(2);
subplot(2,2,1)      
plot(time,y(:,1),"b-",time,theta1_desired,"r--",'linewidth',2);   
lgd = legend('Theta-1','Theta1-desired');
lgd.FontSize = 14;
title('Time vs Theta-1')
xlabel("Time in Seconds")
ylabel("Angle in radians")

subplot(2,2,2)      
plot(time,y(:,2),"b-",time,theta2_desired,"r--",'linewidth',2);   
lgd = legend('Theta-2','Theta2-desired');
lgd.FontSize = 14;
title('Time vs Theta-2')
xlabel("Time in Seconds")
ylabel("Angle in radians")

subplot(2,2,3)      
plot(time,y(:,3),"b-",time,theta1_dot_desired,"r--",'linewidth',2);   
lgd = legend('Theta-1-dot','Theta1-dot-desired');
lgd.FontSize = 14;
title('Time vs Theta-1-dot')
xlabel("Time in Seconds")
ylabel("Velocity in radians per second")

subplot(2,2,4)      
plot(time,y(:,4),"b-",time,theta2_dot_desired,"r--",'linewidth',2);   
lgd = legend('Theta-2-dot','Theta2-dot-desired');
lgd.FontSize = 14;
title('Time vs Theta-2-dot')
xlabel("Time in Seconds")
ylabel("Velocity in radians per second")

figure(3);
subplot(2,1,1);
plot(time,u1,'linewidth',2);
title("Joint-1 Torque")
lgd = legend("Joint-1 Torque");
lgd.FontSize = 14;
xlabel("Time in Seconds")
ylabel("Torque in N-m")
subplot(2,1,2);
plot(time,u2,'linewidth',2);
title("Joint-2 Torque")
lgd = legend("Joint-2 Torque");
lgd.FontSize = 14;
xlabel("Time in Seconds")
ylabel("Torque in N-m")

figure(4)
plot(time,[e1,e2,e1_dot,e2_dot],'linewidth',2);
title('Time vs Errors')
lgd = legend('Theta1-Error', 'Theta2-Error', 'Theta1-dot-Error', 'Theta2-dot-Error');
lgd.FontSize = 14;


%% Robust control law ODE 

function [dx,tau] = robust(t,x,Kp,Kd,P,B,phi,v,rho_0)

    q1 = x(1);
    q2 = x(2);
    dq1 = x(3);
    dq2 = x(4);
    
    if abs(q1) > 2*pi
        q1 = mod(q1, 2*pi);
    end 

    if abs(q2) > 2*pi
        q2 = mod(q2, 2*pi);

    end

    %Actual physical parameters 
    Mmat  = [(9*cos(q2))/10 + 1573/1000  (9*cos(q2))/20 + 573/2000;
             (9*cos(q2))/20 + 573/2000                   573/2000];
    Cmat =[-(9*dq2*sin(q2))/20, -(9*sin(q2)*(dq1 + dq2))/20;
        (9*dq1*sin(q2))/20,                           0];
    Gmat =[- (8829*sin(q1 + q2))/2000 - (28449*sin(q1))/2000;
                        -(8829*sin(q1 + q2))/2000];
    %Nominal parameters
    
    Mmat_hat =[(27*cos(q2))/40 + 4719/4000, (27*cos(q2))/80 + 1719/8000;
        (27*cos(q2))/80 + 1719/8000,                   1719/8000];
    Cmat_hat =[-(27*dq2*sin(q2))/80, -(27*sin(q2)*(dq1 + dq2))/80;
            (27*dq1*sin(q2))/80,                            0];
    Gmat_hat =[- (26487*sin(q1 + q2))/8000 - (85347*sin(q1))/8000;
                        -(26487*sin(q1 + q2))/8000];
        
    Cmat_tilde = Cmat_hat-Cmat;
    Gmat_tilde = Gmat_hat-Gmat;

    uncertainity = (inv(Mmat)*Mmat_hat-eye(size(Mmat)))*v+inv(Mmat)*(Cmat_tilde*[dq1;dq2]+Gmat_tilde);
%     disp(uncertainity);
    rho = norm(uncertainity)+rho_0; 
    disp(rho)
%     rho=5;
    disp("___________")
    
    qd =[(pi*t^3)/500 - (3*pi*t^2)/100 + pi;
        (pi*t^3)/1000 - (3*pi*t^2)/200 + pi/2];
    dqd = [(3*pi*t^2)/500 - (3*pi*t)/50;
        (3*pi*t^2)/1000 - (3*pi*t)/100];
    
    ddqd = [(3*pi*t)/250 - (3*pi)/50;
            (3*pi*t)/500 - (3*pi)/100];
        
    xref = [qd;dqd];
    

    e = x-xref;
    
    K = [Kp Kd];
    if phi>0
        if norm(e'*P*B)>phi
            vr = -rho*(e'*P*B)/(norm(e'*P*B));
        else
            vr = -rho*(e'*P*B)/phi;
        end
    else
        if norm(e'*P*B)~=0
            vr = -rho*(e'*P*B)/(norm(e'*P*B));
        else
            vr=0;
        end
    end
    v = ddqd- K*e +vr';
%     disp(v)
    tau = Mmat_hat*v + Cmat_hat*[dq1;dq2]+Gmat_hat;
%     disp(vr);
    
    ddq = inv(Mmat)*(tau-Cmat*[dq1;dq2]-Gmat);
    
    dx(1,1) = dq1;
    dx(2,1) = dq2;
    dx(3,1) = ddq(1);
    dx(4,1) = ddq(2);
    
end 
    
%% Cubic polynomial solver 

function constants=cubic(t0,tf,theta_0,theta_f,theta_dot_0,theta_dot_f)
    A = [1  t0 t0^2 t0^3 ;
        0  1  2*t0 3*t0^2;
        1  tf tf^2 tf^3  ;
        0  1  2*tf 3*tf^2];
    constants=inv(A)*[theta_0;theta_dot_0;theta_f;theta_dot_f];
    
end
