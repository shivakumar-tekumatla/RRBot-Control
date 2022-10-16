syms t 
%% 4.2.1.a)Cubic polynomical generation 
t0=0;tf=10;
theta1_0=pi;theta1_f=0;theta1_dot_0=0;theta1_dot_f=0; %Initial and final conditions for first joint
theta2_0=pi/2;theta2_f=0;theta2_dot_0=0;theta2_dot_f=0;  %Initial and final conditions for second joint

theta1_desired = [ 1  t t^2 t^3]*cubic(t0,tf,theta1_0,theta1_f,theta1_dot_0,theta1_dot_f);
theta2_desired = [ 1  t t^2 t^3]*cubic(t0,tf,theta2_0,theta2_f,theta2_dot_0,theta2_dot_f);
theta1_dot_desired = jacobian(theta1_desired,t);
theta2_dot_desired = jacobian(theta2_desired,t);
theta1_ddot_desired = jacobian(theta1_dot_desired,t);
theta2_ddot_desired = jacobian(theta2_dot_desired,t);

%% Cubic polynomial solver Function
function constants=cubic(t0,tf,theta_0,theta_f,theta_dot_0,theta_dot_f)
    A = [1  t0 t0^2 t0^3 ;
        0  1  2*t0 3*t0^2;
        1  tf tf^2 tf^3  ;
        0  1  2*tf 3*tf^2];
    constants=inv(A)*[theta_0;theta_dot_0;theta_f;theta_dot_f];
    
end