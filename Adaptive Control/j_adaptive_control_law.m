%% Adaptive control law
n= 2; %Number of generalized coordinates

A = [zeros(n,n) eye(n);
    zeros(n,n) zeros(n,n)];
B = [zeros(n,n);
    eye(n)];

eigen_values= [-3,-3,-4,-4]; 
K = place(A,B,eigen_values);
Kp = K(1:2,1:2);
Kd = K(1:2,3:4);

Acl = [zeros(n,n) eye(n);
        -Kp        -Kd];
    
 Q = eye(size(A));
 
 P = lyap(Acl',Q);
 
 gamma = eye(5).*0.5;