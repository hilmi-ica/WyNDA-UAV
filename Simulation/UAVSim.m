% Clear workspace and load data
close all
clear all

load("trajSim.mat")
load("ctrlSeqSim.mat")

% Define system parameters
g = 9.81; m = 5; l = 0.5; k = 1.5e-6; b = 1e-7;
Ix = 5e-3; Iy = 5e-3; Iz = 0.01;
params = [(Iy-Iz)/Ix;1/Ix;(Iz-Ix)/Iy;1/Iy;(Ix-Iy)/Iz;1/Iz;1/m];
dt = 0.1;  % Time step
n = 12;    % State dimension
r = 7;     % Number of parameters
rng(0)     % Set random seed

% Initialize states and estimation variables
xact = [0 0 0 0 0 0 0 0 0 0.5 0 0]';
xhat = [0 0 0 0 0 0 0 0 0 0.5 0 0]';
yact = xact;
xhatArray = [];
yactArray = [];

thetahat = zeros(r,1);
thetahatArray = [];

% Initialize filter parameters
a = 0.95;
lambda = 0.97;
Q = 0.001*eye(n);
R = 0.1*eye(n);
P = 0.01*eye(n);
S = 50*eye(r);
Upsilon = zeros(n,r);

% Main loop for parameter estimation
for i = 1:length(uHistory)
    xhatArray = [xhatArray xhat];
    yactArray = [yactArray yact];
    thetahatArray = [thetahatArray thetahat];

    % Process noise
    wn = [-0.02+0.04*rand(3,1); -0.08+0.16*rand(3,1);
        -0.05+0.1*rand(3,1); -0.05+0.1*rand(3,1)];
    
    % Update true state
    uact = uHistory(i,:)';
    xact = xact + QuadMod(xact, uact) * dt;
    yact = eye(n) * xact + wn;
    
    % Compute system dynamics
    f = [yact(4)+yact(2)*yact(6); yact(5)-yact(1)*yact(6); 
        yact(6)+yact(1)*yact(5); 0; 0; 0;
        yact(6)*yact(8)-yact(5)*yact(9)-g*yact(2); 
        yact(4)*yact(9)-yact(6)*yact(7)+g*yact(1);
        g-yact(4)*yact(8)+yact(5)*yact(7); 
        yact(7)+yact(2)*yact(9)-yact(3)*yact(8);
        yact(8)-yact(1)*yact(9)+yact(3)*yact(7); 
        yact(9)+yact(1)*yact(8)-yact(2)*yact(7)];
    
    % Compute wide-array approximation matrix
    Phi = [zeros(3,7); yact(5)*yact(6) uact(2) 0 0 0 0 0;
        0 0 yact(4)*yact(6) uact(3) 0 0 0;
        0 0 0 0 yact(4)*yact(5) uact(4) 0; 
        zeros(2,7); 0 0 0 0 0 0 uact(1); zeros(3,7)];

    % State covariance update
    P = P + Q;
    Sigma = P + R;
    K = P * inv(Sigma);
    P = (eye(n) - K) * P;

    % Parameter covariance update
    Omega = Upsilon + Phi;
    Upsilon = (eye(n) - K) * Upsilon + (eye(n) - K) * Phi;
    Lambda = inv(lambda * Sigma + Omega * S * Omega');
    Gamma = S * Omega' * Lambda;
    S = inv(lambda) * S - inv(lambda) * S * Omega' * Lambda * Omega * S;
    
    % Estimation update
    ytilde = yact - xhat;
    thetahat = thetahat + Gamma * ytilde;
    xhat = xhat + f * dt + Phi * thetahat + K * ytilde + Upsilon * Gamma * ytilde;
end

% Time vector
t = 0:dt:(i-1)*dt;

% Plot trajectories
figure(1);
plot3(x,y,z,'k:','LineWidth',2); hold on;
plot3(yactArray(10,:),yactArray(11,:),yactArray(12,:),'LineWidth',2);
plot3(xhatArray(10,:),xhatArray(11,:),xhatArray(12,:),'LineWidth',2, ...
    'Color','#77AC30','LineStyle','--');
plotBoxSim(yactArray(:,101)); plotBoxSim(yactArray(:,201)); 
xlabel('$x\;(m)$','Interpreter','latex'); 
ylabel('$y\;(m)$','Interpreter','latex'); 
zlabel('$z\;(m)$','Interpreter','latex'); grid on; axis equal;
legend('Desired Trajectory','Actual Trajectory','Estimated Trajectory', ...
    'Location','south','Orientation','horizontal');

% Plot parameter estimation results
figure(2);
subplot(4,2,1); plot(t,(Iy-Iz)/Ix*ones(1,length(t)),'k','LineWidth',2)
hold on; plot(t,thetahatArray(1,:)/dt,'--','LineWidth',2,'Color',"#77AC30"); 
grid minor; xlim tight; set(gca,'FontSize',14); 
ylabel('$\frac{I_{yy}-I_{zz}}{I_{xx}}$','Interpreter','latex'); hold off
subplot(4,2,2); plot(t,1/Ix*ones(1,length(t)),'k','LineWidth',2)
hold on; plot(t,thetahatArray(2,:)/dt,'--','LineWidth',2,'Color',"#77AC30"); 
grid minor; xlim tight; set(gca,'FontSize',14);
ylabel('$\frac{1}{I_{xx}}$','Interpreter','latex'); hold off
subplot(4,2,3); plot(t,(Iz-Ix)/Iy*ones(1,length(t)),'k','LineWidth',2)
hold on; plot(t,thetahatArray(3,:)/dt,'--','LineWidth',2,'Color',"#77AC30");
grid minor; xlim tight; set(gca,'FontSize',14);
ylabel('$\frac{I_{zz}-I_{xx}}{I_{yy}}$','Interpreter','latex'); hold off
subplot(4,2,4); plot(t,1/Iy*ones(1,length(t)),'k','LineWidth',2)
hold on; plot(t,thetahatArray(4,:)/dt,'--','LineWidth',2,'Color',"#77AC30"); 
grid minor; xlim tight; set(gca,'FontSize',14);
ylabel('$\frac{1}{I_{yy}}$','Interpreter','latex'); hold off
subplot(4,2,5); plot(t,(Ix-Iy)/Iz*ones(1,length(t)),'k','LineWidth',2)
hold on; plot(t,thetahatArray(5,:)/dt,'--','LineWidth',2,'Color',"#77AC30"); 
grid minor; xlim tight; set(gca,'FontSize',14);
ylabel('$\frac{I_{xx}-I_{yy}}{I_{zz}}$','Interpreter','latex'); hold off
subplot(4,2,6); plot(t,1/Iz*ones(1,length(t)),'k','LineWidth',2)
hold on; plot(t,thetahatArray(6,:)/dt,'--','LineWidth',2,'Color',"#77AC30"); 
grid minor; xlim tight; ylim([-100 100]); set(gca,'FontSize',14);
ylabel('$\frac{1}{I_{zz}}$','Interpreter','latex'); hold off
subplot(4,2,[7,8]); plot(t,1/m*ones(1,length(t)),'k','LineWidth',2)
hold on; plot(t,thetahatArray(7,:)/dt,'--','LineWidth',2,'Color',"#77AC30"); 
grid minor; xlim tight; xlabel('$t\;(s)$','Interpreter','latex'); 
legend({'\textbf{Parameters}','\boldmath$\theta$'},'Interpreter' ...
    ,'latex','Location','southeast','Orientation','horizontal')
ylabel('$\frac{1}{m}$','Interpreter','latex'); hold off
set(gca,'FontSize',14);
