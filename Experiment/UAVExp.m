% Clear workspace and close figures
close all
clear all

% Load experiment log data
log = ulogreader("dataLogExp.ulg");
log.DropoutIntervals; % Check for data dropouts

% Read position, attitude, and angular velocity data
tab_xyz = readTopicMsgs(log,'TopicNames',{'vehicle_local_position'});
tab_rpy = readTopicMsgs(log,'TopicNames',{'vehicle_attitude'});
tab_pqr = readTopicMsgs(log,'TopicNames',{'vehicle_angular_velocity'});

data_xyz = tab_xyz.TopicMessages{:,:};
data_rpy = tab_rpy.TopicMessages{:,:}; 
data_rpy = data_rpy(2:2:end,:); % Downsample attitude data

data_pqr = tab_pqr.TopicMessages{:,:};
data_pqr = data_pqr(3:5:end,:); % Downsample angular velocity data
data_pqr = data_pqr.xyz;

tab_motor = readTopicMsgs(log,'TopicNames',{'actuator_outputs'});
tab_thrust = readTopicMsgs(log,'TopicNames',{'vehicle_thrust_setpoint'});
tab_torque = readTopicMsgs(log,'TopicNames',{'vehicle_thrust_setpoint'});

data_motor = tab_motor.TopicMessages{:,:};
data_thrust = tab_thrust.TopicMessages{:,:};
data_torque = tab_torque.TopicMessages{:,:};

% Extract UAV state variables
xact = data_xyz.x; yact = data_xyz.y; zact = data_xyz.z;
uact = data_xyz.vx; vact = data_xyz.vy; wact = data_xyz.vz;
[psiact, thetaact, phiact] = quat2angle(data_rpy.q); % Convert quaternion to Euler angles
pact = data_pqr(:,1); qact = data_pqr(:,2); ract = data_pqr(:,3);

% Extract actuator data
w1 = data_motor.output(:,1); w2 = data_motor.output(:,2);
w3 = data_motor.output(:,3); w4 = data_motor.output(:,4);
Tz = data_thrust.xyz(:,3); 
tphi = data_torque.xyz(:,1); 
ttheta = data_torque.xyz(:,2); 
tpsi = data_torque.xyz(:,3);

% Define constants
g = 9.81; % Gravity (m/s^2)
dt = 0.1; % Time step (s)
n = 12; % Number of state variables
r = 7;  % Number of unknown parameters

% Initialize state estimation variables
xhat = [phiact(1) thetaact(1) psiact(1) pact(1) qact(1) ract(1) ...
    uact(1) vact(1) wact(1) xact(1) yact(1) zact(1)]';
xhatArray = [];
thetahat = zeros(r,1);
thetahatArray = [];

% Define tuning parameters for estimation
a = 0.999;
lambda = 0.999;
Q = 5e3 * eye(n);
R = 1e-2 * eye(n);
P = 1e-2 * eye(n);
S = 1e-3 * eye(r);
Upsilon = zeros(n,r);

% Recursive estimation loop
for i = 1:length(psiact)
    xhatArray = [xhatArray xhat];
    thetahatArray = [thetahatArray thetahat];
    
    % Measurement and control inputs
    y = [phiact(i) thetaact(i) psiact(i) pact(i) qact(i) ract(i) ...
        uact(i) vact(i) wact(i) xact(i) yact(i) zact(i)]';
    u = [w1(i) w2(i) w3(i) w4(i)]';
    
    % System dynamics (simplified UAV equations of motion)
    f = [y(4)+y(2)*y(6); y(5)-y(1)*y(6); y(6)+y(1)*y(5); 0; 0; 0;
        y(6)*y(8)-y(5)*y(9)-g*y(2); y(4)*y(9)-y(6)*y(7)+g*y(1);
        g-y(4)*y(8)+y(5)*y(7); y(7)+y(2)*y(9)-y(3)*y(8);
        y(8)-y(1)*y(9)+y(3)*y(7); y(9)+y(1)*y(8)-y(2)*y(7)];
    
    % Approximation matrices
    Phi = [zeros(3,7); y(5)*y(6) -u(2)^2+u(4)^2 0 0 0 0 0;
        0 0 y(4)*y(6) -u(1)^2+u(3)^2 0 0 0;
        0 0 0 0 y(4)*y(5) -u(1)^2+u(2)^2-u(3)^2+u(4)^2 0; zeros(2,7);
        0 0 0 0 0 0 u(1)^2+u(2)^2+u(3)^2+u(4)^2; zeros(3,7)];
    
    % Update state covariances
    P = P + Q;
    Sigma = P + R;
    K = P / Sigma;
    P = (eye(n) - K) * P;
    
    % Update parameter covariances
    Omega = Upsilon + Phi;
    Upsilon = (eye(n) - K) * Upsilon + (eye(n) - K) * Phi;
    Lambda = inv(lambda * Sigma + Omega * S * Omega');
    Gamma = S * Omega' * Lambda;
    S = inv(lambda) * S - inv(lambda) * S * Omega' * Lambda * Omega * S;
    
    % Update estimations
    ytilde = y - xhat;
    Xi = K + Upsilon * Gamma; Psi = eye(n) - Xi;
    Q = a * Q + (1 - a) * Xi * ytilde * ytilde' * Xi';
    R = a * R + (1 - a) * Psi * ytilde * ytilde' * Psi';
    thetahat = thetahat + Gamma * ytilde;
    xhat = xhat + f * dt + Phi * thetahat + K * ytilde + Upsilon * Gamma * ytilde;
end

% Load and scale 3D UAV model for visualization
uav = stlread("x500Frame.stl");
vertices = uav.Points * 0.003;
faces = uav.ConnectivityList;
idx = round(length(xact) * 2/3) - 190;
uav_pos = [xact(idx) yact(idx) zact(idx)];
uav_ori = [phiact(idx) + deg2rad(90) thetaact(idx) psiact(idx)];

% Compute rotation matrix for UAV orientation
Rz = [cos(uav_ori(3)) -sin(uav_ori(3)) 0; sin(uav_ori(3)) cos(uav_ori(3)) 0; 0 0 1];
Ry = [cos(uav_ori(2)) 0 sin(uav_ori(2)); 0 1 0; -sin(uav_ori(2)) 0 cos(uav_ori(2))];
Rx = [1 0 0; 0 cos(uav_ori(1)) -sin(uav_ori(1)); 0 sin(uav_ori(1)) cos(uav_ori(1))];
R = Rz * Ry * Rx;
vertices = (R * vertices')' + uav_pos;

% Plot actual vs estimated trajectory
figure; plot3(xact,yact,zact,'k','LineWidth',4); hold on;
plot3(xhatArray(10,:),xhatArray(11,:),xhatArray(12,:),'b--','LineWidth',4);
patch('Faces',faces,'Vertices',vertices, ...
    'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.3 0.3 0.3]);
text(xact(1)+0.3,yact(1)-0.1,zact(1)-0.3, ...
    '$t_0$','FontSize',32,'Interpreter','latex')
text(xact(end)+0.3,yact(end)+0.1,zact(end)-0.3, ...
    '$t_f$','FontSize',32,'Interpreter','latex')
xlabel('$x\;(m)$','Interpreter','latex'); 
ylabel('$y\;(m)$','Interpreter','latex'); 
zlabel('$z\;(m)$','Interpreter','latex'); hold off;
grid on; xlim tight; ylim tight; set(gca,'FontSize',32); axis equal
legend('Actual Trajectory','Estimated Trajectory', 'FontSize', 32, ...
    'Location','northeast','Orientation','vertical');
view(740,20);

% % Uncomment for estimation error evolution (takes significant time)
% for j = 1:length(thetahatArray)
%     xest = [phiact(1) thetaact(1) psiact(1) pact(1) qact(1) ract(1) ...
%             uact(1) vact(1) wact(1) xact(1) yact(1) zact(1)]';
% 
%     for k = 1:length(psiact)
%         y = [phiact(k) thetaact(k) psiact(k) pact(k) qact(k) ract(k) ...
%             uact(k) vact(k) wact(k) xact(k) yact(k) zact(k)]';
%         u = [w1(k) w2(k) w3(k) w4(k)]';
% 
%         f = [y(4)+y(2)*y(6); y(5)-y(1)*y(6); 
%             y(6)+y(1)*y(5); 0; 0; 0;
%             y(6)*y(8)-y(5)*y(9)-g*y(2); 
%             y(4)*y(9)-y(6)*y(7)+g*y(1);
%             g-y(4)*y(8)+y(5)*y(7); 
%             y(7)+y(2)*y(9)-y(3)*y(8);
%             y(8)-y(1)*y(9)+y(3)*y(7); 
%             y(9)+y(1)*y(8)-y(2)*y(7)];
%         Phi = [zeros(3,7); y(5)*y(6) -u(2)^2+u(4)^2 0 0 0 0 0;
%             0 0 y(4)*y(6) -u(1)^2+u(3)^2 0 0 0;
%             0 0 0 0 y(4)*y(5) -u(1)^2+u(2)^2-u(3)^2+u(4)^2 0; 
%             zeros(2,7); 0 0 0 0 0 0 u(1)^2+u(2)^2+u(3)^2+u(4)^2; 
%             zeros(3,7)];
% 
%         xest = y + f * dt + Phi * thetahatArray(:,j);
%         re(:,k) = xest - y;
%     end
% 
%     nrmse(j) = norm(sqrt(mean(re.^2,2)));
%     nmae(j) = norm(mean(abs(re),2));
% end
% 
% fh = figure(2);
% fh.Position = [1500 500 750 250];
% % fh.WindowState = 'maximized';
% plot(nrmse,'LineWidth',4,'Color','#7E2F8E'); hold on; 
% plot(nmae,'LineWidth',4,'Color','#A2142F'); 
% grid on; xlim tight; set(gca,'FontSize',32);
% legend('RMSE Evolution','MAE Evolution');
% ylabel(['$\|\vec{e}(\vec{x},\hat{\vec{x}})\|$'], ...
%     'Interpreter','latex'); 
% xlabel('Iteration'); hold off;
