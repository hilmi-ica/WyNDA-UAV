% Clear workspace
clear all; 
close all;

% Load control data
load("ctrlSeqSim.mat");

% Define system P\parameters
g = 9.81; m = 5; l = 0.5; k = 1.5e-6; b = 1e-7;
Ix = 5e-3; Iy = 5e-3; Iz = 0.01;
true_coeffs = [(Iy-Iz)/Ix; 1/Ix; ...
               (Iz-Ix)/Iy; 1/Iy; ...
               (Ix-Iy)/Iz; 1/Iz; 1/m];

% Simulation settings
dt = 0.1;     
n = 12; 
rng(0);         

% Initialize states
xact = [0 0 0 0 0 0 0 0 0 0.5 0 0]';
yact = xact;

% Simulate system and collect data
for i = 1:length(uHistory)
    % Measurement noise 
    wn = [-0.02+0.04*sum(rand(3,100),2)/100; 
          -0.08+0.16*sum(rand(3,100),2)/100;
          -0.05+0.1*sum(rand(3,100),2)/100; 
          -0.05+0.1*sum(rand(3,100),2)/100];
    
    uact(:,i) = uHistory(i,:)';
    xact = xact + QuadMod(xact, uact(:,i)) * dt;
    yact(:,i) = xact + wn; 
    
    % True model dynamics (used to compute residuals)
    f(:,i) = [ yact(4,i) + yact(2,i)*yact(6,i);
               yact(5,i) - yact(1,i)*yact(6,i);
               yact(6,i) + yact(1,i)*yact(5,i);
               0; 0; 0;
               yact(6,i)*yact(8,i) - yact(5,i)*yact(9,i) - g*yact(2,i);
               yact(4,i)*yact(9,i) - yact(6,i)*yact(7,i) + g*yact(1,i);
               g - yact(4,i)*yact(8,i) + yact(5,i)*yact(7,i);
               yact(7,i) + yact(2,i)*yact(9,i) - yact(3,i)*yact(8,i);
               yact(8,i) - yact(1,i)*yact(9,i) + yact(3,i)*yact(7,i);
               yact(9,i) + yact(1,i)*yact(8,i) - yact(2,i)*yact(7,i) ];
end

% Estimate derivatives (Version 2: First order forward difference)
% dy = diff(yact,1,2) / dt;
% res = dy - f(:,1:end-1);  % Residuals for parameter estimation
% y = yact(:,1:end-1);      % Aligned measurements
% u = uact(:,1:end-1);      % Aligned inputs

% Estimate derivatives (Version 1: Fourth order central difference)
for i = 3:length(yact)-3
    dy(:,i-2) = (1/(12*dt)) * (-yact(:,i+2) + 8*yact(:,i+1) ...
                               - 8*yact(:,i-1) + yact(:,i-2));
end
res = dy - f(:,3:end-3);  % Residuals for parameter estimation
y = yact(:,3:end-3);      % Aligned measurements
u = uact(:,3:end-3);      % Aligned inputs

% Construct library
Theta = [ y(5,:) .* y(6,:);
          u(2,:);
          y(4,:) .* y(6,:);
          u(3,:);
          y(4,:) .* y(5,:);
          u(4,:);
          u(1,:) ]';

% Define mask for active terms for parameter estimation
% Target states: x4, x5, x6, x9
mask = false(7,4);
mask(1:2,1) = true; 
mask(3:4,2) = true;
mask(5:6,3) = true;
mask(7,4) = true; 

% Perform sparse regression
lambda = 0;  % No thresholding (given masking constraint)
Xi = sparsifyDynamics(Theta, res([4,5,6,9],:)', lambda, mask);

% Extract estimated parameters
estimated_coeffs = [Xi(1:2,1)', Xi(3:4,2)', Xi(5:6,3)', Xi(7,4)];

% Display comparison
fprintf('True Coefficients: % .4f  % .4f  % .4f  % .4f  % .4f  % .4f  % .4f\n', true_coeffs);
fprintf('Estimated Coefficients: % .4f  % .4f  % .4f  % .4f  % .4f  % .4f  % .4f\n', estimated_coeffs);

% Custom SINDy with Structure Constraints
function Xi = sparsifyDynamics(Theta, dXdt, lambda, mask)
    % Copyright 2015, All Rights Reserved
    % Code by Steven L. Brunton
    % For Paper, "Discovering Governing Equations from Data: 
    %        Sparse Identification of Nonlinear Dynamical Systems"
    % by S. L. Brunton, J. L. Proctor, and J. N. Kutz

    % mask: [nTerms x nStates] logical matrix (true = parameter active)
    
    Xi = Theta \ dXdt;  % Initial least-squares guess
    
    for k = 1:10
        smallinds = abs(Xi) < lambda;
        Xi(smallinds) = 0;
        
        for i = 1:size(dXdt, 2)
            activeTerms = find(mask(:, i));
            if isempty(activeTerms)
                Xi(:, i) = 0;
            else
                Xi(:, i) = 0;
                Xi(activeTerms, i) = Theta(:, activeTerms) \ dXdt(:, i);
            end
        end
    end
end
