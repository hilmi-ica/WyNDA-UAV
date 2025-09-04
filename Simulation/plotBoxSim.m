function plotBoxSim(state)
    % PLOTQUADCOPTER Visualize the quadcopter as a gray cube with local axes
    
    % Extract position and orientation from state
    roll = state(1); 
    pitch = state(2);
    yaw = state(3); 
    pos = state(10:12);
    
    % Cube dimensions 
    cubeSize = 0.5;
    halfSize = cubeSize / 2;
    
    % Define cube vertices in the body frame
    vertices = [
        -halfSize, -halfSize, -halfSize;
         halfSize, -halfSize, -halfSize;
         halfSize,  halfSize, -halfSize;
        -halfSize,  halfSize, -halfSize;
        -halfSize, -halfSize,  halfSize;
         halfSize, -halfSize,  halfSize;
         halfSize,  halfSize,  halfSize;
        -halfSize,  halfSize,  halfSize
    ];
    
    % Cube faces
    faces = [
        1, 2, 3, 4; 
        5, 6, 7, 8; 
        1, 2, 6, 5;
        2, 3, 7, 6; 
        3, 4, 8, 7;
        4, 1, 5, 8  
    ];
    
    % Rotation matrix from roll, pitch, yaw
    Rz = [cos(yaw), -sin(yaw), 0; sin(yaw), cos(yaw), 0; 0, 0, 1];
    Ry = [cos(pitch), 0, sin(pitch); 0, 1, 0; -sin(pitch), 0, cos(pitch)];
    Rx = [1, 0, 0; 0, cos(roll), -sin(roll); 0, sin(roll), cos(roll)];
    Rfix = [0 1 0; -1 0 0; 0 0 1];
    R = Rfix * Rz * Ry * Rx; 

    % Rotate and translate cube vertices to the global frame
    verticesGlobal = (R * vertices')' + pos';
    
    % Plot cube
    patch('Vertices', verticesGlobal, 'Faces', faces, ...
          'FaceColor', [0.7, 0.7, 0.7], 'EdgeColor', 'k', 'FaceAlpha', 0.8);
    
    % Plot axes of the quadcopter
    axisLength = 2 * cubeSize; 
    origin = pos'; 
    xAxis = origin + R(:, 1)' * axisLength; 
    yAxis = origin + R(:, 2)' * axisLength; 
    zAxis = origin + R(:, 3)' * axisLength; 
    
    % Plot X-axis
    quiver3(origin(1), origin(2), origin(3), ...
            xAxis(1)-origin(1), xAxis(2)-origin(2), xAxis(3)-origin(3), ...
            0, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
    text(xAxis(1), xAxis(2), xAxis(3), 'X', 'Color', 'r', 'FontSize', 32);
    
    % Plot Y-axis
    quiver3(origin(1), origin(2), origin(3), ...
            yAxis(1)-origin(1), yAxis(2)-origin(2), yAxis(3)-origin(3), ...
            0, 'g', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
    text(yAxis(1), yAxis(2), yAxis(3), 'Y', 'Color', 'g', 'FontSize', 32);
    
    % Plot Z-axis
    quiver3(origin(1), origin(2), origin(3), ...
            zAxis(1)-origin(1), zAxis(2)-origin(2), zAxis(3)-origin(3), ...
            0, 'b', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
    text(zAxis(1), zAxis(2), zAxis(3), 'Z', 'Color', 'b', 'FontSize', 32);

end
