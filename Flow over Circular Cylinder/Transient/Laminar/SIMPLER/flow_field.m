clear; 
clc; 
close all;

nx = 200;
ny = 100;
cx = 5.0;  
cy = 5.0;  
R  = 0.5;  

filename = 'cylinder_transient_SIMPLER.txt'; 
video_filename = 'cylinder_transient_SIMPLER.mp4';

if ~isfile(filename)
    error('File "%s" not found. Check the directory.', filename);
end

fprintf('Reading massive data file (this may take a moment)...\n');
rawText = fileread(filename);
blocks = split(rawText, 'ZONE');

if isempty(strtrim(blocks{1}))
    blocks(1) = [];
end

num_frames = length(blocks);
fprintf('Found %d time steps to animate.\n', num_frames);

v = VideoWriter(video_filename, 'MPEG-4');
v.FrameRate = 15; 
open(v);


fig = figure('Name', 'Cylinder Flow Animation', 'Position', [50, 100, 1600, 400]);


theta_cyl = linspace(0, 2*pi, 100);
xc_cyl = cx + R * cos(theta_cyl);
yc_cyl = cy + R * sin(theta_cyl);

fprintf('Generating animation...\n');

for k = 1:num_frames
    
    currentBlock = blocks{k};
    newlineIdx = find(currentBlock == newline, 1, 'first');
    if ~isempty(newlineIdx)
        currentBlock = currentBlock(newlineIdx:end);
    end
    
    data = str2num(currentBlock);
    
    if isempty(data)
        continue; 
    end
    
    X = reshape(data(:, 1), [nx, ny])';
    Y = reshape(data(:, 2), [nx, ny])';
    U = reshape(data(:, 3), [nx, ny])';
    V = reshape(data(:, 4), [nx, ny])';
    P = reshape(data(:, 5), [nx, ny])';
    
    VelMag = sqrt(U.^2 + V.^2);
    
    inside_cyl = (X - cx).^2 + (Y - cy).^2 <= R^2;
    VelMag(inside_cyl) = NaN;
    P(inside_cyl) = NaN;
    
    % Subplot 1: Velocity Magnitude
    subplot(1, 3, 1);
    contourf(X, Y, VelMag, 30, 'LineStyle', 'none');
    hold on;
    fill(xc_cyl, yc_cyl, [0.4 0.4 0.4], 'EdgeColor', 'k', 'LineWidth', 1.5); 
    hold off;
    colormap(gca, 'jet');
    if k == 1, cb1 = colorbar; cb1.Label.String = 'Velocity (m/s)'; end
    axis equal tight;
    title(sprintf('Velocity Magnitude | Frame: %d', k));
    
    % Subplot 2: Flow Streamlines
    subplot(1, 3, 2);
    cla; 

    streamslice(X, Y, U, V, 1.5); 
    hold on;
    fill(xc_cyl, yc_cyl, [0.4 0.4 0.4], 'EdgeColor', 'k', 'LineWidth', 1.5);
    axis equal tight;
    title('Flow Streamlines');
    
    % Subplot 3: Pressure Field
    subplot(1, 3, 3);
    contourf(X, Y, P, 30, 'LineStyle', 'none');
    hold on;
    fill(xc_cyl, yc_cyl, [0.4 0.4 0.4], 'EdgeColor', 'k', 'LineWidth', 1.5); 
    hold off;
    colormap(gca, 'jet');
    if k == 1, cb2 = colorbar; cb2.Label.String = 'Relative Pressure (Pa)'; end
    axis equal tight;
    title('Pressure Field');
    
    sgtitle('Transient Flow Over a Circular Cylinder (Re = 100)');
    
    drawnow;
    
    frame = getframe(fig);
    writeVideo(v, frame);
    
    if mod(k, 10) == 0
        fprintf('Processed frame %d / %d...\n', k, num_frames);
    end
end

close(v);