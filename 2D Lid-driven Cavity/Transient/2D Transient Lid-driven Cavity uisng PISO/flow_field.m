clear; 
clc; 
close all;

nx = 100;
ny = 100;

filename = 'cavity_transient_PISO.txt';
video_filename = 'Transient_Cavity_PISO.mp4';

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

fig = figure('Name', 'Transient Flow Animation', 'Position', [100, 100, 1400, 400]);
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
    
    % Subplot 1: Velocity Magnitude
    subplot(1, 3, 1);
    contourf(X, Y, VelMag, 30, 'LineStyle', 'none');
    hold on;
    plot([0 1 1 0 0], [0 0 1 1 0], '-k', 'LineWidth', 2);
    hold off;
    colormap(gca, 'jet');
    if k == 1, cb1 = colorbar; cb1.Label.String = 'Velocity (m/s)'; end
    axis equal tight;
    title(sprintf('Velocity Magnitude | Frame: %d', k));
    
    % Subplot 2: Flow Streamlines
    subplot(1, 3, 2);
    cla; 
    streamslice(X, Y, U, V, 2); 
    hold on;
    plot([0 1 1 0 0], [0 0 1 1 0], '-k', 'LineWidth', 2);
    axis equal tight;
    title('Flow Streamlines');
    
    % Subplot 3: Pressure Field
    subplot(1, 3, 3);
    contourf(X, Y, P, 30, 'LineStyle', 'none');
    hold on;
    plot([0 1 1 0 0], [0 0 1 1 0], '-k', 'LineWidth', 2);
    hold off;
    colormap(gca, 'jet');
    if k == 1, cb2 = colorbar; cb2.Label.String = 'Relative Pressure (Pa)'; end
    axis equal tight;
    title('Pressure Field');
    
    sgtitle('Transient PISO - Lid Driven Cavity');
    
    drawnow;
    
    frame = getframe(fig);
    writeVideo(v, frame);
    
    if mod(k, 10) == 0
        fprintf('Processed frame %d / %d...\n', k, num_frames);
    end
end

close(v);
fprintf('Animation complete! Video saved as: %s\n', video_filename);
