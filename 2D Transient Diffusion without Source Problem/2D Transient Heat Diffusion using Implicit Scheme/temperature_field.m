clear; 
clc; 
close all;

% Parameters 
nx = 300;
ny = 300;
Lx = 1.0;
Ly = 1.0;
dt = 0.05;

% Output frequency from Fortran code 
fortran_output_freq = 50;
total_steps = 2000; 

% Load the raw data from Fortran
disp('Loading data from T_out.txt...');
if exist('T_out.txt', 'file') ~= 2
    error('File T_out.txt not found. Make sure it is in the same folder.');
end
rawData = load('T_out.txt');

% Calculate how many frames are in the file
total_rows = size(rawData, 1);
num_frames = total_rows / ny;

% Calculate target frames 
target_steps = 50:50:total_steps;
target_frames = target_steps / fortran_output_freq + 1;
% Ensure target_frames are integers
target_frames = round(target_frames);


% Set up single figure and axes once before the loop for reuse
f = figure('Name', '2D Unsteady Heat Diffusion', 'Position', [100, 100, 600, 500]);
ha = axes(f);
x_coords = linspace(0, Lx, nx);
y_coords = linspace(0, Ly, ny);

% Loop through each target frame and save image
for k = 1:length(target_frames)
    current_frame = target_frames(k);
    current_step = target_steps(k);
    
    % Calculate the rows we need to extract
    start_row = (current_frame - 1) * ny + 1;
    end_row = current_frame * ny;
    
    % If the requested rows don't exist in the file, stop gracefully.
    if end_row > size(rawData, 1)
        fprintf('\nWarning: Data for step %d (frame %d) was not found in T_out.txt.\n', current_step, current_frame);
        disp('Stopping image generation based on available data.');
        break; % Exits the loop safely
    end

    fprintf('Processing step %d (frame %d)...\n', current_step, current_frame);

    % Extract the specific ny-by-nx block for this frame
    T_frame = rawData(start_row:end_row, :);

    % Plot frame on single axes explicitly to update existing plot
    imagesc(ha, x_coords, y_coords, T_frame);
    axis(ha, 'xy');              
    colormap(ha, 'jet');         
    colorbar(ha);
    caxis(ha, [0 100]);          

    % Set title and labels
    title(ha, sprintf('Time Step: %d | Time: %.2f s', current_step, current_step*dt));
    xlabel(ha, 'X (m)');
    ylabel(ha, 'Y (m)');

    % Save figure
    image_name = sprintf('T_step_%04d.png', current_step);
    saveas(f, image_name);
    fprintf('Saved image %s\n', image_name);
end