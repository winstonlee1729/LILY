clear; 
clc;

nx = 300;
ny = 300;
Lx = 1.0;
Ly = 1.0;
output_freq = 25;

disp('Loading data from phi_out.txt... (This may take a moment)');
data = load('phi_out.txt');

% Calculate the total number of frames in the text file
num_frames = size(data, 1) / ny;
disp(['Found ', num2str(num_frames), ' frames. Generating video...']);

% Initialize the Video Writer
video_filename = 'phi_transient_video.mp4';
v = VideoWriter(video_filename, 'MPEG-4');
v.FrameRate = 10; 
open(v);

% Create grid coordinates for the cell centers
dx = Lx / nx;
dy = Ly / ny;
x = (dx/2) : dx : (Lx - dx/2);
y = (dy/2) : dy : (Ly - dy/2);
[X, Y] = meshgrid(x, y);

% Setup figure window
fig = figure;  

% Start loop
for k = 1:num_frames
    
    % Extract the specific 300x300 block for the current frame
    row_start = (k - 1) * ny + 1;
    row_end = k * ny;
    phi_frame = data(row_start:row_end, :);

    % Plot the filled contour
    contourf(X, Y, phi_frame, 100, 'LineColor', 'none');
    
    % Apply colormap
    colormap('jet');
 
    clim([0 100]); 
    cb = colorbar;
    cb.Label.String = '\phi';
    cb.Label.FontSize = 14;
    cb.Label.Rotation = 0; 
    
    % Dynamic Title updating with frame number
    if k == 1
        current_step = 1;
    else
        current_step = (k-1) * output_freq;
    end
    title(sprintf('\\phi Flow Field - Time Step: %d', current_step), 'FontSize', 14);
    xlabel('x', 'FontSize', 12);
    ylabel('y', 'FontSize', 12);
    axis equal tight;
    
    drawnow;
    frame = getframe(fig);
    writeVideo(v, frame);
end

% Close the video file
close(v);