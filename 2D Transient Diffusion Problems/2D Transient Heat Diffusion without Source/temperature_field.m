clear; 
clc;

nx = 300;
ny = 300;
Lx = 1.0;
Ly = 1.0;
output_freq = 50;

disp('Loading data...');
data = load('T_out.txt');
num_frames = size(data, 1) / ny;
disp(['Found ', num2str(num_frames), ' frames. Generating video...']);

% Initialization
video_filename = 'Transient_Heat_Diffusion.mp4';
v = VideoWriter(video_filename, 'MPEG-4');
v.FrameRate = 15; 
open(v);

dx = Lx / nx;
dy = Ly / ny;
x = (dx/2) : dx : (Lx - dx/2);
y = (dy/2) : dy : (Ly - dy/2);
[X, Y] = meshgrid(x, y);


fig = figure('Position', [100, 100, 800, 600]);

for k = 1:num_frames
    
    row_start = (k - 1) * ny + 1;
    row_end = k * ny;
    phi_frame = data(row_start:row_end, :);

    contourf(X, Y, phi_frame, 100, 'LineColor', 'none');
    
    colormap('jet'); 
    clim([0 100]); 
    
    cb = colorbar;
    cb.Label.String = 'Temperature';
    cb.Label.FontSize = 14;
    
    if k == 1
        current_step = 1;
    else
        current_step = (k-1) * output_freq;
    end
    title(sprintf('Transient Heat Diffusion - Time Step: %d', current_step), 'FontSize', 14);
    xlabel('x', 'FontSize', 12);
    ylabel('y', 'FontSize', 12);
    axis equal tight;
    
    drawnow;
    frame = getframe(fig);
    writeVideo(v, frame);
end

close(v);
