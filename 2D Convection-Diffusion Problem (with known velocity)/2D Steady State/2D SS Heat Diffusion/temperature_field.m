clear; 
clc;

% Load the file
data = load('flow_field.txt');

% Extract the columns (x, y, phi)
x_1D = data(:, 1);
y_1D = data(:, 2);
phi_1D = data(:, 3);

% Define grid dimensions 
nx = 100;
ny = 100;

% Reshape the 1D arrays into 2D matrices
X = reshape(x_1D, [ny, nx]);
Y = reshape(y_1D, [ny, nx]);
PHI = reshape(phi_1D, [ny, nx]);

% Plot the filled contour
figure;
contourf(X, Y, PHI, 100, 'LineColor', 'none'); 
cb = colorbar;
ylabel(cb, 'Temperature T')
colormap('jet');

% Formatting the plot
title('2D SS Heat Diffusion without Source', 'FontSize', 14);
xlabel('x', 'FontSize', 12);
ylabel('y', 'FontSize', 12);
axis equal; 
axis tight;