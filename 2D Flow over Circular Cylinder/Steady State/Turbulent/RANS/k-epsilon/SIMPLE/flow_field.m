clear; 
clc;

nx = 200;
ny = 100;
cx = 5.0;
cy = 5.0;
R = 0.5;

disp('Loading data...');
data = load('SS_Cylinder_RANS.txt');

X   = reshape(data(:,1), [nx, ny])';
Y   = reshape(data(:,2), [nx, ny])';
U   = reshape(data(:,3), [nx, ny])';
V   = reshape(data(:,4), [nx, ny])';
P   = reshape(data(:,5), [nx, ny])';
Mut = reshape(data(:,6), [nx, ny])'; 

Vel_Mag = sqrt(U.^2 + V.^2);

theta_cyl = linspace(0, 2*pi, 100);
x_cyl = cx + R*cos(theta_cyl);
y_cyl = cy + R*sin(theta_cyl);

figure('Position', [100, 100, 1000, 800]);
% Subplot 1: Velocity Magnitude 
subplot(2,1,1);
contourf(X, Y, Vel_Mag, 100, 'LineColor', 'none');
colormap(gca, 'jet');
cb1 = colorbar;
cb1.Label.String = 'V (m/s)';
cb1.Label.FontSize = 12;

hold on;
fill(x_cyl, y_cyl, 'k', 'EdgeColor', 'none'); % Draw cylinder
hold off;

title('RANS k-\epsilon: Steady State Velocity Field', 'FontSize', 14);
xlabel('x (m)', 'FontSize', 12); 
ylabel('y (m)', 'FontSize', 12);
axis equal tight;

% Subplot 2: Eddy Viscosity (Turbulence) 
subplot(2,1,2);
contourf(X, Y, Mut, 100, 'LineColor', 'none');
colormap(gca, 'hot'); 
cb2 = colorbar;
cb2.Label.String = '\mu_t (Turbulent Eddy Viscosity)';
cb2.Label.FontSize = 12;

hold on;
fill(x_cyl, y_cyl, 'w', 'EdgeColor', 'k'); % Draw cylinder in white
hold off;

title('RANS k-\epsilon: Turbulent Eddy Viscosity Field', 'FontSize', 14);
xlabel('x (m)', 'FontSize', 12); 
ylabel('y (m)', 'FontSize', 12);
axis equal tight;