clear; 
clc;

nx = 200;
ny = 100;
cx = 5.0;
cy = 5.0;
R = 0.5;

disp('Loading data...');
data = load('SS_Cylinder_SIMPLE.txt');

X = reshape(data(:,1), [nx, ny])';
Y = reshape(data(:,2), [nx, ny])';
U = reshape(data(:,3), [nx, ny])';
V = reshape(data(:,4), [nx, ny])';
P = reshape(data(:,5), [nx, ny])';
Vel_Mag = sqrt(U.^2 + V.^2);


figure('Position', [100, 100, 1000, 800]);

% Subplot 1: Velocity Magnitude 
subplot(2,1,1);
contourf(X, Y, Vel_Mag, 50, 'LineColor', 'none');
colormap(gca, 'jet');
cb1 = colorbar;
cb1.Label.String = 'V (m/s)';
cb1.Label.FontSize = 12;

hold on;
theta = linspace(0, 2*pi, 100);
fill(cx + R*cos(theta), cy + R*sin(theta), 'k', 'EdgeColor', 'none'); 
hold off;

title('Velocity Field', 'FontSize', 14);
xlabel('x (m)', 'FontSize', 12); 
ylabel('y (m)', 'FontSize', 12);
axis equal tight;

% Subplot 2: Pressure Field 
subplot(2,1,2);
contourf(X, Y, P, 50, 'LineColor', 'none');
colormap(gca, 'jet');
cb2 = colorbar;
cb2.Label.String = 'Static Pressure (Pa)';
cb2.Label.FontSize = 12;

hold on;
fill(cx + R*cos(theta), cy + R*sin(theta), 'k', 'EdgeColor', 'none'); 
hold off;

title('Pressure Field', 'FontSize', 14);
xlabel('x (m)', 'FontSize', 12); 
ylabel('y (m)', 'FontSize', 12);
axis equal tight;