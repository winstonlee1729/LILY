clear; clc;
data = load('cavity_SIMPLE_omp.txt');

% Change grid size when they were altered in Fortran solver
nx = 100; ny = 100;

X = reshape(data(:, 1), nx, ny)';
Y = reshape(data(:, 2), nx, ny)';
U = reshape(data(:, 3), nx, ny)';
V = reshape(data(:, 4), nx, ny)';

% Calculate velocity magnitude: |V| = sqrt(u^2 + v^2)
VelMag = sqrt(U.^2 + V.^2);

% Plot 1: Velocity magnitude field 
figure('Color', 'w', 'Position', [100 100 800 600]);
contourf(X, Y, VelMag, 50, 'LineStyle', 'none');
colormap(jet);
cb = colorbar;
ylabel(cb, 'Velocity Magnitude (m/s)');
hold on;

% Add vectors for direction
skip = 2;
quiver(X(1:skip:end, 1:skip:end), Y(1:skip:end, 1:skip:end), ...
       U(1:skip:end, 1:skip:end), V(1:skip:end, 1:skip:end), 1.2, 'w');

title('Steady-State Velocity Profile (Lid-Driven Cavity)');
xlabel('x (m)'); ylabel('y (m)');
axis equal tight;

% Plot 2: Benchmark mid-line profiles
figure('Color', 'w', 'Position', [950 100 500 600]);

% U-velocity along vertical mid-line (x = 0.5)
subplot(2,1,1);
mid_x = round(nx/2);
plot(U(:, mid_x), Y(:, mid_x), 'b-o', 'LineWidth', 1.5);
grid on;
title('u-velocity along vertical mid-line (x=0.5)');
xlabel('u (m/s)'); ylabel('y (m)');

% V-velocity along horizontal mid-line (y = 0.5)
subplot(2,1,2);
mid_y = round(ny/2);
plot(X(mid_y, :), V(mid_y, :), 'r-s', 'LineWidth', 1.5);
grid on;
title('v-velocity along horizontal mid-line (y=0.5)');
xlabel('x (m)'); ylabel('v (m/s)');