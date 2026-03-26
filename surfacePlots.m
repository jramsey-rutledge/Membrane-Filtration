clear; clc; close all;

%% --- USER INPUT ---
Nx = 64;          % must match simulation
Ny = 64;
Lx = 0.05;            % change per case
Ly = 0.002;

%% --- Load data ---
C_full = load(['C_',num2str(Nx),'.dat']);

Nx2 = Nx + 2;
Ny2 = Ny + 2;

% reshape (C++ is column-major flattened as i + j*Nx2)
C_full = reshape(C_full, [Nx2, Ny2])';

% update boundaries
C_plot = C_full;
% C_plot(1, :) = 0.5 * (C_full(1, :) + C_full(2, :));
% C_plot(end, :) = 0.5 * (C_full(end, :) + C_full(end-1, :));
% C_plot(:, 1) = 0.5 * (C_full(:, 1) + C_full(:, 2));
% C_plot(:, end) = 0.5 * (C_full(:, end) + C_full(:, end-1));
% 
% % update corners
% C_plot(1,1) = 2*C_plot(2,2) - C_plot(3,3);
% C_plot(1,end) = 2*C_plot(2,end-1) - C_plot(3,end-2);
% C_plot(end,1) = 2*C_plot(end-1,2) - C_plot(end-2,3);
% C_plot(end,end) = 2*C_plot(end-1,end-1) - C_plot(end-2,end-2);

% only plot interior points
C_plot = C_plot(2:end-1, 2:end-1);
%% --- Create mesh ---
%x_full = linspace(0, Lx, Nx2);
%y_full = linspace(0, Ly, Ny2);
%[X, Y] = meshgrid(x_full, y_full);

dx = Lx/Nx; 
dy = Ly/Ny; 
x = 0.5*dx + (0:Nx-1)*dx; 
y = 0.5*dy + (0:Ny-1)*dy; 

[X,Y] = meshgrid(x,y);

%% --- Create Concentration Plot --
C = C_plot;
%% 3-D Plot
figure;
surf(X, Y, C);   
shading interp;
colorbar;

xlabel('x');
ylabel('y');
zlabel('C');

title('Concentration Field C(x,y)');

%% 2-D Boundary Plot (C vs x at y ≈ 0)
figure;

C_bottom = C(1,:);   

plot(X(1,:), C_bottom, 'b-', 'LineWidth', 2);
xlabel('x');
ylabel('C(x,0)');
title('Boundary Profile at y ≈ 0');
ylim([1 2]);

grid on;