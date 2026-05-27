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

% only plot interior points
C_plot = C_plot(2:end-1, 2:end-1);
%% --- Create mesh ---
xc = load('x.dat');
yc = load('y.dat');
x = xc(2:end-1);
y = yc(2:end-1);

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