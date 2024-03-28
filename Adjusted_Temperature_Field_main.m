clear all;
clear variables;


%% define temperature data
addpath('Data\')
% addpath('functions')

Dat=load('TempField.mat');
X=double(Dat.phi);
Y=double(Dat.lambda);
T=Dat.T;
%% Polynomial function

% set degree of the polynomials p=1,2,3,4,5
p = 5;


% transformation of coordinates in the range of [0, 1]

    [x] = scaleCoord(X);               % function scaleCoord
    [y] = scaleCoord(Y);
    % to stable the values/grid;

% Linear Least Squares Adjustment
% compute design matrices in x and y
   Ax = evalPolynomial(x, p);                  % function evalPolynomial
   Ay = evalPolynomial(y, p);

l = T(:);

%% solve linear GMM

% combine matrices
A = kron(Ay, Ax);

[xS, lS] = linearGMM(A, l);                 % function linearGMM
L_mat=reshape(lS, size(T));

V_cap = lS - l; % residuals
Residuals=reshape(V_cap, size(T));
   

%% Harmonic Analysis 
% set degree of the Harmonic q = 30,40,50,60.....
  q=60;
  x_H = X;
  y_H = Y;

% Define the period T
Tx = 60;
Ty = 360;

% Define the frequency
% fx=@(t) sin(t*2*pi/Tx);
fx=2*pi/Tx;
% fy=@(t) cos(t*2*pi/Tx);
fy=2*pi/Ty;

%Too low frequency check
% fx=0.005;   % fx=2*pi/2*Tx; 
% fy=0.005;   % fy=2*pi/2*Ty;


%% Compute the cosine and sine terms
cos_ax = cos(bsxfun(@times, fx*x_H, 1:q));
sin_ax = sin(bsxfun(@times, fx*x_H, 1:q));
cos_ay = cos(bsxfun(@times, fy*y_H, 1:q));
sin_ay = sin(bsxfun(@times, fy*y_H, 1:q));


% Linear Least Squares Adjustment

    A_x = zeros(length(x_H), 2*q+1);
    A_y = zeros(length(y_H), 2*q+1);
    A_x(:, 1) = 1;
    A_y(:, 1) = 1;

j = 2;
for i = 1:q
    A_x(:, j:j+1) = [sin_ax(:, i), cos_ax(:, i)];
    A_y(:, j:j+1) = [sin_ay(:, i), cos_ay(:, i)];
    j = j + 2;
end


Nx = inv(A_x'*A_x);
Ny =  inv(A_y'*A_y);

X_tilda_x = Nx*A_x';
X_tilda_y = A_y*Ny;

x_S = X_tilda_x*T*X_tilda_y;
l_S = A_x*x_S*A_y';

V_tilda = l_S - T; % residuals
L_tilda = reshape(l_S, size(T));


%% create figure for Polynomial Function
figure(); 
% colored plot of the temperature data
subplot(1,2,1);
pcolor(Y, X, L_mat)
hold on; 
shading interp;
% This files contains data for plotting coastlines:
load('coast.mat');
plot(long, lat, 'k', 'LineWidth', 1.5);
set(gca, 'xlim', [-180 180], 'ylim', [0 90], 'DataAspectRatio', [1 1 1], ...
  'Xtick', -180:30:180,"YTick",0:20:90); 
% add labels and title
xlabel('longitude');
ylabel('latitude ');
title(['Adjusted Temperature field using polynomial p= ',num2str(p)]);
% add colorbar
colorbar
colormap('jet') % change colormap 
hold off 

subplot(1,2,2); 
% colored plot of the temperature data residual
pcolor(Y, X, Residuals)
hold on; 
shading interp;
% This files contains data for plotting coastlines:
load('coast.mat');
plot(long, lat, 'k', 'LineWidth', 1.5);
set(gca, 'xlim', [-180 180], 'ylim', [0 90], 'DataAspectRatio', [1 1 1], ...
  'Xtick', -180:30:180,"YTick",0:20:90); 
% add labels and title
xlabel('longitude');
ylabel('latitude ');
title(['Adjusted Temperature field of Residual using polynomial p= ',num2str(p)]);
% add colorbar
colorbar
colormap('jet') % change colormap 
hold off 


%% create fihure for Harmonic Analysis 

figure(); 
% colored plot of the temperature data
subplot(2,2,1);
pcolor(Y, X, L_tilda)
hold on; 
shading interp;
% This files contains data for plotting coastlines:
load('coast.mat');
plot(long, lat, 'k', 'LineWidth', 1.5);
set(gca, 'xlim', [-180 180], 'ylim', [0 90], 'DataAspectRatio', [1 1 1], ...
  'Xtick', -180:30:180,"YTick",0:20:90); 
% add labels and title
xlabel('longitude');
ylabel('latitude ');
title(['Adjusted Temperature field using Harmonic q= ',num2str(q)]);
% add colorbar
colorbar
colormap('jet') % change colormap 
hold off 

subplot(2,2,2); 
% colored plot of the temperature data residual
pcolor(Y, X, V_tilda)
hold on; 
shading interp;
% This files contains data for plotting coastlines:
load('coast.mat');
plot(long, lat, 'k', 'LineWidth', 1.5);
set(gca, 'xlim', [-180 180], 'ylim', [0 90], 'DataAspectRatio', [1 1 1], ...
  'Xtick', -180:30:180,"YTick",0:20:90); 
% add labels and title
xlabel('longitude');
ylabel('latitude ');
title(['Adjusted Temperature field of Residual using Harmonic q= ',num2str(q)]);
% add colorbar
colorbar
colormap('jet') % change colormap 
hold off 

