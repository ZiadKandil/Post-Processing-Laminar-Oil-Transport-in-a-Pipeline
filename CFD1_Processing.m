close all
clear all
clc

% set(0, 'DefaultFigureWindowStyle', 'docked')

% Input Data

D = 0.15;
Wb = 0.45;
rho = 910.0;
nu = 3.5e-4;
L = 3.0;
mu = nu * rho;
R = D / 2;

% Calculate Reynold's number

Re = D * Wb / nu;
fprintf('Reynolds Number = %.2f\n', Re)

% Claculate entrance length based on 2 formulas
Le1 = 0.05 * Re * D;
Le2 = 0.057 * Re * D;

fprintf('Entrance Length based on first formula = %.2f\n', Le1)
fprintf('Entrance Length based on second formula = %.2f\n', Le2)

% import mesh and mesh data

mesh = load('main.mat');
NZ = mesh.NZ;
NY = mesh.NY;
Y = mesh.Y_C;
dz = L / NZ;
z = dz / 2 : dz : L - dz / 2;       % Cells centers
W = mesh.W1;
W = squeeze(W);
P = mesh.P1;
P = squeeze(P);
DWDY = mesh.DWDY;
DWDY = squeeze(DWDY);
tau_wall = - mu * DWDY(30, :);

%  Velocity and tau_wall at 0.5m, 1.0m, 1.45m, 1.65m, 2.0m to verify fully developed flow

[~, i] = min(abs(z - 0.5));
W_500 = W(:, i);
tau_wall_500 = tau_wall(i);

[~, i] = min(abs(z - 1.0));
W_1000 = W(:, i);
tau_wall_1000 = tau_wall(i);

[~, i] = min(abs(z - 1.45));
W_1450 = W(:, i);
tau_wall_1450 = tau_wall(i);

[~, i] = min(abs(z - 1.65));
W_1650 = W(:, i);
tau_wall_1650 = tau_wall(i);

[~, i] = min(abs(z - 2.0));
W_2000 = W(:, i);
tau_wall_2000 = tau_wall(i);


figure('Color','w'); hold on; box on
colors = lines(5);

plot(W_500, Y, '-','LineWidth', 2.5, 'Color', colors(1,:), 'Marker','none')
plot(W_1000, Y, '--d','LineWidth', 2.5, 'Color', colors(2,:), 'MarkerSize',5, 'MarkerFaceColor', colors(2,:))
plot(W_1450, Y, '--o','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(3,:))
plot(W_1650, Y, '--s','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(4,:))
plot(W_2000, Y, '--','LineWidth', 2.5, 'Color', colors(5,:), 'MarkerSize',5, 'MarkerFaceColor', colors(5,:))

xlabel('Velocity W (m/s)','FontSize',12,'FontWeight','bold')
ylabel('Pipe Radius R (m)','FontSize',12,'FontWeight','bold')
legend({'0.5 m','1.0 m','1.45 m','1.65 m', '2.0 m'}, 'Location','best', 'FontSize',11)

figure('Color','w'); hold on; box on
colors = lines(5);

plot(0.5,  tau_wall_500,  'd','LineWidth',2.5,'MarkerSize',6,...
     'Color',colors(1,:),'MarkerFaceColor',colors(1,:))
plot(1.0,  tau_wall_1000, 'o','LineWidth',2.5,'MarkerSize',6,...
     'Color',colors(2,:),'MarkerFaceColor',colors(2,:))
plot(1.45, tau_wall_1450, 's','LineWidth',2.5,'MarkerSize',6,...
     'Color',colors(3,:),'MarkerFaceColor',colors(3,:))
plot(1.65, tau_wall_1650, '^','LineWidth',2.5,'MarkerSize',6,...
     'Color',colors(4,:),'MarkerFaceColor',colors(4,:))
plot(2.0,  tau_wall_2000, 'v','LineWidth',2.5,'MarkerSize',6,...
     'Color',colors(5,:),'MarkerFaceColor',colors(5,:))

xlim([0 3.0])
ylim([6 10])
grid on

xlabel('Axial Distance (m)','FontSize',12,'FontWeight','bold')
ylabel('Wall Shear Stress (Pa)','FontSize',12,'FontWeight','bold')
legend({'0.5 m','1.0 m','1.45 m','1.65 m', '2.0 m'}, 'Location','best', 'FontSize',11)


% Grid Independence study

mesh1 = load('mesh1.mat');
mesh2 = load('mesh2.mat');
mesh3 = load('mesh3.mat');

W1 = mesh1.W1;
W1 = squeeze(W1);
P1 = mesh1.P1;
P1 = squeeze(P1);
DWDY_1 = mesh1.DWDY;
DWDY_1 = squeeze(DWDY_1);

W2 = mesh2.W1;
W2 = squeeze(W2);
P2 = mesh2.P1;
P2 = squeeze(P2);
DWDY_2 = mesh2.DWDY;
DWDY_2 = squeeze(DWDY_2);

W3 = mesh3.W1;
W3 = squeeze(W3);
P3 = mesh3.P1;
P3 = squeeze(P3);
DWDY_3 = mesh3.DWDY;
DWDY_3 = squeeze(DWDY_3);

% Extracting values of velocity and pressure at the center of the pipe and
% 2.0 meters from inlet

[~, i] = min(abs(z - 2.0));
Wloc = W(1, i);
ploc = P(1, i);

NZ1 = mesh1.NZ;
NY1 = mesh1.NY;
Y1 = mesh1.Y_C;
dz1 = L / NZ1;
z1 = dz1 / 2 : dz1 : L - dz1 / 2; 

[~, i] = min(abs(z1 - 2.0));
Wloc1 = W1(1, i);
ploc1 = P1(1, i);
W1_2000 = W1(:, i);
tau_wall1 = - mu * DWDY_1(NY1, :);
tau_wall1_2000 = tau_wall1(i);

NZ2 = mesh2.NZ;
NY2 = mesh2.NY;
Y2 = mesh2.Y_C;
dz2 = L / NZ2;
z2 = dz2 / 2 : dz2 : L - dz2 / 2; 

[~, i] = min(abs(z2 - 2.0));
Wloc2 = W2(1, i);
ploc2 = P2(1, i);
W2_2000 = W2(:, i);
tau_wall2 = - mu * DWDY_2(NY2, :);
tau_wall2_2000 = tau_wall2(i);

NZ3 = mesh3.NZ;
NY3 = mesh3.NY;
Y3 = mesh3.Y_C;
dz3 = L / NZ3;
z3 = dz3 / 2 : dz3 : L - dz3 / 2; 

[~, i] = min(abs(z3 - 2.0));
Wloc3 = W3(1, i);
ploc3 = P3(1, i);
W3_2000 = W3(:, i);
tau_wall3 = - mu * DWDY_3(NY3, :);
tau_wall3_2000 = tau_wall3(i);


figure('Color','w');
tiledlayout(1,3,'Padding','compact','TileSpacing','compact')

nexttile
hold on; box on
colors = lines(4);

plot(NZ1, ploc1, 'o', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', colors(1,:), ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2)

plot(NZ2, ploc2, 'x', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'LineWidth', 2, ...
    'Color', colors(2,:))

plot(NZ, ploc, 's', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', colors(3,:), ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2)

plot(NZ3, ploc3, 'd', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', colors(4,:), ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2)

xlim([50 500])
ylim([20 500])
xlabel('Number of Elements', 'FontSize', 12, 'FontWeight','bold')
ylabel('Pressure (Pa)', 'FontSize', 12, 'FontWeight','bold')
legend({'Mesh 1','Mesh 2','Mesh 3','Mesh 4'}, 'Location','best', 'FontSize',11)
grid on

nexttile
hold on; box on

plot(NZ1, Wloc1, 'o', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', colors(1,:), ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2)

plot(NZ2, Wloc2, 'x', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'LineWidth', 2, ...
    'Color', colors(2,:))

plot(NZ, Wloc, 's', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', colors(3,:), ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2)

plot(NZ3, Wloc3, 'd', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', colors(4,:), ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2)

xlim([50 500])
ylim([0.75 1.1])
xlabel('Number of Elements', 'FontSize', 12, 'FontWeight','bold')
ylabel('Velocity (m/s)', 'FontSize', 12, 'FontWeight','bold')
legend({'Mesh 1','Mesh 2','Mesh 3','Mesh 4'}, 'Location','best', 'FontSize',11)
grid on

nexttile
hold on; box on

plot(NZ1, tau_wall1_2000, 'o', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', colors(1,:), ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2)

plot(NZ2, tau_wall2_2000, 'x', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'LineWidth', 2, ...
    'Color', colors(2,:))

plot(NZ, tau_wall_2000, 's', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', colors(3,:), ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2)

plot(NZ3, tau_wall3_2000, 'd', ...
    'LineStyle','none', ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', colors(4,:), ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.2)

xlim([50 500])
ylim([6 9])
xlabel('Number of Elements', 'FontSize', 12, 'FontWeight','bold')
ylabel('Wall Shear Stress (Pa)', 'FontSize', 12, 'FontWeight','bold')
legend({'Mesh 1','Mesh 2','Mesh 3','Mesh 4'}, 'Location','best', 'FontSize',11)
grid on


% Velocity Profile at 2.0 meters from inlet

figure('Color','w'); hold on; box on
colors = lines(4);

plot(W1_2000, Y1, '-','LineWidth', 2.5, 'Color', colors(1,:), 'Marker','none')
plot(W2_2000, Y2, '--d','LineWidth', 2.5, 'Color', colors(2,:), 'MarkerSize',5, 'MarkerFaceColor', colors(2,:))
plot(W_2000, Y, '--o','LineWidth', 2.5, 'Color', colors(3,:), 'MarkerSize',5, 'MarkerFaceColor', colors(3,:))
plot(W3_2000, Y3, '--s','LineWidth', 2.5, 'Color', colors(4,:), 'MarkerSize',5, 'MarkerFaceColor', colors(4,:))

xlabel('Velocity W (m/s)','FontSize',12,'FontWeight','bold')
ylabel('Pipe Radius R (m)','FontSize',12,'FontWeight','bold')
legend({'Mesh 1','Mesh 2','Mesh 3','Mesh 4'}, 'Location','best', 'FontSize',11)


% Validation against analytical solution

% Simulation Values at 2.0m

[~, i] = min(abs(z - 2.0));
tau = - mu * DWDY(:,i);
[~, i] = min(abs(z - 2.0));
p_2000 = P(1, i);
[~, i] = min(abs(z - 2.5));
p_2500 = P(1, i);
DPDZ = (p_2500 - p_2000) / 0.5;

% Analytical values

W_analytical = 2 * Wb * (1 - (Y.^2) ./ (R^2))';
DPDZ_analytical = - 8 * mu * Wb / (R^2);
tau_analytical = (4 * mu * Wb / (R^2)) * Y';

fprintf('CFD Pressure Gradient = %.2f\n', DPDZ)
fprintf('Analytical Pressure Gradient = %.2f\n', DPDZ_analytical)


figure('Color','w'); hold on; box on
colors = lines(2);

plot(W_2000, Y, '-','LineWidth', 2.5, 'Color', colors(1,:), 'Marker','none')
plot(W_analytical, Y, '--d','LineWidth', 2.5, 'Color', colors(2,:), 'MarkerSize',5, 'MarkerFaceColor', colors(2,:))

xlabel('Velocity W (m/s)','FontSize',12,'FontWeight','bold')
ylabel('Pipe Radius R (m)','FontSize',12,'FontWeight','bold')
legend({'CFD','Analytical'}, 'Location','best', 'FontSize',11)

figure('Color','w'); hold on; box on
colors = lines(2);

plot(tau, Y, '-','LineWidth', 2.5, 'Color', colors(1,:), 'Marker','none')
plot(tau_analytical, Y, '--d','LineWidth', 2.5, 'Color', colors(2,:), 'MarkerSize',5, 'MarkerFaceColor', colors(2,:))

xlabel('Shear Stress (Pa)','FontSize',12,'FontWeight','bold')
ylabel('Pipe Radius R (m)','FontSize',12,'FontWeight','bold')
legend({'CFD','Analytical'}, 'Location','best', 'FontSize',11)

figure('Color','w'); hold on; box on
colors = lines(1);

plot(z, P(1,:), '-','LineWidth', 2.5, 'Color', colors(1,:), 'Marker','none')

xlabel('Axial Distance (m)','FontSize',12,'FontWeight','bold')
ylabel('Pressure (Pa)','FontSize',12,'FontWeight','bold')

% Friction Factor

% CFD Results
f = (8 * tau(30)) / (Wb^2 * rho);

% Analytical
f_analytical = 64 / Re;

fprintf('CFD Friction Factor = %.2f\n', f)
fprintf('Analytical Friction Factor = %.2f\n', f_analytical)

% Different Reynolds numbers

bulk_vel = [0.3, 1.0, 2.0, 3.0];
bulk_Re = bulk_vel * D / nu;

% Analytical friction factors

f_analytical = 64 ./ bulk_Re;

% CFD Results

Re_03 = load('WB_03.mat');
Re_1 = load('WB_1.mat');
Re_2 = load('WB_2.mat');
Re_3 = load('WB_3.mat');

DWDY03 = squeeze(Re_03.DWDY);
DWDY1 = squeeze(Re_1.DWDY);
DWDY2 = squeeze(Re_2.DWDY);
DWDY3 = squeeze(Re_3.DWDY);

[~, i] = min(abs(z - 2.0));

tau03 = - mu * DWDY03(:, i);
tau1 = - mu * DWDY1(:, i);
tau2 = - mu * DWDY2(:, i);
tau3 = - mu * DWDY3(:, i);

f03 = (8 * tau03(30)) / ((0.3)^2 * rho);
f1 = (8 * tau1(30)) / ((1)^2 * rho);
f2 = (8 * tau2(30)) / ((2)^2 * rho);
f3 = (8 * tau3(30)) / ((3)^2 * rho);

f_cfd = [f03, f1, f2, f3];

figure('Color','w'); hold on; box on
colors = lines(2);

plot(bulk_Re, f_cfd, '-','LineWidth', 2.5, 'Color', colors(1,:), 'Marker','none')
plot(bulk_Re, f_analytical, '--d','LineWidth', 2.5, 'Color', colors(2,:), 'MarkerSize',5, 'MarkerFaceColor', colors(2,:))

xlabel('Reynold Number Re','FontSize',12,'FontWeight','bold')
ylabel('friction Factor','FontSize',12,'FontWeight','bold')
legend({'CFD','Analytical'}, 'Location','best', 'FontSize',11)