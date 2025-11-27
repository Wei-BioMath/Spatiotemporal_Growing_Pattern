addpath('/Volumes/T7 Shield/Quail_scMultiome/')
%% Growing-domain RD system (5 variables) with stretching advection (MOL)
clear; clc;

set(0, 'DefaultAxesFontSize', 18);             % For axis tick labels
set(0, 'DefaultAxesTitleFontSizeMultiplier', 1.2);  % Relative to axis font size
set(0, 'DefaultAxesLabelFontSizeMultiplier', 1.1);  % For xlabel/ylabel
set(0, 'DefaultLegendFontSize', 18);           % For legend text


rng(12345)
%% ---- Growth law L(t) and parameters
scale_factor = 10;
scale_factor_x = 2/5.49;
L0    = 5.49*scale_factor_x;                 % initial length
growR = 0.0155323;               % exponential growth rate
t0 =  50.64511;
Lfun  = @(t) L0*exp(growR*(t-t0));
dLdt  = @(t) growR*L0*exp(growR*(t-t0));

%% ---- Spatial grid on fixed computational domain xi in [0,1]
N   = 100;                    % grid points
xi  = linspace(0,1,N).';      % column
dxi = xi(2)-xi(1);

%% ---- Time span
tspan = linspace(0,100,10000).';

%% ---- Model parameters
Du = 1e-2*scale_factor_x^2;  Dv = 1e-1*3*1*scale_factor_x^2;  Dw = 5e-2*scale_factor_x^2*0.1;
De = 0.01*Dw*scale_factor_x^2; 
Ds = 1*0.01*2*0.01*Dw*500*scale_factor_x^2;
a = 0.1; b = 0.1; c = 0.4; e = 0.1; s = 0.1;

%% ---- Steady state (reaction-only) to seed ICs
F = @(X)[ ...
    a - X(1) + X(1)^2/X(2); ...
    b - X(2) + X(1)^2; ...
    c - X(3) + X(4)/(0.1+0.9*1.5*X(5)); ...
    e - X(4) + X(1); ...
    (s - X(5) + X(1))*0.001];
x0_guess = [1;1;1;1;1];
steady = fsolve(F, x0_guess, optimset('Display','off'));
U0 = steady(1); V0 = steady(2); W0 = steady(3); Y0 = steady(4); Z0 = 0.15; % your choice for Z

%% ---- Initial conditions on xi-grid (small noise + two bumps)
ic_bump = zeros(N,1);
mid_idx = N/2 + 1; % for even N
offset = round(N/4);
ic_bump(mid_idx - offset) = 0.1;
ic_bump(mid_idx + offset) = 0.1;

U_init = U0 + 0.01*randn(N,1) + ic_bump;
V_init = V0 + 0.01*randn(N,1);
W_init = W0 + 0.01*randn(N,1);
Y_init = Y0 + 0.01*randn(N,1);
Z_init = Z0 + 0.01*randn(N,1);

y0 = [U_init; V_init; W_init; Y_init; Z_init];   % stack into 5N vector

%% ---- ODE options
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

%% ---- Integrate with ode15s (MOL)
[tout, Y] = ode15s(@(t,y) mol_rhs_non_conserve(t,y,xi,dxi,Lfun,dLdt,Du,Dv,Dw,De,Ds,a,b,c,e,s,scale_factor), tspan, y0, opts);

%% ---- Reshape solution for convenience (time x space x var)
NT = numel(tout);
U = zeros(NT,N,5);
for k = 1:NT
    Yk = reshape(Y(k,:), [N,5]);  % columns: U,V,W,Y,Z
    U(k,:,:) = Yk;
end

%% ---- Build centered physical coordinates Xreal(t,xi) = (xi-0.5)*L(t)
[XI, TT] = meshgrid(xi, tout);
LTT = arrayfun(Lfun, TT);
Xreal_centered = (XI - 0.5) .* LTT;


n = 256;  % resolution
blueMap = [linspace(1,0,n)', linspace(1,0,n)', ones(n,1)];

lightBlue = [0.6 0.8 1];   % RGB for a pale blue
blueMap = [linspace(1,lightBlue(1),n)', ...
           linspace(1,lightBlue(2),n)', ...
           linspace(1,lightBlue(3),n)'];

mediumBlue = [0.3 0.5 1];   % brighter but not too dark
blueMap = [linspace(1,mediumBlue(1),n)', ...
           linspace(1,mediumBlue(2),n)', ...
           linspace(1,mediumBlue(3),n)'];

%% ---- Plot ASIP dynamics over space and time
figure; 
surf(Xreal_centered, TT-t0, squeeze(U(:,:,3))); 
colormap(blueMap); 
shading interp; view(2); % heatmap view
xlabel('x (centered)'); 
ylabel('t'); 
title('ASIP on growing domain');
colorbar;caxis([0 5.5]);
set(gca,'YDir','reverse')

ylim([-50 40]);   % desired t-range

% --- add dashed lines at specific time points ---
hold on
times_to_plot = [0, 24, 36];   % time points (already shifted by -t0 in Y axis)
xVals = get(gca,'XLim');       % current x-axis range
for tt = times_to_plot
    plot(xVals, [tt tt], 'k--', 'LineWidth', 1.5);  % dashed horizontal line
end
hold off

grid off;  % remove background grid
box on;    % keep axis box

set(gcf, 'Units', 'inches');fig_width = 8; fig_height = 7;
set(gcf, 'Position', [1, 1, fig_width, fig_height]);


%% ---- Plot ASIP snapshots at E4, E5, E5.5
t_values = [0, 24, 36] + t0;
t_indices = arrayfun(@(t) find(abs(tout - t) == min(abs(tout - t)), 1), t_values);

Lfun(t_values)/Lfun(t_values(1))
% --- precompute profiles and color limits
profiles2D = cell(1,length(t_indices));
x_coords = cell(1,length(t_indices));
cmin = inf; cmax = -inf;

for i = 1:length(t_indices)
    profile = squeeze(U(t_indices(i), :, 3));       % 1 x Nx
    profiles2D{i} = repmat(profile, 50, 1);        % 50 x Nx
    
    x_at_t = Xreal_centered(t_indices(i), :);      % 1 x Nx
    x_coords{i} = x_at_t;
    
    cmin = min(cmin, min(profile(:)));
    cmax = max(cmax, max(profile(:)));
end

% --- determine relative widths for each subplot
x_mins = cellfun(@min, x_coords);
x_maxs = cellfun(@max, x_coords);
rel_widths = (x_maxs - x_mins) / sum(x_maxs - x_mins);

% --- grid size (number of columns)
grid_cols = 30;
col_edges = round([0 cumsum(rel_widths)*grid_cols]);
labels = {'E4','E5','E5.5'};
fontsize = 18;

% --- create figure and axes
figure;
ax_all = gobjects(1,length(t_indices));  % store axes handles

for i = 1:length(t_indices)
    cols = col_edges(i)+1 : col_edges(i+1);   
    ax_all(i) = subplot(1, grid_cols, cols);  % single row grid
    
    imagesc(x_coords{i}, 1:50, profiles2D{i},[0 5.5]);
colormap(blueMap);
    set(ax_all(i),'YDir','normal');
     yticks(ax_all(i), []);
     xlabel('X','FontSize',fontsize);
    ylabel('Y','FontSize',fontsize);
    title(labels{i},'FontSize',fontsize);
    set(ax_all(i),'FontSize',fontsize);  
    if i==3
        colorbar;caxis([0 5.5]);  % unified color scale
    end 
end

% --- adjust each subplot to leave space for colorbar
for i = 1:length(ax_all)
    pos = get(ax_all(i),'Position');
    pos(3) = pos(3)*0.85;  % shrink width to fit colorbar
    set(ax_all(i),'Position',pos);
end

set(gcf, 'Units', 'inches');fig_width = 16; fig_height = 6;
set(gcf, 'Position', [1, 1, fig_width, fig_height]);
