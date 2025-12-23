clear;

%% ------------------------------------------------------------------------
%  Thermophysical properties of the PCM
%  All calculations are performed in SI units
% -------------------------------------------------------------------------

T_m   = -10;    % Phase change temperature [°C]
T_0   = -18;    % Boundary temperature at x = 0 [°C]
T_i   =  27;    % Initial uniform temperature [°C]

L     = 330e3;  % Latent heat [J/kg]

k_s   = 4.25;   % Thermal conductivity (solid) [W/mK]
k_l   = 0.60;   % Thermal conductivity (liquid) [W/mK]

rho_s = 1057;   % Density (solid) [kg/m^3]
rho_l = 1125;   % Density (liquid) [kg/m^3]

C_p_s = 1.9e3;  % Specific heat (solid) [J/kgK]
C_p_l = 3.4e3;  % Specific heat (liquid) [J/kgK]

% Thermal diffusivities
alfa_s = k_s / (rho_s * C_p_s);
alfa_l = k_l / (rho_l * C_p_l);

%% ------------------------------------------------------------------------
%  Lambda (λ) calculation for the Stefan problem
% -------------------------------------------------------------------------

func = @(lambda) ...
    (exp(-lambda^2) / erf(lambda)) + ...
    (k_l / k_s) * sqrt(alfa_s / alfa_l) * ...
    ((T_m - T_i) / (T_m - T_0)) * ...
    (exp(-lambda^2 * (alfa_s / alfa_l)) / ...
     erfc(lambda * sqrt(alfa_s / alfa_l))) - ...
    (lambda * L * sqrt(pi)) / (C_p_s * (T_m - T_0));

lambda = fzero(func, [0.0001, 5]);

%% ------------------------------------------------------------------------
%  Spatial and temporal discretization
% -------------------------------------------------------------------------

dx = 0.001;          % Spatial step [m]
dt = 1;              % Time step [s]

x  = 0:dx:1;         % Spatial domain [m]
t  = 1:dt:100000;    % Time vector [s]

T  = zeros(1, length(x));  % Temperature field

%% ------------------------------------------------------------------------
%  Import numerical (ANSYS Fluent) data
% -------------------------------------------------------------------------

data = importdata('1D_PCM\Numeric Data\1D_PCM-1000');
x_n  = data.data(:,2);

for k = 1:19
    data = importdata(['\1D_PCM\Numeric Data\1D_PCM-' num2str(k*1000)]);
    T_n(:,k) = data.data(:,4) - 273.15;   % Convert K → °C
end

%% ------------------------------------------------------------------------
%  Time marching and analytical temperature calculation
% -------------------------------------------------------------------------

counter = 0;

for ti = 1:length(t)

    % Solid–liquid interface position
    s_t = 2 * lambda * sqrt(alfa_s * t(ti));

    % Temperature distribution
    for xi = 1:length(x)
        if x(xi) < s_t
            T(xi) = ((T_m - T_0) * ...
                     erf(x(xi) / (2 * sqrt(alfa_s * t(ti))))) / ...
                     erf(lambda) + T_0;
        else
            T(xi) = ((T_m - T_i) * ...
                     erfc(x(xi) / (2 * sqrt(alfa_l * t(ti))))) / ...
                     erfc(lambda * sqrt(alfa_s / alfa_l)) + T_i;
        end
    end

    %% --------------------------------------------------------------------
    %  Plot results every 1000 seconds
    % --------------------------------------------------------------------

    if mod(t(ti), 1000) == 0
        counter = counter + 1;

        plot(x*100, T, '-r', 'LineWidth', 2.8);
        hold on

        idx = 1:5:length(x_n);
        plot(x_n(idx)*100, T_n(idx,counter), 'bo', ...
             'MarkerSize',5, 'MarkerFaceColor','none', 'LineWidth',1.2);

        xline(s_t*100, 'k-', 'LineWidth', 2);

        xlabel('Length [cm]');
        ylabel('Temperature [°C]');
        title('1D PCM Solidification: Analytical vs Numerical');

        legend('Analytical Solution', ...
               'Numerical Solution', ...
               'Solid–Liquid Interface', ...
               'Location','east');

        xlim([0 30]);

        % Time annotation (hours & minutes)
        time_hr  = floor(t(ti)/3600);
        time_min = floor(mod(t(ti),3600)/60);
        time_str = sprintf('t = %d h %02d min', time_hr, time_min);

        text(0.65, 0.70, time_str, ...
             'Units','normalized', ...
             'FontSize',12, ...
             'FontWeight','bold', ...
             'BackgroundColor','w', ...
             'EdgeColor','k');

        hold off
        pause(0.2)
    end
end
