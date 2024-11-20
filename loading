% Advanced Wind and Ice Loading Analysis for Transmission Lines with Bar Plots
clc; clear;

% Fixed Parameters
span = 200;               % Span between towers (m)
diameter = 0.03;          % Diameter of the cable (m)
line_density = 2700;      % Line material density (kg/m^3)
g = 9.81;                 % Gravity (m/s^2)
rho_air = 1.225;          % Air density (kg/m^3)
drag_coefficient = 1.2;   % Drag coefficient for wind
elastic_modulus = 2e11;   % Elastic modulus of cable material (Pa)

% Variable Parameters
ice_thickness_values = linspace(0, 0.03, 5);  % Ice thickness (m)
wind_speed_values = linspace(10, 50, 5);      % Wind speed (m/s)

% Calculations
x = linspace(0, span, 100); % Horizontal positions along the span (m)

% Preallocate for plotting
tension_curves = zeros(length(ice_thickness_values), length(x));
sag_curves = zeros(length(wind_speed_values), length(x));
max_tension = zeros(size(ice_thickness_values)); % For bar plot
max_sag = zeros(size(wind_speed_values));        % For bar plot

% Constants
line_area = pi * (diameter / 2)^2;  % Cross-sectional area of cable (m^2)
line_weight = line_density * g * line_area; % Line weight per unit length (N/m)

% Tension Analysis: Varying Ice Thickness
for i = 1:length(ice_thickness_values)
    ice_thickness = ice_thickness_values(i);
    ice_weight = pi * diameter * ice_thickness * 900 * g; % Ice weight per unit length (N/m)
    effective_diameter = diameter + 2 * ice_thickness;    % Effective cable diameter (m)
    
    % Total weight per unit length (N/m)
    total_weight = line_weight + ice_weight;
    
    for j = 1:length(x)
        % Formula for vertical tension (due to wind and ice loading)
        vertical_tension = total_weight * abs(x(j) - span / 2);
        
        % Horizontal tension (constant due to line weight)
        horizontal_tension = sqrt(vertical_tension^2 + (tension_curves(1, 1))^2);
        
        % Storing tension values
        tension_curves(i, j) = horizontal_tension;
    end
    % Store maximum tension for bar plot
    max_tension(i) = max(tension_curves(i, :));
end

% Sag Analysis: Varying Wind Speed
for i = 1:length(wind_speed_values)
    wind_speed = wind_speed_values(i);
    wind_force = 0.5 * rho_air * wind_speed^2 * drag_coefficient * diameter; % Wind force (N/m)
    total_load = sqrt(line_weight^2 + wind_force^2); % Combined load (N/m)
    
    % Formula for sag calculation (using catenary approximation)
    for j = 1:length(x)
        sag_curves(i, j) = (total_load / (2 * tension_curves(1, 1))) * x(j) * (span - x(j));
    end
    % Store maximum sag for bar plot
    max_sag(i) = max(sag_curves(i, :));
end

% Plot Tension Curves (Varying Ice Thickness)
figure;
hold on;
for i = 1:length(ice_thickness_values)
    plot(x, tension_curves(i, :), 'LineWidth', 1.5);
end
grid on;
title('Tension Distribution for Varying Ice Thickness');
xlabel('Horizontal Position (m)');
ylabel('Tension (N)');
legend(arrayfun(@(t) sprintf('Ice Thickness: %.2f m', t), ice_thickness_values, 'UniformOutput', false), 'Location', 'Best');
hold off;

% Plot Sag Curves (Varying Wind Speed)
figure;
hold on;
for i = 1:length(wind_speed_values)
    plot(x, sag_curves(i, :), 'LineWidth', 1.5);
end
grid on;
title('Sag Profile for Varying Wind Speeds');
xlabel('Horizontal Position (m)');
ylabel('Sag (m)');
legend(arrayfun(@(w) sprintf('Wind Speed: %.1f m/s', w), wind_speed_values, 'UniformOutput', false), 'Location', 'Best');
hold off;

% Bar Plot: Maximum Tension for Different Ice Thicknesses
figure;
bar(ice_thickness_values, max_tension, 'FaceColor', [0.2, 0.6, 0.8]);
grid on;
title('Maximum Tension for Varying Ice Thickness');
xlabel('Ice Thickness (m)');
ylabel('Maximum Tension (N)');
xticks(ice_thickness_values);
xticklabels(arrayfun(@(t) sprintf('%.2f', t), ice_thickness_values, 'UniformOutput', false));

% Bar Plot: Maximum Sag for Different Wind Speeds
figure;
bar(wind_speed_values, max_sag, 'FaceColor', [0.8, 0.4, 0.4]);
grid on;
title('Maximum Sag for Varying Wind Speeds');
xlabel('Wind Speed (m/s)');
ylabel('Maximum Sag (m)');
xticks(wind_speed_values);
xticklabels(arrayfun(@(w) sprintf('%.1f', w), wind_speed_values, 'UniformOutput', false));

% Print Results
fprintf('--- Summary of Results ---\n');
fprintf('Tension Analysis:\n');
for i = 1:length(ice_thickness_values)
    fprintf('Ice Thickness: %.2f m -> Max Tension: %.2f N\n', ice_thickness_values(i), max_tension(i));
end
fprintf('\nSag Analysis:\n');
for i = 1:length(wind_speed_values)
    fprintf('Wind Speed: %.1f m/s -> Max Sag: %.2f m\n', wind_speed_values(i), max_sag(i));
end 

