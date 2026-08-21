%% SOGRAND Memory Interval Optimization
close all;
clear all;

set_plot_style();

fig = figure();
fig.Units = 'centimeters';
fig.Position = [1, 1, 6.2, 4.5];
fig.Color = 'w';

result = load("data/results_paper/250_10_15_sogrand_relay_mem_interval_opt.mat").result;
x = result.sweep_list;
y = result.sweep_list_2;

imagesc(x, y, log10(result.bler_deg_list'));
set(gca, 'YDir', 'normal'); 
colormap(parula);
cb = colorbar;
cb.Label.String = "Logical BLER (log_{10})";

xlabel('Center \gamma_{c}');
ylabel('Width \gamma_{w}');

yticks(0:0.4:1.3);
xticks(-0.2:0.1:0.2);

cb.Ticks = [-4, -3, -2];
cb.TickLabels = {'10^{-4}', '10^{-3}', '10^{-2}'};

cb.Ticks = [-4, -3, -2];
cb.TickLabels = {'10^{-4}', '10^{-3}', '10^{-2}'};

[minVal, minIdx] = min(result.bler_deg_list(:));
[row, col] = ind2sub(size(result.bler_deg_list), minIdx);

hold on;
plot(x(row), y(col), 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'none', 'LineWidth', 1.5);

exportgraphics(gcf, 'relay_sweep_interval_sogrand.pdf', 'ContentType', 'vector');

%% BP Memory Interval Optimization
fig = figure();
fig.Units = 'centimeters';
fig.Position = [1, 1, 6.2, 4.5];
fig.Color = 'w';

result = load("data/results_paper/250_10_15_bp_relay_mem_interval_opt.mat");
x = result.center_list;
y = result.width_list;

bler_log10  = log10(result.bler_deg_list');

imagesc(x, y, log10(result.bler_deg_list'));
set(gca, 'YDir', 'normal'); 

colormap(parula);
cb = colorbar;
cb.Label.String = "Logical BLER (log_{10})";

xlabel('Center \gamma_{c}');
ylabel('Width \gamma_{w}');

cb.Ticks = [-4, -3, -2];
cb.TickLabels = {'10^{-4}', '10^{-3}', '10^{-2}'};

[minVal, minIdx] = min(result.bler_deg_list(:));
[row, col] = ind2sub(size(result.bler_deg_list), minIdx);

hold on;
plot(x(row), y(col), 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
exportgraphics(gcf, 'relay_sweep_interval_bp.pdf', 'ContentType', 'vector');

% Gamm0 Optimization
fig = figure();
fig.Units = 'centimeters';
fig.Position = [1, 1, 6.2, 4.5];
fig.Color = 'w';

result_sogrand = load("data/results_paper/250_10_15_sogrand_relay_gamma0_opt.mat").result;
result_bp = load("data/results_paper/250_10_15_bp_relay_gamma0_opt.mat");

hold on;
semilogy(result_bp.gamma0_list, result_bp.bler_deg_list, "-", "DisplayName", "BP", "Color", get_color(1));
semilogy(result_sogrand.sweep_list, result_sogrand.bler_deg_list, "-", "DisplayName", "SOGRAND", "Color", get_color(2));

ylim([1e-4, 1e-2]);
xlim([0, 0.25]);
grid();
set(gca, 'YScale', 'log');
xlabel("Uniform Memory Strength \gamma");
ylabel("Logical BLER P_{L}");
legend("Location", "southeast");
xline(0.07, "--", "DisplayName", "Optimum", "LineWidth", 1.2);

exportgraphics(gcf, 'relay_sweep_gamma0.pdf', 'ContentType', 'vector');