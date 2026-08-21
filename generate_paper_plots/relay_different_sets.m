clear all
close all;

set_plot_style();

fig = figure();
fig.Units = 'centimeters';
fig.Position = [1, 1, 8.8, 7]; 
fig.Color = 'w';

result_bp_5 = load("data/results_paper/250_10_15_bp_relay_set5.mat");
result_sogrand_5 = load("data/results_paper/250_10_15_sogrand_relay_set5.mat").result;

result_bp_20 = load("data/results_paper/250_10_15_bp_relay_set20.mat");
result_sogrand_20 = load("data/results_paper/250_10_15_sogrand_relay_set20.mat").result;

result_bp_60 = load("data/results_paper/250_10_15_bp_relay_iter3620.mat");
result_sogrand_60 = load("data/results_paper/250_10_15_sogrand_relay_iter3620.mat").result;

hold on;
loglog(result_bp_5.p_list, result_bp_5.bler_deg_list, "-s", 'MarkerFaceColor', 'w', "DisplayName", "BP+Relay (R=5)", "Color", get_color(1));
loglog(result_bp_20.p_list, result_bp_20.bler_deg_list, "-s", 'MarkerFaceColor', 'w', "DisplayName", "BP+Relay (R=20)", "Color", get_color(3));
loglog(result_bp_60.p_list, result_bp_60.bler_deg_list, "-s", 'MarkerFaceColor', 'w', "DisplayName", "BP+Relay (R=60)", "Color", get_color(4));

loglog(result_sogrand_5.p_list, result_sogrand_5.bler_deg_list, '-o', 'MarkerFaceColor', get_color(1), "DisplayName", "SOGRAND+Relay (R=5)", "Color", get_color(1));
loglog(result_sogrand_20.p_list, result_sogrand_20.bler_deg_list, '-o', 'MarkerFaceColor', get_color(3), "DisplayName", "SOGRAND+Relay (R=20)", "Color", get_color(3));
loglog(result_sogrand_60.p_list, result_sogrand_60.bler_deg_list, '-o', 'MarkerFaceColor', get_color(4), "DisplayName", "SOGRAND+Relay (R=60)", "Color", get_color(4));

legend("Location", "northoutside", "NumColumns", 2);

ylim([1e-7, 1]);
xlim([1e-2, 10^(-0.8)]);
grid();
set(gca, 'XScale', 'log', 'YScale', 'log');
ylabel("Logical BLER  P_{L}");
xlabel("Physical Error Rate p");

exportgraphics(gcf, 'bler_curve_relay_different_sets.pdf', 'ContentType', 'vector');
