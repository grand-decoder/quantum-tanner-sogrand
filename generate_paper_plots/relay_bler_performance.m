clear all
close all;

set_plot_style();

fig = figure(1);
fig.Units = 'centimeters';
fig.Position = [1, 1, 8.5, 6.5]; 
fig.Color = 'w';

result_sogrand = load("data/results_paper/250_10_15_sogrand_iter3620.mat").result;
result_sogrand_osd = load("data/results_paper/250_10_15_sogrand_osd_iter3620.mat").result;
result_sogrand_relay = load("data/results_paper/250_10_15_sogrand_relay_iter3620.mat").result;

result_bp = load("data/results_paper/250_10_15_bp_3620iter.mat");
result_bp_osd = load("data/results_paper/250_10_15_bp_osd_iter3620.mat"); 
result_bp_relay = load("data/results_paper/250_10_15_bp_relay_iter3620.mat");

hold on;
plot(result_bp.p_list, result_bp.bler_deg_list, '-s', 'MarkerFaceColor', 'w', "DisplayName", "BP", "Color", get_color(1));
plot(result_sogrand.p_list, result_sogrand.bler_deg_list, '-o', 'MarkerFaceColor', get_color(1), 'MarkerFaceColor', get_color(1), "DisplayName", "SOGRAND", "Color", get_color(1));

plot(result_bp_osd.p_list, result_bp_osd.bler_deg_list, '-s', 'MarkerFaceColor', 'w', "DisplayName", "BP+OSD-9", "Color", get_color(3));
plot(result_sogrand_osd.p_list, result_sogrand_osd.bler_deg_list, '-o', 'MarkerFaceColor', get_color(1), 'MarkerFaceColor', get_color(3), "DisplayName", "SOGRAND+OSD-9", "Color", get_color(3));

plot(result_bp_relay.p_list, result_bp_relay.bler_deg_list, "-s", 'MarkerFaceColor', 'w', "DisplayName", "BP+Relay", "Color", get_color(4));
plot(result_sogrand_relay.p_list, result_sogrand_relay.bler_deg_list, '-o', 'MarkerFaceColor', get_color(4), "DisplayName", "SOGRAND+Relay", "Color", get_color(4));
ylim([1e-7, 1]);
xlim([1e-2, 10^(-0.8)]);
grid();
set(gca, 'XScale', 'log', 'YScale', 'log');
lgd = legend("Location", "northoutside", "NumColumns", 2, 'Orientation', 'horizontal');

ylabel("Logical BLER  P_{L}");
xlabel("Physical Error Rate p");
exportgraphics(gcf, 'bler_curve_relay.pdf', 'ContentType', 'vector');