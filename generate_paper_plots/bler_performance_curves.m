close all;

set_plot_style();

% [[250, 10, 15]]
fig = figure();
fig.Units = 'centimeters';
fig.Position = [1, 1, 8.5, 8]; 
fig.Color = 'w';

result_250_bp = load("data\results_paper\250_10_15_bp.mat");
result_250_bp_osd_9 = load("data\results_paper\250_10_15_bp_osd_order9.mat");
result_250_bp_osd_0 = load("data\results_paper\250_10_15_bp_osd_order0.mat");
result_250_bp_qbp = readtable('data\results_paper\250_10_15_bp_xz\bler_curve.txt', 'Delimiter', ' ', 'FileType', 'text');

result_250_sogrand = load("data\results_paper\250_10_15_sogrand.mat").result;
result_250_sogrand_osd_9 = load("data\results_paper\250_10_15_sogrand_osd_order9.mat").result;
result_250_sogrand_osd_0 = load("data\results_paper\250_10_15_sogrand_osd_order0.mat").result;
result_250_sogrand_qbp = load('data\results_paper\250_10_15_sogrand_xz.mat').result;

hold on;
loglog(result_250_bp.p_list, result_250_bp.bler_deg_list, "-s", 'MarkerFaceColor', 'w', "DisplayName", "BP", "Color", get_color(1));
loglog(result_250_bp_osd_0.p_list, result_250_bp_osd_0.bler_deg_list, "-s", 'MarkerFaceColor', 'w', "DisplayName", "BP+OSD-0", "Color", get_color(2));
loglog(result_250_bp_osd_9.p_list, result_250_bp_osd_9.bler_deg_list, "-s", 'MarkerFaceColor', 'w', "DisplayName", "BP+OSD-9", "Color", get_color(3));
loglog(result_250_bp_qbp.p_list, result_250_bp_qbp.bler_deg_list, "-s", 'MarkerFaceColor', 'w', "DisplayName", "BP+XZ", "Color", get_color(4));
loglog(result_250_sogrand.p_list, result_250_sogrand.bler_deg_list, '-o', 'MarkerFaceColor', get_color(1), "DisplayName", "SOGRAND", "Color", get_color(1));
loglog(result_250_sogrand_osd_0.p_list, result_250_sogrand_osd_0.bler_deg_list, '-o', 'MarkerFaceColor', get_color(2), "DisplayName", "SOGRAND+OSD-0", "Color", get_color(2));
loglog(result_250_sogrand_osd_9.p_list, result_250_sogrand_osd_9.bler_deg_list, '-o', 'MarkerFaceColor', get_color(3), "DisplayName", "SOGRAND+OSD-9", "Color", get_color(3));
loglog(result_250_sogrand_qbp.p_list, result_250_sogrand_qbp.bler_deg_list, '-o', 'MarkerFaceColor', get_color(4), "DisplayName", "SOGRAND+XZ", "Color", get_color(4));

set(gca, 'XScale', 'log', 'YScale', 'log');
hold off;
grid on;

lgd = legend('show', 'Location','northoutside', "NumColumns", 2);
lgd.ItemTokenSize = [20, 6];

ylim([1e-7, 1]);
ylabel("Logical BLER  P_{L}");
xlabel("Physical Error Rate p");
exportgraphics(gcf, 'bler_curve_250_10_15_all.pdf', 'ContentType', 'vector');

% [[200, 10, 10]]
fig = figure();
fig.Units = 'centimeters';
fig.Position = [1, 1, 8.5, 8]; 
fig.Color = 'w';

result_200_bp = load("data\results_paper\200_10_10_bp.mat");
result_200_bp_osd_9 = load("data\results_paper\200_10_10_bp_osd_order9.mat");
result_200_bp_osd_0 = load("data\results_paper\200_10_10_bp_osd_order0.mat");
result_200_bp_qbp = readtable('data\results_paper\200_10_10_bp_xz\bler_curve.txt', 'Delimiter', ' ', 'FileType', 'text');

result_200_sogrand = load("data\results_paper\200_10_10_sogrand.mat").result;
result_200_sogrand_osd_9 = load("data\results_paper\200_10_10_sogrand_osd_order9.mat").result;
result_200_sogrand_osd_0 = load("data\results_paper\200_10_10_sogrand_osd_order0.mat").result;
result_200_sogrand_qbp = load('data\results_paper\200_10_10_sogrand_xz.mat').result;

hold on;
loglog(result_200_bp.p_list, result_200_bp.bler_deg_list, "-s", 'MarkerFaceColor', 'w', "DisplayName", "BP", "Color", get_color(1));
loglog(result_200_bp_osd_0.p_list, result_200_bp_osd_0.bler_deg_list, "-s", 'MarkerFaceColor', 'w', "DisplayName", "BP+OSD-0", "Color", get_color(2));
loglog(result_200_bp_osd_9.p_list, result_200_bp_osd_9.bler_deg_list, "-s", 'MarkerFaceColor', 'w', "DisplayName", "BP+OSD-9", "Color", get_color(3));
loglog(result_200_bp_qbp.p_list, result_200_bp_qbp.bler_deg_list, "-s", 'MarkerFaceColor', 'w', "DisplayName", "BP+XZ", "Color", get_color(4));

loglog(result_200_sogrand.p_list, result_200_sogrand.bler_deg_list, '-o', 'MarkerFaceColor', get_color(1), "DisplayName", "SOGRAND", "Color", get_color(1));
loglog(result_200_sogrand_osd_0.p_list, result_200_sogrand_osd_0.bler_deg_list, '-o', 'MarkerFaceColor', get_color(2), "DisplayName", "SOGRAND+OSD-0", "Color", get_color(2));
loglog(result_200_sogrand_osd_9.p_list, result_200_sogrand_osd_9.bler_deg_list, '-o', 'MarkerFaceColor', get_color(3), "DisplayName", "SOGRAND+OSD-9", "Color", get_color(3));
loglog(result_200_sogrand_qbp.p_list, result_200_sogrand_qbp.bler_deg_list, '-o', 'MarkerFaceColor', get_color(4), "DisplayName", "SOGRAND+XZ", "Color", get_color(4));

set(gca, 'XScale', 'log', 'YScale', 'log');
hold off;
grid on;

lgd = legend('show', 'Location','northoutside', "NumColumns", 2);
lgd.ItemTokenSize = [20, 6];

ylim([1e-7, 1]);
ylabel("Logical BLER  P_{L}");
xlabel("Physical Error Rate p");
exportgraphics(gcf, 'bler_curve_200_10_10_all.pdf', 'ContentType', 'vector');