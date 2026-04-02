clear all;
close all;

set_plot_style();


result_no_osd_no_corr = load("data/results_paper/200_10_10_sogrand_osd0_corr0.mat").result;
result_osd_no_corr = load("data/results_paper/200_10_10_sogrand_osd1_corr0.mat").result;
result_osd_corr = load("data/results_paper/200_10_10_sogrand_osd1_corr1.mat").result;

result_bp_no_osd = load("data/results_paper/200_10_10_bp_osdFalse_iter100.mat");
result_bp_osd = load("data/results_paper/200_10_10_bp_osdTrue_iter100.mat");

fig = figure();
fig.Units = 'centimeters';
fig.Position = [1, 1, 6.2, 6.5]; 
fig.Color = 'w';

hold on;
loglog(result_bp_no_osd.p_list, result_bp_no_osd.bler_deg_list, "-s", 'MarkerFaceColor', 'w', "DisplayName", "BP", "Color", get_color(1));
loglog(result_bp_osd.p_list, result_bp_osd.bler_deg_list, "-s", 'MarkerFaceColor', 'w', "DisplayName", "BP+OSD", "Color", get_color(2));
plot(NaN, NaN, 'LineStyle', 'none', 'DisplayName', ' ', 'HandleVisibility', 'on');

loglog(result_no_osd_no_corr.p_list, result_no_osd_no_corr.bler_deg_list, '-o', 'MarkerFaceColor', get_color(1), "DisplayName", "SOGRAND", "Color", get_color(1));
loglog(result_osd_no_corr.p_list, result_osd_no_corr.bler_deg_list, '-o', 'MarkerFaceColor', get_color(2), "DisplayName", "SOGRAND+OSD", "Color", get_color(2));
loglog(result_osd_corr.p_list, result_osd_corr.bler_deg_list, '-o', 'MarkerFaceColor', get_color(3), "DisplayName", "SOGRAND+OSD*", "Color", get_color(3));



set(gca, 'XScale', 'log', 'YScale', 'log');
grid();
ylabel("Logical BLER P_{L}");
xlabel("Physical Error Rate p");
lgd = legend('show', 'Location','northoutside', "NumColumns", 2);
lgd.ItemTokenSize = [20, 6];

ylim([1e-7, 1]);
xlim([1e-3, 10^(-0.8)]);
xticks([1e-3, 1e-2, 1e-1]);
exportgraphics(gcf, 'bler_curve_200_10_10_small.pdf', 'ContentType', 'vector');
