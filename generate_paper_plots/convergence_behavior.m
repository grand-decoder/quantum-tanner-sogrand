clear all
close all;

result_sogrand_no_osd_no_corr = load("data/results_paper/250_10_15_sogrand_convergence.mat").result;
result_bp_no_osd = load("data/results_paper/250_10_15_bp_convergence.mat");
result_bp_no_osd.bler_deg_per_iteration_list = result_bp_no_osd.bler_deg_list';

set_plot_style();

fig = figure;
fig.Units = 'centimeters';
fig.Position = [1, 1, 8.5, 5.5]; 
fig.Color = 'w';

ax1 = gca;

hold on;

[~, ~, iteration_list_pseudo, threshold_list_pseudo] = create_curves(result_bp_no_osd);
semilogy(ax1, iteration_list_pseudo, threshold_list_pseudo, '-', 'MarkerFaceColor', 'w', "DisplayName", "BP", "Color", get_color(1));
[~, ~, iteration_list_pseudo, threshold_list_pseudo] = create_curves(result_sogrand_no_osd_no_corr);
semilogy(ax1, iteration_list_pseudo, threshold_list_pseudo, '-', 'MarkerFaceColor', get_color(1), "DisplayName", "SOGRAND", "Color", get_color(2));

grid(ax1, 'on');
xlim(ax1, [0, 50]);
ylim(ax1, [1e-2, 10^(-0.7)]);
yticks(ax1, [1e-2, 1e-1]);
xlabel(ax1, 'Max. Iterations N_{iter}');
ylabel(ax1, 'Pseudothreshold p_{th}');
lgd = legend(ax1, 'Location', 'southeast');
set(ax1, 'YScale', 'log', 'XScale', 'linear');

exportgraphics(gcf, 'pseudo_threshold_over_iteration.pdf', 'ContentType', 'vector');

fig = figure;
fig.Units = 'centimeters';
fig.Position = [1, 1, 8.5, 5.5]; 
fig.Color = 'w';

ax2 = gca;

p_idx = 13;
hold on;
semilogy(ax2, 1:size(result_bp_no_osd.bler_deg_per_iteration_list, 2), ...
    result_bp_no_osd.bler_deg_per_iteration_list(p_idx, :), ...
    '-', 'MarkerFaceColor', 'w', "Color", get_color(1), "DisplayName", "BP");

semilogy(ax2, 1:size(result_sogrand_no_osd_no_corr.bler_deg_per_iteration_list, 2), ...
    result_sogrand_no_osd_no_corr.bler_deg_per_iteration_list(p_idx, :), ...
    '-', 'MarkerFaceColor', get_color(1), "Color", get_color(2), "DisplayName", "SOGRAND");

set(ax2, 'YScale', 'log', 'XScale', 'linear');
grid(ax2, 'on');
ylim(ax2, [1e-7, 1]);
yticks(ax2, [1e-7, 1e-4, 1e-1]);
xlim(ax2, [0, 200]);
xlabel(ax2, 'Max. Iterations N_{iter}');
ylabel(ax2, 'P_{L} at p=10^{-1.8}');
lgd = legend(ax2, 'Location', 'southwest');

assert(abs(result_sogrand_no_osd_no_corr.p_list(p_idx) - 10^(-1.8)) < eps);
assert(abs(result_bp_no_osd.p_list(p_idx) - 10^(-1.8)) < eps);

xZoom = [0 8];
yZoom = [1e-3 1];

mainAx = gca;

rectangle('Position', [xZoom(1), yZoom(1), diff(xZoom), diff(yZoom)], ...
          'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.75)

insetAx = axes('Position', [0.62 0.63 0.27 0.27]);
box on;
hold on;
semilogy(1:size(result_bp_no_osd.bler_deg_per_iteration_list, 2), result_bp_no_osd.bler_deg_per_iteration_list(p_idx, :), "-x", "Color", get_color(1));
semilogy(1:size(result_sogrand_no_osd_no_corr.bler_deg_per_iteration_list, 2), result_sogrand_no_osd_no_corr.bler_deg_per_iteration_list(p_idx, :), "-x", "Color", get_color(2));
set(insetAx, 'YScale', 'log', 'XScale', 'linear');
xlim(xZoom)
ylim(yZoom)
grid on
insetAx.XMinorGrid = 'off';
insetAx.YMinorGrid = 'off';
xticks([0, 4, 8]);
yticks([1e-3, 1e-2, 1e-1]);

pos1 = dataToFigureLog(mainAx, [xZoom(2), yZoom(2)]);
pos2 = dataToFigureLog(mainAx, [xZoom(1), yZoom(1)]);

insetPos = insetAx.Position;
insetTopLeft    = [insetPos(1), insetPos(2)+insetPos(4)];
insetBottomLeft = [insetPos(1), insetPos(2)];

annotation('line', [pos1(1) insetTopLeft(1)],    [pos1(2) insetTopLeft(2)],    'Color', [0.5 0.5 0.5], 'LineWidth', 0.75)
annotation('line', [pos2(1) insetBottomLeft(1)], [pos2(2) insetBottomLeft(2)], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.75)

exportgraphics(gcf, 'bler_at_1eminus1_8_over_iteration.pdf', 'ContentType', 'vector');

function figPos = dataToFigureLog(ax, dataPoint)
    axPos = ax.Position;
    xl = ax.XLim;
    yl = ax.YLim;
    xFrac = (dataPoint(1) - xl(1)) / diff(xl);
    yFrac = (log10(dataPoint(2)) - log10(yl(1))) / (log10(yl(2)) - log10(yl(1)));
    figPos = [axPos(1) + xFrac*axPos(3), axPos(2) + yFrac*axPos(4)];
end

function [iteration_list_1e_2, threshold_list_1e_2, iteration_list_pseudo, threshold_list_pseudo] = create_curves(result)
    cut_indices_right = 7;
    bler_threshold = 1e-2;
    threshold_line_y = [bler_threshold, bler_threshold];
    threshold_line_x = [result.p_list(1), result.p_list(end)];
    
    [iteration_list_1e_2, threshold_list_1e_2] = extract_thresholds(threshold_line_x, threshold_line_y, result.p_list, result.bler_deg_per_iteration_list, cut_indices_right, false);
    
    k = 10;
    bler_uncoded = 1 - (1 - result.p_list').^k;
    
    [iteration_list_pseudo, threshold_list_pseudo] = extract_thresholds(result.p_list, bler_uncoded, result.p_list, result.bler_deg_per_iteration_list, cut_indices_right, false);
end


function [iteration_list, threshold_list] = extract_thresholds(threshold_line_x, threshold_line_y, p_list, bler_deg_per_iteration_list, cut_indices_right, plot_results)
    iteration_list = [];
    threshold_list = [];    
    
    p_list_cut = p_list(1:end-cut_indices_right);
    
    if plot_results
        figure;
        hold on;
        plot(log10(threshold_line_x), log10(threshold_line_y), "-", "DisplayName", "u");
    end
    for i=1:size(bler_deg_per_iteration_list, 2)
        bler_deg = bler_deg_per_iteration_list(1:end-cut_indices_right, i);
        [x0,y0] = polyxpoly(log10(threshold_line_x), log10(threshold_line_y), log10(p_list_cut), log10(bler_deg));

        if plot_results
            plot(x0, y0, "x");
            plot(log10(p_list_cut), log10(bler_deg), "-", "DisplayName", string(i));
        end
        
        if length(x0) ~= 1
            continue;
            % error("not one intersection");
        end
    
        p_threshold = 10^x0;
    
        iteration_list = [iteration_list, i];
        threshold_list = [threshold_list, p_threshold];
    end
    
    if plot_results
        grid();
        legend('Location', 'southeast');
        xlabel("Physical Error Rate p");
        ylabel("BLER");
    end
end