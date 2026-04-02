clear all
close all;

% SOGRAND
result_sogrand_no_osd_no_corr = load("data/results_paper/250_10_15_sogrand_osd0_corr0_threshold.mat").result;


% BP
result_bp_no_osd = load("data/results_paper/250_10_15_bp_osdFalse_iter100_threshold.mat");
result_bp_no_osd.bler_deg_per_iteration_list = result_bp_no_osd.bler_deg_list';

set_plot_style();

% Create the figure
fig = figure;
fig.Units = 'centimeters';
fig.Position = [1, 1, 5, 6.5]; 
fig.Color = 'w';

% Create a 2x1 layout
t = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Top Plot (Now Pseudothreshold) ---
ax1 = nexttile;
hold on;


[iteration_list_1e_2, threshold_list_1e_2, iteration_list_pseudo, threshold_list_pseudo] = create_curves(result_bp_no_osd);
semilogy(ax1, iteration_list_pseudo, threshold_list_pseudo, '-s', 'MarkerFaceColor', 'w', "DisplayName", "BP", "Color", get_color(1));
[iteration_list_1e_2, threshold_list_1e_2, iteration_list_pseudo, threshold_list_pseudo] = create_curves(result_sogrand_no_osd_no_corr);
semilogy(ax1, iteration_list_pseudo, threshold_list_pseudo, '-o', 'MarkerFaceColor', get_color(1), "DisplayName", "SOGRAND", "Color", get_color(1));

set(ax1, 'YScale', 'log', 'XScale', 'linear');
grid(ax1, 'on');
xlim(ax1, [0, 12]);
xticklabels(ax1, {}); % Hide x-labels for the top plot
ylabel(ax1, 'Pseudothreshold p_{th}');

ylim(ax1, [1e-2, 10^(-0.7)]);
yticks(ax1, [1e-2, 1e-1]);

% --- Bottom Plot (Now BLER) ---
ax2 = nexttile;
hold on;
semilogy(ax2, 1:size(result_bp_no_osd.bler_deg_per_iteration_list, 2), ...
    result_bp_no_osd.bler_deg_per_iteration_list(1, :), ...
    '-s', 'MarkerFaceColor', 'w', "Color", get_color(1), "DisplayName", "BP");

semilogy(ax2, 1:size(result_sogrand_no_osd_no_corr.bler_deg_per_iteration_list, 2), ...
    result_sogrand_no_osd_no_corr.bler_deg_per_iteration_list(1, :), ...
    '-o', 'MarkerFaceColor', get_color(1), "Color", get_color(1), "DisplayName", "SOGRAND");

set(ax2, 'YScale', 'log', 'XScale', 'linear');
grid(ax2, 'on');
ylim(ax2, [1e-7, 1]);
yticks(ax2, [1e-7, 1e-4, 1e-1]);
xlim(ax2, [0, 10]);
xlabel(ax2, 'Max. Iterations N_{iter}');
ylabel(ax2, 'P_{L} at p=10^{-3}');
% legend(ax2, 'Location', 'best');
lgd = legend(ax2, 'Location', 'northeast');
currentPos = lgd.Position;
lgd.Position = [currentPos(1)+0.02, currentPos(2)+0.02, currentPos(3), currentPos(4)];

% --- Link the axes ---
linkaxes([ax1, ax2], 'x');

exportgraphics(gcf, 'threshold_over_iteration_new_small.pdf', 'ContentType', 'vector');

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