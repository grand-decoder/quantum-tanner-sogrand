function color = get_color(idx)
    % Define the palette
    colors = [0.90, 0.62, 0.00;  % Orange
              0.33, 0.70, 0.91;  % Sky Blue
              0.00, 0.62, 0.45;  % Green
              0.83, 0.36, 0.00]; % Vermillion

    color = colors(mod(idx - 1, size(colors, 1)) + 1, :);
end