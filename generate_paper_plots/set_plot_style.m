function set_publication_style()
    % --- Global Publication Defaults ---
    set(groot, 'DefaultLineLineWidth', 1.2);      % Thicker lines
    set(groot, 'DefaultLineMarkerSize', 6);       % Visible markers
    set(groot, 'DefaultAxesFontName', 'Arial');   % Clean font
    set(groot, 'DefaultAxesFontSize', 9);         % Professional scale
    set(groot, 'DefaultAxesBox', 'on');           % Contained plots

    set(groot, 'DefaultTextInterpreter', 'tex');
    set(groot, 'DefaultLegendInterpreter', 'tex');
end