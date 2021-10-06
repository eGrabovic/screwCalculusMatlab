function plotFrame(ax, T, NameValueArgs)
    arguments
        ax
        T
        NameValueArgs.LineWidth (1,1) {mustBeNumeric} = 1
    end
    LineWidth = NameValueArgs.LineWidth;
    pos = T(1:3,4);
    xaxis = [pos, pos + T(1:3,1)];
    yaxis = [pos, pos + T(1:3,2)];
    zaxis = [pos, pos + T(1:3,3)];
    line(ax, xaxis(1,:), xaxis(2,:), xaxis(3,:), 'color', 'r', 'LineWidth', LineWidth);
    line(ax, yaxis(1,:), yaxis(2,:), yaxis(3,:), 'color', 'g', 'LineWidth', LineWidth);
    line(ax, zaxis(1,:), zaxis(2,:), zaxis(3,:), 'color', 'b', 'LineWidth', LineWidth);
end