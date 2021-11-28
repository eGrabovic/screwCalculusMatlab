function projection = proj4(pt4D, angle) % default value for angle: 0 
%     Chop({pt4D((1)),pt4D((2)),Cos(angle)*pt4D((3)) + Sin(angle)*pt4D((4))})//N
% Chop should substitute values close to 0 with actual 0
% // N just returns the numerical value of the expression
    if ~exist(angle, 'var')
        angle = 0;
    end
    projection = [pt4D(1),...
                  pt4D(2),...
                  cos(angle)*pt4D(3) + Sin(angle)*pt4D(4)];
end
