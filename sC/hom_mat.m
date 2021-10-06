function out = hom_mat(rotm,disp)
            %
            % Data matrice di rotazione e vettore di traslazione riporta la
            % matrice di rototraslazione omogenea associata
            %
if isa(rotm, 'ADvar') || isa(disp, 'ADvar')
    b = ADvar([0 0 0 1], [0 0 0 0]);
else
b = [0 0 0 1];
end

out = [rotm disp; b];
            
end