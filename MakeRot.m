function rot = MakeRot(n, angle)
%Warning: these are the reverse of conventional argument order:*)
%change order once all is consistent to (angle, nhat).*)

    c = cos(angle);
    s = sin(angle);
    cm = 1 - c;

    rot = [[c+cm*n(1)^2,          cm*n(2)*n(1)-s*n(3),  cm*n(3)*n(1)+s*n(2)]
           [cm*n(1)*n(2)+s*n(3),  c+cm*n(2)^2,          cm*n(3)*n(2)-s*n(1)]
           [cm*n(1)*n(3)-s*n(2),  cm*n(2)*n(3)+s*n(1),  c+cm*n(3)^2]];
end
