function points = createParallelepiped(edge, len)
 % TODO

 % TODO: why as method?
    a = edge;

    if length(len) == 2
        h1 = len(1);
        h2 = len(2);

    else
        h1 = 0;
        h2 = len;
    end

    P1 = [ a/2; -a/2; -h1]; 
    P2 = [ a/2;  a/2; -h1]; 
    P3 = [ a/2;  a/2;  h2]; 
    P4 = [ a/2; -a/2;  h2];
    P5 = [-a/2;  a/2; -h1]; 
    P6 = [-a/2;  a/2;  h2]; 
    P7 = [-a/2; -a/2; -h1]; 
    P8 = [-a/2; -a/2;  h2];

    points = [P1,P2,P3,P4,P5,P6,P7,P8];
    
end