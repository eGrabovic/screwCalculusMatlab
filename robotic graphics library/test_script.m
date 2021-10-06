%% animation scripting

a = 1;
alpha  = linspace(0, 2*pi, 360);
j1 = Joint('r', 2.5, 15, sC.hom_mat(sC.rotX(pi/2),[0; 0; 0]));
j1.plotJoint();
view(45,45);
counter = 0;

while a == 1
    
    if counter == 360
        counter = 0;
    end
    
    counter = counter + 1;
    matrix = sC.hom_mat(sC.rotX(alpha(counter)),[alpha(counter);alpha(counter).^2 - alpha(counter);3]);
    j1.changeFrameMat(matrix)
    pause(0.01)
    drawnow
    
end