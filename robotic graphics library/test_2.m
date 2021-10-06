k = 1;
j1 = Joint(1, 2.5, 15, TrotX(pi/2));
j2 = Joint(0, 2.5, [-20, 10],Ttz(30)*TrotX(pi/2));
j3 = Joint(1, 4, 5, Ttz(30));
j4 = Joint(1, 2.5, 8, Ttx(20));

gyroscope = SerialChain('local', eye(4), {[j1, j2, j3, j4]});

gyroscope.plotSerial()

theta = linspace(0, 2*pi, 360);
dist = sin(10.*theta).*4;
i = 0;

while true
    
    i = i + 1;
    gyroscope.updateJoints([theta(i), dist(i), 10*theta(i), 0]);
    
    if i == 360
        i = 0;
    end
    pause(0.05)
end