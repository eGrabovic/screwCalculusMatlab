% script cradle machine graphics

% upper body joints  type : 1 REVOLUTE; 0 PRISMATIC
jt(1) = Joint(1, 2.5, 8, eye(4));
jt(2) = Joint(0, 2.5, 8, [rotY(pi/2), [0;0;0]; 0,0,0,1]);
jt(3) = Joint(1, 2.5, 8, [rotY(pi/2), [0;0;0]; 0,0,0,1]);
jt(4) = Joint(1, 2.5, 8, [rotY(-pi/2), [0;0;0]; 0,0,0,1]);
% lower body joints
jg(1) = Joint(0, 2.5, [-8, 14], [rotX(pi/2), [0;0;0]; 0,0,0,1]);
jg(2) = Joint(0, 2.5, [-8, 14], [rotX(-pi/2), [0;0;0]; 0,0,0,1]);
jg(3) = Joint(1, 2.5, 8, [rotX(pi/2), [0;0;0]; 0,0,0,1]);
jg(4) = Joint(0, 2.5, [-8, 14], [rotY(pi/2), [0;0;0]; 0,0,0,1]);
jg(5) = Joint(1, 2.5, 8, eye(4));
% upperbody serial
upperCradleMachine = SerialChain('local', [rotX(pi/2), [0;0;0]; 0,0,0,1],{jt}); % [rotX(-pi/2), [0;0;0]; 0,0,0,1]
%lower body serial
lowerCradleMachine = SerialChain('local', [rotX(pi/2), [0;0;0]; 0,0,0,1],{jg});
% plot 
ax = upperCradleMachine.plotSerial();
lowerCradleMachine.plotSerial(ax)
% update joint values
upperCradleMachine.updateJoints([pi/4, 50, 0, 0])
lowerCradleMachine.updateJoints([40, 40, pi/10, 40, pi/4])
% extra axis
plot3([0,60], [0, 0], [0,0], 'color', 'k', 'LineStyle', '--')
plot3([0,0], [-60, 60], [0,0], 'color', 'k', 'LineStyle', '--')
plot3([0,0], [0, 0], [-60, 60], 'color', 'k', 'LineStyle', '--')
rootPointAxis = [lowerCradleMachine.links{3}.XData(1); lowerCradleMachine.links{3}.YData(1); lowerCradleMachine.links{3}.ZData(1)];
plot3([rootPointAxis(1),rootPointAxis(1) + 60], ...
    [rootPointAxis(2),rootPointAxis(2)],...
    [rootPointAxis(3),rootPointAxis(3)],...
    'color', 'k', 'LineStyle', '--')
view(-45,45)