function [] = plotkite(kiteshape, CoP, CoG)
%PLOTKITE Plots the shape of a simple diamond kite along with CoP and CoG*
%   Input dimensions of the cross beams of your kite** and gaze in awe as my
%   function quickly plots the shape! Mostly used to clean up code in main
%   script. 
%
%   *only input the x-coordinate of CoP/CoG
%   **in the form of a 2x3 double: [y1,y2,y3; x1,x2,x3]

%Governing geometry
L = kiteshape(2,3);
cline = [0 0;0 L];

%plotting@
figure()
    %plot the kite and a dashed centerline
    plot(kiteshape(2,:),kiteshape(1,:),'Color','blue'); axis equal
    hold on; plot(cline(2,:),cline(1,:), '--', 'Color','black'); 
    plot(kiteshape(2,:),-kiteshape(1,:),'Color','blue');
    plot([L/4,L/4],[-kiteshape(1,2),kiteshape(1,2)],'--','Color','black')

    %plot CoP and CoG
    plot(CoP,0, '.',"Color",'red','MarkerSize',20)
    plot(CoG,0,'.','Color','green','MarkerSize',20)
    legend('','','','', 'Center of Pressure','Center of Gravity')

    xlabel('Length (m)'); ylabel('Width (m)')
    title('Visualization of Proposed Kite with CoP and CoG')
    hold off
end

