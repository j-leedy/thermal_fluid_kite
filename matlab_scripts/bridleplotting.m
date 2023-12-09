function [BiP] = bridleplotting(L, CoG, CoP,r1v)
%BRIDLEPLOTTING will amaze you by plotting your kite's bridle point and
%cords
%   By simply inputting the length, center of gravity, center of pressure,
%   and eithe reaction vectors, this function will quickly and easily plot
%   a detailed side view of your kite. It will include a labeled CoP, CoG,
%   reaction vectors, bridle point, and bridle cord.

BiP = [CoG-r1v(1);0-r1v(2);0]; %location of bridal point with respect to kite
BiX = BiP(1);
BiY = BiP(2);

%PLOT!
figure();
    plot([0,L],[0,0],'Color','Black') %Kite side view
    axis equal; hold on 
    plot(CoP,0,'.','Color','Red','MarkerSize',20) %CoP
    plot(CoG,0,'.','Color','Green','MarkerSize',20) %CoG
    plot(BiX,BiY,'.','Color','Blue','MarkerSize',20) %bridle pt
    plot([BiX,CoP],[BiY,0],'--','Color','Red') %r1
    plot([BiX,CoG],[BiY,0],'--','Color','Green') %r2
    plot([0,BiX],[0,BiY],'Color','Blue') %bridle cords
    plot([BiX,L],[BiY,0],'Color','Blue')
    legend('','CoP','CoG','Bridle Point','','','Bridle Cords','','Location','Southeast')
    xlabel('Length (m)'); ylabel('Height (m)')
hold off

end