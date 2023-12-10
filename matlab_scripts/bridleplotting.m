function [] = bridleplotting(L, CoG, CoP,r1v,r2v,BiV)
%BRIDLEPLOTTING will amaze you by plotting your kite's bridle point and
%cords
%   By simply inputting the length, center of gravity, center of pressure,
%   and eithe reaction vectors, this function will quickly and easily plot
%   a detailed side view of your kite. It will include a labeled CoP, CoG,
%   reaction vectors, bridle point, and bridle cord.

%PLOT!
figure();
    plot([0,L],[0,0],'Color','Black') %Kite side view
    axis equal; hold on 
    plot(CoP,0,'.','Color','Red','MarkerSize',20) %CoP
    plot(CoG,0,'.','Color','Green','MarkerSize',20) %CoG
    plot(BiV(1),BiV(2),'.','Color','Blue','MarkerSize',20) %bridle pt
    %plot([BiX,CoP],[BiY,0],'--','Color','Red') %r1
    %plot([BiX,CoG],[BiY,0],'--','Color','Green') %r2
    quiver(BiV(1),BiV(2),r1v(1),r1v(2),'--','Color','Red')
    quiver(BiV(1),BiV(2),r2v(1),r2v(2),'--','Color','Green')
    plot([0,BiV(1)],[0,BiV(2)],'Color','Blue') %bridle cords
    plot([BiV(1),L],[BiV(2),0],'Color','Blue')
    legend('','CoP','CoG','Bridle Point','','','Bridle Cords','','Location','Southeast')
    xlabel('Length (m)'); ylabel('Height (m)')
hold off

end