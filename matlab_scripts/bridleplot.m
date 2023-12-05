function [] = bridleplot(L,x1,x2,CoP,CoG,r1)
%BRIDLEPLOT plots the bridle points of your kite with ease!!
%   Given inputs of the length of your (simple, diamond-shaped) kite, the
%   offsets of your bridle points (x1,x2), and the primary reaction vector 
%   (r1) magnitude, this function will amaze you with its ability to
%   quickly and easily plot bridle points, cords, and the CoP and CoG of
%   your kite (includes color coding!!!). Primarily created to make my 
%   main script cleaner.

figure()
    plot([0,L],[0,0],"Color","Black"); axis equal; hold on
    plot((L-x1),0,'.','MarkerSize',20,'Color','Blue') %near bridal point
    plot(x1,0,'.','MarkerSize',20,'Color','Blue') %far bridal point
    plot([CoP,L-x1],[-r1,0],'--','Color','Black') %bridal cord 1
    plot([x2, CoP],[0,-r1],'--','Color','Black') %bridal cord 2
    plot(CoG,0,'.','MarkerSize',20,'Color','Red') %CoP
    plot(CoP,0,'.','MarkerSize',20,'Color','Green') %CoG
    plot([CoP,CoP],[-r1,0],'--','Color','Green') %r1
    plot([CoP,CoG],[-r1,0],'--','color','Red') %r2
    legend('','Bridal Point','','Bridal Cords','','CoP','CoG','','')
    xlabel('Length(m)'); ylabel('Height (m)')
    title('Side View of Kite with Bridle Points and Cords Visualized')
hold off

end