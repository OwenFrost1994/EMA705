function []=meshplotfunc(nElements,nNodes,nNoEl,connArray,coorNoUD,order,W,H,W1,plotn)
figure(plotn)
% rectangle('Position',[(W-W1)/2,0,W1,H],'LineWidth',1,'linestyle','-','edgecolor','b');
% hold on;
% x=[(W-W1)/2 (W+W1)/2 (W+W1)/2 (W-W1)/2];
% y=[0 0 H H];
% fill(x,y,'y')
% hold on;
for i=1:1:nElements
    if order==1
        for j=1:1:nNoEl
            if j<4
                A=j;
                B=j+1;
            else
                A=j;
                B=j-3;
            end
            plot([coorNoUD(connArray(i,A),1),coorNoUD(connArray(i,B),1)],[coorNoUD(connArray(i,A),2),coorNoUD(connArray(i,B),2)],'b-','Markersize',8,'Markerface','white','linewidth',2.0);
            hold on;
        end
    else
        for j=1:1:nNoEl
            if j<=4
                A=j;
                B=j+4;
            else
                if j==8
                    A=j;
                    B=j-7;
                else
                    A=j;
                    B=j-3;
                end
            end
            plot([coorNoUD(connArray(i,A),1),coorNoUD(connArray(i,B),1)],[coorNoUD(connArray(i,A),2),coorNoUD(connArray(i,B),2)],'b-','Markersize',8,'Markerface','white','linewidth',2.0);
            hold on;
        end
    end
end
for i=1:1:nNodes
    plot(coorNoUD(i,1),coorNoUD(i,2),'ro','Markersize',8,'Markerface','white','linewidth',3.0);
    hold on;
end
grid on;
xlabel('x(m)')
ylabel('y(m)')
axis([-0.01 0.04 0 0.04])
set(gca, 'FontName','Times New Roman','FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Undeformed Mesh');
end