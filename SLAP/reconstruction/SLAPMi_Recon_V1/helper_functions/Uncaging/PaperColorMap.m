function c = PaperColorMap()
% c = hot(1000);
c = jet(1000);
%c(1:125,3) = (c(1:125,3)-c(1,3))./(c(125,3)-c(1,3));
c = [ zeros(125,1) zeros(125,1) linspace(0,0.5, 125)' ;c];
% c(:,3) = c(:,3)/10;
% figure
% hold on;
% plot(c(:,1),'r','linewidth',3);
% plot(c(:,2),'g','linewidth',3);
% plot(c(:,3),'b','linewidth',3);
% c = [0 0 0 ; 1 0 0; 1 1 0; 1 1 1];
% c = interp1(0:3,c, linspace(0,3,1000));
end

