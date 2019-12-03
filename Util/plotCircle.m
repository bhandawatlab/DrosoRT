function h = plotCircle(center,r,N,color)

THETA=linspace(0,2*pi,N);
RHO=ones(1,N)*r;
[X,Y] = pol2cart(THETA,RHO);
X=X+center(1);
Y=Y+center(2);
h=plot(X,Y,'Color',color,'LineStyle','-','Linewidth',1);