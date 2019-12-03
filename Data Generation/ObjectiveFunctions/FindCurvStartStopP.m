function [Cstarts,Cstops,CWST] = FindCurvStartStopP(curv,param)
%This function defines the start and stops of fly turns.
%
%9/7/2018 define start/ stops as periods of near 0 curvature with small
%changes in curvature as well (diff)
%9/12 Instead of searching for straight walks will look for a way to
%distinguish between curevd walks and sharp turns. trying first by looking
%for periods of high curvature and high change ing curvature and labeling
%as sharp turns
% curv = curv.Curv.raw.curv;
% Cndx = C.Curv.ndx;
dcurv = diff(curv);

%close
clear Sndx

cwgt = param(1);
dwgt = param(2);
sm = 1;%ceil(param(3)*16);
Idsc = param(3)*400;
Idsd = param(4)*400;
% OdsT = param(6);
% OdsW = param(7);
% fsm = round(param(8)*8);
OdsT = param(5);
OdsW = param(6);
fsm = 5;

cFlag = zeros(1,length(curv));
dFlag = zeros(1,length(curv));
%find periods of high curvature and change in curvature to locate sharp
%turns
cFlag(abs(curv)>pi/Idsc) = 1;
dFlag(abs(dcurv)>pi/Idsd) = 1;

TFlag = cFlag*cwgt + dFlag*dwgt;
sTFlag = smooth(TFlag,sm);

Tndx = find(sTFlag>(1+OdsT)*mean(sTFlag));
Wndx = find(sTFlag<(1-OdsW)*mean(sTFlag));

CWST = zeros(1,length(curv));
N1 = zeros(1,length(curv));
N2 = zeros(1,length(curv));

CWST(Tndx)=1;
CWST(Wndx)=-1;

% N1(sTFlag>=(1-OdsW))=1;  
% N2(sTFlag<=(1+OdsT))=1;
% 
% N = and(N1,N2);
Nndx = find(CWST==0);

if ~isempty(Tndx) && ~isempty(Wndx)
    for c = Nndx
        if meanmin(abs(c-Tndx),fsm)<meanmin(abs(c-Wndx),fsm)
            CWST(c) = 1;
        else
            CWST(c) = -1;
        end
    end

    Cstarts = find(diff(CWST)>0)+1;
    Cstops = find(diff(CWST)<0);
    if min(Cstops)<min(Cstarts)
        Cstarts=[1 Cstarts];
    end
    if length(Cstarts)>length(Cstops)
        Cstarts=Cstarts(1:end-1);
    end
else
    Cstarts = []; Cstops = [];
end

% plot(s.Center.x,s.Center.y)
% hold on
% plot(s.Center.x(CWST<0),s.Center.y(CWST<0),'*g')
% plot(s.Center.x(CWST>0),s.Center.y(CWST>0),'*r')
% % plot(s.Center.x(Cstarts),s.Center.y(Cstarts),'*bl')
% % plot(s.Center.x(Cstops),s.Center.y(Cstops),'*c')
% legend('Path','Straight','Turn')%,'Starts','Stops')
% % plot(s.Center.xUS(ndx+1),s.Center.yUS(ndx+1),'ob')   
% hold off
end