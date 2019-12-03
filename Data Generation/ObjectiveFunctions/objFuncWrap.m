function [xTmp] = objFuncWrap(obj,gs)
% set up the problem with initial estimates, upper and lower bounds  of 
% parameters
problem = createOptimProblem('fmincon','objective',obj,...
    'x0',rand(6,1).*0.6+0.1,'lb',zeros(6,1)+eps,'ub',ones(6,1));
% run the optimization
xTmp = run(gs,problem);
end