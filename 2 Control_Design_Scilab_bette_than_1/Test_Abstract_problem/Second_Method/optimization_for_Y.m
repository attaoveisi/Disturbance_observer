function [x,fval,exitflag,output,population,score] = optimization_for_Y(nvars,lb,ub,PopulationSize_Data,EliteCount_Data,CrossoverFraction_Data,MigrationFraction_Data,Generations_Data,StallGenLimit_Data,TolFun_Data,TolCon_Data)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = gaoptimset;
%% Modify options setting
options = gaoptimset(options,'PopulationSize', PopulationSize_Data);
options = gaoptimset(options,'EliteCount', EliteCount_Data);
options = gaoptimset(options,'CrossoverFraction', CrossoverFraction_Data);
options = gaoptimset(options,'MigrationFraction', MigrationFraction_Data);
options = gaoptimset(options,'Generations', Generations_Data);
options = gaoptimset(options,'StallGenLimit', StallGenLimit_Data);
options = gaoptimset(options,'TolFun', TolFun_Data);
options = gaoptimset(options,'TolCon', TolCon_Data);
options = gaoptimset(options,'MutationFcn', @mutationadaptfeasible);
options = gaoptimset(options,'Display', 'off');
options = gaoptimset(options,'PlotFcns', { @gaplotbestf });
[x,fval,exitflag,output,population,score] = ...
ga(@finding_Y,nvars,[],[],[],[],lb,ub,[],[],options);
