clear, close all
nlobj_new_3x = nlmpc(4,1,'MV',4,'MD',[1 2 3]);

Ts=900;
nlobj_new_3x.Ts=Ts;
nlobj_new_3x.PredictionHorizon=10;
nlobj_new_3x.ControlHorizon=1;
nlobj_new_3x.MV(1).ScaleFactor=850;
nlobj_new_3x.MV(1).Min=0;
nlobj_new_3x.MV(1).Max=850;

nlobj_new_3x.Model.StateFcn='ZonePcmStateFcnCT_3x';
nlobj_new_3x.Model.OutputFcn='ZonePcmOutputFcn';


nlobj_new_3x.Optimization.CustomCostFcn = 'ZonePcmCostFunction';
nlobj_new_3x.Optimization.ReplaceStandardCost = true;
% nlobj_new.Optimization.SolverOptions.Algorithm='interior-point';
% nlobj_new.Optimization.SolverOptions.Algorithm='active-set';
% nlobj_new.Optimization.SolverOptions.Algorithm='sqp-legacy';
% nlobj_new.Optimization.UseSuboptimalSolution=true;







