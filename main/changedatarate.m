%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Copyright (C) 2011- Flávio Luiz C. Ribeiro
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tnew Xnew] = changedatarate(t,X,dt)
    %modificacao da taxa de amostragem:
    tfinal = t(size(t,1));
    %ts=timeseries(X,t,'Name','speed');
    %ts=transpose(ts);
    %res_ts=resample(ts,0:dt:tfinal);
    tnew = 0:dt:tfinal;
    Xnew = interp1(t,X,tnew);
    %tnew = res_ts.time;
    %Xnew = reshape(res_ts.Data, size(X,2),size(tnew,1))';
end