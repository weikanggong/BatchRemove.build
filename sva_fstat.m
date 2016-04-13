function [fstats,p]=sva_fstat(data,mod,mod0)
% A function for quickly calculating f statistic and p-values
%
% This function does simple linear algebra to calculate f-statistics
% for each row of a data matrix comparing the nested models
% defined by the design matrices for the alternative (mod) and and null (mod0) cases.
% The columns of mod0 must be a subset of the columns of mod.
% Input:
%    data: m-by-n, m feature and n samples
%    mod: n-by-p, The model matrix being used to fit the data.
%    mod0: n-by-p1, The null model being compared when fitting the data.
% Output:
%    fstats: A vector of F-statistic one for each row of dat.
%        p: p-value.
%
n=size(data,2);
df1=size(mod,2);
df0=size(mod0,2);
Id=eye(n);

resid = data * (Id - mod*((mod'*mod)\mod'));
rss1 = sum(resid.*resid,2);
clear resid

resid0 = data * (Id - mod0*((mod0'*mod0)\mod0'));
rss0 = sum(resid0.*resid0,2);
clear resid0

fstats = ((rss0 - rss1)/(df1-df0))./(rss1/(n-df1));

p=fcdf(fstats,df1-df0,n-df1,'upper');







end
