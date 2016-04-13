function [sv,pprob_gam]=sva_2stepsva(data, mod, n_sv)
%      [sv,nsv,pprob_gam]=sva_TwoStepSVA(dat, mod, n_sv)
%      input: data: m-by-n, m voxels, n samples
%             mod: n-by-p, n samples, p covariates
%      output:sv: n-by-nsv matrix: The estimated surrogate variables, one in each column
%                 pprob_gam: A vector of the posterior probabilities each gene is affected by heterogeneity
if fix(n_sv)==n_sv && length(n_sv)==1
    [m,n]=size(data);
    H=mod*((mod'*mod)\mod');
    res=data-(H*data')';
    [~,~,V] = svd(res);
    res_sv = V(:, 1:n_sv);
    use_var=zeros(m,n_sv);
    pp=zeros(m,n_sv);
    for i=1:n_sv
        modx = [ones(n,1),res_sv(:,i)];
        modx0 = ones(n,1);
        [~,pp(:,i)]=sva_fstat(data,modx,modx0);
        use_var(:,i)=mafdr(pp(:,i))<0.01;
    end
    for i=size(use_var,2):1
        if sum(use_var(:,i))<n
            use_var(:,i)=[];
            n_sv=n_sv-1;
            if n_sv<=0
                break;
            end
        end
    end
    if n_sv>0
        sv=zeros(n,n_sv);
        data1=mean_center(data')';
        for i=1:n_sv
            [~,~,V]=svd(data1(logical(use_var(:,i)),:));
            maxcor=0;
            for j=1:(n-1)
                rr=abs(corr(V(:,j),res_sv(:,i)));
                if rr>maxcor
                    maxcor=rr;
                    sv(:,i)=V(:,j);
                end
            end
        end
        pprob_gam=use_var*ones(n_sv,1)>0;
    end
    
end
end