function nsv=sva_numsv(data,mod,method)
%        nsv=sva_numsv(data,mod,method)
%      input: data: m-by-n, m voxels, n samples
%             mod: n-by-p, n samples, p covariates
%             method: 'be' or 'leek'
%      output:nsv, number of estimated surrogate variables

if strcmp(method,'be')
    [m,n]=size(data);
    H=mod*((mod'*mod)\mod');
    %or H=mod*inv(mod'*mod)*mod';
    res=data-(H*data')';
    [~,S,~] = svd(res);
    ndf=min(m,n)-ceil(sum(diag(H)));
    S=diag(S);
    dstat=S(1:ndf).^2./sum(S(1:ndf).^2);
    dstat0=zeros(20,ndf);
    for i=1:20
        res0=res(:,randperm(n));
        res0=res0-(H*res0')';
        [~,s0,~]=svd(res0);
        s0=diag(s0);
        dstat0(i,:)=s0(1:ndf).^2./sum(s0(1:ndf).^2);
    end
    psv=ones(1,n);
    for i=1:ndf
        psv(i) = mean(dstat0(:, i) >= dstat(i));
    end
    for i=2:ndf
        psv(i) = max(psv((i - 1)), psv(i));
    end
    nsv = sum(psv <= 0.1);
end

if strcmp(method,'leek')
    dims=size(data);
    a=0:(2/99):2;
    n=floor(dims(1)/10);
    rhat=zeros(100,10);
    P=diag(1*ones(1,dims(2)))-mod*((mod'*mod)\mod');
    for j=1:10
        dats=data(1:(j*n),:);
        [~,D]=eig(dats'*dats);
        ee=sort(diag(D),'descend');
        sigbar = ee(dims(2))/(j * n);
        R = dats * P;
        wm = (1/(j * n)) .* R' * R - P .* sigbar;
        [~,D]=eig(wm);
        ee=sort(diag(D),'descend');
        v = [ones(1, 100), zeros(1, dims(2))];
        [~,ind]=sort([a * (j * n)^(-1/3) * dims(2), ee'], 'descend');
        v=v(ind);
        u = 1:length(v);
        w = 1:100;
        rhat(:, j) = wrev((u(v == 1) - w));
    end
    ss=var(rhat,0,2)';
    [~,bumpstart] = max(ss > (2 * ss(1)));
    [~,start] = max([(1e+05)*ones(1, bumpstart), ss((bumpstart + 1):100)] < 0.5 * ss(1));
    [~,finish]= max(ss .* [zeros(1, start), ones(1, 100 - start)] > ss(1));
    if finish==1
        finish=100;
    end
    x=rhat(start:finish,10);
    tbl=tabulate(x);
    [~,ind]=sort(-tbl(:,2));
    tbl=tbl(ind,:);
    nsv=tbl(1,1);
end



end