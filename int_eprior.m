function  adjust=int_eprior(sdat,g_hat,d_hat)
  %       adjust=int_eprior(sdat,g_hat,d_hat)
  %input: sdat: n-by-p data matrx: n=number of voxels, p=number of sample in a
  %             single batch
  %       g_hat: n-by-1 vector, prior gamma  of n voxels
  %       d_hat: n-by-1 vector, prior delta  of n voxels
  %output: adjust: 2-by-n matrix: posterior estimated gamma star and delta star for a
  %                               batch with n voxels
       
       r = size(sdat,1);
       g_star=zeros(r,1);
       d_star=zeros(r,1);
       for i=1:r
           g = g_hat; g(i)=[];
           d = d_hat; d(i)=[];		
           x = sdat(i,~isnan(sdat(i,:)));
           n = length(x);
           j = ones(n,1);
           dat = repmat(x,length(g),1);
           resid2 = (dat-repmat(g,1,size(dat,2))).^2;
           sum2 = resid2*j;
           LH = 1./(2.*pi.*d).^(n/2).*exp(-sum2./(2*d));
           LH(isnan(LH))=0;
           g_star(i) = sum(g.*LH)./sum(LH);
           d_star(i) = sum(d.*LH)./sum(LH);
 
       end
       
       adjust = cat(2,g_star,d_star)';



end