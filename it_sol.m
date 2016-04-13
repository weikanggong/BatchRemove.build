  %       adjust=it.sol(sdat,g_hat,d_hat,g_bar,t2,a,b)
  %input: sdat: n-by-p data matrx: n=number of voxels, p=number of sample in a
  %             single batch
  %       g_hat: n-by-1 vector, prior gamma(normal) of n voxels
  %       d_hat: n-by-1 vector, prior delta(inverse gamma) of n voxels
  %       g_bar: a number:mean of g_hat 
  %       t2: a number: variance of g_hat
  %       a: ¡®lammba¡¯ parameter of inverse gamma distribution of a batch
  %       b: 'theta' parameter of inverse gamma distribution of a batch
  %output: adjust: 2-by-n matrix: posterior estimated gamma star and delta star for a
  %                               batch with n voxels
function adjust=it_sol(sdat,g_hat,d_hat,g_bar,t22,a,b)
 
  conv=0.0001;
  n = sum(~isnan(sdat),2);
  g_old = g_hat;
  d_old = d_hat;
  change = 1;
  count =0;
  %EM
  while change>conv
    %calculate posterior gamma_star of a batch with n voxels
    g_new = (t22*n.*g_hat+d_old.*g_bar)./(t22*n+d_old);
    
    %calculate posterior delta_star of a batch with n voxels
    sum2 = sum((sdat-repmat(g_new,1,size(sdat,2))).^2,2);
    d_new = (0.5*sum2+b)./(n/2+a-1);
    
    change = max(max(abs(g_new-g_old)./g_old,abs(d_new-d_old)./d_old));
    g_old = g_new;
    d_old = d_new;
    count = count+1;
  end
  %disp(['This batch took ', num2str(count), ' iterations until convergence']);
  adjust= cat(2,g_new, d_new)';

end
  



