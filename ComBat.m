function bayesdata=ComBat(dat, batch, mod , par_prior) 
%       bayesdata=ComBat(dat, batch, mod , par_prior) 
%input: dat: m-by-n matrix: m voxels , n samples
%       batch: n-by-1 vector:  n sample batch variable
%       mod:  n-by-p matrix:  n samples, p known covariates
%       par_prior : 1 =parametric prior  0=non-parametrc prior (support parallel computing)
%output: bayesdata: adjusted data
  if size(batch,2) > 1 
      error('This version of ComBat only allows one batch variable');
  end
      batch1 = unique(batch);
      batchmod=dummyvar(batch);
      ind=(mean(batchmod==0)~=1);
      batchmod=batchmod(:,ind);
      
      n_batch =length(batch1);
      disp(['Found ', num2str(n_batch), ' batches']);
      n_batches=sum(batchmod);
      n_array=sum(n_batches);
      
      design=cat(2,batchmod,mod);
      
      
      
      n_covariates=size(mod,2);
      disp(['Adjusting for ', num2str(n_covariates), ' covariate(s) or covariate level(s)']);
      
    if (rank(design) < size(design,2)) 
        if size(design,2) == (n_batch + 1) 
             error('The covariate is confounded with batch! Remove the covariate and rerun ComBat');
        end
        if size(design,2) > (n_batch + 1)
           if rank(mod) < size(mod,2) 
              error('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded');
        
           else 
              error('At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat');
           end
        end
    end
      
    NAs=sum(sum(isnan(dat)));  
    if (NAs>0) 
       error(['Found ',num2str(NAs),' Missing Data Values']);
    end  
      
    disp('Standardizing Data across voxels/genes/FCs'); 
    
    if NAs==0
        B_hat=inv(design'*design)*design'*dat';
    end
    
    grand_mean= (n_batches./n_array) * B_hat(1:n_batch,:);
      
    if NAs==0
         var_pooled = ((dat - (design * B_hat)').^2) * 1/n_array*ones(n_array,1);
                                                          
    end
    
    stand_mean = grand_mean' * ones(1, n_array);
    
    if isempty(design)==0
         tmp=design;
         tmp(:, 1:n_batch)= 0;
         stand_mean = stand_mean + (tmp * B_hat)';
        
    end
    
    s_data=(dat - stand_mean)./(sqrt(var_pooled)*ones(1, n_array));
                                                          
    disp('Fitting L/S model and finding priors');
    
    batch_design =design(:, 1:n_batch);
    if NAs==0 
       gamma_hat = inv(batch_design'*batch_design)*batch_design'*s_data'; 
    end
    
    delta_hat=[];
    for i=1:n_batch
        ind1=batch==batch1(i);
        delta_hat= cat(1,delta_hat, var(s_data(:, ind1)'));
    end
    gamma_bar = mean(gamma_hat,2);
    t2 = var(gamma_hat');
    
    m1=mean(delta_hat,2);
    v1=var(delta_hat');
    a_prior=(2*v1'+m1.^2)./v1';
    b_prior=(m1.*v1'+m1.^3)./v1';
    
   gamma_star=zeros(n_batch,size(dat,1));
   delta_star=zeros(n_batch,size(dat,1));
   if par_prior==1
      disp('Finding parametric adjustments');
      for i=1:n_batch
         ind1=batch==batch1(i);
         %caluculate posterior parameter gamma_star and delta_star
         temp =it_sol(s_data(:,ind1),gamma_hat(i,:)',delta_hat(i, :)',gamma_bar(i),t2(i),a_prior(i),b_prior(i)); 
         gamma_star(i,:) = temp(1,:);
         delta_star(i,:) = temp(2,:);
      end
   end
   
   if par_prior==0
       disp('Finding non-parametric adjustments');
       for i=1:n_batch
          disp(['calculating the ',num2str(i),'th batches']);
          ind1=batch==batch1(i);
          temp = int_eprior(s_data(:,ind1),gamma_hat(i,:)',delta_hat(i, :)');
          gamma_star(i,:) = temp(1,:);
          delta_star(i,:) = temp(2,:);
          
       end
   end
   
   disp('Adjusting the Data');
   bayesdata = s_data;
   j=1;
   
   for i=1:n_batch
       ind1=batch==batch1(i);
       bayesdata(:, ind1) = (bayesdata(:, ind1) - (batch_design(ind1,:) * gamma_star)')./(sqrt(delta_star(j,:))' * ones(1,n_batches(j)));
       j = j + 1 ;                                                                                                                                                    
   end
   
   bayesdata = (bayesdata .* (sqrt(var_pooled) * ones(1, n_array))) + stand_mean;
                                                        
    
end

  
  
  