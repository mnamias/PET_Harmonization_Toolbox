function [CRCmax, CRCmean, SUVmax, SUVmean, CV, error, fail] = get_CRCs_islands_EARL(imagen, data, maskx1, maskx2, masky1, masky2, I, bkg_voi, EARL_mode, half_life);

    t1 = data(1).sphere_reference_time;
    t2 = data(1).img_time;
    h1 = str2num(t1(1:2));
    h2 = str2num(t2(1:2));
    m1 = str2num(t1(3:4));
    m2 = str2num(t2(3:4));
    s1 = str2num(t1(5:6));
    s2 = str2num(t2(5:6));
    
    dt = (h2-h1)*60+(m2-m1)+(s2-s1)/60; %% Elapsed time in minutes
    disp(['dt: ' num2str(dt)  ' mins']);
    
    Fd = exp(-log(2)*dt/half_life); %%  decay factor
    %Fd = exp(-log(2)*dt/67.629); F18 decay factor
    
     %data(1).sphere_A
%data(1).sphere_vol
    true_sphere_A = data(1).sphere_A/data(1).sphere_vol * Fd *1e6; %% expected sphere activity concentration, in Bq/ml
    true_sphere_A/1000
 %   pause(1)
    
    if(EARL_mode == 1)
   CRCmax_min = [0.34 0.59 0.73 0.83 0.91 0.95];
   CRCmax_max = [0.57 0.85 1.01 1.09 1.13 1.16];
   CRCmax_target = (CRCmax_min+CRCmax_max)/2;
   
   CRCmean_min = [0.27 0.44 0.57 0.63 0.72 0.76];
   CRCmean_max = [0.43 0.60 0.73 0.78 0.85 0.89];
   CRCmean_target = (CRCmean_min+CRCmean_max)/2;
   CVmax = 0.15;
    elseif (EARL_mode == 2)
        
   CRCmax_min = [0.52 0.85 1.00 1.01 1.01 1.05];
   CRCmax_max = [0.88 1.22 1.38 1.32 1.26 1.29];
   CRCmax_target = (CRCmax_min+CRCmax_max)/2;
   
   CRCmean_min = [0.39 0.63 0.76 0.80 0.82 0.85];
   CRCmean_max = [0.61 0.86 0.97 0.99 0.97 1.00];
   CRCmean_target = (CRCmean_min+CRCmean_max)/2;
   CVmax = 0.15;
        
    elseif (EARL_mode == 3)
        
  CRCmax_min = [0.34 0.59 0.73 0.83 0.91 0.95];
   CRCmax_max = [0.57 1.85 1.01 1.09 1.13 1.16];
   CRCmax_target = [0.5 1.0 1.0 1.0 1.0 1.0];
   
   CRCmean_min = [0.27 0.44 0.57 0.63 0.72 0.76];
   CRCmean_max = [0.43 0.60 0.73 0.78 0.85 0.89];
   CRCmean_target = (CRCmean_min+CRCmean_max)/2;
   
    end
    
   
   diameters = [10 13 17 22 28 37];
   

                       % title('Estimating CRCs...');

                       bkg = mean(bkg_voi(:));
                       CV = std(bkg_voi(:))/mean(bkg_voi(:));

                       for i = 1:6
%                            size(imagen)
%                            masky1(i)
%                            masky2(i)
%                            maskx1(i)
%                            maskx1(i)
%                            I-21
%                            I+21
                        voi =  imagen(masky1(i):masky2(i),maskx1(i):maskx2(i),I-15:I+15);
                        maxi = max(voi(:));
                        SUVmax(i) = maxi;
                        thr = (maxi+bkg)/2;
                        
                        mask = voi>thr;
                     %   mask = bwareaopen(mask,4);
                        SUVmean(i) = mean(voi(mask==1));
                        
                        CRCmax(i) = SUVmax(i)/true_sphere_A;
                        
                        CRCmean(i) =  SUVmean(i)/true_sphere_A;
                       end

%                        figure
%                        plot(diameters,CRCmax,'k +', diameters,CRCmax_min,'g-.' , diameters,CRCmax_max,'g-.', diameters, CRCmax_target, 'r- ' )
%                        xlim([5 40]);
%                        xlabel('Sphere Diameter [mm]');
%                        ylabel('CRCmax');
% %                       title(['Ratio: ' num2str(ratio)])
                       
%                        figure
%                        plot(diameters,CRCmean,'k +', diameters,CRCmean_min,'g-.' , diameters,CRCmean_max,'g-.', diameters, CRCmean_target, 'r- ' )
%                        xlim([5 40]);
%                        xlabel('Sphere Diameter [mm]');
%                        ylabel('CRCmean');
             %          title(['Ratio: ' num2str(ratio)])
                     
                       error = sum((CRCmax-CRCmax_target).^2./CRCmax_target)+ sum((CRCmean-CRCmean_target).^2./CRCmean_target);
%                        dummy = CRCmax(CRCmax>CRCmax_target);
%                        error = error + sum(dummy(dummy>0))*2;
%                        dummy = CRCmax(CRCmax<CRCmax_target);
%                        error = error + sum(dummy(dummy>0))*2;
%                        dummy = CRCmean(CRCmean<CRCmean_target);
%                        error = error + sum(dummy(dummy>0))*2;
%                        dummy = CRCmean(CRCmean>CRCmean_target);
%                        error = error + sum(dummy(dummy>0))*2;
                       if(CV > CVmax)
                           error = error*10;
                       end
                       
                         fail1 = CRCmax(CRCmax>CRCmax_max);
                         fail1 = sum(fail1>0);
                         fail2 = CRCmax(CRCmax<CRCmax_min);
                         fail2 = sum(fail2>0);
                         fail3 = CRCmean(CRCmean<CRCmean_min);
                         fail3 = sum(fail3>0);
                         fail4 = CRCmean(CRCmean>CRCmean_max);
                         fail4 = sum(fail4>0);
                         fail = fail1+fail2+fail3+fail4;
                         
                         %error = error*(1+fail*10)
                         error = error*(1+fail*0.5)
                         
                         CRCmax = CRCmax';
                         CRCmean = CRCmean';
                         SUVmax = SUVmax';
                         SUVmean = SUVmean';