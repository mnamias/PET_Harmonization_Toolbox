function [CRCmax, CRCmean, error, fail] = get_CRCs_quick(spheres, mean_bkg, sphere_signal, showfig);
% Name of code: get_CRCs_quick.m
% Version: 1.0 (April 20, 2018)
% Level of code: function. This function is from the PET harmonization toolbox GUI
% Copyright (c) 2018, Mauro Namías mnamias@gmail.com
% All rights reserved.
% This code is intended to accompany the paper:
% Namías et. al
% A novel approach to quantitative harmonization in PET
% PMB (2018)
%
% function [CRCmax, CRCmean, error, fail] = get_CRCs_quick(spheres, mean_bkg, sphere_signal, showfig);
%
% Estimates CRCmax and CRCmean values and the error between the simulated spheres CRC´s and the target CRC values. 
%
% Inputs: 
%   spheres: (Nmax x 6) MATLAB structure with the simulated spheres 
%   sphere_signal: ideal sphere uptake (default: 9.75 for EARL harmonization).
%   mean_bkg: background signal (default: 1)
%   showfig: if equal to 1, will show figures. 
%
% Outputs:
% CRCmax: CRCmax values for the NEMA spheres
% CRCmean: CRCmean values for the NEMA spheres
% error: RMS error between target and measured CRCs
% fail: number of CRCs outside EARL tolerances
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

[a,b] = size(spheres);
    
    
   CRCmax_min = [0.34 0.59 0.73 0.83 0.91 0.95];
   CRCmax_max = [0.57 0.85 1.01 1.09 1.13 1.16];
   CRCmax_target = (CRCmax_min+CRCmax_max)/2;
   
   CRCmean_min = [0.27 0.44 0.57 0.63 0.72 0.76];
   CRCmean_max = [0.43 0.60 0.73 0.78 0.85 0.89];
   CRCmean_target = (CRCmean_min+CRCmean_max)/2;
   
   diameters = [10 13 17 22 28 37];
   

                       for i = 1:a
                           for j = 1:b
                               
                               
                        voi = spheres(i,j).sphere(:);
                        maxi = max(voi(:));
                        SUVmax(i,j) = maxi;
                        thr = (maxi+mean_bkg)/2;
                        
                        mask = voi>thr;
                     
                        SUVmean(i,j) = mean(voi(mask==1));
                        
                        CRCmax(i,j) = SUVmax(i,j)/sphere_signal;
                        
                        CRCmean(i,j) =  SUVmean(i,j)/sphere_signal;
                        error(i,j) = (CRCmax(i,j)-CRCmax_target(j)).^2/CRCmax_target(j)+(CRCmean(i,j)-CRCmean_target(j)).^2/CRCmean_target(j) ;
                       end
                       end
                       
    
                       
                      if(showfig) 
                      figure
                      plot(diameters(1:6),CRCmax_max,'k --')
                      hold on
                      plot(diameters(1:6),CRCmax_min,'k --')
                      
                      CRCmax_avg = mean(CRCmax,1);
                      CRCmean_avg = mean(CRCmean,1);

                      
                      for i = 1:a
                      plot(diameters(1:6),CRCmax(i,:), 'k .')
                      hold on
                      end
                 %     plot(diameters(1:6),CRCmax_avg, 'w +')
                       
                      xlabel('Sphere diameter (mm)')
                      ylabel('CRCmax')
                      axis([4 40 0.1 1.6])
                       
                      %%%%%%%%%%%%%%%%%%%%%%
                      figure
                      
                      h = boxplot(CRCmax, diameters(1:6),'outliersize',1)
                      set(findobj(get(h(1), 'parent'), 'type', 'text'), 'fontsize', 11, 'FontName', 'Century Gothic');
                      
                      hold on
                      plot(1:6,CRCmax_max,'k --')
                      plot(1:6,CRCmax_min,'k --')
                      xlabel('Sphere diameter (mm)')
                      ylabel('CRCmax')
                      axis([0.5 6.5 0.1 1.6])
                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                      
                      figure
                      plot(diameters(1:6),CRCmean_max,'k --')
                      hold on
                      plot(diameters(1:6),CRCmean_min,'k --')
                      
                      for i = 1:a
                      plot(diameters(1:6),CRCmean(i,:), 'k .')
                      hold on
                      end
                 %     plot(diameters(1:6),CRCmean_avg, 'w +')

                       xlabel('Sphere diameter (mm)')
                      ylabel('CRCmean')
                     axis([4 40 0.1 1.2])
                     
                     %%%%%%%%%%%%
                     
                     figure
                      
                    h =   boxplot(CRCmean, diameters(1:6),'outliersize',1)
%                    set(findobj(get(h(1), 'parent'), 'type', 'text'), 'fontsize', 14);
                      set(findobj(get(h(1), 'parent'), 'type', 'text'), 'fontsize', 11, 'FontName', 'Century Gothic');

                      hold on
                      plot(1:6,CRCmean_max,'k --')
                      plot(1:6,CRCmean_min,'k --')
                      xlabel('Sphere diameter (mm)')
                      ylabel('CRCmean')
                      axis([0 7 0.1 1.2])
                      end
                      
                  %     error = sum((CRCmax-CRCmax_target).^2./CRCmax_target)+ sum((CRCmean-CRCmean_target).^2./CRCmean_target);
%                 
                        CRCmax_max = repmat(CRCmax_max,[a 1]);
                        CRCmax_min = repmat(CRCmax_min,[a 1]);
                        CRCmean_max = repmat(CRCmean_max,[a 1]);
                        CRCmean_min = repmat(CRCmean_min,[a 1]);
                        
                        
                         fail1 = CRCmax(CRCmax>CRCmax_max);
                         fail1 = sum(fail1>0);
                         fail2 = CRCmax(CRCmax<CRCmax_min);
                         fail2 = sum(fail2>0);
                         fail3 = CRCmean(CRCmean<CRCmean_min);
                         fail3 = sum(fail3>0);
                         fail4 = CRCmean(CRCmean>CRCmean_max);
                         fail4 = sum(fail4>0);
                         fail = fail1+fail2+fail3+fail4;
                         
                         error = sqrt(sum(error(:).^2));
                         error = error*(1+fail*10); % amplify the error if there are CRC values outside range
                         
                       