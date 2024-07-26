%% Plot MLT Anomaly versus the Decomposition of ML Warming from Surface Heat Flux Anomaly

%  1981-2023:
%  Monthly MLT from IAP and SST from ERA5 (as a proxy of MLT) 
%  Interpolate monthly heat flux data to 1 degree grids
%  Monthly heat flux climatology during 1981-2010



%% ########################################################################
%  ########################################################################
%% 1. Monthly MLT from IAP for 1981-2023
%  Monthly MLT 1981-2022 
clc;clear
addpath /srv/ccrc/OceanME/z5195509/MATLAB_Functions
time_ann=(1981:2022)';

% #########################################################################
count_yr=1;
for year=1981:2022
    % #####################################################################
    % Mixed-Layer Temp from IAP
    % Monthly MLD from IAP
    disp(['   MLT from IAPv3, Year#',num2str(year)])
      cd('/srv/ccrc/OceanME/z5195509/Data/IAP_Global_Ocean_Gridded_Product/Global_GSW_Monthly_V2_Double/MLD')
      % dRh0=0.125 ########################################################
      load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld','lon_IAP','lat_IAP')
      lat_IAP=lat_IAP(:,1);
      mld0(:,:,1:12)=squeeze(mld(:,:,1,1:12));
      clear mld
      
      lon_IAP0(1:340,1)=lon_IAP(21:360,1);
      lon_IAP0(341:360,1)=lon_IAP(1:20,1)+360; clear lon_IAP
      lon_IAP=lon_IAP0; clear lon_IAP0
      
      mld00(1:340,:,:)=mld0(21:360,:,:);
      mld00(341:360,:,:)=mld0(1:20,:,:); clear mld0
      mld_monthly=mld00; clear mld00
      % figure;imagesc(mld_monthly(:,:,1))
      % figure;plot(1:12,squeeze(mld_monthly(240,35,1:12)))
    % #####################################################################
      
      
    % #####################################################################
    % Monthly MLT from IAP
      cd('/srv/ccrc/OceanME/z5195509/Data/IAP_Global_Ocean_Gridded_Product/Global_GSW_Monthly_V2_Double/')
      load(['CT_depth_Monthly_gswTSinterp_Glob_10itvl_IAP_C17_',num2str(year),'.mat'],'CT_depth','depth10')
      CT_depth0(1:340,1:180,1:196,1:12)=CT_depth(21:360,1:180,1:196,1:12);
      CT_depth0(341:360,1:180,1:196,1:12)=CT_depth(1:20,1:180,1:196,1:12);
      clear CT_depth
      depth10=depth10(1:196,1);

      depth_MML=nan(length(lon_IAP),length(lat_IAP),length(depth10),12);
      for month=1:12
          for i=1:length(lon_IAP)
              for j=1:length(lat_IAP)
                  k_MML=find(abs(depth10(:,1)-mld_monthly(i,j,month))==min(abs(depth10(:,1)-mld_monthly(i,j,month))));
                  depth_MML(i,j,1:k_MML,month)=1;
                  clear k_MML
              end
          end
      end
      clear i j month
      % mld_10m=squeeze(nansum(depth_MML(:,:,2:end,:),3)).*10+squeeze(nansum(depth_MML(:,:,1,:),3)).*5;
      % figure;imagesc(mld_10m(:,:,1));load color2; colormap(gca,color2); caxis([0 600]);
      % figure;imagesc(mld0(:,:,1));load color2; colormap(gca,color2); caxis([0 600]);
      
      CT_depth0(isnan(depth_MML))=NaN;
      mlt_monthly(:,:,:)=squeeze(nanmean(CT_depth0,3));
      % figure;imagesc(MLT(:,:,1,2));
      clear CT_depth0 
      clear depth_MML depth10 mld_monthly mld_10m
    % #####################################################################

    
      cd('/srv/ccrc/OceanME/z5195509/Temporary_Cal')
      save(['mlt_monthly_Rh0_0125_GSW_IAP_C17_',num2str(year),'.mat'],...
            'mlt_monthly','lon_IAP','lat_IAP')
      clear  mlt_monthly
      count_yr=count_yr+1;
end
% #########################################################################
% #########################################################################


      

% #########################################################################  
% #########################################################################
% Monthly MLT 2023 
clc;clear

% #########################################################################
    % 2.2 Mixed-Layer Temp from IAP in 2023
    % 2023 Monthly MLD - from IAP
    disp(['MLD, MLT, and ML_dens from IAPv3...'])
      cd('/Users/z5195509/Documents/Data/IAP_Ocean/Global_GSW_Monthly/MLD')
      % dRh0=0.125 ########################################################
      load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_2023.mat'],'mld','lon_IAP','lat_IAP')
      mld(:,:,:,10:12)=NaN;
      % ###################################################################
      
      lat_IAP=lat_IAP(:,1);
      mld0(:,:,:)=squeeze(mld(:,:,1,1:12));
      clear mld
      
      lon_IAP0(1:340,1)=lon_IAP(21:360,1);
      lon_IAP0(341:360,1)=lon_IAP(1:20,1)+360; clear lon_IAP
      lon_IAP=lon_IAP0; clear lon_IAP0
      
      mld00(1:340,:,:)=mld0(21:360,:,:);
      mld00(341:360,:,:)=mld0(1:20,:,:); clear mld0
      % mld00(isnan(basin_dep_mon))=NaN;
      mld0=mld00; clear mld00
      mld_2023(:,:,:)=mld0;


    % 2023 Monthly MLT - from IAP
      cd('/Users/z5195509/Documents/Data/IAP_Ocean/Global_GSW_Monthly/')
      load(['CT_depth_Monthly_gswTSinterp_Glob_10itvl_IAP_C17_2023.mat'],'CT_depth','depth10')
      CT_depth0(1:340,1:180,1:196,1:12)=CT_depth(21:360,1:180,1:196,1:12);
      CT_depth0(341:360,1:180,1:196,1:12)=CT_depth(1:20,1:180,1:196,1:12);
      clear CT_depth
      depth10=depth10(1:196,1);

      depth_MML=nan(length(lon_IAP),length(lat_IAP),length(depth10),12);
      for month=1:12
          for i=1:length(lon_IAP)
              for j=1:length(lat_IAP)
                  k_MML=find(abs(depth10(:,1)-mld0(i,j,month))==min(abs(depth10(:,1)-mld0(i,j,month))));
                  depth_MML(i,j,1:k_MML,month)=1;
                  clear k_MML
              end
          end
      end
      clear i j month
      % mld_10m=squeeze(nansum(depth_MML(:,:,2:end,:),3)).*10+squeeze(nansum(depth_MML(:,:,1,:),3)).*5;
      % figure;imagesc(mld_10m(:,:,1));load color2; colormap(gca,color2); caxis([0 600]);
      % figure;imagesc(mld0(:,:,1));load color2; colormap(gca,color2); caxis([0 600]);
      
      CT_depth0(isnan(depth_MML))=NaN;
      MLT(:,:,:)=squeeze(nanmean(CT_depth0,3));
      % figure;imagesc(MLT(:,:,1,2));
      clear CT_depth0 
      clear depth_MML depth10 mld0 mld_10m
      
      mlt_monthly(:,:,1:9)=MLT(:,:,1:9); clear MLT

      
% #########################################################################
    % 2.3 Mixed-Layer Temp from IAPv4 in 2023
    % 2023 Monthly MLD - from IAPv4
    disp(['MLD, MLT, and ML_dens from IAPv4...'])
      cd('/Users/z5195509/Documents/Data/IAP_Ocean/Global_GSW_Monthly/MLD')
      % dRh0=0.125 ########################################################
      load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAPv4_Temp_2023.mat'],'mld','lon_IAP','lat_IAP')
      % ###################################################################
      
      lat_IAP=lat_IAP(:,1);
      mld0(:,:,:)=squeeze(mld(:,:,1,1:12));
      clear mld
      
      lon_IAP0(1:340,1)=lon_IAP(21:360,1);
      lon_IAP0(341:360,1)=lon_IAP(1:20,1)+360; clear lon_IAP
      lon_IAP=lon_IAP0; clear lon_IAP0
      
      mld00(1:340,:,:)=mld0(21:360,:,:);
      mld00(341:360,:,:)=mld0(1:20,:,:); clear mld0
      % mld00(isnan(basin_dep_mon))=NaN;
      mld0=mld00; clear mld00
      mld_2023(:,:,:,3)=mld0;
      
      
    % 2023 Monthly MLT - from IAPv4
      cd('/Users/z5195509/Documents/Data/IAP_Ocean/Global_GSW_Monthly/')
      load(['CT_depth_Monthly_gswTSinterp_Glob_10itvl_IAPv4_Temp_2023.mat'],'CT_depth','depth10')
      CT_depth0(1:340,1:180,1:196,1:12)=CT_depth(21:360,1:180,1:196,1:12);
      CT_depth0(341:360,1:180,1:196,1:12)=CT_depth(1:20,1:180,1:196,1:12);
      clear CT_depth
      depth10=depth10(1:196,1);

      depth_MML=nan(length(lon_IAP),length(lat_IAP),length(depth10),12);
      for month=1:12
          for i=1:length(lon_IAP)
              for j=1:length(lat_IAP)
                  k_MML=find(abs(depth10(:,1)-mld0(i,j,month))==min(abs(depth10(:,1)-mld0(i,j,month))));
                  depth_MML(i,j,1:k_MML,month)=1;
                  clear k_MML
              end
          end
      end
      clear i j month
      % mld_10m=squeeze(nansum(depth_MML(:,:,2:end,:),3)).*10+squeeze(nansum(depth_MML(:,:,1,:),3)).*5;
      % figure;imagesc(mld_10m(:,:,1));load color2; colormap(gca,color2); caxis([0 600]);
      % figure;imagesc(mld0(:,:,1));load color2; colormap(gca,color2); caxis([0 600]);
      
      CT_depth0(isnan(depth_MML))=NaN;
      MLT(:,:,:)=squeeze(nanmean(CT_depth0,3));
      % figure;imagesc(MLT(:,:,1,3));
      clear CT_depth0 
      clear depth_MML depth10 mld0 mld_10m
      
      mlt_monthly(:,:,10:12)=MLT(:,:,10:12); clear MLT
      
      cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')
      save(['mlt_monthly_Rh0_0125_GSW_IAP_C17_2023.mat'],...
            'mlt_monthly','lon_IAP','lat_IAP')
      clear  mlt_monthly

% #########################################################################
% #########################################################################



%% ########################################################################
%% ########################################################################
%% 2. Interpolating monthly heat flux to 1 degree grids
clc;clear   
% #########################################################################
for year=1981:2022
    disp([' Heat flux into 1deg Year#',num2str(year)])
    % #####################################################################
    % Monthly MLD from IAP
    disp(['   MLD from IAPv3, Year#',num2str(year)])
      cd('/Users/z5195509/Documents/Data/IAP_Ocean/Global_GSW_Monthly/MLD')
      % dRh0=0.125 ########################################################
      load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld','lon_IAP','lat_IAP')
      lat_IAP=lat_IAP(:,1);
      mld0(:,:,1:12)=squeeze(mld(:,:,1,1:12));
      clear mld
      
      lon_IAP0(1:340,1)=lon_IAP(21:360,1);
      lon_IAP0(341:360,1)=lon_IAP(1:20,1)+360; clear lon_IAP
      lon_IAP=lon_IAP0; clear lon_IAP0
      
      mld00(1:340,:,:)=mld0(21:360,:,:);
      mld00(341:360,:,:)=mld0(1:20,:,:); clear mld0
      mld_monthly=mld00; clear mld00
      % figure;imagesc(mld_monthly(:,:,1))
      % figure;plot(1:12,squeeze(mld_monthly(240,35,1:12)))
    % #####################################################################

    
    % #####################################################################
    % Monthly heat flux terms
      cd('/Users/z5195509/Documents/Data/ERA-5/Monthly_Averaged_Reanalysis_on_Single_Levels/Heat_Flux/')
      load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'.mat'],...
            'latent','longwave','sensible','shortwave','lon','lat')   
      % 20E-380E
      lon0(1:1360,1)=lon(81:1440,1);
      lon0(1361:1440,1)=lon(1:80,1)+360; clear lon
      lon_ERA5=lon0; clear lon0
      lat_ERA5=lat;  clear lat
      
      heatflux0(1:1360,:,:)=latent(81:1440,:,:);
      heatflux0(1361:1440,:,:)=latent(1:80,:,:); clear latent
      latent=heatflux0; clear heatflux0
      
      heatflux0(1:1360,:,:)=sensible(81:1440,:,:);
      heatflux0(1361:1440,:,:)=sensible(1:80,:,:); clear sensible
      sensible=heatflux0; clear heatflux0
            
      heatflux0(1:1360,:,:)=longwave(81:1440,:,:);
      heatflux0(1361:1440,:,:)=longwave(1:80,:,:); clear longwave
      longwave=heatflux0; clear heatflux0

      heatflux0(1:1360,:,:)=shortwave(81:1440,:,:);
      heatflux0(1361:1440,:,:)=shortwave(1:80,:,:); clear shortwave
      shortwave=heatflux0; clear heatflux0 
      % figure;imagesc(shortwave(240:360,91:150,1))
      % figure;imagesc(longwave(:,:,1))
    % #####################################################################
      
    
    % #####################################################################
    % Interpolate into 1-deg grids
      disp('   Interpolate 1/4-deg Qnet to 1 degree...')
      [lo,la]=meshgrid((lon_IAP)',(lat_IAP)');
      for month=1:size(mld_monthly,3)
           disp(['    Mon#',num2str(month),' Year#',num2str(year),' interpt Qnet to 1-deg'])
           shortwave_1deg(:,:,month)=griddata(lon_ERA5,lat_ERA5,squeeze(shortwave(:,:,month))',lo',la','linear');
           longwave_1deg(:,:,month) =griddata(lon_ERA5,lat_ERA5,squeeze(longwave(:,:,month))',lo',la','linear');
           sensible_1deg(:,:,month) =griddata(lon_ERA5,lat_ERA5,squeeze(sensible(:,:,month))',lo',la','linear');
           latent_1deg(:,:,month)   =griddata(lon_ERA5,lat_ERA5,squeeze(latent(:,:,month))',lo',la','linear');
      end
      clear lo la month shortwave longwave sensible latent
      % figure;imagesc(shortwave_1deg(:,:,1)) % figure;imagesc(shortwave(:,:,1))
      shortwave= shortwave_1deg; clear shortwave_1deg
      longwave = longwave_1deg;  clear longwave_1deg
      sensible = sensible_1deg;  clear sensible_1deg
      latent   = latent_1deg;    clear latent_1deg
      
      mld_monthly(mld_monthly==0)=NaN;
      shortwave(isnan(mld_monthly))=NaN;
      longwave(isnan(mld_monthly))=NaN;
      sensible(isnan(mld_monthly))=NaN;
      latent(isnan(mld_monthly))=NaN;
      clear mld_monthly
      % figure;imagesc(shortwave(:,:,1))
    % #####################################################################

      cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg/')
      save(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'_1deg.mat'],...
            'shortwave','latent','longwave','sensible','lon_IAP','lat_IAP','-v7.3')
      clear  shortwave latent longwave sensible
end
clear 
% #########################################################################



% #########################################################################
clc;clear   
% #########################################################################
for year=2023
    disp([' Heat flux into 1deg Year#',num2str(year)])
    % #####################################################################
    % Monthly MLD from IAP
    disp(['   MLD from IAPv3, Year#',num2str(year)])
      cd('/Users/z5195509/Documents/Data/IAP_Ocean/Global_GSW_Monthly/MLD')
      % dRh0=0.125 ########################################################
      load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld','lon_IAP','lat_IAP')
      lat_IAP=lat_IAP(:,1);
      mld0(:,:,1:9)=squeeze(mld(:,:,1,1:9));
      clear mld
      load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAPv4_Temp_2023.mat'],'mld')
      mld0(:,:,10:12)=squeeze(mld(:,:,1,10:12));
      clear mld
      
      lon_IAP0(1:340,1)=lon_IAP(21:360,1);
      lon_IAP0(341:360,1)=lon_IAP(1:20,1)+360; clear lon_IAP
      lon_IAP=lon_IAP0; clear lon_IAP0
      
      mld00(1:340,:,:)=mld0(21:360,:,:);
      mld00(341:360,:,:)=mld0(1:20,:,:); clear mld0
      mld_monthly=mld00; clear mld00
      % figure;imagesc(mld_monthly(:,:,1))
      % figure;plot(1:12,squeeze(mld_monthly(240,35,1:12)))
    % #####################################################################

    
    % #####################################################################
    % Monthly heat flux terms
      cd('/Users/z5195509/Documents/Data/ERA-5/Monthly_Averaged_Reanalysis_on_Single_Levels/Heat_Flux/')
      load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'.mat'],...
            'latent','longwave','sensible','shortwave','lon','lat')   
      % 20E-380E
      lon0(1:1360,1)=lon(81:1440,1);
      lon0(1361:1440,1)=lon(1:80,1)+360; clear lon
      lon_ERA5=lon0; clear lon0
      lat_ERA5=lat;  clear lat
      
      heatflux0(1:1360,:,:)=latent(81:1440,:,:);
      heatflux0(1361:1440,:,:)=latent(1:80,:,:); clear latent
      latent=heatflux0; clear heatflux0
      
      heatflux0(1:1360,:,:)=sensible(81:1440,:,:);
      heatflux0(1361:1440,:,:)=sensible(1:80,:,:); clear sensible
      sensible=heatflux0; clear heatflux0
            
      heatflux0(1:1360,:,:)=longwave(81:1440,:,:);
      heatflux0(1361:1440,:,:)=longwave(1:80,:,:); clear longwave
      longwave=heatflux0; clear heatflux0

      heatflux0(1:1360,:,:)=shortwave(81:1440,:,:);
      heatflux0(1361:1440,:,:)=shortwave(1:80,:,:); clear shortwave
      shortwave=heatflux0; clear heatflux0 
      % figure;imagesc(shortwave(240:360,91:150,1))
      % figure;imagesc(longwave(:,:,1))
    % #####################################################################
      
    
    % #####################################################################
    % Interpolate into 1-deg grids
      disp('   Interpolate 1/4-deg Qnet to 1 degree...')
      [lo,la]=meshgrid((lon_IAP)',(lat_IAP)');
      for month=1:size(mld_monthly,3)
           disp(['    Mon#',num2str(month),' Year#',num2str(year),' interpt Qnet to 1-deg'])
           shortwave_1deg(:,:,month)=griddata(lon_ERA5,lat_ERA5,squeeze(shortwave(:,:,month))',lo',la','linear');
           longwave_1deg(:,:,month) =griddata(lon_ERA5,lat_ERA5,squeeze(longwave(:,:,month))',lo',la','linear');
           sensible_1deg(:,:,month) =griddata(lon_ERA5,lat_ERA5,squeeze(sensible(:,:,month))',lo',la','linear');
           latent_1deg(:,:,month)   =griddata(lon_ERA5,lat_ERA5,squeeze(latent(:,:,month))',lo',la','linear');
      end
      clear lo la month shortwave longwave sensible latent
      % figure;imagesc(shortwave_1deg(:,:,1)) % figure;imagesc(shortwave(:,:,1))
      shortwave= shortwave_1deg; clear shortwave_1deg
      longwave = longwave_1deg;  clear longwave_1deg
      sensible = sensible_1deg;  clear sensible_1deg
      latent   = latent_1deg;    clear latent_1deg
      
      mld_monthly(mld_monthly==0)=NaN;
      shortwave(isnan(mld_monthly))=NaN;
      longwave(isnan(mld_monthly))=NaN;
      sensible(isnan(mld_monthly))=NaN;
      latent(isnan(mld_monthly))=NaN;
      clear mld_monthly
      % figure;imagesc(shortwave(:,:,1))
    % #####################################################################

      cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg/')
      save(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'_1deg.mat'],...
            'shortwave','latent','longwave','sensible','lon_IAP','lat_IAP','-v7.3')
      clear  shortwave latent longwave sensible
end
clear 

% #########################################################################
% #########################################################################  



%% 3. ML Heat Budget Using Monthly MLD and Qnet - V2 Accumulation of heat flux anomalies
% #########################################################################
% 
% % #########################################################################
% % 3-1. 1981-2010 climatology of monthly heat flux
% clc;clear
% time_ann=(1981:2010)';
%      
% % #########################################################################
% % #########################################################################
% count_yr=1;
% for year=1981:2010
%     disp([' Heat flux climtology Year#',num2str(year)])
%     % #####################################################################
%     % Monthly MLD from IAP
%     disp(['   MLD from IAPv3, Year#',num2str(year)])
%       cd('/Users/z5195509/Documents/Data/IAP_Ocean/Global_GSW_Monthly/MLD')
%       % dRh0=0.125 ########################################################
%       load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld','lon_IAP','lat_IAP')
%       lat_IAP=lat_IAP(:,1);
%       mld0(:,:,1:12)=squeeze(mld(:,:,1,1:12));
%       clear mld
%       
%       lon_IAP0(1:340,1)=lon_IAP(21:360,1);
%       lon_IAP0(341:360,1)=lon_IAP(1:20,1)+360; clear lon_IAP
%       lon_IAP=lon_IAP0; clear lon_IAP0
%       
%       mld00(1:340,:,:)=mld0(21:360,:,:);
%       mld00(341:360,:,:)=mld0(1:20,:,:); clear mld0
%       mld_monthly=mld00; clear mld00
%       % figure;imagesc(mld_monthly(:,:,1)); caxis([0 100])
%       % figure;plot(1:12,squeeze(mld_monthly(240,35,1:12)))
%     % #####################################################################
% 
%         
%     % Monthly heat flux terms
%       cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg/')
%       load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'_1deg.mat'],...
%             'shortwave','latent','longwave','sensible')
%       % figure;imagesc(shortwave(:,:,1))
%       latent_clim(:,:,:,count_yr)    = latent;
%       longwave_clim(:,:,:,count_yr)  = longwave;
%       sensible_clim(:,:,:,count_yr)  = sensible;
%       shortwave_clim(:,:,:,count_yr) = shortwave;
%       Qnet_clim(:,:,:,count_yr)      = latent+longwave+sensible+shortwave;
%       clear latent longwave sensible
%       
%       % ###################################################################
%       % Remove the shortwave penatration
%         % The Vertical Redistribution of SWR, at MLD=mld_2023_daily
%         R=0.58; h1=0.35; h2=23; 
%         z=mld_monthly; 
%         F_PS77=R.*exp(-z./h1)+(1-R).*exp(-z./h2);
%         
%         shortwave_mld(:,:,:,count_yr)= - shortwave.*F_PS77;
%         clear R h1 h2 z F_PS77 mld_monthly
%         % figure;imagesc(shortwave_mld(240:360,91:150,1))
% %         % Qnet and Shortwave radiation after removing the shortwave penatration at mixed-layer base
% %         Qnet_clim(:,:,:,count_yr)      = Qnet     -shortwave_mld;
% %         shortwave_clim(:,:,:,count_yr) = shortwave-shortwave_mld;
% %         clear shortwave_mld Qnet shortwave
%       % ###################################################################    
%       % ###################################################################
%       clear Qnet shortwave
%       count_yr=count_yr+1;
% end
% clear year count_yr
% 
% Qnet_clim     = nanmean(Qnet_clim,4);
% shortwave_clim= nanmean(shortwave_clim,4);
% shortwave_mld = nanmean(shortwave_mld,4);
% latent_clim   = nanmean(latent_clim,4);
% longwave_clim = nanmean(longwave_clim,4);
% sensible_clim = nanmean(sensible_clim,4);
% % figure;imagesc(Qnet_clim(:,:,1))
% % figure;imagesc(shortwave_mld(:,:,1))
% 
% cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg/')
% save(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_19812010clim_1deg.mat'],...
%       'Qnet_clim','shortwave_clim','shortwave_mld','latent_clim','longwave_clim','sensible_clim',...
%       'lon_IAP','lat_IAP','-v7.3')
% % #########################################################################
% % #########################################################################
% 
%   
% 
% 
% %% #########################################################################
% %  #########################################################################
% % 3-2. Accumulation of heat flux anomalie in 1981-2023
% clc;clear
% time_ann=(1981:2010)';
% 
% % #########################################################################
% %   % ACCESS Om2 0.25
% %     cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/access-om2-025_daily_sst_fields/')
% %     load('basin_mask_NA_ACCESS-OM2-025_era5_iaf.mat','basin_mask_NA','lon_025','lat_025')
% %     
% %     cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
% %     load(['SST_NA_100W20E_Eq70N_ERA5_Daily_Averaged_Reanalysis_on_Single_Levels_2022.mat'],'sst_daily_NA','lon','lat')
% %        [lo,la]=meshgrid((lon)', (lat)');
% %        basin_mask_NA=griddata(lon_025,lat_025,basin_mask_NA',lo',la','nearest');
% %        clear lo la lon_025 lat_025
% %        sst_daily_NA=nanmean(sst_daily_NA,3);
% %        basin_mask_NA(isnan(sst_daily_NA))=NaN; clear sst_daily_NA
% %     
% %        [Sxy,~,~]=function_Cgrid_Area_Distance((lon)',(lat)');
% %        Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
% %        Sxy_NA=Sxy;
% %      save(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','Sxy_NA','lon','lat')
% 
%      % cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
%      load(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','lon','lat')
%        [lo,la]=meshgrid((260:380)', (0.5:59.5)');
%        basin_mask_NA=griddata(lon,lat,basin_mask_NA',lo',la','nearest');
%        clear lo la lon lat
%        [Sxy,~,~]=function_Cgrid_Area_Distance((260:380)',(0.5:59.5)');
%        Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
%        Sxy_NA=Sxy; clear Sxy
% % #########################################################################
% 
% 
% % #########################################################################
% % Net heat flux climatology in 1981-2010 
% cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg/')
% load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_19812010clim_1deg.mat'],...
%       'Qnet_clim','shortwave_clim','shortwave_mld','latent_clim','longwave_clim','sensible_clim',...
%       'lon_IAP','lat_IAP')
%   shortwave_mld_clim=shortwave_mld; clear shortwave_mld
%   Sxy_NA(isnan(shortwave_mld_clim(240:360,91:150,1)))=NaN;
%   
%   % figure;imagesc(Qnet_clim(240:360,91:150,1))
%   % figure;imagesc(shortwave_mld_clim(240:360,91:150,1))
%   % figure;imagesc(Sxy_NA)
% % #########################################################################
%      
% 
% % #########################################################################
% count_yr=1;
% for year=1981:2023
%     disp([' Predicting NA daily SST Year#',num2str(year)])
%     % #####################################################################
%     % Monthly MLD from IAP
%     disp(['   MLD from IAPv3, Year#',num2str(year)])
%       cd('/Users/z5195509/Documents/Data/IAP_Ocean/Global_GSW_Monthly/MLD')
%       % dRh0=0.125 ########################################################
%       if year<2023
%           load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld','lon_IAP','lat_IAP')
%           lat_IAP=lat_IAP(:,1);
%           mld0(:,:,1:12)=squeeze(mld(:,:,1,1:12));
%           clear mld
%       else
%           % dRh0=0.125 ########################################################
%           load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld','lon_IAP','lat_IAP')
%           lat_IAP=lat_IAP(:,1);
%           mld0(:,:,1:9)=squeeze(mld(:,:,1,1:9));
%           clear mld
%           load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAPv4_Temp_2023.mat'],'mld')
%           mld0(:,:,10:12)=squeeze(mld(:,:,1,10:12));
%           clear mld
%       end
%       
%       lon_IAP0(1:340,1)=lon_IAP(21:360,1);
%       lon_IAP0(341:360,1)=lon_IAP(1:20,1)+360; clear lon_IAP
%       lon_IAP=lon_IAP0; clear lon_IAP0
%       
%       mld00(1:340,:,:)=mld0(21:360,:,:);
%       mld00(341:360,:,:)=mld0(1:20,:,:); clear mld0
%       mld_monthly=mld00; clear mld00
%       % figure;imagesc(mld_monthly(:,:,1))
%       % figure;plot(1:12,squeeze(mld_monthly(240,35,1:12)))
%     % #####################################################################
%     
%         
%     % Monthly heat flux terms
%       cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg/')
%       load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'_1deg.mat'],...
%             'shortwave','latent','longwave','sensible') 
%       Qnet=latent+longwave+sensible+shortwave;
%       % figure;imagesc(Qnet(240:360,91:150,1))
%         
%       % ##################################################################
%       % Remove the shortwave penatration
%         % The Vertical Redistribution of SWR, at MLD=mld_2023_daily
%         R=0.58; h1=0.35; h2=23; 
%         z=mld_monthly;
%         F_PS77=R.*exp(-z./h1)+(1-R).*exp(-z./h2);
%         
%         shortwave_mld = - shortwave.*F_PS77;
%         clear R h1 h2 z F_PS77 
%         % figure;imagesc(shortwave_mld(240:360,91:150,1))
% %         % Qnet and Shortwave radiation after removing the shortwave penatration at mixed-layer base
% %         Qnet_mld  = Qnet     -shortwave_mld;
% %         shortwave = shortwave-shortwave_mld;
% %         clear shortwave_mld Qnet
%         
% %       Qnet_mld      = Qnet;
% %       shortwave = shortwave;
% %       clear Qnet
%       % ##################################################################   
%        
%       
%       % ##################################################################   
%       % Heat flux anomaly relative to 1981-2010 climatology 
%         Qnet      = Qnet      - Qnet_clim;
%         shortwave = shortwave - shortwave_clim;
%         latent    = latent    - latent_clim;
%         sensible  = sensible  - sensible_clim;
%         longwave  = longwave  - longwave_clim;
%         
%         shortwave_mld = shortwave_mld - shortwave_mld_clim;
%         Qnet_mld  = Qnet + shortwave_mld;
%       % ##################################################################
%       
%        
%       % ##################################################################
%       % Mixed layer warming due to Qnet
%         disp(['   Projected ML Warming in Year#',num2str(year)])
%         % Heat flux term in heat budget
%         MLT_Qnet_mon(:,:,:)  = Qnet./mld_monthly./3992./1027.*(24*60*60*30);  % K/month
%         clear Qnet
%         MLT_Qswr_mon(:,:,:)  = shortwave./mld_monthly./3992./1027.*(24*60*60*30); % K/month
%         clear shortwave
%         MLT_Qlat_mon(:,:,:)  = latent./mld_monthly./3992./1027.*(24*60*60*30);    % K/month
%         clear latent
%         MLT_Qlon_mon(:,:,:)  = longwave./mld_monthly./3992./1027.*(24*60*60*30);  % K/month
%         clear longwave
%         MLT_Qsen_mon(:,:,:)  = sensible./mld_monthly./3992./1027.*(24*60*60*30);  % K/month
%         clear sensible
%         
%         MLT_Qnet_mld_mon(:,:,:)  = Qnet_mld./mld_monthly./3992./1027.*(24*60*60*30);  % K/month
%         clear Qnet_mld
%         MLT_Qswr_mld_mon(:,:,:)  = shortwave_mld./mld_monthly./3992./1027.*(24*60*60*30); % K/month
%         clear shortwave_mld
%         clear mld_monthly
%         % figure;imagesc(MLT_Qnet_mon(240:360,91:150,1)); caxis([-3 3])
% 
%         
%          % ################################################################
%          % NA (100W-20E,0-60N) averaged MLT tendency by heat flux terms 
%          dMLT_Qnet_mon0=MLT_Qnet_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qswr_mon0=MLT_Qswr_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qlat_mon0=MLT_Qlat_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qlon_mon0=MLT_Qlon_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qsen_mon0=MLT_Qsen_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qnet_mld_mon0=MLT_Qnet_mld_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          dMLT_Qswr_mld_mon0=MLT_Qswr_mld_mon(240:360,91:150,:).*Sxy_NA(:,:);
%          % figure;imagesc(dMLT_Qnet_mon0(:,:,1)); % figure;imagesc(dMLT_Qswr_mld_mon0(:,:,1));
%          % figure;imagesc(Sxy_NA(:,:,1)); 
%          Sxy_NA0=Sxy_NA;
%          for month=2:12
%             Sxy_NA0(:,:,month)=Sxy_NA0(:,:,1); 
%          end
%          Sxy_NA0(isnan(dMLT_Qnet_mon0))=NaN; 
%          dMLT_Qnet_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qswr_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qlat_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qlon_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qsen_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qnet_mld_mon0(isnan(Sxy_NA0))=NaN; 
%          dMLT_Qswr_mld_mon0(isnan(Sxy_NA0))=NaN; 
%          Sxy_NA0(isnan(dMLT_Qnet_mon0))=NaN; 
%          % figure;imagesc(dMLT_Qnet_mon0(:,:,1)); 
%          % figure;imagesc(Sxy_NA0(:,:,1)); 
%          
%          dMLT_Qnet_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qnet_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qswr_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qswr_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qlat_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qlat_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qlon_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qlon_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qsen_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qsen_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qnet_mld_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qnet_mld_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          dMLT_Qswr_mld_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qswr_mld_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          clear dMLT_*_mon0 Sxy_NA0
%       % ###################################################################
%          
%         
%       % ###################################################################
%       % Predicted MLT = accummulated / intgrated heat flux anomalies
%         project_sst_Qnet=nan(size(MLT_Qnet_mon));
%         project_sst_Qswr=nan(size(MLT_Qnet_mon));
%         project_sst_Qlat=nan(size(MLT_Qnet_mon));
%         project_sst_Qlon=nan(size(MLT_Qnet_mon));
%         project_sst_Qsen=nan(size(MLT_Qnet_mon));
%         project_sst_Qnet_mld=nan(size(MLT_Qnet_mon));
%         project_sst_Qswr_mld=nan(size(MLT_Qnet_mon));
%         
%         
%         % Jan-15
%         cd('/Users/z5195509/Documents/Data/ERA-5/Monthly_Averaged_Reanalysis_on_Single_Levels/SST_SAT/')
%         load(['SST_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'.mat'],...
%               'sst')
%         % 20E-380E
%         heatflux0(1:1360,:,1)=sst(81:1440,:,1);
%         heatflux0(1361:1440,:,1)=sst(1:80,:,1); clear sst
%         sst_daily0=heatflux0; clear heatflux0
%         
%            [lo,la]=meshgrid((lon_IAP)',(lat_IAP)');
%            sst_daily(:,:,1)=griddata((20:0.25:379.75)',(-90:0.25:90)',squeeze(sst_daily0(:,:,1))',lo',la','linear');
%         
%         sst_daily(isnan(MLT_Qnet_mon(:,:,1)))=NaN;
%          % figure;imagesc(sst_daily(240:360,91:150,1)); 
%          % figure;imagesc(MLT_Qnet_mon(240:360,91:150,1)); 
%          % figure;imagesc(Sxy_NA(:,1:end,1)); 
%         
%         
%         % SST Jan-15
%         project_sst_Qnet(:,:,1)=sst_daily(:,:,1);
%         project_sst_Qswr(:,:,1)=sst_daily(:,:,1);
%         project_sst_Qlat(:,:,1)=sst_daily(:,:,1);
%         project_sst_Qlon(:,:,1)=sst_daily(:,:,1);
%         project_sst_Qsen(:,:,1)=sst_daily(:,:,1);
%         project_sst_Qnet_mld(:,:,1)=sst_daily(:,:,1);
%         project_sst_Qswr_mld(:,:,1)=sst_daily(:,:,1);
%         
%         % Projected SST from Feb-15
%         for month=2:12
%             disp(['   Projeted monthly MLT mon#',num2str(month),' year#',num2str(year)])
%             project_sst_Qnet(:,:,month)=sst_daily(:,:,1)+squeeze(nansum(MLT_Qnet_mon(:,:,2:month),3));
%             project_sst_Qswr(:,:,month)=sst_daily(:,:,1)+squeeze(nansum(MLT_Qswr_mon(:,:,2:month),3));
%             project_sst_Qlat(:,:,month)=sst_daily(:,:,1)+squeeze(nansum(MLT_Qlat_mon(:,:,2:month),3));
%             project_sst_Qlon(:,:,month)=sst_daily(:,:,1)+squeeze(nansum(MLT_Qlon_mon(:,:,2:month),3));
%             project_sst_Qsen(:,:,month)=sst_daily(:,:,1)+squeeze(nansum(MLT_Qsen_mon(:,:,2:month),3));
%             project_sst_Qnet_mld(:,:,month)=sst_daily(:,:,1)+squeeze(nansum(MLT_Qnet_mld_mon(:,:,2:month),3));
%             project_sst_Qswr_mld(:,:,month)=sst_daily(:,:,1)+squeeze(nansum(MLT_Qswr_mld_mon(:,:,2:month),3));
%         end
%         clear sst_daily MLT_*_mon
%         project_sst_Qnet(project_sst_Qnet==0)=NaN;
%         project_sst_Qswr(project_sst_Qswr==0)=NaN;
%         project_sst_Qlat(project_sst_Qlat==0)=NaN;
%         project_sst_Qlon(project_sst_Qlon==0)=NaN;
%         project_sst_Qnet_mld(project_sst_Qnet_mld==0)=NaN;
%         project_sst_Qswr_mld(project_sst_Qswr_mld==0)=NaN;
%         
%         % figure;imagesc(project_sst_Qnet(240:360,91:150,1)); caxis([-3 3])
%         
% 
%         % NA (100W-20E,0-60N) averaged projected SST 
%          project_sst_Qnet=project_sst_Qnet(240:360,91:150,:).*Sxy_NA(:,1:end);
%          project_sst_Qswr=project_sst_Qswr(240:360,91:150,:).*Sxy_NA(:,1:end);
%          project_sst_Qlat=project_sst_Qlat(240:360,91:150,:).*Sxy_NA(:,1:end);
%          project_sst_Qlon=project_sst_Qlon(240:360,91:150,:).*Sxy_NA(:,1:end);
%          project_sst_Qsen=project_sst_Qsen(240:360,91:150,:).*Sxy_NA(:,1:end);
%          project_sst_Qnet_mld=project_sst_Qnet_mld(240:360,91:150,:).*Sxy_NA(:,1:end);
%          project_sst_Qswr_mld=project_sst_Qswr_mld(240:360,91:150,:).*Sxy_NA(:,1:end);
%          % figure;imagesc(project_sst_Qnet(:,:,1)); 
%          % figure;imagesc(project_sst_Qswr_mld(:,:,1));
%          % figure;imagesc(Sxy_NA(:,1:end,1)); 
%          Sxy_NA0=Sxy_NA;
%          for month=2:12
%             Sxy_NA0(:,:,month)=Sxy_NA0(:,:,1); 
%          end
%          Sxy_NA0(isnan(project_sst_Qnet))=NaN; 
%          project_sst_Qnet(isnan(Sxy_NA0))=NaN; 
%          project_sst_Qswr(isnan(Sxy_NA0))=NaN; 
%          project_sst_Qlat(isnan(Sxy_NA0))=NaN; 
%          project_sst_Qlon(isnan(Sxy_NA0))=NaN; 
%          project_sst_Qsen(isnan(Sxy_NA0))=NaN; 
%          project_sst_Qnet_mld(isnan(Sxy_NA0))=NaN; 
%          project_sst_Qswr_mld(isnan(Sxy_NA0))=NaN; 
%          Sxy_NA0(isnan(project_sst_Qnet))=NaN; 
%          % figure;imagesc(project_sst_Qnet(:,:,1));  caxis([-1e11 1e11])
%          % figure;imagesc(Sxy_NA0(:,1:end,1)); caxis([0 1e10])
%          
%          proj_NA_SST_mon_Qnet(1:12,count_yr)=squeeze(nansum(nansum(project_sst_Qnet(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          proj_NA_SST_mon_Qswr(1:12,count_yr)=squeeze(nansum(nansum(project_sst_Qswr(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          proj_NA_SST_mon_Qlat(1:12,count_yr)=squeeze(nansum(nansum(project_sst_Qlat(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          proj_NA_SST_mon_Qlon(1:12,count_yr)=squeeze(nansum(nansum(project_sst_Qlon(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          proj_NA_SST_mon_Qsen(1:12,count_yr)=squeeze(nansum(nansum(project_sst_Qsen(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          proj_NA_SST_mon_Qnet_mld(1:12,count_yr)=squeeze(nansum(nansum(project_sst_Qnet_mld(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          proj_NA_SST_mon_Qswr_mld(1:12,count_yr)=squeeze(nansum(nansum(project_sst_Qswr_mld(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
%          clear project_sst_* Sxy_NA0
%          % ################################################################
%          
%          
%          % ################################################################
%          % Saving every year
%          NA_SST_mon_Qnet(1:12,1)=squeeze(proj_NA_SST_mon_Qnet(1:12,count_yr));
%          NA_SST_mon_Qswr(1:12,1)=squeeze(proj_NA_SST_mon_Qswr(1:12,count_yr));
%          NA_SST_mon_Qlat(1:12,1)=squeeze(proj_NA_SST_mon_Qlat(1:12,count_yr));
%          NA_SST_mon_Qlon(1:12,1)=squeeze(proj_NA_SST_mon_Qlon(1:12,count_yr));
%          NA_SST_mon_Qsen(1:12,1)=squeeze(proj_NA_SST_mon_Qsen(1:12,count_yr));
%          NA_SST_mon_Qnet_mld(1:12,1)=squeeze(proj_NA_SST_mon_Qnet_mld(1:12,count_yr));
%          NA_SST_mon_Qswr_mld(1:12,1)=squeeze(proj_NA_SST_mon_Qswr_mld(1:12,count_yr));
%          
%          cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
%          save(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Anomaly_Accumulated_',num2str(year),'.mat'],...
%               'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen','NA_SST_mon_Qnet_mld','NA_SST_mon_Qswr_mld',...
%               'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon','dMLT_Qnet_mld_mon','dMLT_Qswr_mld_mon',...
%               'basin_mask_NA','Sxy_NA','lon_IAP','lat_IAP')
%           clear NA_SST_mon* dMLT_*_mon
%          % ################################################################
%       count_yr=count_yr+1;
% end
% clear *_dail_clim
% % #########################################################################
% % #########################################################################  




%% #######################################################################
%% 4. Figure# Anomaly of Monthly MLT, SST, and MLT decomposition in 2023
% % Part 1: Twitter plot for daily SSTA
% clc;clear
% time_ann=(1981:2022)';
% 
% % #########################################################################
% % #########################################################################
% % Mixed layer warming by air-sea flues
% count_yr=1;
% for year=1981:2023
%     cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
%     load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Anomaly_Accumulated_',num2str(year),'.mat'],...
%           'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen','NA_SST_mon_Qnet_mld','NA_SST_mon_Qswr_mld',...
%           'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon','dMLT_Qnet_mld_mon','dMLT_Qswr_mld_mon',...
%           'basin_mask_NA','Sxy_NA','lon_IAP','lat_IAP')
%     NA_MLT_Qnet(:,count_yr)=NA_SST_mon_Qnet; 
%     NA_MLT_Qswr(:,count_yr)=NA_SST_mon_Qswr; 
%     NA_MLT_Qlat(:,count_yr)=NA_SST_mon_Qlat; 
%     NA_MLT_Qlon(:,count_yr)=NA_SST_mon_Qlon; 
%     NA_MLT_Qsen(:,count_yr)=NA_SST_mon_Qsen; 
%     NA_MLT_Qswr_mld(:,count_yr)=NA_SST_mon_Qswr_mld; 
%     clear NA_SST_mon*  
%     count_yr=count_yr+1;
% end
% % Climatology
%     NA_MLT_Qnet_clim=nanmean(NA_MLT_Qnet(:,1:30),2);
%     NA_MLT_Qswr_clim=nanmean(NA_MLT_Qswr(:,1:30),2);
%     NA_MLT_Qlat_clim=nanmean(NA_MLT_Qlat(:,1:30),2);
%     NA_MLT_Qlon_clim=nanmean(NA_MLT_Qlon(:,1:30),2);
%     NA_MLT_Qsen_clim=nanmean(NA_MLT_Qsen(:,1:30),2);
%     NA_MLT_Qswr_mld_clim=nanmean(NA_MLT_Qswr_mld(:,1:30),2);
% 
%     
% % Anomaly relative to the clim-monthly mean
%     count=1;
%     for year=1:size(NA_MLT_Qnet,2)
%         NA_MLT_Qnet_ano(:,count)=NA_MLT_Qnet(:,count)-NA_MLT_Qnet_clim(1,1);
%         NA_MLT_Qswr_ano(:,count)=NA_MLT_Qswr(:,count)-NA_MLT_Qswr_clim(1,1);
%         NA_MLT_Qlat_ano(:,count)=NA_MLT_Qlat(:,count)-NA_MLT_Qlat_clim(1,1);
%         NA_MLT_Qlon_ano(:,count)=NA_MLT_Qlon(:,count)-NA_MLT_Qlon_clim(1,1);
%         NA_MLT_Qsen_ano(:,count)=NA_MLT_Qsen(:,count)-NA_MLT_Qsen_clim(1,1);
%         NA_MLT_Qswr_mld_ano(:,count)=NA_MLT_Qswr_mld(:,count)-NA_MLT_Qswr_mld_clim(1,1);
%         count=count+1;
%     end
%     %clear *_clim
% % #########################################################################
% % #########################################################################
% 
% 
% 
% 
% % #########################################################################
% % SSTa from ERA5 reanalysis
% % #########################################################################
% %   % ACCESS Om2 0.25
% %     cd('/Users/z5195509/Documents/MATLAB/Function_MATLAB')
% %     basin_mask=ncread('Basin_mask_ERA5-global-land-sea-mask.nc','lsm');
% %         basin_mask=fliplr(basin_mask);
% %         basin_mask(basin_mask>0.5)=NaN; % land
% %         basin_mask(basin_mask<=0.5)=2;  % ocean
% %     lon_mask=ncread('Basin_mask_ERA5-global-land-sea-mask.nc','longitude');
% %     lat_mask=ncread('Basin_mask_ERA5-global-land-sea-mask.nc','latitude');
% %         lat_mask=flipud(lat_mask);
% %         lon_mask(1441:1520,1)=lon_mask(1:80,1)+360;
% %         basin_mask(1441:1520,:)=basin_mask(1:80,:);
% %         
% %         % Basin mask in NA
% %         lat_mask=lat_mask(361:641,1);
% %         lon_mask=lon_mask(1041:1520,1);
% %         basin_mask_NA=basin_mask(1041:1520,361:641);
% %        
% %         [Sxy,~,~]=function_Cgrid_Area_Distance((260:0.25:379.75)',(0:0.25:70)');
% %         Sxy(isnan(basin_mask_NA))=NaN;
% 
%     cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/access-om2-025_daily_NA_SST_fields/')
%     % load('basin_mask_NA_ACCESS-OM2-025_era5_iaf_no_medi.mat','basin_mask_NA','lon_025','lat_025')
%     load('basin_mask_NA_ACCESS-OM2-025_era5_iaf.mat','basin_mask_NA','lon_025','lat_025')
%     
%     cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
%     load(['SST_NA_100W20E_Eq70N_ERA5_Daily_Averaged_Reanalysis_on_Single_Levels_2022.mat'],'sst_daily_NA','lon','lat')
%        [lo,la]=meshgrid((lon)', (lat)');
%        basin_mask_NA=griddata(lon_025,lat_025,basin_mask_NA',lo',la','nearest');
%        clear lo la lon_025 lat_025
%        sst_daily_NA=nanmean(sst_daily_NA,3);
%        basin_mask_NA(isnan(sst_daily_NA))=NaN; clear sst_daily_NA
%     
%        [Sxy,~,~]=function_Cgrid_Area_Distance((lon)',(lat)');
%        Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
% % #########################################################################
% 
% 
% % #########################################################################
%     % Daily SST during 1981-2010
%     cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
%     sst_daily_19812022=nan(366,length(1981:2022)');
%     count=1;
%     for year=1981:2022
%           disp(['Daily SST Year# ',num2str(year)])
%           load(['SST_NA_100W20E_Eq70N_ERA5_Daily_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'.mat'],'sst_daily_NA','lon','lat')
% 
%           for days=1:size(sst_daily_NA,3)
%               sst_daily_NA0=squeeze(sst_daily_NA(:,:,days));
%               sst_daily_NA0(isnan(basin_mask_NA))=NaN;
%               sst_daily_NA(:,:,days)=sst_daily_NA0;
%               clear sst_daily_NA0
%           end
%           clear days
%           sst_daily_NA=sst_daily_NA.*Sxy;
%           
%           % NA 0-60N
%           sst_daily_19812022(1:size(sst_daily_NA,3),count)=squeeze(nansum(nansum(sst_daily_NA(:,1:241,:),2),1))./squeeze(nansum(nansum(Sxy(:,1:241),2),1));
%           clear sst_daily_NA
%           count=count+1;
%     end
%     clear year count
% 
%     
%     % 2023 Daily SST Anomaly relative to the clim-mean
%           disp(['Daily SST Year# 2023'])
%           load('SST_NA_100W20E_Eq70N_ERA5_Daily_Averaged_Reanalysis_on_Single_Levels_2023.mat','sst_daily')
%           sst_daily_NA=sst_daily; clear sst_daily
% 
%           for days=1:size(sst_daily_NA,3)
%               sst_daily_NA0=squeeze(sst_daily_NA(:,:,days));
%               sst_daily_NA0(isnan(basin_mask_NA))=NaN;
%               sst_daily_NA(:,:,days)=sst_daily_NA0;
%               clear sst_daily_NA0
%           end
%           clear days
%           sst_daily_NA=sst_daily_NA.*Sxy;
%           
%           % NA 0-60N
%           sst_daily_2023(1:365,1)=squeeze(nansum(nansum(sst_daily_NA(:,1:241,1:365),2),1))./squeeze(nansum(nansum(Sxy(:,1:241),2),1));
%           sst_daily_2023(366,1)=NaN;
%           clear sst_daily_NA
% 
%           
%     % Anomaly is relative to the clim-daily mean
%     sst_daily_19812010_clim=squeeze(nanmean(sst_daily_19812022(:,1:30),2));
%     count=1;
%     for year=1:size(sst_daily_19812022,2)
%         sst_daily_19812022_ano(:,count)=sst_daily_19812022(:,count)-sst_daily_19812010_clim;
%         count=count+1;
%     end
%     
%     sst_daily_2023_ano(:,1)=sst_daily_2023(:,1)-sst_daily_19812010_clim; 
% % #########################################################################
% 
% 
% 
% 
% % #########################################################################
% % % MLTa from IAP data
% % % #########################################################################
% %   % ACCESS Om2 0.25
%      % cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
%      load(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','Sxy_NA','lon','lat')
%        [lo,la]=meshgrid((260:380)', (0.5:59.5)');
%        basin_mask_NA=griddata(lon,lat,basin_mask_NA',lo',la','nearest');
%        clear lo la lon_025 lat_025
%        [Sxy,~,~]=function_Cgrid_Area_Distance((260:380)',(0.5:59.5)');
%        Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
%        Sxy_NA=Sxy; clear Sxy
% % #########################################################################
% 
% % #########################################################################
%     % Monthly MLT during 1981-2023
%     cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023')
%     mlt_month_19812023=nan(12,length(1981:2023)');
%     count=1;
%     for year=1981:2023
%           disp(['Monthly MLT Year# ',num2str(year)])
%           load(['mlt_monthly_Rh0_0125_GSW_IAP_C17_',num2str(year),'.mat'],'mlt_monthly','lon_IAP','lat_IAP')
%           mlt_monthly_NA=mlt_monthly(240:360,91:150,:);
%           clear mlt_monthly
% 
%           for month=1:size(mlt_monthly_NA,3)
%               mlt_monthly_NA0=squeeze(mlt_monthly_NA(:,:,month));
%               mlt_monthly_NA0(isnan(basin_mask_NA))=NaN;
%               mlt_monthly_NA(:,:,month)=mlt_monthly_NA0;
%               clear mlt_monthly_NA0
%           end
%           clear days
%           mlt_monthly_NA=mlt_monthly_NA.*Sxy_NA;
%           
%           % NA 0-60N
%           mlt_month_19812023(1:size(mlt_monthly_NA,3),count)=squeeze(nansum(nansum(mlt_monthly_NA(:,:,:),2),1))./squeeze(nansum(nansum(Sxy_NA(:,:),2),1));
%           clear mlt_monthly_NA
%           count=count+1;
%     end
%     clear year count
% 
%     % Anomaly is relative to the clim-daily mean
%     mlt_month_19812010_clim=squeeze(nanmean(mlt_month_19812023(:,1:30),2));
%     count=1;
%     for year=1:size(mlt_month_19812023,2)
%         mlt_month_19812023_ano(:,count)=mlt_month_19812023(:,count)-mlt_month_19812010_clim;
%         count=count+1;
%     end
% % #########################################################################
% 
% 
% 
% 
% %% ########################################################################
% %  Plotting 4-1: Daily SSTA from ERA5 and MLTA by Heat Flux
% clc;
% figure('Color',[1 1 1]);  %create a new figure of white color background
% ixs = 0.060; ixe = 0.450;  ixd = 0.10; ixw = (1-ixs-ixe-0*ixd)/1;
% iys = 0.050; iye = 0.050;  iyd = 0.08; iyw = (1-iys-iye-1*iyd)/2;
% 
% %          [left            bottom      width height]
% pos{101}  = [ixs          iys+0.1*iyw+1*iyd   ixw 1.7*iyw]; 
% pos{102}  = [ixs          iys+0.0*iyw+0*iyd   ixw 1*iyw]; 
% 
% 
% % dc1                 = [0    , 0.447, 0.741]; % Blue-ish.
% % dc6                 = [0.301, 0.745, 0.933]; % Light Blue-ish.
% % dc2                 = [0.850, 0.325, 0.098]; % Orange/Red-ish.
% % dc7                 = [0.635, 0.078, 0.184]; % Red/Brown-ish.
% % dc5                 = [0.466, 0.674, 0.188]; % Green-ish.
% % dc3                 = [0.929, 0.694, 0.125]; % Yellow/Gold-ish.
% % dc4                 = [0.494, 0.184, 0.556]; % Purple-ish.
% 
% 
% clear color color0 
% color=cbrewer('seq', 'Blues', 60,'pchip');
% color(:,:)=color(60:-1:1,:);
% 
%     
% subplot('position',pos{101})
% %     for year=1:42
% %         line00=plot(1:366,sst_daily_19812022_ano(:,year));
% %         set(line00,'color',color(year,:),'LineWidth',1.5,'linestyle','-'); 
% %         hold on
% %     end
% %     
% %     line00=plot(1:366,sst_daily_19812022_ano(:,5),'LineWidth',2);
% % 
% %     hold on
%     line2023=plot(1:366,sst_daily_2023_ano(:,1));
%     set(line2023,'color',[0.3,0.3,0.3],'LineWidth',3,'linestyle','-'); 
%     
%     hold on
%     line2023MLT3=plot(16:30:270,NA_MLT_Qlat_ano(1:9,43));
%     set(line2023MLT3,'color',[0.466, 0.674, 0.188],'LineWidth',3,'linestyle','-'); 
%     hold on
%     line2023MLT4=plot(16:30:270,NA_MLT_Qlon_ano(1:9,43));
%     set(line2023MLT4,'color',[0.5,0.5,0.5],'LineWidth',3,'linestyle','-'); 
%     hold on
%     line2023MLT5=plot(16:30:270,NA_MLT_Qsen_ano(1:9,43));
%     set(line2023MLT5,'color',[0.301, 0.745, 0.933],'LineWidth',3,'linestyle','-'); 
% %     hold on
% %     line2023MLT6=plot(16:30:270,NA_MLT_Qswr_mld_ano(1:9,43));
% %     set(line2023MLT6,'color',[0.494, 0.184, 0.556],'LineWidth',3,'linestyle','--'); 
%     hold on
%     line2023MLT2=plot(16:30:270,NA_MLT_Qswr_ano(1:9,43));
%     set(line2023MLT2,'color',[0.494, 0.184, 0.556],'LineWidth',3,'linestyle','-'); 
%     hold on
%     line2023MLT1=plot(16:30:270,NA_MLT_Qnet_ano(1:9,43));
%     set(line2023MLT1,'color',[0.929, 0.694, 0.125],'LineWidth',6,'linestyle','-'); 
%     
%     hold on
%     line2023MLT0=plot(16:30:270,mlt_month_19812023_ano(1:9,43));
%     set(line2023MLT0,'color',[0.9, 0.3, 0.3],'LineWidth',6,'linestyle','-'); 
% 
%     
% %     leg101=legend([line2023 line2023MLT0 line2023MLT1 line2023MLT2 line2023MLT6 line2023MLT5 line2023MLT4 line2023MLT3],...
% %                '2023 SSTA, ERA5','2023 MLTA','2023 MLTA, Qnet','2023 MLTA, shortwave','2023 MLTA, shortwave_{mlb}','2023 MLTA, sensible','2023 MLTA, longwave','2023 MLTA, latent heat',...
% %                'Location','northwest','NumColumns',1);
%     leg101=legend([line2023 line2023MLT0 line2023MLT1 line2023MLT2 line2023MLT5 line2023MLT4 line2023MLT3],...
%                '2023 SSTA, ERA5','2023 MLTA','2023 MLTA, Qnet','2023 MLTA, shortwave','2023 MLTA, sensible','2023 MLTA, longwave','2023 MLTA, latent heat',...
%                'Location','northwest','NumColumns',1);
%     set(leg101,'fontsize',20)
%     hold on
%     title(leg101,'North Atlantic SST and MLT Anomaly','fontsize',20')
%     legend('boxoff')
% 
%     
%     set(gca,'Ylim',[0 1.5],'ycolor','k') 
%     set(gca,'YTick',-1.5:0.5:4)
%     set(gca,'YTickLabel',{'-1.5','-1.0','-0.5','0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0'},'fontsize',16)
%     set(gca,'Xlim',[1 366]) 
%     set(gca,'XTick',1:30:360)
%     set(gca,'XTickLabel',{'           Jan','           Feb','           Mar','           Apr','           May','           Jun',...
%                           '           Jul','           Aug','           Sep','           Oct','           Nov','            Dec'},'fontsize',20)
% 
%     grid on
%     set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
%     ylabel(['Difference [ ^oC ]'],'fontsize',20,'color','k','FontWeight','normal')
% 
%     title('a. North Atlantic SST and predicted MLT anomaly','fontsize',24,'color','k','FontWeight','bold')
% 
%     
%         
% cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')




%% ########################################################################
% %  Plotting 4-2: Bar Charts for Monthly MLTA and ML Warming Decomposition by Monthly Heat Flux
% clc;clear
% time_ann=(1981:2022)';
% 
% % #########################################################################
% % #########################################################################
% % Mixed layer warming by air-sea flues in 2023 - accumuated monthly heat flux anomalies
% count_yr=1;
% for year=2023
%     cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
%     load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Anomaly_Accumulated_',num2str(year),'.mat'],...
%           'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon','dMLT_Qnet_mld_mon','dMLT_Qswr_mld_mon',...
%           'basin_mask_NA','Sxy_NA','lon_IAP','lat_IAP')
%     clear NA_SST_mon*  
%     count_yr=count_yr+1;
% end
% % #########################################################################
% % #########################################################################
% 
% 
% 
% % #########################################################################
% % % MLTa from IAP data
% % % #########################################################################
% %   % ACCESS Om2 0.25
%      % cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
%      load(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','Sxy_NA','lon','lat')
%        [lo,la]=meshgrid((260:380)', (0.5:59.5)');
%        basin_mask_NA=griddata(lon,lat,basin_mask_NA',lo',la','nearest');
%        clear lo la lon_025 lat_025
%        [Sxy,~,~]=function_Cgrid_Area_Distance((260:380)',(0.5:59.5)');
%        Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
%        Sxy_NA=Sxy; clear Sxy
% % #########################################################################
% 
% % #########################################################################
%     % Monthly MLT during 1981-2023
%     cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023')
%     mlt_month_19812023=nan(12,length(1981:2023)');
%     count=1;
%     for year=1981:2023
%           disp(['Monthly MLT Year# ',num2str(year)])
%           load(['mlt_monthly_Rh0_0125_GSW_IAP_C17_',num2str(year),'.mat'],'mlt_monthly','lon_IAP','lat_IAP')
%           mlt_monthly_NA=mlt_monthly(240:360,91:150,1:12);
%           clear mlt_monthly
% 
%           for month=1:size(mlt_monthly_NA,3)
%               mlt_monthly_NA0=squeeze(mlt_monthly_NA(:,:,month));
%               mlt_monthly_NA0(isnan(basin_mask_NA))=NaN;
%               mlt_monthly_NA(:,:,month)=mlt_monthly_NA0;
%               clear mlt_monthly_NA0
%           end
%           clear days
%           mlt_monthly_NA=mlt_monthly_NA.*Sxy_NA;
%           
%           % NA 0-60N
%           mlt_month_19812023(1:size(mlt_monthly_NA,3),count)=squeeze(nansum(nansum(mlt_monthly_NA(:,:,:),2),1))./squeeze(nansum(nansum(Sxy_NA(:,:),2),1));
%           clear mlt_monthly_NA
%           count=count+1;
%     end
%     clear year count
% 
%     % Anomaly is relative to the clim-daily mean
%     mlt_month_19812010_clim=squeeze(nanmean(mlt_month_19812023(:,1:30),2));
%     count=1;
%     for year=1:size(mlt_month_19812023,2)
%         mlt_month_19812023_ano(:,count)=mlt_month_19812023(:,count)-mlt_month_19812010_clim;
%         count=count+1;
%     end
%     dmlt_month_19812023_ano=mlt_month_19812023_ano(2:12,:)-mlt_month_19812023_ano(1:11,:);
% % #########################################################################
% 
% 
% 
% % #########################################################################
% figure('Color',[1 1 1]);  %create a new figure of white color background
% ixs = 0.060; ixe = 0.450;  ixd = 0.10; ixw = (1-ixs-ixe-0*ixd)/1;
% iys = 0.050; iye = 0.050;  iyd = 0.08; iyw = (1-iys-iye-1*iyd)/2;
% 
% %          [left            bottom      width height]
% pos{101}  = [ixs          iys+0.1*iyw+1*iyd   ixw 1.7*iyw]; 
% pos{102}  = [ixs          iys+0.0*iyw+0*iyd   ixw 1*iyw]; 
% 
% 
% % dc1                 = [0    , 0.447, 0.741]; % Blue-ish.
% % dc6                 = [0.301, 0.745, 0.933]; % Light Blue-ish.
% % dc2                 = [0.850, 0.325, 0.098]; % Orange/Red-ish.
% % dc7                 = [0.635, 0.078, 0.184]; % Red/Brown-ish.
% % dc5                 = [0.466, 0.674, 0.188]; % Green-ish.
% % dc3                 = [0.929, 0.694, 0.125]; % Yellow/Gold-ish.
% % dc4                 = [0.494, 0.184, 0.556]; % Purple-ish.
% 
% 
% clear color color0 
% color=cbrewer('seq', 'Blues', 60,'pchip');
% color(:,:)=color(60:-1:1,:);
% 
%     
% bar_dmlt(:,1)=dmlt_month_19812023_ano(1:8,43);
% bar_dmlt(:,2)=(dMLT_Qnet_mon(1:8,1)+dMLT_Qnet_mon(2:9,1))./2;
% bar_dmlt(:,3)=(dMLT_Qswr_mon(1:8,1)+dMLT_Qswr_mon(2:9,1))./2;
% bar_dmlt(:,4)=(dMLT_Qlat_mon(1:8,1)+dMLT_Qlat_mon(2:9,1))./2;
% bar_dmlt(:,5)=(dMLT_Qlon_mon(1:8,1)+dMLT_Qlon_mon(2:9,1))./2;
% bar_dmlt(:,6)=(dMLT_Qsen_mon(1:8,1)+dMLT_Qsen_mon(2:9,1))./2;
% bar_dmlt(:,7)=bar_dmlt(:,1)-bar_dmlt(:,2);
% 
% subplot('position',pos{101})
%     b=bar(4:8,bar_dmlt(4:8,:));
%     hold on
%     % set 3 display names for the 3 handles
%     set(b, {'DisplayName'}, {'MLT tendency','Qnet','SWR','Latent','Longwave','Sensible','Residual'}')
%     % Legend will show names for each color
%     legend() 
%     set(gca,'Ylim',[-0.2 0.2],'ycolor','k') 
%     set(gca,'YTick',-0.2:0.1:0.2)
%     set(gca,'YTickLabel',{'-0.2','-0.1','0','0.1','0.2'},'fontsize',16)
%     set(gca,'Xlim',[3.5 8.5]) 
%     set(gca,'XTick',4:1:8)
%     set(gca,'XTickLabel',{'Apr-May','May-Jun','Jun-Jul',...
%                           'Jul-Aug','Aug-Sep'},'fontsize',20)
% 
%     grid on
%     set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
%     ylabel(['[ ^oC per month ]'],'fontsize',20,'color','k','FontWeight','normal')
% 
% 
%         
% cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')




%% #########################################################################
%  #########################################################################
% 5. Accumulation of heat flux in 1981-2023
clc;clear
time_ann=(1981:2010)';

% #########################################################################
%   % ACCESS Om2 0.25
%     cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/access-om2-025_daily_sst_fields/')
%     load('basin_mask_NA_ACCESS-OM2-025_era5_iaf.mat','basin_mask_NA','lon_025','lat_025')
%     
%     cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
%     load(['SST_NA_100W20E_Eq70N_ERA5_Daily_Averaged_Reanalysis_on_Single_Levels_2022.mat'],'sst_daily_NA','lon','lat')
%        [lo,la]=meshgrid((lon)', (lat)');
%        basin_mask_NA=griddata(lon_025,lat_025,basin_mask_NA',lo',la','nearest');
%        clear lo la lon_025 lat_025
%        sst_daily_NA=nanmean(sst_daily_NA,3);
%        basin_mask_NA(isnan(sst_daily_NA))=NaN; clear sst_daily_NA
%     
%        [Sxy,~,~]=function_Cgrid_Area_Distance((lon)',(lat)');
%        Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
%        Sxy_NA=Sxy;
%      save(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','Sxy_NA','lon','lat')

     % cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
     load(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','lon','lat')
       [lo,la]=meshgrid((260:380)', (0.5:59.5)');
       basin_mask_NA=griddata(lon,lat,basin_mask_NA',lo',la','nearest');
       clear lo la lon lat
       [Sxy,~,~]=function_Cgrid_Area_Distance((260:380)',(0.5:59.5)');
       Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
       Sxy_NA=Sxy; clear Sxy
% #########################################################################


% #########################################################################
      % second_2_month=60*60*24*30;
      second_2_month(1,1)=60*60*24*31;
      second_2_month(2,1)=60*60*24*28.25;
      second_2_month(3,1)=60*60*24*31;
      second_2_month(4,1)=60*60*24*30;
      second_2_month(5,1)=60*60*24*31;
      second_2_month(6,1)=60*60*24*30;
      second_2_month(7,1)=60*60*24*31;
      second_2_month(8,1)=60*60*24*31;
      second_2_month(9,1)=60*60*24*30;
      second_2_month(10,1)=60*60*24*31;
      second_2_month(11,1)=60*60*24*30;
      second_2_month(12,1)=60*60*24*31;
      
      for month=1:12
          second_2_month0(:,:,month)=second_2_month(month,1).*ones(length(lon_025),length(lat_025));
      end
      second_2_month=second_2_month0; clear second_2_month0
% #########################################################################
      

% #########################################################################
count_yr=1;
for year=1981:2023
    disp([' Predicting NA daily SST Year#',num2str(year)])
    % #####################################################################
    % Monthly MLD from IAP
    disp(['   MLD from IAPv3, Year#',num2str(year)])
      cd('/Users/z5195509/Documents/Data/IAP_Ocean/Global_GSW_Monthly/MLD')
      % dRh0=0.125 ########################################################
      if year<2023
          load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld','lon_IAP','lat_IAP')
          lat_IAP=lat_IAP(:,1);
          mld0(:,:,1:12)=squeeze(mld(:,:,1,1:12));
          clear mld
      else
          % dRh0=0.125 ########################################################
          load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld','lon_IAP','lat_IAP')
          lat_IAP=lat_IAP(:,1);
          mld0(:,:,1:9)=squeeze(mld(:,:,1,1:9));
          clear mld
          load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAPv4_Temp_2023.mat'],'mld')
          mld0(:,:,10:12)=squeeze(mld(:,:,1,10:12));
          clear mld
      end
      
      lon_IAP0(1:340,1)=lon_IAP(21:360,1);
      lon_IAP0(341:360,1)=lon_IAP(1:20,1)+360; clear lon_IAP
      lon_IAP=lon_IAP0; clear lon_IAP0
      
      mld00(1:340,:,:)=mld0(21:360,:,:);
      mld00(341:360,:,:)=mld0(1:20,:,:); clear mld0
      mld_monthly=mld00; clear mld00
      % figure;imagesc(mld_monthly(:,:,1))
      % figure;plot(1:12,squeeze(mld_monthly(240,35,1:12)))
    % #####################################################################
    
        
    % Monthly heat flux terms
      cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg/')
      load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'_1deg.mat'],...
            'shortwave','latent','longwave','sensible') 
      Qnet=latent+longwave+sensible+shortwave;
      % figure;imagesc(Qnet(240:360,91:150,1))
        
      % ##################################################################
      % Remove the shortwave penatration
        % The Vertical Redistribution of SWR, at MLD=mld_2023_daily
        R=0.58; h1=0.35; h2=23; 
        z=mld_monthly;
        F_PS77=R.*exp(-z./h1)+(1-R).*exp(-z./h2);
        
        shortwave_mld = - shortwave.*F_PS77;
        clear R h1 h2 z F_PS77 
        % figure;imagesc(shortwave_mld(240:360,91:150,1))
        Qnet  = Qnet + shortwave_mld;
        shortwave = shortwave + shortwave_mld;
        clear shortwave_mld
      % ##################################################################   
       
      
      % ##################################################################   
      % Heat flux anomaly relative to 1981-2010 climatology 
%         Qnet      = Qnet      - Qnet_clim;
%         shortwave = shortwave - shortwave_clim;
%         latent    = latent    - latent_clim;
%         sensible  = sensible  - sensible_clim;
%         longwave  = longwave  - longwave_clim;
      % ##################################################################
      
       
      % ##################################################################
      % Mixed layer warming due to Qnet
        disp(['   Projected ML Warming in Year#',num2str(year)])
        % Heat flux term in heat budget
        MLT_Qnet_mon(:,:,:)  = Qnet./mld_monthly./3992./1027.*(24*60*60*30);  % K/month
        clear Qnet
        MLT_Qswr_mon(:,:,:)  = shortwave./mld_monthly./3992./1027.*(24*60*60*30); % K/month
        clear shortwave
        MLT_Qlat_mon(:,:,:)  = latent./mld_monthly./3992./1027.*(24*60*60*30);    % K/month
        clear latent
        MLT_Qlon_mon(:,:,:)  = longwave./mld_monthly./3992./1027.*(24*60*60*30);  % K/month
        clear longwave
        MLT_Qsen_mon(:,:,:)  = sensible./mld_monthly./3992./1027.*(24*60*60*30);  % K/month
        clear sensible
        clear mld_monthly
        % figure;imagesc(MLT_Qnet_mon(240:360,91:150,1)); caxis([-3 3])

        
         % ################################################################
         % NA (100W-20E,0-60N) averaged MLT tendency by heat flux terms 
         dMLT_Qnet_mon0=MLT_Qnet_mon(240:360,91:150,:).*Sxy_NA(:,:);
         dMLT_Qswr_mon0=MLT_Qswr_mon(240:360,91:150,:).*Sxy_NA(:,:);
         dMLT_Qlat_mon0=MLT_Qlat_mon(240:360,91:150,:).*Sxy_NA(:,:);
         dMLT_Qlon_mon0=MLT_Qlon_mon(240:360,91:150,:).*Sxy_NA(:,:);
         dMLT_Qsen_mon0=MLT_Qsen_mon(240:360,91:150,:).*Sxy_NA(:,:);
         % figure;imagesc(dMLT_Qnet_mon0(:,:,1)); % figure;imagesc(dMLT_Qswr_mld_mon0(:,:,1));
         % figure;imagesc(Sxy_NA(:,:,1)); 
         Sxy_NA0=Sxy_NA;
         for month=2:12
            Sxy_NA0(:,:,month)=Sxy_NA0(:,:,1); 
         end
         Sxy_NA0(isnan(dMLT_Qnet_mon0))=NaN; 
         dMLT_Qnet_mon0(isnan(Sxy_NA0))=NaN; 
         dMLT_Qswr_mon0(isnan(Sxy_NA0))=NaN; 
         dMLT_Qlat_mon0(isnan(Sxy_NA0))=NaN; 
         dMLT_Qlon_mon0(isnan(Sxy_NA0))=NaN; 
         dMLT_Qsen_mon0(isnan(Sxy_NA0))=NaN; 
         Sxy_NA0(isnan(dMLT_Qnet_mon0))=NaN; 
         % figure;imagesc(dMLT_Qnet_mon0(:,:,1)); 
         % figure;imagesc(Sxy_NA0(:,:,1)); 
         
         dMLT_Qnet_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qnet_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
         dMLT_Qswr_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qswr_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
         dMLT_Qlat_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qlat_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
         dMLT_Qlon_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qlon_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
         dMLT_Qsen_mon(1:12,1)=squeeze(nansum(nansum(dMLT_Qsen_mon0(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
         clear dMLT_*_mon0 Sxy_NA0
      % ###################################################################
         
      
        
      % ###################################################################
      % Predicted MLT = accummulated / intgrated heat flux anomalies
        project_sst_Qnet=nan(size(MLT_Qnet_mon));
        project_sst_Qswr=nan(size(MLT_Qnet_mon));
        project_sst_Qlat=nan(size(MLT_Qnet_mon));
        project_sst_Qlon=nan(size(MLT_Qnet_mon));
        project_sst_Qsen=nan(size(MLT_Qnet_mon));
        
        
        % ################################################################
%         % from SST
%         % Jan-1
%         cd('/Users/z5195509/Documents/Data/ERA-5/Monthly_Averaged_Reanalysis_on_Single_Levels/SST_SAT/')
%         load(['SST_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'.mat'],...
%               'sst')
%         % 20E-380E
%         heatflux0(1:1360,:,1)=sst(81:1440,:,1);
%         heatflux0(1361:1440,:,1)=sst(1:80,:,1); clear sst
%         sst_daily0=heatflux0; clear heatflux0
%         
%            [lo,la]=meshgrid((lon_IAP)',(lat_IAP)');
%            sst_daily(:,:,1)=griddata((20:0.25:379.75)',(-90:0.25:90)',squeeze(sst_daily0(:,:,1))',lo',la','linear');
%         
%         sst_daily(isnan(MLT_Qnet_mon(:,:,1)))=NaN;
%          % figure;imagesc(sst_daily(240:360,91:150,1)); 
%          % figure;imagesc(MLT_Qnet_mon(240:360,91:150,1)); 
%          % figure;imagesc(Sxy_NA(:,1:end,1)); 
        % ################################################################
         
        
        % ################################################################
        % from MLT Jan 15
        cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023')
        disp(['Monthly MLT Year# ',num2str(year)])
        load(['mlt_monthly_Rh0_0125_GSW_IAP_C17_',num2str(year),'.mat'],'mlt_monthly','lon_IAP','lat_IAP')
        sst_daily=mlt_monthly(:,:,1);
        clear mlt_monthly
        sst_daily(isnan(MLT_Qnet_mon(:,:,1)))=NaN;
        % ################################################################
        
        
        % ################################################################
        % Jan-15
        project_sst_Qnet(:,:,1)=sst_daily(:,:,1);
        project_sst_Qswr(:,:,1)=sst_daily(:,:,1);
        project_sst_Qlat(:,:,1)=sst_daily(:,:,1);
        project_sst_Qlon(:,:,1)=sst_daily(:,:,1);
        project_sst_Qsen(:,:,1)=sst_daily(:,:,1);
        
        % Projected SST from Feb-15
        for month=2:12
            disp(['   Projeted monthly MLT mon#',num2str(month),' year#',num2str(year)])
            project_sst_Qnet(:,:,month)=sst_daily(:,:,1)+squeeze(nansum(MLT_Qnet_mon(:,:,2:month),3));
            project_sst_Qswr(:,:,month)=sst_daily(:,:,1)+squeeze(nansum(MLT_Qswr_mon(:,:,2:month),3));
            project_sst_Qlat(:,:,month)=sst_daily(:,:,1)+squeeze(nansum(MLT_Qlat_mon(:,:,2:month),3));
            project_sst_Qlon(:,:,month)=sst_daily(:,:,1)+squeeze(nansum(MLT_Qlon_mon(:,:,2:month),3));
            project_sst_Qsen(:,:,month)=sst_daily(:,:,1)+squeeze(nansum(MLT_Qsen_mon(:,:,2:month),3));
        end
        clear sst_daily MLT_*_mon
        project_sst_Qnet(project_sst_Qnet==0)=NaN;
        project_sst_Qswr(project_sst_Qswr==0)=NaN;
        project_sst_Qlat(project_sst_Qlat==0)=NaN;
        project_sst_Qlon(project_sst_Qlon==0)=NaN;
        % figure;imagesc(project_sst_Qnet(240:360,91:150,1)); caxis([-3 3])
        

        % NA (100W-20E,0-60N) averaged projected SST 
         project_sst_Qnet=project_sst_Qnet(240:360,91:150,:).*Sxy_NA(:,1:end);
         project_sst_Qswr=project_sst_Qswr(240:360,91:150,:).*Sxy_NA(:,1:end);
         project_sst_Qlat=project_sst_Qlat(240:360,91:150,:).*Sxy_NA(:,1:end);
         project_sst_Qlon=project_sst_Qlon(240:360,91:150,:).*Sxy_NA(:,1:end);
         project_sst_Qsen=project_sst_Qsen(240:360,91:150,:).*Sxy_NA(:,1:end);
         % figure;imagesc(project_sst_Qnet(:,:,1)); 
         % figure;imagesc(Sxy_NA(:,1:end,1)); 
         Sxy_NA0=Sxy_NA;
         for month=2:12
            Sxy_NA0(:,:,month)=Sxy_NA0(:,:,1); 
         end
         Sxy_NA0(isnan(project_sst_Qnet))=NaN; 
         project_sst_Qnet(isnan(Sxy_NA0))=NaN; 
         project_sst_Qswr(isnan(Sxy_NA0))=NaN; 
         project_sst_Qlat(isnan(Sxy_NA0))=NaN; 
         project_sst_Qlon(isnan(Sxy_NA0))=NaN; 
         project_sst_Qsen(isnan(Sxy_NA0))=NaN; 
         Sxy_NA0(isnan(project_sst_Qnet))=NaN; 
         % figure;imagesc(project_sst_Qnet(:,:,1));  caxis([-1e11 1e11])
         % figure;imagesc(Sxy_NA0(:,1:end,1)); caxis([0 1e10])
         
         proj_NA_SST_mon_Qnet(1:12,count_yr)=squeeze(nansum(nansum(project_sst_Qnet(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
         proj_NA_SST_mon_Qswr(1:12,count_yr)=squeeze(nansum(nansum(project_sst_Qswr(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
         proj_NA_SST_mon_Qlat(1:12,count_yr)=squeeze(nansum(nansum(project_sst_Qlat(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
         proj_NA_SST_mon_Qlon(1:12,count_yr)=squeeze(nansum(nansum(project_sst_Qlon(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
         proj_NA_SST_mon_Qsen(1:12,count_yr)=squeeze(nansum(nansum(project_sst_Qsen(:,1:end,:),2),1))./squeeze(nansum(nansum(Sxy_NA0(:,1:end,:),2),1));
         clear project_sst_* Sxy_NA0
         % ################################################################
         
         
         % ################################################################
         % Saving every year
         NA_SST_mon_Qnet(1:12,1)=squeeze(proj_NA_SST_mon_Qnet(1:12,count_yr));
         NA_SST_mon_Qswr(1:12,1)=squeeze(proj_NA_SST_mon_Qswr(1:12,count_yr));
         NA_SST_mon_Qlat(1:12,1)=squeeze(proj_NA_SST_mon_Qlat(1:12,count_yr));
         NA_SST_mon_Qlon(1:12,1)=squeeze(proj_NA_SST_mon_Qlon(1:12,count_yr));
         NA_SST_mon_Qsen(1:12,1)=squeeze(proj_NA_SST_mon_Qsen(1:12,count_yr));
         
         cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
         save(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Accumulated_',num2str(year),'.mat'],...
              'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen',...
              'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon',...
              'basin_mask_NA','Sxy_NA','lon_IAP','lat_IAP')
          clear NA_SST_mon* dMLT_*_mon
         % ################################################################
      count_yr=count_yr+1;
end
clear *_dail_clim
% #########################################################################
% #########################################################################  



%% 5-1. Figure# Anomaly of Monthly MLT, SST, and MLT decomposition in 2023
% Part 1: Twitter plot for daily SSTA
clc;clear
time_ann=(1981:2022)';

% #########################################################################
% #########################################################################
% Mixed layer warming by air-sea flues
count_yr=1;
for year=1981:2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
    load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Accumulated_',num2str(year),'.mat'],...
          'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen',...
          'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon',...
          'basin_mask_NA','Sxy_NA','lon_IAP','lat_IAP')
    NA_MLT_Qnet(:,count_yr)=NA_SST_mon_Qnet; 
    NA_MLT_Qswr(:,count_yr)=NA_SST_mon_Qswr; 
    NA_MLT_Qlat(:,count_yr)=NA_SST_mon_Qlat; 
    NA_MLT_Qlon(:,count_yr)=NA_SST_mon_Qlon; 
    NA_MLT_Qsen(:,count_yr)=NA_SST_mon_Qsen; 
    clear NA_SST_mon*  
    count_yr=count_yr+1;
end
% Climatology
    NA_MLT_Qnet_clim=nanmean(NA_MLT_Qnet(:,1:30),2);
    NA_MLT_Qswr_clim=nanmean(NA_MLT_Qswr(:,1:30),2);
    NA_MLT_Qlat_clim=nanmean(NA_MLT_Qlat(:,1:30),2);
    NA_MLT_Qlon_clim=nanmean(NA_MLT_Qlon(:,1:30),2);
    NA_MLT_Qsen_clim=nanmean(NA_MLT_Qsen(:,1:30),2);

    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(NA_MLT_Qnet,2)
        NA_MLT_Qnet_ano(:,count)=NA_MLT_Qnet(:,count)-NA_MLT_Qnet_clim(:,1);
        NA_MLT_Qswr_ano(:,count)=NA_MLT_Qswr(:,count)-NA_MLT_Qswr_clim(:,1);
        NA_MLT_Qlat_ano(:,count)=NA_MLT_Qlat(:,count)-NA_MLT_Qlat_clim(:,1);
        NA_MLT_Qlon_ano(:,count)=NA_MLT_Qlon(:,count)-NA_MLT_Qlon_clim(:,1);
        NA_MLT_Qsen_ano(:,count)=NA_MLT_Qsen(:,count)-NA_MLT_Qsen_clim(:,1);
        count=count+1;
    end
    %clear *_clim
% #########################################################################
% #########################################################################




% #########################################################################
% SSTa from ERA5 reanalysis
% #########################################################################
%   % ACCESS Om2 0.25
    cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/access-om2-025_daily_NA_SST_fields/')
    % load('basin_mask_NA_ACCESS-OM2-025_era5_iaf_no_medi.mat','basin_mask_NA','lon_025','lat_025')
    load('basin_mask_NA_ACCESS-OM2-025_era5_iaf.mat','basin_mask_NA','lon_025','lat_025')
    
    cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
    load(['SST_NA_100W20E_Eq70N_ERA5_Daily_Averaged_Reanalysis_on_Single_Levels_2022.mat'],'sst_daily_NA','lon','lat')
       [lo,la]=meshgrid((lon)', (lat)');
       basin_mask_NA=griddata(lon_025,lat_025,basin_mask_NA',lo',la','nearest');
       clear lo la lon_025 lat_025
       sst_daily_NA=nanmean(sst_daily_NA,3);
       basin_mask_NA(isnan(sst_daily_NA))=NaN; clear sst_daily_NA
    
       [Sxy,~,~]=function_Cgrid_Area_Distance((lon)',(lat)');
       Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
% #########################################################################


% #########################################################################
    % Daily SST during 1981-2010
    cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
    sst_daily_19812022=nan(366,length(1981:2022)');
    count=1;
    for year=1981:2022
          disp(['Daily SST Year# ',num2str(year)])
          load(['SST_NA_100W20E_Eq70N_ERA5_Daily_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'.mat'],'sst_daily_NA','lon','lat')

          for days=1:size(sst_daily_NA,3)
              sst_daily_NA0=squeeze(sst_daily_NA(:,:,days));
              sst_daily_NA0(isnan(basin_mask_NA))=NaN;
              sst_daily_NA(:,:,days)=sst_daily_NA0;
              clear sst_daily_NA0
          end
          clear days
          sst_daily_NA=sst_daily_NA.*Sxy;
          
          % NA 0-60N
          sst_daily_19812022(1:size(sst_daily_NA,3),count)=squeeze(nansum(nansum(sst_daily_NA(:,1:241,:),2),1))./squeeze(nansum(nansum(Sxy(:,1:241),2),1));
          clear sst_daily_NA
          count=count+1;
    end
    clear year count

    
    % 2023 Daily SST Anomaly relative to the clim-mean
          disp(['Daily SST Year# 2023'])
          load('SST_NA_100W20E_Eq70N_ERA5_Daily_Averaged_Reanalysis_on_Single_Levels_2023.mat','sst_daily')
          sst_daily_NA=sst_daily; clear sst_daily

          for days=1:size(sst_daily_NA,3)
              sst_daily_NA0=squeeze(sst_daily_NA(:,:,days));
              sst_daily_NA0(isnan(basin_mask_NA))=NaN;
              sst_daily_NA(:,:,days)=sst_daily_NA0;
              clear sst_daily_NA0
          end
          clear days
          sst_daily_NA=sst_daily_NA.*Sxy;
          
          % NA 0-60N
          sst_daily_2023(1:365,1)=squeeze(nansum(nansum(sst_daily_NA(:,1:241,1:365),2),1))./squeeze(nansum(nansum(Sxy(:,1:241),2),1));
          sst_daily_2023(366,1)=NaN;
          clear sst_daily_NA

          
    % Anomaly is relative to the clim-daily mean
    sst_daily_19812010_clim=squeeze(nanmean(sst_daily_19812022(:,1:30),2));
    count=1;
    for year=1:size(sst_daily_19812022,2)
        sst_daily_19812022_ano(:,count)=sst_daily_19812022(:,count)-sst_daily_19812010_clim;
        count=count+1;
    end
    
    sst_daily_2023_ano(:,1)=sst_daily_2023(:,1)-sst_daily_19812010_clim; 
% #########################################################################




% #########################################################################
% % MLTa from IAP data
% % #########################################################################
%   % ACCESS Om2 0.25
     % cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
     load(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','Sxy_NA','lon','lat')
       [lo,la]=meshgrid((260:380)', (0.5:59.5)');
       basin_mask_NA=griddata(lon,lat,basin_mask_NA',lo',la','nearest');
       clear lo la lon_025 lat_025
       [Sxy,~,~]=function_Cgrid_Area_Distance((260:380)',(0.5:59.5)');
       Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
       Sxy_NA=Sxy; clear Sxy
% #########################################################################

% #########################################################################
    % Monthly MLT during 1981-2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023')
    mlt_month_19812023=nan(12,length(1981:2023)');
    count=1;
    for year=1981:2023
          disp(['Monthly MLT Year# ',num2str(year)])
          load(['mlt_monthly_Rh0_0125_GSW_IAP_C17_',num2str(year),'.mat'],'mlt_monthly','lon_IAP','lat_IAP')
          mlt_monthly_NA=mlt_monthly(240:360,91:150,:);
          clear mlt_monthly

          for month=1:size(mlt_monthly_NA,3)
              mlt_monthly_NA0=squeeze(mlt_monthly_NA(:,:,month));
              mlt_monthly_NA0(isnan(basin_mask_NA))=NaN;
              mlt_monthly_NA(:,:,month)=mlt_monthly_NA0;
              clear mlt_monthly_NA0
          end
          clear days
          mlt_monthly_NA=mlt_monthly_NA.*Sxy_NA;
          
          % NA 0-60N
          mlt_month_19812023(1:size(mlt_monthly_NA,3),count)=squeeze(nansum(nansum(mlt_monthly_NA(:,:,:),2),1))./squeeze(nansum(nansum(Sxy_NA(:,:),2),1));
          clear mlt_monthly_NA
          count=count+1;
    end
    clear year count

    % Anomaly is relative to the clim-daily mean
    mlt_month_19812010_clim=squeeze(nanmean(mlt_month_19812023(:,1:30),2));
    count=1;
    for year=1:size(mlt_month_19812023,2)
        mlt_month_19812023_ano(:,count)=mlt_month_19812023(:,count)-mlt_month_19812010_clim;
        count=count+1;
    end
% #########################################################################




%% ########################################################################
%  Plotting 5-1: Daily SSTA from ERA5 and MLTA by Heat Flux
clc;
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.060; ixe = 0.450;  ixd = 0.10; ixw = (1-ixs-ixe-0*ixd)/1;
iys = 0.050; iye = 0.050;  iyd = 0.08; iyw = (1-iys-iye-1*iyd)/2;

%          [left            bottom      width height]
pos{101}  = [ixs          iys+0.1*iyw+1*iyd   ixw 1.7*iyw]; 
pos{102}  = [ixs          iys+0.0*iyw+0*iyd   ixw 1*iyw]; 


% dc1                 = [0    , 0.447, 0.741]; % Blue-ish.
% dc6                 = [0.301, 0.745, 0.933]; % Light Blue-ish.
% dc2                 = [0.850, 0.325, 0.098]; % Orange/Red-ish.
% dc7                 = [0.635, 0.078, 0.184]; % Red/Brown-ish.
% dc5                 = [0.466, 0.674, 0.188]; % Green-ish.
% dc3                 = [0.929, 0.694, 0.125]; % Yellow/Gold-ish.
% dc4                 = [0.494, 0.184, 0.556]; % Purple-ish.


clear color color0 
color=cbrewer('seq', 'Blues', 60,'pchip');
color(:,:)=color(60:-1:1,:);

    
subplot('position',pos{101})
%     for year=1:42
%         line00=plot(1:366,sst_daily_19812022_ano(:,year));
%         set(line00,'color',color(year,:),'LineWidth',1.5,'linestyle','-'); 
%         hold on
%     end
%     
%     line00=plot(1:366,sst_daily_19812022_ano(:,5),'LineWidth',2);
% 
%     hold on
    line2023=plot(1:366,sst_daily_2023_ano(:,1));
    set(line2023,'color',[0.3,0.3,0.3],'LineWidth',3,'linestyle','-'); 
    
    hold on
    line2023MLT3=plot(16:30:270,NA_MLT_Qlat_ano(1:9,43));
    set(line2023MLT3,'color',[0.466, 0.674, 0.188],'LineWidth',3,'linestyle','-'); 
    hold on
    line2023MLT4=plot(16:30:270,NA_MLT_Qlon_ano(1:9,43));
    set(line2023MLT4,'color',[0.5,0.5,0.5],'LineWidth',3,'linestyle','-'); 
    hold on
    line2023MLT5=plot(16:30:270,NA_MLT_Qsen_ano(1:9,43));
    set(line2023MLT5,'color',[0.301, 0.745, 0.933],'LineWidth',3,'linestyle','-'); 
    hold on
    line2023MLT2=plot(16:30:270,NA_MLT_Qswr_ano(1:9,43));
    set(line2023MLT2,'color',[0.494, 0.184, 0.556],'LineWidth',3,'linestyle','-'); 
    hold on
    line2023MLT1=plot(16:30:270,NA_MLT_Qnet_ano(1:9,43));
    set(line2023MLT1,'color',[0.929, 0.694, 0.125],'LineWidth',6,'linestyle','-'); 
    
    hold on
    line2023MLT0=plot(16:30:270,mlt_month_19812023_ano(1:9,43));
    set(line2023MLT0,'color',[0.9, 0.3, 0.3],'LineWidth',6,'linestyle','-'); 

    
    leg101=legend([line2023 line2023MLT0 line2023MLT1 line2023MLT2 line2023MLT5 line2023MLT4 line2023MLT3],...
               '2023 SSTA, ERA5','2023 MLTA','2023 MLTA, Qnet','2023 MLTA, shortwave','2023 MLTA, sensible','2023 MLTA, longwave','2023 MLTA, latent heat',...
               'Location','northwest','NumColumns',1);
    set(leg101,'fontsize',20)
    hold on
    title(leg101,'North Atlantic SST and MLT Anomaly','fontsize',20')
    legend('boxoff')

    
    set(gca,'Ylim',[-0.0 1.5],'ycolor','k') 
    set(gca,'YTick',-1.5:0.5:4)
    set(gca,'YTickLabel',{'-1.5','-1.0','-0.5','0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0'},'fontsize',16)
    set(gca,'Xlim',[1 366]) 
    set(gca,'XTick',1:30:360)
    set(gca,'XTickLabel',{'           Jan','           Feb','           Mar','           Apr','           May','           Jun',...
                          '           Jul','           Aug','           Sep','           Oct','           Nov','            Dec'},'fontsize',20)

    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
    ylabel(['Difference [ ^oC ]'],'fontsize',20,'color','k','FontWeight','normal')

    % title('a. North Atlantic SST and predicted MLT anomaly','fontsize',24,'color','k','FontWeight','bold')

    
        
cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')




%% ########################################################################
%  Plotting 5-2: Bar Charts for Monthly MLTA and ML Warming Decomposition by Monthly Heat Flux
clc;clear
time_ann=(1981:2022)';

% #########################################################################
% #########################################################################
% Mixed layer warming by air-sea flues in 2023 - accumuated monthly heat flux 
count_yr=1;
for year=1981:2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
    load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Accumulated_',num2str(year),'.mat'],...
          'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen',...
          'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon')
    dMLT_Qnet_mon0(:,count_yr)=dMLT_Qnet_mon;
    dMLT_Qswr_mon0(:,count_yr)=dMLT_Qswr_mon;
    dMLT_Qlat_mon0(:,count_yr)=dMLT_Qlat_mon;
    dMLT_Qlon_mon0(:,count_yr)=dMLT_Qlon_mon;
    dMLT_Qsen_mon0(:,count_yr)=dMLT_Qsen_mon;
    clear NA_SST_mon*  
    count_yr=count_yr+1;
end
% Climatology
    NA_MLT_Qnet_clim=nanmean(dMLT_Qnet_mon0(:,1:30),2);
    NA_MLT_Qswr_clim=nanmean(dMLT_Qswr_mon0(:,1:30),2);
    NA_MLT_Qlat_clim=nanmean(dMLT_Qlat_mon0(:,1:30),2);
    NA_MLT_Qlon_clim=nanmean(dMLT_Qlon_mon0(:,1:30),2);
    NA_MLT_Qsen_clim=nanmean(dMLT_Qsen_mon0(:,1:30),2);

    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(dMLT_Qnet_mon0,2)
        NA_MLT_Qnet_ano(:,count)=dMLT_Qnet_mon0(:,count)-NA_MLT_Qnet_clim(:,1);
        NA_MLT_Qswr_ano(:,count)=dMLT_Qswr_mon0(:,count)-NA_MLT_Qswr_clim(:,1);
        NA_MLT_Qlat_ano(:,count)=dMLT_Qlat_mon0(:,count)-NA_MLT_Qlat_clim(:,1);
        NA_MLT_Qlon_ano(:,count)=dMLT_Qlon_mon0(:,count)-NA_MLT_Qlon_clim(:,1);
        NA_MLT_Qsen_ano(:,count)=dMLT_Qsen_mon0(:,count)-NA_MLT_Qsen_clim(:,1);
        count=count+1;
    end
    %clear *_clim
% #########################################################################
% #########################################################################



% #########################################################################
% % MLTa from IAP data
% % #########################################################################
%   % ACCESS Om2 0.25
     % cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
     load(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','Sxy_NA','lon','lat')
       [lo,la]=meshgrid((260:380)', (0.5:59.5)');
       basin_mask_NA=griddata(lon,lat,basin_mask_NA',lo',la','nearest');
       clear lo la lon_025 lat_025
       [Sxy,~,~]=function_Cgrid_Area_Distance((260:380)',(0.5:59.5)');
       Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
       Sxy_NA=Sxy; clear Sxy
% #########################################################################

% #########################################################################
    % Monthly MLT during 1981-2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023')
    mlt_month_19812023=nan(12,length(1981:2023)');
    count=1;
    for year=1981:2023
          disp(['Monthly MLT Year# ',num2str(year)])
          load(['mlt_monthly_Rh0_0125_GSW_IAP_C17_',num2str(year),'.mat'],'mlt_monthly','lon_IAP','lat_IAP')
          mlt_monthly_NA=mlt_monthly(240:360,91:150,1:12);
          clear mlt_monthly

          for month=1:size(mlt_monthly_NA,3)
              mlt_monthly_NA0=squeeze(mlt_monthly_NA(:,:,month));
              mlt_monthly_NA0(isnan(basin_mask_NA))=NaN;
              mlt_monthly_NA(:,:,month)=mlt_monthly_NA0;
              clear mlt_monthly_NA0
          end
          clear days
          mlt_monthly_NA=mlt_monthly_NA.*Sxy_NA;
          
          % NA 0-60N
          mlt_month_19812023(1:size(mlt_monthly_NA,3),count)=squeeze(nansum(nansum(mlt_monthly_NA(:,:,:),2),1))./squeeze(nansum(nansum(Sxy_NA(:,:),2),1));
          clear mlt_monthly_NA
          count=count+1;
    end
    clear year count

    % Anomaly is relative to the clim-daily mean
    mlt_month_19812010_clim=squeeze(nanmean(mlt_month_19812023(:,1:30),2));
    count=1;
    for year=1:size(mlt_month_19812023,2)
        mlt_month_19812023_ano(:,count)=mlt_month_19812023(:,count)-mlt_month_19812010_clim;
        count=count+1;
    end
    dmlt_month_19812023_ano=mlt_month_19812023_ano(2:12,:)-mlt_month_19812023_ano(1:11,:);
% #########################################################################



% #########################################################################
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.220; ixe = 0.220;  ixd = 0.10; ixw = (1-ixs-ixe-0*ixd)/1;
iys = 0.150; iye = 0.100;  iyd = 0.10; iyw = (1-iys-iye-1*iyd)/2;

%          [left            bottom      width height]
pos{101}  = [ixs          iys+1*iyw+1*iyd   ixw 1*iyw]; 

% dc1                 = [0    , 0.447, 0.741]; % Blue-ish.
% dc6                 = [0.301, 0.745, 0.933]; % Light Blue-ish.
% dc2                 = [0.850, 0.325, 0.098]; % Orange/Red-ish.
% dc7                 = [0.635, 0.078, 0.184]; % Red/Brown-ish.
% dc5                 = [0.466, 0.674, 0.188]; % Green-ish.
% dc3                 = [0.929, 0.694, 0.125]; % Yellow/Gold-ish.
% dc4                 = [0.494, 0.184, 0.556]; % Purple-ish.


clear color color0 
color=cbrewer('seq', 'Blues', 60,'pchip');
color(:,:)=color(60:-1:1,:);

% Centered at end-month
% bar_dmlt(:,1)=dmlt_month_19812023_ano(1:8,43);
% bar_dmlt(:,2)=(NA_MLT_Qnet_ano(1:8,43)+NA_MLT_Qnet_ano(2:9,43))./2;
% bar_dmlt(:,3)=(NA_MLT_Qswr_ano(1:8,43)+NA_MLT_Qswr_ano(2:9,43))./2;
% bar_dmlt(:,4)=(NA_MLT_Qlat_ano(1:8,43)+NA_MLT_Qlat_ano(2:9,43))./2;
% bar_dmlt(:,5)=(NA_MLT_Qlon_ano(1:8,43)+NA_MLT_Qlon_ano(2:9,43))./2;
% bar_dmlt(:,6)=(NA_MLT_Qsen_ano(1:8,43)+NA_MLT_Qsen_ano(2:9,43))./2;
% bar_dmlt(:,7)=bar_dmlt(:,1)-bar_dmlt(:,2);

% Centered at mid-month
bar_dmlt(2:8,1)=(dmlt_month_19812023_ano(1:7,43)+dmlt_month_19812023_ano(2:8,43))./2;
bar_dmlt(1:9,2)=NA_MLT_Qnet_ano(1:9,43);
bar_dmlt(1:9,3)=NA_MLT_Qswr_ano(1:9,43);
bar_dmlt(1:9,4)=NA_MLT_Qlat_ano(1:9,43);
bar_dmlt(1:9,5)=NA_MLT_Qlon_ano(1:9,43);
bar_dmlt(1:9,6)=NA_MLT_Qsen_ano(1:9,43);
bar_dmlt(1:9,7)=bar_dmlt(:,1)-bar_dmlt(:,2);

subplot('position',pos{101}) 
   h0=bar(5:8,bar_dmlt(5:8,:));
       set(h0,'BarWidth',0.94); 
       set(h0(1),'FaceColor',[0.850, 0.325, 0.098],'EdgeColor',[0.850, 0.325, 0.098])
       hold on
       set(h0(2),'FaceColor',[0.959, 0.494, 0.225],'EdgeColor',[0.959, 0.494, 0.225])
       hold on
       set(h0(3),'FaceColor',[0.929, 0.694, 0.125],'EdgeColor',[0.929, 0.694, 0.125])
       hold on
       set(h0(4),'FaceColor',[0    , 0.447, 0.641],'EdgeColor',[0    , 0.447, 0.641])
       hold on
       set(h0(5),'FaceColor',[0.101, 0.645, 0.833],'EdgeColor',[0.101, 0.645, 0.833])
       hold on
       set(h0(6),'FaceColor',[0.400, 0.850, 0.933],'EdgeColor',[0.400, 0.850, 0.933])
       hold on
       set(h0(7),'FaceColor',[0.594, 0.284, 0.556],'EdgeColor',[0.594, 0.284, 0.556])

    hold on
    legend([h0(1) h0(2) h0(3) h0(4) h0(5) h0(6) h0(7)],...
           'MLT tendency','Qnet','Shortwave','Latent','Longwave','Sensible','Residual',...
           'Location','northeast','Orientation','vertical','NumColumns',2)
    hold on
    set(legend,'fontsize',20)
    hold on
    % title(leg4,'Monthly SST Anomaly','fontsize',20')
    legend('boxoff')
    
    % Legend will show names for each color
    set(gca,'Ylim',[-0.45 0.8],'ycolor','k') 
    set(gca,'YTick',-0.4:0.4:1.2)
    set(gca,'YTickLabel',{'-0.4','0','0.4','0.8','1.2'},'fontsize',24)
    set(gca,'Xlim',[4.5 8.5]) 
    set(gca,'XTick',4.5:0.5:8.5)
    set(gca,'XTickLabel',{[],'May',[],'Jun',[],'Jul',[],'Aug',[]},'fontsize',20)
%     set(gca,'XTickLabel',{[],'May-Jun',[],'Jun-Jul',[],'Jul-Aug',[],'Aug-Sep',[]},'fontsize',20)
    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
    ylabel(['[ ^oC per month ]'],'fontsize',24,'color','k','FontWeight','normal')
    
        
cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')




%% ########################################################################
%  Plotting 5-2: Bar Charts for Monthly MLTA and ML Warming Decomposition by Monthly Heat Flux
% #########################################################################
%  V2: Adding Bars of MLD/Qnet Climatology
% #########################################################################
clc;clear
time_ann=(1981:2022)';


% #########################################################################
% #########################################################################
% 1. MLTa from IAP data
% #########################################################################
%   % ACCESS OM2 0.25
     % cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
     load(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','Sxy_NA','lon','lat')
       [lo,la]=meshgrid((260:380)', (0.5:59.5)');
       basin_mask_NA=griddata(lon,lat,basin_mask_NA',lo',la','nearest');
       clear lo la lon_025 lat_025
       [Sxy,~,~]=function_Cgrid_Area_Distance((260:380)',(0.5:59.5)');
       Sxy(isnan(basin_mask_NA))=NaN; 
       Sxy_NA=Sxy; clear Sxy % figure;imagesc(basin_mask_NA)
       
       
%      % cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
%      load(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','Sxy_NA','lon','lat')
%        [lo,la]=meshgrid((260:0.25:379.75)', (0:0.25:70)');
%        basin_mask_NA=griddata(lon,lat,basin_mask_NA',lo',la','nearest');
% 
%         figure('Color',[1 1 1]);  %create a new figure of white color background
%         ixs = 0.250; ixe = 0.250;  ixd = 0.10; ixw = (1-ixs-ixe-0*ixd)/1;
%         iys = 0.200; iye = 0.200;  iyd = 0.10; iyw = (1-iys-iye-0*iyd)/1;
% 
%         %          [left            bottom      width height]
%         pos{101}  = [ixs          iys+0*iyw+0*iyd   ixw 1*iyw]; 
% 
%         subplot('position',pos{101})
%         m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
%         m_pcolor(lon,lat,basin_mask_NA(:,:)');
%         shading flat
%         
%         hold on       
%         m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
%                     'xtick',[20:30:380],'xticklabels',[],...
%                     'yaxislocation','left','ytick',[-60:20:60]);
%         m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 


        
        
       
% #########################################################################

% #########################################################################
    % Monthly MLT during 1981-2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023')
    mlt_month_19812023=nan(12,length(1981:2023)');
    count=1;
    for year=1981:2023
          disp(['Monthly MLT Year# ',num2str(year)])
          load(['mlt_monthly_Rh0_0125_GSW_IAP_C17_',num2str(year),'.mat'],'mlt_monthly','lon_IAP','lat_IAP')
          mlt_monthly_NA=mlt_monthly(240:360,91:150,1:12);
          clear mlt_monthly

          for month=1:size(mlt_monthly_NA,3)
              mlt_monthly_NA0=squeeze(mlt_monthly_NA(:,:,month));
              mlt_monthly_NA0(isnan(basin_mask_NA))=NaN;
              mlt_monthly_NA(:,:,month)=mlt_monthly_NA0;
              clear mlt_monthly_NA0
          end
          clear month 
          mlt_monthly_NA=mlt_monthly_NA.*Sxy_NA;
          
          % NA 0-60N
          mlt_month_19812023(1:size(mlt_monthly_NA,3),count)=squeeze(nansum(nansum(mlt_monthly_NA(:,:,:),2),1))./squeeze(nansum(nansum(Sxy_NA(:,:),2),1));
          clear mlt_monthly_NA
          count=count+1;
    end
    clear year count Sxy_NA lon lat

    % Anomaly is relative to the clim-daily mean
    mlt_month_19812010_clim=squeeze(nanmean(mlt_month_19812023(:,1:30),2));
    count=1;
    for year=1:size(mlt_month_19812023,2)
        mlt_month_19812023_ano(:,count)=mlt_month_19812023(:,count)-mlt_month_19812010_clim;
        count=count+1;
    end
    dmlt_month_19812023_ano=mlt_month_19812023_ano(2:12,:)-mlt_month_19812023_ano(1:11,:);
    clear mlt_month* 
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
% 2.1 Mixed layer warming by air-sea flues in 2023 - accumuated monthly heat flux 
count_yr=1;
for year=1981:2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
    load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Accumulated_',num2str(year),'.mat'],...
          'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen',...
          'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon')
    dMLT_Qnet_mon0(:,count_yr)=dMLT_Qnet_mon;
    dMLT_Qswr_mon0(:,count_yr)=dMLT_Qswr_mon;
    dMLT_Qlat_mon0(:,count_yr)=dMLT_Qlat_mon;
    dMLT_Qlon_mon0(:,count_yr)=dMLT_Qlon_mon;
    dMLT_Qsen_mon0(:,count_yr)=dMLT_Qsen_mon;
    clear NA_SST_mon*  dMLT_*_mon
    count_yr=count_yr+1;
end
% Climatology
    NA_MLT_Qnet_clim=nanmean(dMLT_Qnet_mon0(:,1:30),2);
    NA_MLT_Qswr_clim=nanmean(dMLT_Qswr_mon0(:,1:30),2);
    NA_MLT_Qlat_clim=nanmean(dMLT_Qlat_mon0(:,1:30),2);
    NA_MLT_Qlon_clim=nanmean(dMLT_Qlon_mon0(:,1:30),2);
    NA_MLT_Qsen_clim=nanmean(dMLT_Qsen_mon0(:,1:30),2);

    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(dMLT_Qnet_mon0,2)
        NA_MLT_Qnet_ano(:,count)=dMLT_Qnet_mon0(:,count)-NA_MLT_Qnet_clim(:,1);
        NA_MLT_Qswr_ano(:,count)=dMLT_Qswr_mon0(:,count)-NA_MLT_Qswr_clim(:,1);
        NA_MLT_Qlat_ano(:,count)=dMLT_Qlat_mon0(:,count)-NA_MLT_Qlat_clim(:,1);
        NA_MLT_Qlon_ano(:,count)=dMLT_Qlon_mon0(:,count)-NA_MLT_Qlon_clim(:,1);
        NA_MLT_Qsen_ano(:,count)=dMLT_Qsen_mon0(:,count)-NA_MLT_Qsen_clim(:,1);
        count=count+1;
    end
    clear NA_MLT_*_clim dMLT_*_mon0
% #########################################################################
% #########################################################################


% #########################################################################
% #########################################################################
% 2.2 Mixed layer warming by air-sea flues in 2023 - test with MLD clim values
count_yr=1;
for year=1981:2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
    % Test with MLD clim values
    load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Accumulated_',num2str(year),'_2_Clim_MLD.mat'],...
          'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen',...
          'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon')
    dMLT_Qnet_mon0(:,count_yr)=dMLT_Qnet_mon;
    dMLT_Qswr_mon0(:,count_yr)=dMLT_Qswr_mon;
    dMLT_Qlat_mon0(:,count_yr)=dMLT_Qlat_mon;
    dMLT_Qlon_mon0(:,count_yr)=dMLT_Qlon_mon;
    dMLT_Qsen_mon0(:,count_yr)=dMLT_Qsen_mon;
    clear NA_SST_mon*  dMLT_*_mon
    count_yr=count_yr+1;
end
% Climatology
    NA_MLT_Qnet_clim=nanmean(dMLT_Qnet_mon0(:,1:30),2);
    NA_MLT_Qswr_clim=nanmean(dMLT_Qswr_mon0(:,1:30),2);
    NA_MLT_Qlat_clim=nanmean(dMLT_Qlat_mon0(:,1:30),2);
    NA_MLT_Qlon_clim=nanmean(dMLT_Qlon_mon0(:,1:30),2);
    NA_MLT_Qsen_clim=nanmean(dMLT_Qsen_mon0(:,1:30),2);

    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(dMLT_Qnet_mon0,2)
        NA_MLT_Qnet_ano_MLDclim(:,count)=dMLT_Qnet_mon0(:,count)-NA_MLT_Qnet_clim(:,1);
        NA_MLT_Qswr_ano_MLDclim(:,count)=dMLT_Qswr_mon0(:,count)-NA_MLT_Qswr_clim(:,1);
        NA_MLT_Qlat_ano_MLDclim(:,count)=dMLT_Qlat_mon0(:,count)-NA_MLT_Qlat_clim(:,1);
        NA_MLT_Qlon_ano_MLDclim(:,count)=dMLT_Qlon_mon0(:,count)-NA_MLT_Qlon_clim(:,1);
        NA_MLT_Qsen_ano_MLDclim(:,count)=dMLT_Qsen_mon0(:,count)-NA_MLT_Qsen_clim(:,1);
        count=count+1;
    end
    clear NA_MLT_*_clim dMLT_*_mon0
% #########################################################################
% #########################################################################


% #########################################################################
% #########################################################################
% 2.3 Mixed layer warming by air-sea flues in 2023 - test with Qnet clim values
count_yr=1;
for year=1981:2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
    % Test with Qnet clim values
    load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Accumulated_',num2str(year),'_3_Clim_Qnet.mat'],...
          'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen',...
          'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon')
    dMLT_Qnet_mon0(:,count_yr)=dMLT_Qnet_mon;
    dMLT_Qswr_mon0(:,count_yr)=dMLT_Qswr_mon;
    dMLT_Qlat_mon0(:,count_yr)=dMLT_Qlat_mon;
    dMLT_Qlon_mon0(:,count_yr)=dMLT_Qlon_mon;
    dMLT_Qsen_mon0(:,count_yr)=dMLT_Qsen_mon;
    clear NA_SST_mon*  dMLT_*_mon
    count_yr=count_yr+1;
end
% Climatology
    NA_MLT_Qnet_clim=nanmean(dMLT_Qnet_mon0(:,1:30),2);
    NA_MLT_Qswr_clim=nanmean(dMLT_Qswr_mon0(:,1:30),2);
    NA_MLT_Qlat_clim=nanmean(dMLT_Qlat_mon0(:,1:30),2);
    NA_MLT_Qlon_clim=nanmean(dMLT_Qlon_mon0(:,1:30),2);
    NA_MLT_Qsen_clim=nanmean(dMLT_Qsen_mon0(:,1:30),2);

    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(dMLT_Qnet_mon0,2)
        NA_MLT_Qnet_ano_Qclim(:,count)=dMLT_Qnet_mon0(:,count)-NA_MLT_Qnet_clim(:,1);
        NA_MLT_Qswr_ano_Qclim(:,count)=dMLT_Qswr_mon0(:,count)-NA_MLT_Qswr_clim(:,1);
        NA_MLT_Qlat_ano_Qclim(:,count)=dMLT_Qlat_mon0(:,count)-NA_MLT_Qlat_clim(:,1);
        NA_MLT_Qlon_ano_Qclim(:,count)=dMLT_Qlon_mon0(:,count)-NA_MLT_Qlon_clim(:,1);
        NA_MLT_Qsen_ano_Qclim(:,count)=dMLT_Qsen_mon0(:,count)-NA_MLT_Qsen_clim(:,1);
        count=count+1;
    end
    clear NA_MLT_*_clim dMLT_*_mon0
% #########################################################################
% #########################################################################


% #########################################################################
% #########################################################################
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.230; ixe = 0.230;  ixd = 0.10; ixw = (1-ixs-ixe-0*ixd)/1;
iys = 0.150; iye = 0.100;  iyd = 0.10; iyw = (1-iys-iye-1*iyd)/2;

%          [left            bottom      width height]
pos{101}  = [ixs          iys+1*iyw+1*iyd   ixw 1*iyw]; 

% dc1                 = [0    , 0.447, 0.741]; % Blue-ish.
% dc6                 = [0.301, 0.745, 0.933]; % Light Blue-ish.
% dc2                 = [0.850, 0.325, 0.098]; % Orange/Red-ish.
% dc7                 = [0.635, 0.078, 0.184]; % Red/Brown-ish.
% dc5                 = [0.466, 0.674, 0.188]; % Green-ish.
% dc3                 = [0.929, 0.694, 0.125]; % Yellow/Gold-ish.
% dc4                 = [0.494, 0.184, 0.556]; % Purple-ish.


clear color color0 
color=cbrewer('seq', 'Blues', 60,'pchip');
color(:,:)=color(60:-1:1,:);

% Centered at end-month
% bar_dmlt(:,1)=dmlt_month_19812023_ano(1:8,43);
% bar_dmlt(:,2)=(NA_MLT_Qnet_ano(1:8,43)+NA_MLT_Qnet_ano(2:9,43))./2;
% bar_dmlt(:,3)=(NA_MLT_Qswr_ano(1:8,43)+NA_MLT_Qswr_ano(2:9,43))./2;
% bar_dmlt(:,4)=(NA_MLT_Qlat_ano(1:8,43)+NA_MLT_Qlat_ano(2:9,43))./2;
% bar_dmlt(:,5)=(NA_MLT_Qlon_ano(1:8,43)+NA_MLT_Qlon_ano(2:9,43))./2;
% bar_dmlt(:,6)=(NA_MLT_Qsen_ano(1:8,43)+NA_MLT_Qsen_ano(2:9,43))./2;
% bar_dmlt(:,7)=bar_dmlt(:,1)-bar_dmlt(:,2);

% Centered at mid-month
bar_dmlt(2:8,1)=(dmlt_month_19812023_ano(1:7,43)+dmlt_month_19812023_ano(2:8,43))./2;
bar_dmlt(1:9,2)=NA_MLT_Qnet_ano(1:9,43);
bar_dmlt(1:9,3)=NA_MLT_Qswr_ano(1:9,43);
bar_dmlt(1:9,4)=NA_MLT_Qlat_ano(1:9,43);
bar_dmlt(1:9,5)=NA_MLT_Qlon_ano(1:9,43);
bar_dmlt(1:9,6)=NA_MLT_Qsen_ano(1:9,43);
bar_dmlt(1:9,7)=bar_dmlt(:,1)-bar_dmlt(:,2);


subplot('position',pos{101}) 
   h0=bar(5:8,bar_dmlt(5:8,:));
       set(h0,'BarWidth',0.94); 
       set(h0(1),'FaceColor',[0.850, 0.325, 0.098],'EdgeColor',[0.850, 0.325, 0.098])
       hold on
       set(h0(2),'FaceColor',[0.959, 0.494, 0.225],'EdgeColor',[0.959, 0.494, 0.225])
       hold on
       set(h0(3),'FaceColor',[0.929, 0.694, 0.125],'EdgeColor',[0.929, 0.694, 0.125])
       hold on
       set(h0(4),'FaceColor',[0    , 0.447, 0.641],'EdgeColor',[0    , 0.447, 0.641])
       hold on
       set(h0(5),'FaceColor',[0.101, 0.645, 0.833],'EdgeColor',[0.101, 0.645, 0.833])
       hold on
       set(h0(6),'FaceColor',[0.400, 0.850, 0.933],'EdgeColor',[0.400, 0.850, 0.933])
       hold on
       set(h0(7),'FaceColor',[0.594, 0.284, 0.556],'EdgeColor',[0.594, 0.284, 0.556])

       
       % ##################################################################
        % Bars and lines from MLD/Qnet climatology
        hold on
        bar_dmlt_Clim(1:9,1:14)=NaN;
        bar_dmlt_Clim(1:9,5)=NA_MLT_Qswr_ano_Qclim(1:9,43);
        bar_dmlt_Clim(1:9,6)=NA_MLT_Qswr_ano_MLDclim(1:9,43);
        h0_Clim=bar(5:8,bar_dmlt_Clim(5:8,:));
           set(h0_Clim,'BarWidth',0.36); 
           hold on
           set(h0_Clim(5),'FaceColor',[0.9, 0.2, 0.2],'EdgeColor',[0.9, 0.2, 0.2])
           hold on
           set(h0_Clim(6),'FaceColor',[0.2, 0.2, 0.2],'EdgeColor',[0.2, 0.2, 0.2])
           
        hold on
        plot((5.00:0.001:5.15),0.72*ones(length((5.05:0.001:5.20)')),'-','color',[0.929, 0.694, 0.125],'linewidth',18)
        hold on
        plot((5.00:0.001:5.15),0.59*ones(length((5.05:0.001:5.20)')),'-','color',[0.929, 0.694, 0.125],'linewidth',18)
        hold on
        plot((5.00:0.001:5.15),0.72*ones(length((5.05:0.001:5.20)')),'-','color',[0.9, 0.2, 0.2],'linewidth',4)
        hold on
        plot((5.00:0.001:5.15),0.59*ones(length((5.05:0.001:5.20)')),'-','color',[0.2, 0.2, 0.2],'linewidth',4)
        
        text(5.17, 0.71,'$\mathbf{\overline{Q_{sw}} / MLD^\prime}$', 'Interpreter', 'latex','fontsize',19,'FontName', 'Aptos')
        text(5.17, 0.58,'$\mathbf{{Q_{sw}}^\prime / \overline{MLD}}$', 'Interpreter', 'latex','fontsize',19,'FontName', 'Aptos')
       % ##################################################################
       
       
       
    % #####################################################################
    hold on
    leg1=legend([h0(1) h0(2) h0(3) h0(4) h0(5) h0(6) h0(7)],...
           'MLT tendency','Qnet','Shortwave','Latent','Longwave','Sensible','Residual',...
           'Location','northeast','Orientation','vertical','NumColumns',2);
    hold on
    set(leg1,'fontsize',19)
    hold on
    % title(leg4,'Monthly SST Anomaly','fontsize',20')
    legend('boxoff')
    
    % Second legend 
%     ah1=axes('position',get(gca,'position'),'visible','off');
%     hold on
%     leg2=legend(ah1,[h0_Clim(5) h0_Clim(6)],...
%            'Qnet_clim','MLD',...
%            'Location','northwest','Orientation','vertical','NumColumns',1);
%     hold on
%     set(leg2,'fontsize',18)
%     hold on
%     % title(leg4,'Monthly SST Anomaly','fontsize',20')
%     legend('boxoff')
    % #####################################################################

    
    % Legend will show names for each color
    set(gca,'Ylim',[-0.44 0.8],'ycolor','k') 
    set(gca,'YTick',-0.4:0.4:1.2)
    set(gca,'YTickLabel',{'-0.4','0','0.4','0.8','1.2'},'fontsize',22)
    set(gca,'Xlim',[4.5 8.5]) 
    set(gca,'XTick',5:1:8)
    set(gca,'XTickLabel',{'May','Jun','Jul','Aug'},'fontsize',22)
%     set(gca,'XTick',4.5:0.5:8.5)
%     set(gca,'XTickLabel',{[],'May',[],'Jun',[],'Jul',[],'Aug',[]},'fontsize',22)
%     set(gca,'XTickLabel',{[],'May-Jun',[],'Jun-Jul',[],'Jul-Aug',[],'Aug-Sep',[]},'fontsize',20)
    grid off
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
    ylabel(['[ ^oC per month ]'],'fontsize',22,'color','k','FontWeight','normal')
    
        
cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')



% #########################################################################
% #########################################################################
% Add all heat flux terms
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.230; ixe = 0.230;  ixd = 0.10; ixw = (1-ixs-ixe-0*ixd)/1;
iys = 0.150; iye = 0.100;  iyd = 0.10; iyw = (1-iys-iye-1*iyd)/2;

%          [left            bottom      width height]
pos{101}  = [ixs          iys+1*iyw+1*iyd   ixw 1*iyw]; 

% dc1                 = [0    , 0.447, 0.741]; % Blue-ish.
% dc6                 = [0.301, 0.745, 0.933]; % Light Blue-ish.
% dc2                 = [0.850, 0.325, 0.098]; % Orange/Red-ish.
% dc7                 = [0.635, 0.078, 0.184]; % Red/Brown-ish.
% dc5                 = [0.466, 0.674, 0.188]; % Green-ish.
% dc3                 = [0.929, 0.694, 0.125]; % Yellow/Gold-ish.
% dc4                 = [0.494, 0.184, 0.556]; % Purple-ish.


clear color color0 
color=cbrewer('seq', 'Blues', 60,'pchip');
color(:,:)=color(60:-1:1,:);

% Centered at end-month
% bar_dmlt(:,1)=dmlt_month_19812023_ano(1:8,43);
% bar_dmlt(:,2)=(NA_MLT_Qnet_ano(1:8,43)+NA_MLT_Qnet_ano(2:9,43))./2;
% bar_dmlt(:,3)=(NA_MLT_Qswr_ano(1:8,43)+NA_MLT_Qswr_ano(2:9,43))./2;
% bar_dmlt(:,4)=(NA_MLT_Qlat_ano(1:8,43)+NA_MLT_Qlat_ano(2:9,43))./2;
% bar_dmlt(:,5)=(NA_MLT_Qlon_ano(1:8,43)+NA_MLT_Qlon_ano(2:9,43))./2;
% bar_dmlt(:,6)=(NA_MLT_Qsen_ano(1:8,43)+NA_MLT_Qsen_ano(2:9,43))./2;
% bar_dmlt(:,7)=bar_dmlt(:,1)-bar_dmlt(:,2);

% Centered at mid-month
bar_dmlt(2:8,1)=(dmlt_month_19812023_ano(1:7,43)+dmlt_month_19812023_ano(2:8,43))./2;
bar_dmlt(1:9,2)=NA_MLT_Qnet_ano(1:9,43);
bar_dmlt(1:9,3)=NA_MLT_Qswr_ano(1:9,43);
bar_dmlt(1:9,4)=NA_MLT_Qlat_ano(1:9,43);
bar_dmlt(1:9,5)=NA_MLT_Qlon_ano(1:9,43);
bar_dmlt(1:9,6)=NA_MLT_Qsen_ano(1:9,43);
bar_dmlt(1:9,7)=bar_dmlt(:,1)-bar_dmlt(:,2);


subplot('position',pos{101}) 
   h0=bar(5:8,bar_dmlt(5:8,:));
       set(h0,'BarWidth',0.94); 
       set(h0(1),'FaceColor',[0.850, 0.325, 0.098],'EdgeColor',[0.850, 0.325, 0.098])
       hold on
       set(h0(2),'FaceColor',[0.959, 0.494, 0.225],'EdgeColor',[0.959, 0.494, 0.225])
       hold on
       set(h0(3),'FaceColor',[0.929, 0.694, 0.125],'EdgeColor',[0.929, 0.694, 0.125])
       hold on
       set(h0(4),'FaceColor',[0    , 0.447, 0.641],'EdgeColor',[0    , 0.447, 0.641])
       hold on
       set(h0(5),'FaceColor',[0.101, 0.645, 0.833],'EdgeColor',[0.101, 0.645, 0.833])
       hold on
       set(h0(6),'FaceColor',[0.400, 0.850, 0.933],'EdgeColor',[0.400, 0.850, 0.933])
       hold on
       set(h0(7),'FaceColor',[0.594, 0.284, 0.556],'EdgeColor',[0.594, 0.284, 0.556])

       
       % ##################################################################
        % Bars and lines from MLD/Qnet climatology
        hold on
        bar_dmlt_Clim(1:9,1:14)=NaN;
        bar_dmlt_Clim(1:9,3)=NA_MLT_Qnet_ano_Qclim(1:9,43);
        bar_dmlt_Clim(1:9,4)=NA_MLT_Qnet_ano_MLDclim(1:9,43);
        
        bar_dmlt_Clim(1:9,5)=NA_MLT_Qswr_ano_Qclim(1:9,43);
        bar_dmlt_Clim(1:9,6)=NA_MLT_Qswr_ano_MLDclim(1:9,43);
        
        bar_dmlt_Clim(1:9,7)=NA_MLT_Qlat_ano_Qclim(1:9,43);
        bar_dmlt_Clim(1:9,8)=NA_MLT_Qlat_ano_MLDclim(1:9,43);
        
        bar_dmlt_Clim(1:9,9) =NA_MLT_Qlon_ano_Qclim(1:9,43);
        bar_dmlt_Clim(1:9,10)=NA_MLT_Qlon_ano_MLDclim(1:9,43);
        
        bar_dmlt_Clim(1:9,11)=NA_MLT_Qsen_ano_Qclim(1:9,43);
        bar_dmlt_Clim(1:9,12)=NA_MLT_Qsen_ano_MLDclim(1:9,43);
        
        h0_Clim=bar(5:8,bar_dmlt_Clim(5:8,:));
           set(h0_Clim,'BarWidth',0.25); 
           hold on
           set(h0_Clim(3),'FaceColor',[0.9, 0.2, 0.2],'EdgeColor',[0.9, 0.2, 0.2])
           hold on
           set(h0_Clim(4),'FaceColor',[0.2, 0.2, 0.2],'EdgeColor',[0.2, 0.2, 0.2])
           hold on
           set(h0_Clim(5),'FaceColor',[0.9, 0.2, 0.2],'EdgeColor',[0.9, 0.2, 0.2])
           hold on
           set(h0_Clim(6),'FaceColor',[0.2, 0.2, 0.2],'EdgeColor',[0.2, 0.2, 0.2])
           hold on
           set(h0_Clim(7),'FaceColor',[0.9, 0.2, 0.2],'EdgeColor',[0.9, 0.2, 0.2])
           hold on
           set(h0_Clim(8),'FaceColor',[0.2, 0.2, 0.2],'EdgeColor',[0.2, 0.2, 0.2])
           hold on
           set(h0_Clim(9),'FaceColor',[0.9, 0.2, 0.2],'EdgeColor',[0.9, 0.2, 0.2])
           hold on
           set(h0_Clim(10),'FaceColor',[0.2, 0.2, 0.2],'EdgeColor',[0.2, 0.2, 0.2])
           hold on
           set(h0_Clim(11),'FaceColor',[0.9, 0.2, 0.2],'EdgeColor',[0.9, 0.2, 0.2])
           hold on
           set(h0_Clim(12),'FaceColor',[0.2, 0.2, 0.2],'EdgeColor',[0.2, 0.2, 0.2])
           
        hold on
        plot((5.00:0.001:5.15),0.72*ones(length((5.05:0.001:5.20)')),'-','color',[0.7, 0.7, 0.7],'linewidth',18)
        hold on
        plot((5.00:0.001:5.15),0.59*ones(length((5.05:0.001:5.20)')),'-','color',[0.7, 0.7, 0.7],'linewidth',18)
        hold on
        plot((5.00:0.001:5.15),0.72*ones(length((5.05:0.001:5.20)')),'-','color',[0.9, 0.2, 0.2],'linewidth',4)
        hold on
        plot((5.00:0.001:5.15),0.59*ones(length((5.05:0.001:5.20)')),'-','color',[0.2, 0.2, 0.2],'linewidth',4)
        
        text(5.17, 0.71,'$\mathbf{\overline{Q} / MLD^\prime}$', 'Interpreter', 'latex','fontsize',19,'FontName', 'Aptos','FontWeight','bold')
        text(5.17, 0.58,'$\mathbf{{Q}^\prime / \overline{MLD}}$', 'Interpreter', 'latex','fontsize',19,'FontName', 'Aptos','FontWeight','bold')
       % ##################################################################
       
       
       
    % #####################################################################
    hold on
    leg1=legend([h0(1) h0(2) h0(3) h0(4) h0(5) h0(6) h0(7)],...
           'MLT tendency','Qnet','Shortwave','Latent','Longwave','Sensible','Residual',...
           'Location','northeast','Orientation','vertical','NumColumns',2);
    hold on
    set(leg1,'fontsize',19)
    hold on
    % title(leg4,'Monthly SST Anomaly','fontsize',20')
    legend('boxoff')
    
    % Second legend 
%     ah1=axes('position',get(gca,'position'),'visible','off');
%     hold on
%     leg2=legend(ah1,[h0_Clim(5) h0_Clim(6)],...
%            'Qnet_clim','MLD',...
%            'Location','northwest','Orientation','vertical','NumColumns',1);
%     hold on
%     set(leg2,'fontsize',18)
%     hold on
%     % title(leg4,'Monthly SST Anomaly','fontsize',20')
%     legend('boxoff')
    % #####################################################################

    
    % Legend will show names for each color
    set(gca,'Ylim',[-0.44 0.8],'ycolor','k') 
    set(gca,'YTick',-0.4:0.4:1.2)
    set(gca,'YTickLabel',{'-0.4','0','0.4','0.8','1.2'},'fontsize',22)
    set(gca,'Xlim',[4.5 8.5]) 
    set(gca,'XTick',5:1:8)
    set(gca,'XTickLabel',{'May','Jun','Jul','Aug'},'fontsize',22)
%     set(gca,'XTick',4.5:0.5:8.5)
%     set(gca,'XTickLabel',{[],'May',[],'Jun',[],'Jul',[],'Aug',[]},'fontsize',22)
%     set(gca,'XTickLabel',{[],'May-Jun',[],'Jun-Jul',[],'Jul-Aug',[],'Aug-Sep',[]},'fontsize',20)
    grid off
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
    ylabel(['[ ^oC per month ]'],'fontsize',22,'color','k','FontWeight','normal')
    
        
cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')



%% ########################################################################
%  Plotting 5-2: Bar Charts for Monthly MLTA and ML Warming Decomposition by Monthly Heat Flux
% #########################################################################
%  V3: Adding Bars of MLD/Qsw/Qsw_h Climatology, new three-way decomposition
% #########################################################################
clc;clear
time_ann=(1981:2022)';


% #########################################################################
% #########################################################################
% 1. MLTa from IAP data
% #########################################################################
%   % ACCESS OM2 0.25
     % cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
     load(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','Sxy_NA','lon','lat')
       [lo,la]=meshgrid((260:380)', (0.5:59.5)');
       basin_mask_NA=griddata(lon,lat,basin_mask_NA',lo',la','nearest');
       clear lo la lon_025 lat_025
       [Sxy,~,~]=function_Cgrid_Area_Distance((260:380)',(0.5:59.5)');
       Sxy(isnan(basin_mask_NA))=NaN; 
       Sxy_NA=Sxy; clear Sxy % figure;imagesc(basin_mask_NA)
% #########################################################################



% #########################################################################
    % Monthly MLT during 1981-2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023')
    mlt_month_19812023=nan(12,length(1981:2023)');
    count=1;
    for year=1981:2023
          disp(['Monthly MLT Year# ',num2str(year)])
          load(['mlt_monthly_Rh0_0125_GSW_IAP_C17_',num2str(year),'.mat'],'mlt_monthly','lon_IAP','lat_IAP')
          mlt_monthly_NA=mlt_monthly(240:360,91:150,1:12);
          clear mlt_monthly

          for month=1:size(mlt_monthly_NA,3)
              mlt_monthly_NA0=squeeze(mlt_monthly_NA(:,:,month));
              mlt_monthly_NA0(isnan(basin_mask_NA))=NaN;
              mlt_monthly_NA(:,:,month)=mlt_monthly_NA0;
              clear mlt_monthly_NA0
          end
          clear month 
          mlt_monthly_NA=mlt_monthly_NA.*Sxy_NA;
          
          % NA 0-60N
          mlt_month_19812023(1:size(mlt_monthly_NA,3),count)=squeeze(nansum(nansum(mlt_monthly_NA(:,:,:),2),1))./squeeze(nansum(nansum(Sxy_NA(:,:),2),1));
          clear mlt_monthly_NA
          count=count+1;
    end
    clear year count Sxy_NA lon lat

    % Anomaly is relative to the clim-daily mean
    mlt_month_19812010_clim=squeeze(nanmean(mlt_month_19812023(:,1:30),2));
    count=1;
    for year=1:size(mlt_month_19812023,2)
        mlt_month_19812023_ano(:,count)=mlt_month_19812023(:,count)-mlt_month_19812010_clim;
        count=count+1;
    end
    dmlt_month_19812023_ano=mlt_month_19812023_ano(2:12,:)-mlt_month_19812023_ano(1:11,:);
    clear mlt_month* 
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
% 2.1 MLT budhet in 1981-2023
count_yr=1;
for year=1981:2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
    load(['MLT_budget_NA_monthly_ERA5_IAP_V2_1deg_V3_3way_decomposition_',num2str(year),'.mat'],...
          'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen',...
          'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon')
    dMLT_Qnet_mon0(:,count_yr)=dMLT_Qnet_mon;
    dMLT_Qswr_mon0(:,count_yr)=dMLT_Qswr_mon;
    dMLT_Qlat_mon0(:,count_yr)=dMLT_Qlat_mon;
    dMLT_Qlon_mon0(:,count_yr)=dMLT_Qlon_mon;
    dMLT_Qsen_mon0(:,count_yr)=dMLT_Qsen_mon;
    clear NA_SST_mon*  dMLT_*_mon
    count_yr=count_yr+1;
end
% Climatology
    NA_MLT_Qnet_clim=nanmean(dMLT_Qnet_mon0(:,1:30),2);
    NA_MLT_Qswr_clim=nanmean(dMLT_Qswr_mon0(:,1:30),2);
    NA_MLT_Qlat_clim=nanmean(dMLT_Qlat_mon0(:,1:30),2);
    NA_MLT_Qlon_clim=nanmean(dMLT_Qlon_mon0(:,1:30),2);
    NA_MLT_Qsen_clim=nanmean(dMLT_Qsen_mon0(:,1:30),2);

    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(dMLT_Qnet_mon0,2)
        NA_MLT_Qnet_ano(:,count)=dMLT_Qnet_mon0(:,count)-NA_MLT_Qnet_clim(:,1);
        NA_MLT_Qswr_ano(:,count)=dMLT_Qswr_mon0(:,count)-NA_MLT_Qswr_clim(:,1);
        NA_MLT_Qlat_ano(:,count)=dMLT_Qlat_mon0(:,count)-NA_MLT_Qlat_clim(:,1);
        NA_MLT_Qlon_ano(:,count)=dMLT_Qlon_mon0(:,count)-NA_MLT_Qlon_clim(:,1);
        NA_MLT_Qsen_ano(:,count)=dMLT_Qsen_mon0(:,count)-NA_MLT_Qsen_clim(:,1);
        count=count+1;
    end
    clear NA_MLT_*_clim dMLT_*_mon0
% #########################################################################
% #########################################################################


% #########################################################################
% #########################################################################
% 2.2 MLT Budget in 1981-2023 - using MLD and SWR climatology
count_yr=1;
for year=1981:2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
    % Test with MLD and SWR climatology
    load(['MLT_budget_NA_monthly_ERA5_IAP_V2_1deg_V3_3way_decomposition_',num2str(year),'_2_Clim_MLD_SWR.mat'],...
          'dMLT_Qswr_mon')
    dMLT_Qswr_mon0(:,count_yr)=dMLT_Qswr_mon;
    clear NA_SST_mon*  dMLT_*_mon
    count_yr=count_yr+1;
end
% Climatology
    NA_MLT_Qswr_clim=nanmean(dMLT_Qswr_mon0(:,1:30),2);
    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(dMLT_Qswr_mon0,2)
        NA_MLT_Qswr_ano_2_clim_MLD_SWR(:,count)=dMLT_Qswr_mon0(:,count)-NA_MLT_Qswr_clim(:,1);
        count=count+1;
    end
    clear NA_MLT_*_clim dMLT_*_mon0
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
% 2.3 MLT Budget in 1981-2023 - using MLD and Qsw_H climatology
count_yr=1;
for year=1981:2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
    % Test with MLD and Qsw_H climatology
    load(['MLT_budget_NA_monthly_ERA5_IAP_V2_1deg_V3_3way_decomposition_',num2str(year),'_3_Clim_MLD_Qsw_H.mat'],...
          'dMLT_Qswr_mon')
    dMLT_Qswr_mon0(:,count_yr)=dMLT_Qswr_mon;
    clear NA_SST_mon*  dMLT_*_mon
    count_yr=count_yr+1;
end
% Climatology
    NA_MLT_Qswr_clim=nanmean(dMLT_Qswr_mon0(:,1:30),2);
    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(dMLT_Qswr_mon0,2)
        NA_MLT_Qswr_ano_3_clim_MLD_Qsw_H(:,count)=dMLT_Qswr_mon0(:,count)-NA_MLT_Qswr_clim(:,1);
        count=count+1;
    end
    clear NA_MLT_*_clim dMLT_*_mon0
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
% 2.4 MLT Budget in 1981-2023 - using SWR and Qsw_H climatology
count_yr=1;
for year=1981:2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
    % Test with MLD and Qsw_H climatology
    load(['MLT_budget_NA_monthly_ERA5_IAP_V2_1deg_V3_3way_decomposition_',num2str(year),'_4_Clim_SWR_Qsw_H.mat'],...
          'dMLT_Qswr_mon')
    dMLT_Qswr_mon0(:,count_yr)=dMLT_Qswr_mon;
    clear NA_SST_mon*  dMLT_*_mon
    count_yr=count_yr+1;
end
% Climatology
    NA_MLT_Qswr_clim=nanmean(dMLT_Qswr_mon0(:,1:30),2);
    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(dMLT_Qswr_mon0,2)
        NA_MLT_Qswr_ano_4_clim_SWR_Qsw_H(:,count)=dMLT_Qswr_mon0(:,count)-NA_MLT_Qswr_clim(:,1);
        count=count+1;
    end
    clear NA_MLT_*_clim dMLT_*_mon0
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.230; ixe = 0.230;  ixd = 0.10; ixw = (1-ixs-ixe-0*ixd)/1;
iys = 0.150; iye = 0.100;  iyd = 0.10; iyw = (1-iys-iye-1*iyd)/2;

%          [left            bottom      width height]
pos{101}  = [ixs          iys+1*iyw+1*iyd   ixw 1*iyw]; 

% dc1                 = [0    , 0.447, 0.741]; % Blue-ish.
% dc6                 = [0.301, 0.745, 0.933]; % Light Blue-ish.
% dc2                 = [0.850, 0.325, 0.098]; % Orange/Red-ish.
% dc7                 = [0.635, 0.078, 0.184]; % Red/Brown-ish.
% dc5                 = [0.466, 0.674, 0.188]; % Green-ish.
% dc3                 = [0.929, 0.694, 0.125]; % Yellow/Gold-ish.
% dc4                 = [0.494, 0.184, 0.556]; % Purple-ish.


clear color color0 
color=cbrewer('seq', 'Blues', 60,'pchip');
color(:,:)=color(60:-1:1,:);

% Centered at mid-month
bar_dmlt(2:8,1)=(dmlt_month_19812023_ano(1:7,43)+dmlt_month_19812023_ano(2:8,43))./2;
bar_dmlt(1:9,2)=NA_MLT_Qnet_ano(1:9,43);
bar_dmlt(1:9,3)=NA_MLT_Qswr_ano(1:9,43);
bar_dmlt(1:9,4)=NA_MLT_Qlat_ano(1:9,43);
bar_dmlt(1:9,5)=NA_MLT_Qlon_ano(1:9,43);
bar_dmlt(1:9,6)=NA_MLT_Qsen_ano(1:9,43);
bar_dmlt(1:9,7)=bar_dmlt(:,1)-bar_dmlt(:,2);


subplot('position',pos{101}) 
   h0=bar(5:8,bar_dmlt(5:8,:));
       set(h0,'BarWidth',0.94); 
       set(h0(1),'FaceColor',[0.850, 0.325, 0.098],'EdgeColor',[0.850, 0.325, 0.098])
       hold on
       set(h0(2),'FaceColor',[0.959, 0.494, 0.225],'EdgeColor',[0.959, 0.494, 0.225])
       hold on
       set(h0(3),'FaceColor',[0.929, 0.694, 0.125],'EdgeColor',[0.929, 0.694, 0.125])
       hold on
       set(h0(4),'FaceColor',[0    , 0.447, 0.641],'EdgeColor',[0    , 0.447, 0.641])
       hold on
       set(h0(5),'FaceColor',[0.101, 0.645, 0.833],'EdgeColor',[0.101, 0.645, 0.833])
       hold on
       set(h0(6),'FaceColor',[0.400, 0.850, 0.933],'EdgeColor',[0.400, 0.850, 0.933])
       hold on
       set(h0(7),'FaceColor',[0.594, 0.284, 0.556],'EdgeColor',[0.594, 0.284, 0.556])

       
       % ##################################################################
        % Bars and lines from MLD/Qsw_H/SWR climatology
        hold on
        bar_dmlt_Clim(1:9,1:21)=NaN;
        % SWR and Q_h cimatology
        bar_dmlt_Clim(1:9,7)=NA_MLT_Qswr_ano_4_clim_SWR_Qsw_H(1:9,43);       
        % MLD and Q_h cimatology
        bar_dmlt_Clim(1:9,8)=NA_MLT_Qswr_ano_3_clim_MLD_Qsw_H(1:9,43);
        % MLD and Q_swr cimatology
        bar_dmlt_Clim(1:9,9)=NA_MLT_Qswr_ano_2_clim_MLD_SWR(1:9,43);  
        
        h0_Clim=bar(5:8,bar_dmlt_Clim(5:8,:));
           set(h0_Clim,'BarWidth',0.68); 
           hold on
           set(h0_Clim(7),'FaceColor',[0.9, 0.2, 0.2],'EdgeColor',[0.9, 0.2, 0.2])
           hold on
           set(h0_Clim(8),'FaceColor',[0.2, 0.3, 0.2],'EdgeColor',[0.2, 0.3, 0.2])
           hold on
           set(h0_Clim(9),'FaceColor',[0.2, 0.8, 0.6],'EdgeColor',[0.2, 0.8, 0.6])
           
%         hold on
%         plot((5.00:0.001:5.15),0.72*ones(length((5.05:0.001:5.20)')),'-','color',[0.929, 0.694, 0.125],'linewidth',18)
%         hold on
%         plot((5.00:0.001:5.15),0.59*ones(length((5.05:0.001:5.20)')),'-','color',[0.929, 0.694, 0.125],'linewidth',18)
%         hold on
%         plot((5.00:0.001:5.15),0.72*ones(length((5.05:0.001:5.20)')),'-','color',[0.9, 0.2, 0.2],'linewidth',4)
%         hold on
%         plot((5.00:0.001:5.15),0.59*ones(length((5.05:0.001:5.20)')),'-','color',[0.2, 0.2, 0.2],'linewidth',4)
%         
%         text(5.17, 0.71,'$\mathbf{\overline{Q_{sw}} / MLD^\prime}$', 'Interpreter', 'latex','fontsize',19,'FontName', 'Aptos')
%         text(5.17, 0.58,'$\mathbf{{Q_{sw}}^\prime / \overline{MLD}}$', 'Interpreter', 'latex','fontsize',19,'FontName', 'Aptos')
%        % ##################################################################

           hold on
           plot((6.10:0.001:6.23),0.70*ones(length((6.10:0.001:6.23)')),'-','color',[0.929, 0.694, 0.125],'linewidth',18)
           hold on
           plot((6.10:0.001:6.23),0.59*ones(length((6.10:0.001:6.23)')),'-','color',[0.929, 0.694, 0.125],'linewidth',18)
           hold on
           plot((6.10:0.001:6.23),0.48*ones(length((6.10:0.001:6.23)')),'-','color',[0.929, 0.694, 0.125],'linewidth',18)  
           hold on
           plot((6.10:0.001:6.23),0.70*ones(length((6.10:0.001:6.23)')),'-','color',[0.9, 0.2, 0.2],'linewidth',4)
           hold on
           plot((6.10:0.001:6.23),0.59*ones(length((6.10:0.001:6.23)')),'-','color',[0.2, 0.2, 0.2],'linewidth',4)
           hold on
           plot((6.10:0.001:6.23),0.48*ones(length((6.10:0.001:6.23)')),'-','color',[0.2, 0.8, 0.6],'linewidth',4)  

           text(6.25, 0.70,'$\mathbf{MLD^\prime}$', 'Interpreter', 'latex','fontsize',18,'FontName', 'Aptos')
           text(6.25, 0.59,'$\mathbf{{Q_{sw}}^\prime}$', 'Interpreter', 'latex','fontsize',18,'FontName', 'Aptos')
           text(6.25, 0.48,'$\mathbf{{Q_{sw,H}}^\prime}$', 'Interpreter', 'latex','fontsize',18,'FontName', 'Aptos')
       % ##################################################################      
       
       
    % #####################################################################
    hold on
    leg1=legend([h0(1) h0(2) h0(3) h0(4) h0(5) h0(6) h0(7)],...
           'MLT tendency','Qnet','Shortwave','Latent','Longwave','Sensible','Residual',...
           'Location','northeast','Orientation','vertical','NumColumns',2);
    hold on
    set(leg1,'fontsize',18)
    hold on
    % title(leg4,'Monthly SST Anomaly','fontsize',20')
    legend('boxoff')
    
    % Second legend 
%     ah1=axes('position',get(gca,'position'),'visible','off');
%     hold on
%     leg2=legend(ah1,[h0_Clim(5) h0_Clim(6)],...
%            'Qnet_clim','MLD',...
%            'Location','northwest','Orientation','vertical','NumColumns',1);
%     hold on
%     set(leg2,'fontsize',18)
%     hold on
%     % title(leg4,'Monthly SST Anomaly','fontsize',20')
%     legend('boxoff')
    % #####################################################################

    
%     % Legend will show names for each color
%     set(gca,'Ylim',[-0.44 0.8],'ycolor','k') 
%     set(gca,'YTick',-0.4:0.4:1.2)
%     set(gca,'YTickLabel',{'-0.4','0','0.4','0.8','1.2'},'fontsize',22)
%     set(gca,'Xlim',[4.5 8.5]) 
%     set(gca,'XTick',5:1:8)
%     set(gca,'XTickLabel',{'May','Jun','Jul','Aug'},'fontsize',22)
% %     set(gca,'XTick',4.5:0.5:8.5)
% %     set(gca,'XTickLabel',{[],'May',[],'Jun',[],'Jul',[],'Aug',[]},'fontsize',22)
% %     set(gca,'XTickLabel',{[],'May-Jun',[],'Jun-Jul',[],'Jul-Aug',[],'Aug-Sep',[]},'fontsize',20)

    % Legend will show names for each color
    set(gca,'Ylim',[-0.44 0.8],'ycolor','k') 
    set(gca,'YTick',-0.4:0.4:1.6)
    set(gca,'YTickLabel',{'-0.4','0','0.4','0.8','1.2','1.6'},'fontsize',21)
    set(gca,'Xlim',[4.5 8.5]) 
    set(gca,'XTick',4.5:1:8.5)
    set(gca,'XTickLabel',{[],[],[],[]},'fontsize',21)
%     set(gca,'XTick',5:1:8)
%     set(gca,'XTickLabel',{'May','Jun','Jul','Aug'},'fontsize',20)
    text(4.915,-0.54,'May','fontsize',21,'color','k','FontWeight','normal')
    text(5.930,-0.54,'Jun','fontsize',21,'color','k','FontWeight','normal')
    text(6.945,-0.54,'Jul','fontsize',21,'color','k','FontWeight','normal')
    text(7.920,-0.54,'Aug','fontsize',21,'color','k','FontWeight','normal')

    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
    ylabel(['[ ^oC per month ]'],'fontsize',21,'color','k','FontWeight','normal')
    

cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')




% #########################################################################
% #########################################################################
%% 5.3 Combined 5-1 and 5-2
% 5-1. Figure# Anomaly of Monthly MLT, SST, and MLT decomposition in 2023
% Part 1: Twitter plot for daily SSTA
clc;clear
time_ann=(1981:2022)';

% #########################################################################
% #########################################################################
% Mixed layer warming by air-sea flues
count_yr=1;
for year=1981:2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
    load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Accumulated_',num2str(year),'.mat'],...
          'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen',...
          'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon',...
          'basin_mask_NA','Sxy_NA','lon_IAP','lat_IAP')
    NA_MLT_Qnet(:,count_yr)=NA_SST_mon_Qnet; 
    NA_MLT_Qswr(:,count_yr)=NA_SST_mon_Qswr; 
    NA_MLT_Qlat(:,count_yr)=NA_SST_mon_Qlat; 
    NA_MLT_Qlon(:,count_yr)=NA_SST_mon_Qlon; 
    NA_MLT_Qsen(:,count_yr)=NA_SST_mon_Qsen; 
    clear NA_SST_mon*  
    count_yr=count_yr+1;
end
% Climatology
    NA_MLT_Qnet_clim=nanmean(NA_MLT_Qnet(:,1:30),2);
    NA_MLT_Qswr_clim=nanmean(NA_MLT_Qswr(:,1:30),2);
    NA_MLT_Qlat_clim=nanmean(NA_MLT_Qlat(:,1:30),2);
    NA_MLT_Qlon_clim=nanmean(NA_MLT_Qlon(:,1:30),2);
    NA_MLT_Qsen_clim=nanmean(NA_MLT_Qsen(:,1:30),2);

    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(NA_MLT_Qnet,2)
        NA_MLT_Qnet_ano(:,count)=NA_MLT_Qnet(:,count)-NA_MLT_Qnet_clim(:,1);
        NA_MLT_Qswr_ano(:,count)=NA_MLT_Qswr(:,count)-NA_MLT_Qswr_clim(:,1);
        NA_MLT_Qlat_ano(:,count)=NA_MLT_Qlat(:,count)-NA_MLT_Qlat_clim(:,1);
        NA_MLT_Qlon_ano(:,count)=NA_MLT_Qlon(:,count)-NA_MLT_Qlon_clim(:,1);
        NA_MLT_Qsen_ano(:,count)=NA_MLT_Qsen(:,count)-NA_MLT_Qsen_clim(:,1);
        count=count+1;
    end
    %clear *_clim
% #########################################################################
% #########################################################################




% #########################################################################
% SSTa from ERA5 reanalysis
% #########################################################################
%   % ACCESS Om2 0.25
    cd('/Users/z5195509/Documents/Data/ACCESS-OM2-025/access-om2-025_daily_NA_SST_fields/')
    % load('basin_mask_NA_ACCESS-OM2-025_era5_iaf_no_medi.mat','basin_mask_NA','lon_025','lat_025')
    load('basin_mask_NA_ACCESS-OM2-025_era5_iaf.mat','basin_mask_NA','lon_025','lat_025')
    
    cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
    load(['SST_NA_100W20E_Eq70N_ERA5_Daily_Averaged_Reanalysis_on_Single_Levels_2022.mat'],'sst_daily_NA','lon','lat')
       [lo,la]=meshgrid((lon)', (lat)');
       basin_mask_NA=griddata(lon_025,lat_025,basin_mask_NA',lo',la','nearest');
       clear lo la lon_025 lat_025
       sst_daily_NA=nanmean(sst_daily_NA,3);
       basin_mask_NA(isnan(sst_daily_NA))=NaN; clear sst_daily_NA
    
       [Sxy,~,~]=function_Cgrid_Area_Distance((lon)',(lat)');
       Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
% #########################################################################


% #########################################################################
    % Daily SST during 1981-2010
    cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
    sst_daily_19812022=nan(366,length(1981:2022)');
    count=1;
    for year=1981:2022
          disp(['Daily SST Year# ',num2str(year)])
          load(['SST_NA_100W20E_Eq70N_ERA5_Daily_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'.mat'],'sst_daily_NA','lon','lat')

          for days=1:size(sst_daily_NA,3)
              sst_daily_NA0=squeeze(sst_daily_NA(:,:,days));
              sst_daily_NA0(isnan(basin_mask_NA))=NaN;
              sst_daily_NA(:,:,days)=sst_daily_NA0;
              clear sst_daily_NA0
          end
          clear days
          sst_daily_NA=sst_daily_NA.*Sxy;
          
          % NA 0-60N
          sst_daily_19812022(1:size(sst_daily_NA,3),count)=squeeze(nansum(nansum(sst_daily_NA(:,1:241,:),2),1))./squeeze(nansum(nansum(Sxy(:,1:241),2),1));
          clear sst_daily_NA
          count=count+1;
    end
    clear year count

    
    % 2023 Daily SST Anomaly relative to the clim-mean
          disp(['Daily SST Year# 2023'])
          load('SST_NA_100W20E_Eq70N_ERA5_Daily_Averaged_Reanalysis_on_Single_Levels_2023.mat','sst_daily')
          sst_daily_NA=sst_daily; clear sst_daily

          for days=1:size(sst_daily_NA,3)
              sst_daily_NA0=squeeze(sst_daily_NA(:,:,days));
              sst_daily_NA0(isnan(basin_mask_NA))=NaN;
              sst_daily_NA(:,:,days)=sst_daily_NA0;
              clear sst_daily_NA0
          end
          clear days
          sst_daily_NA=sst_daily_NA.*Sxy;
          
          % NA 0-60N
          sst_daily_2023(1:365,1)=squeeze(nansum(nansum(sst_daily_NA(:,1:241,1:365),2),1))./squeeze(nansum(nansum(Sxy(:,1:241),2),1));
          sst_daily_2023(366,1)=NaN;
          clear sst_daily_NA

          
    % Anomaly is relative to the clim-daily mean
    sst_daily_19812010_clim=squeeze(nanmean(sst_daily_19812022(:,1:30),2));
    count=1;
    for year=1:size(sst_daily_19812022,2)
        sst_daily_19812022_ano(:,count)=sst_daily_19812022(:,count)-sst_daily_19812010_clim;
        count=count+1;
    end
    
    sst_daily_2023_ano(:,1)=sst_daily_2023(:,1)-sst_daily_19812010_clim; 
% #########################################################################




% #########################################################################
% % MLTa from IAP data
% % #########################################################################
%   % ACCESS Om2 0.25
     % cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
     load(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','Sxy_NA','lon','lat')
       [lo,la]=meshgrid((260:380)', (0.5:59.5)');
       basin_mask_NA=griddata(lon,lat,basin_mask_NA',lo',la','nearest');
       clear lo la lon_025 lat_025
       [Sxy,~,~]=function_Cgrid_Area_Distance((260:380)',(0.5:59.5)');
       Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
       Sxy_NA=Sxy; clear Sxy
% #########################################################################

% #########################################################################
    % Monthly MLT during 1981-2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023')
    mlt_month_19812023=nan(12,length(1981:2023)');
    count=1;
    for year=1981:2023
          disp(['Monthly MLT Year# ',num2str(year)])
          load(['mlt_monthly_Rh0_0125_GSW_IAP_C17_',num2str(year),'.mat'],'mlt_monthly','lon_IAP','lat_IAP')
          mlt_monthly_NA=mlt_monthly(240:360,91:150,:);
          clear mlt_monthly

          for month=1:size(mlt_monthly_NA,3)
              mlt_monthly_NA0=squeeze(mlt_monthly_NA(:,:,month));
              mlt_monthly_NA0(isnan(basin_mask_NA))=NaN;
              mlt_monthly_NA(:,:,month)=mlt_monthly_NA0;
              clear mlt_monthly_NA0
          end
          clear days
          mlt_monthly_NA=mlt_monthly_NA.*Sxy_NA;
          
          % NA 0-60N
          mlt_month_19812023(1:size(mlt_monthly_NA,3),count)=squeeze(nansum(nansum(mlt_monthly_NA(:,:,:),2),1))./squeeze(nansum(nansum(Sxy_NA(:,:),2),1));
          clear mlt_monthly_NA
          count=count+1;
    end
    clear year count

    % Anomaly is relative to the clim-daily mean
    mlt_month_19812010_clim=squeeze(nanmean(mlt_month_19812023(:,1:30),2));
    count=1;
    for year=1:size(mlt_month_19812023,2)
        mlt_month_19812023_ano(:,count)=mlt_month_19812023(:,count)-mlt_month_19812010_clim;
        count=count+1;
    end
% #########################################################################




% ########################################################################
%  Plotting 5-1: Daily SSTA from ERA5 and MLTA by Heat Flux
clc;
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.220; ixe = 0.220;  ixd = 0.10; ixw = (1-ixs-ixe-0*ixd)/1;
iys = 0.050; iye = 0.100;  iyd = 0.10; iyw = (1-iys-iye-1*iyd)/2;

%          [left            bottom      width height]
pos{102}  = [ixs          iys+0*iyw+0.1*iyd   ixw 1.35*iyw]; 


% dc1                 = [0    , 0.447, 0.741]; % Blue-ish.
% dc6                 = [0.301, 0.745, 0.933]; % Light Blue-ish.
% dc2                 = [0.850, 0.325, 0.098]; % Orange/Red-ish.
% dc7                 = [0.635, 0.078, 0.184]; % Red/Brown-ish.
% dc5                 = [0.466, 0.674, 0.188]; % Green-ish.
% dc3                 = [0.929, 0.694, 0.125]; % Yellow/Gold-ish.
% dc4                 = [0.494, 0.184, 0.556]; % Purple-ish.


clear color color0 
color=cbrewer('seq', 'Blues', 60,'pchip');
color(:,:)=color(60:-1:1,:);

    
subplot('position',pos{102})
    line2023=plot(1:366,sst_daily_2023_ano(:,1));
    set(line2023,'color',[0.3,0.3,0.3],'LineWidth',3,'linestyle','-'); 
    
    hold on
    line2023MLT3=plot(16:30:270,NA_MLT_Qlat_ano(1:9,43));
    set(line2023MLT3,'color',[0    , 0.447, 0.641],'LineWidth',3,'linestyle','-'); 
    hold on
    line2023MLT4=plot(16:30:270,NA_MLT_Qlon_ano(1:9,43));
    set(line2023MLT4,'color',[0.101, 0.645, 0.833],'LineWidth',3,'linestyle','-'); 
    hold on
    line2023MLT5=plot(16:30:270,NA_MLT_Qsen_ano(1:9,43));
    set(line2023MLT5,'color',[0.400, 0.850, 0.933],'LineWidth',3,'linestyle','-'); 
    hold on
    line2023MLT2=plot(16:30:270,NA_MLT_Qswr_ano(1:9,43));
    set(line2023MLT2,'color',[0.929, 0.694, 0.125],'LineWidth',3,'linestyle','-'); 
    hold on
    line2023MLT1=plot(16:30:270,NA_MLT_Qnet_ano(1:9,43));
    set(line2023MLT1,'color',[0.959, 0.494, 0.225],'LineWidth',6,'linestyle','-'); 
    
    hold on
    line2023MLT0=plot(16:30:270,mlt_month_19812023_ano(1:9,43));
    set(line2023MLT0,'color',[0.850, 0.325, 0.098],'LineWidth',6,'linestyle','-'); 

    
    leg101=legend([line2023 line2023MLT0 line2023MLT1 line2023MLT2 line2023MLT5 line2023MLT4 line2023MLT3],...
               '2023 SSTA, ERA5','2023 MLTA','2023 MLTA, Qnet','2023 MLTA, shortwave','2023 MLTA, sensible','2023 MLTA, longwave','2023 MLTA, latent heat',...
               'Location','northwest','NumColumns',1);
    set(leg101,'fontsize',20)
    hold on
    title(leg101,'North Atlantic SST and MLT Anomaly','fontsize',20')
    legend('boxoff')

    
    set(gca,'Ylim',[-0.0 1.5],'ycolor','k') 
    set(gca,'YTick',-1.5:0.5:4)
    set(gca,'YTickLabel',{'-1.5','-1.0','-0.5','0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0'},'fontsize',24)
    set(gca,'Xlim',[1 366]) 
    set(gca,'XTick',1:30:360)
    set(gca,'XTickLabel',{'           Jan','           Feb','           Mar','           Apr','           May','           Jun',...
                          '           Jul','           Aug','           Sep','           Oct','           Nov','            Dec'},'fontsize',24)

    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
    ylabel(['Difference [ ^oC ]'],'fontsize',24,'color','k','FontWeight','normal')

    % title('a. North Atlantic SST and predicted MLT anomaly','fontsize',24,'color','k','FontWeight','bold')
    text(-40,1.5,'b.','fontsize',24,'color','k','FontWeight','bold')

        
cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')




% ########################################################################
% ########################################################################
% Plotting 5-2: Bar Charts for Monthly MLTA and ML Warming Decomposition by Monthly Heat Flux
clc;clear
time_ann=(1981:2022)';

% #########################################################################
% #########################################################################
% Mixed layer warming by air-sea flues in 2023 - accumuated monthly heat flux 
count_yr=1;
for year=1981:2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
    load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Accumulated_',num2str(year),'.mat'],...
          'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen',...
          'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon')
    dMLT_Qnet_mon0(:,count_yr)=dMLT_Qnet_mon;
    dMLT_Qswr_mon0(:,count_yr)=dMLT_Qswr_mon;
    dMLT_Qlat_mon0(:,count_yr)=dMLT_Qlat_mon;
    dMLT_Qlon_mon0(:,count_yr)=dMLT_Qlon_mon;
    dMLT_Qsen_mon0(:,count_yr)=dMLT_Qsen_mon;
    clear NA_SST_mon*  
    count_yr=count_yr+1;
end
% Climatology
    NA_MLT_Qnet_clim=nanmean(dMLT_Qnet_mon0(:,1:30),2);
    NA_MLT_Qswr_clim=nanmean(dMLT_Qswr_mon0(:,1:30),2);
    NA_MLT_Qlat_clim=nanmean(dMLT_Qlat_mon0(:,1:30),2);
    NA_MLT_Qlon_clim=nanmean(dMLT_Qlon_mon0(:,1:30),2);
    NA_MLT_Qsen_clim=nanmean(dMLT_Qsen_mon0(:,1:30),2);

    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(dMLT_Qnet_mon0,2)
        NA_MLT_Qnet_ano(:,count)=dMLT_Qnet_mon0(:,count)-NA_MLT_Qnet_clim(:,1);
        NA_MLT_Qswr_ano(:,count)=dMLT_Qswr_mon0(:,count)-NA_MLT_Qswr_clim(:,1);
        NA_MLT_Qlat_ano(:,count)=dMLT_Qlat_mon0(:,count)-NA_MLT_Qlat_clim(:,1);
        NA_MLT_Qlon_ano(:,count)=dMLT_Qlon_mon0(:,count)-NA_MLT_Qlon_clim(:,1);
        NA_MLT_Qsen_ano(:,count)=dMLT_Qsen_mon0(:,count)-NA_MLT_Qsen_clim(:,1);
        count=count+1;
    end
    %clear *_clim
% #########################################################################
% #########################################################################



% #########################################################################
% % MLTa from IAP data
% % #########################################################################
%   % ACCESS Om2 0.25
     % cd('/Users/z5195509/Documents/Data/ERA-5/Daily_Averaged_Reanalysis_Single_Levels/Daily_Averaged_SST/')
     load(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','Sxy_NA','lon','lat')
       [lo,la]=meshgrid((260:380)', (0.5:59.5)');
       basin_mask_NA=griddata(lon,lat,basin_mask_NA',lo',la','nearest');
       clear lo la lon_025 lat_025
       [Sxy,~,~]=function_Cgrid_Area_Distance((260:380)',(0.5:59.5)');
       Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
       Sxy_NA=Sxy; clear Sxy
% #########################################################################

% #########################################################################
    % Monthly MLT during 1981-2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023')
    mlt_month_19812023=nan(12,length(1981:2023)');
    count=1;
    for year=1981:2023
          disp(['Monthly MLT Year# ',num2str(year)])
          load(['mlt_monthly_Rh0_0125_GSW_IAP_C17_',num2str(year),'.mat'],'mlt_monthly','lon_IAP','lat_IAP')
          mlt_monthly_NA=mlt_monthly(240:360,91:150,1:12);
          clear mlt_monthly

          for month=1:size(mlt_monthly_NA,3)
              mlt_monthly_NA0=squeeze(mlt_monthly_NA(:,:,month));
              mlt_monthly_NA0(isnan(basin_mask_NA))=NaN;
              mlt_monthly_NA(:,:,month)=mlt_monthly_NA0;
              clear mlt_monthly_NA0
          end
          clear days
          mlt_monthly_NA=mlt_monthly_NA.*Sxy_NA;
          
          % NA 0-60N
          mlt_month_19812023(1:size(mlt_monthly_NA,3),count)=squeeze(nansum(nansum(mlt_monthly_NA(:,:,:),2),1))./squeeze(nansum(nansum(Sxy_NA(:,:),2),1));
          clear mlt_monthly_NA
          count=count+1;
    end
    clear year count

    % Anomaly is relative to the clim-daily mean
    mlt_month_19812010_clim=squeeze(nanmean(mlt_month_19812023(:,1:30),2));
    count=1;
    for year=1:size(mlt_month_19812023,2)
        mlt_month_19812023_ano(:,count)=mlt_month_19812023(:,count)-mlt_month_19812010_clim;
        count=count+1;
    end
    dmlt_month_19812023_ano=mlt_month_19812023_ano(2:12,:)-mlt_month_19812023_ano(1:11,:);
% #########################################################################



% #########################################################################
ixs = 0.220; ixe = 0.220;  ixd = 0.10; ixw = (1-ixs-ixe-0*ixd)/1;
iys = 0.150; iye = 0.100;  iyd = 0.10; iyw = (1-iys-iye-1*iyd)/2;

%          [left            bottom      width height]
pos{101}  = [ixs          iys+1.2*iyw+1*iyd   ixw 1*iyw]; 

% dc1                 = [0    , 0.447, 0.741]; % Blue-ish.
% dc6                 = [0.301, 0.745, 0.933]; % Light Blue-ish.
% dc2                 = [0.850, 0.325, 0.098]; % Orange/Red-ish.
% dc7                 = [0.635, 0.078, 0.184]; % Red/Brown-ish.
% dc5                 = [0.466, 0.674, 0.188]; % Green-ish.
% dc3                 = [0.929, 0.694, 0.125]; % Yellow/Gold-ish.
% dc4                 = [0.494, 0.184, 0.556]; % Purple-ish.


clear color color0 
color=cbrewer('seq', 'Blues', 60,'pchip');
color(:,:)=color(60:-1:1,:);

    
% bar_dmlt(:,1)=dmlt_month_19812023_ano(1:8,43);
% bar_dmlt(:,2)=(NA_MLT_Qnet_ano(1:8,43)+NA_MLT_Qnet_ano(2:9,43))./2;
% bar_dmlt(:,3)=(NA_MLT_Qswr_ano(1:8,43)+NA_MLT_Qswr_ano(2:9,43))./2;
% bar_dmlt(:,4)=(NA_MLT_Qlat_ano(1:8,43)+NA_MLT_Qlat_ano(2:9,43))./2;
% bar_dmlt(:,5)=(NA_MLT_Qlon_ano(1:8,43)+NA_MLT_Qlon_ano(2:9,43))./2;
% bar_dmlt(:,6)=(NA_MLT_Qsen_ano(1:8,43)+NA_MLT_Qsen_ano(2:9,43))./2;
% bar_dmlt(:,7)=bar_dmlt(:,1)-bar_dmlt(:,2);

bar_dmlt(2:8,1)=(dmlt_month_19812023_ano(1:7,43)+dmlt_month_19812023_ano(2:8,43))./2;
bar_dmlt(1:9,2)=NA_MLT_Qnet_ano(1:9,43);
bar_dmlt(1:9,3)=NA_MLT_Qswr_ano(1:9,43);
bar_dmlt(1:9,4)=NA_MLT_Qlat_ano(1:9,43);
bar_dmlt(1:9,5)=NA_MLT_Qlon_ano(1:9,43);
bar_dmlt(1:9,6)=NA_MLT_Qsen_ano(1:9,43);
bar_dmlt(1:9,7)=bar_dmlt(:,1)-bar_dmlt(:,2);

subplot('position',pos{101}) 
   h0=bar(5:8,bar_dmlt(5:8,:));
       set(h0,'BarWidth',0.94); 
       set(h0(1),'FaceColor',[0.850, 0.325, 0.098],'EdgeColor',[0.850, 0.325, 0.098])
       hold on
       set(h0(2),'FaceColor',[0.959, 0.494, 0.225],'EdgeColor',[0.959, 0.494, 0.225])
       hold on
       set(h0(3),'FaceColor',[0.929, 0.694, 0.125],'EdgeColor',[0.929, 0.694, 0.125])
       hold on
       set(h0(4),'FaceColor',[0    , 0.447, 0.641],'EdgeColor',[0    , 0.447, 0.641])
       hold on
       set(h0(5),'FaceColor',[0.101, 0.645, 0.833],'EdgeColor',[0.101, 0.645, 0.833])
       hold on
       set(h0(6),'FaceColor',[0.400, 0.850, 0.933],'EdgeColor',[0.400, 0.850, 0.933])
       hold on
       set(h0(7),'FaceColor',[0.594, 0.284, 0.556],'EdgeColor',[0.594, 0.284, 0.556])

    hold on
    legend([h0(1) h0(2) h0(3) h0(4) h0(5) h0(6) h0(7)],...
           'MLT tendency','Qnet','Shortwave','Latent','Longwave','Sensible','Residual',...
           'Location','northeast','Orientation','vertical','NumColumns',2)
    hold on
    set(legend,'fontsize',20)
    hold on
    % title(leg4,'Monthly SST Anomaly','fontsize',20')
    legend('boxoff')
    
    % Legend will show names for each color
    set(gca,'Ylim',[-0.45 0.8],'ycolor','k') 
    set(gca,'YTick',-0.4:0.4:1.2)
    set(gca,'YTickLabel',{'-0.4','0','0.4','0.8','1.2'},'fontsize',24)
    set(gca,'Xlim',[4.5 8.5]) 
    set(gca,'XTick',4.5:0.5:8.5)
    set(gca,'XTickLabel',{[],'May',[],'Jun',[],'Jul',[],'Aug',[]},'fontsize',24)
%     set(gca,'XTickLabel',{[],'May-Jun',[],'Jun-Jul',[],'Jul-Aug',[],'Aug-Sep',[]},'fontsize',20)
    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
    ylabel(['[ ^oC per month ]'],'fontsize',24,'color','k','FontWeight','normal')
    
    text(4.0,0.8,'a.','fontsize',24,'color','k','FontWeight','bold')
        
    
cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')




%% #########################################################################
%  #########################################################################
%  6. Monthly MLT heat budget in 1981-2023 - for consistency with bar charts
clc;clear
time_ann=(1981:2023)';

% #########################################################################
%   % ACCESS OM2 0.25
%      load(['basin_mask_NA_100W20E_Eq70N_ERA5_from_ACCESS-OM2-025_era5_iaf.mat'],'basin_mask_NA','lon','lat')
%        [lo,la]=meshgrid((260:380)', (0.5:59.5)');
%        basin_mask_NA=griddata(lon,lat,basin_mask_NA',lo',la','nearest');
%        clear lo la lon lat
%        [Sxy,~,~]=function_Cgrid_Area_Distance((260:380)',(0.5:59.5)');
%        Sxy(isnan(basin_mask_NA))=NaN; % figure;imagesc(Sxy)
%        Sxy_NA=Sxy; clear Sxy
% #########################################################################


% #########################################################################
count_yr=1;
for year=1981:2023
    disp([' Cal monthly ML heat budget Year#',num2str(year)])
    % #####################################################################
    % Monthly MLD from IAP
      cd('/Users/z5195509/Documents/Data/IAP_Ocean/Global_GSW_Monthly/MLD')
      % dRh0=0.125 ########################################################
      if year<2023
          load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld','lon_IAP','lat_IAP')
          lat_IAP=lat_IAP(:,1);
          mld0(:,:,1:12)=squeeze(mld(:,:,1,1:12));
          clear mld
      else
          % dRh0=0.125 ########################################################
          load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld','lon_IAP','lat_IAP')
          lat_IAP=lat_IAP(:,1);
          mld0(:,:,1:9)=squeeze(mld(:,:,1,1:9));
          clear mld
          load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAPv4_Temp_2023.mat'],'mld')
          mld0(:,:,10:12)=squeeze(mld(:,:,1,10:12));
          clear mld
      end
      
      lon_IAP0(1:340,1)=lon_IAP(21:360,1);
      lon_IAP0(341:360,1)=lon_IAP(1:20,1)+360; clear lon_IAP
      lon_IAP=lon_IAP0; clear lon_IAP0
      
      mld00(1:340,:,:)=mld0(21:360,:,:);
      mld00(341:360,:,:)=mld0(1:20,:,:); clear mld0
      mld_monthly=mld00; clear mld00
      % figure;imagesc(mld_monthly(:,:,1))
    % #####################################################################
    
        
    % Monthly heat flux terms
      cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg/')
      load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'_1deg.mat'],...
            'shortwave','latent','longwave','sensible') 
      Qnet=latent+longwave+sensible+shortwave;
      % figure;imagesc(Qnet(240:360,91:150,1))
        
      % ##################################################################
      % Remove the shortwave penatration
        % The Vertical Redistribution of SWR, at MLD=mld_2023_daily
        R=0.58; h1=0.35; h2=23; 
        z=mld_monthly;
        F_PS77=R.*exp(-z./h1)+(1-R).*exp(-z./h2);
        
        shortwave_mld = - shortwave.*F_PS77;
        clear R h1 h2 z F_PS77 
        % figure;imagesc(shortwave_mld(240:360,91:150,1))
        Qnet  = Qnet + shortwave_mld;
        shortwave = shortwave + shortwave_mld;
        clear shortwave_mld
      % ##################################################################   
      
       
      % ##################################################################
      % Mixed layer warming due to Qnet
        % Heat flux term in heat budget
        MLT_Qnet_mon(:,:,:,count_yr)  = Qnet./mld_monthly./3992./1027.*(24*60*60*30);      % K/month
        clear Qnet
        MLT_Qswr_mon(:,:,:,count_yr)  = shortwave./mld_monthly./3992./1027.*(24*60*60*30); % K/month
        clear shortwave
        MLT_Qlat_mon(:,:,:,count_yr)  = latent./mld_monthly./3992./1027.*(24*60*60*30);    % K/month
        clear latent
        MLT_Qlon_mon(:,:,:,count_yr)  = longwave./mld_monthly./3992./1027.*(24*60*60*30);  % K/month
        clear longwave
        MLT_Qsen_mon(:,:,:,count_yr)  = sensible./mld_monthly./3992./1027.*(24*60*60*30);  % K/month
        clear sensible
        clear mld_monthly
        % figure;imagesc(MLT_Qnet_mon(240:360,91:150,5,43)); caxis([-3 3])
      % ###################################################################

      count_yr=count_yr+1;
end
clear year count_yr

% Climatology
    NA_MLT_Qnet_clim=nanmean(MLT_Qnet_mon(:,:,:,1:30),4);
    NA_MLT_Qswr_clim=nanmean(MLT_Qswr_mon(:,:,:,1:30),4);
    NA_MLT_Qlat_clim=nanmean(MLT_Qlat_mon(:,:,:,1:30),4);
    NA_MLT_Qlon_clim=nanmean(MLT_Qlon_mon(:,:,:,1:30),4);
    NA_MLT_Qsen_clim=nanmean(MLT_Qsen_mon(:,:,:,1:30),4);

    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(MLT_Qnet_mon,4)
        NA_MLT_Qnet_ano(:,:,:,count)=MLT_Qnet_mon(:,:,:,count)-NA_MLT_Qnet_clim;
        NA_MLT_Qswr_ano(:,:,:,count)=MLT_Qswr_mon(:,:,:,count)-NA_MLT_Qswr_clim;
        NA_MLT_Qlat_ano(:,:,:,count)=MLT_Qlat_mon(:,:,:,count)-NA_MLT_Qlat_clim;
        NA_MLT_Qlon_ano(:,:,:,count)=MLT_Qlon_mon(:,:,:,count)-NA_MLT_Qlon_clim;
        NA_MLT_Qsen_ano(:,:,:,count)=MLT_Qsen_mon(:,:,:,count)-NA_MLT_Qsen_clim;
        count=count+1;
    end
    clear *_clim *_mon
    % figure;imagesc(NA_MLT_Qnet_ano(:,:,5,43)); caxis([-3 3])
% #########################################################################
% #########################################################################


% #########################################################################
% MLTa from IAP data
    % Monthly MLT during 1981-2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023')
    count=1;
    for year=1981:2023
          disp(['Monthly MLT Year# ',num2str(year)])
          load(['mlt_monthly_Rh0_0125_GSW_IAP_C17_',num2str(year),'.mat'],'mlt_monthly','lon_IAP','lat_IAP')
          mlt_monthly_NA(:,:,:,count)=mlt_monthly(:,:,1:12);
          clear mlt_monthly
          count=count+1;
    end
    clear year count

    % Anomaly is relative to the clim-daily mean
    mlt_monthly_NA_clim=squeeze(nanmean(mlt_monthly_NA(:,:,:,1:30),4));
    count=1;
    for year=1:size(mlt_monthly_NA,4)
        mlt_monthly_NA_ano(:,:,:,count)=mlt_monthly_NA(:,:,:,count)-mlt_monthly_NA_clim;
        count=count+1;
    end
    clear mlt_monthly_NA_clim mlt_monthly_NA year count
    dmlt_month_19812023_ano=mlt_monthly_NA_ano(:,:,2:12,:)-mlt_monthly_NA_ano(:,:,1:11,:);
    clear mlt_monthly_NA_ano
    % Centered at mid-month
    bar_dmlt(:,:,2:8,:)=(dmlt_month_19812023_ano(:,:,1:7,:)+dmlt_month_19812023_ano(:,:,2:8,:))./2;
% #########################################################################
% #########################################################################



%% ########################################################################
% Plotting 6-1: Heat Budget Anomaly in May, June, July 2023
% #########################################################################
% #########################################################################  
clc;
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.110; ixe = 0.110;  ixd = 0.04; ixw = (1-ixs-ixe-2*ixd)/3;
iys = 0.000; iye = 0.050;  iyd = 0.06; iyw = (1-iys-iye-2*iyd)/3;

%          [left            bottom      width height]
pos{11}  = [1*ixs              iys+2*iyw+2*iyd   ixw 1.0*iyw]; 
pos{12}  = [1*ixs+1*ixw+1*ixd  iys+2*iyw+2*iyd   ixw 1.0*iyw]; 
pos{13}  = [ixs+2*ixw+2*ixd    iys+2*iyw+2*iyd   ixw 1.0*iyw]; 

pos{21}  = [1*ixs              iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{22}  = [1*ixs+1*ixw+1*ixd  iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{23}  = [ixs+2*ixw+2*ixd    iys+1*iyw+1*iyd   ixw 1.0*iyw]; 

pos{31}  = [1*ixs              iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
pos{32}  = [1*ixs+1*ixw+1*ixd  iys+0*iyw+0*iyd   ixw 1.0*iyw]; 
pos{33}  = [ixs+2*ixw+2*ixd    iys+0*iyw+0*iyd   ixw 1.0*iyw]; 


% dc1                 = [0    , 0.447, 0.741]; % Blue-ish.
% dc6                 = [0.301, 0.745, 0.933]; % Light Blue-ish.
% dc2                 = [0.850, 0.325, 0.098]; % Orange/Red-ish.
% dc7                 = [0.635, 0.078, 0.184]; % Red/Brown-ish.
% dc5                 = [0.466, 0.674, 0.188]; % Green-ish.
% dc3                 = [0.929, 0.694, 0.125]; % Yellow/Gold-ish.
% dc4                 = [0.494, 0.184, 0.556]; % Purple-ish.


clear color color0 
color=cbrewer('div', 'RdBu', 14,'pchip');
color0(1:12,:)=color(13:-1:2,:);
% color0(8,:)=(color0(7,:)+color0(8,:))./2;   
color0(6,:)=(color0(6,:)+color0(5,:))./2;    


% #########################################################################    
% 1. MLT tendency term in Jun-July 2023
subplot('position',pos{11})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,squeeze(bar_dmlt(:,:,6,43))');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',21,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.45 0.45 0.45]); 
        hold on;
        colormap(gca,color0)
        caxis([-1.8 1.8]);
        title('a. MLT tendency (June 2023)','fontsize',22,'FontWeight','bold')
        
         
subplot('position',pos{21})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,squeeze(bar_dmlt(:,:,7,43))');
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',21,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.45 0.45 0.45]); 
        hold on;
        colormap(gca,color0)
        caxis([-1.8 1.8]);
        title('d. MLT tendency (July 2023)','fontsize',22,'FontWeight','bold')
        
        
        
% #########################################################################    
% 2. Heat flux term in Jun-July 2023    
subplot('position',pos{12})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,squeeze(NA_MLT_Qnet_ano(:,:,6,43))');
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',21,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.45 0.45 0.45]); 
        hold on;
        colormap(gca,color0)
        caxis([-1.8 1.8]);
        title('b. Heat flux term (June 2023)','fontsize',22,'FontWeight','bold')
        

        
subplot('position',pos{22})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,squeeze(NA_MLT_Qnet_ano(:,:,7,43))');
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',21,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.45 0.45 0.45]); 
        hold on;
        colormap(gca,color0)
        caxis([-1.8 1.8]);
        title('e. Heat flux term (July 2023)','fontsize',22,'FontWeight','bold')


        

% #########################################################################    
% 3. Residual term in Jun-Aug 2023
subplot('position',pos{13})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,squeeze(bar_dmlt(:,:,6,43))'-squeeze(NA_MLT_Qnet_ano(:,:,6,43))');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',21,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.45 0.45 0.45]); 
        hold on;
        colormap(gca,color0)
        caxis([-1.8 1.8]);
        title('c. Residual','fontsize',22,'FontWeight','bold')
        
        
        
subplot('position',pos{23})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,squeeze(bar_dmlt(:,:,7,43))'-squeeze(NA_MLT_Qnet_ano(:,:,7,43))');
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',21,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.45 0.45 0.45]); 
        hold on;
        colormap(gca,color0)
        caxis([-1.8 1.8]);
        title('f. Residual','fontsize',22,'FontWeight','bold')
        
        
        hBar2 = colorbar('EastOutside','vertical');
        get(hBar2, 'Position');
        set(hBar2, 'Position', [ixs+3*ixw+2*ixd+0.016 iys+1.5*iyw+1*iyd 0.012 1*iyw+1*iyd]);
        set(hBar2, 'ytick',-1.8:0.3:1.8,'yticklabel',{'<-1.8',[],'-1.2',[],'-0.6',[],'0',[],'0.6',[],'1.2',[],'>1.8'},'fontsize',22,'FontName','Arial','LineWidth',1.2,'TickLength',0.06);
        ylabel(hBar2, '[ ^oC/month ]','rotation',90);

        
cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')





%% #########################################################################
% Plotting 6-2: Surface Heat Flux Anomaly in June July 2023
clc;
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.335; ixe = 0.335;  ixd = 0.045; ixw = (1-ixs-ixe-1*ixd)/2;
iys = 0.035; iye = 0.035;  iyd = 0.036; iyw = (1-iys-iye-3*iyd)/4;
ixw./iyw

%          [left            bottom      width height]
pos{11}  = [ixs+0*ixw+0*ixd    iys+3*iyw+3*iyd   ixw 1.0*iyw]; 
pos{12}  = [ixs+0*ixw+0*ixd    iys+2*iyw+2*iyd   ixw 1.0*iyw]; 
pos{13}  = [ixs+0*ixw+0*ixd    iys+1*iyw+1*iyd   ixw 1.0*iyw];
pos{14}  = [ixs+0*ixw+0*ixd    iys+0*iyw+0*iyd   ixw 1.0*iyw];

pos{21}  = [ixs+1*ixw+1*ixd    iys+3*iyw+3*iyd   ixw 1.0*iyw]; 
pos{22}  = [ixs+1*ixw+1*ixd    iys+2*iyw+2*iyd   ixw 1.0*iyw]; 
pos{23}  = [ixs+1*ixw+1*ixd    iys+1*iyw+1*iyd   ixw 1.0*iyw]; 
pos{24}  = [ixs+1*ixw+1*ixd    iys+0*iyw+0*iyd   ixw 1.0*iyw]; 


% dc1                 = [0    , 0.447, 0.741]; % Blue-ish.
% dc6                 = [0.301, 0.745, 0.933]; % Light Blue-ish.
% dc2                 = [0.850, 0.325, 0.098]; % Orange/Red-ish.
% dc7                 = [0.635, 0.078, 0.184]; % Red/Brown-ish.
% dc5                 = [0.466, 0.674, 0.188]; % Green-ish.
% dc3                 = [0.929, 0.694, 0.125]; % Yellow/Gold-ish.
% dc4                 = [0.494, 0.184, 0.556]; % Purple-ish.


clear color color0 
color=cbrewer('div', 'RdBu', 14,'pchip');
color0(1:12,:)=color(13:-1:2,:);
% color0(8,:)=(color0(7,:)+color0(8,:))./2;   
color0(6,:)=(color0(6,:)+color0(5,:))./2;    


% #########################################################################       
subplot('position',pos{11})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,NA_MLT_Qswr_ano(:,:,6,43)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('a. Net shortwave (June 2023)','fontsize',19,'FontWeight','bold')
        
        
         
subplot('position',pos{12})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,NA_MLT_Qlon_ano(:,:,6,43)');
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('c. Longwave radiation','fontsize',19,'FontWeight','bold')

        

subplot('position',pos{13})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,NA_MLT_Qsen_ano(:,:,6,43)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('e. Sensible heat','fontsize',19,'FontWeight','bold')
        
        
        
subplot('position',pos{14})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,NA_MLT_Qlat_ano(:,:,6,43)');
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('g. Latent heat','fontsize',19,'FontWeight','bold')
        
        


% #########################################################################       
subplot('position',pos{21})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,NA_MLT_Qswr_ano(:,:,7,43)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('b. Net shortwave (July 2023)','fontsize',19,'FontWeight','bold')
        
        
         
subplot('position',pos{22})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,NA_MLT_Qlon_ano(:,:,7,43)');
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('d. Longwave Radiation','fontsize',19,'FontWeight','bold')

        

subplot('position',pos{23})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,NA_MLT_Qsen_ano(:,:,7,43)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('f. Sensible heat','fontsize',19,'FontWeight','bold')
        
        
        
subplot('position',pos{24})
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,NA_MLT_Qlat_ano(:,:,7,43)');
        shading flat
        hold on
        
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',19,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('h. Latent heat','fontsize',19,'FontWeight','bold')
        
        
        hBar1 = colorbar('EastOutside','vertical');
        get(hBar1, 'Position');
        set(hBar1, 'Position', [ixs+2*ixw+1*ixd+0.016 iys+1.2*iyw+1*iyd 0.010 1.6*iyw+1*iyd]);
        set(hBar1, 'ytick',-2.4:0.4:2.4,'yticklabel',{'<-2.4',[],'-1.6',[],'-0.8',[],'0',[],'0.8',[],'1.6',[],'>2.4'},'fontsize',19,'FontName','Arial','LineWidth',1.2,'TickLength',0.055);
        ylabel(hBar1, '[ ^oC/month ]','rotation',90);


cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')



