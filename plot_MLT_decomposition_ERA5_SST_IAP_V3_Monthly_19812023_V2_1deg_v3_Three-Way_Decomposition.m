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



%% #########################################################################
%  #########################################################################
% 2.1 MLT Budget in 1981-2023
clc;clear
time_ann=(1981:2010)';

% #########################################################################
  % ACCESS Om2 0.25
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
          second_2_month0(:,:,month)=second_2_month(month,1).*ones(360,180);
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
      % Mixed layer warming due to Qnet
        disp(['   Projected ML Warming in Year#',num2str(year)])
        % Heat flux term in heat budget
        MLT_Qnet_mon(:,:,:)  = Qnet./mld_monthly./3992./1027.*second_2_month;  % K/month
        clear Qnet
        MLT_Qswr_mon(:,:,:)  = shortwave./mld_monthly./3992./1027.*second_2_month; % K/month
        clear shortwave
        MLT_Qlat_mon(:,:,:)  = latent./mld_monthly./3992./1027.*second_2_month;    % K/month
        clear latent
        MLT_Qlon_mon(:,:,:)  = longwave./mld_monthly./3992./1027.*second_2_month;  % K/month
        clear longwave
        MLT_Qsen_mon(:,:,:)  = sensible./mld_monthly./3992./1027.*second_2_month;  % K/month
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
        clear sst_daily %MLT_*_mon
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
         save(['MLT_budget_NA_monthly_ERA5_IAP_V2_1deg_V3_3way_decomposition_',num2str(year),'.mat'],...
              'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen',...
              'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon',...
              'MLT_Qnet_mon','MLT_Qswr_mon','MLT_Qlat_mon','MLT_Qlon_mon','MLT_Qsen_mon',...
              'basin_mask_NA','Sxy_NA','lon_IAP','lat_IAP')
          clear NA_SST_mon* dMLT_*_mon
         % ################################################################
      count_yr=count_yr+1;
end
clear *_dail_clim
% #########################################################################
% #########################################################################  



%% 2.2 MLT Budget in 1981-2023 - using MLD and SWR climatology
clc;clear
time_ann=(1981:2010)';

% #########################################################################
  % ACCESS Om2 0.25
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
          second_2_month0(:,:,month)=second_2_month(month,1).*ones(360,180);
      end
      second_2_month=second_2_month0; clear second_2_month0
% #########################################################################
      

% #########################################################################
  % Monthly MLD Climatology from IAP
    count=1;
    for year=1981:2010
        disp(['   MLD from IAPv3, Year#',num2str(year)])
          cd('/Users/z5195509/Documents/Data/IAP_Ocean/Global_GSW_Monthly/MLD')
          % dRh0=0.125 ########################################################
          load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld','lon_IAP','lat_IAP')
          lat_IAP=lat_IAP(:,1);
          if year==2023
             mld(:,:,:,10:12)=NaN;
          end
          mld0(:,:,1:12)=squeeze(mld(:,:,1,1:12));
          clear mld

          lon_IAP0(1:340,1)=lon_IAP(21:360,1);
          lon_IAP0(341:360,1)=lon_IAP(1:20,1)+360; clear lon_IAP
          lon_IAP=lon_IAP0; clear lon_IAP0

          mld00(1:340,:,:)=mld0(21:360,:,:);
          mld00(341:360,:,:)=mld0(1:20,:,:); clear mld0
          mld_monthly(:,:,:,count)=mld00; clear mld00
          count=count+1;
    end
    clear year count
    mld_monthly=nanmean(mld_monthly,4);
      % figure;imagesc(mld_monthly(:,:,1))
      % figure;plot(1:12,squeeze(mld_monthly(240,35,1:12)))
% #########################################################################


% #########################################################################
  % Monthly SWR climatology
      cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg/')
      count=1;
      for year_qnet=1981:2010
          disp(['    heat fluxes Year#',num2str(year_qnet)])
          load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year_qnet),'_1deg.mat'],...
                'shortwave') 
          shortwave0(:,:,:,count)=shortwave;
          clear shortwave
          count=count+1;
      end
      clear year_qnet count
      SWR_clim=nanmean(shortwave0,4); clear shortwave0
      % figure;imagesc(Qnet(240:360,91:150,1))    
% #########################################################################
      

% #########################################################################
count_yr=1;
for year=1981:2023
    disp([' Predicting NA daily SST Year#',num2str(year)])
    % #####################################################################
    % Monthly MLD from IAP
      cd('/Users/z5195509/Documents/Data/IAP_Ocean/Global_GSW_Monthly/MLD')
      % dRh0=0.125 ####################################################
      load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld')
      if year==2023
         mld(:,:,:,10:12)=NaN;
      end
      mld0(:,:,1:12)=squeeze(mld(:,:,1,1:12));
      clear mld
      mld00(1:340,:,:)=mld0(21:360,:,:);
      mld00(341:360,:,:)=mld0(1:20,:,:); clear mld0
      mld0(:,:,:)=mld00; clear mld00
          
          
    % ##################################################################### 
    % Monthly heat flux terms
      cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg/')
      load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'_1deg.mat'],...
            'shortwave','latent','longwave','sensible') 
      Qnet=latent+longwave+sensible+SWR_clim;
      % figure;imagesc(Qnet(240:360,91:150,1))
 
      
      % ##################################################################
      % Remove the shortwave penatration
        % The Vertical Redistribution of SWR, at MLD=mld_2023_daily
        R=0.58; h1=0.35; h2=23; 
        z=mld0;
        F_PS77=R.*exp(-z./h1)+(1-R).*exp(-z./h2);
        
        shortwave_mld = - shortwave.*F_PS77;
        clear R h1 h2 z F_PS77 
        % figure;imagesc(shortwave_mld(240:360,91:150,1))
        Qnet  = Qnet + shortwave_mld;
        shortwave = SWR_clim + shortwave_mld;
        clear shortwave_mld
      % ##################################################################   
       

      % ##################################################################
      % Mixed layer warming due to Qnet
        disp(['   Projected ML Warming in Year#',num2str(year)])
        % Heat flux term in heat budget
        MLT_Qnet_mon(:,:,:)  = Qnet./mld_monthly./3992./1027.*second_2_month;      % K/month
        clear Qnet
        MLT_Qswr_mon(:,:,:)  = shortwave./mld_monthly./3992./1027.*second_2_month; % K/month
        clear shortwave
        MLT_Qlat_mon(:,:,:)  = latent./mld_monthly./3992./1027.*second_2_month;    % K/month
        clear latent
        MLT_Qlon_mon(:,:,:)  = longwave./mld_monthly./3992./1027.*second_2_month;  % K/month
        clear longwave
        MLT_Qsen_mon(:,:,:)  = sensible./mld_monthly./3992./1027.*second_2_month;  % K/month
        clear sensible
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
        clear sst_daily %MLT_*_mon
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
         save(['MLT_budget_NA_monthly_ERA5_IAP_V2_1deg_V3_3way_decomposition_',num2str(year),'_2_Clim_MLD_SWR.mat'],...
              'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen',...
              'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon',...
              'MLT_Qnet_mon','MLT_Qswr_mon','MLT_Qlat_mon','MLT_Qlon_mon','MLT_Qsen_mon',...
              'basin_mask_NA','Sxy_NA','lon_IAP','lat_IAP')
          clear NA_SST_mon* dMLT_*_mon
         % ################################################################
      count_yr=count_yr+1;
end
clear *_dail_clim
% #########################################################################
% #########################################################################  



%% 2.3 MLT Budget in 1981-2023 - using MLD and Qsw_H climatology
clc;clear
time_ann=(1981:2010)';

% #########################################################################
  % ACCESS Om2 0.25
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
          second_2_month0(:,:,month)=second_2_month(month,1).*ones(360,180);
      end
      second_2_month=second_2_month0; clear second_2_month0 month
% #########################################################################
      

% #########################################################################
  % Monthly MLD Climatology from IAP
    count=1;
    for year=1981:2010
        disp(['   MLD from IAPv3, Year#',num2str(year)])
          cd('/Users/z5195509/Documents/Data/IAP_Ocean/Global_GSW_Monthly/MLD')
          % dRh0=0.125 ########################################################
          load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld','lon_IAP','lat_IAP')
          lat_IAP=lat_IAP(:,1);
          if year==2023
             mld(:,:,:,10:12)=NaN;
          end
          mld0(:,:,1:12)=squeeze(mld(:,:,1,1:12));
          clear mld

          lon_IAP0(1:340,1)=lon_IAP(21:360,1);
          lon_IAP0(341:360,1)=lon_IAP(1:20,1)+360; clear lon_IAP
          lon_IAP=lon_IAP0; clear lon_IAP0

          mld00(1:340,:,:)=mld0(21:360,:,:);
          mld00(341:360,:,:)=mld0(1:20,:,:); clear mld0
          mld_monthly(:,:,:,count)=mld00; clear mld00
          count=count+1;
    end
    clear year count
    mld_monthly=nanmean(mld_monthly,4);
      % figure;imagesc(mld_monthly(:,:,1))
% #########################################################################


% #########################################################################
  % Monthly Qswr_H climatology
    count=1;
    for year=1981:2010
        disp(['   Qsw and MLD, Year#',num2str(year)])
          cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg/')
          load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'_1deg.mat'],...
                'shortwave') 

          cd('/Users/z5195509/Documents/Data/IAP_Ocean/Global_GSW_Monthly/MLD')
          % dRh0=0.125 ####################################################
          load(['mld_Rh0_0125_GSW_SA_CT_Interp_IAP_C17_',num2str(year),'.mat'],'mld')
          mld0(:,:,1:12)=squeeze(mld(:,:,1,1:12));
          clear mld
          mld00(1:340,:,:)=mld0(21:360,:,:);
          mld00(341:360,:,:)=mld0(1:20,:,:); clear mld0
          mld0(:,:,:)=mld00; clear mld00
          
          % ###############################################################
          % Remove the shortwave penatration
            % The Vertical Redistribution of SWR, at MLD=mld_2023_daily
            R=0.58; h1=0.35; h2=23; 
            z=mld0; clear mld0
            F_PS77=R.*exp(-z./h1)+(1-R).*exp(-z./h2);

            shortwave_mld = - shortwave.*F_PS77;
            clear R h1 h2 z F_PS77 shortwave
            % figure;imagesc(shortwave_mld(240:360,91:150,1))
            SWR_mld(:,:,:,count)=shortwave_mld;
            clear shortwave_mld
          % ###############################################################
          count=count+1;
    end
    clear year count 
    SWR_mld=nanmean(SWR_mld,4);
% #########################################################################
      

% #########################################################################
count_yr=1;
for year=1981:2023
    disp([' Predicting NA daily SST Year#',num2str(year)])
    % ##################################################################### 
    % Monthly heat flux terms
      cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg/')
      load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year),'_1deg.mat'],...
            'shortwave','latent','longwave','sensible') 
      Qnet=latent+longwave+sensible+shortwave;
      % figure;imagesc(Qnet(240:360,91:150,1))
        
      % ##################################################################
%       % Remove the shortwave penatration
%         % The Vertical Redistribution of SWR, at MLD=mld_2023_daily
%         R=0.58; h1=0.35; h2=23; 
%         z=mld_monthly;
%         F_PS77=R.*exp(-z./h1)+(1-R).*exp(-z./h2);
%         
%         shortwave_mld = - shortwave.*F_PS77;
%         clear R h1 h2 z F_PS77 
%         % figure;imagesc(shortwave_mld(240:360,91:150,1))
        Qnet  = Qnet + SWR_mld;
        shortwave = shortwave + SWR_mld;
      % ##################################################################   
       

      % ##################################################################
      % Mixed layer warming due to Qnet
        disp(['   Projected ML Warming in Year#',num2str(year)])
        % Heat flux term in heat budget
        MLT_Qnet_mon(:,:,:)  = Qnet./mld_monthly./3992./1027.*second_2_month;      % K/month
        clear Qnet
        MLT_Qswr_mon(:,:,:)  = shortwave./mld_monthly./3992./1027.*second_2_month; % K/month
        clear shortwave
        MLT_Qlat_mon(:,:,:)  = latent./mld_monthly./3992./1027.*second_2_month;    % K/month
        clear latent
        MLT_Qlon_mon(:,:,:)  = longwave./mld_monthly./3992./1027.*second_2_month;  % K/month
        clear longwave
        MLT_Qsen_mon(:,:,:)  = sensible./mld_monthly./3992./1027.*second_2_month;  % K/month
        clear sensible
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
        clear sst_daily %MLT_*_mon
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
         save(['MLT_budget_NA_monthly_ERA5_IAP_V2_1deg_V3_3way_decomposition_',num2str(year),'_3_Clim_MLD_Qsw_H.mat'],...
              'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen',...
              'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon',...
              'MLT_Qnet_mon','MLT_Qswr_mon','MLT_Qlat_mon','MLT_Qlon_mon','MLT_Qsen_mon',...
              'basin_mask_NA','Sxy_NA','lon_IAP','lat_IAP')
          clear NA_SST_mon* dMLT_*_mon
         % ################################################################
      count_yr=count_yr+1;
end
clear *_dail_clim
% #########################################################################
% #########################################################################  



%% 2.4 MLT Budget in 1981-2023 - using SWR and Qsw_H climatology
clc;clear
time_ann=(1981:2010)';

% #########################################################################
  % ACCESS Om2 0.25
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
          second_2_month0(:,:,month)=second_2_month(month,1).*ones(360,180);
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
      count=1;
      for year_qnet=1981:2010
          disp(['    heat fluxes Year#',num2str(year_qnet)])
          load(['heatflux_W_m2_ERA5_Monthly_Averaged_Reanalysis_on_Single_Levels_',num2str(year_qnet),'_1deg.mat'],...
                'shortwave','latent','longwave','sensible') 
          shortwave0(:,:,:,count)=shortwave;
          latent0(:,:,:,count)=latent;
          longwave0(:,:,:,count)=longwave;
          sensible0(:,:,:,count)=sensible;
          Qnet0(:,:,:,count)=latent+longwave+sensible+shortwave;
          clear latent longwave sensible shortwave
          count=count+1;
      end
      clear year_qnet count
      Qnet=nanmean(Qnet0,4); clear Qnet0
      shortwave=nanmean(shortwave0,4); clear shortwave0
      latent=nanmean(latent0,4); clear latent0
      longwave=nanmean(longwave0,4); clear longwave0
      sensible=nanmean(sensible0,4); clear sensible0
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
        disp(['   Projected ML Warming in Year#',num2str(year)])
        % Heat flux term in heat budget
        MLT_Qnet_mon(:,:,:)  = Qnet./mld_monthly./3992./1027.*second_2_month;  % K/month
        clear Qnet
        MLT_Qswr_mon(:,:,:)  = shortwave./mld_monthly./3992./1027.*second_2_month; % K/month
        clear shortwave
        MLT_Qlat_mon(:,:,:)  = latent./mld_monthly./3992./1027.*second_2_month;    % K/month
        clear latent
        MLT_Qlon_mon(:,:,:)  = longwave./mld_monthly./3992./1027.*second_2_month;  % K/month
        clear longwave
        MLT_Qsen_mon(:,:,:)  = sensible./mld_monthly./3992./1027.*second_2_month;  % K/month
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
        clear sst_daily %MLT_*_mon
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
         save(['MLT_budget_NA_monthly_ERA5_IAP_V2_1deg_V3_3way_decomposition_',num2str(year),'_4_Clim_SWR_Qsw_H.mat'],...
              'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen',...
              'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon',...
              'MLT_Qnet_mon','MLT_Qswr_mon','MLT_Qlat_mon','MLT_Qlon_mon','MLT_Qsen_mon',...
              'basin_mask_NA','Sxy_NA','lon_IAP','lat_IAP')
          clear NA_SST_mon* dMLT_*_mon
         % ################################################################
      count_yr=count_yr+1;
end
clear *_dail_clim
% #########################################################################
% #########################################################################  




%% ########################################################################
%  Plotting 5-4: Bar Charts for Monthly MLTA and ML Warming Decomposition by Monthly Heat Flux
clc;clear
time_ann=(1981:2022)';


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
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.250; ixe = 0.250;  ixd = 0.10; ixw = (1-ixs-ixe-0*ixd)/1;
iys = 0.045; iye = 0.025;  iyd = 0.05; iyw = (1-iys-iye-2*iyd)/3;

%          [left            bottom      width height]
pos{101}  = [ixs          iys+2*iyw+2*iyd   ixw 1*iyw]; 
pos{102}  = [ixs          iys+1*iyw+1*iyd   ixw 1*iyw]; 
pos{103}  = [ixs          iys+0*iyw+0*iyd   ixw 1*iyw]; 

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
% Mixed layer warming by air-sea flues in 2023 - accumuated monthly heat flux 
count_yr=1;
for year=1981:2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
    load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Accumulated_',num2str(year),'.mat'],...
          'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen',...
          'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon')
%     % test with MLD clim values
%     load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Accumulated_',num2str(year),'_2_Clim_MLD.mat'],...
%           'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen',...
%           'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon')
%     % test with Qnet clim values
%     load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Accumulated_',num2str(year),'_3_Clim_Qnet.mat'],...
%           'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen',...
%           'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon')
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
           'MLT tendency','Surface flux term','Shortwave','Latent','Longwave','Sensible','Residual',...
           'Location','northeast','Orientation','vertical','NumColumns',2)
    hold on
    set(legend,'fontsize',16)
    hold on
    % title(leg4,'Monthly SST Anomaly','fontsize',20')
    legend('boxoff')
    
    % Legend will show names for each color
    set(gca,'Ylim',[-0.45 0.8],'ycolor','k') 
    set(gca,'YTick',-0.4:0.4:1.2)
    set(gca,'YTickLabel',{'-0.4','0','0.4','0.8','1.2'},'fontsize',22)
    set(gca,'Xlim',[4.5 8.5]) 
    set(gca,'XTick',4.5:0.5:8.5)
    set(gca,'XTickLabel',{[],'May',[],'Jun',[],'Jul',[],'Aug',[]},'fontsize',22)
%     set(gca,'XTickLabel',{[],'May-Jun',[],'Jun-Jul',[],'Jul-Aug',[],'Aug-Sep',[]},'fontsize',22)
    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
    ylabel(['[ ^oC per month ]'],'fontsize',22,'color','k','FontWeight','normal')
    
    text(4.6,0.7,'a. Temperature budget anomalies','fontsize',22,'color','k','FontWeight','bold')
    
    


% #########################################################################
% #########################################################################
%     % Monthly MLT during 1981-2023 - using MLD climatology
%     cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023')
%     mlt_month_19812023=nan(12,length(1981:2023)');
%     count=1;
%     for year=1981:2023
%           disp(['Monthly MLT Year# ',num2str(year)])
%           load(['mlt_monthly_Rh0_0125_GSW_IAP_C17_',num2str(year),'_MLD_Climatology.mat'],'mlt_monthly','lon_IAP','lat_IAP')
%           if year==2023
%               mlt_monthly(:,:,10:12)=NaN;
%           end
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
% #########################################################################


% #########################################################################
% Mixed layer warming by air-sea flues in 2023 - accumuated monthly heat flux 
clear dMLT* NA_MLT*
count_yr=1;
for year=1981:2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
    % test with MLD clim values - 1981-2010
    load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Accumulated_',num2str(year),'_2_Clim_MLD.mat'],...
          'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen',...
          'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon')
%     % test with MLD clim values - 1981-2023
%     load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Accumulated_',num2str(year),'_2_Clim_MLD_19812023.mat'],...
%           'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen',...
%           'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon')
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



% Centered at mid-month
bar_dmlt(2:8,1)=(dmlt_month_19812023_ano(1:7,43)+dmlt_month_19812023_ano(2:8,43))./2;
bar_dmlt(1:9,2)=NA_MLT_Qnet_ano(1:9,43);
bar_dmlt(1:9,3)=NA_MLT_Qswr_ano(1:9,43);
bar_dmlt(1:9,4)=NA_MLT_Qlat_ano(1:9,43);
bar_dmlt(1:9,5)=NA_MLT_Qlon_ano(1:9,43);
bar_dmlt(1:9,6)=NA_MLT_Qsen_ano(1:9,43);
bar_dmlt(1:9,7)=bar_dmlt(:,1)-bar_dmlt(:,2);

subplot('position',pos{102}) 
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
           'MLT tendency','Surface flux term','Shortwave','Latent','Longwave','Sensible','Residual',...
           'Location','northeast','Orientation','vertical','NumColumns',2)
    hold on
    set(legend,'fontsize',16)
    hold on
    % title(leg4,'Monthly SST Anomaly','fontsize',20')
    legend('boxoff')
    
    % Legend will show names for each color
    set(gca,'Ylim',[-0.45 0.8],'ycolor','k') 
    set(gca,'YTick',-0.4:0.4:1.2)
    set(gca,'YTickLabel',{'-0.4','0','0.4','0.8','1.2'},'fontsize',22)
    set(gca,'Xlim',[4.5 8.5]) 
    set(gca,'XTick',4.5:0.5:8.5)
    set(gca,'XTickLabel',{[],'May',[],'Jun',[],'Jul',[],'Aug',[]},'fontsize',22)
%     set(gca,'XTickLabel',{[],'May-Jun',[],'Jun-Jul',[],'Jul-Aug',[],'Aug-Sep',[]},'fontsize',22)
    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
    ylabel(['[ ^oC per month ]'],'fontsize',22,'color','k','FontWeight','normal')
        
    text(4.6,0.7,'b. MLD climatology 1981-2010','fontsize',22,'color','k','FontWeight','bold')
%     text(4.6,0.7,'b. MLD climatology 1981-2023','fontsize',22,'color','k','FontWeight','bold')  
    
    
% #########################################################################
% #########################################################################
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
% #########################################################################


% #########################################################################
% Mixed layer warming by air-sea flues in 2023 - accumuated monthly heat flux 
clear dMLT* NA_MLT*
count_yr=1;
for year=1981:2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
    % test with Qnet clim values - 1981-2010
    load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Accumulated_',num2str(year),'_3_Clim_Qnet.mat'],...
          'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen',...
          'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon')
%     % test with Qnet clim values - 1981-2023
%     load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Accumulated_',num2str(year),'_3_Clim_Qnet_19812023.mat'],...
%           'NA_SST_mon_Qnet','NA_SST_mon_Qswr','NA_SST_mon_Qlat','NA_SST_mon_Qlon','NA_SST_mon_Qsen',...
%           'dMLT_Qnet_mon','dMLT_Qswr_mon','dMLT_Qlat_mon','dMLT_Qlon_mon','dMLT_Qsen_mon')
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



% Centered at mid-month
bar_dmlt(2:8,1)=(dmlt_month_19812023_ano(1:7,43)+dmlt_month_19812023_ano(2:8,43))./2;
bar_dmlt(1:9,2)=NA_MLT_Qnet_ano(1:9,43);
bar_dmlt(1:9,3)=NA_MLT_Qswr_ano(1:9,43);
bar_dmlt(1:9,4)=NA_MLT_Qlat_ano(1:9,43);
bar_dmlt(1:9,5)=NA_MLT_Qlon_ano(1:9,43);
bar_dmlt(1:9,6)=NA_MLT_Qsen_ano(1:9,43);
bar_dmlt(1:9,7)=bar_dmlt(:,1)-bar_dmlt(:,2);

subplot('position',pos{103}) 
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
           'MLT tendency','Surface flux term','Shortwave','Latent','Longwave','Sensible','Residual',...
           'Location','northeast','Orientation','vertical','NumColumns',2)
    hold on
    set(legend,'fontsize',16)
    hold on
    % title(leg4,'Monthly SST Anomaly','fontsize',20')
    legend('boxoff')
    
    % Legend will show names for each color
    set(gca,'Ylim',[-0.45 0.8],'ycolor','k') 
    set(gca,'YTick',-0.4:0.4:1.2)
    set(gca,'YTickLabel',{'-0.4','0','0.4','0.8','1.2'},'fontsize',22)
    set(gca,'Xlim',[4.5 8.5]) 
    set(gca,'XTick',4.5:0.5:8.5)
    set(gca,'XTickLabel',{[],'May',[],'Jun',[],'Jul',[],'Aug',[]},'fontsize',22)
%     set(gca,'XTickLabel',{[],'May-Jun',[],'Jun-Jul',[],'Jul-Aug',[],'Aug-Sep',[]},'fontsize',22)
    grid on
    set(gca,'GridColor',[.8 .8 .8],'GridAlpha',0.3,'GridLineStyle','-')
    ylabel(['[ ^oC per month ]'],'fontsize',22,'color','k','FontWeight','normal')
      
    text(4.6,0.7,'c. Qnet climatology 1981-2010','fontsize',22,'color','k','FontWeight','bold')
%     text(4.6,0.7,'c. Qnet climatology 1981-2023','fontsize',22,'color','k','FontWeight','bold')  
    
        
cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')




%% ########################################################################
%  Plotting 5-5: Maps for ML Warming by Heat Flux, May-Aug 2023
clc;clear
time_ann=(1981:2022)';


% #########################################################################
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.040; ixe = 0.085;  ixd = 0.03; ixw = (1-ixs-ixe-3*ixd)/4;
iys = 0.060; iye = 0.060;  iyd = 0.05; iyw = (1-iys-iye-2*iyd)/3;

ixw./iyw

%          [left            bottom      width height]
pos{101}  = [ixs          iys+2*iyw+2*iyd   ixw 1*iyw]; 
pos{201}  = [ixs          iys+1*iyw+1*iyd   ixw 1*iyw]; 
pos{301}  = [ixs          iys+0*iyw+0*iyd   ixw 1*iyw]; 

pos{102}  = [ixs+1*ixw+1*ixd   iys+2*iyw+2*iyd   ixw 1*iyw]; 
pos{202}  = [ixs+1*ixw+1*ixd   iys+1*iyw+1*iyd   ixw 1*iyw]; 
pos{302}  = [ixs+1*ixw+1*ixd   iys+0*iyw+0*iyd   ixw 1*iyw]; 

pos{103}  = [ixs+2*ixw+2*ixd   iys+2*iyw+2*iyd   ixw 1*iyw]; 
pos{203}  = [ixs+2*ixw+2*ixd   iys+1*iyw+1*iyd   ixw 1*iyw]; 
pos{303}  = [ixs+2*ixw+2*ixd   iys+0*iyw+0*iyd   ixw 1*iyw]; 

pos{104}  = [ixs+3*ixw+3*ixd   iys+2*iyw+2*iyd   ixw 1*iyw]; 
pos{204}  = [ixs+3*ixw+3*ixd   iys+1*iyw+1*iyd   ixw 1*iyw]; 
pos{304}  = [ixs+3*ixw+3*ixd   iys+0*iyw+0*iyd   ixw 1*iyw]; 

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
% #########################################################################
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
    mlt_month_19812010_clim=squeeze(nanmean(mlt_monthly_NA(:,:,:,1:30),4));
    count=1;
    for year=1:size(mlt_monthly_NA,4)
        mlt_month_19812023_ano(:,:,:,count)=mlt_monthly_NA(:,:,:,count)-mlt_month_19812010_clim;
        count=count+1;
    end
    dmlt_month_19812023_ano=mlt_month_19812023_ano(:,:,2:12,:)-mlt_month_19812023_ano(:,:,1:11,:);
    
    % Centered at mid-month
    dmlt_2023(:,:,2:8,1)=(dmlt_month_19812023_ano(:,:,1:7,43)+dmlt_month_19812023_ano(:,:,2:8,43))./2;
    clear mlt* dmlt_month_19812023_ano
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
% Mixed layer warming by air-sea flues in 2023 - accumuated monthly heat flux 
count_yr=1;
for year=1981:2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
    load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Accumulated_',num2str(year),'.mat'],...
          'MLT_Qnet_mon','MLT_Qswr_mon','MLT_Qlat_mon','MLT_Qlon_mon','MLT_Qsen_mon')
    MLT_Qnet_mon0(:,:,:,count_yr)=MLT_Qnet_mon;
    MLT_Qswr_mon0(:,:,:,count_yr)=MLT_Qswr_mon;
    MLT_Qlat_mon0(:,:,:,count_yr)=MLT_Qlat_mon;
    MLT_Qlon_mon0(:,:,:,count_yr)=MLT_Qlon_mon;
    MLT_Qsen_mon0(:,:,:,count_yr)=MLT_Qsen_mon;
    clear MLT_*_mon
    count_yr=count_yr+1;
end
clear count_yr year
% Climatology
    NA_MLT_Qnet_clim=nanmean(MLT_Qnet_mon0(:,:,:,1:30),4);
    NA_MLT_Qswr_clim=nanmean(MLT_Qswr_mon0(:,:,:,1:30),4);
    NA_MLT_Qlat_clim=nanmean(MLT_Qlat_mon0(:,:,:,1:30),4);
    NA_MLT_Qlon_clim=nanmean(MLT_Qlon_mon0(:,:,:,1:30),4);
    NA_MLT_Qsen_clim=nanmean(MLT_Qsen_mon0(:,:,:,1:30),4);
    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(MLT_Qnet_mon0,4)
        NA_MLT_Qnet_ano(:,:,:,count)=MLT_Qnet_mon0(:,:,:,count)-NA_MLT_Qnet_clim(:,:,:,1);
        NA_MLT_Qswr_ano(:,:,:,count)=MLT_Qswr_mon0(:,:,:,count)-NA_MLT_Qswr_clim(:,:,:,1);
        NA_MLT_Qlat_ano(:,:,:,count)=MLT_Qlat_mon0(:,:,:,count)-NA_MLT_Qlat_clim(:,:,:,1);
        NA_MLT_Qlon_ano(:,:,:,count)=MLT_Qlon_mon0(:,:,:,count)-NA_MLT_Qlon_clim(:,:,:,1);
        NA_MLT_Qsen_ano(:,:,:,count)=MLT_Qsen_mon0(:,:,:,count)-NA_MLT_Qsen_clim(:,:,:,1);
        count=count+1;
    end
    %clear *_clim
    dmlt_2023_Qnet=NA_MLT_Qnet_ano(:,:,1:12,43);
% #########################################################################



% #########################################################################       
subplot('position',pos{101}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,5)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('a. Heat flux term (May 2023)','fontsize',20,'FontWeight','bold')
        
subplot('position',pos{102}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,6)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('b. Heat flux term (Jun 2023)','fontsize',20,'FontWeight','bold')
        
subplot('position',pos{103}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,7)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('c. Heat flux term (Jul 2023)','fontsize',20,'FontWeight','bold')
        
subplot('position',pos{104}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,8)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('d. Heat flux term (Aug 2023)','fontsize',20,'FontWeight','bold')

        clear dmlt_2023_Qnet NA_MLT* MLT*
% #########################################################################
% #########################################################################  



% #########################################################################
% #########################################################################
% Mixed layer warming by air-sea flues in 2023 - accumuated monthly heat flux 
clear dMLT* NA_MLT* 
count_yr=1;
for year=1981:2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
    % test with MLD clim values
    load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Accumulated_',num2str(year),'_2_Clim_MLD.mat'],...
          'MLT_Qnet_mon','MLT_Qswr_mon','MLT_Qlat_mon','MLT_Qlon_mon','MLT_Qsen_mon')
    MLT_Qnet_mon0(:,:,:,count_yr)=MLT_Qnet_mon;
    MLT_Qswr_mon0(:,:,:,count_yr)=MLT_Qswr_mon;
    MLT_Qlat_mon0(:,:,:,count_yr)=MLT_Qlat_mon;
    MLT_Qlon_mon0(:,:,:,count_yr)=MLT_Qlon_mon;
    MLT_Qsen_mon0(:,:,:,count_yr)=MLT_Qsen_mon;
    clear MLT_*_mon
    count_yr=count_yr+1;
end
clear count_yr year
% Climatology
    NA_MLT_Qnet_clim=nanmean(MLT_Qnet_mon0(:,:,:,1:30),4);
    NA_MLT_Qswr_clim=nanmean(MLT_Qswr_mon0(:,:,:,1:30),4);
    NA_MLT_Qlat_clim=nanmean(MLT_Qlat_mon0(:,:,:,1:30),4);
    NA_MLT_Qlon_clim=nanmean(MLT_Qlon_mon0(:,:,:,1:30),4);
    NA_MLT_Qsen_clim=nanmean(MLT_Qsen_mon0(:,:,:,1:30),4);
    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(MLT_Qnet_mon0,4)
        NA_MLT_Qnet_ano(:,:,:,count)=MLT_Qnet_mon0(:,:,:,count)-NA_MLT_Qnet_clim(:,:,:,1);
        NA_MLT_Qswr_ano(:,:,:,count)=MLT_Qswr_mon0(:,:,:,count)-NA_MLT_Qswr_clim(:,:,:,1);
        NA_MLT_Qlat_ano(:,:,:,count)=MLT_Qlat_mon0(:,:,:,count)-NA_MLT_Qlat_clim(:,:,:,1);
        NA_MLT_Qlon_ano(:,:,:,count)=MLT_Qlon_mon0(:,:,:,count)-NA_MLT_Qlon_clim(:,:,:,1);
        NA_MLT_Qsen_ano(:,:,:,count)=MLT_Qsen_mon0(:,:,:,count)-NA_MLT_Qsen_clim(:,:,:,1);
        count=count+1;
    end
    %clear *_clim
    dmlt_2023_Qnet=NA_MLT_Qnet_ano(:,:,1:12,43);
% #########################################################################



% #########################################################################       
subplot('position',pos{201}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,5)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('e. Heat flux term, MLD clim.','fontsize',20,'FontWeight','bold')
        
subplot('position',pos{202}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,6)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('f. Heat flux term, MLD clim.','fontsize',20,'FontWeight','bold')
        
subplot('position',pos{203}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,7)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('g. Heat flux term, MLD clim.','fontsize',20,'FontWeight','bold')
        
subplot('position',pos{204}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,8)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('h. Heat flux term, MLD clim.','fontsize',20,'FontWeight','bold')

        clear dmlt_2023_Qnet NA_MLT* MLT*
% #########################################################################
% #########################################################################


    
    
% #########################################################################
% #########################################################################
% Mixed layer warming by air-sea flues in 2023 - accumuated monthly heat flux 
clear dMLT* NA_MLT* 
count_yr=1;
for year=1981:2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
    % test with Qnet clim values
    load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Accumulated_',num2str(year),'_3_Clim_Qnet.mat'],...
          'MLT_Qnet_mon','MLT_Qswr_mon','MLT_Qlat_mon','MLT_Qlon_mon','MLT_Qsen_mon')
    MLT_Qnet_mon0(:,:,:,count_yr)=MLT_Qnet_mon;
    MLT_Qswr_mon0(:,:,:,count_yr)=MLT_Qswr_mon;
    MLT_Qlat_mon0(:,:,:,count_yr)=MLT_Qlat_mon;
    MLT_Qlon_mon0(:,:,:,count_yr)=MLT_Qlon_mon;
    MLT_Qsen_mon0(:,:,:,count_yr)=MLT_Qsen_mon;
    clear MLT_*_mon
    count_yr=count_yr+1;
end
clear count_yr year
% Climatology
    NA_MLT_Qnet_clim=nanmean(MLT_Qnet_mon0(:,:,:,1:30),4);
    NA_MLT_Qswr_clim=nanmean(MLT_Qswr_mon0(:,:,:,1:30),4);
    NA_MLT_Qlat_clim=nanmean(MLT_Qlat_mon0(:,:,:,1:30),4);
    NA_MLT_Qlon_clim=nanmean(MLT_Qlon_mon0(:,:,:,1:30),4);
    NA_MLT_Qsen_clim=nanmean(MLT_Qsen_mon0(:,:,:,1:30),4);
    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(MLT_Qnet_mon0,4)
        NA_MLT_Qnet_ano(:,:,:,count)=MLT_Qnet_mon0(:,:,:,count)-NA_MLT_Qnet_clim(:,:,:,1);
        NA_MLT_Qswr_ano(:,:,:,count)=MLT_Qswr_mon0(:,:,:,count)-NA_MLT_Qswr_clim(:,:,:,1);
        NA_MLT_Qlat_ano(:,:,:,count)=MLT_Qlat_mon0(:,:,:,count)-NA_MLT_Qlat_clim(:,:,:,1);
        NA_MLT_Qlon_ano(:,:,:,count)=MLT_Qlon_mon0(:,:,:,count)-NA_MLT_Qlon_clim(:,:,:,1);
        NA_MLT_Qsen_ano(:,:,:,count)=MLT_Qsen_mon0(:,:,:,count)-NA_MLT_Qsen_clim(:,:,:,1);
        count=count+1;
    end
    %clear *_clim
    dmlt_2023_Qnet=NA_MLT_Qnet_ano(:,:,1:12,43);
% #########################################################################



% #########################################################################       
subplot('position',pos{301}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,5)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('i. Heat flux term, Qnet clim.','fontsize',20,'FontWeight','bold')
        
subplot('position',pos{302}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,6)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('j. Heat flux term, Qnet clim.','fontsize',20,'FontWeight','bold')
        
subplot('position',pos{303}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,7)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('k. Heat flux term, Qnet clim.','fontsize',20,'FontWeight','bold')
        
subplot('position',pos{304}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,8)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('l. Heat flux term, Qnet clim.','fontsize',20,'FontWeight','bold')

        clear dmlt_2023_Qnet NA_MLT* MLT*
% #########################################################################
    

hBar1 = colorbar('EastOutside','vertical');
get(hBar1, 'Position');
set(hBar1, 'Position', [ixs+4*ixw+3*ixd+0.012 iys+1.0*iyw+0*iyd 0.012 1.0*iyw+2*iyd]);
set(hBar1, 'ytick',-2.4:0.4:2.4,'yticklabel',{'<-2.4',[],'-1.6',[],'-0.8',[],'0',[],'0.8',[],'1.6',[],'>2.4'},'fontsize',20,'FontName','Arial','LineWidth',1.2,'TickLength',0.055);
ylabel(hBar1, '[ ^oC/month ]','rotation',90);
        

cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')




%% ########################################################################
%  Plotting 5-6: Maps for ML Warming by SWR, May-Aug 2023
clc;clear
time_ann=(1981:2022)';


% #########################################################################
figure('Color',[1 1 1]);  %create a new figure of white color background
ixs = 0.040; ixe = 0.085;  ixd = 0.03; ixw = (1-ixs-ixe-3*ixd)/4;
iys = 0.060; iye = 0.060;  iyd = 0.05; iyw = (1-iys-iye-2*iyd)/3;

ixw./iyw

%          [left            bottom      width height]
pos{101}  = [ixs          iys+2*iyw+2*iyd   ixw 1*iyw]; 
pos{201}  = [ixs          iys+1*iyw+1*iyd   ixw 1*iyw]; 
pos{301}  = [ixs          iys+0*iyw+0*iyd   ixw 1*iyw]; 

pos{102}  = [ixs+1*ixw+1*ixd   iys+2*iyw+2*iyd   ixw 1*iyw]; 
pos{202}  = [ixs+1*ixw+1*ixd   iys+1*iyw+1*iyd   ixw 1*iyw]; 
pos{302}  = [ixs+1*ixw+1*ixd   iys+0*iyw+0*iyd   ixw 1*iyw]; 

pos{103}  = [ixs+2*ixw+2*ixd   iys+2*iyw+2*iyd   ixw 1*iyw]; 
pos{203}  = [ixs+2*ixw+2*ixd   iys+1*iyw+1*iyd   ixw 1*iyw]; 
pos{303}  = [ixs+2*ixw+2*ixd   iys+0*iyw+0*iyd   ixw 1*iyw]; 

pos{104}  = [ixs+3*ixw+3*ixd   iys+2*iyw+2*iyd   ixw 1*iyw]; 
pos{204}  = [ixs+3*ixw+3*ixd   iys+1*iyw+1*iyd   ixw 1*iyw]; 
pos{304}  = [ixs+3*ixw+3*ixd   iys+0*iyw+0*iyd   ixw 1*iyw]; 

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
% #########################################################################
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
    mlt_month_19812010_clim=squeeze(nanmean(mlt_monthly_NA(:,:,:,1:30),4));
    count=1;
    for year=1:size(mlt_monthly_NA,4)
        mlt_month_19812023_ano(:,:,:,count)=mlt_monthly_NA(:,:,:,count)-mlt_month_19812010_clim;
        count=count+1;
    end
    dmlt_month_19812023_ano=mlt_month_19812023_ano(:,:,2:12,:)-mlt_month_19812023_ano(:,:,1:11,:);
    
    % Centered at mid-month
    dmlt_2023(:,:,2:8,1)=(dmlt_month_19812023_ano(:,:,1:7,43)+dmlt_month_19812023_ano(:,:,2:8,43))./2;
    clear mlt* dmlt_month_19812023_ano
% #########################################################################
% #########################################################################



% #########################################################################
% #########################################################################
% Mixed layer warming by air-sea flues in 2023 - accumuated monthly heat flux 
count_yr=1;
for year=1981:2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
    load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Accumulated_',num2str(year),'.mat'],...
          'MLT_Qnet_mon','MLT_Qswr_mon','MLT_Qlat_mon','MLT_Qlon_mon','MLT_Qsen_mon')
    MLT_Qnet_mon0(:,:,:,count_yr)=MLT_Qnet_mon;
    MLT_Qswr_mon0(:,:,:,count_yr)=MLT_Qswr_mon;
    MLT_Qlat_mon0(:,:,:,count_yr)=MLT_Qlat_mon;
    MLT_Qlon_mon0(:,:,:,count_yr)=MLT_Qlon_mon;
    MLT_Qsen_mon0(:,:,:,count_yr)=MLT_Qsen_mon;
    clear MLT_*_mon
    count_yr=count_yr+1;
end
clear count_yr year
% Climatology
    NA_MLT_Qnet_clim=nanmean(MLT_Qnet_mon0(:,:,:,1:30),4);
    NA_MLT_Qswr_clim=nanmean(MLT_Qswr_mon0(:,:,:,1:30),4);
    NA_MLT_Qlat_clim=nanmean(MLT_Qlat_mon0(:,:,:,1:30),4);
    NA_MLT_Qlon_clim=nanmean(MLT_Qlon_mon0(:,:,:,1:30),4);
    NA_MLT_Qsen_clim=nanmean(MLT_Qsen_mon0(:,:,:,1:30),4);
    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(MLT_Qnet_mon0,4)
        NA_MLT_Qnet_ano(:,:,:,count)=MLT_Qnet_mon0(:,:,:,count)-NA_MLT_Qnet_clim(:,:,:,1);
        NA_MLT_Qswr_ano(:,:,:,count)=MLT_Qswr_mon0(:,:,:,count)-NA_MLT_Qswr_clim(:,:,:,1);
        NA_MLT_Qlat_ano(:,:,:,count)=MLT_Qlat_mon0(:,:,:,count)-NA_MLT_Qlat_clim(:,:,:,1);
        NA_MLT_Qlon_ano(:,:,:,count)=MLT_Qlon_mon0(:,:,:,count)-NA_MLT_Qlon_clim(:,:,:,1);
        NA_MLT_Qsen_ano(:,:,:,count)=MLT_Qsen_mon0(:,:,:,count)-NA_MLT_Qsen_clim(:,:,:,1);
        count=count+1;
    end
    %clear *_clim
    dmlt_2023_Qnet=NA_MLT_Qswr_ano(:,:,1:12,43);
% #########################################################################



% #########################################################################       
subplot('position',pos{101}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,5)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('a. SWR term (May 2023)','fontsize',20,'FontWeight','bold')
        
subplot('position',pos{102}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,6)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('b. SWR term (Jun 2023)','fontsize',20,'FontWeight','bold')
        
subplot('position',pos{103}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,7)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('c. SWR term (Jul 2023)','fontsize',20,'FontWeight','bold')
        
subplot('position',pos{104}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,8)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('d. SWR term (Aug 2023)','fontsize',20,'FontWeight','bold')

        clear dmlt_2023_Qnet NA_MLT* MLT*
% #########################################################################
% #########################################################################  



% #########################################################################
% #########################################################################
% Mixed layer warming by air-sea flues in 2023 - accumuated monthly heat flux 
clear dMLT* NA_MLT* 
count_yr=1;
for year=1981:2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
    % test with MLD clim values
    load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Accumulated_',num2str(year),'_2_Clim_MLD.mat'],...
          'MLT_Qnet_mon','MLT_Qswr_mon','MLT_Qlat_mon','MLT_Qlon_mon','MLT_Qsen_mon')
    MLT_Qnet_mon0(:,:,:,count_yr)=MLT_Qnet_mon;
    MLT_Qswr_mon0(:,:,:,count_yr)=MLT_Qswr_mon;
    MLT_Qlat_mon0(:,:,:,count_yr)=MLT_Qlat_mon;
    MLT_Qlon_mon0(:,:,:,count_yr)=MLT_Qlon_mon;
    MLT_Qsen_mon0(:,:,:,count_yr)=MLT_Qsen_mon;
    clear MLT_*_mon
    count_yr=count_yr+1;
end
clear count_yr year
% Climatology
    NA_MLT_Qnet_clim=nanmean(MLT_Qnet_mon0(:,:,:,1:30),4);
    NA_MLT_Qswr_clim=nanmean(MLT_Qswr_mon0(:,:,:,1:30),4);
    NA_MLT_Qlat_clim=nanmean(MLT_Qlat_mon0(:,:,:,1:30),4);
    NA_MLT_Qlon_clim=nanmean(MLT_Qlon_mon0(:,:,:,1:30),4);
    NA_MLT_Qsen_clim=nanmean(MLT_Qsen_mon0(:,:,:,1:30),4);
    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(MLT_Qnet_mon0,4)
        NA_MLT_Qnet_ano(:,:,:,count)=MLT_Qnet_mon0(:,:,:,count)-NA_MLT_Qnet_clim(:,:,:,1);
        NA_MLT_Qswr_ano(:,:,:,count)=MLT_Qswr_mon0(:,:,:,count)-NA_MLT_Qswr_clim(:,:,:,1);
        NA_MLT_Qlat_ano(:,:,:,count)=MLT_Qlat_mon0(:,:,:,count)-NA_MLT_Qlat_clim(:,:,:,1);
        NA_MLT_Qlon_ano(:,:,:,count)=MLT_Qlon_mon0(:,:,:,count)-NA_MLT_Qlon_clim(:,:,:,1);
        NA_MLT_Qsen_ano(:,:,:,count)=MLT_Qsen_mon0(:,:,:,count)-NA_MLT_Qsen_clim(:,:,:,1);
        count=count+1;
    end
    %clear *_clim
    dmlt_2023_Qnet=NA_MLT_Qswr_ano(:,:,1:12,43);
% #########################################################################



% #########################################################################       
subplot('position',pos{201}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,5)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('e. SWR term, MLD clim.','fontsize',20,'FontWeight','bold')
        
subplot('position',pos{202}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,6)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('f. SWR term, MLD clim.','fontsize',20,'FontWeight','bold')
        
subplot('position',pos{203}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,7)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('g. SWR term, MLD clim.','fontsize',20,'FontWeight','bold')
        
subplot('position',pos{204}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,8)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],'xticklabels',[],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('h. SWR term, MLD clim.','fontsize',20,'FontWeight','bold')

        clear dmlt_2023_Qnet NA_MLT* MLT*
% #########################################################################
% #########################################################################


    
    
% #########################################################################
% #########################################################################
% Mixed layer warming by air-sea flues in 2023 - accumuated monthly heat flux 
clear dMLT* NA_MLT* 
count_yr=1;
for year=1981:2023
    cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/plot_MLT_decomposition_ERA5_SST_IAP_V3_Monthly_19812023_V2_1deg')
    % test with Qnet clim values
    load(['Predicted_NA_monthly_MLT_ERA5_IAP_V3_Accumulated_',num2str(year),'_3_Clim_Qnet.mat'],...
          'MLT_Qnet_mon','MLT_Qswr_mon','MLT_Qlat_mon','MLT_Qlon_mon','MLT_Qsen_mon')
    MLT_Qnet_mon0(:,:,:,count_yr)=MLT_Qnet_mon;
    MLT_Qswr_mon0(:,:,:,count_yr)=MLT_Qswr_mon;
    MLT_Qlat_mon0(:,:,:,count_yr)=MLT_Qlat_mon;
    MLT_Qlon_mon0(:,:,:,count_yr)=MLT_Qlon_mon;
    MLT_Qsen_mon0(:,:,:,count_yr)=MLT_Qsen_mon;
    clear MLT_*_mon
    count_yr=count_yr+1;
end
clear count_yr year
% Climatology
    NA_MLT_Qnet_clim=nanmean(MLT_Qnet_mon0(:,:,:,1:30),4);
    NA_MLT_Qswr_clim=nanmean(MLT_Qswr_mon0(:,:,:,1:30),4);
    NA_MLT_Qlat_clim=nanmean(MLT_Qlat_mon0(:,:,:,1:30),4);
    NA_MLT_Qlon_clim=nanmean(MLT_Qlon_mon0(:,:,:,1:30),4);
    NA_MLT_Qsen_clim=nanmean(MLT_Qsen_mon0(:,:,:,1:30),4);
    
% Anomaly relative to the clim-monthly mean
    count=1;
    for year=1:size(MLT_Qnet_mon0,4)
        NA_MLT_Qnet_ano(:,:,:,count)=MLT_Qnet_mon0(:,:,:,count)-NA_MLT_Qnet_clim(:,:,:,1);
        NA_MLT_Qswr_ano(:,:,:,count)=MLT_Qswr_mon0(:,:,:,count)-NA_MLT_Qswr_clim(:,:,:,1);
        NA_MLT_Qlat_ano(:,:,:,count)=MLT_Qlat_mon0(:,:,:,count)-NA_MLT_Qlat_clim(:,:,:,1);
        NA_MLT_Qlon_ano(:,:,:,count)=MLT_Qlon_mon0(:,:,:,count)-NA_MLT_Qlon_clim(:,:,:,1);
        NA_MLT_Qsen_ano(:,:,:,count)=MLT_Qsen_mon0(:,:,:,count)-NA_MLT_Qsen_clim(:,:,:,1);
        count=count+1;
    end
    %clear *_clim
    dmlt_2023_Qnet=NA_MLT_Qswr_ano(:,:,1:12,43);
% #########################################################################



% #########################################################################       
subplot('position',pos{301}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,5)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('i. SWR term, SWR clim.','fontsize',20,'FontWeight','bold')
        
subplot('position',pos{302}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,6)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('j. SWR term, SWR clim.','fontsize',20,'FontWeight','bold')
        
subplot('position',pos{303}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,7)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('k. SWR term, SWR clim.','fontsize',20,'FontWeight','bold')
        
subplot('position',pos{304}) 
        m_proj('miller','lon',[260.5 379.5],'lat',[0.5 69.5]);  
        m_pcolor(lon_IAP,lat_IAP,dmlt_2023_Qnet(:,:,8)');
        shading flat
        
        hold on       
        m_grid('box','on','tickdir','in','linestyle','none','ticklen',0.01,'fontsize',20,...
                    'xtick',[20:30:380],...
                    'yaxislocation','left','ytick',[-60:20:60],'yticklabels',[]);
        m_coast('patch',[0.65 0.65 0.65],'edgecolor',[0.4 0.4 0.4]); 
        hold on;
        colormap(gca,color0)
        caxis([-2.4 2.4]);
        title('l. SWR term, SWR clim.','fontsize',20,'FontWeight','bold')

        clear dmlt_2023_Qnet NA_MLT* MLT*
% #########################################################################
    

hBar1 = colorbar('EastOutside','vertical');
get(hBar1, 'Position');
set(hBar1, 'Position', [ixs+4*ixw+3*ixd+0.012 iys+1.0*iyw+0*iyd 0.012 1.0*iyw+2*iyd]);
set(hBar1, 'ytick',-2.4:0.4:2.4,'yticklabel',{'<-2.4',[],'-1.6',[],'-0.8',[],'0',[],'0.8',[],'1.6',[],'>2.4'},'fontsize',20,'FontName','Arial','LineWidth',1.2,'TickLength',0.055);
ylabel(hBar1, '[ ^oC/month ]','rotation',90);
        

cd('/Users/z5195509/Documents/6_NA_MHW_MLD/2_1_MLT_Decomposition/')


