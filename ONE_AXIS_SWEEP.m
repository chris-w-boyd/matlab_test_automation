%RUN THE SWEEP SETUP PROGRAM BEFORE RUNNING THE FIRST SECTION OF THIS
%PROGRAM 
date = datetime('now','Format','MM_dd_yyyy_h_mm_a');

%====================================================
if(forward_and_reverse == 0)
    if(sweep_which_axis == 0)
        sigout_ch = 1;
        [frequencies180,r180_meas,theta180] = main_ZI_sweep(sigout_ch,amplitude(1),num_points,TC_settle,settle_time,start_f,stop_f,samp_avg,scan_order,aux1_out,aux2_out,aux3_out ...
            ,settings_filename,filter_slope,filter_BW,demod_samp_rate,diff_sense,sense_range);
    elseif(sweep_which_axis == 1)
        sigout_ch = 2;
        [frequencies225,r225_meas,theta225] = main_ZI_sweep(sigout_ch,amplitude(2),num_points,TC_settle,settle_time,start_f,stop_f,samp_avg,scan_order,aux1_out,aux2_out,aux3_out ...
             ,settings_filename,filter_slope,filter_BW,demod_samp_rate,diff_sense,sense_range);    
    elseif(sweep_which_axis == 2)
        sigout_ch = 1;
        [frequencies180,r180_meas,theta180] = main_ZI_sweep(sigout_ch,amplitude(1),num_points,TC_settle,settle_time,start_f,stop_f,samp_avg,scan_order,aux1_out,aux2_out,aux3_out ...
             ,settings_filename,filter_slope,filter_BW,demod_samp_rate,diff_sense,sense_range);
        pause(pause_time);
        sigout_ch = 2;
        [frequencies225,r225_meas,theta225] = main_ZI_sweep(sigout_ch,amplitude(2),num_points,TC_settle,settle_time,start_f,stop_f,samp_avg,scan_order,aux1_out,aux2_out,aux3_out ...
             ,settings_filename,filter_slope,filter_BW,demod_samp_rate,diff_sense,sense_range);
    end   
else %do a sweep up and then down in frequency 
    %forward sweep (drive axis 1)
    sigout_ch = 1;
    frequencies180 = nan(2,num_points);
    r180_meas = nan(4,num_points);
    theta180 = nan(4,num_points);
    [frequencies180(1,:),r180_meas(1:2,:),theta180(1:2,:)] = main_ZI_sweep(sigout_ch,amplitude(1),num_points,TC_settle,settle_time,start_f,stop_f,samp_avg,0);
    pause(pause_time)
    %reverse sweep (drive axis 1)
    [frequencies180(2,:),r180_meas(3:4,:),theta180(3:4,:)] = main_ZI_sweep(sigout_ch,amplitude(1),num_points,TC_settle,settle_time,start_f,stop_f,samp_avg,3);
    
    %pause(10)
    %{
    sigout_ch = 2;
    frequencies225 = nan(2,num_points);
    r225 = nan(4,num_points);
    theta225 = nan(4,num_points);
    %forward sweep (drive axis 2)
    [frequencies225(1,:),r225(1:2,:),theta225(1:2,:)] = main_ZI_sweep(sigout_ch,amplitude,num_points,TC_settle,settle_time,start_f,stop_f,samp_avg,0);
    pause(10)
    %reverse sweep (drive axis 2) 
    [frequencies225(2,:),r225(3:4,:),theta225(3:4,:)] = main_ZI_sweep(sigout_ch,amplitude,num_points,TC_settle,settle_time,start_f,stop_f,samp_avg,3);
    %}
end
%% plot data / find resonance frequencies and amplitudes
%============================================
%compensate for 10x attentuator on the board
ise180 = evalin( 'base', 'exist(''r180_meas'',''var'') == 1' );
ise225 = evalin( 'base', 'exist(''r225_meas'',''var'') == 1' );
if(ise180)
    r180 = 10.*r180_meas;
end

if(ise225)
    r225 = 10.*r225_meas;    
end
    
%}
%=============================================
if(forward_and_reverse == 0)
    if(ise180 && (~ise225))
        plot_fr_data(frequencies180,r180,theta180,'-',1,0,2,0,1);
    elseif((~ise180) && (ise225))
        plot_fr_data(frequencies225,r225,theta225,'-',2,0,2,0,1);
    elseif((ise180) && (ise225))
        plot_fr_data(frequencies180,r180,theta180,'-',1,0,3,0,1);
        plot_fr_data(frequencies225,r225,theta225,'-',2,0,4,0,1);
        plot_fr_data(frequencies180,[r180(1,:);r225(1,:)],[theta180(1,:);theta225(1,:)],'-',3,0,5,0,0);
        plot_fr_data(frequencies180,[r180(2,:);r225(2,:)],[theta180(2,:);theta225(2,:)],'-',4,0,6,0,0);
       
        %combined plot (no phase)
        plot_fr_data(frequencies180, [r180(1,:);r225(1,:)],[theta180(1,:);theta225(1,:)],'-',3,0,7,0,0);
        plot_fr_data(frequencies180, [r180(2,:);r225(2,:)],[theta180(2,:);theta225(2,:)],'-',4,0,7,0,0);
    
        plot_fr_data(frequencies180, [r180(1,:);r225(1,:)],[theta180(1,:);theta225(1,:)],'-',3,0,8,1,0);
        plot_fr_data(frequencies180, [r180(2,:);r225(2,:)],[theta180(2,:);theta225(2,:)],'-',4,0,8,1,0);
        
    end
else
    %{
    plot_fr_data(frequencies180(1,:),[r180(1,:);r225(1,:)],[theta180(1,:);theta225(1,:)],'-',3,1,3,0);
    plot_fr_data(frequencies180(2,:),[r180(3,:);r225(3,:)],[theta180(3,:);theta225(3,:)],'-',5,0,3,0);

    plot_fr_data(frequencies180(1,:),[r180(2,:);r225(2,:)],[theta180(2,:);theta225(2,:)],'-',4,1,4,0);
    plot_fr_data(frequencies180(2,:),[r180(4,:);r225(4,:)],[theta180(4,:);theta225(4,:)],'-',6,0,4,0);

    plot_fr_data(frequencies180(1,:),[r180(1,:);r225(1,:)],[theta180(1,:);theta225(1,:)],'-',3,0,5,1);
    plot_fr_data(frequencies180(2,:),[r180(3,:);r225(3,:)],[theta180(3,:);theta225(3,:)],'-',5,0,5,1);

    plot_fr_data(frequencies180(1,:),[r180(2,:);r225(2,:)],[theta180(2,:);theta225(2,:)],'-',4,0,6,1);
    plot_fr_data(frequencies180(2,:),[r180(4,:);r225(4,:)],[theta180(4,:);theta225(4,:)],'-',6,0,6,1);  
    %}
end

if( (sweep_which_axis == 0) || (sweep_which_axis == 2))
    df = mean(diff(frequencies180(1,:)))*1000; %resolution in mHz
else
    df = mean(diff(frequencies225(1,:)))*1000; %resolution in mHz
end

numPeaks = 4;
minPeakProminence = 2;%10; %10 is a good starting point
%text_offset = 0.05; %0.5 is good for 20-30Hz sweep span , use smaller offset as span decreases
text_offset = 0.025;

if(forward_and_reverse == 0)
    if(scan_order == 0)
        %(7/16/19) I only modified code for case of a single forward scan
        %find resonance peaks for in drive180/sense0 (in-axis)
        
        if(ise180)
            [pks1,locs1,w1,p1] = findpeaks(20*log10(r180(1,:)),frequencies180,'Npeaks',numPeaks,'SortStr','descend',...
            'MinPeakProminence',minPeakProminence);
        
            %find resonance peaks for drive180/sense45 (cross-axis)
            [pks3,locs3,w3,p3] = findpeaks(20*log10(r180(2,:)),frequencies180,'Npeaks',numPeaks,'SortStr','descend',...
            'MinPeakProminence',minPeakProminence);
        end
    
        if(ise225)
            %find resonance peaks for in drive225/sense45 (in-axis)
            [pks2,locs2,w2,p2] = findpeaks(20*log10(r225(1,:)),frequencies225,'Npeaks',numPeaks,'SortStr','descend',...
            'MinPeakProminence',minPeakProminence);
            
            %find resonance peaks for drive225/sense0 (cross-axis)
            [pks4,locs4,w4,p4] = findpeaks(20*log10(r225(2,:)),frequencies225,'Npeaks',numPeaks,'SortStr','descend',...
            'MinPeakProminence',minPeakProminence);
        end
        
    elseif(scan_order == 3)
        %find resonance peaks for in drive180/sense0 (in-axis)
        [pks1,locs1,w1,p1] = findpeaks(20*log10(fliplr(r180(1,:))),fliplr(frequencies180),'Npeaks',numPeaks,'SortStr','descend',...
        'MinPeakProminence',minPeakProminence);

        %find resonance peaks for in drive225/sense45 (in-axis)
        [pks2,locs2,w2,p2] = findpeaks(20*log10(fliplr(r225(1,:))),fliplr(frequencies225),'Npeaks',numPeaks,'SortStr','descend',...
        'MinPeakProminence',minPeakProminence);

        %find resonance peaks for drive180/sense45 (cross-axis)
        [pks3,locs3,w3,p3] = findpeaks(20*log10(fliplr(r180(2,:))),fliplr(frequencies180),'Npeaks',numPeaks,'SortStr','descend',...
        'MinPeakProminence',minPeakProminence);

        %find resonance peaks for drive225/sense0 (cross-axis)
        [pks4,locs4,w4,p4] = findpeaks(20*log10(fliplr(r225(2,:))),fliplr(frequencies225),'Npeaks',numPeaks,'SortStr','descend',...
        'MinPeakProminence',minPeakProminence);     
    else
       %DO NOTHING  
    end
else
    %find forward scan peaks 
    %find resonance peaks for in drive180/sense0 (in-axis)
        [pks1,locs1,w1,p1] = findpeaks(20*log10(r180(1,:)),frequencies180(1,:),'Npeaks',numPeaks,'SortStr','descend',...
        'MinPeakProminence',minPeakProminence);

        %find resonance peaks for in drive225/sense45 (in-axis)
        [pks2,locs2,w2,p2] = findpeaks(20*log10(r225(1,:)),frequencies225(1,:),'Npeaks',numPeaks,'SortStr','descend',...
        'MinPeakProminence',minPeakProminence);

        %find resonance peaks for drive180/sense45 (cross-axis)
        [pks3,locs3,w3,p3] = findpeaks(20*log10(r180(2,:)),frequencies180(1,:),'Npeaks',numPeaks,'SortStr','descend',...
        'MinPeakProminence',minPeakProminence);

        %find resonance peaks for drive225/sense0 (cross-axis)
        [pks4,locs4,w4,p4] = findpeaks(20*log10(r225(2,:)),frequencies225(1,:),'Npeaks',numPeaks,'SortStr','descend',...
        'MinPeakProminence',minPeakProminence);
    
     %find reverse peaks
     %find resonance peaks for in drive180/sense0 (in-axis)
        [pks5,locs5,w5,p5] = findpeaks(20*log10(fliplr(r180(3,:))),fliplr(frequencies180(2,:)),'Npeaks',numPeaks,'SortStr','descend',...
        'MinPeakProminence',minPeakProminence);

        %find resonance peaks for in drive225/sense45 (in-axis)
        [pks6,locs6,w6,p6] = findpeaks(20*log10(fliplr(r225(3,:))),fliplr(frequencies225(2,:)),'Npeaks',numPeaks,'SortStr','descend',...
        'MinPeakProminence',minPeakProminence);

        %find resonance peaks for drive180/sense45 (cross-axis)
        [pks7,locs7,w7,p7] = findpeaks(20*log10(fliplr(r180(4,:))),fliplr(frequencies180(2,:)),'Npeaks',numPeaks,'SortStr','descend',...
        'MinPeakProminence',minPeakProminence);

        %find resonance peaks for drive225/sense0 (cross-axis)
        [pks8,locs8,w8,p8] = findpeaks(20*log10(fliplr(r225(4,:))),fliplr(frequencies225(2,:)),'Npeaks',numPeaks,'SortStr','descend',...
        'MinPeakProminence',minPeakProminence);     
end 

%place markers on the plots to identify peaks 

 if(ise180 && (~ise225))

    figure(2)
    hold on
    plot(locs1,pks1,'v','MarkerFaceColor','g');
    text(locs1+text_offset,pks1,num2str((1:numel(pks1))'))

    plot(locs3,pks3,'v','MarkerFaceColor','m');
    text(locs3+text_offset,pks3,num2str((1:numel(pks3))'))

 elseif((~ise180) && (ise225))
     
     figure(2)
     hold on
     plot(locs2,pks2,'v','MarkerFaceColor','b');
     text(locs2+text_offset,pks2,num2str((1:numel(pks2))'))

     plot(locs4,pks4,'v','MarkerFaceColor','r');
     text(locs4+text_offset,pks4,num2str((1:numel(pks4))'))
     
     
 elseif((ise180) && (ise225))
     
     figure(3)
     hold on
     plot(locs1,pks1,'v','MarkerFaceColor','g');
     text(locs1+text_offset,pks1,num2str((1:numel(pks1))'))

     plot(locs3,pks3,'v','MarkerFaceColor','m');
     text(locs3+text_offset,pks3,num2str((1:numel(pks3))'))
     
     figure(4)
     hold on
     plot(locs2,pks2,'v','MarkerFaceColor','b');
     text(locs2+text_offset,pks2,num2str((1:numel(pks2))'))

     plot(locs4,pks4,'v','MarkerFaceColor','r');
     text(locs4+text_offset,pks4,num2str((1:numel(pks4))'))
     
     figure(5)
     hold on
     plot(locs1,pks1,'v','MarkerFaceColor','g');
     text(locs1+text_offset,pks1,num2str((1:numel(pks1))'))
     plot(locs2,pks2,'v','MarkerFaceColor','b');
     text(locs2+text_offset,pks2,num2str((1:numel(pks2))'))
     
     figure(6)
     hold on
     plot(locs3,pks3,'v','MarkerFaceColor','m');
     text(locs3+text_offset,pks3,num2str((1:numel(pks3))'))
     plot(locs4,pks4,'v','MarkerFaceColor','r');
     text(locs4+text_offset,pks4,num2str((1:numel(pks4))'))
     
     figure(7)
     hold on
     plot(locs1(1),pks1(1),'v','MarkerFaceColor','g');
     text(locs1(1)+text_offset,pks1(1),num2str(1))
     plot(locs2(1),pks2(1),'v','MarkerFaceColor','b');
     text(locs2(1)+text_offset,pks2(1),num2str(1))
     
     plot(locs3(1),pks3(1),'v','MarkerFaceColor','m');
     text(locs3(1)+text_offset,pks3(1),num2str(1));
     plot(locs4(1),pks4(1),'v','MarkerFaceColor','r');
     text(locs4(1)+text_offset,pks4(1),num2str(1))
     
     
 end
    
  
    
%{
figure(3)
hold on
plot(locs1,pks1,'v','MarkerFaceColor','g');
text(locs1+text_offset,pks1,num2str((1:numel(pks1))'))

plot(locs2,pks2,'v','MarkerFaceColor','b');
text(locs2+text_offset,pks2,num2str((1:numel(pks2))'))

if(forward_and_reverse == 1) 
    plot(locs5,pks5,'v','MarkerFaceColor','c');
    text(locs5+text_offset,pks5,num2str((1:numel(pks5))'))

    plot(locs6,pks6,'v','MarkerFaceColor','k');
    text(locs6+text_offset,pks6,num2str((1:numel(pks6))'))
end

figure(4)
hold on
plot(locs3,pks3,'v','MarkerFaceColor','m');
text(locs3+text_offset,pks3,num2str((1:numel(pks3))'))

plot(locs4,pks4,'v','MarkerFaceColor','r');
text(locs4+text_offset,pks4,num2str((1:numel(pks4))'))

if(forward_and_reverse == 1) 
    plot(locs7,pks7,'v','MarkerFaceColor','y');
    text(locs7+text_offset,pks7,num2str((1:numel(pks7))'))

    plot(locs8,pks8,'v','MarkerFaceColor','k');
    text(locs8+text_offset,pks8,num2str((1:numel(pks8))'))
end
%}

fh = fopen(file_name,'a+');

%compute deltaf
 if(ise180 && (~ise225))
      in_axis_delta_f = (locs1(1)-locs3(1)); %Hz
 elseif((~ise180) && (ise225))
      in_axis_delta_f = (locs2(1)-locs4(1)); %Hz=
 elseif((ise180) && (ise225))
     in_axis_delta_f = (locs1(1)-locs2(1)); %Hz
 end
     
    


axis_name = cell(1,4);
axis_name{1,1} = 'DRIVE_180/SENSE_0 (IN-AXIS)';
axis_name{1,2} = 'DRIVE_225/SENSE_45 (IN-AXIS)';
axis_name{1,3} = 'DRIVE_180/SENSE_45 (CROSS-AXIS)';
axis_name{1,4} = 'DRIVE_225/SENSE_0 (CROSS-AXIS)';

fprintf('==============================================================\n');
fprintf(fh,'==============================================================\n');
fprintf('STEP# %d\n',step_num);
fprintf(fh,'STEP# %d\n',step_num);
fprintf(strcat(char(date),'\n'));
fprintf(fh,strcat(char(date),'\n'));
fprintf('START_F(Hz): %g STOP_F(Hz): %g NUMPOINTS: %g\n',start_f,stop_f,num_points);
fprintf(fh,'START_F(Hz): %g STOP_F(Hz): %g NUMPOINTS: %g\n',start_f,stop_f,num_points);
fprintf('DRIVE_AMP (V): AX1 %g AX2 %g SETTLE_TIME(s): %g\n',amplitude(1),amplitude(2),settle_time);
fprintf(fh,'DRIVE_AMP (V): AX1 %g AX2 %g SETTLE_TIME(s): %g\n',amplitude(1),amplitude(2),settle_time);
fprintf('TC_SETTLE: %g SAMPLE_AVG: %g\n',TC_settle,samp_avg);
fprintf(fh,'TC_SETTLE: %g SAMPLE_AVG: %g\n',TC_settle,samp_avg);
fprintf('FILTER_SLOPE(dB/Oct): %d FILTER_BW(Hz): %g DEMOD_STREAM_RATE(Hz): %g\n',filter_slope,filter_BW,demod_samp_rate);
fprintf(fh,'FILTER_SLOPE(dB/Oct): %d FILTER_BW(Hz): %g DEMOD_STREAM_RATE(Hz): %g\n',filter_slope,filter_BW,demod_samp_rate);
fprintf(strcat(comment1,'\n'));
fprintf(fh,strcat(comment2,'\n'));
fprintf(strcat(comment2,'\n'));
fprintf(fh,strcat(comment1,'\n'));
fprintf('IN-AXIS delta f(Hz): %g IN-AXIS delta_f(mHz): %g \n',in_axis_delta_f,in_axis_delta_f*1000);
fprintf(fh,'IN-AXIS delta f(Hz): %g IN-AXIS delta_f(mHz): %g \n',in_axis_delta_f,in_axis_delta_f*1000);
fprintf('Frequency Resolution(mHz): %f\n',df);
fprintf(fh,'Frequency Resolution(mHz): %f\n',df);

if(forward_and_reverse == 1)
    for(i=1:4)
        fprintf('%s\n',axis_name{1,i});
        fprintf(fh,'%s\n',axis_name{1,i});
        if(i==1)
            pks_f = pks1;
            locs_f = locs1;
            p_f = p1;
            
            pks_r = pks5;
            locs_r = locs5;
            p_r = p5;
        
        elseif(i==2)
            pks_f = pks2;
            locs_f = locs2;
            p_f = p2;
            
            pks_r = pks6;
            locs_r = locs6;
            p_r = p6;
        
        elseif(i==3)
            pks_f = pks3;
            locs_f = locs3;
            p_f = p3;
            
            pks_r = pks7;
            locs_r = locs7;
            p_r = p7;
        elseif(i==4)
            pks_f = pks4;
            locs_f = locs4;
            p_f = p4;
            
            pks_r = pks8;
            locs_r = locs8;
            p_r = p8;
        end
        
        fprintf('FORWARD SWEEP\n');
        fprintf(fh,'FORWARD SWEEP\n');
        
        for(ind=1:numel(pks_f))
            fprintf('Peak: %d Freq(Hz): %f Amp(dBV): %g Prom: %g\n',ind,locs_f(ind),pks_f(ind),p_f(ind));
            fprintf(fh,'Peak: %d Freq(Hz): %f Amp(dBV): %g Prom: %g\n',ind,locs_f(ind),pks_f(ind),p_f(ind)); 
        end
        
        fprintf('REVERSE SWEEP\n');
        fprintf(fh,'REVERSE SWEEP\n');
        
        for(ind=1:numel(pks_r))
            fprintf('Peak: %d Freq(Hz): %f Amp(dBV): %g Prom: %g\n',ind,locs_r(ind),pks_r(ind),p_r(ind));
            fprintf(fh,'Peak: %d Freq(Hz): %f Amp(dBV): %g Prom: %g\n',ind,locs_r(ind),pks_r(ind),p_r(ind)); 
        end   
    end
else  %(forward_and_reverse == 0) (forward sweep) 
    
    if((ise180) && (~ise225))
        test_cond = [1 3];
    elseif((~ise180) && (ise225))
        test_cond = [2 4];
    elseif((ise180) && (ise225))
        test_cond = [1 3 2 4];
    end
 
    for(i=test_cond) %modified on 7/16/19 to only do drive 180 (sense0/45)
        fprintf('%s\n',axis_name{1,i});
        fprintf(fh,'%s\n',axis_name{1,i});
        if(i==1)
            pks = pks1;
            locs = locs1;
            p = p1;
        elseif(i==2)
            pks = pks2;
            locs = locs2;
            p = p2;
        elseif(i==3)
            pks = pks3;
            locs = locs3;
            p = p3;
        elseif(i==4)
            pks = pks4;
            locs = locs4;
            p = p4;
        end
    
        for(ind=1:numel(pks))
            fprintf('Peak: %d Freq(Hz): %f Amp(dBV): %g Prom: %g\n',ind,locs(ind),pks(ind),p(ind));
            fprintf(fh,'Peak: %d Freq(Hz): %f Amp(dBV): %g Prom: %g\n',ind,locs(ind),pks(ind),p(ind)); 
        end
    end
end
fclose(fh);
%}
%===================================================================================================================

%%  Calculate delta f and delta amplitude (used to measure quadrature error) 

%date = datetime('now','Format','MM-dd-yyyy-h:mm-a');
%save figures!
savepath1 = 'F:\Users\cwboyd\Desktop\DATA_2021\May_Sept_21\gain_phase_logs\plots\7_13_21\';
savepath2 = 'F:\Users\cwboyd\Desktop\DATA_2021\May_Sept_21\gain_phase_logs\raw_data\7_13_21\';

str1 = sprintf('step%d',step_num);

if((ise180) && (~ise225))

    name1 = strcat(str1,'_',char(date),'_GP_drive180');     %sprintf('step%d_gain_inaxis.png',step_num);
    name2 = strcat(str1,'_',char(date),'_GAIN_drive180');   %sprintf('step%d_gain_crossaxis.png',step_num);

    saveas(1,strcat(savepath1,name1),'png');
    saveas(2,strcat(savepath1,name2),'png');

    name7 = strcat(str1,'_',char(date),'_raw_data','.mat');

    save(strcat(savepath2,name7),'frequencies180','r180','theta180');
    
elseif((~ise180) && (ise225))
    
    name1 = strcat(str1,'_',char(date),'_GP_drive225');     %sprintf('step%d_gain_inaxis.png',step_num);
    name2 = strcat(str1,'_',char(date),'_GAIN_drive225');  %sprintf('step%d_gain_crossaxis.png',step_num);
    
    saveas(1,strcat(savepath1,name1),'png');
    saveas(2,strcat(savepath1,name2),'png');

    name7 = strcat(str1,'_',char(date),'_raw_data','.mat');

    save(strcat(savepath2,name7),'frequencies225','r225','theta225');
    
elseif((ise180) && (ise225))
    
    name1 = strcat(str1,'_',char(date),'_GP_drive180');     %sprintf('step%d_gain_inaxis.png',step_num);
    name2 = strcat(str1,'_',char(date),'_GAIN_drive180');  %sprintf('step%d_gain_crossaxis.png',step_num);
    name3 = strcat(str1,'_',char(date),'_GP_drive225');     %sprintf('step%d_gain_inaxis.png',step_num);
    name4 = strcat(str1,'_',char(date),'_GAIN_drive225');  %sprintf('step%d_gain_crossaxis.png',step_num);
    name5 = strcat(str1,'_',char(date),'_GAIN_INAXIS');  %sprintf('step%d_gain_crossaxis.png',step_num);
    name6 = strcat(str1,'_',char(date),'_GAIN_CROSSAXIS');  %sprintf('step%d_gain_crossaxis.png',step_num);
    name7 = strcat(str1,'_',char(date),'_GAIN_4CHANNEL');
    name8 = strcat(str1,'_',char(date),'_GP_4CHANNEL');
    
    saveas(1,strcat(savepath1,name1),'png');
    saveas(2,strcat(savepath1,name3),'png');
    saveas(3,strcat(savepath1,name2),'png');
    saveas(4,strcat(savepath1,name4),'png');
    saveas(5,strcat(savepath1,name5),'png');
    saveas(6,strcat(savepath1,name6),'png');
    saveas(7,strcat(savepath1,name7),'png');
    saveas(8,strcat(savepath1,name8),'png');

    name9 = strcat(str1,'_',char(date),'_raw_data','.mat');

    save(strcat(savepath2,name9),'frequencies180','r180','theta180','frequencies225','r225','theta225');
   
end
    
%name3 = strcat(str1,'_GP_inaxis_',char(date));  %sprintf('step%d_GP_inaxis.png',step_num);
%name4 = strcat(str1,'_GP_crossaxis_',char(date));  %sprintf('step%d_GP_crossaxis.png',step_num);
%name5 = strcat(str1,'_GP_drive180_',char(date));  %sprintf('step%d_GP_drive180.png',step_num);
%name6 = strcat(str1,'_GP_drive225_',char(date));  %sprintf('step%d_GP_drive225.png',step_num);


%saveas(5,strcat(savepath1,name3),'png');
%saveas(6,strcat(savepath1,name4),'png');
%saveas(1,strcat(savepath1,name5),'png');
%saveas(2,strcat(savepath1,name6),'png');




%{
npks1 = input('Number valid in-axis peaks:');
npks2 = input('Number valid cross-axis peaks:');
%}

%{
%assuming peak prominence identifies peaks in correct order 
delta_f = abs(locs1(1)-locs2(1));
delta_f_prim_0_45 = abs(locs1(1)-locs2(2));
delta_amp1 = pks1(1)-pks2(1);
delta_amp2 = pks1(1)-pks2(2);
delta_phase = (theta(1,pk1_index(1))-theta(2,pk2_index(1)))*(180/pi);
phase0deg = theta(1,pk1_index(1))*(180/pi);
phase45deg = theta(2,pk1_index(1))*(180/pi);
fprintf('Calculated delta f and delta amp\n');
fprintf('delta f(Hz): %g delta_f(mHz) %g \n',delta_f,delta_f*1000);
fprintf('delta amp1(dB): %g delta amp2(dB) %g \n',delta_amp1,delta_amp2);
fprintf('delta_f_prim_0_45(mHz): %g \n',delta_f_prim_0_45*1000);
fprintf('Phase difference at peak freq (deg): %g \n',delta_phase);
fprintf('Phase(0sense): %g Phase(45sense): %g \n',phase0deg,phase45deg);
fprintf('amp45(f1): %g\n',20*log10(r(2,pk1_index(1))));
%}
















