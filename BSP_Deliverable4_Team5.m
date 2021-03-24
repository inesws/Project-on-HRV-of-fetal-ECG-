%CTG signal: Healthy - 34.78+/-0.53 ; IUGR - 32.27 +/- 2.79
load Healthy_1.mat ;Healthy_1 = data;
load Healthy_2.mat ;Healthy_2 = data;
load Healthy_3.mat ;Healthy_3 = data;
load Healthy_4.mat ;Healthy_4 = data;
load Healthy_5.mat ;Healthy_5 = data;

load IUGR_1.mat;IUGR_1 = data;
load IUGR_2.mat;IUGR_2 = data;
load IUGR_3.mat;IUGR_3 = data;
load IUGR_4.mat;IUGR_4 = data;
load IUGR_5.mat;IUGR_5 = data;

%Plot all the raw data

% figure
% subplot(3,2,1), plot(Healthy_1), title('Healthy\_1'), xlabel('beat number');
% ylabel('FHR bpm');
% subplot(3,2,2), plot(Healthy_2), title('Healthy\_2'), xlabel('beat number');
% ylabel('FHR bpm');
% subplot(3,2,3), plot(Healthy_3), title('Healthy\_3'), xlabel('beat number');
% ylabel('FHR bpm');
% subplot(3,2,4), plot(Healthy_4), title('Healthy\_4'), xlabel('beat number');
% ylabel('FHR bpm');
% subplot(3,2,5), plot(Healthy_5), title('Healthy\_5'), xlabel('beat number');
% ylabel('FHR bpm');
% figure
% subplot(3,2,1), plot(IUGR_1), title('IUGR\_1'), xlabel('beat number');
% ylabel('FHR bpm');
% subplot(3,2,2), plot(IUGR_2), title('IUGR\_2'), xlabel('beat number');
% ylabel('FHR bpm');
% subplot(3,2,3), plot(IUGR_3), title('IUGR\_3'), xlabel('beat number');
% ylabel('FHR bpm');
% subplot(3,2,4), plot(IUGR_4), title('IUGR\_4'), xlabel('beat number');
% ylabel('FHR bpm');
% subplot(3,2,5), plot(IUGR_5), title('IUGR\_5'), xlabel('beat number');
% ylabel('FHR bpm');

% %Compare healthy_1 with IUGR_1
% figure
% subplot(2,1,1), plot(Healthy_1), title('Healthy\_1'), xlabel('beat number');
% ylabel('FHR bpm');
% subplot(2,1,2), plot(IUGR_1), title('IUGR\_1'), xlabel('beat number');
% ylabel('FHR bpm');

%Sample at 2Hz 
fs=2;
t_Healthy_1 = ( 0:1/fs:length(Healthy_1)/fs-1/fs); 
t_Healthy_2 = 0: 1/fs :length(Healthy_2)/fs - 1/fs ;
t_Healthy_3 = 0:1/fs :length(Healthy_3)/fs - 1/fs ;
t_Healthy_4 =0:1/fs :length(Healthy_4)/fs - 1/fs ;
t_Healthy_5 =0:1/fs :length(Healthy_5)/fs - 1/fs;

t_IUGR_1 = 0:1/fs: length(IUGR_1)/fs - 1/fs ;
t_IUGR_2 = 0:1/fs:length(IUGR_2)/fs - 1/fs ;
t_IUGR_3 = 0:1/fs:length(IUGR_3)/fs - 1/fs ;
t_IUGR_4 =0:1/fs:length(IUGR_4)/fs - 1/fs;
t_IUGR_5 =0:1/fs:length(IUGR_5)/fs - 1/fs ;


figure;
subplot(2,3,1); plot(t_Healthy_1,Healthy_1);title('Healthy\_1 CTG'); xlabel('Time (s)'); ylabel('FHR (bpm)');
subplot(2,3,2); plot(t_Healthy_2,Healthy_2);title('Healthy\_2 CTG'); xlabel('Time (s)'); ylabel('FHR (bpm)');
subplot(2,3,3); plot(t_Healthy_3,Healthy_3);title('Healthy\_3 CTG'); xlabel('Time (s)'); ylabel('FHR (bpm)');
subplot(2,3,4); plot(t_Healthy_4,Healthy_4);title('Healthy\_4 CTG'); xlabel('Time (s)'); ylabel('FHR (bpm)');
subplot(2,3,5); plot(t_Healthy_5,Healthy_5);title('Healthy\_5 CTG'); xlabel('Time (s)'); ylabel('FHR (bpm)');

figure;
subplot(2,3,1); plot(t_IUGR_1,IUGR_1);title('IUGR\_1 CTG');xlabel('Time (s)'); ylabel('FHR (bpm)');
subplot(2,3,2); plot(t_IUGR_2,IUGR_2);title('IUGR\_2 CTG');xlabel('Time (s)'); ylabel('FHR (bpm)');
subplot(2,3,3); plot(t_IUGR_3,IUGR_3);title('IUGR\_3 CTG');xlabel('Time (s)'); ylabel('FHR (bpm)');
subplot(2,3,4); plot(t_IUGR_4,IUGR_4);title('IUGR\_4 CTG');xlabel('Time (s)'); ylabel('FHR (bpm)');
subplot(2,3,5); plot(t_IUGR_5,IUGR_5);title('IUGR\_5 CTG');xlabel('Time (s)'); ylabel('FHR (bpm)');

%% 1)PREPROCESSING, SUBDIVIDING IN SUBINTERVALS HEALTHY SIGNALS

%a)Detecting bad quality intervals (size = windo)
%Bad quality intervals: Defined here as intervals where inside a window there are too much
%missing data as a whole or it appears consecutive above a threshold
%Algorithm: Deleted (by putting samples to NaN) if more than 5percent of missing data or spikes
%or more than 8 consecutive bad quality points are identify. A window to
%analyse the subinterval is choosen according to the feature to analyse
%forward.
%windo=360 samples (180s) was choosen for frequency domain features
%windo=120(60s) was choosen to analyse time domain features

%b)Analysing the remaining intervals 
%Identification of outliers/artifacts and correction of those
%with a moving average proceeding 
%Algorithm: Using an artifact detection method from Van Geijn, H. P., Jongsma, H. W., De Haan, J., & Eskes, T. K. A. B. (1980). Analysis of heart rate and beat-to-beat variability: Interval difference index. American Journal of Obstetrics and Gynecology, 138(3), 246???252. https://doi.org/10.1016/0002-9378(80)90242-2
%The identified points are then substituted with a moving average
%proceeding that uses a window of 5 previous and forward samples.

windo = 360;


%HEALTHY 1 

%Cutting intervals by putting those to NaN values
Healthy_1_corr = Healthy_1;
count = 0;
percent = 0;
index = 1;
i = 1;
consec_8=0;
while i <= length(Healthy_1) - (windo + length(Healthy_1)-floor(length(Healthy_1)/windo)*windo) +1
    for k = 0:(windo-1)
        if Healthy_1(i+k) <= 50 | Healthy_1(i+k) >= 200
            percent = percent +1;
            if i+k == index +1
                count = count +1;
            else 
                if count >= 8 
                    Healthy_1_corr(i+k-count-1:i+k-1) = NaN([1 count+1]);
                    consec_8=1; %add to save if any time inside a window there were 8 consect;
                    count = 0;                    
                else
                    count=0;
                end
            end
            index = k + i;
        end
    end
    if (percent >= floor(5*windo/100) || consec_8==1)
        Healthy_1_corr(i:i+windo-1) = NaN([1 windo]);
        percent = 0;
        consec_8=0;
    end
    i = i+windo;
end

% figure
% plot(Healthy_1_corr)

%Remaing "for analysis" intervals: correct missing data + artifact 
count=0;
consec=0;
j=0;
out_ind_H1=zeros(1,1);
for i=1:length(Healthy_1)-3
    if(~isnan(Healthy_1_corr(i:i+3)))
        for j=i:i+2
            if 1/(Healthy_1_corr(j+1)/60) < 2*(1/(Healthy_1_corr(j)/60)) - 0.3 && 1/(Healthy_1_corr(j+1)/60) > 0.57*1/(Healthy_1_corr(j)/60) + 0.43 * 0.3
                consec=consec+1;
            end
        end
        if(consec==3)
            consec=0;
        else
            if(Healthy_1_corr(i)==0)
                count=count+1;
                out_ind_H1(count)=i;
            end
            count=count+1;
            out_ind_H1(count)=i+1;
        end
    end
    consec=0;
    
end

% figure
% plot(Healthy_1_corr), hold on, stem(out_ind_H1, Healthy_1_corr(out_ind_H1), '-.or')

%Put all artifacts and missin data to NaN so it's ignore on the next step
%to correct the values throught moving average 
for i=1:length(out_ind_H1)
    Healthy_1_corr(out_ind_H1(i))= NaN ;
end
% figure
% plot(Healthy_1_corr);

%Interpolate with movavg
count=0;
for i=1:length(out_ind_H1)

        movavg=movmean(Healthy_1_corr, [5,5], 'omitnan');
        Healthy_1_corr(out_ind_H1(i))=movavg(out_ind_H1(i));
end

% figure
% plot(Healthy_1_corr), hold on, stem(out_ind_H1, Healthy_1_corr(out_ind_H1), '-.or')


%Creating Surbinterval of size = windo
Healthy_1_inter = zeros(windo,floor(length(Healthy_1)/windo));
t_Healthy_1_inter = zeros(windo,floor(length(Healthy_1)/windo));
for i = 1:floor(length(Healthy_1)/windo)
    Healthy_1_inter(:,i) = Healthy_1_corr((i-1)*windo+1:i*windo);
    t_Healthy_1_inter(:,i) = t_Healthy_1((i-1)*windo+1:i*windo);
end
% figure;
% for i = 1:floor(length(Healthy_1)/windo)
%     plot(t_Healthy_1_inter(:,1), Healthy_1_inter(:,i));hold on;
% end

%HEALTHY 2 
%Cutting intervals with too missing data so moving average is not perform
%in those
Healthy_2_corr = Healthy_2;
count = 0;
percent = 0;
index = 1;
i = 1;
consec_8=0;
while i <= length(Healthy_2) - (windo + length(Healthy_2)-floor(length(Healthy_2)/windo)*windo) +1
    for k = 0:(windo-1)
        if Healthy_2(i+k) <= 50 | Healthy_2(i+k) >= 200
            percent = percent +1;
            if i+k == index +1
                count = count +1;
            else 
                if count >= 8 
                    Healthy_2_corr(i+k-count-1:i+k-1) = NaN([1 count+1]);
                    consec_8=1; %add to save if any time inside a window there were 5 consect;
                    count = 0;                    
                else
                    count=0;
                end
            end
            index = k + i;
        end
    end
    if (percent >= floor(5*windo/100) || consec_8==1)
        Healthy_2_corr(i:i+windo-1) = NaN([1 windo]);
        percent = 0;
        consec_8=0;
    end
    i = i+windo;
end
% figure
% plot(Healthy_2_corr)

%Remaing "for analysis" intervals: correct missing data + artifact 
count=0;
consec=0;
j=0;
out_ind_H2=zeros(1,1);
for i=1:length(Healthy_2)-3
    if(~isnan(Healthy_2_corr(i:i+3)))
        for j=i:i+2
            if((1/(Healthy_2_corr(j+1)/60)) < (2*(1/(Healthy_2_corr(j)/60))) - 0.3 && (1/(Healthy_2_corr(j+1)/60))>0.57*(1/(Healthy_2_corr(j)/60)) +0.43*0.3)
                consec=consec+1;
            end
        end
        if(consec==3)
            consec=0;
        else
            if(Healthy_2_corr(i)==0)
                count=count+1;
                out_ind_H2(count)=i;
            end
            count=count+1;
            out_ind_H2(count)=i+1;
        end
    end
    consec=0;
    
end

% figure
% plot(Healthy_2_corr), hold on, stem(out_ind_H2, Healthy_2_corr(out_ind_H2), '-.or')

%Put all artifacts and missin data to NaN so it's ignore on the next step
%to correct the values throught moving average 
for i=1:length(out_ind_H2)
    Healthy_2_corr(out_ind_H2(i))=NaN;
end
% figure
% plot(Healthy_2_corr)
% 
%Interpolate with movavg
count=0;
for i=1:length(out_ind_H2)

        movavg=movmean(Healthy_2_corr, [5,5], 'omitnan');
        Healthy_2_corr(out_ind_H2(i))=movavg(out_ind_H2(i));
end

% figure
% plot(Healthy_2_corr), hold on, stem(out_ind_H2, Healthy_2_corr(out_ind_H2), '-.or')


%Creating Surbinterval of size = windo
Healthy_2_inter = zeros(windo,floor(length(Healthy_2)/windo));
t_Healthy_2_inter = zeros(windo,floor(length(Healthy_2)/windo));
for i = 1:floor(length(Healthy_2)/windo)
    Healthy_2_inter(:,i) = Healthy_2_corr((i-1)*windo+1:i*windo);
    t_Healthy_2_inter(:,i) = t_Healthy_2((i-1)*windo+1:i*windo);
end
% figure;
% for i = 1:floor(length(Healthy_2)/windo)
%     plot(t_Healthy_2_inter(:,1), Healthy_2_inter(:,i));hold on;
% end


%HEALTHY 3
%Cutting intervals with too missing data so moving average is not perform
%in those
Healthy_3_corr = Healthy_3;
count = 0;
percent = 0;
index = 1;
i = 1;
consec_8=0;
while i <= length(Healthy_3) - (windo + length(Healthy_3)-floor(length(Healthy_3)/windo)*windo) +1
    
    for k = 0:(windo-1)
        if Healthy_3(i+k) <= 50 | Healthy_3(i+k) >= 200
            percent = percent +1;
            if i+k == index +1
                count = count +1;
            else 
                if count >= 8 
                    %Healthy_3(i+k-count-1:i+k-1) = NaN([1 count+1]);
                    consec_8=1; %add to save if any time inside a window there were 5 consect;
                    count = 0;                    
                else
                    count=0;
                end
            end
            index = k + i;
        end
    end
    if (percent >= floor(5*windo/100) || consec_8==1)
        Healthy_3_corr(i:i+windo-1) = NaN([1 windo]);
        percent = 0;
        consec_8=0;
    end
    i = i+windo;
end
% figure
% plot(Healthy_3)

%Remaing "for analysis" intervals: correct missing data + artifact 
count=0;
consec=0;
j=0;
out_ind_H3=zeros(1,1);
for i=1:length(Healthy_3)-3
    if(~isnan(Healthy_3_corr(i:i+3)))
        for j=i:i+2
            if((1/(Healthy_3_corr(j+1)/60)) < (2*(1/(Healthy_3_corr(j)/60))) - 0.3 && (1/(Healthy_3_corr(j+1)/60))>0.57*(1/(Healthy_3_corr(j)/60)) +0.43*0.3)
                consec=consec+1;
            end
        end
        if(consec==3)
            consec=0;
        else
            if(Healthy_3_corr(i)==0)
                count=count+1;
                out_ind_H3(count)=i;
            end
            count=count+1;
            out_ind_H3(count)=i+1;
        end
    end
    consec=0;
    
end

% figure
% plot(Healthy_3), hold on, stem(out_ind_H3, Healthy_3(out_ind_H3), '-.or')

%Put all artifacts and missin data to NaN so it's ignore on the next step
%to correct the values throught moving average 
for i=1:length(out_ind_H3)
    Healthy_3_corr(out_ind_H3(i))=NaN;
end
% figure
% plot(Healthy_3)

%Interpolate with movavg
count=0;
for i=1:length(out_ind_H3)

        movavg=movmean(Healthy_3_corr, [5,5], 'omitnan');
        Healthy_3_corr(out_ind_H3(i))=movavg(out_ind_H3(i));
end

% figure
% plot(Healthy_3), hold on, stem(out_ind_H3, Healthy_3(out_ind_H3), '-.or')

%Creating Surbinterval of size = windo
Healthy_3_inter = zeros(windo,floor(length(Healthy_3)/windo));
t_Healthy_3_inter = zeros(windo,floor(length(Healthy_3)/windo));
for i = 1:floor(length(Healthy_3)/windo)
    Healthy_3_inter(:,i) = Healthy_3_corr((i-1)*windo+1:i*windo);
    t_Healthy_3_inter(:,i) = t_Healthy_3((i-1)*windo+1:i*windo);
end
% figure;
% for i = 1:floor(length(Healthy_3)/windo)
%     plot(t_Healthy_3_inter(:,1), Healthy_3_inter(:,i));hold on;
% end

%HEALTHY 4
%Cutting intervals with too missing data so moving average is not perform
%in those
Healthy_4_corr = Healthy_4;
count = 0;
percent = 0;
index = 1;
i = 1;
consec_8=0;
while i <= length(Healthy_4) - (windo + length(Healthy_4)-floor(length(Healthy_4)/windo)*windo) +1
    
    for k = 0:(windo-1)
        if Healthy_4_corr(i+k) <= 50 | Healthy_4_corr(i+k) >= 200
            percent = percent +1;
            if i+k == index +1
                count = count +1;
            else 
                if count >= 8 
                    %Healthy_3(i+k-count-1:i+k-1) = NaN([1 count+1]);
                    consec_8=1; %add to save if any time inside a window there were 5 consect;
                    count = 0;                    
                else
                    count=0;
                end
            end
            index = k + i;
        end
    end
    if (percent >= floor(5*windo/100) || consec_8==1)
        Healthy_4_corr(i:i+windo-1) = NaN([1 windo]);
        percent = 0;
        consec_8=0;
    end
    i = i+windo;
end
% figure
% plot(Healthy_4)

%Remaing "for analysis" intervals: correct missing data + artifact 
count=0;
consec=0;
j=0;
out_ind_H4=zeros(1,1);
for i=1:length(Healthy_4)-3
    if(~isnan(Healthy_4_corr(i:i+3)))
        for j=i:i+2
            if((1/(Healthy_4_corr(j+1)/60)) < (2*(1/(Healthy_4_corr(j)/60))) - 0.3 && (1/(Healthy_4_corr(j+1)/60))>0.57*(1/(Healthy_4_corr(j)/60)) +0.43*0.3)
                consec=consec+1;
            end
        end
        if(consec==3)
            consec=0;
        else
            if(Healthy_4_corr(i)==0)
                count=count+1;
                out_ind_H4(count)=i;
            end
            count=count+1;
            out_ind_H4(count)=i+1;
        end
    end
    consec=0;
    
end

% figure
% plot(Healthy_4), hold on, stem(out_ind_H4, Healthy_4(out_ind_H4), '-.or')

%Put all artifacts and missin data to NaN so it's ignore on the next step
%to correct the values throught moving average 
for i=1:length(out_ind_H4)
    Healthy_4_corr(out_ind_H4(i))=NaN;
end
% figure
% plot(Healthy_4)

%Interpolate with movavg
count=0;
for i=1:length(out_ind_H4)
        movavg=movmean(Healthy_4_corr, [5,5], 'omitnan');
        Healthy_4_corr(out_ind_H4(i))=movavg(out_ind_H4(i));
end

% figure
% plot(Healthy_4), hold on, stem(out_ind_H4, Healthy_4(out_ind_H4), '-.or')

%Creating Surbinterval of size = windo
Healthy_4_inter = zeros(windo,floor(length(Healthy_4)/windo));
t_Healthy_4_inter = zeros(windo,floor(length(Healthy_4)/windo));
for i = 1:floor(length(Healthy_4)/windo)
    Healthy_4_inter(:,i) = Healthy_4_corr((i-1)*windo+1:i*windo);
    t_Healthy_4_inter(:,i) = t_Healthy_4((i-1)*windo+1:i*windo);
end
% figure;
% for i = 1:floor(length(Healthy_4)/windo)
%     plot(t_Healthy_4_inter(:,1), Healthy_4_inter(:,i));hold on;
% end


%HEALTHY 5
%Cutting intervals with too missing data so moving average is not perform
%in those
Healthy_5_corr = Healthy_5;
count = 0;
percent = 0;
index = 1;
i = 1;
consec_8=0;
while i <= length(Healthy_5) - (windo + length(Healthy_5)-floor(length(Healthy_5)/windo)*windo) +1
   
    for k = 0:(windo-1)
        if Healthy_5(i+k) <= 50 | Healthy_5(i+k) >= 200
            percent = percent +1;
            if i+k == index +1
                count = count +1;
            else 
                if count >= 8 
                    %Healthy_3(i+k-count-1:i+k-1) = NaN([1 count+1]);
                    consec_8=1; %add to save if any time inside a window there were 5 consect;
                    count = 0;                    
                else
                    count=0;
                end
            end
            index = k + i;
        end
    end
    if (percent >= floor(5*windo/100) || consec_8==1)
        Healthy_5_corr(i:i+windo-1) = NaN([1 windo]);
        percent = 0;
        consec_8=0;
    end
    i = i+windo;
end
% figure
% plot(Healthy_5)

%Remaing "for analysis" intervals: correct missing data + artifact 
count=0;
consec=0;
j=0;
out_ind_H5=zeros(1,1);
for i=1:length(Healthy_5)-3
    if(~isnan(Healthy_5_corr(i:i+3)))
        for j=i:i+2
            if((1/(Healthy_5_corr(j+1)/60)) < (2*(1/(Healthy_5_corr(j)/60))) - 0.3 && (1/(Healthy_5_corr(j+1)/60))>0.57*(1/(Healthy_5_corr(j)/60)) +0.43*0.3)
                consec=consec+1;
            end
        end
        if(consec==3)
            consec=0;
        else
            if(Healthy_5_corr(i)==0)
                count=count+1;
                out_ind_H5(count)=i;
            end
            count=count+1;
            out_ind_H5(count)=i+1;
        end
    end
    consec=0;
    
end

% figure
% plot(Healthy_5), hold on, stem(out_ind_H5, Healthy_5(out_ind_H5), '-.or')

%Put all artifacts and missin data to NaN so it's ignore on the next step
%to correct the values throught moving average 
for i=1:length(out_ind_H5)
    Healthy_5_corr(out_ind_H5(i))=NaN;
end
% figure
% plot(Healthy_5)

%Interpolate with movavg
count=0;
for i=1:length(out_ind_H5)
        movavg=movmean(Healthy_5_corr, [5,5], 'omitnan');
        Healthy_5_corr(out_ind_H5(i))=movavg(out_ind_H5(i));
end
% 
% figure
% plot(Healthy_5), hold on, stem(out_ind_H5, Healthy_5(out_ind_H5), '-.or')

%Creating Surbinterval of size = windo
Healthy_5_inter = zeros(windo,floor(length(Healthy_5)/windo));
t_Healthy_5_inter = zeros(windo,floor(length(Healthy_5)/windo));
for i = 1:floor(length(Healthy_5)/windo)
    Healthy_5_inter(:,i) = Healthy_5_corr((i-1)*windo+1:i*windo);
    t_Healthy_5_inter(:,i) = t_Healthy_5((i-1)*windo+1:i*windo);
end
% figure;
% for i = 1:floor(length(Healthy_5)/windo)
%     plot(t_Healthy_5_inter(:,1), Healthy_5_inter(:,i));hold on;
% end

%% IUGR Pre-processing 
%IUGR 1
%Cutting intervals with too missing data so moving average is not perform
%in those
IUGR_1_corr = IUGR_1;
count = 0;
percent = 0;
index = 1;
i = 1;
consec_8=0;
while i <= length(IUGR_1) - (windo + length(IUGR_1)-floor(length(IUGR_1)/windo)*windo) +1
    
    for k = 0:(windo-1)
        if IUGR_1(i+k) <= 50 | IUGR_1(i+k) >= 200
            percent = percent +1;
            if i+k == index +1
                count = count +1;
            else 
                if count >= 8 
                    %Healthy_3(i+k-count-1:i+k-1) = NaN([1 count+1]);
                    consec_8=1; %add to save if any time inside a window there were 5 consect;
                    count = 0;                    
                else
                    count=0;
                end
            end
            index = k + i;
        end
    end
    if (percent >= floor(5*windo/100) || consec_8==1)
        IUGR_1_corr(i:i+windo-1) = NaN([1 windo]);
        percent = 0;
        consec_8=0;
    end
    i = i+windo;
end
% figure
% plot(IUGR_1_corr)

%Remaing "for analysis" intervals: correct missing data + artifact 
count=0;
consec=0;
j=0;
out_ind_I1=zeros(1,1);
for i=1:length(IUGR_1)-3
    if(~isnan(IUGR_1_corr(i:i+3)))
        for j=i:i+2
            if((1/(IUGR_1_corr(j+1)/60)) < (2*(1/(IUGR_1_corr(j)/60))) - 0.3 && (1/(IUGR_1_corr(j+1)/60))>0.57*(1/(IUGR_1_corr(j)/60)) +0.43*0.3)
                consec=consec+1;
            end
        end
        if(consec==3)
            consec=0;
        else
            if(IUGR_1_corr(i)==0)
                count=count+1;
                out_ind_I1(count)=i;
            end
            count=count+1;
            out_ind_I1(count)=i+1;
        end
    end
    consec=0;
    
end

if (out_ind_I1 >0)
%     figure
%     plot(IUGR_1), hold on, stem(out_ind_I1, IUGR_1(out_ind_I1), '-.or')
    
    %Put all artifacts and missin data to NaN so it's ignore on the next step
    %to correct the values throught moving average
    for i=1:length(out_ind_I1)
        IUGR_1_corr(out_ind_I1(i))=NaN;
    end
    %figure
    %plot(IUGR_1)
    
    %Interpolate with movavg
    count=0;
    for i=1:length(out_ind_I1)
        movavg=movmean(IUGR_1_corr, [5,5], 'omitnan');
        IUGR_1_corr(out_ind_I1(i))=movavg(out_ind_I1(i));
    end
    
    
    % figure
    % plot(IUGR_1), hold on, stem(out_ind_I1, IUGR_1(out_ind_I1), '-.or')
    
else
    fprintf('No artifacts or missing data found');
end

%Creating Surbinterval of size = windo
IUGR_1_inter = zeros(windo,floor(length(IUGR_1)/windo));
t_IUGR_1_inter = zeros(windo,floor(length(IUGR_1)/windo));
for i = 1:floor(length(IUGR_1)/windo)
    IUGR_1_inter(:,i) = IUGR_1_corr((i-1)*windo+1:i*windo);
    t_IUGR_1_inter(:,i) = t_IUGR_1((i-1)*windo+1:i*windo);
end
% figure;
% for i = 1:floor(length(IUGR_1)/windo)
%     plot(t_IUGR_1_inter(:,1), IUGR_1_inter(:,i));hold on;
% end


%IUGR 2
IUGR_2_corr = IUGR_2;
%Cutting intervals with too missing data so moving average is not perform
%in those
count = 0;
percent = 0;
index = 1;
i = 1;
consec_8=0;
while i <= length(IUGR_2) - (windo + length(IUGR_2)-floor(length(IUGR_2)/windo)*windo) +1
    
    for k = 0:(windo-1)
        if IUGR_2(i+k) <= 50 | IUGR_2(i+k) >= 200
            percent = percent +1;
            if i+k == index +1
                count = count +1;
            else 
                if count >= 8 
                    %Healthy_3(i+k-count-1:i+k-1) = NaN([1 count+1]);
                    consec_8=1; %add to save if any time inside a window there were 5 consect;
                    count = 0;                    
                else
                    count=0;
                end
            end
            index = k + i;
        end
    end
    if (percent >= floor(5*windo/100) || consec_8==1)
        IUGR_2_corr(i:i+windo-1) = NaN([1 windo]);
        percent = 0;
        consec_8=0;
    end
    i = i+windo;
end
% figure
% plot(IUGR_2)

%Remaing "for analysis" intervals: correct missing data + artifact 
count=0;
consec=0;
j=0;
out_ind_I2=zeros(1,1);
for i=1:length(IUGR_2)-3
    if(~isnan(IUGR_2_corr(i:i+3)))
        for j=i:i+2
            if((1/(IUGR_2_corr(j+1)/60)) < (2*(1/(IUGR_2_corr(j)/60))) - 0.3 && (1/(IUGR_2_corr(j+1)/60))>0.57*(1/(IUGR_2_corr(j)/60)) +0.43*0.3)
                consec=consec+1;
            end
        end
        if(consec==3)
            consec=0;
        else
            if(IUGR_2_corr(i)==0)
                count=count+1;
                out_ind_I2(count)=i;
            end
            count=count+1;
            out_ind_I2(count)=i+1;
        end
    end
    consec=0;
    
end

% if (out_ind_I2 >0)
%     figure
%     plot(IUGR_2), hold on, stem(out_ind_I2, IUGR_2(out_ind_I2), '-.or')
% else
%     fprintf('No artifacts or missing data found');
% end

%Put all artifacts and missin data to NaN so it's ignore on the next step
%to correct the values throught moving average 
for i=1:length(out_ind_I2)
    IUGR_2_corr(out_ind_I2(i))=NaN;
end
% figure
% plot(IUGR_2)

%Interpolate with movavg
count=0;
for i=1:length(out_ind_I2)

        movavg=movmean(IUGR_2_corr, [5,5], 'omitnan');
        IUGR_2_corr(out_ind_I2(i))=movavg(out_ind_I2(i));
end

% figure
% plot(IUGR_2), hold on, stem(out_ind_I2, IUGR_2(out_ind_I2), '-.or')

%Creating Surbinterval of size = windo
IUGR_2_inter = zeros(windo,floor(length(IUGR_2)/windo));
t_IUGR_2_inter = zeros(windo,floor(length(IUGR_2)/windo));
for i = 1:floor(length(IUGR_2)/windo)
    IUGR_2_inter(:,i) = IUGR_2_corr((i-1)*windo+1:i*windo);
    t_IUGR_2_inter(:,i) = t_IUGR_2((i-1)*windo+1:i*windo);
end
% figure;
% for i = 1:floor(length(IUGR_2)/windo)
%     plot(t_IUGR_2_inter(:,1), IUGR_2_inter(:,i));hold on;
% end


%IUGR 3
%Cutting intervals with too missing data so moving average is not perform
%in those
IUGR_3_corr = IUGR_3;
count = 0;
percent = 0;
index = 1;
i = 1;
consec_8=0;
while i <= length(IUGR_3) - (windo + length(IUGR_3)-floor(length(IUGR_3)/windo)*windo) +1
    
    for k = 0:(windo-1)
        if IUGR_3(i+k) <= 50 | IUGR_3(i+k) >= 200
            percent = percent +1;
            if i+k == index +1
                count = count +1;
            else 
                if count >= 8 
                    %Healthy_3(i+k-count-1:i+k-1) = NaN([1 count+1]);
                    consec_8=1; %add to save if any time inside a window there were 5 consect;
                    count = 0;                    
                else
                    count=0;
                end
            end
            index = k + i;
        end
    end
    if (percent >= floor(5*windo/100) || consec_8==1)
        IUGR_3_corr(i:i+windo-1) = NaN([1 windo]);
        percent = 0;
        consec_8=0;
    end
    i = i+windo;
end
% figure
% plot(IUGR_3)

%Remaing "for analysis" intervals: correct missing data + artifact 
count=0;
consec=0;
j=0;
out_ind_I3=zeros(1,1);
for i=1:length(IUGR_3)-3
    if(~isnan(IUGR_3_corr(i:i+3)))
        for j=i:i+2
            if((1/(IUGR_3_corr(j+1)/60)) < (2*(1/(IUGR_3_corr(j)/60))) - 0.3 && (1/(IUGR_3_corr(j+1)/60))>0.57*(1/(IUGR_3_corr(j)/60)) +0.43*0.3)
                consec=consec+1;
            end
        end
        if(consec==3)
            consec=0;
        else
            if(IUGR_3_corr(i)==0)
                count=count+1;
                out_ind_I3(count)=i;
            end
            count=count+1;
            out_ind_I3(count)=i+1;
        end
    end
    consec=0;
    
end

% if (out_ind_I3 >0)
%     figure
%     plot(IUGR_3), hold on, stem(out_ind_I3, IUGR_3(out_ind_I3), '-.or')
% else
%     fprintf('No artifacts or missing data found');
% end

%Put all artifacts and missin data to NaN so it's ignore on the next step
%to correct the values throught moving average 
for i=1:length(out_ind_I3)
    IUGR_3_corr(out_ind_I3(i))=NaN;
end
% figure
% plot(IUGR_3)

%Interpolate with movavg
count=0;
for i=1:length(out_ind_I3)

        movavg=movmean(IUGR_3_corr, [5,5], 'omitnan');
        IUGR_3_corr(out_ind_I3(i))=movavg(out_ind_I3(i));
end

% figure
% plot(IUGR_3), hold on, stem(out_ind_I3, IUGR_3(out_ind_I3), '-.or')

%Creating Surbinterval of size = windo
IUGR_3_inter = zeros(windo,floor(length(IUGR_3)/windo));
t_IUGR_3_inter = zeros(windo,floor(length(IUGR_3)/windo));
for i = 1:floor(length(IUGR_3)/windo)
    IUGR_3_inter(:,i) = IUGR_3_corr((i-1)*windo+1:(i)*windo);
    t_IUGR_3_inter(:,i) = t_IUGR_3((i-1)*windo+1:i*windo);
end
% figure;
% for i = 1:floor(length(IUGR_3)/windo)
%     plot(t_IUGR_3_inter(:,1), IUGR_3_inter(:,i));hold on;
% end

%IUGR 4
%Cutting intervals with too missing data so moving average is not perform
%in those
IUGR_4_corr = IUGR_4;
count = 0;
percent = 0;
index = 1;
i = 1;
consec_8=0;
while i <= length(IUGR_4) - (windo + length(IUGR_4)-floor(length(IUGR_4)/windo)*windo) +1
    
    for k = 0:(windo-1)
        if IUGR_4(i+k) <= 50 | IUGR_4(i+k) >= 200
            percent = percent +1;
            if i+k == index +1
                count = count +1;
            else 
                if count >= 8 
                    %Healthy_3(i+k-count-1:i+k-1) = NaN([1 count+1]);
                    consec_8=1; %add to save if any time inside a window there were 5 consect;
                    count = 0;                    
                else
                    count=0;
                end
            end
            index = k + i;
        end
    end
    if (percent >= floor(5*windo/100) || consec_8==1)
        IUGR_4_corr(i:i+windo-1) = NaN([1 windo]);
        percent = 0;
        consec_8=0;
    end
    i = i+windo;
end
% figure
% plot(IUGR_4)

%Remaing "for analysis" intervals: correct missing data + artifact 
count=0;
consec=0;
j=0;
out_ind_I4=zeros(1,1);
for i=1:length(IUGR_4)-3
    if(~isnan(IUGR_4_corr(i:i+3)))
        for j=i:i+2
            if((1/(IUGR_4_corr(j+1)/60)) < (2*(1/(IUGR_4_corr(j)/60))) - 0.3 && (1/(IUGR_4_corr(j+1)/60))>0.57*(1/(IUGR_4_corr(j)/60)) +0.43*0.3)
                consec=consec+1;
            end
        end
        if(consec==3)
            consec=0;
        else
            if(IUGR_4_corr(i)==0)
                count=count+1;
                out_ind_I4(count)=i;
            end
            count=count+1;
            out_ind_I4(count)=i+1;
        end
    end
    consec=0;
    
end

% if (out_ind_I4 >0)
%     figure
%     plot(IUGR_4), hold on, stem(out_ind_I4, IUGR_4(out_ind_I4), '-.or')
% else
%     fprintf('No artifacts or missing data found');
% end

%Put all artifacts and missin data to NaN so it's ignore on the next step
%to correct the values throught moving average 
for i=1:length(out_ind_I4)
    IUGR_4_corr(out_ind_I4(i))=NaN;
end
% figure
% plot(IUGR_4)

%Interpolate with movavg
count=0;
for i=1:length(out_ind_I4)

        movavg=movmean(IUGR_4_corr, [5,5], 'omitnan');
        IUGR_4_corr(out_ind_I4(i))=movavg(out_ind_I4(i));
end

% figure
% plot(IUGR_4), hold on, stem(out_ind_I4, IUGR_4(out_ind_I4), '-.or')

%Creating Surbinterval of size = windo
IUGR_4_inter = zeros(windo,floor(length(IUGR_4)/windo));
t_IUGR_4_inter = zeros(windo,floor(length(IUGR_4)/windo));
for i = 1:floor(length(IUGR_4)/windo)
    IUGR_4_inter(:,i) = IUGR_4_corr((i-1)*windo+1:(i)*windo);
    t_IUGR_4_inter(:,i) = t_IUGR_4((i-1)*windo+1:i*windo);
end
% figure;
% for i = 1:floor(length(IUGR_4)/windo)
%     plot(t_IUGR_4_inter(:,1), IUGR_4_inter(:,i));hold on;
% end

%IUGR 5
%Cutting intervals with too missing data so moving average is not perform
%in those
IUGR_5_corr= IUGR_5;
count = 0;
percent = 0;
index = 1;
i = 1;
consec_8=0;
while i <= length(IUGR_5) - (windo + length(IUGR_5)-floor(length(IUGR_5)/windo)*windo) +1
    
    for k = 0:(windo-1)
        if IUGR_5(i+k) <= 50 | IUGR_5(i+k) >= 200
            percent = percent +1;
            if i+k == index +1
                count = count +1;
            else 
                if count >= 8 
                    %Healthy_3(i+k-count-1:i+k-1) = NaN([1 count+1]);
                    consec_8=1; %add to save if any time inside a window there were 5 consect;
                    count = 0;                    
                else
                    count=0;
                end
            end
            index = k + i;
        end
    end
    if (percent >= floor(5*windo/100) || consec_8==1)
        IUGR_5_corr(i:i+windo-1) = NaN([1 windo]);
        percent = 0;
        consec_8=0;
    end
    i = i+windo;
end
% figure
% plot(IUGR_5)

%Remaing "for analysis" intervals: correct missing data + artifact 
count=0;
consec=0;
j=0;
out_ind_I5=zeros(1,1);
for i=1:length(IUGR_5)-3
    if(~isnan(IUGR_5_corr(i:i+3)))
        for j=i:i+2
            if((1/(IUGR_5_corr(j+1)/60)) < (2*(1/(IUGR_5_corr(j)/60))) - 0.3 && (1/(IUGR_5_corr(j+1)/60))>0.57*(1/(IUGR_5_corr(j)/60)) +0.43*0.3)
                consec=consec+1;
            end
        end
        if(consec==3)
            consec=0;
        else
            if(IUGR_5_corr(i)==0)
                count=count+1;
                out_ind_I5(count)=i;
            end
            count=count+1;
            out_ind_I5(count)=i+1;
        end
    end
    consec=0;
    
end

% if (out_ind_I5 >0)
%     figure
%     plot(IUGR_5), hold on, stem(out_ind_I5, IUGR_5(out_ind_I5), '-.or')
% else
%     fprintf('No artifacts or missing data found');
% end

%Put all artifacts and missin data to NaN so it's ignore on the next step
%to correct the values throught moving average 
for i=1:length(out_ind_I5)
    IUGR_5_corr(out_ind_I5(i))=NaN;
end
% figure
% plot(IUGR_5)

%Interpolate with movavg
count=0;
for i=1:length(out_ind_I5)
        movavg=movmean(IUGR_5_corr, [5,5], 'omitnan');
        IUGR_5_corr(out_ind_I5(i))=movavg(out_ind_I5(i));
end

% figure
% plot(IUGR_5), hold on, stem(out_ind_I5, IUGR_5(out_ind_I5), '-.or')

%Creating Surbinterval of size = windo
IUGR_5_inter = zeros(windo,floor(length(IUGR_5)/windo));
t_IUGR_5_inter = zeros(windo,floor(length(IUGR_5)/windo));
for i = 1:floor(length(IUGR_5)/windo)
    IUGR_5_inter(:,i) = IUGR_5_corr((i-1)*windo+1:(i)*windo);
    t_IUGR_5_inter(:,i) = t_IUGR_5((i-1)*windo+1:i*windo);
end
% figure;
% for i = 1:floor(length(IUGR_5)/windo)
%     plot(t_IUGR_5_inter(:,1), IUGR_5_inter(:,i));hold on;
% end

%% SIGNAL LOSS : 
ind_1 = find(isnan(Healthy_1_corr));
loss_Healthy_1 = length(ind_1)/length(Healthy_1_corr);
sig_loss_H1=['Sig\_Loss=', num2str(loss_Healthy_1)];

ind_2 = find(isnan(Healthy_2_corr));
loss_Healthy_2 = length(ind_2)/length(Healthy_2_corr);
sig_loss_H2=['Sig\_Loss=', num2str(loss_Healthy_2)];

ind_3 = find(isnan(Healthy_3_corr));
loss_Healthy_3 = length(ind_3)/length(Healthy_3_corr);
sig_loss_H3=['Sig\_Loss=', num2str(loss_Healthy_3)];

ind_4 = find(isnan(Healthy_4_corr));
loss_Healthy_4 = length(ind_4)/length(Healthy_4_corr);
sig_loss_H4=['Sig\_Loss=', num2str(loss_Healthy_4)];

ind_5 = find(isnan(Healthy_5_corr));
loss_Healthy_5 = length(ind_5)/length(Healthy_5_corr);
sig_loss_H5=['Sig\_Loss=', num2str(loss_Healthy_5)];

indI_1 = find(isnan(IUGR_1_corr));
loss_IUGR_1 = length(indI_1)/length(IUGR_1_corr);
sig_loss_I1=['Sig\_Loss=', num2str(loss_IUGR_1)];

indI_2 = find(isnan(IUGR_2_corr));
loss_IUGR_2 = length(indI_2)/length(IUGR_2_corr);
sig_loss_I2=['Sig\_Loss=', num2str(loss_IUGR_2)];


indI_3 = find(isnan(IUGR_3_corr));
loss_IUGR_3 = length(indI_3)/length(IUGR_3_corr);
sig_loss_I3=['Sig\_Loss=', num2str(loss_IUGR_3)];


indI_4 = find(isnan(IUGR_4_corr));
loss_IUGR_4 = length(indI_4)/length(IUGR_4_corr);
sig_loss_I4=['Sig\_Loss=', num2str(loss_IUGR_4)];


indI_5 = find(isnan(IUGR_5_corr));
loss_IUGR_5 = length(indI_5)/length(IUGR_5_corr);
sig_loss_I5=['Sig\_Loss=', num2str(loss_IUGR_5)];


%Plot preprocessed signal with corresponding signal loss ratio
figure
subplot(3,2,1), plot(t_Healthy_1, Healthy_1_corr), title('Healthy\_1 corr'), xlabel('Time in sec');
ylabel('FHR bpm'), text(1400, 180, sig_loss_H1,'HorizontalAlignment','left');
subplot(3,2,2), plot(t_Healthy_2, Healthy_2_corr), title('Healthy\_2 corr'), xlabel('Time in sec');
ylabel('FHR bpm'), text(1000, 190, sig_loss_H2,'HorizontalAlignment','left');;
subplot(3,2,3), plot(t_Healthy_3, Healthy_3_corr), title('Healthy\_3 corr'), xlabel('Time in sec');
ylabel('FHR bpm'), text(1400, 190, sig_loss_H3,'HorizontalAlignment','left');
subplot(3,2,4), plot(t_Healthy_4, Healthy_4_corr), title('Healthy\_4 corr'), xlabel('Time in sec');
ylabel('FHR bpm'), text(1000, 180, sig_loss_H4,'HorizontalAlignment','left');
subplot(3,2,5), plot(t_Healthy_5, Healthy_5_corr), title('Healthy\_5 corr'), xlabel('Time in sec');
ylabel('FHR bpm'), text(1400, 110, sig_loss_H5,'HorizontalAlignment','left');

figure
subplot(2,3,1), plot(t_IUGR_1,IUGR_1_corr), xlabel('Time in sec');
ylabel('FHR bpm'), text(1300, 125, sig_loss_I1,'HorizontalAlignment','left');
title('Preprocessed IUGR\_1');
subplot(2,3,2), plot(t_IUGR_2,IUGR_2_corr), xlabel('Time in sec');
ylabel('FHR bpm'), text(1100, 145, sig_loss_I2,'HorizontalAlignment','left');
title('Preprocessed IUGR\_2');
subplot(2,3,3), plot(t_IUGR_3,IUGR_3_corr), xlabel('Time in sec');
ylabel('FHR bpm'), text(1300, 155, sig_loss_I3,'HorizontalAlignment','left');
title('Preprocessed IUGR\_3');
subplot(2,3,4),plot(t_IUGR_4,IUGR_4_corr), xlabel('Time in sec');
ylabel('FHR bpm'), text(1200, 135, sig_loss_I4,'HorizontalAlignment','left');
title('Preprocessed IUGR\_4');
subplot(2,3,5),plot(t_IUGR_5,IUGR_5_corr),xlabel('Time in sec');
ylabel('FHR bpm'), text(1500, 160, sig_loss_I5,'HorizontalAlignment','left');
title('Preprocessed IUGR\_5');

%% Mean Variance Stationarity 
%Healthy_1
wrapping=windo; %defined in beat
num_seg=floor(length(Healthy_1_corr)/wrapping) ;
Healthy_1_cut=Healthy_1_corr(1:num_seg*wrapping);

Healthy_1_reshape=reshape(Healthy_1_cut, [wrapping, length(Healthy_1_cut)/wrapping]);
Healthy_1_mean=mean(Healthy_1_reshape);
Healthy_1_var=var(Healthy_1_reshape);
t_Healthy_1_wrapping=0+wrapping:wrapping:floor(length(Healthy_1)/wrapping)*wrapping;

%Healthy_2
wrapping=windo; %defined in beat
num_seg=floor(length(Healthy_2_corr)/wrapping) ;
Healthy_2_cut=Healthy_2_corr(1:num_seg*wrapping);
Healthy_2_reshape=reshape(Healthy_2_cut, [wrapping, length(Healthy_2_cut)/wrapping]);
Healthy_2_mean=mean(Healthy_2_reshape);
Healthy_2_var=var(Healthy_2_reshape);
t_Healthy_2_wrapping=0+wrapping:wrapping:floor(length(Healthy_2)/wrapping)*wrapping;

%Healthy_3
wrapping=windo; %defined in beat
num_seg=floor(length(Healthy_3_corr)/wrapping) ;
Healthy_3_cut=Healthy_3_corr(1:num_seg*wrapping);
Healthy_3_reshape=reshape(Healthy_3_cut, [wrapping, length(Healthy_3_cut)/wrapping]);
Healthy_3_mean=mean(Healthy_3_reshape);
Healthy_3_var=var(Healthy_3_reshape);
t_Healthy_3_wrapping=0+wrapping:wrapping:floor(length(Healthy_3)/wrapping)*wrapping;
%Healthy_4
wrapping=windo; %defined in beat
num_seg=floor(length(Healthy_4_corr)/wrapping) ;
Healthy_4_cut=Healthy_4_corr(1:num_seg*wrapping);
Healthy_4_reshape=reshape(Healthy_4_cut, [wrapping, length(Healthy_4_cut)/wrapping]);
Healthy_4_mean=mean(Healthy_4_reshape);
Healthy_4_var=var(Healthy_4_reshape);
t_Healthy_4_wrapping=0+wrapping:wrapping:floor(length(Healthy_4)/wrapping)*wrapping;
%Healthy_5
wrapping=windo; %defined in beat
num_seg=floor(length(Healthy_5_corr)/wrapping) ;
Healthy_5_cut=Healthy_5_corr(1:num_seg*wrapping);
Healthy_5_reshape=reshape(Healthy_5_cut, [wrapping, length(Healthy_5_cut)/wrapping]);
Healthy_5_mean=mean(Healthy_5_reshape);
Healthy_5_var=var(Healthy_5_reshape);
t_Healthy_5_wrapping=0+wrapping:wrapping:floor(length(Healthy_5)/wrapping)*wrapping;

%IUGR_1
wrapping=windo; %defined in beat
num_seg=floor(length(IUGR_1_corr)/wrapping) ;
IUGR_1_cut=IUGR_1_corr(1:num_seg*wrapping);

IUGR_1_reshape=reshape(IUGR_1_cut, [wrapping, length(IUGR_1_cut)/wrapping]);
IUGR_1_mean=mean(IUGR_1_reshape);
IUGR_1_var=var(IUGR_1_reshape);
t_IUGR_1_wrapping=0+wrapping:wrapping:floor(length(IUGR_1)/wrapping)*wrapping;

%IUGR_2
wrapping=windo; %defined in beat
num_seg=floor(length(IUGR_2_corr)/wrapping) ;
IUGR_2_cut=IUGR_2_corr(1:num_seg*wrapping);
IUGR_2_reshape=reshape(IUGR_2_cut, [wrapping, length(IUGR_2_cut)/wrapping]);
IUGR_2_mean=mean(IUGR_2_reshape);
IUGR_2_var=var(IUGR_2_reshape);
t_IUGR_2_wrapping=0+wrapping:wrapping:floor(length(IUGR_2)/wrapping)*wrapping;

%IUGR_3
wrapping=windo; %defined in beat
num_seg=floor(length(IUGR_3_corr)/wrapping) ;
IUGR_3_cut=IUGR_3_corr(1:num_seg*wrapping);
IUGR_3_reshape=reshape(IUGR_3_cut, [wrapping, length(IUGR_3_cut)/wrapping]);
IUGR_3_mean=mean(IUGR_3_reshape);
IUGR_3_var=var(IUGR_3_reshape);
t_IUGR_3_wrapping=0+wrapping:wrapping:floor(length(IUGR_3)/wrapping)*wrapping;

%IUGR_5
wrapping=windo; %defined in beat
num_seg=floor(length(IUGR_4_corr)/wrapping) ;
IUGR_4_cut=IUGR_4_corr(1:num_seg*wrapping);
IUGR_4_reshape=reshape(IUGR_4_cut, [wrapping, length(IUGR_4_cut)/wrapping]);
IUGR_4_mean=mean(IUGR_4_reshape);
IUGR_4_var=var(IUGR_4_reshape);
t_IUGR_4_wrapping=0+wrapping:wrapping:floor(length(IUGR_4)/wrapping)*wrapping;
%IUGR_5
wrapping=windo; %defined in beat
num_seg=floor(length(IUGR_5_corr)/wrapping) ;
IUGR_5_cut=IUGR_5_corr(1:num_seg*wrapping);
IUGR_5_reshape=reshape(IUGR_5_cut, [wrapping, length(IUGR_5_cut)/wrapping]);
IUGR_5_mean=mean(IUGR_5_reshape);
IUGR_5_var=var(IUGR_5_reshape);
t_IUGR_5_wrapping=0+wrapping:wrapping:floor(length(IUGR_5)/wrapping)*wrapping;

Healthy_1_av_mean = mean(Healthy_1_mean,'omitnan');
Healthy_1_av_var = mean(Healthy_1_var,'omitnan');
Mean_H1=['Mean =', num2str(Healthy_1_av_mean)];
Var_H1 = ['Var =', num2str(Healthy_1_av_var)];

Healthy_2_av_mean = mean(Healthy_2_mean,'omitnan');
Healthy_2_av_var = mean(Healthy_2_var,'omitnan');
Mean_H2=['Mean =', num2str(Healthy_2_av_mean)];
Var_H2 = ['Var =', num2str(Healthy_2_av_var)];
Healthy_3_av_mean = mean(Healthy_3_mean,'omitnan');
Healthy_3_av_var = mean(Healthy_3_var,'omitnan');
Mean_H3=['Mean =', num2str(Healthy_3_av_mean)];
Var_H3 = ['Var =', num2str(Healthy_3_av_var)];
Healthy_4_av_mean = mean(Healthy_4_mean,'omitnan');
Healthy_4_av_var = mean(Healthy_4_var,'omitnan');
Mean_H4=['Mean =', num2str(Healthy_4_av_mean)];
Var_H4 = ['Var =', num2str(Healthy_4_av_var)];
Healthy_5_av_mean = mean(Healthy_5_mean,'omitnan');
Healthy_5_av_var = mean(Healthy_5_var,'omitnan');
Mean_H5=['Mean =', num2str(Healthy_5_av_mean)];
Var_H5= ['Var =', num2str(Healthy_5_av_var)];
IUGR_1_av_mean = mean(IUGR_1_mean,'omitnan');
IUGR_1_av_var = mean(IUGR_1_var,'omitnan');
Mean_I1=['Mean =', num2str(IUGR_1_av_mean)];
Var_I1 = ['Var =', num2str(IUGR_1_av_var)];
IUGR_2_av_mean = mean(IUGR_2_mean,'omitnan');
IUGR_2_av_var = mean(IUGR_2_var,'omitnan');
Mean_I2=['Mean =', num2str(IUGR_2_av_mean)];
Var_I2 = ['Var =', num2str(IUGR_2_av_var)];
IUGR_3_av_mean = mean(IUGR_3_mean,'omitnan');
IUGR_3_av_var = mean(IUGR_3_var,'omitnan');
Mean_I3=['Mean =', num2str(IUGR_3_av_mean)];
Var_I3 = ['Var =', num2str(IUGR_3_av_var)];
IUGR_4_av_mean = mean(IUGR_4_mean,'omitnan');
IUGR_4_av_var = mean(IUGR_4_var,'omitnan');
Mean_I4=['Mean =', num2str(IUGR_4_av_mean)];
Var_I4 = ['Var =', num2str(IUGR_4_av_var)];
IUGR_5_av_mean = mean(IUGR_5_mean,'omitnan');
IUGR_5_av_var = mean(IUGR_5_var,'omitnan');
Mean_I5=['Mean =', num2str(IUGR_5_av_mean)];
Var_I5= ['Var =', num2str(IUGR_5_av_var)];

Mean_Healthy = (Healthy_1_av_mean+Healthy_2_av_mean+Healthy_3_av_mean+Healthy_4_av_mean+Healthy_5_av_mean)/5;
Var_Healthy = (Healthy_1_av_var+Healthy_2_av_var+Healthy_3_av_var+Healthy_4_av_var+Healthy_5_av_var)/5;


figure;
subplot(2,5,1);plot(t_IUGR_1_wrapping, IUGR_1_mean, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('IUGR\_1 mean [\muV]');
subplot(2,5,6);plot(t_IUGR_1_wrapping, IUGR_1_var, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('IUGR\_1 variance [\muV^2]');
subplot(2,5,2);plot(t_IUGR_2_wrapping, IUGR_2_mean, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('IUGR\_2 mean [\muV]');
subplot(2,5,7);plot(t_IUGR_2_wrapping, IUGR_2_var, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('IUGR\_2 variance [\muV^2]');
subplot(2,5,3);plot(t_IUGR_3_wrapping, IUGR_3_mean, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('IUGR\_3 mean [\muV]');
subplot(2,5,8);plot(t_IUGR_3_wrapping, IUGR_3_var, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('IUGR\_3 variance [\muV^2]');
subplot(2,5,4);plot(t_IUGR_4_wrapping, IUGR_4_mean, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('IUGR\_4 mean [\muV]');
subplot(2,5,9);plot(t_IUGR_4_wrapping, IUGR_4_var, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('IUGR\_4 variance [\muV^2]');
subplot(2,5,5);plot(t_IUGR_5_wrapping, IUGR_5_mean, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('IUGR\_5 mean [\muV]');
subplot(2,5,10);plot(t_IUGR_5_wrapping, IUGR_5_var, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('IUGR\_5 variance [\muV^2]');

figure;
subplot(2,5,1);plot(t_Healthy_1_wrapping, Healthy_1_mean, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Healthy\_1 mean [\muV]');
subplot(2,5,6);plot(t_Healthy_1_wrapping, Healthy_1_var, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Healthy\_1 variance [\muV^2]');
subplot(2,5,2);plot(t_Healthy_2_wrapping, Healthy_2_mean, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Healthy\_2 mean [\muV]');
subplot(2,5,7);plot(t_Healthy_2_wrapping, Healthy_2_var, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Healthy\_2 variance [\muV^2]');
subplot(2,5,3);plot(t_Healthy_3_wrapping, Healthy_3_mean, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Healthy\_3 mean [\muV]');
subplot(2,5,8);plot(t_Healthy_3_wrapping, Healthy_3_var, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Healthy\_3 variance [\muV^2]');
subplot(2,5,4);plot(t_Healthy_4_wrapping, Healthy_4_mean, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Healthy\_4 mean [\muV]');
subplot(2,5,9);plot(t_Healthy_4_wrapping, Healthy_4_var, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Healthy\_4 variance [\muV^2]');
subplot(2,5,5);plot(t_Healthy_5_wrapping, Healthy_5_mean, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Healthy\_5 mean [\muV]');
subplot(2,5,10);plot(t_Healthy_5_wrapping, Healthy_5_var, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Healthy\_5 variance [\muV^2]');


%% Floating Line Analysis and Removal

%H1 
float_H1=zeros(windo-1,floor(length(Healthy_1)/windo));
for h=1:floor(length(Healthy_1)/windo)
    pp=0;
    t=0;
    if(isnan(Healthy_1_inter(:,h)))
        float_H1(:, h)=NaN;
    else
        pp=csaps(t_Healthy_1_inter(:,h), Healthy_1_inter(:,h),0.007);
        t=fnplt(pp);
        count=0;
        for i=1:length(t)
            if(ceil(t(1,i))== t(1,i))
                count=count+1;
                float_H1(count, h)=t(2,i);
            end
        end
    end
end

H1_de_float=zeros(windo-1,floor(length(Healthy_1)/windo));
for i=1:length(Healthy_1_inter(1,:))
% %     figure
%     subplot(2,1,1), plot(t_Healthy_1_inter(1:windo-1,i),float_H1(1:windo-1,i), 'b');
%     hold on, plot(t_Healthy_1_inter(1:windo-1,i),Healthy_1_inter(1:windo-1,i), 'r');
    H1_de_float(:,i)= (Healthy_1_inter(1:windo-1,i))- float_H1(:,i);
%     subplot(2,1,2), plot(t_Healthy_1_inter(1:windo-1,i),H1_de_float(:,i));
%     xlabel('Time in sec'), ylabel('FRHV in bpm');
    
end 

% Removing floating line of H2
float_H2=zeros(windo-1,floor(length(Healthy_2)/windo));
for h=1:floor(length(Healthy_2)/windo)
    pp=0;
    t=0;
    if(isnan(Healthy_2_inter(:,h)))
        float_H2(:, h)=NaN;
    else
        pp=csaps(t_Healthy_2_inter(:,h), Healthy_2_inter(:,h),0.007);
        t=fnplt(pp);
        count=0;
        for i=1:length(t)
            if(ceil(t(1,i))== t(1,i))
                count=count+1;
                float_H2(count, h)=t(2,i);
            end
        end
    end
end

H2_de_float=zeros(windo-1,floor(length(Healthy_2)/windo));
for i=1:length(Healthy_2_inter(1,:))
% %     figure
%     subplot(2,1,1), plot(t_Healthy_2_inter(1:windo-1,i),float_H2(1:windo-1,i), 'b');
%     hold on, plot(t_Healthy_2_inter(1:windo-1,i),Healthy_2_inter(1:windo-1,i), 'r');
    H2_de_float(:,i)= (Healthy_2_inter(1:windo-1,i))- float_H2(:,i);
%     subplot(2,1,2), plot(t_Healthy_2_inter(1:windo-1,i),H2_de_float(:,i));
%     xlabel('Time in sec'), ylabel('FRHV in bpm');
    
end

% Removing floating line of H3
float_H3=zeros(windo-1,floor(length(Healthy_3)/windo));
for h=1:floor(length(Healthy_3)/windo)
    pp=0;
    t=0;
    if(isnan(Healthy_3_inter(:,h)))
        float_H3(:, h)=NaN;
    else
        pp=csaps(t_Healthy_3_inter(:,h), Healthy_3_inter(:,h),0.007);
        t=fnplt(pp);
        count=0;
        for i=1:length(t)
            if(ceil(t(1,i))== t(1,i))
                count=count+1;
                float_H3(count, h)=t(2,i);
            end
        end
    end
end

H3_de_float=zeros(windo-1,floor(length(Healthy_3)/windo));
for i=1:length(Healthy_3_inter(1,:))
% %     figure
%     subplot(2,1,1), plot(t_Healthy_3_inter(1:windo-1,i),float_H3(1:windo-1,i), 'b');
%     hold on, plot(t_Healthy_3_inter(1:windo-1,i),Healthy_3_inter(1:windo-1,i), 'r');
    H3_de_float(:,i)= (Healthy_3_inter(1:windo-1,i))- float_H3(:,i);
%     subplot(2,1,2), plot(t_Healthy_3_inter(1:windo-1,i),H3_de_float(:,i));
%     xlabel('Time in sec'), ylabel('FRHV  H3 in bpm');
    
end

% Removing floating line of H4
float_H4=zeros(windo-1,floor(length(Healthy_4)/windo));
for h=1:floor(length(Healthy_4)/windo)
    pp=0;
    t=0;
    if(isnan(Healthy_4_inter(:,h)))
        float_H4(:, h)=NaN;
    else
        pp=csaps(t_Healthy_4_inter(:,h), Healthy_4_inter(:,h),0.007);
        t=fnplt(pp);
        count=0;
        for i=1:length(t)
            if(ceil(t(1,i))== t(1,i))
                count=count+1;
                float_H4(count, h)=t(2,i);
            end
        end
    end
end

H4_de_float=zeros(windo-1,floor(length(Healthy_4)/windo));
for i=1:length(Healthy_4_inter(1,:))
% %     figure
%     subplot(2,1,1), plot(t_Healthy_4_inter(1:windo-1,i),float_H4(1:windo-1,i), 'b');
%     hold on, plot(t_Healthy_4_inter(1:windo-1,i),Healthy_4_inter(1:windo-1,i), 'r');
    H4_de_float(:,i)= (Healthy_4_inter(1:windo-1,i))- float_H4(:,i);
%     subplot(2,1,2), plot(t_Healthy_4_inter(1:windo-1,i),H4_de_float(:,i));
%     xlabel('Time in sec'), ylabel('FRHV H4 in bpm');
    
end


%H5
float_H5=zeros(windo-1,floor(length(Healthy_5)/windo));
for h=1:floor(length(Healthy_5)/windo)
    pp=0;
    t=0;
    if(isnan(Healthy_5_inter(:,h)))
        float_H5(:, h)=NaN;
    else
        pp=csaps(t_Healthy_5_inter(:,h), Healthy_5_inter(:,h),0.007);
        t=fnplt(pp);
        count=0;
        for i=1:length(t)
            if(ceil(t(1,i))== t(1,i))
                count=count+1;
                float_H5(count, h)=t(2,i);
            end
        end
    end
end

H5_de_float=zeros(windo-1,floor(length(Healthy_5)/windo));
for i=1:length(Healthy_5_inter(1,:))
% %     figure
%     subplot(2,1,1), plot(t_Healthy_5_inter(1:windo-1,i),float_H5(1:windo-1,i), 'b');
%     hold on, plot(t_Healthy_5_inter(1:windo-1,i),Healthy_5_inter(1:windo-1,i), 'r');
    H5_de_float(:,i)= (Healthy_5_inter(1:windo-1,i))- float_H5(:,i);
%     subplot(2,1,2), plot(t_Healthy_5_inter(1:windo-1,i),H5_de_float(:,i));
%     xlabel('Time in sec'), ylabel('FRHV H5 in bpm');
    
end


%Removing floating line of IUGR_1
float_I1=zeros(windo-1,floor(length(IUGR_1)/windo));
for h=1:floor(length(IUGR_1)/windo)
    pp=0;
    t=0;
    if(isnan(IUGR_1_inter(:,h)))
        float_I1(:, h)=NaN;
    else
        pp=csaps(t_IUGR_1_inter(:,h), IUGR_1_inter(:,h),0.007);
        t=fnplt(pp);
        count=0;
        for i=1:length(t)
            if(ceil(t(1,i))== t(1,i))
                count=count+1;
                float_I1(count, h)=t(2,i);
            end
        end
    end
end

I1_de_float=zeros(windo-1,floor(length(IUGR_1)/windo));
for i=1:length(IUGR_1_inter(1,:))
%     figure
%     subplot(2,1,1), plot(t_IUGR_1_inter(1:windo-1,i),float_I1(1:windo-1,i), 'b');
%     hold on, plot(t_IUGR_1_inter(1:windo-1,i),IUGR_1_inter(1:windo-1,i), 'r');
    I1_de_float(:,i)= (IUGR_1_inter(1:windo-1,i))- float_I1(:,i);
%     subplot(2,1,2), plot(t_IUGR_1_inter(1:windo-1,i),I1_de_float(:,i));
%     xlabel('Time in sec'), ylabel('FRHV in bpm');
%     
end

%Removing floating line of IUGR_2
float_I2=zeros(windo-1,floor(length(IUGR_2)/windo));
for h=1:floor(length(IUGR_2)/windo)
    pp=0;
    t=0;
    if(isnan(IUGR_2_inter(:,h)))
        float_I2(:, h)=NaN;
    else
        pp=csaps(t_IUGR_2_inter(:,h), IUGR_2_inter(:,h),0.007);
        t=fnplt(pp);
        count=0;
        for i=1:length(t)
            if(ceil(t(1,i))== t(1,i))
                count=count+1;
                float_I2(count, h)=t(2,i);
            end
        end
    end
end

I2_de_float=zeros(windo-1,floor(length(IUGR_2)/windo));
for i=1:length(IUGR_2_inter(1,:))
%     figure
%     subplot(2,1,1), plot(t_IUGR_2_inter(1:windo-1,i),float_I2(1:windo-1,i), 'b');
%     hold on, plot(t_IUGR_2_inter(1:windo-1,i),IUGR_2_inter(1:windo-1,i), 'r');
    I2_de_float(:,i)= (IUGR_2_inter(1:windo-1,i))- float_I2(:,i);
%     subplot(2,1,2), plot(t_IUGR_2_inter(1:windo-1,i),I2_de_float(:,i));
%     xlabel('Time in sec'), ylabel('FRHV in bpm');
    
end


%Removing floating line of IUGR_3
float_I3=zeros(windo-1,floor(length(IUGR_3)/windo));
for h=1:floor(length(IUGR_3)/windo)
    pp=0;
    t=0;
    if(isnan(IUGR_3_inter(:,h)))
        float_I3(:, h)=NaN;
    else
        pp=csaps(t_IUGR_3_inter(:,h), IUGR_3_inter(:,h),0.007);
        t=fnplt(pp);
        count=0;
        for i=1:length(t)
            if(ceil(t(1,i))== t(1,i))
                count=count+1;
                float_I3(count, h)=t(2,i);
            end
        end
    end
end

I3_de_float=zeros(windo-1,floor(length(IUGR_3)/windo));
for i=1:length(IUGR_3_inter(1,:))
%     figure
%     subplot(2,1,1), plot(t_IUGR_3_inter(1:windo-1,i),float_I3(1:windo-1,i), 'b');
%     hold on, plot(t_IUGR_3_inter(1:windo-1,i),IUGR_3_inter(1:windo-1,i), 'r');
    I3_de_float(:,i)= (IUGR_3_inter(1:windo-1,i))- float_I3(:,i);
%     subplot(2,1,2), plot(t_IUGR_3_inter(1:windo-1,i),I3_de_float(:,i));
%     xlabel('Time in sec'), ylabel('FRHV IUGR 3 in bpm');
    
end

%Removing floating line of IUGR_4
float_I4=zeros(windo-1,floor(length(IUGR_4)/windo));
for h=1:floor(length(IUGR_4)/windo)
    pp=0;
    t=0;
    if(isnan(IUGR_4_inter(:,h)))
        float_I4(:, h)=NaN;
    else
        pp=csaps(t_IUGR_4_inter(:,h), IUGR_4_inter(:,h),0.007);
        t=fnplt(pp);
        count=0;
        for i=1:length(t)
            if(ceil(t(1,i))== t(1,i))
                count=count+1;
                float_I4(count, h)=t(2,i);
            end
        end
    end
end

I4_de_float=zeros(windo-1,floor(length(IUGR_4)/windo));
for i=1:length(IUGR_4_inter(1,:))
%     figure
%     subplot(2,1,1), plot(t_IUGR_4_inter(1:windo-1,i),float_I4(1:windo-1,i), 'b');
%     hold on, plot(t_IUGR_4_inter(1:windo-1,i),IUGR_4_inter(1:windo-1,i), 'r');
    I4_de_float(:,i)= (IUGR_4_inter(1:windo-1,i))- float_I4(:,i);
%     subplot(2,1,2), plot(t_IUGR_4_inter(1:windo-1,i),I4_de_float(:,i));
%     xlabel('Time in sec'), ylabel('FRHV  I4 in bpm');
    
end

%Removing floating line of IUGR_5
float_I5=zeros(windo-1,floor(length(IUGR_5)/windo));
for h=1:floor(length(IUGR_5)/windo)
    pp=0;
    t=0;
    if(isnan(IUGR_5_inter(:,h)))
        float_I5(:, h)=NaN;
    else
        pp=csaps(t_IUGR_5_inter(:,h), IUGR_5_inter(:,h),0.007);
        t=fnplt(pp);
        count=0;
        for i=1:length(t)
            if(ceil(t(1,i))== t(1,i))
                count=count+1;
                float_I5(count, h)=t(2,i);
            end
        end
    end
end

I5_de_float=zeros(windo-1,floor(length(IUGR_5)/windo));
for i=1:length(IUGR_5_inter(1,:))
%     figure
%     subplot(2,1,1), plot(t_IUGR_5_inter(1:windo-1,i),float_I5(1:windo-1,i), 'b');
%     hold on, plot(t_IUGR_5_inter(1:windo-1,i),IUGR_5_inter(1:windo-1,i), 'r');
    I5_de_float(:,i)= (IUGR_5_inter(1:windo-1,i))- float_I5(:,i);
%     subplot(2,1,2), plot(t_IUGR_5_inter(1:windo-1,i),I5_de_float(:,i));
%     xlabel('Time in sec'), ylabel('FRHV I5 in bpm');
    
end

%Plot an example of the floating line with the signal
figure
plot(t_Healthy_1_inter(1:windo-1,:), float_H1(:,:),'LineWidth', 1.90, 'Color', 'b'), xlabel('Time (s)'), ylabel('bpm'),title('Floating Line superimposed to FHR signal'); 
hold on, plot(t_Healthy_1_inter(1:windo-1,:),Healthy_1_inter(1:windo-1,:), 'r'), hold off;

% Mean and Variance of defloated signals
H1_defloat_mean = mean(H1_de_float,'omitnan');
H1_defloat_var = var(H1_de_float,'omitnan');
Mean_H1_defloat = mean(H1_defloat_mean);
Var_H1_defloat = var(H1_defloat_var);
Mean_H1=['Mean =', num2str(Mean_H1_defloat)];
Var_H1 = ['Var =', num2str(Var_H1_defloat)];

H2_defloat_mean = mean(H2_de_float,'omitnan');
H2_defloat_var = var(H2_de_float,'omitnan');
Mean_H2_defloat = mean(H2_defloat_mean);
Var_H2_defloat = var(H2_defloat_var);
Mean_H2=['Mean =', num2str(Mean_H2_defloat)];
Var_H2 = ['Var =', num2str(Var_H2_defloat)];

H3_defloat_mean = mean(H3_de_float,'omitnan');
H3_defloat_var = var(H3_de_float,'omitnan');
Mean_H3_defloat = mean(H3_defloat_mean);
Var_H3_defloat = var(H3_defloat_var);
Mean_H3=['Mean =', num2str(Mean_H3_defloat)];
Var_H3 = ['Var =', num2str(Var_H3_defloat)];

H4_defloat_mean = mean(H4_de_float,'omitnan');
H4_defloat_var = var(H4_de_float,'omitnan');
Mean_H4_defloat = mean(H4_defloat_mean);
Var_H4_defloat = var(H4_defloat_var);
Mean_H4=['Mean =', num2str(Mean_H4_defloat)];
Var_H4 = ['Var =', num2str(Var_H4_defloat)];

H5_defloat_mean = mean(H5_de_float,'omitnan');
H5_defloat_var = var(H5_de_float,'omitnan');
Mean_H5_defloat = mean(H5_defloat_mean);
Var_H5_defloat = var(H5_defloat_var);
Mean_H5=['Mean =', num2str(Mean_H5_defloat)];
Var_H5 = ['Var =', num2str(Var_H5_defloat)];

I1_defloat_mean = mean(I1_de_float,'omitnan');
I1_defloat_var = var(I1_de_float,'omitnan');
Mean_I1_defloat = mean(I1_defloat_mean);
Var_I1_defloat = var(I1_defloat_var);
Mean_I1=['Mean =', num2str(Mean_I1_defloat)];
Var_I1 = ['Var =', num2str(Var_I1_defloat)];

I2_defloat_mean = mean(I2_de_float,'omitnan');
I2_defloat_var = var(I2_de_float,'omitnan');
Mean_I2_defloat = mean(I2_defloat_mean);
Var_I2_defloat = var(I2_defloat_var);
Mean_I2=['Mean =', num2str(Mean_I2_defloat)];
Var_I2 = ['Var =', num2str(Var_I2_defloat)];

I3_defloat_mean = mean(I3_de_float,'omitnan');
I3_defloat_var = var(I3_de_float,'omitnan');
Mean_I3_defloat = mean(I3_defloat_mean);
Var_I3_defloat = var(I3_defloat_var);
Mean_I3=['Mean =', num2str(Mean_I3_defloat)];
Var_I3 = ['Var =', num2str(Var_I3_defloat)];

I4_defloat_mean = mean(I4_de_float,'omitnan');
I4_defloat_var = var(I4_de_float,'omitnan');
Mean_I4_defloat = mean(I4_defloat_mean);
Var_I4_defloat = var(I4_defloat_var);
Mean_I4=['Mean =', num2str(Mean_I4_defloat)];
Var_I4 = ['Var =', num2str(Var_I4_defloat)];

I5_defloat_mean = mean(I5_de_float,'omitnan');
I5_defloat_var = var(I5_de_float,'omitnan');
Mean_I5_defloat = mean(I5_defloat_mean);
Var_I5_defloat = var(I5_defloat_var);
Mean_I5=['Mean =', num2str(Mean_I5_defloat)];
Var_I5 = ['Var =', num2str(Var_I5_defloat)];
%Plot all the signals without the floating line
figure
subplot(2,3,1), plot(t_Healthy_1_inter(1:windo-1,:),H1_de_float(:,:)), xlabel('Time in sec'), ylabel('FRHV of H1 in bpm');
subplot(2,3,2), plot(t_Healthy_2_inter(1:windo-1,:),H2_de_float(:,:)), xlabel('Time in sec'), ylabel('FRHV OF H2 in bpm');
subplot(2,3,3), plot(t_Healthy_3_inter(1:windo-1,:),H3_de_float(:,:)), xlabel('Time in sec'), ylabel('FRHV of H3 in bpm');
subplot(2,3,4), plot(t_Healthy_4_inter(1:windo-1,:),H4_de_float(:,:)), xlabel('Time in sec'), ylabel('FRHV of H4 in bpm');
subplot(2,3,3), plot(t_Healthy_5_inter(1:windo-1,:),H5_de_float(:,:)), xlabel('Time in sec'), ylabel('FRHV of H5 in bpm');

figure
subplot(2,3,1), plot(t_IUGR_1_inter(1:windo-1,:),I1_de_float(:,:)), xlabel('Time in sec'), ylabel('FRHV IUGR1 in bpm');
subplot(2,3,2), plot(t_IUGR_2_inter(1:windo-1,:),I2_de_float(:,:)), xlabel('Time in sec'), ylabel('FRHV IUGR2 in bpm');
subplot(2,3,3), plot(t_IUGR_3_inter(1:windo-1,:),I3_de_float(:,:)), xlabel('Time in sec'), ylabel('FRHV IUGR3 in bpm');
subplot(2,3,4), plot(t_IUGR_4_inter(1:windo-1,:),I4_de_float(:,:)), xlabel('Time in sec'), ylabel('FRHV IUGR4 in bpm');
subplot(2,3,5), plot(t_IUGR_5_inter(1:windo-1,:),I5_de_float(:,:)), xlabel('Time in sec'), ylabel('FRHV IUGR5 in bpm');

%Plot Mean and Variance in Time of the floating time
figure;
subplot(2,5,1);plot(t_IUGR_1_wrapping, I1_defloat_mean, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Mean [bpm]'),ylim([-1 1]); title('Defloated IUGR 1');
subplot(2,5,6);plot(t_IUGR_1_wrapping, I1_defloat_var, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Var [bpm^2]');ylim([0 20]);
subplot(2,5,2);plot(t_IUGR_2_wrapping, I2_defloat_mean, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Mean [bpm]');ylim([-1 1]);title('Defloated IUGR 2');
subplot(2,5,7);plot(t_IUGR_2_wrapping, I2_defloat_var, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Var [bpm^2]');ylim([0 20]);
subplot(2,5,3);plot(t_IUGR_3_wrapping, I3_defloat_mean, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Mean [bpm]');ylim([-1 1]);title('Defloated IUGR 3');
subplot(2,5,8);plot(t_IUGR_3_wrapping, I3_defloat_var, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Var [bpm^2]');ylim([0 20]);
subplot(2,5,4);plot(t_IUGR_4_wrapping, I4_defloat_mean, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Mean [bpm]');ylim([-1 1]);title('Defloated IUGR 4');
subplot(2,5,9);plot(t_IUGR_4_wrapping, I4_defloat_var, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Var [bpm^2]');ylim([0 20]);
subplot(2,5,5);plot(t_IUGR_5_wrapping, I5_defloat_mean, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Mean [bpm]');ylim([-1 1]);title('Defloated IUGR 5');
subplot(2,5,10);plot(t_IUGR_5_wrapping, I5_defloat_var, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Var [bpm^2]');ylim([0 20]);

figure;
subplot(2,5,1);plot(t_Healthy_1_wrapping, H1_defloat_mean, '-ok', 'MarkerFaceColor' ,'r'),
xlabel('Beat'), ylabel('Mean [bpm]'),ylim([-1 1]), title('Defloated Healhy 1');
subplot(2,5,6);plot(t_Healthy_1_wrapping, H1_defloat_var, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Var [bpm^2]');ylim([0 20]),
subplot(2,5,2);plot(t_Healthy_2_wrapping, H2_defloat_mean, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Mean [bpm]'),ylim([-1 1]),title('Defloated Healhy 2');
subplot(2,5,7);plot(t_Healthy_2_wrapping, H2_defloat_var, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Var [bpm^2]'); ylim([0 20])
subplot(2,5,3);plot(t_Healthy_3_wrapping, H3_defloat_mean, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Mean [bpm]'),ylim([-1 1]),title('Defloated Healhy 3');
subplot(2,5,8);plot(t_Healthy_3_wrapping, H3_defloat_var, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Var [bpm^2]'),ylim([0 20]);
subplot(2,5,4);plot(t_Healthy_4_wrapping, H4_defloat_mean, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Mean [bpm]'),ylim([-1 1]),title('Defloated Healhy 4');
subplot(2,5,9);plot(t_Healthy_4_wrapping, H4_defloat_var, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Var [bpm^2]'),ylim([0 20]);
subplot(2,5,5);plot(t_Healthy_5_wrapping, H5_defloat_mean, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Mean [bpm]'),ylim([-1 1]);title('Defloated Healhy 5');
subplot(2,5,10);plot(t_Healthy_5_wrapping, H5_defloat_var, '-ok', 'MarkerFaceColor' ,'r');
xlabel('Beat'); ylabel('Var[bpm^2]'),ylim([0 20]);



IUGR_defloated_means=[Mean_I1_defloat;Mean_I2_defloat;Mean_I3_defloat;Mean_I4_defloat;Mean_I5_defloat];
Healthy_defloated_means=[Mean_H1_defloat;Mean_H2_defloat;Mean_H3_defloat;Mean_H4_defloat;Mean_H5_defloat];
IUGR_defloated_vars=[Var_I1_defloat;Var_I2_defloat;Var_I3_defloat;Var_I4_defloat;Var_I5_defloat];
Healthy_defloated_vars=[Var_H1_defloat;Var_H2_defloat;Var_H3_defloat;Var_H4_defloat;Var_H5_defloat];

T =table(IUGR_defloated_means,Healthy_defloated_means,IUGR_defloated_vars,Healthy_defloated_vars,'VariableNames', {'Mean_IUGR_defloated','Mean_Healthy_defloated','Var_IUGR_defloated','Var_Healthy_defloated'});

%% PSD calculation

%Frequencty Bands
VLF = [0 0.03]; LF = [0.03 0.15]; MF = [0.15 0.5]; HF = [0.5 1];
labels= {'VLF' 'LF', 'MF', 'HF'};


%Healthy_1 PSD
PSDH1_YW_mean = 0; count = 0;
H1_de_float = detrend(H1_de_float);
for i = 1:length(H1_de_float(1,:))
    if isnan(H1_de_float(:,i))
    else
        count = count + 1;
        [PSDH1_YW, fH1_YW] = pyulear(H1_de_float(:,i),15,256,2,'onesided');
    end
    PSDH1_YW_mean = PSDH1_YW_mean + PSDH1_YW;
    PSDH1_YW = 0;
end
PSDH1_YW_mean = PSDH1_YW_mean./count;

%PARAMETERS
f_VLF = and(ge(fH1_YW,VLF(1)),le(fH1_YW,VLF(2)));
f_LF = and(ge(fH1_YW,LF(1)),le(fH1_YW,LF(2)));
f_MF = and(ge(fH1_YW,MF(1)),le(fH1_YW,MF(2)));
f_HF = and(ge(fH1_YW,HF(1)),le(fH1_YW,HF(2)));

VLF_H1_YW=trapz(PSDH1_YW_mean(f_VLF));
LF_H1_YW=trapz(PSDH1_YW_mean(f_LF));
MF_H1_YW=trapz(PSDH1_YW_mean(f_MF));
HF_H1_YW=trapz(PSDH1_YW_mean(f_HF));
freq_range_H1 = [VLF_H1_YW LF_H1_YW MF_H1_YW HF_H1_YW];
%VLF_H1_YW

LFsuHFandMF_H1=LF_H1_YW/(HF_H1_YW + MF_H1_YW);

%Healthy_2 PSD
PSDH2_YW_mean = 0; count = 0;
H2_de_float = detrend(H2_de_float);
for i = 1:length(H2_de_float(1,:))
    if isnan(H2_de_float(:,i))
    else
        count = count + 1;
        [PSDH2_YW, fH2_YW] = pyulear(H2_de_float(:,i),15,256,2,'onesided');
    end
    PSDH2_YW_mean = PSDH2_YW_mean + PSDH2_YW;
    PSDH2_YW = 0;
end
PSDH2_YW_mean = PSDH2_YW_mean./count;

%PARAMETERS
f_VLF = and(ge(fH2_YW,VLF(1)),le(fH2_YW,VLF(2)));
f_LF = and(ge(fH2_YW,LF(1)),le(fH2_YW,LF(2)));
f_MF = and(ge(fH2_YW,MF(1)),le(fH2_YW,MF(2)));
f_HF = and(ge(fH2_YW,HF(1)),le(fH2_YW,HF(2)));

VLF_H2_YW=trapz(PSDH2_YW_mean(f_VLF));
LF_H2_YW=trapz(PSDH2_YW_mean(f_LF));
MF_H2_YW=trapz(PSDH2_YW_mean(f_MF));
HF_H2_YW=trapz(PSDH2_YW_mean(f_HF));
freq_range_H2 = [VLF_H2_YW LF_H2_YW MF_H2_YW HF_H2_YW];

LFsuHFandMF_H2=LF_H2_YW/(HF_H2_YW + MF_H2_YW);

%Healthy_3 PSD
PSDH3_YW_mean = 0; count = 0;
H3_de_float = detrend(H3_de_float);
for i = 1:length(H3_de_float(1,:))
    if isnan(H3_de_float(:,i))
    else
        count = count + 1;
        [PSDH3_YW, fH3_YW] = pyulear(H3_de_float(:,i),15,256,2,'onesided');
    end
    PSDH3_YW_mean = PSDH3_YW_mean + PSDH3_YW;
    PSDH3_YW = 0;
end
PSDH3_YW_mean = PSDH3_YW_mean./count;

%PARAMETERS
f_VLF = and(ge(fH3_YW,VLF(1)),le(fH3_YW,VLF(2)));
f_LF = and(ge(fH3_YW,LF(1)),le(fH3_YW,LF(2)));
f_MF = and(ge(fH3_YW,MF(1)),le(fH3_YW,MF(2)));
f_HF = and(ge(fH3_YW,HF(1)),le(fH3_YW,HF(2)));

VLF_H3_YW=trapz(PSDH3_YW_mean(f_VLF));
LF_H3_YW=trapz(PSDH3_YW_mean(f_LF));
MF_H3_YW=trapz(PSDH3_YW_mean(f_MF));
HF_H3_YW=trapz(PSDH3_YW_mean(f_HF));
freq_range_H3 = [VLF_H3_YW LF_H3_YW MF_H3_YW HF_H3_YW];

LFsuHFandMF_H3=LF_H3_YW/(HF_H3_YW + MF_H3_YW);

%Healthy_4 PSD
PSDH4_YW_mean = 0; count = 0;
H4_de_float = detrend(H4_de_float);
for i = 1:length(H4_de_float(1,:))
    if isnan(H4_de_float(:,i))
    else
        count = count + 1;
        [PSDH4_YW, fH4_YW] = pyulear(H4_de_float(:,i),15,256,2,'onesided');
    end
    PSDH4_YW_mean = PSDH4_YW_mean + PSDH4_YW;
    PSDH4_YW = 0;
end
PSDH4_YW_mean = PSDH4_YW_mean./count;

%PARAMETERS
f_VLF = and(ge(fH4_YW,VLF(1)),le(fH4_YW,VLF(2)));
f_LF = and(ge(fH4_YW,LF(1)),le(fH4_YW,LF(2)));
f_MF = and(ge(fH4_YW,MF(1)),le(fH4_YW,MF(2)));
f_HF = and(ge(fH4_YW,HF(1)),le(fH4_YW,HF(2)));

VLF_H4_YW=trapz(PSDH4_YW_mean(f_VLF));
LF_H4_YW=trapz(PSDH4_YW_mean(f_LF));
MF_H4_YW=trapz(PSDH4_YW_mean(f_MF));
HF_H4_YW=trapz(PSDH4_YW_mean(f_HF));
freq_range_H4 = [VLF_H4_YW LF_H4_YW MF_H4_YW HF_H4_YW];

LFsuHFandMF_H4=LF_H4_YW/(HF_H4_YW + MF_H4_YW);

%Healthy_5 PSD
PSDH5_YW_mean = 0; count = 0;
H5_de_float = detrend(H5_de_float);
for i = 1:length(H5_de_float(1,:))
    if isnan(H5_de_float(:,i))
    else
        count = count + 1;
        [PSDH5_YW, fH5_YW] = pyulear(H5_de_float(:,i),15,256,2,'onesided');
    end
    PSDH5_YW_mean = PSDH5_YW_mean + PSDH5_YW;
    PSDH5_YW = 0;
end
PSDH5_YW_mean = PSDH5_YW_mean./count;

%PARAMETERS
f_VLF = and(ge(fH5_YW,VLF(1)),le(fH5_YW,VLF(2)));
f_LF = and(ge(fH5_YW,LF(1)),le(fH5_YW,LF(2)));
f_MF = and(ge(fH5_YW,MF(1)),le(fH5_YW,MF(2)));
f_HF = and(ge(fH5_YW,HF(1)),le(fH5_YW,HF(2)));

VLF_H5_YW=trapz(PSDH5_YW_mean(f_VLF));
LF_H5_YW=trapz(PSDH5_YW_mean(f_LF));
MF_H5_YW=trapz(PSDH5_YW_mean(f_MF));
HF_H5_YW=trapz(PSDH5_YW_mean(f_HF));
freq_range_H5 = [VLF_H5_YW LF_H5_YW MF_H5_YW HF_H5_YW];

LFsuHFandMF_H5=LF_H5_YW/(HF_H5_YW + MF_H5_YW);

%IUGR_1 PSD
PSDI1_YW_mean = 0; count = 0;
I1_de_float = detrend(I1_de_float);
for i = 1:length(I1_de_float(1,:))
    if isnan(I1_de_float(:,i))
    else
        count = count + 1;
        [PSDI1_YW, fI1_YW] = pyulear(I1_de_float(:,i),15,256,2,'onesided');
    end
    PSDI1_YW_mean = PSDI1_YW_mean + PSDI1_YW;
    PSDI1_YW = 0;
end
PSDI1_YW_mean = PSDI1_YW_mean./count;

%PARAMETERS
f_VLF = and(ge(fI1_YW,VLF(1)),le(fI1_YW,VLF(2)));
f_LF = and(ge(fI1_YW,LF(1)),le(fI1_YW,LF(2)));
f_MF = and(ge(fI1_YW,MF(1)),le(fI1_YW,MF(2)));
f_HF = and(ge(fI1_YW,HF(1)),le(fI1_YW,HF(2)));

VLF_I1_YW=trapz(PSDI1_YW_mean(f_VLF));
LF_I1_YW=trapz(PSDI1_YW_mean(f_LF));
MF_I1_YW=trapz(PSDI1_YW_mean(f_MF));
HF_I1_YW=trapz(PSDI1_YW_mean(f_HF));
freq_range_I1 = [VLF_I1_YW LF_I1_YW MF_I1_YW HF_I1_YW];

LFsuHFandMF_I1=LF_I1_YW/(HF_I1_YW + MF_I1_YW);

%IUGR_2 PSD
PSDI2_YW_mean = 0; count = 0;
I2_de_float = detrend(I2_de_float);
for i = 1:length(I2_de_float(1,:))
    if isnan(I2_de_float(:,i))
    else
        count = count + 1;
        [PSDI2_YW, fI2_YW] = pyulear(I2_de_float(:,i),15,256,2,'onesided');
    end
    PSDI2_YW_mean = PSDI2_YW_mean + PSDI2_YW;
    PSDI2_YW = 0;
end
PSDI2_YW_mean = PSDI2_YW_mean./count;

%PARAMETERS
f_VLF = and(ge(fI2_YW,VLF(1)),le(fI2_YW,VLF(2)));
f_LF = and(ge(fI2_YW,LF(1)),le(fI2_YW,LF(2)));
f_MF = and(ge(fI2_YW,MF(1)),le(fI2_YW,MF(2)));
f_HF = and(ge(fI2_YW,HF(1)),le(fI2_YW,HF(2)));

VLF_I2_YW=trapz(PSDI2_YW_mean(f_VLF));
LF_I2_YW=trapz(PSDI2_YW_mean(f_LF));
MF_I2_YW=trapz(PSDI2_YW_mean(f_MF));
HF_I2_YW=trapz(PSDI2_YW_mean(f_HF));
freq_range_I2 = [VLF_I2_YW LF_I2_YW MF_I2_YW HF_I2_YW];

LFsuHFandMF_I2=LF_I2_YW/(HF_I2_YW + MF_I2_YW);

%IUGR_3 PSD
PSDI3_YW_mean = 0; count = 0;
I3_de_float = detrend(I3_de_float);
for i = 1:length(I3_de_float(1,:))
    if isnan(I3_de_float(:,i))
    else
        count = count + 1;
        [PSDI3_YW, fI3_YW] = pyulear(I3_de_float(:,i),15,256,2,'onesided');
    end
    PSDI3_YW_mean = PSDI3_YW_mean + PSDI3_YW;
    PSDI3_YW = 0;
end
PSDI3_YW_mean = PSDI3_YW_mean./count;

%PARAMETERS
f_VLF = and(ge(fI1_YW,VLF(1)),le(fI1_YW,VLF(2)));
f_LF = and(ge(fI1_YW,LF(1)),le(fI1_YW,LF(2)));
f_MF = and(ge(fI1_YW,MF(1)),le(fI1_YW,MF(2)));
f_HF = and(ge(fI1_YW,HF(1)),le(fI1_YW,HF(2)));

VLF_I3_YW=trapz(PSDI3_YW_mean(f_VLF));
LF_I3_YW=trapz(PSDI3_YW_mean(f_LF));
MF_I3_YW=trapz(PSDI3_YW_mean(f_MF));
HF_I3_YW=trapz(PSDI3_YW_mean(f_HF));
freq_range_I3 = [VLF_I3_YW LF_I3_YW MF_I3_YW HF_I3_YW];

LFsuHFandMF_I3=LF_I3_YW/(HF_I3_YW + MF_I3_YW);

%IUGR 4 PSD
PSDI4_YW_mean = 0; count = 0;
I4_de_float = detrend(I4_de_float);
for i = 1:length(I4_de_float(1,:))
    if isnan(I4_de_float(:,i))
    else
        count = count + 1;
        [PSDI4_YW, fI4_YW] = pyulear(I4_de_float(:,i),15,256,2,'onesided');
    end
    PSDI4_YW_mean = PSDI4_YW_mean + PSDI4_YW;
    PSDI4_YW = 0;
end
PSDI4_YW_mean = PSDI4_YW_mean./count;

%PARAMETERS
f_VLF = and(ge(fI1_YW,VLF(1)),le(fI1_YW,VLF(2)));
f_LF = and(ge(fI1_YW,LF(1)),le(fI1_YW,LF(2)));
f_MF = and(ge(fI1_YW,MF(1)),le(fI1_YW,MF(2)));
f_HF = and(ge(fI1_YW,HF(1)),le(fI1_YW,HF(2)));

VLF_I4_YW=trapz(PSDI4_YW_mean(f_VLF));
LF_I4_YW=trapz(PSDI4_YW_mean(f_LF));
MF_I4_YW=trapz(PSDI4_YW_mean(f_MF));
HF_I4_YW=trapz(PSDI4_YW_mean(f_HF));
freq_range_I4 = [VLF_I4_YW LF_I4_YW MF_I4_YW HF_I4_YW];

LFsuHFandMF_I4=LF_I4_YW/(HF_I4_YW + MF_I4_YW);

%IUGR 5 PSD
PSDI5_YW_mean = 0; count = 0; PSDI5_YW=0;
I5_de_float = detrend(I5_de_float);
for i = 1:length(I5_de_float(1,:))
    if isnan(I5_de_float(:,i))
    else
        count = count + 1;
        [PSDI5_YW, fI5_YW] = pyulear(I5_de_float(:,i),15,256,2);
    end
    PSDI5_YW_mean = PSDI5_YW_mean + PSDI5_YW;
    PSDI5_YW = 0;
end
PSDI5_YW_mean = PSDI5_YW_mean./count;

%IUGR 5 PARAMETERS 
f_VLF = and(ge(fI5_YW,VLF(1)),le(fI5_YW,VLF(2)));
f_LF = and(ge(fI5_YW,LF(1)),le(fI5_YW,LF(2)));
f_MF = and(ge(fI5_YW,MF(1)),le(fI5_YW,MF(2)));
f_HF = and(ge(fI5_YW,HF(1)),le(fI5_YW,HF(2)));

VLF_I5_YW=trapz(PSDI5_YW_mean(f_VLF));
LF_I5_YW=trapz(PSDI5_YW_mean(f_LF));
MF_I5_YW=trapz(PSDI5_YW_mean(f_MF));
HF_I5_YW=trapz(PSDI5_YW_mean(f_HF));
freq_range_I5 = [VLF_I5_YW LF_I5_YW MF_I5_YW HF_I5_YW];

LFsuHFandMF_I5=LF_I5_YW/(HF_I5_YW + MF_I5_YW);

%Plot all results

figure;
subplot(5,2,1),plot(fH1_YW, PSDH1_YW_mean), title('PSD Healty1'), xlabel('Frequency (Hz)'); ylabel('bps^2/Hz');
str= ['LF/(HF+MF)=',num2str(LFsuHFandMF_H1)]; text(0.6,0.8*max(PSDH1_YW_mean),str,'HorizontalAlignment','left');
subplot(5,2,2), pie(freq_range_H1), legend(labels,'Location','best');
subplot(5,2,3),plot(fH2_YW, PSDH2_YW_mean), title('PSD Healty2'), xlabel('Frequency (Hz)'); ylabel('bps^2/Hz');
str= ['LF/(HF+MF)=',num2str(LFsuHFandMF_H2)]; text(0.6,0.8*max(PSDH2_YW_mean),str,'HorizontalAlignment','left');
subplot(5,2,4),pie(freq_range_H2);
subplot(5,2,5),plot(fH3_YW, PSDH3_YW_mean), title('PSD Healty3'), xlabel('Frequency (Hz)'); ylabel('bps^2/Hz');
str= ['LF/(HF+MF)=',num2str(LFsuHFandMF_H3)]; text(0.6,0.8*max(PSDH3_YW_mean),str,'HorizontalAlignment','left');
subplot(5,2,6),pie(freq_range_H3);
subplot(5,2,7),plot(fH4_YW, PSDH4_YW_mean), title('PSD Healty4'), xlabel('Frequency (Hz)'); ylabel('bps^2/Hz');
str= ['LF/(HF+MF)=',num2str(LFsuHFandMF_H4)]; text(0.6,0.8*max(PSDH4_YW_mean),str,'HorizontalAlignment','left');
subplot(5,2,8),pie(freq_range_H4);
subplot(5,2,9),plot(fH5_YW, PSDH5_YW_mean), title('PSD Healty5'), xlabel('Frequency (Hz)'); ylabel('bps^2/Hz');
str= ['LF/(HF+MF)=',num2str(LFsuHFandMF_H5)]; text(0.6,0.8*max(PSDH5_YW_mean),str,'HorizontalAlignment','left');
subplot(5,2,10),pie(freq_range_H5);

figure;
subplot(5,2,1),plot(fI1_YW, PSDI1_YW_mean), title('PSD IUGR 1'), xlabel('Frequency (Hz)'); ylabel('bps^2/Hz');
str= ['LF/(HF+MF)=',num2str(LFsuHFandMF_I1)]; text(0.6,0.8*max(PSDI1_YW_mean),str,'HorizontalAlignment','left');
subplot(5,2,2), pie(freq_range_I1), legend(labels,'Location','best');
subplot(5,2,3),plot(fI2_YW, PSDI2_YW_mean), title('PSD IUGR 2'), xlabel('Frequency (Hz)'); ylabel('bps^2/Hz');
str= ['LF/(HF+MF)=',num2str(LFsuHFandMF_I2)]; text(0.6,0.8*max(PSDI2_YW_mean),str,'HorizontalAlignment','left');
subplot(5,2,4),pie(freq_range_I2);
subplot(5,2,5),plot(fI3_YW, PSDI3_YW_mean), title('PSD IUGR 3'), xlabel('Frequency (Hz)'); ylabel('bps^2/Hz');
str= ['LF/(HF+MF)=',num2str(LFsuHFandMF_I3)]; text(0.6,0.8*max(PSDI3_YW_mean),str,'HorizontalAlignment','left');
subplot(5,2,6),pie(freq_range_I3);
subplot(5,2,7),plot(fI4_YW, PSDI4_YW_mean), title('PSD IUGR4'), xlabel('Frequency (Hz)'); ylabel('bps^2/Hz');
str= ['LF/(HF+MF)=',num2str(LFsuHFandMF_I4)]; text(0.6,0.8*max(PSDI4_YW_mean),str,'HorizontalAlignment','left');
subplot(5,2,8),pie(freq_range_I4);
subplot(5,2,9),plot(fI5_YW, PSDI5_YW_mean), title('PSD IUGR5'), xlabel('Frequency (Hz)'); ylabel('bps^2/Hz');
str= ['LF/(HF+MF)=',num2str(LFsuHFandMF_I5)]; text(0.6,0.8*max(PSDI5_YW_mean),str,'HorizontalAlignment','left');
subplot(5,2,10),pie(freq_range_I5);

%Computing a Table with Results
SVB_Healthy= [ LFsuHFandMF_H1; LFsuHFandMF_H2; LFsuHFandMF_H3; LFsuHFandMF_H4; LFsuHFandMF_H5];
SVB_IUGR= [ LFsuHFandMF_I1; LFsuHFandMF_I2; LFsuHFandMF_I3; LFsuHFandMF_I4; LFsuHFandMF_I5];
SVB_mean_std= [ mean(SVB_Healthy) std(SVB_Healthy); mean(SVB_IUGR) std(SVB_IUGR)];
%SVB_I_mean_std=[mean(SVB_IUGR) std(SVB_IUGR)];
Row = {'Healthy';'IUGR'};
T_SVB=table(SVB_mean_std,'RowNames', Row);


%Complete table with svb mean 

%PSD FOR FLOATING LINE - assessing accelerations and decelerations
%H1 PSD
PSD_f_H1_YW_mean = 0; count = 0; PSD_f_H1_YW=0;
float_H1 = detrend(float_H1);
for i = 1:length(float_H1(1,:))
    if isnan(float_H1(:,i))
    else
        count = count + 1;
        [PSD_f_H1_YW, f_f_H1_YW] = pyulear(float_H1(:,i),15,256,2);
    end
    PSD_f_H1_YW_mean = PSD_f_H1_YW_mean + PSD_f_H1_YW;
   PSD_f_H1_YW = 0;
end
PSD_f_H1_YW_mean= PSD_f_H1_YW_mean./count;

%VLF 
f_f_VLF= and(ge(f_f_H1_YW,VLF(1)),le(f_f_H1_YW,VLF(2)));

VLF_f_H1_YW=trapz(PSD_f_H1_YW_mean(f_f_VLF));
freq_range_H1(1)=freq_range_H1(1)+ VLF_f_H1_YW;

%H2 PSD
PSD_f_H2_YW_mean = 0; count = 0; PSD_f_H2_YW=0;
float_H2 = detrend(float_H2);
for i = 1:length(float_H2(1,:))
    if isnan(float_H2(:,i))
    else
        count = count + 1;
        [PSD_f_H2_YW, f_f_H2_YW] = pyulear(float_H2(:,i),15,256,2);
    end
    PSD_f_H2_YW_mean = PSD_f_H2_YW_mean + PSD_f_H2_YW;
   PSD_f_H2_YW = 0;
end
PSD_f_H2_YW_mean= PSD_f_H2_YW_mean./count;

%VLF 
f_f_VLF= and(ge(f_f_H2_YW,VLF(1)),le(f_f_H2_YW,VLF(2)));

VLF_f_H2_YW=trapz(PSD_f_H2_YW_mean(f_f_VLF));

freq_range_H2(1)=freq_range_H2(1)+ VLF_f_H2_YW;

%H3 PSD
PSD_f_H3_YW_mean = 0; count = 0; PSD_f_H3_YW=0;
float_H3 = detrend(float_H3);
for i = 1:length(float_H3(1,:))
    if isnan(float_H3(:,i))
    else
        count = count + 1;
        [PSD_f_H3_YW, f_f_H3_YW] = pyulear(float_H3(:,i),15,256,2);
    end
    PSD_f_H3_YW_mean = PSD_f_H3_YW_mean + PSD_f_H3_YW;
   PSD_f_H3_YW = 0;
end
PSD_f_H3_YW_mean= PSD_f_H3_YW_mean./count;

%VLF 
f_f_VLF= and(ge(f_f_H3_YW,VLF(1)),le(f_f_H3_YW,VLF(2)));

VLF_f_H3_YW=trapz(PSD_f_H3_YW_mean(f_f_VLF));

freq_range_H3(1)=freq_range_H3(1)+ VLF_f_H3_YW;

%H4 PSD
PSD_f_H4_YW_mean = 0; count = 0; PSD_f_H4_YW=0;
float_H4 = detrend(float_H4);
for i = 1:length(float_H4(1,:))
    if isnan(float_H4(:,i))
    else
        count = count + 1;
        [PSD_f_H4_YW, f_f_H4_YW] = pyulear(float_H4(:,i),15,256,2);
    end
    PSD_f_H4_YW_mean = PSD_f_H4_YW_mean + PSD_f_H4_YW;
   PSD_f_H4_YW = 0;
end
PSD_f_H4_YW_mean= PSD_f_H4_YW_mean./count;

%VLF 
f_f_VLF= and(ge(f_f_H4_YW,VLF(1)),le(f_f_H4_YW,VLF(2)));

VLF_f_H4_YW=trapz(PSD_f_H4_YW_mean(f_f_VLF));

freq_range_H4(1)=freq_range_H4(1)+ VLF_f_H4_YW;

%H5 PSD
PSD_f_H5_YW_mean = 0; count = 0; PSD_f_H5_YW=0;
float_H5 = detrend(float_H5);
for i = 1:length(float_H5(1,:))
    if isnan(float_H5(:,i))
    else
        count = count + 1;
        [PSD_f_H5_YW, f_f_H5_YW] = pyulear(float_H5(:,i),15,256,2);
    end
    PSD_f_H5_YW_mean = PSD_f_H5_YW_mean + PSD_f_H5_YW;
   PSD_f_H5_YW = 0;
end
PSD_f_H5_YW_mean= PSD_f_H5_YW_mean./count;

%VLF 
f_f_VLF= and(ge(f_f_H5_YW,VLF(1)),le(f_f_H5_YW,VLF(2)));

VLF_f_H5_YW=trapz(PSD_f_H5_YW_mean(f_f_VLF));
freq_range_H5(1)=freq_range_H5(1)+ VLF_f_H5_YW;

%IUGR 1 PSD 

PSD_f_I1_YW_mean = 0; count = 0; PSD_f_I1_YW=0;
float_I1 = detrend(float_I1);
for i = 1:length(float_I1(1,:))
    if isnan(float_I1(:,i))
    else
        count = count + 1;
        [PSD_f_I1_YW, f_f_I1_YW] = pyulear(float_I1(:,i),15,256,2);
    end
    PSD_f_I1_YW_mean = PSD_f_I1_YW_mean + PSD_f_I1_YW;
   PSD_f_I1_YW = 0;
end
PSD_f_I1_YW_mean= PSD_f_I1_YW_mean./count;

%VLF 
f_f_VLF = and(ge(f_f_I1_YW,VLF(1)),le(f_f_I1_YW,VLF(2)));

VLF_f_I1_YW=trapz(PSD_f_I1_YW_mean(f_f_VLF));
freq_range_I1(1)=freq_range_I1(1)+ VLF_f_I1_YW;

%IUGR 2 PSD
PSD_f_I2_YW_mean = 0; count = 0; PSD_f_I2_YW=0;
float_I2 = detrend(float_I2);
for i = 1:length(float_I2(1,:))
    if isnan(float_I2(:,i))
    else
        count = count + 1;
        [PSD_f_I2_YW, f_f_I2_YW] = pyulear(float_I2(:,i),15,256,2);
    end
    PSD_f_I2_YW_mean = PSD_f_I2_YW_mean + PSD_f_I2_YW;
   PSD_f_I2_YW = 0;
end
PSD_f_I2_YW_mean= PSD_f_I2_YW_mean./count;

%VLF 
f_f_VLF = and(ge(f_f_I2_YW,VLF(1)),le(f_f_I2_YW,VLF(2)));

VLF_f_I2_YW=trapz(PSD_f_I2_YW_mean(f_f_VLF));
freq_range_I2(1)=freq_range_I2(1)+ VLF_f_I2_YW;

%IUGR 3 PSD
PSD_f_I3_YW_mean = 0; count = 0; PSD_f_I3_YW=0;
float_I3 = detrend(float_I3);
for i = 1:length(float_I3(1,:))
    if isnan(float_I3(:,i))
    else
        count = count + 1;
        [PSD_f_I3_YW, f_f_I3_YW] = pyulear(float_I3(:,i),15,256,2);
    end
    PSD_f_I3_YW_mean = PSD_f_I3_YW_mean + PSD_f_I3_YW;
   PSD_f_I3_YW = 0;
end
PSD_f_I3_YW_mean= PSD_f_I3_YW_mean./count;

%VLF 
f_f_VLF = and(ge(f_f_I3_YW,VLF(1)),le(f_f_I3_YW,VLF(2)));

VLF_f_I3_YW=trapz(PSD_f_I3_YW_mean(f_f_VLF));
freq_range_I3(1)=freq_range_I3(1)+ VLF_f_I3_YW;

%IUGR 4 PSD
PSD_f_I4_YW_mean = 0; count = 0; PSD_f_I4_YW=0;
float_I4 = detrend(float_I4);
for i = 1:length(float_I4(1,:))
    if isnan(float_I4(:,i))
    else
        count = count + 1;
        [PSD_f_I4_YW, f_f_I4_YW] = pyulear(float_I4(:,i),15,256,2);
    end
    PSD_f_I4_YW_mean = PSD_f_I4_YW_mean + PSD_f_I4_YW;
   PSD_f_I4_YW = 0;
end
PSD_f_I4_YW_mean= PSD_f_I4_YW_mean./count;

%VLF 
f_f_VLF = and(ge(f_f_I4_YW,VLF(1)),le(f_f_I4_YW,VLF(2)));

VLF_f_I4_YW=trapz(PSD_f_I4_YW_mean(f_f_VLF));
freq_range_I4(1)=freq_range_I4(1)+ VLF_f_I4_YW;

%IUGR 5 PSD
PSD_f_I5_YW_mean = 0; count = 0; PSD_f_I5_YW=0;
float_I5 = detrend(float_I5);
for i = 1:length(float_I5(1,:))
    if isnan(float_I5(:,i))
    else
        count = count + 1;
        [PSD_f_I5_YW, f_f_I5_YW] = pyulear(float_I5(:,i),15,256,2);
    end
    PSD_f_I5_YW_mean = PSD_f_I5_YW_mean + PSD_f_I5_YW;
   PSD_f_I5_YW = 0;
end
PSD_f_I5_YW_mean= PSD_f_I5_YW_mean./count;

%VLF 
f_f_VLF = and(ge(f_f_I5_YW,VLF(1)),le(f_f_I5_YW,VLF(2)));

VLF_f_I5_YW=trapz(PSD_f_I5_YW_mean(f_f_VLF));
freq_range_I5(1)=freq_range_I5(1)+ VLF_f_I5_YW;


% %Example of Floating line for H1 and IUGR1
% figure
% subplot(2,1,1), plot(t_Healthy_1_inter(1:windo-1,:), float_H1), title('Floating line of H1 and IUGR1');
% subplot(2,1,2), plot(t_IUGR_3_inter(1:windo-1,:), float_I3), ylim([-20 40]);
% 

figure;
subplot(5,2,1),plot(f_f_H1_YW, PSD_f_H1_YW_mean), title('PSD floatH1'); 
xlabel('Frequency (Hz)'); ylabel('bpm^2/Hz'), xlim([0 0.1]);
str=['VLF=',num2str(VLF_f_H1_YW)], text(0.06,0.8*max(PSD_f_H1_YW_mean),str,'HorizontalAlignment','left');
subplot(5,2,3),plot(f_f_H2_YW, PSD_f_H2_YW_mean), title('PSD floatH2'), xlabel('Frequency (Hz)'); ylabel('bpm^2/Hz');
str= ['VLF=',num2str(VLF_f_H2_YW)]; text(0.06,0.8*max(PSD_f_H2_YW_mean),str,'HorizontalAlignment','left');
subplot(5,2,5),plot(f_f_H3_YW, PSD_f_H3_YW_mean), title('PSD floatH3'), xlabel('Frequency (Hz)'); ylabel('bpm^2/Hz');
str= ['VLF(bpm2)=',num2str(VLF_f_H3_YW)]; text(0.6,0.8*max(PSD_f_H3_YW_mean),str,'HorizontalAlignment','left');
subplot(5,2,7),plot(f_f_H4_YW, PSD_f_H4_YW_mean), title('PSD floatH4'), xlabel('Frequency (Hz)'); ylabel('bpm^2/Hz');
str= ['VLF=',num2str(VLF_f_H4_YW)]; text(0.6,0.8*max(PSD_f_H4_YW_mean),str,'HorizontalAlignment','left');
subplot(5,2,9),plot(f_f_H5_YW, PSD_f_H5_YW_mean), title('PSD floatH5'), xlabel('Frequency (Hz)'); ylabel('bpm^2/Hz');
str= ['VLF=',num2str(VLF_f_H5_YW)]; text(0.6,0.8*max(PSD_f_H5_YW_mean),str,'HorizontalAlignment','left');

subplot(5,2,2),plot(f_f_I1_YW, PSD_f_I1_YW_mean), title('PSD floatI1'), xlabel('Frequency (Hz)'); ylabel('bpm^2/Hz');
str= ['VLF=',num2str(VLF_f_I1_YW)]; text(0.6,0.8*max(PSD_f_I1_YW_mean),str,'HorizontalAlignment','left');
subplot(5,2,4),plot(f_f_I2_YW, PSD_f_I2_YW_mean), title('PSD floatI2'), xlabel('Frequency (Hz)'); ylabel('bpm^2/Hz');
str= ['VLF=',num2str(VLF_f_I2_YW)]; text(0.6,0.8*max(PSD_f_I2_YW_mean),str,'HorizontalAlignment','left');
subplot(5,2,6),plot(f_f_I3_YW, PSD_f_I3_YW_mean), title('PSD floatI3'), xlabel('Frequency (Hz)'); ylabel('bpm^2/Hz');
str= ['VLF=',num2str(VLF_f_I3_YW)]; text(0.6,0.8*max(PSD_f_I3_YW_mean),str,'HorizontalAlignment','left');
subplot(5,2,8),plot(f_f_I4_YW, PSD_f_I4_YW_mean), title('PSD floatI4'), xlabel('Frequency (Hz)'); ylabel('bpm^2/Hz');
str= ['VLF=',num2str(VLF_f_I4_YW)]; text(0.6,0.8*max(PSD_f_I4_YW_mean),str,'HorizontalAlignment','left');
subplot(5,2,10),plot(f_f_I5_YW, PSD_f_I5_YW_mean), title('PSD floatI5'), xlabel('Frequency (Hz)'); ylabel('bpm^2/Hz');
str= ['VLF=',num2str(VLF_f_I5_YW)]; text(0.6,0.8*max(PSD_f_I5_YW_mean),str,'HorizontalAlignment','left');

%Do table 
Row = {'Healthy';'IUGR'};
%Var = {'LF_mean_std','MF_mean_std)','HF_mean_std'};
F_H = [freq_range_H1; freq_range_H2; freq_range_H3; freq_range_H4; freq_range_H5];
F_I = [freq_range_I1; freq_range_I2; freq_range_I3; freq_range_I4; freq_range_I5];
VLF = [mean(F_H(:,1)) std(F_H(:,1)); mean(F_I(:,1)) std(F_I(:,1))];
LF = [mean(F_H(:,2)) std(F_H(:,2)); mean(F_I(:,2)) std(F_I(:,2))];
MF = [mean(F_H(:,3)) std(F_H(:,3)); mean(F_I(:,3)) std(F_I(:,3))];
HF = [mean(F_H(:,4)) std(F_H(:,4)); mean(F_I(:,4)) std(F_I(:,4))];
T = table(VLF, LF,MF,HF,'RowNames',Row);
T.Properties.VariableNames{'VLF'}= 'VLF_mean_std_bpm2';
T.Properties.VariableNames{'LF'}= 'LF_mean_std_bpm2';
T.Properties.VariableNames{'MF'}='MF_mean_std_bpm2';
T.Properties.VariableNames{'HF'}='HF_mean_std_bpm2';
T.Properties.VariableUnits = {'bpm^2' 'bpm^2' 'bpm^2' 'bpm^2'};

% % Get the table in string form.
% T = evalc('disp(T)');
% % Use TeX Markup for bold formatting and underscores.
% T= strrep(T,'<strong>','\bf');
% T= strrep(T,'</strong>','\rm');
% T = strrep(T,'_','\_');
% % Get a fixed-width font.
% FixedWidth = get(0,'FixedWidthFontName');
% % Output the table using the annotation command.
% figure
% annotation(gcf,'Textbox','String',T,'Interpreter','Tex',...
%     'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);
%___________________________________________________________________________________________________________-
%Counting accelerations 
% 
% figure
% plot(t_Healthy_1_inter(1:windo-1, :), float_H1(:,:), 'r');

%% Apen Calculation
%H1
c=1;
for i = 1:length(Healthy_1_inter(1,:))
    if isnan(Healthy_1_inter(:,i))
    else
           ApEn_H(c) = ApEn(1,0.1,Healthy_1_inter(:,i),0);
            c = c+1
    end
    
end
% figure;
% plot(ApEn_H)
% ylim([-1 2])
% title('Healthy1')
ApEn_H1=mean(ApEn_H);

%H2
c=1;
for i = 1:length(Healthy_2_inter(1,:))
    if isnan(Healthy_2_inter(:,i))
    else
        ApEn_H(c) = ApEn(1,0.2,Healthy_2_inter(:,i),0);
        c = c+1
    end
    
end
% figure;
% plot(ApEn_H)
% ylim([-1 2])
% title('Healthy2')
ApEn_H2=mean(ApEn_H);


%H3
c=1;
for i = 1:length(Healthy_3_inter(1,:))
    if isnan(Healthy_3_inter(:,i))
    else
        ApEn_H(c) = ApEn(1,0.2,Healthy_3_inter(:,i),0);
        c = c+1
    end
    
end
% figure;
% plot(ApEn_H)
% ylim([-1 2])
% title('Healthy3')
ApEn_H3=mean(ApEn_H);

%H4
c=1;
for i = 1:length(Healthy_4_inter(1,:))
    if isnan(Healthy_4_inter(:,i))
    else
        ApEn_H(c) = ApEn(1,0.2,Healthy_4_inter(:,i),0);
        c = c+1
    end
    
end
% figure;
% plot(ApEn_H)
% ylim([-1 2])
% title('Healthy4')
ApEn_H4=mean(ApEn_H);

%H5
c=1;
for i = 1:length(Healthy_5_inter(1,:))
    if isnan(Healthy_5_inter(:,i))
    else
        ApEn_H(c) = ApEn(1,0.2,Healthy_5_inter(:,i),0);
        c = c+1
    end
    
end
% figure;
% plot(ApEn_H)
% ylim([-1 2])
% title('Healthy5')
ApEn_H5=mean(ApEn_H);

%I1
c=1;
for i = 1:length(IUGR_1_inter(1,:))
    if isnan(IUGR_1_inter(:,i))
    else
        ApEn_H(c) = ApEn(1,0.2,IUGR_1_inter(:,i),0);
        c = c+1
    end
    
end
% figure;
% plot(ApEn_H)
% ylim([-1 2])
% title('IUGR1')
ApEn_I1=mean(ApEn_H);


%I2
c=1;
for i = 1:length(IUGR_2_inter(1,:))
    if isnan(IUGR_2_inter(:,i))
    else
        ApEn_H(c) = ApEn(1,0.2,IUGR_2_inter(:,i),0);
        c = c+1
    end
    
end
% figure;
% plot(ApEn_H)
% ylim([-1 2])
% title('IUGR2')
ApEn_I2=mean(ApEn_H);

%I3
c=1;
for i = 1:length(IUGR_3_inter(1,:))
    if isnan(IUGR_3_inter(:,i))
    else
        ApEn_H(c) = ApEn(1,0.2,IUGR_3_inter(:,i),0);
        c = c+1
    end
    
end
% figure;
% plot(ApEn_H)
% ylim([-1 2])
% title('IUGR3')
ApEn_I3=mean(ApEn_H);

%I4
c=1;
for i = 1:length(IUGR_4_inter(1,:))
    if isnan(IUGR_4_inter(:,i))
    else
        ApEn_H(c) = ApEn(1,0.2,IUGR_4_inter(:,i),0);
        c = c+1
    end
    
end
% figure;
% plot(ApEn_H)
% ylim([-1 2])
% title('IUGR4')
ApEn_I4=mean(ApEn_H);


%I5
c=1;
for i = 1:length(IUGR_5_inter(1,:))
    if isnan(IUGR_5_inter(:,i))
    else
        ApEn_H(c) = ApEn(1,0.2,IUGR_5_inter(:,i),0);
        c = c+1
    end
    
end
% figure;
% plot(ApEn_H)
% ylim([-1 2])
% title('IUGR5')
ApEn_I5=mean(ApEn_H);

Names_H={'Healhy_1';'Healthy_2';'Healhy_3';'Healhy_4'; 'Healhy_5'; 'IUGR_1'; 'IUGR_2'; 'IUGR_3'; 'IUGR_4'; 'IUGR_5'};
ApEn_values=[ ApEn_H1; ApEn_H2; ApEn_H3; ApEn_H4; ApEn_H5; ApEn_I1; ApEn_I2; ApEn_I3; ApEn_I4; ApEn_I5];
%IUGR={'IUGR_1';'IUGR_2';'IUGR_3';'IUGR_4'; 'IUGR_5'};
T_ApEn=table(ApEn_values,'RowNames', Names_H);

% Get the table in string form.
TString = evalc('disp(T)');
% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
figure
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
    'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

ApEn_mean_std=[mean(ApEn_values(1:5)) std(ApEn_values(1:5));  mean(ApEn_values(6:10)) std(ApEn_values(6:10))];
T_ApEn_mean=table(ApEn_mean_std, 'RowNames', Row);

%% Pointcare Plot Healthy
RRH1 = 1000.* (1./(Healthy_1_corr/60));
RRH2 = 1000.* (1./(Healthy_2_corr/60));
RRH3 = 1000.* (1./(Healthy_3_corr/60));
RRH4 = 1000.* (1./(Healthy_4_corr/60));
RRH5 = 1000.* (1./(Healthy_5_corr/60));
co=cos(pi/4);
si=sin(pi/4);
alpha= linspace(0,2*pi,100);
figure
subplot(3,2,1), plot(RRH1(1:end-1),RRH1(2:end),'.'), xlabel('RR(n) [ms]'),
xlim([300 600]),ylim([300 600]),
ylabel('RR(n+1) [ms]'),title('Pointcare Plot Healthy_1');
xp= RRH1(1:end-1);
xp(end)=[];
xm= RRH1(2:end);
xm(1)=[];
PC_SD1= std(xp-xm, 'omitnan')/sqrt(2);
PC_SD2= std(xp+xm, 'omitnan')/sqrt(2);
xpcenter=mean(xp,'omitnan');
xmcenter=mean(xm,'omitnan');
hold on 
plot(si*PC_SD2*cos(alpha)+co*PC_SD1*sin(alpha)+xpcenter,co*PC_SD2*cos(alpha)-si*PC_SD1*sin(alpha)+xmcenter,xpcenter,xmcenter,'*');
hold on
plot([xpcenter xpcenter+PC_SD1/sqrt(2)],[xmcenter xmcenter-PC_SD1/sqrt(2)]);
hold on
plot([xpcenter xpcenter+PC_SD2/sqrt(2)],[xmcenter xmcenter+PC_SD2/sqrt(2)]);
PC_SD1_H1 = PC_SD1;
PC_SD2_H1 = PC_SD2;

subplot(3,2,2), plot(RRH2(1:end-1),RRH2(2:end),'.'), xlabel('RR(n) [ms]'),
xlim([300 600]),ylim([300 600]),
ylabel('RR(n+1) [ms]'),title('Pointcare Plot Healthy_2');
xp= RRH2(1:end-1);
xp(end)=[];
xm= RRH2(2:end);
xm(1)=[];
PC_SD1= std(xp-xm, 'omitnan')/sqrt(2);
PC_SD2= std(xp+xm, 'omitnan')/sqrt(2);
xpcenter=mean(xp,'omitnan');
xmcenter=mean(xm,'omitnan');
hold on 
plot(si*PC_SD2*cos(alpha)+co*PC_SD1*sin(alpha)+xpcenter,co*PC_SD2*cos(alpha)-si*PC_SD1*sin(alpha)+xmcenter,xpcenter,xmcenter,'*');
hold on
plot([xpcenter xpcenter+PC_SD1/sqrt(2)],[xmcenter xmcenter-PC_SD1/sqrt(2)]);
hold on
plot([xpcenter xpcenter+PC_SD2/sqrt(2)],[xmcenter xmcenter+PC_SD2/sqrt(2)]);
PC_SD1_H2 = PC_SD1;
PC_SD2_H2 = PC_SD2;


subplot(3,2,3), plot(RRH3(1:end-1),RRH3(2:end),'.'), xlabel('RR(n) [ms]'),
xlim([300 600]),ylim([300 600]),
ylabel('RR(n+1) [ms]'),title('Pointcare Plot Healthy_3');
xp= RRH3(1:end-1);
xp(end)=[];
xm= RRH3(2:end);
xm(1)=[];
PC_SD1= std(xp-xm, 'omitnan')/sqrt(2);
PC_SD2= std(xp+xm, 'omitnan')/sqrt(2);
xpcenter=mean(xp,'omitnan');
xmcenter=mean(xm,'omitnan');
hold on 
plot(si*PC_SD2*cos(alpha)+co*PC_SD1*sin(alpha)+xpcenter,co*PC_SD2*cos(alpha)-si*PC_SD1*sin(alpha)+xmcenter,xpcenter,xmcenter,'*');
hold on
plot([xpcenter xpcenter+PC_SD1/sqrt(2)],[xmcenter xmcenter-PC_SD1/sqrt(2)]);
hold on
plot([xpcenter xpcenter+PC_SD2/sqrt(2)],[xmcenter xmcenter+PC_SD2/sqrt(2)]);
PC_SD1_H3 = PC_SD1;
PC_SD2_H3 = PC_SD2;


subplot(3,2,4), plot(RRH4(1:end-1),RRH4(2:end),'.'), xlabel('RR(n) [ms]'),
xlim([300 600]),ylim([300 600]),
ylabel('RR(n+1) [ms]'),title('Pointcare Plot Healthy_4');
xp= RRH4(1:end-1);
xp(end)=[];
xm= RRH4(2:end);
xm(1)=[];
PC_SD1= std(xp-xm, 'omitnan')/sqrt(2);
PC_SD2= std(xp+xm, 'omitnan')/sqrt(2);
xpcenter=mean(xp,'omitnan');
xmcenter=mean(xm,'omitnan');
hold on 
plot(si*PC_SD2*cos(alpha)+co*PC_SD1*sin(alpha)+xpcenter,co*PC_SD2*cos(alpha)-si*PC_SD1*sin(alpha)+xmcenter,xpcenter,xmcenter,'*');
hold on
plot([xpcenter xpcenter+PC_SD1/sqrt(2)],[xmcenter xmcenter-PC_SD1/sqrt(2)]);
hold on
plot([xpcenter xpcenter+PC_SD2/sqrt(2)],[xmcenter xmcenter+PC_SD2/sqrt(2)]);
PC_SD1_H4 = PC_SD1;
PC_SD2_H4 = PC_SD2;


subplot(3,2,5), plot(RRH5(1:end-1),RRH5(2:end),'.'), xlabel('RR(n) [ms]'),
xlim([300 600]),ylim([300 600]),
ylabel('RR(n+1) [ms]'),title('Pointcare Plot Healthy_5');
xp= RRH5(1:end-1);
xp(end)=[];
xm= RRH5(2:end);
xm(1)=[];
PC_SD1= std(xp-xm, 'omitnan')/sqrt(2);
PC_SD2= std(xp+xm, 'omitnan')/sqrt(2);
xpcenter=mean(xp,'omitnan');
xmcenter=mean(xm,'omitnan');
hold on 
plot(si*PC_SD2*cos(alpha)+co*PC_SD1*sin(alpha)+xpcenter,co*PC_SD2*cos(alpha)-si*PC_SD1*sin(alpha)+xmcenter,xpcenter,xmcenter,'*');
hold on
plot([xpcenter xpcenter+PC_SD1/sqrt(2)],[xmcenter xmcenter-PC_SD1/sqrt(2)]);
hold on
plot([xpcenter xpcenter+PC_SD2/sqrt(2)],[xmcenter xmcenter+PC_SD2/sqrt(2)]);
PC_SD1_H5 = PC_SD1;
PC_SD2_H5 = PC_SD2;


%% Pointcare Plot IUGR
RRIUGR1 = 1000.* (1./(IUGR_1_corr/60));
RRIUGR2 = 1000.* (1./(IUGR_2_corr/60));
RRIUGR3 = 1000.* (1./(IUGR_3_corr/60));
RRIUGR4 = 1000.* (1./(IUGR_4_corr/60));
RRIUGR5 = 1000.* (1./(IUGR_5_corr/60));

figure
subplot(3,2,1), plot(RRIUGR1(1:end-1),RRIUGR1(2:end),'.'), xlabel('RR(n) [ms]'),
xlim([300 600]),ylim([300 600]),
ylabel('RR(n+1) [ms]'),title('Pointcare Plot IUGR_1');
xp= RRIUGR1(1:end-1);
xp(end)=[];
xm= RRIUGR1(2:end);
xm(1)=[];
PC_SD1= std(xp-xm, 'omitnan')/sqrt(2);
PC_SD2= std(xp+xm, 'omitnan')/sqrt(2);
xpcenter=mean(xp,'omitnan');
xmcenter=mean(xm,'omitnan');
hold on 
plot(si*PC_SD2*cos(alpha)+co*PC_SD1*sin(alpha)+xpcenter,co*PC_SD2*cos(alpha)-si*PC_SD1*sin(alpha)+xmcenter,xpcenter,xmcenter,'*');
hold on
plot([xpcenter xpcenter+PC_SD1/sqrt(2)],[xmcenter xmcenter-PC_SD1/sqrt(2)]);
hold on
plot([xpcenter xpcenter+PC_SD2/sqrt(2)],[xmcenter xmcenter+PC_SD2/sqrt(2)]);
PC_SD1_I1 = PC_SD1;
PC_SD2_I1 = PC_SD2;

subplot(3,2,2), plot(RRIUGR2(1:end-1),RRIUGR2(2:end),'.'), xlabel('RR(n) [ms]'),
xlim([300 600]),ylim([300 600]),
ylabel('RR(n+1) [ms]'),title('Pointcare Plot IUGR_2');
xp= RRIUGR2(1:end-1);
xp(end)=[];
xm= RRIUGR2(2:end);
xm(1)=[];
PC_SD1= std(xp-xm, 'omitnan')/sqrt(2);
PC_SD2= std(xp+xm, 'omitnan')/sqrt(2);
xpcenter=mean(xp,'omitnan');
xmcenter=mean(xm,'omitnan');
hold on 
plot(si*PC_SD2*cos(alpha)+co*PC_SD1*sin(alpha)+xpcenter,co*PC_SD2*cos(alpha)-si*PC_SD1*sin(alpha)+xmcenter,xpcenter,xmcenter,'*');
hold on
plot([xpcenter xpcenter+PC_SD1/sqrt(2)],[xmcenter xmcenter-PC_SD1/sqrt(2)]);
hold on
plot([xpcenter xpcenter+PC_SD2/sqrt(2)],[xmcenter xmcenter+PC_SD2/sqrt(2)]);
PC_SD1_I2 = PC_SD1;
PC_SD2_I2 = PC_SD2;

subplot(3,2,3), plot(RRIUGR3(1:end-1),RRIUGR3(2:end),'.'), xlabel('RR(n) [ms]'),
xlim([300 600]),ylim([300 600]),
ylabel('RR(n+1) [ms]'),title('Pointcare Plot IUGR_3');
xp= RRIUGR3(1:end-1);
xp(end)=[];
xm= RRIUGR3(2:end);
xm(1)=[];
PC_SD1= std(xp-xm, 'omitnan')/sqrt(2);
PC_SD2= std(xp+xm, 'omitnan')/sqrt(2);
xpcenter=mean(xp,'omitnan');
xmcenter=mean(xm,'omitnan');
hold on 
plot(si*PC_SD2*cos(alpha)+co*PC_SD1*sin(alpha)+xpcenter,co*PC_SD2*cos(alpha)-si*PC_SD1*sin(alpha)+xmcenter,xpcenter,xmcenter,'*');
hold on
plot([xpcenter xpcenter+PC_SD1/sqrt(2)],[xmcenter xmcenter-PC_SD1/sqrt(2)]);
hold on
plot([xpcenter xpcenter+PC_SD2/sqrt(2)],[xmcenter xmcenter+PC_SD2/sqrt(2)]);
PC_SD1_I3 = PC_SD1;
PC_SD2_I3 = PC_SD2;

subplot(3,2,4), plot(RRIUGR4(1:end-1),RRIUGR4(2:end),'.'), xlabel('RR(n) [ms]'),
xlim([300 600]),ylim([300 600]),
ylabel('RR(n+1) [ms]'),title('Pointcare Plot IUGR_4');
xp= RRIUGR4(1:end-1);
xp(end)=[];
xm= RRIUGR4(2:end);
xm(1)=[];
PC_SD1= std(xp-xm, 'omitnan')/sqrt(2);
PC_SD2= std(xp+xm, 'omitnan')/sqrt(2);
xpcenter=mean(xp,'omitnan');
xmcenter=mean(xm,'omitnan');
hold on 
plot(si*PC_SD2*cos(alpha)+co*PC_SD1*sin(alpha)+xpcenter,co*PC_SD2*cos(alpha)-si*PC_SD1*sin(alpha)+xmcenter,xpcenter,xmcenter,'*');
hold on
plot([xpcenter xpcenter+PC_SD1/sqrt(2)],[xmcenter xmcenter-PC_SD1/sqrt(2)]);
hold on
plot([xpcenter xpcenter+PC_SD2/sqrt(2)],[xmcenter xmcenter+PC_SD2/sqrt(2)]);
PC_SD1_I4 = PC_SD1;
PC_SD2_I4 = PC_SD2;

subplot(3,2,5), plot(RRIUGR5(1:end-1),RRIUGR5(2:end),'.'), xlabel('RR(n) [ms]'),
xlim([300 600]),ylim([300 600]),
ylabel('RR(n+1) [ms]'),title('Pointcare Plot IUGR_5');
xp= RRIUGR5(1:end-1);
xp(end)=[];
xm= RRIUGR5(2:end);
xm(1)=[];
PC_SD1= std(xp-xm, 'omitnan')/sqrt(2);
PC_SD2= std(xp+xm, 'omitnan')/sqrt(2);
xpcenter=mean(xp,'omitnan');
xmcenter=mean(xm,'omitnan');
hold on 
plot(si*PC_SD2*cos(alpha)+co*PC_SD1*sin(alpha)+xpcenter,co*PC_SD2*cos(alpha)-si*PC_SD1*sin(alpha)+xmcenter,xpcenter,xmcenter,'*');
hold on
plot([xpcenter xpcenter+PC_SD1/sqrt(2)],[xmcenter xmcenter-PC_SD1/sqrt(2)]);
hold on
plot([xpcenter xpcenter+PC_SD2/sqrt(2)],[xmcenter xmcenter+PC_SD2/sqrt(2)]);
PC_SD1_I5 = PC_SD1;
PC_SD2_I5 = PC_SD2;

%% Calculate baseline accordind to Mantel Healthy 1
%Plot RR intervals in time
RRH1 = 1000.* (1./(Healthy_1_corr/60));
%figure; plot(RRH1); figure; plot(t_Healthy_1,RRH1);

% Calculate the mean av_1 of 2.5s window
RRH1B = RRH1;
time = 5; %2.5sec 
av_1 = zeros(1,floor(length(RRH1)/time));
for i=1:floor(length(RRH1)/time)
        av_1(i)= mean(RRH1((i-1)*time+1:i*time), 'omitnan');
end
t_av = 0:time:length(RRH1) - (length(RRH1) - time*floor(length(RRH1)/time))-1;
%figure;plot(t_av,av_1);hold on; plot(RRH1);

% frequency distribution of RR intervals with 1ms BW
[n,edges] = histcounts(av_1,'BinWidth',1);
figure;histogram(av_1,'BinWidth',1), xlabel('bpm'), ylabel('Counts');

% Calculating the peak P as define in Mantel algorithm for the next steop
for i = length(n)-1:-1:2
    if n(i)>n(i+1) & n(i)>=n(i-1) & n(i)>=n(i-2) & n(i)>=n(i-3) & n(i)>=n(i-4) & n(i)>=n(i-5) & sum(n(i:end)) > sum(n)/8
        PH1 = edges(i)
        break
    else 
        continue
    end
end
% 1st Filtering step of RRH1B 
RRH1B(1)=PH1;
for i = 1:length(RRH1B)
    if abs(RRH1B(i) - PH1)<=60
        RRH1B(1)= 0.95*RRH1B(1)+0.05*RRH1B(i);
    end
end
for i = 2:length(RRH1B)
    if abs(RRH1B(i) - PH1)<=60
        RRH1B(i)= 0.95*RRH1B(i-1)+0.05*RRH1B(i);
    else
        RRH1B(i) = RRH1B(i-1);
    end
end
for i = length(RRH1B)-1:1
    RRH1B(i)= 0.95*RRH1B(i+1)+0.05*RRH1B(i);
end
%Ploting the resulting RRH1B baseline in RR intervals and in FHR1
%figure; plot(RRH1B);hold on; plot(RRH1);title('RRH1B 1st filtering');
FHR1 = 1000.*60.*(1./RRH1B);
% 1st Trimming at U = L = 20 bpm,comparison between Health1 and FHR1, change
%in RRH1B and FHR1
U = 20;
L=20;
%figure;plot(FHR1+L);hold on;plot(FHR1-U); 
%hold on; plot(Healthy_1); title('FHR1 before 1st trimming');
for i = 1:length(RRH1)
    if  FHR1(i) - L <= Healthy_1_corr(i) & Healthy_1_corr(i)<= FHR1(i) +U
        FHR1(i)=Healthy_1_corr(i);
        RRH1B(i)=RRH1(i);
    end
end
%figure;plot(FHR1);hold on; plot(Healthy_1_corr); title('FHR1 after 1st trimming');

% 2nd filtering
RRH1B(1)=PH1;
for i = 1:length(RRH1B)
    if abs(RRH1B(i) - PH1)<=60
        RRH1B(1)= 0.95*RRH1B(1)+0.05*RRH1B(i);
    end
end
for i = 2:length(RRH1B)
    if abs(RRH1B(i) - PH1)<=60
        RRH1B(i)= 0.95*RRH1B(i-1)+0.05*RRH1B(i);
    else
        RRH1B(i) = RRH1B(i-1);
    end
end
for i = length(RRH1B)-1:1
    RRH1B(i)= 0.95*RRH1B(i+1)+0.05*RRH1B(i);
end

FHR1 = 1000.*60.*(1./RRH1B);
%figure;plot(FHR1);hold on;plot(Healthy_1_corr);title('FHR1 2nd filtering');

% 2nd Trimming
U = 15;
L=15;
%figure;plot(FHR1+L);hold on;plot(FHR1-U)
%hold on; plot(Healthy_1); title('FHR1 before 2snd trimming');
for i = 1:length(RRH1)
    if  FHR1(i)- L <= Healthy_1_corr(i) & Healthy_1_corr(i)<= FHR1(i) + L
        FHR1(i)=Healthy_1_corr(i);
        RRH1B(i)=RRH1(i);
    end
end
%figure;plot(FHR1);hold on; plot(Healthy_1_corr);title('FHR1 after 2nd trimming');

% 3rd filtering
RRH1B(1)=PH1;
for i = 1:length(RRH1B)
    if abs(RRH1B(i) - PH1)<=60
        RRH1B(1)= 0.95*RRH1B(1)+0.05*RRH1B(i);
    end
end
for i = 2:length(RRH1B)
    if abs(RRH1B(i) - PH1)<=60
        RRH1B(i)= 0.95*RRH1B(i-1)+0.05*RRH1B(i);
    else
        RRH1B(i) = RRH1B(i-1);
    end
end
for i = length(RRH1B)-1:1
    RRH1B(i)= 0.95*RRH1B(i+1)+0.05*RRH1B(i);
end

FHR1 = 1000.*60.*(1./RRH1B);
%figure;plot(FHR1);hold on;plot(Healthy_1_corr);title('FHR1 3dr filtering');

% 3nd Trimming
U = 10;
L=10;
%figure;plot(FHR1+L);hold on;plot(FHR1-U)
%hold on; plot(Healthy_1); title('FHR1 before 3rd trimming');
for i = 1:length(RRH1)
    if  FHR1(i)- L <= Healthy_1_corr(i) & Healthy_1_corr(i)<= FHR1(i) + L
        FHR1(i)=Healthy_1_corr(i);
        RRH1B(i)=RRH1(i);
    end
end
%figure;plot(FHR1);hold on; plot(Healthy_1_corr);title('FHR1 after 3rd trimming');

% 4rd filtering
RRH1B(1)=PH1;
for i = 1:length(RRH1B)
    if abs(RRH1B(i) - PH1)<=60
        RRH1B(1)= 0.95*RRH1B(1)+0.05*RRH1B(i);
    end
end
for i = 2:length(RRH1B)
    if abs(RRH1B(i) - PH1)<=60
        RRH1B(i)= 0.95*RRH1B(i-1)+0.05*RRH1B(i);
    else
        RRH1B(i) = RRH1B(i-1);
    end
end
for i = length(RRH1B)-1:1
    RRH1B(i)= 0.95*RRH1B(i+1)+0.05*RRH1B(i);
end

FHR1 = 1000.*60.*(1./RRH1B);
%figure;plot(FHR1);hold on;plot(Healthy_1_corr);title('FHR1 4th filtering');

% 4nd Trimming
U = 5;
L=5;
%figure;plot(FHR1+L);hold on;plot(FHR1-U)
%hold on; plot(Healthy_1); title('FHR1 before 4th trimming');
for i = 1:length(RRH1)
    if  FHR1(i)- L <= Healthy_1_corr(i) & Healthy_1_corr(i)<= FHR1(i) + L
        FHR1(i)=Healthy_1_corr(i);
        RRH1B(i)=RRH1(i);
    end
end
%figure;plot(FHR1);hold on; plot(Healthy_1_corr);title('FHR1 after 4th trimming');

% 5th filtering
RRH1B(1)=PH1;
for i = 1:length(RRH1B)
    if abs(RRH1B(i) - PH1)<=60
        RRH1B(1)= 0.95*RRH1B(1)+0.05*RRH1B(i);
    end
end
for i = 2:length(RRH1B)
    if abs(RRH1B(i) - PH1)<=60
        RRH1B(i)= 0.95*RRH1B(i-1)+0.05*RRH1B(i);
    else
        RRH1B(i) = RRH1B(i-1);
    end
end
for i = length(RRH1B)-1:1
    RRH1B(i)= 0.95*RRH1B(i+1)+0.05*RRH1B(i);
end

FHR1 = 1000.*60.*(1./RRH1B);
figure;plot(FHR1);hold on; plot(Healthy_1_corr); title('FHR1 and Baseline'),legend('Baseline Healthy 1','Healthy 1');


%% From FHR to RR intervals , calculate baseline accordind to Mantel Healthy 2
%Plot RR intervals in time
RRH2 = 1000.* (1./(Healthy_2_corr/60));
%figure; plot(RRH2);

% Calculate the mean av_1 of 2.5s window
RRH2B = RRH2;
time = 5; %2.5sec 
av_1 = zeros(1,floor(length(RRH2)/time));
for i=1:floor(length(RRH2)/time)
        av_2(i)= mean(RRH2((i-1)*time+1:i*time));
end
t_av = 0:time:length(RRH2) - (length(RRH2) - time*floor(length(RRH2)/time))-1;
%figure;plot(t_av,av_2); figure; plot(RRH2B);


% frequency distribution of RR intervals with 1ms BW
[n,edges] = histcounts(av_2,'BinWidth',1);
% figure;histogram(av_2,'BinWidth',1);

% Calculating the peak P as define in Mantel algorithm for the next steop
for i = length(n)-1:-1:2
    
    if n(i+1) < n(i) & n(i)>=n(i-1) & n(i)>=n(i-2) & n(i)>=n(i-3) & n(i)>=n(i-4) & n(i)>=n(i-5) & sum(n(i:end)) > sum(n)/8
        PH2 = edges(i)
        break
    else 
        continue
    end
end

% 1st Filtering step of RRH1B
RRH2B(1)=PH2;
i=1;
while  i <= length(RRH2B)
    if (~isnan((RRH2B(i)))) & abs(RRH2B(i) - PH2)<=60
        RRH2B(1)= 0.95*RRH2B(1)+0.05*RRH2B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRH2B)
    if (~isnan((RRH2B(i))))
        if  abs(RRH2B(i) - PH2)<=60
            if isnan((RRH2B(i-1)))
                RRH2B(i)=RRH2B(i);
            else
                RRH2B(i)= 0.95*RRH2B(i-1)+0.05*RRH2B(i);
            end
        else
            if isnan((RRH2B(i-1)))
                RRH2B(i)=RRH2B(i);
            else
                RRH2B(i)= RRH2B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRH2B)-1;
while i >= 1
    if (~isnan((RRH2B(i))))
        if isnan((RRH2B(i+1)))
            RRH2B(i)= RRH2B(i);
        else
            RRH2B(i)= 0.95*RRH2B(i+1)+0.05*RRH2B(i);
        end
    end
    i = i-1;
end

FHR2 = 1000.*60.*(1./RRH2B);
%figure;plot(FHR2);hold on; plot(Healthy_2_corr); title('FHR2 1st filtering');
%figure;plot(RRH2B);

% 1st Trimming at U = L = 20 bpm,comparison between Health1 and FHR1, change
%in RRH1B and FHR1
U = 20;
L=20;
%figure;plot(FHR2+L);hold on;plot(FHR2-U); 
%hold on; plot(Healthy_2_corr); title('FHR1 before 1st trimming');
for i = 1:length(RRH2)
    if  FHR2(i) - L <= Healthy_2_corr(i) & Healthy_2_corr(i)<= FHR2(i) +U
        FHR2(i)=Healthy_2_corr(i);
        RRH2B(i)=RRH2(i);
    end
end
%figure;plot(FHR2);hold on; plot(Healthy_2_corr); title('FHR1 after 1st trimming');

% 2nd filtering
RRH2B(1)=PH2;
i=1;
while  i <= length(RRH2B)
    if (~isnan((RRH2B(i)))) & abs(RRH2B(i) - PH2)<=60
        RRH2B(1)= 0.95*RRH2B(1)+0.05*RRH2B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRH2B)
    if (~isnan((RRH2B(i))))
        if  abs(RRH2B(i) - PH2)<=60
            if isnan((RRH2B(i-1)))
                RRH2B(i)=RRH2B(i);
            else
                RRH2B(i)= 0.95*RRH2B(i-1)+0.05*RRH2B(i);
            end
        else
            if isnan((RRH2B(i-1)))
                RRH2B(i)=RRH2B(i);
            else
                RRH2B(i)= RRH2B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRH2B)-1;
while i >= 1
    if (~isnan((RRH2B(i))))
        if isnan((RRH2B(i+1)))
            RRH2B(i)= RRH2B(i);
        else
            RRH2B(i)= 0.95*RRH2B(i+1)+0.05*RRH2B(i);
        end
    end
    i = i-1;
end

FHR2 = 1000.*60.*(1./RRH2B);
%figure;plot(FHR2);hold on; plot(Healthy_2_corr); title('FHR2 2nd filtering');

% 2nd Trimming
U = 15;
L=15;
%figure;plot(FHR2+L);hold on;plot(FHR2-U)
%hold on; plot(Healthy_2_corr); title('FHR2 before 2nd trimming');
for i = 1:length(RRH2)
    if  FHR2(i) - L <= Healthy_2_corr(i) & Healthy_2_corr(i)<= FHR2(i) +U
        FHR2(i)=Healthy_2_corr(i);
        RRH2B(i)=RRH2(i);
    end
end
%figure;plot(FHR2);hold on; plot(Healthy_2_corr); title('FHR1 after 2nd trimming');

% 3rd filtering
RRH2B(1)=PH2;
i=1;
while  i <= length(RRH2B)
    if (~isnan((RRH2B(i)))) & abs(RRH2B(i) - PH2)<=60
        RRH2B(1)= 0.95*RRH2B(1)+0.05*RRH2B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRH2B)
    if (~isnan((RRH2B(i))))
        if  abs(RRH2B(i) - PH2)<=60
            if isnan((RRH2B(i-1)))
                RRH2B(i)=RRH2B(i);
            else
                RRH2B(i)= 0.95*RRH2B(i-1)+0.05*RRH2B(i);
            end
        else
            if isnan((RRH2B(i-1)))
                RRH2B(i)=RRH2B(i);
            else
                RRH2B(i)= RRH2B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRH2B)-1;
while i >= 1
    if (~isnan((RRH2B(i))))
        if isnan((RRH2B(i+1)))
            RRH2B(i)= RRH2B(i);
        else
            RRH2B(i)= 0.95*RRH2B(i+1)+0.05*RRH2B(i);
        end
    end
    i = i-1;
end

FHR2 = 1000.*60.*(1./RRH2B);
%figure;plot(FHR2);hold on; plot(Healthy_2_corr); title('FHR2 3rd filtering');

% 3nd Trimming
U = 10;
L=10;
%figure;plot(FHR2+L);hold on;plot(FHR2-U)
%hold on; plot(Healthy_2_corr); title('FHR2 before 2nd trimming');
for i = 1:length(RRH2)
    if  FHR2(i) - L <= Healthy_2_corr(i) & Healthy_2_corr(i)<= FHR2(i) +U
        FHR2(i)=Healthy_2_corr(i);
        RRH2B(i)=RRH2(i);
    end
end
%figure;plot(FHR2);hold on; plot(Healthy_2_corr); title('FHR2 after 3rd trimming');

% 4rd filtering
RRH2B(1)=PH2;
i=1;
while  i <= length(RRH2B)
    if (~isnan((RRH2B(i)))) & abs(RRH2B(i) - PH2)<=60
        RRH2B(1)= 0.95*RRH2B(1)+0.05*RRH2B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRH2B)
    if (~isnan((RRH2B(i))))
        if  abs(RRH2B(i) - PH2)<=60
            if isnan((RRH2B(i-1)))
                RRH2B(i)=RRH2B(i);
            else
                RRH2B(i)= 0.95*RRH2B(i-1)+0.05*RRH2B(i);
            end
        else
            if isnan((RRH2B(i-1)))
                RRH2B(i)=RRH2B(i);
            else
                RRH2B(i)= RRH2B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRH2B)-1;
while i >= 1
    if (~isnan((RRH2B(i))))
        if isnan((RRH2B(i+1)))
            RRH2B(i)= RRH2B(i);
        else
            RRH2B(i)= 0.95*RRH2B(i+1)+0.05*RRH2B(i);
        end
    end
    i = i-1;
end

FHR2 = 1000.*60.*(1./RRH2B);
%figure;plot(FHR2);hold on; plot(Healthy_2_corr); title('FHR2 4th filtering');

% 4nd Trimming of RHB2
U = 5;
L=5;

%figure;plot(FHR2+L);hold on;plot(FHR2-U)
%hold on; plot(Healthy_2_corr); title('FHR2 before 2nd trimming');
for i = 1:length(RRH2)
    if  FHR2(i) - L <= Healthy_2_corr(i) & Healthy_2_corr(i)<= FHR2(i) +U
        FHR2(i)=Healthy_2_corr(i);
        RRH2B(i)=RRH2(i);
    end
end
%figure;plot(FHR2);hold on; plot(Healthy_2_corr); title('FHR2 after 4th trimming');


% 5th filtering

RRH2B(1)=PH2;
i=1;
while  i <= length(RRH2B)
    if (~isnan((RRH2B(i)))) & abs(RRH2B(i) - PH2)<=60
        RRH2B(1)= 0.95*RRH2B(1)+0.05*RRH2B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRH2B)
    if (~isnan((RRH2B(i))))
        if  abs(RRH2B(i) - PH2)<=60
            if isnan((RRH2B(i-1)))
                RRH2B(i)=RRH2B(i);
            else
                RRH2B(i)= 0.95*RRH2B(i-1)+0.05*RRH2B(i);
            end
        else
            if isnan((RRH2B(i-1)))
                RRH2B(i)=RRH2B(i);
            else
                RRH2B(i)= RRH2B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRH2B)-1;
while i >= 1
    if (~isnan((RRH2B(i))))
        if isnan((RRH2B(i+1)))
            RRH2B(i)= RRH2B(i);
        else
            RRH2B(i)= 0.95*RRH2B(i+1)+0.05*RRH2B(i);
        end
    end
    i = i-1;
end

FHR2 = 1000.*60.*(1./RRH2B);
figure;plot(FHR2);hold on; plot(Healthy_2_corr); title('Baseline of Healthy 2'),legend('Baseline','Healthy 2');

%% Calculate baseline accordind to Mantel Healthy 3
%Plot RR intervals in time
RRH3 = 1000.* (1./(Healthy_3_corr/60));
%figure; plot(RRH3); figure; plot(t_Healthy_3,RRH3);

% Calculate the mean av_3 of 2.5s window
RRH3B = RRH3;
time = 5; %2.5sec 
av_3 = zeros(1,floor(length(RRH3)/time));
for i=1:floor(length(RRH3)/time)
        av_3(i)= mean(RRH3((i-1)*time+1:i*time), 'omitnan');
end
t_av = 0:time:length(RRH3) - (length(RRH3) - time*floor(length(RRH3)/time))-1;
%figure;plot(t_av,av_3);hold on; plot(RRH3);

% frequency distribution of RR intervals with 1ms BW
[n,edges] = histcounts(av_3,'BinWidth',1);
% figure;histogram(av_3,'BinWidth',1);

% Calculating the peak P as define in Mantel algorithm for the next steop
for i = length(n)-1:-1:2

    if n(i+1) < n(i) & n(i)>=n(i-1) & n(i)>=n(i-2) & n(i)>=n(i-3) & n(i)>=n(i-4) & n(i)>=n(i-5) & sum(n(i:end)) > sum(n)/5
        PH3 = edges(i)
        break
    else 
        continue
    end
end

% 1st Filtering step of RRH3B
RRH3B(1)=PH3;
for i = 1:length(RRH3B)
    if abs(RRH3B(i) - PH3)<=60
        RRH3B(1)= 0.95*RRH3B(1)+0.05*RRH3B(i);
    end
end
for i = 2:length(RRH3B)
    if abs(RRH3B(i) - PH3)<=60
        RRH3B(i)= 0.95*RRH3B(i-1)+0.05*RRH3B(i);
    else
        RRH3B(i) = RRH3B(i-1);
    end
end
for i = length(RRH3B)-1:1
    RRH3B(i)= 0.95*RRH3B(i+1)+0.05*RRH3B(i);
end
%Ploting the resulting RRH3B baseline in RR intervals and in FHR3
%figure; plot(RRH3B);hold on; plot(RRH3);title('RRH3B 1st filtering');
FHR3 = 1000.*60.*(1./RRH3B);

% 1st Trimming at U = L = 20 bpm,comparison between Health3 and FHR3, change
%in RRH3B and FHR3
U = 20;
L=20;
%figure;plot(FHR3+L);hold on;plot(FHR3-U); 
%hold on; plot(Healthy_3); title('FHR3 before 3st trimming');
for i = 1:length(RRH3)
    if  FHR3(i) - L <= Healthy_3_corr(i) & Healthy_3_corr(i)<= FHR3(i) +U
        FHR3(i)=Healthy_3_corr(i);
        RRH3B(i)=RRH3(i);
    end
end
%figure;plot(FHR3);hold on; plot(Healthy_3_corr); title('FHR3 after 1st trimming');

% 2nd filtering
RRH3B(1)=PH3;
for i = 1:length(RRH3B)
    if abs(RRH3B(i) - PH3)<=60
        RRH3B(1)= 0.95*RRH3B(1)+0.05*RRH3B(i);
    end
end
for i = 2:length(RRH3B)
    if abs(RRH3B(i) - PH3)<=60
        RRH3B(i)= 0.95*RRH3B(i-1)+0.05*RRH3B(i);
    else
        RRH3B(i) = RRH3B(i-1);
    end
end
for i = length(RRH3B)-1:1
    RRH3B(i)= 0.95*RRH3B(i+1)+0.05*RRH3B(i);
end

FHR3 = 1000.*60.*(1./RRH3B);
%figure;plot(FHR3);hold on;plot(Healthy_3_corr);title('FHR3 2nd filtering');

% 2nd Trimming
U = 15;
L=15;
%figure;plot(FHR3+L);hold on;plot(FHR3-U)
%hold on; plot(Healthy_3); title('FHR3 before 2snd trimming');
for i = 1:length(RRH3)
    if  FHR3(i)- L <= Healthy_3_corr(i) & Healthy_3_corr(i)<= FHR3(i) + L
        FHR3(i)=Healthy_3_corr(i);
        RRH3B(i)=RRH3(i);
    end
end
%figure;plot(FHR3);hold on; plot(Healthy_3_corr);title('FHR3 after 2nd trimming');

% 3rd filtering
RRH3B(1)=PH3;
for i = 1:length(RRH3B)
    if abs(RRH3B(i) - PH3)<=60
        RRH3B(1)= 0.95*RRH3B(1)+0.05*RRH3B(i);
    end
end
for i = 2:length(RRH3B)
    if abs(RRH3B(i) - PH3)<=60
        RRH3B(i)= 0.95*RRH3B(i-1)+0.05*RRH3B(i);
    else
        RRH3B(i) = RRH3B(i-1);
    end
end
for i = length(RRH3B)-1:1
    RRH3B(i)= 0.95*RRH3B(i+3)+0.05*RRH3B(i);
end

FHR3 = 1000.*60.*(1./RRH3B);
%figure;plot(FHR3);hold on;plot(Healthy_3_corr);title('FHR3 3dr filtering');

% 3nd Trimming
U = 10;
L=10;
%figure;plot(FHR1+L);hold on;plot(FHR1-U)
%hold on; plot(Healthy_1); title('FHR1 before 3rd trimming');
for i = 1:length(RRH3)
    if  FHR3(i)- L <= Healthy_3_corr(i) & Healthy_3_corr(i)<= FHR3(i) + L
        FHR3(i)=Healthy_3_corr(i);
        RRH3B(i)=RRH3(i);
    end
end
%figure;plot(FHR3);hold on; plot(Healthy_3_corr);title('FHR3 after 3rd trimming');


% 4rd filtering
RRH3B(1)=PH3;
for i = 1:length(RRH3B)
    if abs(RRH3B(i) - PH3)<=60
        RRH3B(1)= 0.95*RRH3B(1)+0.05*RRH3B(i);
    end
end
for i = 2:length(RRH3B)
    if abs(RRH3B(i) - PH3)<=60
        RRH3B(i)= 0.95*RRH3B(i-1)+0.05*RRH3B(i);
    else
        RRH3B(i) = RRH3B(i-1);
    end
end
for i = length(RRH3B)-1:1
    RRH3B(i)= 0.95*RRH3B(i+1)+0.05*RRH3B(i);
end

FHR3 = 1000.*60.*(1./RRH3B);
%figure;plot(FHR3);hold on;plot(Healthy_3_corr);title('FHR3 4th filtering');

% 4nd Trimming
U = 5;
L=5;
%figure;plot(FHR3+L);hold on;plot(FHR3-U)
%hold on; plot(Healthy_3); title('FHR3 before 4th trimming');
for i = 1:length(RRH3)
    if  FHR3(i)- L <= Healthy_3_corr(i) & Healthy_3_corr(i)<= FHR3(i) + L
        FHR3(i)=Healthy_3_corr(i);
        RRH3B(i)=RRH3(i);
    end
end
%figure;plot(FHR3);hold on; plot(Healthy_3_corr);title('FHR3 after 4th trimming');

% 5th filtering
RRH3B(1)=PH3;
for i = 1:length(RRH3B)
    if abs(RRH3B(i) - PH3)<=60
        RRH3B(1)= 0.95*RRH3B(1)+0.05*RRH3B(i);
    end
end
for i = 2:length(RRH3B)
    if abs(RRH3B(i) - PH3)<=60
        RRH3B(i)= 0.95*RRH3B(i-1)+0.05*RRH3B(i);
    else
        RRH3B(i) = RRH3B(i-1);
    end
end
for i = length(RRH3B)-1:1
    RRH3B(i)= 0.95*RRH3B(i+1)+0.05*RRH3B(i);
end

FHR3 = 1000.*60.*(1./RRH3B);
figure;plot(FHR3);hold on; plot(Healthy_3_corr); title('Baseline of Healthy 3'),legend('Baseline','Healthy 3');

%% From FHR to RR intervals , calculate baseline accordind to Mantel Healthy 4

%Plot RR intervals in time
RRH4 = 1000.* (1./(Healthy_4_corr/60));
%figure; plot(RRH4);

% Calculate the mean av_1 of 2.4s window
RRH4B = RRH4;
time = 4; %2.4sec 
av_4 = zeros(1,floor(length(RRH4)/time));
for i=1:floor(length(RRH4)/time)
        av_4(i)= mean(RRH4((i-1)*time+1:i*time));
end
t_av = 0:time:length(RRH4) - (length(RRH4) - time*floor(length(RRH4)/time))-1;
%figure;plot(t_av,av_4); figure; plot(RRH4B);

% frequency distribution of RR intervals with 1ms BW
[n,edges] = histcounts(av_4,'BinWidth',1);
% figure;histogram(av_4,'BinWidth',1);
% Calculating the peak P as define in Mantel algorithm for the next steop
for i = length(n)-1:-1:2
    if n(i+1) < n(i) & n(i)>=n(i-1) & n(i)>=n(i-2) & n(i)>=n(i-3) & n(i)>=n(i-4) & n(i)>=n(i-5) & sum(n(i:end)) > sum(n)/8
        PH4 = edges(i)
        break
    else 
        continue
    end
end


% 1st Filtering step of RRH5B 
RRH4B(1)=PH4;
i=1;
while  i <= length(RRH4B)
    if (~isnan((RRH4B(i)))) & abs(RRH4B(i) - PH4)<=60
        RRH4B(1)= 0.95*RRH4B(1)+0.05*RRH4B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRH4B)
    if (~isnan((RRH4B(i))))
        if  abs(RRH4B(i) - PH4)<=60
            if isnan((RRH4B(i-1)))
                RRH4B(i)=RRH4B(i);
            else
                RRH4B(i)= 0.95*RRH4B(i-1)+0.05*RRH4B(i);
            end
        else
            if isnan((RRH4B(i-1)))
                RRH4B(i)=RRH4B(i);
            else
                RRH4B(i)= RRH4B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRH4B)-1;
while i >= 1
    if (~isnan((RRH4B(i))))
        if isnan((RRH4B(i+1)))
            RRH4B(i)= RRH4B(i);
        else
            RRH4B(i)= 0.95*RRH4B(i+1)+0.05*RRH4B(i);
        end
    end
    i = i-1;
end

FHR4 = 1000.*60.*(1./RRH4B);
%figure;plot(FHR4);hold on; plot(Healthy_4_corr); title('FHR4 1st filtering');
%figure;plot(RRH4B);

% 1st Trimming at U = L = 20 bpm of FHR4
U = 20;
L=20;
%figure;plot(FHR4+L);hold on;plot(FHR4-U); 
%hold on; plot(Healthy_4_corr); title('FHR1 before 1st trimming');
for i = 1:length(RRH4)
    if  FHR4(i) - L <= Healthy_4_corr(i) & Healthy_4_corr(i)<= FHR4(i) +U
        FHR4(i)=Healthy_4_corr(i);
        RRH4B(i)=RRH4(i);
    end
end
%figure;plot(FHR4);hold on; plot(Healthy_4_corr); title('FHR1 after 1st trimming');

%  2nd filtering of FHR4
RRH4B(1)=PH4;
i=1;
while  i <= length(RRH4B)
    if (~isnan((RRH4B(i)))) & abs(RRH4B(i) - PH4)<=60
        RRH4B(1)= 0.95*RRH4B(1)+0.05*RRH4B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRH4B)
    if (~isnan((RRH4B(i))))
        if  abs(RRH4B(i) - PH4)<=60
            if isnan((RRH4B(i-1)))
                RRH4B(i)=RRH4B(i);
            else
                RRH4B(i)= 0.95*RRH4B(i-1)+0.05*RRH4B(i);
            end
        else
            if isnan((RRH4B(i-1)))
                RRH4B(i)=RRH4B(i);
            else
                RRH4B(i)= RRH4B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRH4B)-1;
while i >= 1
    if (~isnan((RRH4B(i))))
        if isnan((RRH4B(i+1)))
            RRH4B(i)= RRH4B(i);
        else
            RRH4B(i)= 0.95*RRH4B(i+1)+0.05*RRH4B(i);
        end
    end
    i = i-1;
end

FHR4 = 1000.*60.*(1./RRH4B);
%figure;plot(FHR4);hold on; plot(Healthy_4_corr); title('FHR4 2nd filtering');

% 2nd Trimming of FHR4
U = 15;
L=15;
%figure;plot(FHR4+L);hold on;plot(FHR4-U)
%hold on; plot(Healthy_4_corr); title('FHR4 before 2nd trimming');
for i = 1:length(RRH4)
    if  FHR4(i) - L <= Healthy_4_corr(i) & Healthy_4_corr(i)<= FHR4(i) +U
        FHR4(i)=Healthy_4_corr(i);
        RRH4B(i)=RRH4(i);
    end
end
%figure;plot(FHR4);hold on; plot(Healthy_4_corr); title('FHR1 after 2nd trimming');

%  3rd filtering of FHR4
RRH4B(1)=PH4;
i=1;
while  i <= length(RRH4B)
    if (~isnan((RRH4B(i)))) & abs(RRH4B(i) - PH4)<=60
        RRH4B(1)= 0.95*RRH4B(1)+0.05*RRH4B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRH4B)
    if (~isnan((RRH4B(i))))
        if  abs(RRH4B(i) - PH4)<=60
            if isnan((RRH4B(i-1)))
                RRH4B(i)=RRH4B(i);
            else
                RRH4B(i)= 0.95*RRH4B(i-1)+0.05*RRH4B(i);
            end
        else
            if isnan((RRH4B(i-1)))
                RRH4B(i)=RRH4B(i);
            else
                RRH4B(i)= RRH4B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRH4B)-1;
while i >= 1
    if (~isnan((RRH4B(i))))
        if isnan((RRH4B(i+1)))
            RRH4B(i)= RRH4B(i);
        else
            RRH4B(i)= 0.95*RRH4B(i+1)+0.05*RRH4B(i);
        end
    end
    i = i-1;
end

FHR4 = 1000.*60.*(1./RRH4B);
%figure;plot(FHR4);hold on; plot(Healthy_4_corr); title('FHR4 3rd filtering');

% 3nd Trimming of FHR4
U = 10;
L=10;
%figure;plot(FHR4+L);hold on;plot(FHR4-U)
%hold on; plot(Healthy_4_corr); title('FHR4 before 2nd trimming');
for i = 1:length(RRH4)
    if  FHR4(i) - L <= Healthy_4_corr(i) & Healthy_4_corr(i)<= FHR4(i) +U
        FHR4(i)=Healthy_4_corr(i);
        RRH4B(i)=RRH4(i);
    end
end
%figure;plot(FHR4);hold on; plot(Healthy_4_corr); title('FHR4 after 3rd trimming');

%  5rd filtering of FHR4
RRH4B(1)=PH4;
i=1;
while  i <= length(RRH4B)
    if (~isnan((RRH4B(i)))) & abs(RRH4B(i) - PH4)<=60
        RRH4B(1)= 0.95*RRH4B(1)+0.05*RRH4B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRH4B)
    if (~isnan((RRH4B(i))))
        if  abs(RRH4B(i) - PH4)<=60
            if isnan((RRH4B(i-1)))
                RRH4B(i)=RRH4B(i);
            else
                RRH4B(i)= 0.95*RRH4B(i-1)+0.05*RRH4B(i);
            end
        else
            if isnan((RRH4B(i-1)))
                RRH4B(i)=RRH4B(i);
            else
                RRH4B(i)= RRH4B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRH4B)-1;
while i >= 1
    if (~isnan((RRH4B(i))))
        if isnan((RRH4B(i+1)))
            RRH4B(i)= RRH4B(i);
        else
            RRH4B(i)= 0.95*RRH4B(i+1)+0.05*RRH4B(i);
        end
    end
    i = i-1;
end

FHR4 = 1000.*60.*(1./RRH4B);
%figure;plot(FHR4);hold on; plot(Healthy_4_corr); title('FHR4 4th filtering');

% 5nd Trimming of FHR4
U = 5;
L=5;

%figure;plot(FHR4+L);hold on;plot(FHR4-U)
%hold on; plot(Healthy_4_corr); title('FHR4 before 2nd trimming');
for i = 1:length(RRH4)
    if  FHR4(i) - L <= Healthy_4_corr(i) & Healthy_4_corr(i)<= FHR4(i) +U
        FHR4(i)=Healthy_4_corr(i);
        RRH4B(i)=RRH4(i);
    end
end
%figure;plot(FHR4);hold on; plot(Healthy_4_corr); title('FHR4 after 5th trimming');


% 5th filtering of FHR4

RRH4B(1)=PH4;
i=1;
while  i <= length(RRH4B)
    if (~isnan((RRH4B(i)))) & abs(RRH4B(i) - PH4)<=60
        RRH4B(1)= 0.95*RRH4B(1)+0.05*RRH4B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRH4B)
    if (~isnan((RRH4B(i))))
        if  abs(RRH4B(i) - PH4)<=60
            if isnan((RRH4B(i-1)))
                RRH4B(i)=RRH4B(i);
            else
                RRH4B(i)= 0.95*RRH4B(i-1)+0.05*RRH4B(i);
            end
        else
            if isnan((RRH4B(i-1)))
                RRH4B(i)=RRH4B(i);
            else
                RRH4B(i)= RRH4B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRH4B)-1;
while i >= 1
    if (~isnan((RRH4B(i))))
        if isnan((RRH4B(i+1)))
            RRH4B(i)= RRH4B(i);
        else
            RRH4B(i)= 0.95*RRH4B(i+1)+0.05*RRH4B(i);
        end
    end
    i = i-1;
end

FHR4 = 1000.*60.*(1./RRH4B);
figure;plot(FHR4);hold on; plot(Healthy_4_corr); title('Baseline of Healthy 4'),legend('Baseline','Healthy 4');

%% From FHR to RR intervals , calculate baseline accordind to Mantel Healthy 5
%Plot RR intervals in time
RRH5 = 1000.* (1./(Healthy_5_corr/60));
%figure; plot(RRH5);

% Calculate the mean av_1 of 2.5s window
RRH5B = RRH5;
time = 5; %2.5sec 
av_5 = zeros(1,floor(length(RRH5)/time));
for i=1:floor(length(RRH5)/time)
        av_5(i)= mean(RRH5((i-1)*time+1:i*time));
end
t_av = 0:time:length(RRH5) - (length(RRH5) - time*floor(length(RRH5)/time))-1;
%figure;plot(t_av,av_5); figure; plot(RRH5B);

% frequency distribution of RR intervals with 1ms BW
[n,edges] = histcounts(av_5,'BinWidth',1);
% figure;histogram(av_5,'BinWidth',1);
% Calculating the peak P as define in Mantel algorithm for the next steop
for i = length(n)-1:-1:2
    
    if n(i+1) < n(i) 
        continue
    end
    if n(i)>=n(i-1) & n(i)>=n(i-2) & n(i)>=n(i-3) & n(i)>=n(i-5) & n(i)>=n(i-5) & sum(n(i:end)) > sum(n)/8
        PH5 = edges(i)
        break
    else 
        continue
    end
end
% 1st Filtering step of RRH5B 
RRH5B(1)=PH5;
i=1;
while  i <= length(RRH5B)
    if (~isnan((RRH5B(i)))) & abs(RRH5B(i) - PH5)<=60
        RRH5B(1)= 0.95*RRH5B(1)+0.05*RRH5B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRH5B)
    if (~isnan((RRH5B(i))))
        if  abs(RRH5B(i) - PH5)<=60
            if isnan((RRH5B(i-1)))
                RRH5B(i)=RRH5B(i);
            else
                RRH5B(i)= 0.95*RRH5B(i-1)+0.05*RRH5B(i);
            end
        else
            if isnan((RRH5B(i-1)))
                RRH5B(i)=RRH5B(i);
            else
                RRH5B(i)= RRH5B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRH5B)-1;
while i >= 1
    if (~isnan((RRH5B(i))))
        if isnan((RRH5B(i+1)))
            RRH5B(i)= RRH5B(i);
        else
            RRH5B(i)= 0.95*RRH5B(i+1)+0.05*RRH5B(i);
        end
    end
    i = i-1;
end

FHR5 = 1000.*60.*(1./RRH5B);
%figure;plot(FHR5);hold on; plot(Healthy_5_corr); title('FHR5 1st filtering');
%figure;plot(RRH5B);

% 1st Trimming at U = L = 20 bpm,comparison between Health1 and FHR1, change
%in RRH1B and FHR1
U = 20;
L=20;
%figure;plot(FHR5+L);hold on;plot(FHR5-U); 
%hold on; plot(Healthy_5_corr); title('FHR5 before 1st trimming');
for i = 1:length(RRH5)
    if  FHR5(i) - L <= Healthy_5_corr(i) & Healthy_5_corr(i)<= FHR5(i) +U
        FHR5(i)=Healthy_5_corr(i);
        RRH5B(i)=RRH5(i);
    end
end
%figure;plot(FHR5);hold on; plot(Healthy_5_corr); title('FHR1 after 1st trimming');

%  2nd filtering of FRH5
RRH5B(1)=PH5;
i=1;
while  i <= length(RRH5B)
    if (~isnan((RRH5B(i)))) & abs(RRH5B(i) - PH5)<=60
        RRH5B(1)= 0.95*RRH5B(1)+0.05*RRH5B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRH5B)
    if (~isnan((RRH5B(i))))
        if  abs(RRH5B(i) - PH5)<=60
            if isnan((RRH5B(i-1)))
                RRH5B(i)=RRH5B(i);
            else
                RRH5B(i)= 0.95*RRH5B(i-1)+0.05*RRH5B(i);
            end
        else
            if isnan((RRH5B(i-1)))
                RRH5B(i)=RRH5B(i);
            else
                RRH5B(i)= RRH5B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRH5B)-1;
while i >= 1
    if (~isnan((RRH5B(i))))
        if isnan((RRH5B(i+1)))
            RRH5B(i)= RRH5B(i);
        else
            RRH5B(i)= 0.95*RRH5B(i+1)+0.05*RRH5B(i);
        end
    end
    i = i-1;
end

FHR5 = 1000.*60.*(1./RRH5B);
%figure;plot(FHR5);hold on; plot(Healthy_5_corr); title('FHR5 2nd filtering');

% 2nd Trimming of FRH5
U = 15;
L=15;
%figure;plot(FHR5+L);hold on;plot(FHR5-U)
%hold on; plot(Healthy_5_corr); title('FHR5 before 2nd trimming');
for i = 1:length(RRH5)
    if  FHR5(i) - L <= Healthy_5_corr(i) & Healthy_5_corr(i)<= FHR5(i) +U
        FHR5(i)=Healthy_5_corr(i);
        RRH5B(i)=RRH5(i);
    end
end
%figure;plot(FHR5);hold on; plot(Healthy_5_corr); title('FHR1 after 2nd trimming');

% 3rd filtering of FRH5
RRH5B(1)=PH5;
i=1;
while  i <= length(RRH5B)
    if (~isnan((RRH5B(i)))) & abs(RRH5B(i) - PH5)<=60
        RRH5B(1)= 0.95*RRH5B(1)+0.05*RRH5B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRH5B)
    if (~isnan((RRH5B(i))))
        if  abs(RRH5B(i) - PH5)<=60
            if isnan((RRH5B(i-1)))
                RRH5B(i)=RRH5B(i);
            else
                RRH5B(i)= 0.95*RRH5B(i-1)+0.05*RRH5B(i);
            end
        else
            if isnan((RRH5B(i-1)))
                RRH5B(i)=RRH5B(i);
            else
                RRH5B(i)= RRH5B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRH5B)-1;
while i >= 1
    if (~isnan((RRH5B(i))))
        if isnan((RRH5B(i+1)))
            RRH5B(i)= RRH5B(i);
        else
            RRH5B(i)= 0.95*RRH5B(i+1)+0.05*RRH5B(i);
        end
    end
    i = i-1;
end

FHR5 = 1000.*60.*(1./RRH5B);
%figure;plot(FHR5);hold on; plot(Healthy_5_corr); title('FHR5 3rd filtering');

% 3rd Trimming of FHR5
U = 10;
L=10;
%figure;plot(FHR5+L);hold on;plot(FHR5-U)
%hold on; plot(Healthy_5_corr); title('FHR5 before 2nd trimming');
for i = 1:length(RRH5)
    if  FHR5(i) - L <= Healthy_5_corr(i) & Healthy_5_corr(i)<= FHR5(i) +U
        FHR5(i)=Healthy_5_corr(i);
        RRH5B(i)=RRH5(i);
    end
end
%figure;plot(FHR5);hold on; plot(Healthy_5_corr); title('FHR5 after 3rd trimming');

% 4th filtering of FHR5
RRH5B(1)=PH5;
i=1;
while  i <= length(RRH5B)
    if (~isnan((RRH5B(i)))) & abs(RRH5B(i) - PH5)<=60
        RRH5B(1)= 0.95*RRH5B(1)+0.05*RRH5B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRH5B)
    if (~isnan((RRH5B(i))))
        if  abs(RRH5B(i) - PH5)<=60
            if isnan((RRH5B(i-1)))
                RRH5B(i)=RRH5B(i);
            else
                RRH5B(i)= 0.95*RRH5B(i-1)+0.05*RRH5B(i);
            end
        else
            if isnan((RRH5B(i-1)))
                RRH5B(i)=RRH5B(i);
            else
                RRH5B(i)= RRH5B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRH5B)-1;
while i >= 1
    if (~isnan((RRH5B(i))))
        if isnan((RRH5B(i+1)))
            RRH5B(i)= RRH5B(i);
        else
            RRH5B(i)= 0.95*RRH5B(i+1)+0.05*RRH5B(i);
        end
    end
    i = i-1;
end

FHR5 = 1000.*60.*(1./RRH5B);
%figure;plot(FHR5);hold on; plot(Healthy_5_corr); title('FHR5 5th filtering');

% 4TH Trimming
U = 5;
L=5;

%figure;plot(FHR5+L);hold on;plot(FHR5-U)
%hold on; plot(Healthy_5_corr); title('FHR5 before 2nd trimming');
for i = 1:length(RRH5)
    if  FHR5(i) - L <= Healthy_5_corr(i) & Healthy_5_corr(i)<= FHR5(i) +U
        FHR5(i)=Healthy_5_corr(i);
        RRH5B(i)=RRH5(i);
    end
end


%figure;plot(FHR5);hold on; plot(Healthy_5_corr); title('FHR5 after 5th trimming');


% 5th filtering of FHR5

RRH5B(1)=PH5;
i=1;
while  i <= length(RRH5B)
    if (~isnan((RRH5B(i)))) & abs(RRH5B(i) - PH5)<=60
        RRH5B(1)= 0.95*RRH5B(1)+0.05*RRH5B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRH5B)
    if (~isnan((RRH5B(i))))
        if  abs(RRH5B(i) - PH5)<=60
            if isnan((RRH5B(i-1)))
                RRH5B(i)=RRH5B(i);
            else
                RRH5B(i)= 0.95*RRH5B(i-1)+0.05*RRH5B(i);
            end
        else
            if isnan((RRH5B(i-1)))
                RRH5B(i)=RRH5B(i);
            else
                RRH5B(i)= RRH5B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRH5B)-1;
while i >= 1
    if (~isnan((RRH5B(i))))
        if isnan((RRH5B(i+1)))
            RRH5B(i)= RRH5B(i);
        else
            RRH5B(i)= 0.95*RRH5B(i+1)+0.05*RRH5B(i);
        end
    end
    i = i-1;
end

FHR5 = 1000.*60.*(1./RRH5B);
figure;plot(FHR5);hold on; plot(Healthy_5_corr); title('Baseline of Healthy 5'),legend('Baseline','Healthy 5');


%% Calculate baseline accordind to Mantel IUGR 1
%Plot RR intervals in time
RRI1 = 1000.* (1./(IUGR_1_corr/60));
%figure; plot(RRI1); figure; plot(t_IUGR_1,RRI1);

% Calculate the mean av_1 of 2.5s window
RRI1B = RRI1;
time = 5; %2.5sec 
av_1_I = zeros(1,floor(length(RRI1)/time));
for i=1:floor(length(RRI1)/time)
        av_1_I(i)= mean(RRI1((i-1)*time+1:i*time), 'omitnan');
end
t_av_I = 0:time:length(RRI1) - (length(RRI1) - time*floor(length(RRI1)/time))-1;
%figure;plot(t_av_I,av_1_I);hold on; plot(RRI1);

% frequency distribution of RR intervals with 1ms BW
[n,edges] = histcounts(av_1_I,'BinWidth',1);
% figure;histogram(av_1_I,'BinWidth',1);

% Calculating the peak P as define in Mantel algorithm for the next steop
for i = length(n)-1:-1:2
    if n(i+1)< n(i) & n(i)>=n(i-1) & n(i)>=n(i-2) & n(i)>=n(i-3) & n(i)>=n(i-5) & n(i)>=n(i-5) & sum(n(i:end)) > sum(n)/8
        PI1 = edges(i)
        break
    else 
        continue
    end
end

% 1st Filtering step of RRI1B 
RRI1B(1)=PI1;
for i = 1:length(RRI1B)
    if abs(RRI1B(i) - PI1)<=60
        RRI1B(1)= 0.95*RRI1B(1)+0.05*RRI1B(i);
    end
end
for i = 2:length(RRI1B)
    if abs(RRI1B(i) - PI1)<=60
        RRI1B(i)= 0.95*RRI1B(i-1)+0.05*RRI1B(i);
    else
        RRI1B(i) = RRI1B(i-1);
    end
end
for i = length(RRI1B)-1:1
    RRI1B(i)= 0.95*RRI1B(i+1)+0.05*RRI1B(i);
end


%Ploting the resulting RRI1B baseline in RR intervals and in FIR1
%figure; plot(RRI1B);hold on; plot(RRI1);title('RRI1B 1st filtering');
FHR1_IUGR = 1000.*60.*(1./RRI1B);
% 1st Trimming at U = L = 20 bpm OF IUGR1

U = 20;
L=20;
%figure;plot(FHR1_IUGR+L);hold on;plot(FHR1_IUGR-U); 
%hold on; plot(IUGR_1_corr); title('FHR1_I before 1st trimming');
for i = 1:length(RRI1B)
    if  FHR1_IUGR(i) - L <= IUGR_1_corr(i) & IUGR_1_corr(i)<= FHR1_IUGR(i) +U
        FHR1_IUGR(i)=IUGR_1_corr(i);
        RRI1B(i)=RRI1(i);
    end
end
%figure;plot(FHR1_IUGR);hold on; plot(IUGR_1_corr); title('FHR1_IUGR after 1st trimming');
% 2nd filtering OF IUGR 1
RRI1B(1)=PI1;
for i = 1:length(RRI1B)
    if abs(RRI1B(i) - PI1)<=60
        RRI1B(1)= 0.95*RRI1B(1)+0.05*RRI1B(i);
    end
end
for i = 2:length(RRI1B)
    if abs(RRI1B(i) - PI1)<=60
        RRI1B(i)= 0.95*RRI1B(i-1)+0.05*RRI1B(i);
    else
        RRI1B(i) = RRI1B(i-1);
    end
end
for i = length(RRI1B)-1:1
    RRI1B(i)= 0.95*RRI1B(i+1)+0.05*RRI1B(i);
end

FHR1_IUGR = 1000.*60.*(1./RRI1B);
%figure;plot(FHR1_IUGR);hold on;plot(IUGR_1_corr);title('FHR1_IUGR 2nd filtering');
% 2nd Trimming OF IUGR 1
U = 15;
L=15;
%figure;plot(FHR1+L);hold on;plot(FHR1-U)
%hold on; plot(Healthy_1); title('FHR1 before 2snd trimming');
for i = 1:length(RRI1)
    if  FHR1_IUGR(i)- L <= IUGR_1_corr(i) & IUGR_1_corr(i)<= FHR1_IUGR(i) + L
        FHR1_IUGR(i)=IUGR_1_corr(i);
        RRI1B(i)=RRI1(i);
    end
end
%figure;plot(FHR1_IUGR);hold on; plot(IUGR_1_corr);title('FHR1_IUGR after 2nd trimming');

% 3rd filtering OF IUGR 1
RRI1B(1)=PI1;
for i = 1:length(RRI1B)
    if abs(RRI1B(i) - PI1)<=60
        RRI1B(1)= 0.95*RRI1B(1)+0.05*RRI1B(i);
    end
end
for i = 2:length(RRI1B)
    if abs(RRI1B(i) - PI1)<=60
        RRI1B(i)= 0.95*RRI1B(i-1)+0.05*RRI1B(i);
    else
        RRI1B(i) = RRI1B(i-1);
    end
end
for i = length(RRI1B)-1:1
    RRI1B(i)= 0.95*RRI1B(i+1)+0.05*RRI1B(i);
end

FHR1_IUGR = 1000.*60.*(1./RRI1B);
%figure;plot(FHR1_IUGR);hold on;plot(IUGR_1_corr);title('FHR1_IUGR 3dr filtering');
% 3nd TrimmingOF IUGR 1
U = 10;
L=10;
%figure;plot(FHR1+L);hold on;plot(FHR1-U)
%hold on; plot(Healthy_1); title('FHR1 before 3rd trimming');
for i = 1:length(RRI1)
    if  FHR1_IUGR(i)- L <= IUGR_1_corr(i) & IUGR_1_corr(i)<= FHR1_IUGR(i) + L
        FHR1_IUGR(i)=IUGR_1_corr(i);
        RRI1B(i)=RRI1(i);
    end
end
%figure;plot(FHR1_IUGR);hold on; plot(IUGR_1_corr);title('FHR1_IUGR after 3rd trimming');

% 5rd filtering OF IUGR 1
RRI1B(1)=PI1;
for i = 1:length(RRI1B)
    if abs(RRI1B(i) - PI1)<=60
        RRI1B(1)= 0.95*RRI1B(1)+0.05*RRI1B(i);
    end
end
for i = 2:length(RRI1B)
    if abs(RRI1B(i) - PI1)<=60
        RRI1B(i)= 0.95*RRI1B(i-1)+0.05*RRI1B(i);
    else
        RRI1B(i) = RRI1B(i-1);
    end
end
for i = length(RRI1B)-1:1
    RRI1B(i)= 0.95*RRI1B(i+1)+0.05*RRI1B(i);
end

FHR1_IUGR = 1000.*60.*(1./RRI1B);
%figure;plot(FHR1_IUGR);hold on;plot(IUGR_1_corr);title('FHR1_IUGR 5th filtering');
% 5nd Trimming OF IUGR 1
U = 5;
L=5;
%figure;plot(FHR1+L);hold on;plot(FHR1-U)
%hold on; plot(Healthy_1); title('FHR1 before 5th trimming');
for i = 1:length(RRI1)
    if  FHR1_IUGR(i)- L <= IUGR_1_corr(i) & IUGR_1_corr(i)<= FHR1_IUGR(i) + L
        FHR1_IUGR(i)=IUGR_1_corr(i);
        RRI1B(i)=RRI1(i);
    end
end
%figure;plot(FHR1_IUGR);hold on; plot(IUGR_1_corr);title('FHR1_IUGR after 5th trimming');

% 5th filtering OF IUGR 1
RRI1B(1)=PI1;
for i = 1:length(RRI1B)
    if abs(RRI1B(i) - PI1)<=60
        RRI1B(1)= 0.95*RRI1B(1)+0.05*RRI1B(i);
    end
end
for i = 2:length(RRI1B)
    if abs(RRI1B(i) - PI1)<=60
        RRI1B(i)= 0.95*RRI1B(i-1)+0.05*RRI1B(i);
    else
        RRI1B(i) = RRI1B(i-1);
    end
end
for i = length(RRI1B)-1:1
    RRI1B(i)= 0.95*RRI1B(i+1)+0.05*RRI1B(i);
end

FHR1_IUGR = 1000.*60.*(1./RRI1B);
figure;plot(FHR1_IUGR);hold on; plot(IUGR_1_corr); title('Baseline of IUGR 1'),legend('Baseline','IUGR 1');

%% From FHR to RR intervals , calculate baseline accordind to Mantel IUGR 2
%Plot RR intervals in time
RRI2 = 1000.* (1./(IUGR_2_corr/60));
%figure; plot(RRI2);

%Calculate the mean av_1 of 2.5s window
RRI2B = RRI2;
time = 5; %2.5sec 
av_1 = zeros(1,floor(length(RRI2)/time));
for i=1:floor(length(RRI2)/time)
        av_2_I(i)= mean(RRI2((i-1)*time+1:i*time));
end
t_av_I = 0:time:length(RRI2) - (length(RRI2) - time*floor(length(RRI2)/time))-1;
%figure;plot(t_av_I,av_2_I); figure; plot(RRI2B);

% frequency distribution of RR intervals with 1ms BW
[n,edges] = histcounts(av_2_I,'BinWidth',1);
% figure;histogram(av_2_I,'BinWidth',1);

% Calculating the peak P as define in Mantel algorithm for the next steop
for i = length(n)-1:-1:2
    if n(i+1) < n(i) & n(i)>=n(i-1) & n(i)>=n(i-2) & n(i)>=n(i-3) & n(i)>=n(i-5) & n(i)>=n(i-5)  & sum(n(i:end)) >= sum(n)/8
        PI2 = edges(i)
        break
    end
end
% 1st Filtering step of RRI2B 
RRI2B(1)=PI2;
i=1;
while  i <= length(RRI2B)
    if (~isnan((RRI2B(i)))) & abs(RRI2B(i) - PI2)<=60
        RRI2B(1)= 0.95*RRI2B(1)+0.05*RRI2B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRI2B)
    if (~isnan((RRI2B(i))))
        if  abs(RRI2B(i) - PI2)<=60
            if isnan((RRI2B(i-1)))
                RRI2B(i)=RRI2B(i);
            else
                RRI2B(i)= 0.95*RRI2B(i-1)+0.05*RRI2B(i);
            end
        else
            if isnan((RRI2B(i-1)))
                RRI2B(i)=RRI2B(i);
            else
                RRI2B(i)= RRI2B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRI2B)-1;
while i >= 1
    if (~isnan((RRI2B(i))))
        if isnan((RRI2B(i+1)))
            RRI2B(i)= RRI2B(i);
        else
            RRI2B(i)= 0.95*RRI2B(i+1)+0.05*RRI2B(i);
        end
    end
    i = i-1;
end

FHR2_IUGR = 1000.*60.*(1./RRI2B);
%figure;plot(FHR2_IUGR);hold on; plot(IUGR_2_corr); title('FHR2_IUGR 1st filtering');


% 1st Trimming at U = L = 20 bpm pf UIGR 2
U = 20;
L=20;
for i = 1:length(RRI2)
    if  FHR2_IUGR(i) - L <= IUGR_2_corr(i) & IUGR_2_corr(i)<= FHR2_IUGR(i) +U
        FHR2_IUGR(i)=IUGR_2_corr(i);
        RRI2B(i)=RRI2(i);
    end
end
%figure;plot(FHR2_IUGR);hold on; plot(IUGR_2_corr); title('FHR2_IUGR after 1st trimming');

% 2nd filtering OF IUGR 2
RRI2B(1)=PI2;
i=1;
while  i <= length(RRI2B)
    if (~isnan((RRI2B(i)))) & abs(RRI2B(i) - PI2)<=60
        RRI2B(1)= 0.95*RRI2B(1)+0.05*RRI2B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRI2B)
    if (~isnan((RRI2B(i))))
        if  abs(RRI2B(i) - PI2)<=60
            if isnan((RRI2B(i-1)))
                RRI2B(i)=RRI2B(i);
            else
                RRI2B(i)= 0.95*RRI2B(i-1)+0.05*RRI2B(i);
            end
        else
            if isnan((RRI2B(i-1)))
                RRI2B(i)=RRI2B(i);
            else
                RRI2B(i)= RRI2B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRI2B)-1;
while i >= 1
    if (~isnan((RRI2B(i))))
        if isnan((RRI2B(i+1)))
            RRI2B(i)= RRI2B(i);
        else
            RRI2B(i)= 0.95*RRI2B(i+1)+0.05*RRI2B(i);
        end
    end
    i = i-1;
end

FHR2_IUGR = 1000.*60.*(1./RRI2B);
%figure;plot(FHR2_IUGR);hold on; plot(IUGR_2_corr); title('FHR2_IUGR 2nd filtering');

% 2nd Trimming OF IUGR 2
U = 15;
L=15;

for i = 1:length(RRI2)
    if  FHR2_IUGR(i) - L <= IUGR_2_corr(i) & IUGR_2_corr(i)<= FHR2_IUGR(i) +U
        FHR2_IUGR(i)=IUGR_2_corr(i);
        RRI2B(i)=RRI2(i);
    end
end
%figure;plot(FHR2_IUGR);hold on; plot(IUGR_2_corr); title('FHR2_IUGR after 2nd trimming');

% 3rd filtering OF IUGR 2 
RRI2B(1)=PI2;
i=1;
while  i <= length(RRI2B)
    if (~isnan((RRI2B(i)))) & abs(RRI2B(i) - PI2)<=60
        RRI2B(1)= 0.95*RRI2B(1)+0.05*RRI2B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRI2B)
    if (~isnan((RRI2B(i))))
        if  abs(RRI2B(i) - PI2)<=60
            if isnan((RRI2B(i-1)))
                RRI2B(i)=RRI2B(i);
            else
                RRI2B(i)= 0.95*RRI2B(i-1)+0.05*RRI2B(i);
            end
        else
            if isnan((RRI2B(i-1)))
                RRI2B(i)=RRI2B(i);
            else
                RRI2B(i)= RRI2B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRI2B)-1;
while i >= 1
    if (~isnan((RRI2B(i))))
        if isnan((RRI2B(i+1)))
            RRI2B(i)= RRI2B(i);
        else
            RRI2B(i)= 0.95*RRI2B(i+1)+0.05*RRI2B(i);
        end
    end
    i = i-1;
end

FHR2_IUGR = 1000.*60.*(1./RRI2B);
%figure;plot(FHR2_IUGR);hold on; plot(IUGR_2_corr); title('FHR2_IUGR 3rd filtering');

% 3nd Trimming OF IUGR 2
U = 10;
L=10;

for i = 1:length(RRI2)
    if  FHR2_IUGR(i) - L <= IUGR_2_corr(i) & IUGR_2_corr(i)<= FHR2_IUGR(i) +U
        FHR2_IUGR(i)=IUGR_2_corr(i);
        RRI2B(i)=RRI2(i);
    end
end
%figure;plot(FHR2_IUGR);hold on; plot(IUGR_2_corr); title('FHR2_IUGR after 3rd trimming');

% 4TH filtering OF IUGR 2 
RRI2B(1)=PI2;
i=1;
while  i <= length(RRI2B)
    if (~isnan((RRI2B(i)))) & abs(RRI2B(i) - PI2)<=60
        RRI2B(1)= 0.95*RRI2B(1)+0.05*RRI2B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRI2B)
    if (~isnan((RRI2B(i))))
        if  abs(RRI2B(i) - PI2)<=60
            if isnan((RRI2B(i-1)))
                RRI2B(i)=RRI2B(i);
            else
                RRI2B(i)= 0.95*RRI2B(i-1)+0.05*RRI2B(i);
            end
        else
            if isnan((RRI2B(i-1)))
                RRI2B(i)=RRI2B(i);
            else
                RRI2B(i)= RRI2B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRI2B)-1;
while i >= 1
    if (~isnan((RRI2B(i))))
        if isnan((RRI2B(i+1)))
            RRI2B(i)= RRI2B(i);
        else
            RRI2B(i)= 0.95*RRI2B(i+1)+0.05*RRI2B(i);
        end
    end
    i = i-1;
end

FHR2_IUGR = 1000.*60.*(1./RRI2B);
%figure;plot(FHR2_IUGR);hold on; plot(IUGR_2_corr); title('FHR2_IUGR 5th filtering');

% 5TH Trimming OF IUGR 2
U = 5;
L=5;

for i = 1:length(RRI2)
    if  FHR2_IUGR(i) - L <= IUGR_2_corr(i) & IUGR_2_corr(i)<= FHR2_IUGR(i) +U
        FHR2_IUGR(i)=IUGR_2_corr(i);
        RRI2B(i)=RRI2(i);
    end
end


% 5th filtering OF IUGR 2

RRI2B(1)=PI2;
i=1;
while  i <= length(RRI2B)
    if (~isnan((RRI2B(i)))) & abs(RRI2B(i) - PI2)<=60
        RRI2B(1)= 0.95*RRI2B(1)+0.05*RRI2B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRI2B)
    if (~isnan((RRI2B(i))))
        if  abs(RRI2B(i) - PI2)<=60
            if isnan((RRI2B(i-1)))
                RRI2B(i)=RRI2B(i);
            else
                RRI2B(i)= 0.95*RRI2B(i-1)+0.05*RRI2B(i);
            end
        else
            if isnan((RRI2B(i-1)))
                RRI2B(i)=RRI2B(i);
            else
                RRI2B(i)= RRI2B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRI2B)-1;
while i >= 1
    if (~isnan((RRI2B(i))))
        if isnan((RRI2B(i+1)))
            RRI2B(i)= RRI2B(i);
        else
            RRI2B(i)= 0.95*RRI2B(i+1)+0.05*RRI2B(i);
        end
    end
    i = i-1;
end

FHR2_IUGR = 1000.*60.*(1./RRI2B);

figure;plot(FHR2_IUGR);hold on; plot(IUGR_2_corr); title('Baseline of IUGR 2'),legend('Baseline','IUGR 2');
%% From FHR to RR intervals , calculate baseline accordind to Mantel IUGR 3
%Plot RR intervals in time
RRI3 = 1000.* (1./(IUGR_3_corr/60));
%figure; plot(RRI3);

%Calculate the mean av_1 of 2.5s window
RRI3B = RRI3;
time = 5; %2.5sec 
av_3_I = zeros(1,floor(length(RRI3)/time));
for i=1:floor(length(RRI3)/time)
        av_3_I(i)= mean(RRI3((i-1)*time+1:i*time));
end
t_av_I = 0:time:length(RRI3) - (length(RRI3) - time*floor(length(RRI3)/time))-1;
%figure;plot(t_av_I,av_3_I); figure; plot(RRI3B);

% frequency distribution of RR intervals with 1ms BW
[n,edges] = histcounts(av_3_I,'BinWidth',1);
%figure;histogram(av_3_I,'BinWidth',1);

% Calculating the peak P as define in Mantel algorithm for the next steop
for i = length(n)-1:-1:2
    if n(i+1) < n(i) & n(i)>=n(i-1) & n(i)>=n(i-2) & n(i)>=n(i-3) & n(i)>=n(i-5) & n(i)>=n(i-5)  & sum(n(i:end)) >= sum(n)/8
        PI3 = edges(i)
        break
    end
end
% 1st Filtering step of RRI3B 
RRI3B(1)=PI3;
i=1;
while  i <= length(RRI3B)
    if (~isnan((RRI3B(i)))) & abs(RRI3B(i) - PI3)<=60
        RRI3B(1)= 0.95*RRI3B(1)+0.05*RRI3B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRI3B)
    if (~isnan((RRI3B(i))))
        if  abs(RRI3B(i) - PI3)<=60
            if isnan((RRI3B(i-1)))
                RRI3B(i)=RRI3B(i);
            else
                RRI3B(i)= 0.95*RRI3B(i-1)+0.05*RRI3B(i);
            end
        else
            if isnan((RRI3B(i-1)))
                RRI3B(i)=RRI3B(i);
            else
                RRI3B(i)= RRI3B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRI3B)-1;
while i >= 1
    if (~isnan((RRI3B(i))))
        if isnan((RRI3B(i+1)))
            RRI3B(i)= RRI3B(i);
        else
            RRI3B(i)= 0.95*RRI3B(i+1)+0.05*RRI3B(i);
        end
    end
    i = i-1;
end

FHR3_IUGR = 1000.*60.*(1./RRI3B);
%figure;plot(FHR3_IUGR);hold on; plot(IUGR_3_corr); title('FHR3_IUGR 1st filtering');


% 1st Trimming at U = L = 20 bpm pf UIGR 3
U = 20;
L=20;
for i = 1:length(RRI3)
    if  FHR3_IUGR(i) - L <= IUGR_3_corr(i) & IUGR_3_corr(i)<= FHR3_IUGR(i) +U
        FHR3_IUGR(i)=IUGR_3_corr(i);
        RRI3B(i)=RRI3(i);
    end
end

%figure;plot(FHR3_IUGR);hold on; plot(IUGR_3_corr); title('FHR3_IUGR after 1st trimming');
% 2nd filtering OF IUGR "
RRI3B(1)=PI3;
i=1;
while  i <= length(RRI3B)
    if (~isnan((RRI3B(i)))) & abs(RRI3B(i) - PI3)<=60
        RRI3B(1)= 0.95*RRI3B(1)+0.05*RRI3B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRI3B)
    if (~isnan((RRI3B(i))))
        if  abs(RRI3B(i) - PI3)<=60
            if isnan((RRI3B(i-1)))
                RRI3B(i)=RRI3B(i);
            else
                RRI3B(i)= 0.95*RRI3B(i-1)+0.05*RRI3B(i);
            end
        else
            if isnan((RRI3B(i-1)))
                RRI3B(i)=RRI3B(i);
            else
                RRI3B(i)= RRI3B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRI3B)-1;
while i >= 1
    if (~isnan((RRI3B(i))))
        if isnan((RRI3B(i+1)))
            RRI3B(i)= RRI3B(i);
        else
            RRI3B(i)= 0.95*RRI3B(i+1)+0.05*RRI3B(i);
        end
    end
    i = i-1;
end

FHR3_IUGR = 1000.*60.*(1./RRI3B);
%figure;plot(FHR3_IUGR);hold on; plot(IUGR_3_corr); title('FHR3_IUGR 2nd filtering');

% 2nd Trimming OF IUGR 2
U = 15;
L=15;

for i = 1:length(RRI3)
    if  FHR3_IUGR(i) - L <= IUGR_3_corr(i) & IUGR_3_corr(i)<= FHR3_IUGR(i) +U
        FHR3_IUGR(i)=IUGR_3_corr(i);
        RRI3B(i)=RRI3(i);
    end
end
%figure;plot(FHR3_IUGR);hold on; plot(IUGR_3_corr); title('FHR3_IUGR after 2nd trimming');

% 3rd filtering OF IUGR 2 
RRI3B(1)=PI3;
i=1;
while  i <= length(RRI3B)
    if (~isnan((RRI3B(i)))) & abs(RRI3B(i) - PI3)<=60
        RRI3B(1)= 0.95*RRI3B(1)+0.05*RRI3B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRI3B)
    if (~isnan((RRI3B(i))))
        if  abs(RRI3B(i) - PI3)<=60
            if isnan((RRI3B(i-1)))
                RRI3B(i)=RRI3B(i);
            else
                RRI3B(i)= 0.95*RRI3B(i-1)+0.05*RRI3B(i);
            end
        else
            if isnan((RRI3B(i-1)))
                RRI3B(i)=RRI3B(i);
            else
                RRI3B(i)= RRI3B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRI3B)-1;
while i >= 1
    if (~isnan((RRI3B(i))))
        if isnan((RRI3B(i+1)))
            RRI3B(i)= RRI3B(i);
        else
            RRI3B(i)= 0.95*RRI3B(i+1)+0.05*RRI3B(i);
        end
    end
    i = i-1;
end

FHR3_IUGR = 1000.*60.*(1./RRI3B);
%figure;plot(FHR3_IUGR);hold on; plot(IUGR_3_corr); title('FHR3_IUGR 3rd filtering');

% 3nd Trimming OF IUGR 2
U = 10;
L=10;

for i = 1:length(RRI3)
    if  FHR3_IUGR(i) - L <= IUGR_3_corr(i) & IUGR_3_corr(i)<= FHR3_IUGR(i) +U
        FHR3_IUGR(i)=IUGR_3_corr(i);
        RRI3B(i)=RRI3(i);
    end
end
%figure;plot(FHR3_IUGR);hold on; plot(IUGR_3_corr); title('FHR3_IUGR after 3rd trimming');

% 5rd filtering OF IUGR 2 
RRI3B(1)=PI3;
i=1;
while  i <= length(RRI3B)
    if (~isnan((RRI3B(i)))) & abs(RRI3B(i) - PI3)<=60
        RRI3B(1)= 0.95*RRI3B(1)+0.05*RRI3B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRI3B)
    if (~isnan((RRI3B(i))))
        if  abs(RRI3B(i) - PI3)<=60
            if isnan((RRI3B(i-1)))
                RRI3B(i)=RRI3B(i);
            else
                RRI3B(i)= 0.95*RRI3B(i-1)+0.05*RRI3B(i);
            end
        else
            if isnan((RRI3B(i-1)))
                RRI3B(i)=RRI3B(i);
            else
                RRI3B(i)= RRI3B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRI3B)-1;
while i >= 1
    if (~isnan((RRI3B(i))))
        if isnan((RRI3B(i+1)))
            RRI3B(i)= RRI3B(i);
        else
            RRI3B(i)= 0.95*RRI3B(i+1)+0.05*RRI3B(i);
        end
    end
    i = i-1;
end

FHR3_IUGR = 1000.*60.*(1./RRI3B);
%figure;plot(FHR3_IUGR);hold on; plot(IUGR_3_corr); title('FHR3_IUGR 5th filtering');

% 5TH Trimming OF IUGR 2
U = 5;
L=5;

for i = 1:length(RRI3)
    if  FHR3_IUGR(i) - L <= IUGR_3_corr(i) & IUGR_3_corr(i)<= FHR3_IUGR(i) +U
        FHR3_IUGR(i)=IUGR_3_corr(i);
        RRI3B(i)=RRI3(i);
    end
end


% 5th filtering OF IUGR 2

RRI3B(1)=PI3;
i=1;
while  i <= length(RRI3B)
    if (~isnan((RRI3B(i)))) & abs(RRI3B(i) - PI3)<=60
        RRI3B(1)= 0.95*RRI3B(1)+0.05*RRI3B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRI3B)
    if (~isnan((RRI3B(i))))
        if  abs(RRI3B(i) - PI3)<=60
            if isnan((RRI3B(i-1)))
                RRI3B(i)=RRI3B(i);
            else
                RRI3B(i)= 0.95*RRI3B(i-1)+0.05*RRI3B(i);
            end
        else
            if isnan((RRI3B(i-1)))
                RRI3B(i)=RRI3B(i);
            else
                RRI3B(i)= RRI3B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRI3B)-1;
while i >= 1
    if (~isnan((RRI3B(i))))
        if isnan((RRI3B(i+1)))
            RRI3B(i)= RRI3B(i);
        else
            RRI3B(i)= 0.95*RRI3B(i+1)+0.05*RRI3B(i);
        end
    end
    i = i-1;
end

FHR3_IUGR = 1000.*60.*(1./RRI3B);

figure;plot(FHR3_IUGR);hold on; plot(IUGR_3_corr); title('Baseline of IUGR 3'),legend('Baseline','IUGR 3');

%% From FHR to RR intervals , calculate baseline accordind to Mantel IUGR 4
%Plot RR intervals in time
RRI4 = 1000.* (1./(IUGR_4_corr/60));
%figure; plot(RRI4);

%Calculate the mean av_1 of 2.5s window
RRI4B = RRI4;
time = 5; %2.5sec 
av_4_I = zeros(1,floor(length(RRI4)/time));
for i=1:floor(length(RRI4)/time)
        av_4_I(i)= mean(RRI4((i-1)*time+1:i*time));
end
t_av_I = 0:time:length(RRI4) - (length(RRI4) - time*floor(length(RRI4)/time))-1;
%figure;plot(t_av_I,av_4_I); figure; plot(RRI4B);

% frequency distribution of RR intervals with 1ms BW
[n,edges] = histcounts(av_4_I,'BinWidth',1);
% figure;histogram(av_4_I,'BinWidth',1);

% Calculating the peak P as define in Mantel algorithm for the next steop
for i = length(n)-1:-1:2
    if n(i+1) < n(i) & n(i)>=n(i-1) & n(i)>=n(i-2) & n(i)>=n(i-3) & n(i)>=n(i-4) & n(i)>=n(i-5)  & sum(n(i:end)) >= sum(n)/8
        PI4 = edges(i)
        break
    end
end
% 1st Filtering step of RRI4B 
RRI4B(1)=PI4;
i=1;
while  i <= length(RRI4B)
    if (~isnan((RRI4B(i)))) & abs(RRI4B(i) - PI4)<=60
        RRI4B(1)= 0.95*RRI4B(1)+0.05*RRI4B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRI4B)
    if (~isnan((RRI4B(i))))
        if  abs(RRI4B(i) - PI4)<=60
            if isnan((RRI4B(i-1)))
                RRI4B(i)=RRI4B(i);
            else
                RRI4B(i)= 0.95*RRI4B(i-1)+0.05*RRI4B(i);
            end
        else
            if isnan((RRI4B(i-1)))
                RRI4B(i)=RRI4B(i);
            else
                RRI4B(i)= RRI4B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRI4B)-1;
while i >= 1
    if (~isnan((RRI4B(i))))
        if isnan((RRI4B(i+1)))
            RRI4B(i)= RRI4B(i);
        else
            RRI4B(i)= 0.95*RRI4B(i+1)+0.05*RRI4B(i);
        end
    end
    i = i-1;
end

FHR4_IUGR = 1000.*60.*(1./RRI4B);
%figure;plot(FHR4_IUGR);hold on; plot(IUGR_4_corr); title('FHR4_IUGR 1st filtering');


% 1st Trimming at U = L = 20 bpm pf UIGR 5
U = 20;
L=20;
for i = 1:length(RRI4)
    if  FHR4_IUGR(i) - L <= IUGR_4_corr(i) & IUGR_4_corr(i)<= FHR4_IUGR(i) +U
        FHR4_IUGR(i)=IUGR_4_corr(i);
        RRI4B(i)=RRI4(i);
    end
end

%figure;plot(FHR4_IUGR);hold on; plot(IUGR_4_corr); title('FHR4_IUGR after 1st trimming');
% 2nd filtering OF IUGR "
RRI4B(1)=PI4;
i=1;
while  i <= length(RRI4B)
    if (~isnan((RRI4B(i)))) & abs(RRI4B(i) - PI4)<=60
        RRI4B(1)= 0.95*RRI4B(1)+0.05*RRI4B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRI4B)
    if (~isnan((RRI4B(i))))
        if  abs(RRI4B(i) - PI4)<=60
            if isnan((RRI4B(i-1)))
                RRI4B(i)=RRI4B(i);
            else
                RRI4B(i)= 0.95*RRI4B(i-1)+0.05*RRI4B(i);
            end
        else
            if isnan((RRI4B(i-1)))
                RRI4B(i)=RRI4B(i);
            else
                RRI4B(i)= RRI4B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRI4B)-1;
while i >= 1
    if (~isnan((RRI4B(i))))
        if isnan((RRI4B(i+1)))
            RRI4B(i)= RRI4B(i);
        else
            RRI4B(i)= 0.95*RRI4B(i+1)+0.05*RRI4B(i);
        end
    end
    i = i-1;
end

FHR4_IUGR = 1000.*60.*(1./RRI4B);
%figure;plot(FHR4_IUGR);hold on; plot(IUGR_4_corr); title('FHR4_IUGR 2nd filtering');

% 2nd Trimming OF IUGR 2
U = 15;
L=15;

for i = 1:length(RRI4)
    if  FHR4_IUGR(i) - L <= IUGR_4_corr(i) & IUGR_4_corr(i)<= FHR4_IUGR(i) +U
        FHR4_IUGR(i)=IUGR_4_corr(i);
        RRI4B(i)=RRI4(i);
    end
end
%figure;plot(FHR4_IUGR);hold on; plot(IUGR_4_corr); title('FHR4_IUGR after 2nd trimming');

% 3rd filtering OF IUGR 2 
RRI4B(1)=PI4;
i=1;
while  i <= length(RRI4B)
    if (~isnan((RRI4B(i)))) & abs(RRI4B(i) - PI4)<=60
        RRI4B(1)= 0.95*RRI4B(1)+0.05*RRI4B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRI4B)
    if (~isnan((RRI4B(i))))
        if  abs(RRI4B(i) - PI4)<=60
            if isnan((RRI4B(i-1)))
                RRI4B(i)=RRI4B(i);
            else
                RRI4B(i)= 0.95*RRI4B(i-1)+0.05*RRI4B(i);
            end
        else
            if isnan((RRI4B(i-1)))
                RRI4B(i)=RRI4B(i);
            else
                RRI4B(i)= RRI4B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRI4B)-1;
while i >= 1
    if (~isnan((RRI4B(i))))
        if isnan((RRI4B(i+1)))
            RRI4B(i)= RRI4B(i);
        else
            RRI4B(i)= 0.95*RRI4B(i+1)+0.05*RRI4B(i);
        end
    end
    i = i-1;
end

FHR4_IUGR = 1000.*60.*(1./RRI4B);
%figure;plot(FHR4_IUGR);hold on; plot(IUGR_4_corr); title('FHR4_IUGR 3rd filtering');

% 3nd Trimming OF IUGR 2
U = 10;
L=10;

for i = 1:length(RRI4)
    if  FHR4_IUGR(i) - L <= IUGR_4_corr(i) & IUGR_4_corr(i)<= FHR4_IUGR(i) +U
        FHR4_IUGR(i)=IUGR_4_corr(i);
        RRI4B(i)=RRI4(i);
    end
end
%figure;plot(FHR4_IUGR);hold on; plot(IUGR_4_corr); title('FHR4_IUGR after 3rd trimming');

% 4rd filtering OF IUGR 2 
RRI4B(1)=PI4;
i=1;
while  i <= length(RRI4B)
    if (~isnan((RRI4B(i)))) & abs(RRI4B(i) - PI4)<=60
        RRI4B(1)= 0.95*RRI4B(1)+0.05*RRI4B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRI4B)
    if (~isnan((RRI4B(i))))
        if  abs(RRI4B(i) - PI4)<=60
            if isnan((RRI4B(i-1)))
                RRI4B(i)=RRI4B(i);
            else
                RRI4B(i)= 0.95*RRI4B(i-1)+0.05*RRI4B(i);
            end
        else
            if isnan((RRI4B(i-1)))
                RRI4B(i)=RRI4B(i);
            else
                RRI4B(i)= RRI4B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRI4B)-1;
while i >= 1
    if (~isnan((RRI4B(i))))
        if isnan((RRI4B(i+1)))
            RRI4B(i)= RRI4B(i);
        else
            RRI4B(i)= 0.95*RRI4B(i+1)+0.05*RRI4B(i);
        end
    end
    i = i-1;
end

FHR4_IUGR = 1000.*60.*(1./RRI4B);
%figure;plot(FHR4_IUGR);hold on; plot(IUGR_4_corr); title('FHR4_IUGR 4th filtering');

% 5TH Trimming OF IUGR 2
U = 5;
L=5;

for i = 1:length(RRI4)
    if  FHR4_IUGR(i) - L <= IUGR_4_corr(i) & IUGR_4_corr(i)<= FHR4_IUGR(i) +U
        FHR4_IUGR(i)=IUGR_4_corr(i);
        RRI4B(i)=RRI4(i);
    end
end


% 5th filtering OF IUGR 4

RRI4B(1)=PI4;
i=1;
while  i <= length(RRI4B)
    if (~isnan((RRI4B(i)))) & abs(RRI4B(i) - PI4)<=60
        RRI4B(1)= 0.95*RRI4B(1)+0.05*RRI4B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRI4B)
    if (~isnan((RRI4B(i))))
        if  abs(RRI4B(i) - PI4)<=60
            if isnan((RRI4B(i-1)))
                RRI4B(i)=RRI4B(i);
            else
                RRI4B(i)= 0.95*RRI4B(i-1)+0.05*RRI4B(i);
            end
        else
            if isnan((RRI4B(i-1)))
                RRI4B(i)=RRI4B(i);
            else
                RRI4B(i)= RRI4B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRI4B)-1;
while i >= 1
    if (~isnan((RRI4B(i))))
        if isnan((RRI4B(i+1)))
            RRI4B(i)= RRI4B(i);
        else
            RRI4B(i)= 0.95*RRI4B(i+1)+0.05*RRI4B(i);
        end
    end
    i = i-1;
end

FHR4_IUGR = 1000.*60.*(1./RRI4B);

figure;plot(FHR4_IUGR);hold on; plot(IUGR_4_corr); title('Baseline of IUGR 4'),legend('Baseline','IUGR 4');

%% From FHR to RR intervals , calculate baseline accordind to Mantel IUGR 5
%Plot RR intervals in time
RRI5 = 1000.* (1./(IUGR_5_corr/60));
%figure; plot(RRI5);

%Calculate the mean av_1 of 2.5s window
RRI5B = RRI5;
time = 5; %2.5sec 
av_5_I = zeros(1,floor(length(RRI5)/time));
for i=1:floor(length(RRI5)/time)
        av_5_I(i)= mean(RRI5((i-1)*time+1:i*time));
end
t_av_I = 0:time:length(RRI5) - (length(RRI5) - time*floor(length(RRI5)/time))-1;
%figure;plot(t_av_I,av_5_I); figure; plot(RRI5B);

% frequency distribution of RR intervals with 1ms BW
[n,edges] = histcounts(av_5_I,'BinWidth',1);
figure;histogram(av_5_I,'BinWidth',1), xlabel('bpm'), ylabel('counts');

% Calculating the peak P as define in Mantel algorithm for the next steop
for i = length(n)-1:-1:2
    if n(i+1) < n(i) & n(i)>=n(i-1) & n(i)>=n(i-2) & n(i)>=n(i-3) & n(i)>=n(i-4) & n(i)>=n(i-5)  & sum(n(i:end)) >= sum(n)/8
        PI5 = edges(i)
        break
    end
end
% 1st Filtering step of RRI5B 
RRI5B(1)=PI5;
i=1;
while  i <= length(RRI5B)
    if (~isnan((RRI5B(i)))) & abs(RRI5B(i) - PI5)<=60
        RRI5B(1)= 0.95*RRI5B(1)+0.05*RRI5B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRI5B)
    if (~isnan((RRI5B(i))))
        if  abs(RRI5B(i) - PI5)<=60
            if isnan((RRI5B(i-1)))
                RRI5B(i)=RRI5B(i);
            else
                RRI5B(i)= 0.95*RRI5B(i-1)+0.05*RRI5B(i);
            end
        else
            if isnan((RRI5B(i-1)))
                RRI5B(i)=RRI5B(i);
            else
                RRI5B(i)= RRI5B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRI5B)-1;
while i >= 1
    if (~isnan((RRI5B(i))))
        if isnan((RRI5B(i+1)))
            RRI5B(i)= RRI5B(i);
        else
            RRI5B(i)= 0.95*RRI5B(i+1)+0.05*RRI5B(i);
        end
    end
    i = i-1;
end

FHR5_IUGR = 1000.*60.*(1./RRI5B);
%figure;plot(FHR5_IUGR);hold on; plot(IUGR_5_corr); title('FHR5_IUGR 1st filtering');


% 1st Trimming at U = L = 20 bpm pf UIGR 5
U = 20;
L=20;
for i = 1:length(RRI5)
    if  FHR5_IUGR(i) - L <= IUGR_5_corr(i) & IUGR_5_corr(i)<= FHR5_IUGR(i) +U
        FHR5_IUGR(i)=IUGR_5_corr(i);
        RRI5B(i)=RRI5(i);
    end
end

%figure;plot(FHR5_IUGR);hold on; plot(IUGR_5_corr); title('FHR5_IUGR after 1st trimming');
% 2nd filtering OF IUGR "
RRI5B(1)=PI5;
i=1;
while  i <= length(RRI5B)
    if (~isnan((RRI5B(i)))) & abs(RRI5B(i) - PI5)<=60
        RRI5B(1)= 0.95*RRI5B(1)+0.05*RRI5B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRI5B)
    if (~isnan((RRI5B(i))))
        if  abs(RRI5B(i) - PI5)<=60
            if isnan((RRI5B(i-1)))
                RRI5B(i)=RRI5B(i);
            else
                RRI5B(i)= 0.95*RRI5B(i-1)+0.05*RRI5B(i);
            end
        else
            if isnan((RRI5B(i-1)))
                RRI5B(i)=RRI5B(i);
            else
                RRI5B(i)= RRI5B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRI5B)-1;
while i >= 1
    if (~isnan((RRI5B(i))))
        if isnan((RRI5B(i+1)))
            RRI5B(i)= RRI5B(i);
        else
            RRI5B(i)= 0.95*RRI5B(i+1)+0.05*RRI5B(i);
        end
    end
    i = i-1;
end

FHR5_IUGR = 1000.*60.*(1./RRI5B);
%figure;plot(FHR5_IUGR);hold on; plot(IUGR_5_corr); title('FHR5_IUGR 2nd filtering');

% 2nd Trimming OF IUGR 2
U = 15;
L=15;

for i = 1:length(RRI5)
    if  FHR5_IUGR(i) - L <= IUGR_5_corr(i) & IUGR_5_corr(i)<= FHR5_IUGR(i) +U
        FHR5_IUGR(i)=IUGR_5_corr(i);
        RRI5B(i)=RRI5(i);
    end
end
%figure;plot(FHR5_IUGR);hold on; plot(IUGR_5_corr); title('FHR5_IUGR after 2nd trimming');

% 3rd filtering OF IUGR 5
RRI5B(1)=PI5;
i=1;
while  i <= length(RRI5B)
    if (~isnan((RRI5B(i)))) & abs(RRI5B(i) - PI5)<=60
        RRI5B(1)= 0.95*RRI5B(1)+0.05*RRI5B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRI5B)
    if (~isnan((RRI5B(i))))
        if  abs(RRI5B(i) - PI5)<=60
            if isnan((RRI5B(i-1)))
                RRI5B(i)=RRI5B(i);
            else
                RRI5B(i)= 0.95*RRI5B(i-1)+0.05*RRI5B(i);
            end
        else
            if isnan((RRI5B(i-1)))
                RRI5B(i)=RRI5B(i);
            else
                RRI5B(i)= RRI5B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRI5B)-1;
while i >= 1
    if (~isnan((RRI5B(i))))
        if isnan((RRI5B(i+1)))
            RRI5B(i)= RRI5B(i);
        else
            RRI5B(i)= 0.95*RRI5B(i+1)+0.05*RRI5B(i);
        end
    end
    i = i-1;
end

FHR5_IUGR = 1000.*60.*(1./RRI5B);
%figure;plot(FHR5_IUGR);hold on; plot(IUGR_5_corr); title('FHR5_IUGR 3rd filtering');

% 3nd Trimming OF IUGR 2
U = 10;
L=10;

for i = 1:length(RRI5)
    if  FHR5_IUGR(i) - L <= IUGR_5_corr(i) & IUGR_5_corr(i)<= FHR5_IUGR(i) +U
        FHR5_IUGR(i)=IUGR_5_corr(i);
        RRI5B(i)=RRI5(i);
    end
end
%figure;plot(FHR5_IUGR);hold on; plot(IUGR_5_corr); title('FHR5_IUGR after 3rd trimming');

% 5rd filtering OF IUGR 2 
RRI5B(1)=PI5;
i=1;
while  i <= length(RRI5B)
    if (~isnan((RRI5B(i)))) & abs(RRI5B(i) - PI5)<=60
        RRI5B(1)= 0.95*RRI5B(1)+0.05*RRI5B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRI5B)
    if (~isnan((RRI5B(i))))
        if  abs(RRI5B(i) - PI5)<=60
            if isnan((RRI5B(i-1)))
                RRI5B(i)=RRI5B(i);
            else
                RRI5B(i)= 0.95*RRI5B(i-1)+0.05*RRI5B(i);
            end
        else
            if isnan((RRI5B(i-1)))
                RRI5B(i)=RRI5B(i);
            else
                RRI5B(i)= RRI5B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRI5B)-1;
while i >= 1
    if (~isnan((RRI5B(i))))
        if isnan((RRI5B(i+1)))
            RRI5B(i)= RRI5B(i);
        else
            RRI5B(i)= 0.95*RRI5B(i+1)+0.05*RRI5B(i);
        end
    end
    i = i-1;
end

FHR5_IUGR = 1000.*60.*(1./RRI5B);
%figure;plot(FHR5_IUGR);hold on; plot(IUGR_5_corr); title('FHR5_IUGR 5th filtering');

% 5TH Trimming OF IUGR 5
U = 5;
L=5;

for i = 1:length(RRI5)
    if  FHR5_IUGR(i) - L <= IUGR_5_corr(i) & IUGR_5_corr(i)<= FHR5_IUGR(i) +U
        FHR5_IUGR(i)=IUGR_5_corr(i);
        RRI5B(i)=RRI5(i);
    end
end


% 5th filtering OF IUGR 5

RRI5B(1)=PI5;
i=1;
while  i <= length(RRI5B)
    if (~isnan((RRI5B(i)))) & abs(RRI5B(i) - PI5)<=60
        RRI5B(1)= 0.95*RRI5B(1)+0.05*RRI5B(i);
    end
    i=i+1;
end
i = 2;
while i <=length(RRI5B)
    if (~isnan((RRI5B(i))))
        if  abs(RRI5B(i) - PI5)<=60
            if isnan((RRI5B(i-1)))
                RRI5B(i)=RRI5B(i);
            else
                RRI5B(i)= 0.95*RRI5B(i-1)+0.05*RRI5B(i);
            end
        else
            if isnan((RRI5B(i-1)))
                RRI5B(i)=RRI5B(i);
            else
                RRI5B(i)= RRI5B(i-1);
            end
        end
    end
    i=i+1;
end

i = length(RRI5B)-1;
while i >= 1
    if (~isnan((RRI5B(i))))
        if isnan((RRI5B(i+1)))
            RRI5B(i)= RRI5B(i);
        else
            RRI5B(i)= 0.95*RRI5B(i+1)+0.05*RRI5B(i);
        end
    end
    i = i-1;
end

FHR5_IUGR = 1000.*60.*(1./RRI5B);

figure;plot(FHR5_IUGR);hold on; plot(IUGR_5_corr); title('Baseline of IUGR 5'),legend('Baseline','IUGR 5');


%% Number of accelerations

%ACCELERATIONS HEALTHY_1
% figure, plot(Healthy_1_corr), hold on ,plot(FHR1), hold on, plot(FHR1+15);
nb_acc_H1 = 0;
count = 0;
i=1;
while i <= length(Healthy_1_corr)
    if Healthy_1_corr(i) >= FHR1(i) + 15
        count = count +1;
        j = 1;
        while Healthy_1_corr(i+j) >= FHR1(i+j) + 15
            count = count + 1;
            j = j+1;
        end
        if count >= 30
            nb_acc_H1=nb_acc_H1+1;
        else
            nb_acc_H1=nb_acc_H1;
        end
        i = i+j;
        count = 0;
    end
    i=i+1;
end

%ACCELERATIONS HEALTHY_2
% figure, plot(Healthy_2_corr), hold on ,plot(FHR2), hold on, plot(FHR2+15);
nb_acc_H2 = 0;
count = 0;
i=1;
while i <= length(Healthy_2_corr)
    if Healthy_2_corr(i) >= FHR2(i) + 15
        count = count +1;
        j = 1;
        while Healthy_2_corr(i+j) >= FHR2(i+j) + 15
            count = count + 1;
            j = j+1;
        end
        if count >= 30
            nb_acc_H2=nb_acc_H2+1;
        else
            nb_acc_H2=nb_acc_H2;
        end
        i = i+j;
        count = 0;
    end
    i=i+1;
end


%ACCELERATIONS HEALTHY_3
% figure, plot(Healthy_3_corr), hold on ,plot(FHR3), hold on, plot(FHR3+15);
nb_acc_H3 = 0;
count = 0;
i=1;
while i <= length(Healthy_3_corr)
    if Healthy_3_corr(i) >= FHR3(i) + 15
        count = count +1;
        j = 1;
        while Healthy_3_corr(i+j) >= FHR3(i+j) + 15
            count = count + 1;
            j = j+1;
        end
        if count >= 30
            nb_acc_H3=nb_acc_H3+1;
        else
            nb_acc_H3=nb_acc_H3;
        end
        i = i+j;
        count = 0;
    end
    i=i+1;
end

%ACCELERATIONS HEALTHY_4
% figure, plot(Healthy_4_corr), hold on ,plot(FHR4), hold on, plot(FHR4+15);
nb_acc_H4 = 0;
count = 0;
i=1;
while i <= length(Healthy_4_corr)
    if Healthy_4_corr(i) >= FHR4(i) + 15
        count = count +1;
        j = 1;
        while Healthy_4_corr(i+j) >= FHR4(i+j) + 15
            count = count + 1;
            j = j+1;
        end
        if count >= 30
            nb_acc_H4=nb_acc_H4+1;
        else
            nb_acc_H4=nb_acc_H4;
        end
        i = i+j;
        count = 0;
    end
    i=i+1;
end

%ACCELERATIONS HEALTHY_5
% figure, plot(Healthy_5_corr), hold on ,plot(FHR5), hold on, plot(FHR5+15);
nb_acc_H5 = 0;
count = 0;
i=1;
while i <= length(Healthy_5_corr)
    if Healthy_5_corr(i) >= FHR5(i) + 15
        count = count +1;
        j = 1;
        while Healthy_5_corr(i+j) >= FHR5(i+j) + 15
            count = count + 1
            j = j+1;
        end
        if count >= 30
            nb_acc_H5=nb_acc_H5+1
        else
            nb_acc_H5=nb_acc_H5;
        end
        i = i+j;
        count = 0;
    end
    i=i+1;
end

%ACCELERATIONS IUGR_1

% figure, plot(IUGR_1_corr), hold on ,plot(FHR1_IUGR), hold on, plot(FHR1_IUGR+15);
nb_acc_IUGR_1 = 0;
count = 0;
i=1;
while i <= length(IUGR_1_corr)
    if IUGR_1_corr(i) >= FHR1_IUGR(i) + 15
        count = count +1;
        j = 1;
        while IUGR_1_corr(i+j) >= FHR1_IUGR(i+j) + 15
            count = count + 1;
            j = j+1;
        end
        if count >= 30
            nb_acc_IUGR_1=nb_acc_IUGR_1+1;
        else
            nb_acc_IUGR_1=nb_acc_IUGR_1;
        end
        i = i+j;
        count = 0;
    end
    i=i+1;
end

%ACCELERATIONS IUGR_2

% figure, plot(IUGR_2_corr), hold on ,plot(FHR2_IUGR), hold on, plot(FHR2_IUGR+15);
nb_acc_IUGR_2 = 0;
count = 0;
i=1;
while i <= length(IUGR_2_corr)
    if IUGR_2_corr(i) >= FHR2_IUGR(i) + 15
        count = count +1;
        j = 1;
        while IUGR_2_corr(i+j) >= FHR2_IUGR(i+j) + 15
            count = count + 1;
            j = j+1;
        end
        if count >= 30
            nb_acc_IUGR_2=nb_acc_IUGR_2+1;
        else
            nb_acc_IUGR_2=nb_acc_IUGR_2;
        end
        i = i+j;
        count = 0;
    end
    i=i+1;
end


%ACCELERATIONS IUGR_3
% figure, plot(IUGR_3_corr), hold on ,plot(FHR3_IUGR), hold on, plot(FHR3_IUGR+15);
nb_acc_IUGR_3 = 0;
count = 0;
i=1;
while i <= length(IUGR_3_corr)
    if IUGR_3_corr(i) >= FHR3_IUGR(i) + 15
        count = count +1;
        j = 1;
        while IUGR_3_corr(i+j) >= FHR3_IUGR(i+j) + 15
            count = count + 1;
            j = j+1;
        end
        if count >= 30
            nb_acc_IUGR_3=nb_acc_IUGR_3+1;
        else
            nb_acc_IUGR_3=nb_acc_IUGR_3;
        end
        i = i+j;
        count = 0;
    end
    i=i+1;
end
%ACCELERATIONS IUGR_4

% figure, plot(IUGR_4_corr), hold on ,plot(FHR4_IUGR), hold on, plot(FHR4_IUGR+15);
nb_acc_IUGR_4 = 0;
count = 0;
i=1;
while i <= length(IUGR_4_corr)
    if IUGR_4_corr(i) >= FHR4_IUGR(i) + 15
        count = count +1;
        j = 1;
        while IUGR_4_corr(i+j) >= FHR4_IUGR(i+j) + 15
            count = count + 1;
            j = j+1;
        end
        if count >= 30
            nb_acc_IUGR_4=nb_acc_IUGR_4+1;
        else
            nb_acc_IUGR_4=nb_acc_IUGR_4;
        end
        i = i+j;
        count = 0;
    end
    i=i+1;
end
%ACCELERATIONS IUGR_5

% figure, plot(IUGR_5_corr), hold on ,plot(FHR5_IUGR), hold on, plot(FHR5_IUGR+15);
nb_acc_IUGR_5 = 0;
count = 0;
i=1;
while i <= length(IUGR_5_corr)
    if IUGR_5_corr(i) >= FHR5_IUGR(i) + 15
        count = count +1;
        j = 1;
        while IUGR_5_corr(i+j) >= FHR5_IUGR(i+j) + 15
            count = count + 1;
            j = j+1;
        end
        if count >= 30
            nb_acc_IUGR_5=nb_acc_IUGR_5+1;
        else
            nb_acc_IUGR_5=nb_acc_IUGR_5;
        end
        i = i+j;
        count = 0;
    end
    i=i+1;
end

Accelerations = table([nb_acc_H1;nb_acc_H2;nb_acc_H3;nb_acc_H4;nb_acc_H5],[nb_acc_IUGR_1;nb_acc_IUGR_2;nb_acc_IUGR_3;nb_acc_IUGR_4;nb_acc_IUGR_5],'VariableNames',{'NbAccHealthy','NbAccIUGR'});

%% STV calculation
%Compute STV as SD index with a sliding window of 30 samples with overlaping of 29

%STV for H1
wd=30,STD_H1=zeros(windo-wd,length(H1_de_float(1,:))), count=0;
for i=1:length(H1_de_float(1,:))
    for j=1:windo-30
       count=count+1;
       STD_H1(count,i)=std(H1_de_float(j:j+wd-1,i),'omitnan');
    end
    count=0;
end


%H2
wd=30,STD_H2=zeros(windo-wd,length(H2_de_float(1,:))), count=0;
for i=1:length(H2_de_float(1,:))
    for j=1:windo-30
       count=count+1;
       STD_H2(count,i)=std(H2_de_float(j:j+wd-1,i),'omitnan');
    end
    count=0;
end

%H3
wd=30,STD_H3=zeros(windo-wd,length(H3_de_float(1,:))), count=0;
for i=1:length(H3_de_float(1,:))
    for j=1:windo-30
       count=count+1;
       STD_H3(count,i)=std(H3_de_float(j:j+wd-1,i),'omitnan');
    end
    count=0;
end

%H4
wd=30,STD_H4=zeros(windo-wd,length(H4_de_float(1,:))), count=0;
for i=1:length(H4_de_float(1,:))
    for j=1:windo-30
       count=count+1;
       STD_H4(count,i)=std(H4_de_float(j:j+wd-1,i),'omitnan');
    end
    count=0;
end


%H5
wd=30,STD_H5=zeros(windo-wd,length(H5_de_float(1,:))), count=0;
for i=1:length(H5_de_float(1,:))
    for j=1:windo-30
       count=count+1;
       STD_H5(count,i)=std(H5_de_float(j:j+wd-1,i),'omitnan');
    end
    count=0;
end


%I1
wd=30,STD_I1=zeros(windo-wd,length(I1_de_float(1,:))), count=0;
for i=1:length(I1_de_float(1,:))
    for j=1:windo-30
       count=count+1;
       STD_I1(count,i)=std(I1_de_float(j:j+wd-1,i),'omitnan');
    end
    count=0;
end

%I2
wd=30,STD_I2=zeros(windo-wd,length(I2_de_float(1,:))), count=0;
for i=1:length(I2_de_float(1,:))
    for j=1:windo-30
       count=count+1;
       STD_I2(count,i)=std(I2_de_float(j:j+wd-1,i),'omitnan');
    end
    count=0;
end


%I3
wd=30,STD_I3=zeros(windo-wd,length(I3_de_float(1,:))), count=0;
for i=1:length(I3_de_float(1,:))
    for j=1:windo-30
       count=count+1;
       STD_I3(count,i)=std(I3_de_float(j:j+wd-1,i),'omitnan');
    end
    count=0;
end


%I4
wd=30,STD_I4=zeros(windo-wd,length(I4_de_float(1,:))), count=0;
for i=1:length(I4_de_float(1,:))
    for j=1:windo-30
       count=count+1;
       STD_I4(count,i)=std(I4_de_float(j:j+wd-1,i),'omitnan');
    end
    count=0;
end


%I5
wd=30,STD_I5=zeros(windo-wd,length(I5_de_float(1,:))), count=0;
for i=1:length(I5_de_float(1,:))
    for j=1:windo-30
       count=count+1;
       STD_I5(count,i)=std(I5_de_float(j:j+wd-1,i),'omitnan');
    end
    count=0;
end


STV_as_SD_H1=mean(mean(STD_H1, 'omitnan'), 'omitnan');
fprintf('STV index of Healthy_1 is: %f bpm\n', STV_as_SD_H1);
STV_as_SD_H2=mean(mean(STD_H2, 'omitnan'), 'omitnan');
fprintf('STV index of Healthy_2 is: %f bpm\n', STV_as_SD_H2);
STV_as_SD_H3=mean(mean(STD_H3, 'omitnan'), 'omitnan');
fprintf('STV index of Healthy_3 is: %f bpm\n', STV_as_SD_H3);
STV_as_SD_H4=mean(mean(STD_H4, 'omitnan'), 'omitnan');
fprintf('STV index of Healthy_4 is: %f bpm\n', STV_as_SD_H4);
STV_as_SD_H5=mean(mean(STD_H5, 'omitnan'), 'omitnan');
fprintf('STV index of Healthy_5 is: %f bpm\n', STV_as_SD_H5)

STV_as_SD_I1=mean(mean(STD_I1, 'omitnan'), 'omitnan');
fprintf('STV index of IUGR_1 is: %f bpm\n', STV_as_SD_I1);
STV_as_SD_I2=mean(mean(STD_I2, 'omitnan'), 'omitnan');
fprintf('STV index of IUGR_2 is: %f bpm\n', STV_as_SD_I2);
STV_as_SD_I3=mean(mean(STD_I3, 'omitnan'), 'omitnan');
fprintf('STV index of IUGR_3 is: %f bpm\n', STV_as_SD_I3);
STV_as_SD_I4=mean(mean(STD_I4, 'omitnan'), 'omitnan');
fprintf('STV index of IUGR_4 is: %f bpm\n', STV_as_SD_I4);
STV_as_SD_I5=mean(mean(STD_I5, 'omitnan'), 'omitnan');
fprintf('STV index of IUGR_5 is: %f bpm\n', STV_as_SD_I5);
STV_for_Healthy_bpm=[ STV_as_SD_H1; STV_as_SD_H2; STV_as_SD_H3 ;STV_as_SD_H4; STV_as_SD_H5];
STV_for_IUGR_bpm=[ STV_as_SD_I1; STV_as_SD_I2; STV_as_SD_I3 ;STV_as_SD_I4; STV_as_SD_I5];

Names_H={'Healhy_1';'Healthy_2';'Healhy_3';'Healhy_4'; 'Healhy_5'};
IUGR={'IUGR_1';'IUGR_2';'IUGR_3';'IUGR_4'; 'IUGR_5'};


STV_mean_std=[mean(STV_for_Healthy_bpm) std(STV_for_Healthy_bpm); mean(STV_for_IUGR_bpm) std(STV_for_IUGR_bpm)];
T_STV_mean=table(STV_mean_std, 'RowNames', Row);

