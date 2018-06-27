%  Madison Bombardier


%  This script uses irisFetch to collect a set of data which it then splits
%  into equal sections of a specified length in order to cross correlate
%  them against a single master event using the seismic data analysis 
%  toolbox, GISMO. There are five locations where manual entry is needed.
%  The master event needs to be the same length as all the other sections.
%  Due to the rounding-down of the nwf variable, the length of the master
%  event should always be too short if it is not the same. The script takes
%  care of that.

%  QUALIFIER
%  This script computes the correlation coefficients and lag times of all
%  sections compared to all other sections. Alas, we throw out almost all
%  of that when we create masterCo, which includes only the coefficients
%  from the comparisons between each section with the master event.



%  ***REQUIRED MANUAL INPUT***
%  choose channel tag fields
%  SETUP ChannelTag(network, station, location, channel)
ct = ChannelTag('UW','RCS','*','EHZ');
%ct = ChannelTag('OO','AXEC2','*','HHZ');

%  ***REQUIRED MANUAL INPUT***
%  choose start and end times for your data
%  SETUP irisFetch.Traces(network, station, location, channel, startDate, endDate)
%  SETUP dates in 'YYYY-MM-DD hh:mm:ss' format
tr = irisFetch.Traces(ct.network,ct.station,ct.location,ct.channel,'2014-11-25 20:00:00','2014-11-25 23:59:59');

%  ***REQUIRED MANUAL INPUT***
%  choose the length of sections to compare
sec_length_sec = 10.0;  % length of sections in seconds

data = tr.data;   % from struct file
ST = tr.startTime;   % from struct file in serial date numbers
ET = tr.endTime;   % from struct file in serial date numbers

nwf = floor(size(data,1)/(tr.sampleRate*sec_length_sec));   % number of waveforms
sec_length_samp = floor(size(data,1)/nwf);   % divide data into correct wf-length sample sections
samptime = (ET-ST)/size(data,1);  % ST+samptime = time at second sample --- time between each sample


for i = 1:size(data,1)
    data(i,2) = (ST+(samptime*i)-samptime);  % attach sample time to each sample
end

trigS = [];
w = waveform();

%  SETUP waveform(ChannelTag, samplerate, starttime, data, units)
%  example: data = data(1:3000,1) for w1, data(3001:6000,1) for w2, etc.
for r = 1:nwf   % only loop for as many wf as we want
    if data(r*sec_length_samp,1) < size(data,1)  % make sure we don't go out of bounds
        trigS(r,1) = data((1+(r-1)*sec_length_samp),2);  % grab trigger time for each wf section
        % create a wf for each section of data!
        w(r,1) = waveform(ct,tr.sampleRate,trigS(r,1),data((1+(r-1)*sec_length_samp):(r*sec_length_samp),1),'m/s');
        w(r,1) = w(r,1)-mean(w(r,1));  % recalibrate the data around zero
    else
    end
end



%  ***REQUIRED MANUAL INPUT***
%  choose the start and end times of a single master event that will be compared
%  to all other wf. Must be the same length as sec_length_sec
master = irisFetch.Traces(ct.network,ct.station,ct.location,ct.channel,'2014-11-25 21:20:44.000','2014-11-25 21:20:54.000');

%  if master is a few samples off of the length it should be, find how many
%  more samples it needs, and put zeros in their place
if length(master.data) < sec_length_samp
    mEvent = master.data;  % include all the data
    diff = sec_length_samp - length(master.data);  % find the diff in length
    mEvent = [mEvent
        zeros(diff,1)];  % put zeros in the empty spots at the end
else
    mEvent = master.data;
end

mEvent = mEvent-mean(mEvent);  % recalibrate the data around zero


%  place master at the end of the list
w(nwf+1,1) = waveform(ct,master.sampleRate,master.startTime,mEvent,'m/s');
%  this will mess with the format of your wf object: "with dissimilar
%  fields". it's fine though

cObj = correlation(w);  % create the correlation object with all your wf

%  ***REQUIRED MANUAL INPUT***
%  SETUP c = butter(correlation object, [lofq hifq])
c = butter(cObj,[2 20]);  % bandpass filter so noise doesn't affect correlation coeffs
%  based on my previous analysis, the pulses exist within this fq range

C = xcorr(c);  % cross correlate all samples to get corr coeffs and lag times
C = adjusttrig(C);  % line up samples according to position of highest corr
corr_matrix = get(C,'corr');  % create matrix of all corr coeffs
masterCo = corr_matrix(:,length(w));  % create matrix of coeffs from correlation with the master only


%  ***PLOTS***

figure('Position',[0,400,1500,350])
plot(data(:,1));  % plot the raw data
title('data');

figure('Position',[0,400,1500,350])
plot(w(length(w),1));  % plot master event
title('master event');

figure('Position',[0,400,1500,350])
plot(masterCo);  % plot corr coeffs
title('coefficients')

%  SPECTROGRAM
fSamp = 1/sec_length_sec;
fNy = fSamp/2.0;

window = ceil(1000*fNy);
step = ceil(200*fNy);

%  SETUP spectrogram(data,window,noverlap,nfft,samplerate)
%  create a spectrogram of the corr coeffs to see fq of pulses
figure('Position',[0,400,1500,350])
spectrogram(masterCo,window,window-step,2^nextpow2(window),fSamp,'yaxis')
colormap(jet)
ax = gca;
as = gca;
ax.XMinorTick = 'on';
as.YMinorTick = 'on';



