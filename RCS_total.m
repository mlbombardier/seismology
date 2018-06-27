%  This script uses irisFetch to collect data from IRIS. A FIR filter is used to lowpass
%  the data in order to plot. A spectrogram is made from the unfiltered data.

ct = ChannelTag('UW','RCS','*','EHZ');
tr = irisFetch.Traces(ct.network,ct.station,ct.location,ct.channel,'2014-11-25 20:00:00','2014-11-25 23:59:59');

data = tr.data;
fSamp = tr.sampleRate;
fNy = fSamp/2.0;
fhi = 20;

% create filter parameters
window = ceil(50.0*fNy);
step = ceil(10.0*fNy);
t = linspace(0.0, size(data,1)/fSamp, size(data,1));

% create filter and apply it to the data
k = fir1(256, fhi/fNy, 'low');
dfil = filter(k, 1, data);

% SPECTROGRAM
figure('Position',[0,400,1500,350])
spectrogram(data,window,window-step,2^nextpow2(window),fSamp,'yaxis')
colormap(jet)
at = gca;
as = gca;
at.XMinorTick = 'on';
as.YMinorTick = 'on';

% PLOT
figure('Position',[0,100,1500,350]);
plot(t,dfil,'b-')
ax = gca;
ay = gca;
ax.XLabel.String = 'Time (Seconds)';
ay.YLabel.String = 'Sampler Counts';
ax.XLim = [0 size(t,2)/fSamp];
ay.YLim = [-2300 2300];
ax.XMinorTick = 'on';
ay.YMinorTick = 'on';

