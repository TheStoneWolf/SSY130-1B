M = 5;
dd = dir('dataB/*.mat');
Nfiles = length(dd);

NN = 8000; % Length of Data to analyze and plot
offset = 4000; %Offset into the dataset
m0 = 0;

range = 0.03;
averaging_ele = 15;

fs = 128;

for ff = 1:Nfiles
    load(['dataB/',dd(ff).name]);
    x0 = data(offset+(1:NN),1);

   % p should be here, idk why
    p = 3;
    % Here is a lmfir_diff filter designed (needs to be completed)
    [h1diff, R] = lmfir_diff(@monofun,@monoderfun, p, M, m0);
    h1diff = h1diff(end:-1:1);
    h1diff = h1diff / sum(abs(h1diff));
    y1diff = filter(h1diff,1,x0);
    
    % y1diff = x0; % uncomment this if directly peak-find on ECG trace.
    %subplot(4,2,(ff-1)*2 + 1)
    figure(1)
    plot(y1diff), hold on
    
    % Wavelet Transform
    % continuous wavelet transform
    [c, l] = wavedec(y1diff, 3, 'sym6');
    % layer 3 detail coefficient
    d5 = wrcoef('d', c, l, 'sym6', 3).^2;
    d5 = normalize(d5, "zscore");
    
    % Moving window smoothing
    % Otherwise in some strange case, MATLAB miss peak, wierd
    wsize = round(0.1*fs);
    d5conv=conv(d5, ones(1, wsize)/wsize, 'same');

    % peak finding
    MPH = median(d5conv);
    MPD = round(0.5*fs);
    
    [p, r_indices] = findpeaks(d5conv, 'MinPeakHeight', MPH, ...
        'MinPeakDistance', MPD);
    r_indices = localMin(y1diff, r_indices, round(0.6*fs));

    % Calculate RR intervals
    rr = diff(r_indices);
    % r1 = zeros(length(rr), 1);
    % r2_0 = zeros(length(rr), 1);
    % r2_indices = zeros(length(r_indices), 0);
    % r2_size = 0;

     % Fix outliers
    % aver_array = p(1:averaging_ele);
    % for i = 2:length(p)
    %     p_average = mean(aver_array);
    %     if p_average*(1+range) < p(i) || p_average*(1-range) > p(i)
    %         p(i) = p_average;
    %         new_mean_rr = mean(rr(max(i-6, 1):min(i+5, length(rr)))); 
    %         r1(i-1) = new_mean_rr;
    %         if i ~= length(rr)
    %             rr(i) = new_mean_rr;
    %         end
    %     else
    %         r1(i-1) = rr(i-1);
    %         r2_0(i-1) = rr(i-1); r2_size = r2_size + 1;
    %         r2_indices(i-1) = r_indices(i-1);
    %     end
    %     aver_array(1+mod(i-1, averaging_ele)) = p(i);
    % end

    plot(r_indices,p,'hr')
    hold off

    
    r1 = zeros(length(rr), 1);
    r2 = zeros(length(rr), 1);
    r2_indices = zeros(length(r_indices), 1);
    r2_i = 1;

    aver_array = rr(1:averaging_ele);
    for i = 2:length(rr)
        rr_mean = mean(aver_array);
        aver_array(1+mod(i-1, averaging_ele)) = rr(i);
        if rr(i) > rr_mean*(1+range) || rr(i) < rr_mean*(1-range)
            r1(i) = rr_mean;
        else
            r1(i) = rr(i);
            r2(r2_i) = rr(i);
            r2_indices(r2_i) = r_indices(i);
            r2_i = r2_i+1;
        end
    end
    rr_y = (1./rr)*fs*60;

    figure(2)
    subplot(3,1,1)
    plot(r_indices(1:end-1)/fs, rr_y), axis([0 70 min(rr_y) max(rr_y)]),
    title("Unmodified")
    subplot(3,1,2)
    plot(r_indices(1:end-1)/fs,(1./r1)*fs*60), axis([0 70 min(rr_y) max(rr_y)]),
    title("Replace with moving average")
    subplot(3,1,3)
    indices_used_perc = 100* r2_i ./ length(rr);
    plot(r2_indices(1:r2_i)/fs,(1./r2(1:r2_i))*fs*60), axis([0 70 min(rr_y) max(rr_y)]),
    title("Remove | using " + indices_used_perc + "% which is " + r2_i + "/" + length(rr))
    %break
    pause
end 

function idx=localMin(data, loc, wsize)
    % search for maximum value within a window
    % mid-point is the location
    idx = zeros(size(loc));
    for i=1:size(loc)
        % begin
        idx1 = floor(loc(i)-wsize/2);
        if idx1 < 1
            idx1 = 1;
        end
        % end
        idx2 = ceil(loc(i)+wsize/2);
        if idx2 > length(data)
            idx2 = length(data);
        end
        window = data(idx1:idx2);
        idx(i) = idx1 + find(window == max(window), 1, "first") - 1;
    end
end

function f = monofun(i,m) 
    if i==0
        f = 1;
    elseif i==1
        f = m;
    elseif i>0
        f = m^i;
    else
        error('i must be a positive integer');
    end
end

function fd = monoderfun(i,m) 
    if i==0
        fd = 0;
    elseif i==1
        fd = 1;
    elseif i>1
        fd = i*(m^(i-1));
    else
        error('i must be a positive integer');
    end
end