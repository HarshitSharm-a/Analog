function output = MySolutions(varargin);
%
% ETI220 - Integrated A/D and D/A Converters
%
% Solutions for assignment 1
%
% Created by N/N 2008-11-??
% Last updated by N/N 2008-12-??
%
% Example: (to run exercise 1)
% >> Ass1_Solutions('x',1)
%

% Parameter default values

ex      = 1;             % Run first exercise by default
x       = [];            % Empty signal vector
fin     = 9.97e6;        % Signal frequency
fs      = 81.92e6;       % Sampling frequency
Nx      = 8192;          % FFT length
means   = 16;            % Number of FFTs to average
len     = Nx*means;      % Total signal length
win     = 'rect';        % Desired windowing function
R       = 10;            % Converter resolution
tjit    = 0e-12;         % Std deviation for gaussian jitter to sampling moment
A1      = 1;             % ADC input signal amplitude
Anfl    = 1e-10;         % Noise floor a bit above MATLAB rounding noise
Vref    = 1;             % Reference voltage (single ended; range is from -Vref to Vref)
delta   = 2*Vref/(2^R);  % A quantization step
Arndn   = delta/sqrt(12);% Sets the quantization noise level corresponding
                         % to R bit rectangular quantization noise	   
k2      = 0.000;         % Second order nonlinearity
k3      = 0.000;         % Third order nonlinearity
k4      = 0.000;         % Fourth order nonlinearity
k5      = 0.000;         % Fifth order nonlinearity

jit_gaus = [1e-12 5e-12 15e-12]; % clock jitter 

% Analyse input arguments
index = 1;
while index <= nargin    
    switch (lower(varargin{index}))
    case {'exercise' 'ex' 'x' 'nr' 'ovn' 'number'}
        ex = varargin{index+1};
        index = index+2;
    otherwise
        index=index+1;
    end
end


% Exercise 1
if ex==1
  % Set rectangular window
  win   = 'rect';
  % Generate and sample signal
  x     = sampling('signal','sine','fin',fin,'fs',fs,'ain',A1,'samples',len);
  % Make FFT
  spec  = adcfft('d',x.data,'skip',1,'mean','N',Nx,'w',win);  
  % Plot
  figure('Name','Rectangular Window','NumberTitle','off'); clf;
  plot(0:length(spec)-1,20*log10(abs(spec)),'k-')
  xlabel('Number of bins', 'FontSize',12,'FontWeight','bold','Color','r')
  title('Rectangular Window')
  %Signal peak is at -3.01
end

% Exercise 2
if ex==2
  % Set hann window
  win   = 'hann1';
  % Generate and sample signal
  x     = sampling('signal','sine','fin',fin,'fs',fs,'ain',A1,'samples',len);
  % Make FFT
  spec  = adcfft('d',x.data,'skip',1,'mean','N',Nx,'w',win);  
  % Plot
%   ax1 = subplot(2,1,1);  
  figure('Name','Hann Window','NumberTitle','off'); clf;
  plot(0:length(spec)-1, 20*log10(abs(spec)), 'k-')
  xlabel('Number of bins', 'FontSize',12,'FontWeight','bold','Color','r')
  ylabel('SNQR', 'FontSize',12,'FontWeight','bold','Color','r')
  title('Hann Window')
  
%   ax2 = subplot(2,1,2);
%   plot(ax2, 0:length(spec)-1,20*log10(abs(spec)),'k-')
%   title(ax2, 'Hann Window')
  %Signal peak is at -9.031
end

% Exercise 3
if ex==3
  % Set rectangular window
  win   = 'rect';
  % Generate and sample signal
  x     = sampling('signal','sine','fin',fin,'fs',80e6,'ain',A1,'samples',len);
  % Make FFT
  spec  = adcfft('d',x.data,'skip',1,'mean','N',Nx,'w',win);  
  % Plot
  figure(5); 
  plot(0:length(spec)-1,20*log10(abs(spec)),'k-')
end

% Exercise 4
if ex==4
  % Set rectangular window
  win   = 'rect';
  % Generate and sample signal
  x     = sampling('signal','sine','fin',fin,'fs',80e6,'ain',A1,'samples',len);
  % Make FFT
  spec  = adcfft('d',x.data,'skip',1,'mean','N',Nx,'w',win);  
  % Plot
  figure(5); clf;
  plot(0:length(spec)-1,20*log10(abs(spec)),'k-'); hold on;
  % Set hann window
  win   = 'hann1';
  % Generate and sample signal
  x     = sampling('signal','sine','fin',fin,'fs',80e6,'ain',A1,'samples',len);
  % Make FFT
  spec  = adcfft('d',x.data,'skip',1,'mean','N',Nx,'w',win);  
  % Plot
  % figure(6); 
  plot(0:length(spec)-1,20*log10(abs(spec)),'r-'); hold on;
  % Set hann window
  win   = 'hann2';
  % Generate and sample signal
  x     = sampling('signal','sine','fin',fin,'fs',80e6,'ain',A1,'samples',len);
  % Make FFT
  spec  = adcfft('d',x.data,'skip',1,'mean','N',Nx,'w',win);  
  % Plot
  % figure(6); 
  plot(0:length(spec)-1,20*log10(abs(spec)),'b-')
end

% Exercise 5
if ex==5
  % Set rectangular window
  win   = 'rect';
  % Generate and sample signal with 64 averages
  x     = sampling('signal','sine','fin',fin,'fs',81.92e6,'ain',A1,'samples',Nx);
  % Adding noise
  noise = randn(size(x.data))*(Arndn);
  data_noise = x.data + noise;
  % Make FFT
  spec  = adcfft('d',data_noise,'skip',1,'mean','N',Nx,'w',win);
  % Calculate SNR and SNDR
  perf = adcperf('data',spec,'snr','sndr','w',win);
  msgbox(sprintf('SNR = %d',perf.snr))
  msgbox(sprintf('SNDR = %d',perf.sndr))
  
  % Plot
  clf;
%   ax1 = subplot(2,1,1); 
  plot(0:length(spec)-1,20*log10(abs(spec)),'k-')
  title('FFT with 64 averages')
  
  % Generate and sample signal with zero averages
  x_n     = sampling('signal','sine','fin',fin,'fs',fs,'ain',A1,'samples',Nx);
  % Adding noise
  noise = randn(size(x_n.data))*(Arndn/4);
  data_noise2 = x_n.data + noise;
  % Make FFT
  spec  = adcfft('d',data_noise2,'skip',1,'N',Nx,'w',win);  
  perf = adcperf('data',spec,'snr','sndr','w',win);
  disp(perf.snr)
  disp(perf.sndr)
  % Plot
%   ax2 = subplot(2,1,2); 
%   plot(0:length(spec)-1,20*log10(abs(spec)),'k-')
%   title('FFT with Zero averages')  
end

% Exercise 6
if ex==6
%   % Set hann window
%   win   = 'hann1';
%   % Generate and sample signal
%   x     = sampling('signal','sine','fin',fin,'fs',fs,'ain',A1,'samples',len,'k3',0.01);
%   % Make FFT
%   spec  = adcfft('d',x.data,'skip',1,'mean','N',Nx,'w',win);  
%   % Plot
%   clf;
% %   ax1 = subplot(2,1,1); 
%   plot(0:length(spec)-1, 20*log10(abs(spec)), 'k-')
%   xlabel('Number of bins', 'FontSize',12,'FontWeight','bold','Color','r')
%   title('Hann Window with Fin = 9.97 MHz')
  
  % Set hann window
  win   = 'hann1';
  % Generate and sample signal
  x     = sampling('signal','sine','fin',39.19e6,'fs',fs,'ain',A1,'samples',len,'k3',0.01);
  % Make FFT
  spec  = adcfft('d',x.data,'skip',1,'mean','N',Nx,'w',win);  
  % Plot
%   ax2 = subplot(2,1,2); 
  plot( 0:length(spec)-1, 20*log10(abs(spec)), 'k-')
  xlabel('Number of bins', 'FontSize',12,'FontWeight','bold','Color','r')
  title('Hann Window with Fin = 39.19 MHz')
  
end

% Exercise 7
if ex==7
  % Set hann window
  win   = 'hann1';
  % Cs
  Cs = linspace(0.01e-12,10e-12,50);
  % General 
  kT        = 300*1.38066*1e-23;
  SNR = zeros(1,50); 
  theoSNR = zeros(1,4);
  bit = 8:2:14;
  % Generate and sample signal
  for i = 1:50
      x     = sampling('signal','sine','fin',fin,'fs',fs,'ain',A1,'samples',len,'Cs',Cs(i));
      % Make FFT
      spec  = adcfft('d',x.data,'skip',1,'mean','N',Nx,'w',win);  
      perf = adcperf('data',spec,'snr','w',win);
      SNR(i) = perf.snr;
  end
  % Theoratical SNR for diff bit ADC
  for j = 1:4
      theoSNR(j) = 6.02*bit(j)+1.76;
      disp(theoSNR(j))
  end
  clf;
  plot(Cs,SNR,'b-','DisplayName','SNR vs Cs') 
%   plot(KT_C,SNR)
  hold on
  for i = 1:50
      x     = sampling('signal','sine','fin',fin,'fs',fs,'ain',A1,'samples',len,'Cs',Cs(i));
      KT_C = sqrt(kT/(Cs(i)))*randn(size(x.data));
      % Make FFT
      
      spec  = adcfft('d',x.data+KT_C,'skip',1,'mean','N',Nx,'w',win);  
      perf = adcperf('data',spec,'snr','w',win);
      SNR(i) = perf.snr;
  end
  % Plot
  plot(Cs,SNR,'r-') 
  ylabel('SNR (dB)', 'FontSize',12,'FontWeight','bold','Color','r')
  xlabel('Cs', 'FontSize',12,'FontWeight','bold','Color','r')
  plot(xlim ,[theoSNR(1) theoSNR(1)],'--mo','DisplayName','8-bit ADC')
  plot(xlim ,[theoSNR(2) theoSNR(2)],'--b','DisplayName','10-bit ADC')
  plot(xlim ,[theoSNR(3) theoSNR(3)],'--g','DisplayName','12-bit ADC')
  plot(xlim ,[theoSNR(4) theoSNR(4)],'--v','DisplayName','14-bit ADC')
  title('Hann Window')
  
  hold off
  legend
end

% Exercise 8
if ex==8
  % Set hann window
  win   = 'hann1';
  % General 
  fin = 9.97e6:0.5e6:39.19e6;
  disp(length(fin))
  SNDR_1ps = zeros(1,length(fin)); 
  SNDR_5ps = zeros(1,length(fin)); 
  SNDR_10ps = zeros(1,length(fin)); 
  SNR_1ps = zeros(1,length(fin)); 
  SNR_5ps = zeros(1,length(fin)); 
  SNR_10ps = zeros(1,length(fin));
  % Generate and sample signal and adding clock jitter
  for i = 1:3
      for j = 1:length(fin)
        x     = sampling('signal','sine','fin',fin(j),'fs',fs,'ain',A1,'samples',len,'jit_gaus',jit_gaus(i));
        % Make FFT
        spec  = adcfft('d',x.data,'skip',1,'mean','N',Nx,'w',win);  
        % Calculation of sndr
        perf = adcperf('data',spec,'sndr','w',win);
        if i == 1
         SNDR_1ps(j) = perf.sndr;
         SNR_1ps(j) = -20*log10(2*pi*fin(j)*jit_gaus(i));
        elseif i == 2
         SNDR_5ps(j) = perf.sndr;
         SNR_5ps(j) = -20*log10(2*pi*fin(j)*jit_gaus(i));

        else
         SNDR_10ps(j) = perf.sndr; 
         SNR_10ps(j) = -20*log10(2*pi*fin(j)*jit_gaus(i));

        end
      end
  end
%   plot(fin, SNDR_1ps, 'c-')
%   hold on
%   plot(fin, SNDR_5ps, '-b')
%   plot(fin, SNDR_10ps, '-r')
  
  % Theoretical limitation on the jitter as a functin of input signal
  % frequency. (for graph Jitter degradation of SNR as a function of Input
  % frequency
  % SNR(dBFS) = -20log(2pifin(jitter))
  plot(fin, SNR_1ps, '-g')
  hold on
  plot(fin, SNR_5ps, '-y')
  plot(fin, SNR_10ps, '-c')
end

% Exercise 9
if ex==9
   % Set hann window
   win   = 'hann1';
   fs = 250e6;
   fin = (3277/8192)*fs;
   jit_gaus = 1e-15:.1e-12:20e-12;
   SNDR = zeros(1,length(jit_gaus)); 
   for i = 1:length(jit_gaus)
    x     = sampling('signal','sine','fin',fin,'fs',fs,'ain',A1,'samples',len,'jit_gaus',jit_gaus(i));
    % Make FFT
    spec  = adcfft('d',x.data,'skip',1,'mean','N',Nx,'w',win);
    % Calculation of sndr
    perf = adcperf('data',spec,'sndr','w',win);
    SNDR(i) = perf.sndr;
   end
%    ax1 = subplot(2,1,1); 
   plot(jit_gaus,SNDR,'-r');
   hold on
   fin = 100e6;
   for i = 1:length(jit_gaus)
    x     = sampling('signal','sine','fin',fin,'fs',fs,'ain',A1,'samples',len,'jit_gaus',jit_gaus(i));
    % Make FFT
    spec  = adcfft('d',x.data,'skip',1,'mean','N',Nx,'w',win);
    % Calculation of sndr
    perf = adcperf('data',spec,'sndr','w',win);
    SNDR(i) = perf.sndr;
   end
%    ax2 = subplot(2,1,2);
   plot(jit_gaus,SNDR,'-b');
end

% Exercise 10
if ex==10
  bit = 8;
  Quant_noise = (2/(2^bit))*(1/sqrt(12));
% Set rectangular window
  win   = 'hann1';
  % Generate and sample signal
  x     = sampling('signal','sine','fin',fin,'fs',fs,'ain',A1,'samples',len);
  %Adding noise
  noise = randn(size(x.data))*(Quant_noise);
%   noise = randn(size(x.data))*Arndn;
  data_noise = x.data + noise;  
  % Make FFT
  spec  = adcfft('d',data_noise,'skip',1,'mean','N',Nx,'w',win);  
  %Calculate SNDR
  perf = adcperf('data',spec,'sndr','w',win);
  % Plot
  clf;
  plot(0:length(spec)-1,20*log10(abs(spec)),'k-')
  xlabel('Number of bins', 'FontSize',12,'FontWeight','bold','Color','r')
  title('Hann Window')
  display(perf.sndr)
end

% Exercise 11 part 1 ( Random Noise )
if ex==11 
  % Set rectangular window
  win   = 'hann1';
  % Generate and sample signal
  x     = sampling('signal','sine','fin',fin,'fs',fs,'ain',A1,'samples',len);
  y     = quantization('data',x.data,'R',8,'vref',1); 
  % Make FFT
  spec  = adcfft('d',y,'skip',1,'mean','N',Nx,'w',win);  
  %Calculate SNDR
  perf = adcperf('data',spec,'sndr','w',win);
  % Plot
  figure('Name','Hann Window','NumberTitle','off'); clf;
  plot(0:length(spec)-1,20*log10(abs(spec)),'k-')
  xlabel('Number of bins', 'FontSize',12,'FontWeight','bold','Color','r')
  title('Hann Window')
  display(perf.sndr)
end

% Exercise 11 part 2 ( Quantisation Noise )
if ex==12 
  npow = (Arndn)^2; 
  % General 
  len     = Nx*64;
  % 3 - bit Resolution
  % Set rectangular window
  win   = 'hann1';
  % Generate and sample signal
  x     = sampling('signal','sine','fin',fin,'fs',fs,'ain',A1,'samples',len);
  y     = quantization('data',x.data,'R',8,'vref',1,'npow', npow); 
  % Make FFT
  spec  = adcfft('d',y,'skip',1,'mean','N',Nx,'w',win);  
  %Calculate SNDR
  perf = adcperf('data',spec,'sndr','w',win);
  % Plot
  figure('Name','Hann Window','NumberTitle','off'); clf;
  plot(0:length(spec)-1,20*log10(abs(spec)),'k-')
  xlabel('Number of bins', 'FontSize',12,'FontWeight','bold','Color','r')
  title('Hann Window')
  display(perf.sndr)
end

% Exercise 12
if ex==13
  % 3 - bit Resolution
  res = 3;
  %length = Number of Averages for sin
  leng = 2;
  % Time vector
  t = (1/fs)*(0:len-1);
  disp(length(t)-1)
  s = 2/(length(t)-1);
  ramp =-1:s:1;
%   plot(t,ramp)
  hold on
  % Quantize
  y     = quantization('data',ramp,'res',res,'vref',1); 
  
  Quant_error = ramp - y;
%   plot(t,Quant_error)
  % Generate and sample signal
  t_sine = 0:(1/(150*fin)):(1/fin)*2;

  x = 1*sin(2*pi*fin*t_sine);
  plot (t_sine,x,'-r')
  hold on
  % Quantize
  k     = quantization('data',x,'res',res,'vref',1); 
  error = x - k;
  plot(t_sine,error)
end

% Exercise 14
if ex==14
  bit = 8;
  Quant_noise = (2/(2^(bit+1)))*(1/sqrt(12));
% Set Hann window
  win   = 'hann1';
  % Generate and sample signal
  x     = sampling('signal','sine','fin',fin,'fs',fs,'ain',0.98,'samples',len);
  %Adding quantisation noise with R + 1 bit
  noise = randn(size(x.data))*Quant_noise;
  x_data = x.data + noise;  
  % Make FFT
  spec  = adcfft('d',x_data,'skip',1,'mean','N',Nx,'w',win);  
  %Calculate SNDR
  perf = adcperf('data',spec,'sndr','w',win);
  % Plot
  clf;
  figure('Name','Hann Window','NumberTitle','off'); clf;
  plot(0:length(spec)-1,20*log10(abs(spec)),'k-')
  xlabel('Number of bins', 'FontSize',12,'FontWeight','bold','Color','r')
  title('Hann Window')
  display(perf.sndr)
end

% Exercise 15
if ex==15
  bit = 8;
  Quant_noise = (2/(2^(bit+1)))*(1/sqrt(12));
% Set Hann window
  win   = 'hann1';
  % Generate and sample signal
  x     = sampling('signal','sine','fin',fin,'fs',fs,'ain',1.01,'samples',len);
  %Adding quantisation noise with R + 1 bit
  noise = randn(size(x.data))*Quant_noise;
  x.data = x.data + noise;  
  % Make FFT
  spec  = adcfft('d',x.data,'skip',1,'mean','N',Nx,'w',win);  
  %Calculate SNDR
  perf = adcperf('data',spec,'sndr','w',win);
  % Plot
  clf;
  figure('Name','Hann Window','NumberTitle','off'); clf;
  plot(0:length(spec)-1,20*log10(abs(spec)),'k-')
  xlabel('Number of bins', 'FontSize',12,'FontWeight','bold','Color','r')
  title('Hann Window')
  display(perf.sndr)
end

% Exercise 16
if ex==16
  % Set Hann window
  win   = 'hann1';
  qnpower = 0:1:6;
  % Generate and sample signal
  x     = sampling('signal','sine','fin',fin,'fs',fs,'ain',A1,'samples',len);
  x_data  = x.data + Arndn*randn(size(x.data));
  rand_data = randn(size(x_data));
  for i = 1:length(qnpower)
    rms_noise(i) = qnpower(i)/2;
    data = x_data + Arndn*rand_data*qnpower(i);
    % Make FFT
    spec  = adcfft('d',data,'skip',1,'mean','N',Nx,'w',win);  
    %Calculate SNDR
    perf = adcperf('data',spec,'sndr','w',win);
    sndr(i) = perf.sndr;
%     plot(0:length(spec)-1,20*log10(abs(spec)))
%     hold on
  end
  plot(rms_noise,sndr,'-')
%   for j = 1:7
%       plot(npow,sndr(j))
%       hold on
%       display(sndr(j))
%   end
  %  Plot
%   figure('Name','Hann Window','NumberTitle','off'); clf;
%   plot(0:length(spec)-1,20*log10(abs(spec)),'k-')
%   xlabel('Number of bins', 'FontSize',12,'FontWeight','bold','Color','r')
%   title('Hann Window')
%   display(perf.sndr)
end

% Exercise 17
if ex==17
  % Set rectangular window
  win   = 'rect';
  % data 
  x = dlmread ('3_11_g.txt');
  data = x(:,2);
  % Make FFT
  spec  = adcfft('d',data','skip',1,'mean','N',1024,'w',win);  
  
    %Calculate SNDR
  perf = adcperf('data',spec,'snr','sndr','sfdr','w',win);
  disp(perf.snr)
  disp(perf.sndr)
  disp(perf.sfdr)
  % Plot
  figure('Name','Rectangular Window','NumberTitle','off'); clf;
  plot(0:length(spec)-1,20*log10(abs(spec)),'k-')
  xlabel('Number of bins', 'FontSize',12,'FontWeight','bold','Color','r')
  title('Rectangular Window') 
  
end
 
% Exercise 18
if ex==18
%   t_sine = 0:(1/(4096*fin)):(1/fin)*100;
  t_sine = 1/fs*(0:(8192-1));
  x = 0.98*sin(2*pi*fin*t_sine);
  bit = 8;
  Quant_noise = (2/(2^(bit)))*(1/sqrt(12));
  %Adding quantisation noise with R + 1 bit
  noise = randn(size(x))*Quant_noise;
  x_data = x + noise;  
  hold on
  % Quantize
  k     = quantization('data',x_data,'res',8,'vref',1); 
  error = x - k;
  plot(t_sine,k)
  hold on
  % Make FFT
  spec  = adcfft('d',k,'skip',1,'mean','N',Nx,'w',win);
    perf = adcperf('data',spec,'snr','sndr','sfdr','w',win);
  disp(perf.sndr)
  figure()
  plot(0:length(spec)-1,20*log10(abs(spec)),'-')
end 

% Exercise 18
if ex==19
  fs      = 81.92e6;       % Sampling frequency
   fin     = 9.97e6;        % Signal frequency
   A1      = 1;             % ADC input signal amplitude
   Nx      = 8192;          % FFT length
   means   = 16;            % Number of FFTs to average
   len     = Nx*means;      % Total signal length   
   R       = 10;            % Converter resolution
   Vref    = 1;             % Reference voltage (single ended; range is from -Vref to Vref)
   delta   = 2*Vref/(2^(R));% A quantization step
   Arndn   = delta/sqrt(12);% Sets the quantization noise level corresponding
   % Set hann window
   win     = 'hann1';
   %Sampling
   x
fs      = 81.92e6;       % Sampling frequency
   fin     = 9.97e6;        % Signal frequency
   A1      = 1;             % ADC input signal amplitude
   Nx      = 8192;           % FFT length
   means   = 16;            % Number of FFTs to average
   len     = Nx*means;      % Total signal length   
   R       = 10;            % Converter resolution
   Vref    = 1;             % Reference voltage (single ended; range is from -Vref to Vref)
   delta   = 2*Vref/(2^(R));% A quantization step
   Arndn   = delta/sqrt(12);% Sets the quantization noise level corresponding
   % Set hann window
   win     = 'hann1';
%Sampling
   x       = sampling('signal','sine','fin',fin,'fs',fs,'ain',A1,'samples',len);
   x.data  = x.data + Arndn*randn(size(x.data));
   qnpower = 0:1:6;
   k       = (Arndn);
   z       = randn(size(x.data));
  for i=1:length(qnpower)
   r = qnpower(i)*k*z;
   x.data =x.data + r;
   spec  = adcfft('d',x.data,'skip',1,'mean', means, 'N',Nx,'w',win);
   j{i}= adcperf('data', spec, 'win', win);
   SNDR{i}=j{i}.sndr;
  end
   SNDR=cell2mat(SNDR)
   figure(1); clf;
   plot(qnpower, SNDR, 'k-');
   xlabel('Factor of Quantisation Noise Power');
   ylabel('SNDR (dB)');
   grid on;
end

if ex==20
    a=INLDNL('3_12_d.txt',8,10);
end
