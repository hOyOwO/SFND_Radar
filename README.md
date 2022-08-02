# SFND_Radar

## 1. FMCW Waveform Design
Radar Specifications 
* Frequency of operation = 77GHz
* Max Range = 200m
* Range Resolution = 1 m
* Max Velocity = 100 m/s
* speed of light = 3e8

```matlab

DistanceTarget = 40;
VelocityTarget = 10; %m/s
c = 3e8;
dres = 1; %m
Rmax = 200; %m
RelativeVelocityMax = 100; %m/s
B = c / (2 * dres); %1/s
Tchirp = 5.5 * 2 * Rmax / c; %s
slope = B / Tchirp;
fc= 77e9;  %carrier freq
Nd=128;  % number of chirps
Nr=1024;  %The number of samples on each chirp. 

t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples
% example= linspace (1, 10, 5)
% example = [1, 3.25, 5.5, 7.75, 10]

Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

r_t=zeros(1,length(t));
td=zeros(1,length(t));
```

## 2. Simulation Loop
Simulate Target movement and calculate the beat or mixed signal for every timestamp.
```matlab
for i=1:length(t)         

    r_t(i) = DistanceTarget + VelocityTarget * t(i);
    td(i) = 2 * r_t(i) / c;
    
    Tx(i) = cos(2 * pi * (fc * t(i) + slope*t(i)^2 / 2));
    Rx(i) = cos(2 * pi * (fc * (t(i) - td(i)) + slope*(t(i) - td(i))^2 / 2));
            
    
    Mix(i) = Tx(i).*Rx(i);
    
end
```

## 3. Range FFT (1st FFT)
Implement the Range FFT on the Beat or Mixed Signal and plot the result.

```matlab
Mix = reshape(Mix, [Nr,Nd]);

sig_fft = fft(Mix);
P2 = abs(sig_fft / Nr);
P1 = P2(1:Nr/2+1);

figure ('Name','Range from First FFT')
subplot(2,1,1)

f = 2 * Rmax * (0:(Nr/2)) / Nr;
plot(f, P1); 
xlabel('Range (m)')
ylabel('|P1(f)|')
axis ([0 200 0 1]);
```
<img src = .\img\RangeFFT.jpg>

## 4. RANGE DOPPLER RESPONSE (2D FFT)

```matlab
% Range Doppler Map Generation.

Mix=reshape(Mix,[Nr,Nd]);

sig_fft2 = fft2(Mix,Nr,Nd);
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both dimension
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);
```
<img src = .\img\2DFFT.jpg>

## 5. 2D CFAR
Implement the 2D CFAR process on the output of 2D FFT operation, i.e the Range Doppler Map.
```matlab
%Select the number of Training Cells in both the dimensions.
Tr = 20;
Td = 10;
%Select the number of Guard Cells in both dimensions around the Cell under 
Gr = 10;
Gd = 5;
offset = 5; % offset the threshold by SNR value in dB

%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(Nr/2-(2*Gr+2*Tr+1), Nd-(2*Gd+2*Td+1));


size_T = (2*(Tr+Gr)+1) * (2*(Td+Gd)+1) - (2*Gr+1) * (2*Gd+1);
threshold_CFAR = zeros(Nr/2-(2*Gr+2*Tr+1), Nd-(2*Gd+2*Td+1));
signal_CFAR = zeros(Nr/2, Nd);

for i = 1:(Nr/2-(2*Gr+2*Tr+1))
    for j = 1:(Nd-(2*Gd+2*Td+1))
        noise_T = sum(db2pow(RDM(i:i+2*Tr+2*Gr-1, j:j+2*Td+2*Gd-1)), 'all');
        noise_G = sum(db2pow(RDM(i+Tr:i+Tr+2*Gr-1, j+Td:j+Td+2*Gd-1)), 'all');
        noise_level(i,j) = noise_T - noise_G;
        threshold_CFAR(i,j) = pow2db(noise_level(i,j)/size_T) + offset;

        if RDM(i+Tr+Gr,j+Td+Gd) > threshold_CFAR(i,j)
            signal_CFAR(i+Tr+Gr,j+Td+Gd) = 1;
        end
    end
end

figure,surf(doppler_axis,range_axis,signal_CFAR);
colorbar;
```
<img src = .\img\2DCFAR.jpg>