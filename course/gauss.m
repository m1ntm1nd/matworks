workspace;
format long g;
format compact;



SNR = 2; %approx 3 dB

Len = 100;
step = 0.1;
t = linspace(0, Len, Len/step);
s0 = sin(t); %base signal

subplot(3,1,1); %measured base signal + noise
s1 = awgn(s0, SNR, 'measured');

plot(s0, 'b-');
hold on;
plot(s1, 'r-');

subplot(3,1,2); %base signal + noise
s2 = awgn(s0, SNR);

plot(s0, 'b-');
hold on;
plot(s2, 'g-');

subplot(3,1,3); %base signal + manually generated noise
s3 = zeros(length(s0));

for i=1:length(s0)
    s3(i) = s0(i) + 1/sqrt(SNR) * randn(1,1);
end

plot(s0, 'b-');
hold on;
plot(s3, 'r-');

