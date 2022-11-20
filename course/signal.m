workspace;
format long g;
format compact;

%Lets assume we have M antennas receiving signal
d = 0.5; %half of wave lentgth between antennas
n = 1; %antenna number
phi = pi / 6; %angle between Oy and transmitter
f = 0.1; %Lets assume base freq 0.1 MHz
A = 1; %Amplitude - constant
M = 3; %Antennas amount
tLen = 100; %Transmitting time
step = 0.1;
t = linspace(0, tLen, tLen/step);
SNR = 2;

s = zeros(M, tLen/step+1);

for i = 1:length(t)
    A = cos(t(i));
    for n = 1:M
        s0 = A * exp(1i*2*pi*d*(n-1)*sin(phi)) * exp(1i*2*pi*f*t(i));
        s(n, i+1) = s0;
    end
end


plot(t, s(1, 1:length(t)), 'b-');
hold on;
plot(t, s(2, 1:length(t)), 'g-');
hold on;
plot(t, s(3, 1:length(t)), 'r-');
hold on;



