workspace;
format long g;
format compact;

X = 0 : 0.05*pi : 2*pi;
x1 = 0 : 0.05*pi : 12.25*pi;

b = 0.4;
Y = exp(-X * b); % Get a vector.  No noise in this Y yet.

windowSize = 10; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;

numReplicates = 6;
s1 = repmat(Y, [1, numReplicates]);
s2 = awgn(s1,10);
s3 = filter(b, a, s2);
s4 = matchedFilter(s1);
normK = max(s4) / max(s1);


%s4 = s4(((2*pi)/0.05):length(s4));
dots = length(s4);
x2 = 0 : 0.05*pi : (dots-1)*0.05*pi;
length(x2)
length(s4)

% Ploting
subplot(2, 1, 1);
hold on;
grid on;

plot(x1, s2);
plot(x1, s3, '-');
scatter(x1, s1);
legend({'y = noised signal','y = filtered signal', 'y = real signal'},'Location','southwest');


subplot(2, 1, 2);
hold on;
grid on;


s1 = [zeros(1,40), s1];
length(x2)
length(s1)
plot(x2, s4/normK);
%plot(x1, s2);
%plot(x1, s3);
scatter(x2, s1);
legend({'y = signal detection graphic'},'Location','southwest');

function [s] =  matchedFilter(s0)
    T = 2*pi;
    X = 0 : 0.05*pi : 2*pi;
    b = 0.4;

    h = exp((T+X) * b);
    
    s = conv(s0, h);
end
