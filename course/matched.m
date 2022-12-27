clear all;
a = [+1 +1 +1 -1 0 0 0 0 0 0 0 0]';  % код Баркера
N=length(a); % кол-во временных выборок
M=4; % кол-во элементов решетки
t = (0:1:N-1);

snr0 = 0;
SKO = 1;
SKO_ = SKO;
snr_=snr0;
for snr = -5:5
    snr_dB = 10^(snr/10);% перевод из дБ
    A = sqrt(snr_dB);
    teta = pi/3; % угол прихода сигнала по отношению к нормали к апертуре АР
    fi0 = 0; % начальная фаза
    delta_fi = 2*pi*sin(teta)*0.5; %0.5 = d/lambda, где lambda - длина волны, d - расстояние между элементами,
    m = 1;

    for m = m:(M-1)
        Signal1(m) = exp(1i*(fi0+m*delta_fi));
    end

    Signal = [exp(1i*fi0) Signal1];
    S = A*a*Signal;
    
    Acc = zeros(size(S));
    for number_tests = 1:10000 % количество испытаний
        noise = (randn(N, M)+1i*randn(N, M))/sqrt(2); % добавляет белый гауссов шум к векторному сигналу
        x = S + noise; % складываем сигнал и шум
        y = filter(S(4:-1:1),1,x); % на выходе СФ
        Acc = (Acc+y);
    end
    Acc = Acc / 10000;

    plot(t, abs(Acc(1:end, 1)) / A, 'r|-');
    hold on;

    title("Correlation");
end
