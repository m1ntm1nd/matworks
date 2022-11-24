function fallibility;
N=4; % кол-во временных выборок
M=4; % кол-во элементов решетки
a = [+1 +1 +1 -1]'; % код Баркера
snr0 = 0;
SKO = 1;
SKO_ = SKO;
snr_=snr0;
for snr = snr0:1:25

    snr_dB = 10^(snr/10);% перевод в дБ
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
    number_tests = 1;
    Nol = zeros(size(S)); % нулевая матрица размерности S

    for number_tests = number_tests:1000000 % количество испытаний
        noise = (randn(N, M)+1i*randn(N, M))/sqrt(2); % добавляет белый гауссов шум к векторному сигналу
        x = S + noise; % складываем сигнал и шум
        y = filter(A*S(end:-1:1),1,x); % на выходе СФ
        Nol = (Nol+y)./number_tests; % складывем все результаты с фильтра для всех испытаний
    end

    inaccuracy = Nol(:, end);% выводит последний столбец матрицы Nol, это и есть искомый вектор ошибки
    Signal2 = Signal';
    SKO = mean(sqrt((abs(Signal2-inaccuracy)).^2));
    SKO_ = SKO;
    disp('SKO')
    disp(SKO)
    disp('snr')
    disp(snr)
    hold all;
    plot(snr, SKO, 'b.',[snr snr_], [SKO SKO_], 'b');
    snr_ = snr;
    SKO_ = SKO;
end

grid on
xlabel('ОСШ, дБ');
ylabel('СKO');
end
