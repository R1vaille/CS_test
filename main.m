%   该程序用于验证压缩感知理论（包含了L1最小范数求解和OMP求解）
%   来源：https://blog.csdn.net/Di_Wong/article/details/81191211
%
%
clear all; close all;
%% 产生信号
choice_transform = 1;      % 选择正交基，1为选择DCT变换，0为选择FFT变换，目的是稀疏化
choice_Phi = 0;         %选择测量矩阵，1为部分哈达玛矩阵，0为高斯随机矩阵
%-----------------------利用三角函数生成频域或DCT域离散信号--------------------------
n = 512;
t = [0: n-1];
f = cos(2*pi/256*t) + sin(2*pi/128*t);   % 产生频域稀疏的信号
%-------------------------------信号降采样率-----------------------
n = length(f);
a = 0.2;            %    取原信号的 a%
m = double(int32(a*n));
%--------------------------------------画信号图--------------------------------------
switch choice_transform
    case 1
        ft = dct(f);    % 离散余弦变换
        disp('ft = dct(f)')
    case 0
        ft = fft(f);   
        disp('ft = fft(f)')
end

disp(['信号稀疏度：',num2str(length(find((abs(ft))>0.1)))])
 
figure('name', 'A Tone Time and Frequency Plot');
subplot(2, 1, 1);
plot(f);
xlabel('Time (s)'); 
% ylabel('f(t)');
subplot(2, 1, 2);
 
switch choice_transform
    case 1
        plot(ft)
        disp('plot(ft)')
    case 0
        plot(abs(ft));
        disp('plot(abs(ft))')
end
xlabel('Frequency (Hz)'); 
% ylabel('DCT(f(t))');
%% 产生感知矩阵和稀疏表示矩阵
%--------------------------利用感知矩阵生成测量值---------------------
switch choice_Phi
    case 1
        Phi = PartHadamardMtx(m,n);       % 感知矩阵（测量矩阵）    部分哈达玛矩阵
    case 0
        Phi = sqrt(1/m) * randn(m,n);     % 感知矩阵（测量矩阵）   高斯随机矩阵
end
% Phi =  randn(m,n);    %randn 生成标准正态分布的伪随机数（均值为0，方差为1）
% Phi = rand(m,n);    % rand 生成均匀分布的伪随机数。分布在（0~1）之间
f2 = (Phi * f')';                 % 通过感知矩阵获得测量值
% f2 = f(1:2:n);
switch choice_transform
    case 1
        Psi = dct(eye(n,n));            %离散余弦变换正交基 代码亦可写为Psi = dctmtx(n);
        disp('Psi = dct(eye(n,n));')
    case 0
        Psi = inv(fft(eye(n,n)));     % 傅里叶正变换，频域稀疏正交基（稀疏表示矩阵）
        disp('Psi = inv(fft(eye(n,n)));')
end
A = Phi * Psi;                    % 恢复矩阵 A = Phi * Psi
%%             重建信号
%---------------------使用CVX工具求解L1范数最小值-----------------
cvx_begin;
    variable x(n) complex;
%     variable x(n) ;
    minimize( norm(x,1) );
    subject to
      A*x == f2' ;
cvx_end;
figure;subplot(2,1,2)
switch choice_transform
    case 1
        plot(real(x));
        disp('plot(real(x))')
    case 0
        plot(abs(x));
        disp(' plot(abs(x))')
end
title('Using L1 Norm（Frequency Domain）');
%  ylabel('DCT(f(t))'); xlabel('Frequency (Hz)');
 
switch choice_transform
    case 1
        sig = dct(real(x));
        disp('sig = dct(real(x))')
    case 0
        sig = real(ifft(full(x)));
        disp(' sig = real(ifft(full(x)))')
end
subplot(2,1,1);
plot(f)
hold on;plot(sig);hold off
title('Using L1 Norm (Time Domain)');
% ylabel('f(t)'); xlabel('Time (s)');
legend('Original','Recovery')
%-----------------------------使用OMP算法重建--------------------------
for K = 1:100
    theta = CS_OMP(f2,A,K);
    %     figure;plot(dct(theta));title(['K=',num2str(K)])
    switch choice_transform
        case 1
            re(K) = norm(f'-(dct(theta)));
        case 0
            re(K) = norm(f'-real(ifft(full(theta))));
    end
end
theta = CS_OMP(f2,A,find(re==min(min(re))));
disp(['最佳稀疏度K=',num2str(find(re==min(min(re))))]);
% theta = CS_OMP(f2,A,10);
figure;subplot(2,1,2);
switch choice_transform
    case 1
        plot(theta);
        disp('plot(theta)')
    case 0
        plot(abs(theta));
        disp('plot(abs(theta))')
end
title(['Using OMP(Frequence Domain)  K=',num2str(find(re==min(min(re))))])
switch choice_transform
    case 1
        sig2 = dct(theta);
        disp('sig2 = dct(theta)')
    case 0
        sig2 = real(ifft(full(theta)));
        disp('sig2 = real(ifft(full(theta)))')
end
subplot(2,1,1);plot(f);hold on;
plot(sig2)
hold off;
title(['Using OMP(Time Domain)  K=',num2str(find(re==min(min(re))))]);
legend('Original','Recovery')
%%
