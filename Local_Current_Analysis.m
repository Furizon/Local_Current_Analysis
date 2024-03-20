% 使用格林函数计算局域电流，具体计算方法参见
% J.-E. Yang. 《Topological spin–valley filtering effects based on hybrid silicene-like nanoribbons》
%  New J. Phys., 卷 22, 期 5, 页 053034, 5月 2020,
% doi: 10.1088/1367-2630/ab84b4.
clc;
clear;
%%%%%%%%%%%%%%%%%%%%%全局变量%%%%%%%%%%%%%%%%%%%%%
% 设置最大迭代步数
global SELF_ENG_CONVERGE_MAX_STEPS;
SELF_ENG_CONVERGE_MAX_STEPS = 1000;
% 设置收敛限制
global SELF_ENG_CONVERGE_LIMIT;
SELF_ENG_CONVERGE_LIMIT = 1e-6;
% 小虚部
global YITA;
YITA = 0.0001;



%%%%%%%%%%%%%%%%%%%%%主函数部分%%%%%%%%%%%%%%%%%%%%%%%
nx = 12;
ny = 5;
t = 1.0;
e = 0.01;
% 生成紧束缚哈密顿量
[coordinatesX, coordinatesY, H0, Hv, HvHD] = Tight_Binding_Hamiltonian(nx, ny, t);

% 计算局域电流
[LocalCurrentUpL, LocalCurrentUpR, LocalCurrentDownL, LocalCurrentDownR] = Cal_Local_Current(e, H0, Hv, HvHD, YITA);

Plot_Local_Current(coordinatesX(41:end-40), coordinatesY(41:end-40), LocalCurrentUpR);
%%%%%%%%%%%%%%%%%%%%%迭代求计算电极自能%%%%%%%%%%%%%%%%%%%%%
% ------------------------参数--------------------------%
% ee 能量, 可以为数组
% Hc 器件自相关
% H0 电极自相关
% H1 电极互相关
% Vlc左电极和器件互相关
% Vrc右电极和器件互相关
% yita 小虚部取值
% ------------------------返回--------------------------%
% SigmaX 左/右电极自能
function SigmaX = Cal_Electrode_Self_Eng(e, Hc, H0, H1, Vxc, yita)
    % 电极单元原子数
    N = size(H0, 1);
    % 器件单元原子数
    M = size(Hc, 1);
    global SELF_ENG_CONVERGE_MAX_STEPS;
    global SELF_ENG_CONVERGE_LIMIT;
    
    e = reshape(e, 1, []);
    
    % 用于存储迭代过程中的结果
    alphas = zeros(N, N, SELF_ENG_CONVERGE_MAX_STEPS);
    alpha = zeros(N, N, SELF_ENG_CONVERGE_MAX_STEPS);
    beta = zeros(N, N, SELF_ENG_CONVERGE_MAX_STEPS);
    gamma = zeros(N, N, SELF_ENG_CONVERGE_MAX_STEPS);
    % 设置初始值
    alphas(:, :, 1) = H0;
    alpha(:, :, 1) = H0;
    beta(:, :, 1) = (H1');
    gamma(:, :, 1) = H1;

    eLength = numel(e);
    SigmaX = zeros(M, M, eLength);

    for i = 1:eLength
        % 迭代计算系数
        for n = 2:SELF_ENG_CONVERGE_MAX_STEPS
            alphas(:, :, n) = alphas(:, :, n - 1) + beta(:, :, n - 1) * ((e(i) + 1i * yita) * eye(N) - alpha(:, :, n - 1)) \ gamma(:, :, n - 1);
            alpha(:, :, n) = alpha(:, :, n - 1) + gamma(:, :, n - 1) * ((e(i) + 1i * yita) * eye(N) - alpha(:, :, n - 1)) \ beta(:, :, n - 1) + beta(:, :, n - 1) * ((e(i) + 1i * yita) * eye(N) - alpha(:, :, n - 1)) \ gamma(:, :, n - 1);
            beta(:, :, n) = beta(:, :, n - 1) * ((e(i) + 1i * yita) * eye(N) - alpha(:, :, n - 1)) \ beta(:, :, n - 1);
            gamma(:, :, n) = gamma(:, :, n - 1) * ((e(i) + 1i * yita) * eye(N) - alpha(:, :, n - 1)) \ gamma(:, :, n - 1);
            % 判断收敛，若收敛则提前退出
            if sum(abs(alphas(:, :, n) - alphas(:, :, n - 1)), 'all') < SELF_ENG_CONVERGE_LIMIT
                break;
            end
        end
        % 计算表面格林函数
        surfaceG = inv((e(i) + 1j * yita) * eye(N) - alphas(:, :, n));

        % 计算自能 surfaceG \ Vxc = inv(surfaceG) * Vxc
        SigmaX(:, :, i) = Vxc' * surfaceG * Vxc;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%计算体系局域电流%%%%%%%%%%%%%%%%%%%%
% 同时会保存局域电流为DFTB的形式
% ------------------------参数--------------------------%
% e 能量
% H 总哈密顿矩阵
% Hv
% HvHD
% ------------------------返回--------------------------%
% LocalCurrentUpL 从电极L到R上自旋电流
% LocalCurrentUpR 从电极R到L上自旋电流
% LocalCurrentDownL 从电极L到R下自旋电流
% LocalCurrentDownR 从电极R到L下自旋电流
function [LocalCurrentUpL, LocalCurrentUpR, LocalCurrentDownL, LocalCurrentDownR] = Cal_Local_Current(e, H, Hv, HvHD, yita)
    fprintf("***Calculation of Local Current Start***\n");

    % 构建哈密顿量
    HUp = H + Hv + HvHD;
    HDown = H - Hv + HvHD;
    
    
    % 切割矩阵
    HcUp = HUp(41:end-40, 41:end-40);
    HcDown = HDown(41:end-40, 41:end-40);
    Hl0Up = HUp(21:40, 21:40);
    Hl0Down = HDown(21:40, 21:40);
    Hr0Up = HUp(end-39:end-20, end-39:end-20);
    Hr0Down = HDown(end-39:end-20, end-39:end-20);
    Hl1Up = HUp(1:20, 21:40);
    Hl1Down = HDown(1:20, 21:40);
    Hr1Up = HUp(end-19:end, end-39:end-20);
    Hr1Down = HDown(end-19:end, end-39:end-20);
    VlcUp = HUp(21:40, 41:200);
    VlcDown = HDown(21:40, 41:200);
    VrcUp = HUp(201:220, 41:200);
    VrcDown = HDown(201:220, 41:200);

    M = size(HcUp, 1);

    % 计算上自旋电子左右电极自能
    SigmaLUpRetarded = Cal_Electrode_Self_Eng(e, HcUp, Hl0Up, Hl1Up, VlcUp, yita);
    SigmaRUpRetarded = Cal_Electrode_Self_Eng(e, HcUp, Hr0Up, Hr1Up, VrcUp, yita);
    SigmaLUpAdvanced = Cal_Electrode_Self_Eng(e, HcUp, Hl0Up, Hl1Up, VlcUp, -yita);
    SigmaRUpAdvanced = Cal_Electrode_Self_Eng(e, HcUp, Hr0Up, Hr1Up, VrcUp, -yita);

    % 计算下自旋电子左右电极自能
    SigmaLDownRetarded = Cal_Electrode_Self_Eng(e, HcDown, Hl0Down, Hl1Down, VlcDown, yita);
    SigmaRDownRetarded = Cal_Electrode_Self_Eng(e, HcDown, Hr0Down, Hr1Down, VrcDown, yita);
    SigmaLDownAdvanced = Cal_Electrode_Self_Eng(e, HcDown, Hl0Down, Hl1Down, VlcDown, -yita);
    SigmaRDownAdvanced = Cal_Electrode_Self_Eng(e, HcDown, Hr0Down, Hr1Down, VrcDown, -yita);

    % 线宽函数
    GammaLUp = 1i * (SigmaLUpRetarded - SigmaLUpAdvanced);
    GammaRUp = 1i * (SigmaRUpRetarded - SigmaRUpAdvanced);
    GammaLDown = 1i * (SigmaLDownRetarded - SigmaLDownAdvanced);
    GammaRDown = 1i * (SigmaRDownRetarded - SigmaRDownAdvanced);

    % 计算器件的格林函数 GR(A)cUp R-Retarded A-Advanced
    GRcUp = inv((e + 1i * yita) * eye(M) - HcUp - SigmaLUpRetarded - SigmaRUpRetarded);
    GRcDown = inv((e + 1i * yita) * eye(M) - HcDown - SigmaLDownRetarded - SigmaRDownRetarded);
    GAcUp = inv((e - 1i * yita) * eye(M) - HcUp - SigmaLUpAdvanced - SigmaRUpAdvanced);
    GAcDown = inv((e - 1i * yita) * eye(M) - HcDown - SigmaLDownAdvanced - SigmaRDownAdvanced);

    % 计算格林函数G^< GLess 
    % GLessUp(Down)L(R) = GRcUp(Down) * GammaL(R)Up(Down) * GAcUp(Down) 
    % L(R)指电子来源，例如L指从L电极到R电极
    GLessUpL = GRcUp * GammaLUp * GAcUp;
    GLessDownL = GRcDown * GammaLDown * GAcDown;
    GLessUpR = GRcUp * GammaRUp * GAcUp;
    GLessDownR = GRcDown * GammaRDown * GAcDown;

    % 计算局域电流
    LocalCurrentUpL = calculate(HcUp, GLessUpL);
    LocalCurrentUpR = calculate(HcUp, GLessUpR);
    LocalCurrentDownL = calculate(HcDown, GLessDownL);
    LocalCurrentDownR = calculate(HcDown, GLessDownR);
    fprintf("***Calculation of Local Current Finished***\n");
    %创建文件夹保存局域电流
    if exist('LocalCurrent', 'dir') == 0
        % 如果文件夹不存在，则创建一个
        mkdir('LocalCurrent');
    end
    % 保存局域电流
    fprintf("***Saving Local Current To LocalCurrent Folder***\n");
    Save_Local_Current('.\LocalCurrent\lcurrentUp.txt', LocalCurrentUpL);
    Save_Local_Current('.\LocalCurrent\lcurrentDown.txt', LocalCurrentDownL);
    Save_Local_Current('.\LocalCurrent\rcurrentUp.txt', LocalCurrentUpR);
    Save_Local_Current('.\LocalCurrent\rcurrentDown.txt', LocalCurrentDownR);
end

%%%%%%%%%%%%%%%%%%%%%%%局域电流的具体计算%%%%%%%%%%%%%%%%%%%%
% ------------------------参数--------------------------%
% Hc 器件哈密顿矩阵
% GLess G<
% ------------------------返回--------------------------%
% LocalCurrent 计算的局域电流
function LocalCurrent = calculate(Hc, GLess)
    LocalCurrent = zeros(size(Hc));
    for i = 1:size(Hc, 1)
        for j = 1:size(Hc, 2)
            LocalCurrent(i, j) = -1i * (Hc(j, i) * GLess(i, j) - GLess(j, i) * Hc(i, j));
        end
    end
    LocalCurrent = real(LocalCurrent);
end


%%%%%%%%%%%%%%%%%%%%%%%保存局域电流文件%%%%%%%%%%%%%%%%%%%%
% 这里会保存成和dftb类似的格式，可以用jmol查看
% ------------------------参数--------------------------%
% path 保存路径
% LocalCurrent 要保存的局域电流
% ------------------------返回--------------------------%
function Save_Local_Current(path, LocalCurrent)
    fid = fopen(path, 'w');
    LocalCurrentIndex = abs(LocalCurrent);
    [~, LocalCurrentIndex] = sort(LocalCurrentIndex, 2, 'descend');
    for i = 1:size(LocalCurrent, 1)
        fprintf(fid, '%d  ', i);
        for j = 1:size(LocalCurrent, 2)
            if abs(LocalCurrent(i, LocalCurrentIndex(i, j))) < 1e-50
                break;
            end
            fprintf(fid, '%d  ', LocalCurrentIndex(i, j));
            fprintf(fid, '%f  ', LocalCurrent(i, LocalCurrentIndex(i, j)));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
end

function Plot_Local_Current(coordinatesX, coordinatesY, LocalCurrent)
    Inout1x = [];
    Inout1y = [];
    scatter(coordinatesX, coordinatesY);
    hold on;
    for i  = 1:size(coordinatesX, 2)
        Inout1x(i)=0; Inout1y(i)=0;
        for j = 1:size(coordinatesX, 2)
        
            aph=angle(coordinatesX(j)-coordinatesX(i)+1i*(coordinatesY(j)-coordinatesY(i)));
        
            Inout1x(i)=Inout1x(i) + real(LocalCurrent(i,j))*cos(aph);
            Inout1y(i)=Inout1y(i) + real(LocalCurrent(i,j))*sin(aph);
        
        end
    end
    quiver(coordinatesX, coordinatesY, Inout1x, Inout1y);
end