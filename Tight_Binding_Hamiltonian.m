%%%%%%%%%%%%%%%%%%%%%%%计算体系局域电流%%%%%%%%%%%%%%%%%%%%
% 同时会保存局域电流为DFTB的形式
% ------------------------参数--------------------------%
% nx 纳米带宽度
% ny 纳米带长度
% t Hopping系数
% ------------------------返回--------------------------%
% coordinatesX 原子X坐标
% coordinatesY 原子Y坐标
% H0 基础哈密顿量矩阵
% Hv 自旋轨道耦合项
% HvHD Haldane项 描述电子在不同晶格子之间进行跃迁时的相位差
% H00 磁场作用 上自旋
% H01 磁场作用 下自旋

function [coordinatesX, coordinatesY, H0, Hv, HvHD, H00, H01] = Tight_Binding_Hamiltonian(nx, ny, t)
    % 自旋轨道耦合    
    ldSO = 0.1 * t;
    % 光场项
    ldHD = -0.02 * t;      

    % 原胞
    x0 = [sqrt(3) / 2, 0, 0, sqrt(3) / 2];
    y0 = [-0.5, 0, 1, 1.5];

    % 扩胞成nx*ny的纳米带
    coordinatesX = [];
    coordinatesY = [];
    for i = 1:nx
        for j = 1:ny
            coordinatesX = [coordinatesX, x0 + sqrt(3) * (i - 1)];
            coordinatesY = [coordinatesY, y0 + 3 * (j - 1)];
        end
    end

    
    % 总原子数
    totalAtomNum = nx * ny * 4;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 生成基本哈密顿量矩阵
    H0 = zeros(totalAtomNum);
    Hv = zeros(totalAtomNum);
    HvHD = zeros(totalAtomNum);
    distances = zeros(totalAtomNum);
    for i = 1:totalAtomNum
        for j = 1:totalAtomNum
            % 寻找最近邻和次近邻
            distances(i, j) = sqrt((coordinatesX(i) - coordinatesX(j)) ^ 2 + (coordinatesY(i) - coordinatesY(j)) ^ 2);

            if distances(i, j) > 0.1 && distances(i, j) < 1.1
                H0(i, j) = t;
            elseif distances(i, j) > 1.1 && distances(i, j) < 2.1
                Hv(i, j) = ldSO / (3 * sqrt(3));
                HvHD(i, j) = ldHD / (3 * sqrt(3));
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 自旋轨道耦合(SOC)和Haldane项
    for i = 1:totalAtomNum
        for j = 1:totalAtomNum
            vCoefficient = 0;
            if distances(i, j) > 1.1 && distances(i, j) < 2.1
                for k = 1:totalAtomNum
                    % 寻找一个格点k, i和j都是k的最近邻
                    if distances(i, k) > 0 && distances(i, k) < 1.1 && distances(j, k) > 0 && distances(j, k) < 1.1
                        % 判断相对k来说，j是在i的顺时针方向还是逆时针方向，扩展一维做叉乘来判断
                        crossResult = cross([coordinatesX(j) - coordinatesX(k), coordinatesY(j) - coordinatesY(k) , 0], ...
                            [coordinatesX(i) - coordinatesX(k), coordinatesY(i) - coordinatesY(k) , 0]);
                        if crossResult(3) > 0
                            vCoefficient = 1;
                        else
                            vCoefficient = -1;
                        end
                    end
                end
            end
            Hv(i, j) = Hv(i, j) * vCoefficient * 1i;
            HvHD(i, j) = HvHD(i, j) * vCoefficient * 1i;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AB效应项还有铁磁项
    % 铁磁交换场
    delta = 0.0 * ldSO;   
    % 反铁磁交换场
    deltaAB = 0.0 * ldSO;   
    % AB格点势 同一原胞内不同格点之间相互作用
    dAB = - 0.0 * ldSO;       
    
    deltaR(1:totalAtomNum) = delta;
    deltaABR(1:totalAtomNum) = deltaAB; 
    dABR(1:totalAtomNum) = dAB;
    for i=1:totalAtomNum
        H00(i, i)= deltaR(i) + dABR(i) * (mod(i, 2) - 0.5)*2 + deltaABR(i) * (mod(i, 2) - 0.5) * 2;
        H01(i, i)= -deltaR(i) + dABR(i) * (mod(i, 2) - 0.5)*2 - deltaABR(i) * (mod(i, 2) - 0.5) * 2;
    end
end