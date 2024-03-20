function [coordinatesX, coordinatesY, H0, Hv, HvHD] = tight_binding_hamiltonian(nx, ny, t)
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

    % 生成基本哈密顿量矩阵
    % 总原子数
    totalAtomNum = nx * ny * 4;
    % 仅跳跃能和座位能的哈密顿矩阵
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
            elseif distances(i, j) < 2.1
                Hv(i, j) = ldSO / (3 * sqrt(3));
                HvHD(i, j) = ldHD / (3 * sqrt(3));
            end
        end
    end

    % 自旋轨道耦合(SOC)和Haldane项
    for i = 1:totalAtomNum
        for j = 1:totalAtomNum
            vCoefficient = 0;
            if distances(i, j) > 1.1 && distances(i, j) < 2.1
                for k = 1:totalAtomNum
                    % 寻找一个格点k, i和j都是k的最近邻
                    if distances(i, k) > 0 && distances(i, k) < 1.1 && distances(j, k) > 0 && distances(j, k) < 1.1
                        % 判断相对k来说，j是在i的顺时针方向还是逆时针方向，扩展一维做叉乘
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
            Hv(i, j) = Hv(i, j) * vCoefficient * 1j;
            HvHD(i, j) = HvHD(i, j) * vCoefficient * 1j;
        end
    end
end