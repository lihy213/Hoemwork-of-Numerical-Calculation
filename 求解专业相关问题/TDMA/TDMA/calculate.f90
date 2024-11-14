!***********************************************************
!主计算程序，实现二维传热计算
!***********************************************************
subroutine calculate(Nx, Ny)
    use mesh
    use material
    use boundary_condition
    use initial_condition
    use solution_control
    implicit none
    integer, intent(in) :: Nx, Ny
    integer :: i
    real :: current_time, current_sampling_time, delta_Tau
    real, dimension(0:Nx+1,0:Ny+1) :: T, U, T1, delta_T
    real, dimension(Nx,Ny) :: lambda, aP0, aE, aW, aN, aS, aP, Sv, D, P, Q
    real :: kx(Ny), ky(Nx)
    real :: max_delta_T, max_rate
    real :: tStart, tEnd

    !初始化系数矩阵A
    aP0 = 0
    aE = 0
    aW = 0
    aN = 0
    aS = 0

    !初始化温度矩阵
    T = Tf
    T(1:Nx,1:Ny) = T0
    U = Tf
    T1 = Tf

    !初始化源项矩阵
    Sv = S * delta_X * delta_Y
    Sv(1,:) = Sv(1,:) / 2
    Sv(:,1) = Sv(:,1) / 2
    Sv(:,Ny) = Sv(:,Ny) / 2
    
    open(unit=1, file='result.log')
    !初始化输出
    write(1,"(A4,$)") "Time"
    do i = 1, 5
        write(1,"(A14,I1,$)") "Sample_Point", i
    end do
    write(1,"(A15)") "Sample_Point6"

    !按时间步进行温度场计算
    current_time = 0
    current_sampling_time = sampling_interval
    delta_Tau = 0.0001  !考虑到降温刚开始时的温度变化率很大，初始delta_Tau取一个很小的值
    call CPU_TIME(tStart)
    do while (current_time < time_length)
        !更新当前时间
        current_time = current_time + delta_Tau

        !通过多项式函数计算材料导热系数
        lambda = 9.248 + 1.57E-2 * T(1:Nx,1:Ny)

        !计算系数矩阵A
        !(1)内部节点的计算
        aP0 = rho * Cp * delta_X * delta_Y / (delta_Tau / 2)
        aE(1:Nx-1,:) = 2 / (1 / lambda(1:Nx-1,:) + 1 / lambda(2:Nx,:)) / delta_X * delta_Y
        aW(2:Nx,:) = 2 / (1 / lambda(2:Nx,:) + 1 / lambda(1:Nx-1,:)) / delta_X * delta_Y
        aN(:,1:Ny-1) = 2 / (1 / lambda(:,1:Ny-1) + 1 / lambda(:,2:Ny)) / delta_Y * delta_X
        aS(:,2:Ny) = 2 / (1 / lambda(:,2:Ny) + 1 / lambda(:,1:Ny-1)) / delta_Y * delta_X
        !(2)左侧、上下侧对流边界节点的处理
        aW(1,:) = HTC2 * delta_Y
        aN(:,Ny) = HTC1 * delta_X
        aS(:,1) = HTC2 * delta_X
        !(3)右侧对称边界的处理
        aW(Nx,:) = aW(Nx,:) * 2
        !(4)边界节点的减半处理
        !(4-1)上侧aP0、aE、aW减半
        aP0(:,Ny) = aP0(:,Ny) / 2
        aE(:,Ny) = aE(:,Ny) / 2
        aW(:,Ny) = aW(:,Ny) / 2
        !(4-2)下侧aP0、aE、aW减半
        aP0(:,1) = aP0(:,1) / 2
        aE(:,1) = aE(:,1) / 2
        aW(:,1) = aW(:,1) / 2
        !(4-3)左侧aP0、aN、aS减半
        aP0(1,:) = aP0(1,:) / 2
        aN(1,:) = aN(1,:) / 2
        aS(1,:) = aS(1,:) / 2

        !沿x方向使用TDMA隐式求解
        aP = aP0 + aE + aW
        !(1)定义TDMA中的常数矩阵
        D = (aP0 - aN - aS) * T(1:Nx,1:Ny) + aN * T(1:Nx,2:Ny+1) + aS * T(1:Nx,0:Ny-1) + Sv
        !(2)左侧边界的特殊处理
        D(1,:) = D(1,:) + aW(1,:) * Tf
        !!!!!!
        P(1,:) = aE(1,:) / aP(1,:)
        Q(1,:) = D(1,:) / aP(1,:)
        do i = 2, Nx
            kx = aP(i,:) - aW(i,:) * P(i-1,:)
            P(i,:) = aE(i,:) / kx
            Q(i,:) = (D(i,:) + aW(i,:) * Q(i-1,:)) / kx
        end do
        U(Nx,1:Ny) = Q(Nx,:)
        do i = Nx-1, 1, -1
            U(i,1:Ny) = P(i,:) * U(i+1,1:Ny) + Q(i,:)
        end do

        !沿y方向使用TDMA隐式求解
        aP = aP0 + aN + aS
        !(1)定义TDMA中的常数矩阵
        D = (aP0 - aE - aW) * U(1:Nx,1:Ny) + aE * U(2:Nx+1,1:Ny) + aW * U(0:Nx-1,1:Ny) + Sv
        !(2)上下侧边界的特殊处理
        D(:,1) = D(:,1) + aS(:,1) * Tf
        D(:,Ny) = D(:,Ny) + aN(:,Ny) * Tf
        !!!!!!
        P(:,1) = aN(:,1) / aP(:,1)
        Q(:,1) = D(:,1) / aP(:,1)
        do i = 2, Ny
            ky = aP(:,i) - aS(:,i) * P(:,i-1)
            P(:,i) = aN(:,i) / ky
            Q(:,i) = (D(:,i) + aS(:,i) * Q(:,i-1)) / ky
        end do
        T1(1:Nx,Ny) = Q(:,Ny)
        do i = Ny-1, 1, -1
            T1(1:Nx,i) = P(:,i) * T1(1:Nx,i+1) + Q(:,i)
        end do

        !当达到采样时间时，输出采样点温度
        if (current_time >= current_sampling_time) then
            !输出6个采样点的温度值
            write(1,"(I4,$)") floor(current_time / 60)
            do i = 1, 5
                write(1,"(f14.2,$)") T1(spIndex(i,1), spIndex(i,2))
            end do
            write(1,"(f14.2)") T1(spIndex(6,1), spIndex(6,2))
            !
            !更新下一次需要采样的时间
            current_sampling_time = current_sampling_time + sampling_interval
        end if

        !计算最大冷却速率
        delta_T = abs(T1 - T)
        max_delta_T = maxval(delta_T)
        max_rate = max_delta_T / delta_Tau

        !根据对每个时间步内的最大温度变化率限制，以及下一次采样时间确定下一个时间步长
        delta_Tau = min(max_deltaT_per_step / max_rate, current_sampling_time - current_time)

        !更新温度矩阵T
        T = T1

    end do
    call CPU_TIME(tEnd)
    write(1,"('Calculate complete, time consumption: ', f8.4, ' sec')") (tEnd - tStart)
    close(unit=1)
    write(*,*) "计算结束！请在输出文档中查看结果！*^_^*"

end subroutine calculate