!***********************************************************
!���������ʵ�ֶ�ά���ȼ���
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

    !��ʼ��ϵ������A
    aP0 = 0
    aE = 0
    aW = 0
    aN = 0
    aS = 0

    !��ʼ���¶Ⱦ���
    T = Tf
    T(1:Nx,1:Ny) = T0
    U = Tf
    T1 = Tf

    !��ʼ��Դ�����
    Sv = S * delta_X * delta_Y
    Sv(1,:) = Sv(1,:) / 2
    Sv(:,1) = Sv(:,1) / 2
    Sv(:,Ny) = Sv(:,Ny) / 2
    
    open(unit=1, file='result.log')
    !��ʼ�����
    write(1,"(A4,$)") "Time"
    do i = 1, 5
        write(1,"(A14,I1,$)") "Sample_Point", i
    end do
    write(1,"(A15)") "Sample_Point6"

    !��ʱ�䲽�����¶ȳ�����
    current_time = 0
    current_sampling_time = sampling_interval
    delta_Tau = 0.0001  !���ǵ����¸տ�ʼʱ���¶ȱ仯�ʺܴ󣬳�ʼdelta_Tauȡһ����С��ֵ
    call CPU_TIME(tStart)
    do while (current_time < time_length)
        !���µ�ǰʱ��
        current_time = current_time + delta_Tau

        !ͨ������ʽ����������ϵ���ϵ��
        lambda = 9.248 + 1.57E-2 * T(1:Nx,1:Ny)

        !����ϵ������A
        !(1)�ڲ��ڵ�ļ���
        aP0 = rho * Cp * delta_X * delta_Y / (delta_Tau / 2)
        aE(1:Nx-1,:) = 2 / (1 / lambda(1:Nx-1,:) + 1 / lambda(2:Nx,:)) / delta_X * delta_Y
        aW(2:Nx,:) = 2 / (1 / lambda(2:Nx,:) + 1 / lambda(1:Nx-1,:)) / delta_X * delta_Y
        aN(:,1:Ny-1) = 2 / (1 / lambda(:,1:Ny-1) + 1 / lambda(:,2:Ny)) / delta_Y * delta_X
        aS(:,2:Ny) = 2 / (1 / lambda(:,2:Ny) + 1 / lambda(:,1:Ny-1)) / delta_Y * delta_X
        !(2)��ࡢ���²�����߽�ڵ�Ĵ���
        aW(1,:) = HTC2 * delta_Y
        aN(:,Ny) = HTC1 * delta_X
        aS(:,1) = HTC2 * delta_X
        !(3)�Ҳ�ԳƱ߽�Ĵ���
        aW(Nx,:) = aW(Nx,:) * 2
        !(4)�߽�ڵ�ļ��봦��
        !(4-1)�ϲ�aP0��aE��aW����
        aP0(:,Ny) = aP0(:,Ny) / 2
        aE(:,Ny) = aE(:,Ny) / 2
        aW(:,Ny) = aW(:,Ny) / 2
        !(4-2)�²�aP0��aE��aW����
        aP0(:,1) = aP0(:,1) / 2
        aE(:,1) = aE(:,1) / 2
        aW(:,1) = aW(:,1) / 2
        !(4-3)���aP0��aN��aS����
        aP0(1,:) = aP0(1,:) / 2
        aN(1,:) = aN(1,:) / 2
        aS(1,:) = aS(1,:) / 2

        !��x����ʹ��TDMA��ʽ���
        aP = aP0 + aE + aW
        !(1)����TDMA�еĳ�������
        D = (aP0 - aN - aS) * T(1:Nx,1:Ny) + aN * T(1:Nx,2:Ny+1) + aS * T(1:Nx,0:Ny-1) + Sv
        !(2)���߽�����⴦��
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

        !��y����ʹ��TDMA��ʽ���
        aP = aP0 + aN + aS
        !(1)����TDMA�еĳ�������
        D = (aP0 - aE - aW) * U(1:Nx,1:Ny) + aE * U(2:Nx+1,1:Ny) + aW * U(0:Nx-1,1:Ny) + Sv
        !(2)���²�߽�����⴦��
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

        !���ﵽ����ʱ��ʱ������������¶�
        if (current_time >= current_sampling_time) then
            !���6����������¶�ֵ
            write(1,"(I4,$)") floor(current_time / 60)
            do i = 1, 5
                write(1,"(f14.2,$)") T1(spIndex(i,1), spIndex(i,2))
            end do
            write(1,"(f14.2)") T1(spIndex(6,1), spIndex(6,2))
            !
            !������һ����Ҫ������ʱ��
            current_sampling_time = current_sampling_time + sampling_interval
        end if

        !���������ȴ����
        delta_T = abs(T1 - T)
        max_delta_T = maxval(delta_T)
        max_rate = max_delta_T / delta_Tau

        !���ݶ�ÿ��ʱ�䲽�ڵ�����¶ȱ仯�����ƣ��Լ���һ�β���ʱ��ȷ����һ��ʱ�䲽��
        delta_Tau = min(max_deltaT_per_step / max_rate, current_sampling_time - current_time)

        !�����¶Ⱦ���T
        T = T1

    end do
    call CPU_TIME(tEnd)
    write(1,"('Calculate complete, time consumption: ', f8.4, ' sec')") (tEnd - tStart)
    close(unit=1)
    write(*,*) "�����������������ĵ��в鿴�����*^_^*"

end subroutine calculate