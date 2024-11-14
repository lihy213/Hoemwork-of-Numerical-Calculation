!****************************************************************************
!
!  PROGRAM: TDMA
!
!  �����򣬵��ü��㡢���񡢳�ʼ�������߽��������ӳ���
!
!****************************************************************************
program TDMA
    use mesh
    use material
    use solution_control
    implicit none
    integer :: Nx, Ny
    
    !ȷ��X��Y���������ڵ�����
    Nx = width / delta_X + 1
    Ny = height / delta_Y + 1

    !���ݲ�����λ�û�ò�����������ڵ��е�����
    spIndex(:,1) = floor(sample_points_X / delta_X + 1)
    spIndex(:,2) = floor(sample_points_Y / delta_Y + 1)

    !��ʼ����
    call calculate(Nx, Ny)

end program TDMA
