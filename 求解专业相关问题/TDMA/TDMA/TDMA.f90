!****************************************************************************
!
!  PROGRAM: TDMA
!
!  主程序，调用计算、网格、初始条件、边界条件等子程序
!
!****************************************************************************
program TDMA
    use mesh
    use material
    use solution_control
    implicit none
    integer :: Nx, Ny
    
    !确定X、Y方向的网格节点数量
    Nx = width / delta_X + 1
    Ny = height / delta_Y + 1

    !根据采样点位置获得采样点在网格节点中的索引
    spIndex(:,1) = floor(sample_points_X / delta_X + 1)
    spIndex(:,2) = floor(sample_points_Y / delta_Y + 1)

    !开始计算
    call calculate(Nx, Ny)

end program TDMA
