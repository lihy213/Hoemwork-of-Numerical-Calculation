!***********************************************************
!控制时间步进和采样点等
!***********************************************************
module solution_control
    implicit none
    !物理时间：时间尺度30h，每隔5分钟进行采样点输出
    real, parameter :: time_length = 1800 * 60
    real, parameter :: sampling_interval = 5 * 60
    real, parameter :: max_deltaT_per_step = 10
    real, dimension(6) :: sample_points_X = (/30e-3, 100e-3, 200e-3, 500e-3, 700e-3, 900e-3/)
    real, dimension(6) :: sample_points_Y = (/50e-3, 100e-3, 100e-3, 200e-3, 150e-3, 350e-3/)
    integer :: spIndex(6,2)
end module solution_control