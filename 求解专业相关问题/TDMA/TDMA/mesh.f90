module mesh
    implicit none
    real, parameter :: width = 1000e-3, height = 400e-3     !关于宽度对称处理（2000*400）
    real, parameter :: delta_X = 1e-2, delta_Y = 4e-3       !划分的网格大小
end module mesh