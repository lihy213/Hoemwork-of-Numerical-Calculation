module boundary_condition
    implicit none
    !定义对流换热系数
    real, parameter :: HTC1 = 1000, HTC2 = 500
    !定义流体温度，单位：℃
    real, parameter :: Tf = 500
end module boundary_condition