# AndroidWearDataProcess
    Matlab

    This is the matlab code for data processing.
    The data is from the project AndroidWearSensorLoggor
    AndroidWearSensorLoggor: https://github.com/lylalala/AndroidWearSensorLoggor.git

    environment: MTALAB
--------------------------------------------------------------------------------
    mainForCalibration.m: 校准加速度器，得到校准矩阵
    mainForIntertialNavigation.m: 进行惯导的调用（所调用函数因工程需要涉暂保密）
    mainForPost.m: 获得相对姿态
                    主功能：第1-5部分获得相对姿态
                    辅助功能：第6部分产生惯导所需的数据，供mainForIntertialNavigation使用
                            第7部分可以进行一些画图验证
                            第8部分可以为Project:inmed2014-master产生原始数据