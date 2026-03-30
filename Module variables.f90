Module variables
    Implicit none
    Real(8)::g=9.81
    Real(8)::Pai=3.1416
    Real(8)::C_Air !气相的波速
    Real(8)::C_Waterhammer !水锤波速
    Real(8)::P_ref=101325 !参考大气压力
    Real(8)::Dens_ref=1000 !参考水体密度
    Real(8)::Omiga_ref
    Real(8)::A_ref  !参考面积
    Real(8)::y_ref         !参考水深
    Real(8)::y_dep         !退压阈值水深
    Real(8)::Q_gas  !用来计算封闭气团压力
    Real(8)::H_gas !gauge air pressure head
    Real(8)::Hb !大气压力水头

    Real(8)::D  !圆形管道直径
    Real(8)::W,Wpre !界面移动速度
    Integer::Inter_type !界面类型 1类和2类 1-正界面； 2-负界面三波型；3-负界面四波型
    Real(8)::int_m(1:4) !2类界面左右侧值
    Real(8)::Pg !气团压力
    Real(8)::Air_gamma=1.4d0
    Real(8)::MLeft,MRight
    Real(8)::DiagTime=0.0d0
    Integer::DiagStep=0,DiagIface=0
    Logical::NegTrackReady=.false.
    Logical::NegUseTrackedRightCharRuntime=.true.
    Logical::EnableInterfaceDebugOutput=.true.
    Logical::EnableNegativeShadowDiagnostics=.false.
    Logical::EnableRuntimeStepPrint=.false.
    Real(8)::ShockT1Deltat=0.0d0
    Real(8)::ShockT1Slope=0.0d0
    Real(8)::ShockT1Manning=0.0d0
    Real(8)::ShockT1Friction=0.0d0
    Integer::NegTrackStep=0,NegTrackIface=0
    Real(8)::NegTrackOm1=0.0d0,NegTrackQ1=0.0d0,NegTrackOm2=0.0d0,NegTrackQ2=0.0d0,NegTrackW=0.0d0
    Integer::NegTrackPreparedStep=0
    Logical,Allocatable::NegTrackReadyMap(:),NegTrackPendingReadyMap(:)
    Integer,Allocatable::NegTrackStepMap(:),NegTrackPendingStepMap(:)
    Real(8),Allocatable::NegTrackOm1Map(:),NegTrackQ1Map(:),NegTrackOm2Map(:),NegTrackQ2Map(:),NegTrackWMap(:)
    Real(8),Allocatable::NegTrackPendingOm1Map(:),NegTrackPendingQ1Map(:),NegTrackPendingOm2Map(:),NegTrackPendingQ2Map(:),NegTrackPendingWMap(:)
    Real(8)::tolerance=1.0E-6 !
    Real(8)::Tol1=1.0E-4
    Real(8)::Tol2=1.0E-14
    Real(8)::Af
End module