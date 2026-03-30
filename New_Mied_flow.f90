!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&&&&&A Modified Shock Tracking-Capturing Finite-Volume Approach for Transient Mixed Flows in Stormwater Systems&&&&&&
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&2023-9-12&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Program Tracking_Capturing_Mixed_Flows

    use variables
    implicit none
    
    Interface
        Function Riemann_Solutions(AL,AR,VL,VR) Result(UM)
            Real(8)::UM(2)
	        Real(8)::F,AM,AL,AR,VR,VL,F_Diff,AM2,HL,HR,CL,CR,QM,HM,ThetaM,ThetaM2
	        Real(8)::PHL,PHR,VM,PHAM,NEWPHA,Phi_Circular,ThetaL,ThetaR
	        Real(8)::F_Equation,F_Equation_Diff,Depth,Celerity,Area,Dichotomization    
        End function
        Function Choice(RanSlope,SL,SR,FL,FR,FLL) Result(Choice_Result)
        	Real(8),dimension(2)::FL,FR,FLL,Choice_Result
	        Real(8)::RanSlope,SL,SR   
        End function
        Function Riemann_Interface(OmigaL,OmigaR,Qm_L,Qm_R) Result(UM)
            Real(8)::OmigaL,OmigaR,Qm_L,Qm_R,HL,HR,VL,VR !HP是我们需要求解到的值，中间态的水深
            Real(8)::UM(2)
            Real(8)::APY,Celerity,Depth,Density
            Real(8)::F1,F2,F3,HP1,HP2,HP3,OmigaP1,OmigaP2,OmigaP3
        End function
    Endinterface
    
    Real(8),Allocatable::Dens(:),Q(:)
    Real(8),Allocatable::H(:),A(:),V(:),C(:),UC(:)
    Real(8),Allocatable::UL(:,:),UR(:,:),FL(:,:),FR(:,:),ULL(:,:),FLL(:,:)
    Real(8),Allocatable::U(:,:),F(:,:),SL(:),SR(:)
    Real(8),Allocatable::UM(:,:),F_Riemann(:,:)
    Logical,Allocatable::InterfaceDiag(:)
    Logical,Allocatable::ShockHybridOverride(:)
    Real(8),Allocatable::ShockHybridU(:,:)
    
    Real(8)::L,Deltat,Time,Smax,Cr,Deltax
    Real(8),External::Area,APY,Dichotomization,Density,Celerity,F_Equation,Hydraulic,Depth
    Real(8),External::Phi_Circular,Phi_Circular_Diff,F_Equation_Diff,SOLU
    Real(8)::OriX,OriY,S0,Manning,Friction,Qm1,Qm2,Rhyd,Theta
    Real(8)::flow_depth
    Real(8)::Initial_inflow_charge
    Real(8)::AirPocketLength0,VentHoleRatio,VentHoleArea,VentCd,LeftHeadGauge,LeftPressHead,InitInterfacePos
    Real(8)::Time_step,Time_putout,kkk
    Real(8)::ML,MR,RanSlope
    Real(8)::ChkOm1,ChkQ1,ChkOm2,ChkQ2,ChkU1,ChkU2,ChkH1,ChkH2
    Real(8)::ChkJ1,ChkJ2,ChkRCharL,ChkRMass,ChkRMom,ChkREne,ChkRW
    Real(8)::DenG,MassG_fix,MassG_Ven,PV_fix !管道内气团的压力，密度,封闭气团固定质量；通风条件下气团质量;气团的固定量
    Real(8)::VolumeG !气团体积
    Real(8)::Pg0,DenG0 !初始条件下的气团压力和密度
    Real(8)::ShockVp1,ShockYp1,ShockCp1,ShockVp2,ShockHp2
    Real(8)::ShockVR,ShockHR,ShockCR,ShockVS,ShockHS,ShockCS,ShockDQ,ShockDN
    Real(8)::ShockSF1,ShockSF2,ShockNewH,ShockNewV,ShockNewA,ShockRhyd
    Integer::n,nc,i,j,x,step,InitIfaceCell
    Integer::Right_boundary,Left_boundary
    Integer::Judg !用于判断是否存在气团，若存在，Judg=1；否则，Judg=0.
    Integer::Ventilation !是否为通风条件，若存在，Ventilation=1；否则，Ventilation=0.
    
    Character(50)::str1,str2,str3,str4,str5,str6,str7,str8
    Character(256):: DATALN  
    Open(1,FILE='INPUT_PipeParameter.txt')
    Open(2,FILE='Mixed Flow.txt')
    Open(3,FILE='Gas Pressure.txt')
    Open(99,FILE='interface_debug.dat')
    
    str1='TITLE="Mixed Flow with TPA"'
    str2='VARIABLES= "X", "Y","Phases"'
    str3='ZONE N='
    str4=', E='
    str5=',ZONETYPE=FEQUADRILATERAL'
    str6='DATAPACKING=BLOCK'
    str7='VARLOCATION=([1-2]=NODAL,[3]=CELLCENTERED)'
    
    include'input.inc'
    
    x=1
    step=0
    Orix=0
    Oriy=0
    Time=0
    MassG_fix=0
    MassG_Ven=0
    Theta=-ATAN(S0)
    ShockT1Slope=S0
    ShockT1Manning=Manning
    ShockT1Friction=Friction
    Af=0.25*Pai*D**2
    VentHoleArea=0.25d0*Pai*(VentHoleRatio*D)**2
    VentCd=0.65d0
    y_ref=0.95*D
    A_ref=Area(y_ref)
    Omiga_ref=Dens_ref*A_ref
    Pg=P_ref
    DenG=Pg/C_Air**2
    nc=n+1
    Deltax=L/n
    Hb=P_ref/Dens_ref/g
    Pg0=Pg
    DenG0=DenG
    LeftHeadGauge=Initial_inflow_charge-Hb
    LeftPressHead=max(D+LeftHeadGauge,D+Tol1)
    InitInterfacePos=max(Tol1,min(L-Tol1,L-AirPocketLength0))
    InitIfaceCell=max(0,min(n,int(InitInterfacePos/Deltax)))
    
    Allocate(Dens(0:nc),Q(0:nc))
    Allocate(U(2,0:nc),F(2,0:nc))
    Allocate(H(0:nc),A(0:nc),V(0:nc),C(0:nc),UC(0:nc))
    Allocate(UM(2,nc),F_Riemann(2,nc))
    Allocate(UL(2,nc),UR(2,nc),FL(2,nc),FR(2,nc),ULL(2,nc),FLL(2,nc),SL(nc),SR(nc))
    Allocate(InterfaceDiag(1:nc))
    Allocate(ShockHybridOverride(1:n),ShockHybridU(2,1:n))
    Allocate(NegTrackReadyMap(1:nc),NegTrackPendingReadyMap(1:nc))
    Allocate(NegTrackStepMap(1:nc),NegTrackPendingStepMap(1:nc))
    Allocate(NegTrackOm1Map(1:nc),NegTrackQ1Map(1:nc),NegTrackOm2Map(1:nc),NegTrackQ2Map(1:nc),NegTrackWMap(1:nc))
    Allocate(NegTrackPendingOm1Map(1:nc),NegTrackPendingQ1Map(1:nc),NegTrackPendingOm2Map(1:nc),NegTrackPendingQ2Map(1:nc),NegTrackPendingWMap(1:nc))
    InterfaceDiag=.false.
    ShockHybridOverride=.false.
    ShockHybridU=0.0d0
    NegTrackReadyMap=.false.
    NegTrackPendingReadyMap=.false.
    NegTrackStepMap=0
    NegTrackPendingStepMap=0
    NegTrackOm1Map=0.0d0
    NegTrackQ1Map=0.0d0
    NegTrackOm2Map=0.0d0
    NegTrackQ2Map=0.0d0
    NegTrackWMap=0.0d0
    NegTrackPendingOm1Map=0.0d0
    NegTrackPendingQ1Map=0.0d0
    NegTrackPendingOm2Map=0.0d0
    NegTrackPendingQ2Map=0.0d0
    NegTrackPendingWMap=0.0d0
    NegTrackPreparedStep=0
    NegTrackReady=.false.
    NegTrackStep=0
    NegTrackIface=0
    NegTrackOm1=0.0d0
    NegTrackQ1=0.0d0
    NegTrackOm2=0.0d0
    NegTrackQ2=0.0d0
    NegTrackW=0.0d0
    
    Include 'initial condition.inc'
    
    Write(2,*) str1
    Write(2,*) str2
    Write(3,'(A)') 'Time_s Pg_kPa_minus100'
    if (EnableInterfaceDebugOutput .or. EnableNegativeShadowDiagnostics) then
        Write(99,'(A)') 'PRE x time i type Wpre W U1L U2L U1R U2R HL HR VL VR CR VRmCR int1 int2 int3 int4 UM1 UM2'
        Write(99,'(A)') 'POST x time i U1L U2L U1R U2R U1RR U2RR HL HR HRR VL VR VRR CL CR CRR'
        Write(99,'(A)') 'ALERT x time i HR HRR U1R U2R'
        Write(99,'(A)') 'CHK3 x time i W Wpre Om1 Q1 Om2 Q2 U1 U2 H1 H2 RCharL RMass RMom REne RW'
        Write(99,'(A)') 'M2WTRACK x time i Wraw Wtrack TrackReady RawStatus RawBest RawSamples RawBrackets TrackStatus TrackBest TrackSamples TrackBrackets TrackIface TrackW'
        Write(99,'(A)') 'M2RCHAR x time i CharRaw CharTrack TrackReady RawStatus RawBest RawSamples RawBrackets TrackStatus TrackBest TrackSamples TrackBrackets TrackIface TrackOm2 TrackU2'
        Write(99,'(A)') 'M2RTRY x time i TrackReady TryStatus CharRaw CharTrack Solved TryW UM1 UM2'
    endif
    
    call RunHybridShockFitting()
    goto 900

    !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&检查是否出现堵塞气团&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    If(Ventilation==0) then !封闭气团
        VolumeG=0
        Do i=1,n,1
            If(H(i)>=y_ref) then
                VolumeG=VolumeG
            Else
                VolumeG=VolumeG+Deltax*(A_ref-A(i))
            Endif
        End do
        PV_fix=Pg*VolumeG
        MassG_fix=VolumeG*DenG
        H_gas=Hb
    End if
    !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&检查是否出现堵塞气团&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    Do
    
        Print*,x,Time,Inter_type
        x=x+1
        
        include'output.inc'
        Do i=0,nc,1
            UC(i)=abs(V(i))+C_Waterhammer
        End do
        
        Smax=Maxval(UC)
        Deltat=Cr*Deltax/Smax
        ShockT1Deltat=Deltat
        Time=Time+Deltat
        
        If(Time>2) exit
        !If(((Pg-P_ref)/1000)>420) exit
        !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&检查是否出现堵塞气团&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        If(Ventilation==0) then !封闭气团
            VolumeG=0
            Do i=1,n,1
                If((H(i)-y_ref)>=0) Then
                    VolumeG=VolumeG
                Else
                    VolumeG=VolumeG+Deltax*(A_ref-A(i))
                Endif
            End do
            DenG=MassG_fix/VolumeG
            Pg=PV_fix/VolumeG
        End if
        !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&检查是否出现堵塞气团&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        InterfaceDiag=.false.
        ShockHybridOverride=.false.
        ShockHybridU=0.0d0
        Do i=0,nc,1
            U(1,i)=A(i)*Dens(i)
            U(2,i)=Q(i)*Dens(i)
            F(1,i)=Q(i)*Dens(i)
            F(2,i)=(Q(i)*Dens(i))**2/(A(i)*Dens(i))+APY(A(i)*Dens(i))
        Enddo
        Do i=1,nc,1
            If((H(i-1)-y_ref)>=0.and.(H(i)-y_ref)>=0) Then !压力流
                UM(1,i)=exp(0.5*(log(U(1,i-1)*U(1,i))+(V(i-1)-V(i))/C_Waterhammer))
                UM(2,i)=UM(1,i)*0.5*(V(i-1)+V(i)+C_Waterhammer*(log(U(1,i-1))-log(U(1,i))))
                F_Riemann(1,i)=UM(2,i)
                F_Riemann(2,i)=UM(2,i)**2/UM(1,i)+APY(UM(1,i))
            Else
                If((H(i-1)-y_ref)<0.and.(H(i)-y_ref)<0) Then !自由液面流
                    UM(1,i)=0.5*(U(1,i-1)+U(1,i))+(U(1,i-1)+U(1,i))/(2*(C(i-1)+C(i)))*(V(i-1)-V(i))
                    UM(2,i)=0.5*UM(1,i)*(V(i-1)+V(i))+0.5*UM(1,i)*(F_Equation(UM(1,i),U(1,i))-F_Equation(UM(1,i),U(1,i-1)))
                    UL(:,i)=U(:,i-1)
                    UR(:,i)=U(:,i)
                    FL(:,i)=F(:,i-1)
                    FR(:,i)=F(:,i)
                    If(UM(1,i)>UL(1,i)) Then
                        ML=sqrt((APY(UM(1,i))-APY(UL(1,i)))*UM(1,i)/(UL(1,i)*(UM(1,i)-UL(1,i))))
                    Else
                        ML=C(i-1)
                    Endif
                    If(UM(1,i)>UR(1,i)) Then                
                        MR=sqrt((APY(UM(1,i))-APY(UR(1,i)))*UM(1,i)/(UR(1,i)*(UM(1,i)-UR(1,i))))                
                    Else
                        MR=C(i)
                    Endif
                    SL(i)=V(i-1)-ML
                    SR(i)=V(i)+MR
                    FLL(:,i)=(SR(i)*FL(:,i)-SL(i)*FR(:,i)+SR(i)*SL(i)*(UR(:,i)-UL(:,i)))/(SR(i)-SL(i))
                    RanSlope=0
                    F_Riemann(:,i)=Choice(RanSlope,SL(i),SR(i),FL(:,i),FR(:,i),FLL(:,i))
                Else !转换界面
                    If((H(i-1)-y_ref)>0.and.(H(i)-y_ref)<0) Then !出现流态过渡,上游为压力流，下游为自由液面流
                        DiagStep=x
                        DiagTime=Time
                        DiagIface=i
                        Wpre=(U(2,i)-U(2,i-1))/(U(1,i)-U(1,i-1))
                        UM(:,i)=Riemann_Interface(U(1,i-1),U(1,i),U(2,i-1),U(2,i))
                        F_Riemann(1,i)=UM(2,i)
                        F_Riemann(2,i)=UM(2,i)**2/UM(1,i)+APY(UM(1,i))
                        InterfaceDiag(i)=EnableInterfaceDebugOutput
                        if (EnableInterfaceDebugOutput) then
                            Write(99,*) 'PRE',x,Time,i,Inter_type,Wpre,W, &
                                U(1,i-1),U(2,i-1),U(1,i),U(2,i),H(i-1),H(i), &
                                V(i-1),V(i),C(i),V(i)-C(i),int_m(1),int_m(2), &
                                int_m(3),int_m(4),UM(1,i),UM(2,i)
                        endif
                        if (EnableInterfaceDebugOutput .and. Inter_type==3 .and. int_m(1)>0.0d0 .and. int_m(3)>0.0d0) then
                            ChkOm1=int_m(1)
                            ChkQ1=int_m(2)
                            ChkOm2=int_m(3)
                            ChkQ2=int_m(4)
                            ChkU1=ChkQ1/ChkOm1
                            ChkU2=ChkQ2/ChkOm2
                            ChkH1=Depth(ChkOm1)
                            ChkH2=Depth(ChkOm2)
                            ChkJ1=ChkOm1*(ChkU1-W)
                            ChkJ2=ChkOm2*(ChkU2-W)
                            ChkRCharL=ChkU1-V(i-1)-C_Waterhammer*LOG(U(1,i-1)/ChkOm1)
                            ChkRMass=ChkJ1-ChkJ2
                            ChkRMom=ChkOm1*(ChkU1-W)**2+APY(ChkOm1)-ChkOm2*(ChkU2-W)**2-APY(ChkOm2)
                            ChkREne=ChkOm1*(ChkU1-W)**2+2.0d0*ChkOm1*g*ChkH1- &
                                ChkOm2*(ChkU2-W)**2-2.0d0*ChkOm2*g*ChkH2
                            ChkRW=W-Wpre
                            Write(99,*) 'CHK3',x,Time,i,W,Wpre,ChkOm1,ChkQ1,ChkOm2,ChkQ2, &
                                ChkU1,ChkU2,ChkH1,ChkH2,ChkRCharL,ChkRMass,ChkRMom,ChkREne,ChkRW
                        endif
                        if (int_m(1)>0.0d0 .and. int_m(3)>0.0d0) then
                            ShockVp2=int_m(2)/int_m(1)
                            ShockHp2=Depth(int_m(1))
                            ShockVp1=int_m(4)/int_m(3)
                            ShockYp1=Depth(int_m(3))
                            ShockCp1=Celerity(ShockYp1)
                            if (i>=2) then
                                ShockVR=V(i-1)+Cr*(V(i-2)-V(i-1))
                                ShockHR=H(i-1)+Cr*(H(i-2)-H(i-1))
                                ShockSF1=Friction*ShockVR*abs(ShockVR)/(2.0d0*D)
                                ShockDQ=ShockVR+g*ShockHR/C_Waterhammer+g*(S0-ShockSF1)*Deltat
                                ShockVS=V(i-1)+Cr*(ShockVp2-V(i-1))
                                ShockHS=H(i-1)+Cr*(ShockHp2-H(i-1))
                                ShockSF2=Friction*ShockVS*abs(ShockVS)/(2.0d0*D)
                                ShockDN=ShockVS-g*ShockHS/C_Waterhammer+g*(S0-ShockSF2)*Deltat
                                ShockNewH=(ShockDQ-ShockDN)/(2.0d0*g/C_Waterhammer)
                                ShockNewV=ShockDN+g*ShockNewH/C_Waterhammer
                                if (ShockNewH>y_ref+Tol1) then
                                    ShockHybridU(1,i-1)=A_ref*Density(ShockNewH)
                                    ShockHybridU(2,i-1)=ShockHybridU(1,i-1)*ShockNewV
                                    ShockHybridOverride(i-1)=.true.
                                endif
                            endif
                            if (i<n) then
                                ShockVR=V(i)+Cr*(ShockVp1-V(i))
                                ShockHR=H(i)+Cr*(ShockYp1-H(i))
                                ShockCR=max(Tol1,C(i)+Cr*(ShockCp1-C(i)))
                                if (ShockHR>Tol1 .and. ShockHR<y_ref-Tol1) then
                                    ShockNewA=Area(ShockHR)
                                    ShockRhyd=Hydraulic(ShockNewA)
                                    if (ShockRhyd>Tol1) then
                                        ShockSF1=Manning**2*ShockVR*abs(ShockVR)/(ShockRhyd**(1.33333d0))
                                        ShockDQ=ShockVR+g*ShockHR/ShockCR+g*(S0-ShockSF1)*Deltat
                                        ShockVS=V(i)+Cr*(V(i+1)-V(i))
                                        ShockHS=H(i)+Cr*(H(i+1)-H(i))
                                        ShockCS=max(Tol1,C(i)+Cr*(C(i+1)-C(i)))
                                        if (ShockHS>Tol1 .and. ShockHS<y_ref-Tol1) then
                                            ShockNewA=Area(ShockHS)
                                            ShockRhyd=Hydraulic(ShockNewA)
                                            if (ShockRhyd>Tol1) then
                                                ShockSF2=Manning**2*ShockVS*abs(ShockVS)/(ShockRhyd**(1.33333d0))
                                                ShockDN=ShockVS-g*ShockHS/ShockCS+g*(S0-ShockSF2)*Deltat
                                                ShockNewH=(ShockDQ-ShockDN)/(g/ShockCR+g/ShockCS)
                                                ShockNewV=ShockDN+g*ShockNewH/ShockCS
                                                if (ShockNewH>Tol1 .and. ShockNewH<y_ref-Tol1) then
                                                    ShockNewA=Area(ShockNewH)
                                                    ShockHybridU(1,i)=Dens_ref*ShockNewA
                                                    ShockHybridU(2,i)=ShockHybridU(1,i)*ShockNewV
                                                    ShockHybridOverride(i)=.true.
                                                endif
                                            endif
                                        endif
                                    endif
                                endif
                            endif
                        endif
                
                    Endif
                Endif
            Endif
        Enddo
        Do i=1,n,1
            U(:,i)=U(:,i)-Deltat/Deltax*(F_Riemann(:,i+1)-F_Riemann(:,i))
            If((H(i-1)-y_ref)>=0.and.(H(i)-y_ref)>=0.and.(H(i+1)-y_ref)>=0) Then !有压流
                A(i)=A_ref
                Dens(i)=U(1,i)/A(i)
                Qm1=U(2,i)
                Qm2=Qm1+0.5*Deltat*(g*U(1,i)*(S0-Friction*Qm1*ABS(Qm1)/(U(1,i)**2*2*D*g)))
                U(2,i)=Qm1+Deltat*(g*U(1,i)*(S0-Friction*Qm2*ABS(Qm2)/(U(1,i)**2*2*D*g)))
            Endif
            If((H(i-1)-y_ref)<0.and.(H(i)-y_ref)<0.and.(H(i+1)-y_ref)<0) Then !自由液面流
                Dens(i)=Dens_ref
                A(i)=U(1,i)/Dens(i)
                Qm1=U(2,i)
                Rhyd=Hydraulic(A(i))
                Qm2=Qm1+0.5*Deltat*(g*Dens_ref*A(i)*(S0-Manning**2*Qm1*ABS(Qm1)/(U(1,i)**2*Rhyd**(1.33333))))
                U(2,i)=Qm1+Deltat*(g*Dens_ref*A(i)*(S0-Manning**2*Qm2*ABS(Qm2)/(U(1,i)**2*Rhyd**(1.33333))))
            Endif
        Enddo
        Do i=1,n,1
            if (ShockHybridOverride(i)) then
                U(:,i)=ShockHybridU(:,i)
            endif
        Enddo
        Do i=1,n,1
            V(i)=U(2,i)/U(1,i)
            H(i)=Depth(U(1,i))
            C(i)=Celerity(H(i))
            Dens(i)=Density(H(i))
            A(i)=Area(H(i))
            Q(i)=V(i)*A(i)
        Enddo
        Do i=1,n,1
            if (EnableInterfaceDebugOutput .and. InterfaceDiag(i)) then
                Write(99,*) 'POST',x,Time,i,U(1,i-1),U(2,i-1),U(1,i), &
                    U(2,i),U(1,i+1),U(2,i+1),H(i-1),H(i),H(i+1), &
                    V(i-1),V(i),V(i+1),C(i-1),C(i),C(i+1)
                if (H(i)<=1.2d-2) then
                    Write(99,*) 'ALERT',x,Time,i,H(i),H(i+1),U(1,i),U(2,i)
                endif
            endif
        Enddo
        
        Include'bc.inc'
        
    End do 
    
900 Continue
    Close(1)
    Close(2)
    Close(3)
    Close(99)
    
Contains
    Include 'hybrid_shock_fitting.inc'

End program
