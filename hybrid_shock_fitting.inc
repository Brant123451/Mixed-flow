Subroutine RunHybridShockFitting()
    Implicit None
    Integer::NPShock,NFShock,IfaceIndex,K,RetryLoc,NFActive
    Real(8)::ShockLu,ShockWu,ShockWuNew,ShockDxP,ShockDxF,ShockDxFGeom,ShockDt0,ShockDt1,ShockDt2
    Real(8)::ShockResn,ShockVolumeNew,ShockHGasAbs,ShockMdot,ShockRhoGas,ShockRhoUp,ShockDP,MassGVenNew
    Real(8)::DtTry,VolumeGBase,PgBase,HGasBase,DenGBase,MassGVenBase
    Character(16)::ShockTag
    Logical::ShockOk,StepAccepted
    Real(8),Allocatable::HP(:),VP(:),HNextP(:),VNextP(:)
    Real(8),Allocatable::YF(:),VF(:),YNextF(:),VNextF(:),CF(:)

    Write(99,*) 'HYB_STAGE','ENTER_RUN'
    call InitializeHybridState(NPShock,NFShock,ShockLu,ShockWu,HP,VP,HNextP,VNextP, &
        YF,VF,YNextF,VNextF,CF)
    call ProjectHybridStateToFixedGrid(NPShock,NFShock,ShockLu,HP,VP,YF,VF,CF)

    call ComputeFreeHybridMetrics(NFShock,ShockLu,NFActive,ShockDxFGeom,ShockDxF)
    call UpdateHybridGasVolume(NFShock,NFActive,ShockDxF,YF,VolumeG)
    if (Ventilation==0) then
        H_gas=0.0d0
        Pg=P_ref
    else
        H_gas=0.0d0
        Pg=P_ref
    endif
    DenG=Pg/C_Air**2
    MassG_Ven=max(1.0d-12,DenG*VolumeG)
    Write(3,'(F24.14,1X,F24.14)') Time,Pg/1000.0d0-100.0d0

    do
        call ProjectHybridStateToFixedGrid(NPShock,NFShock,ShockLu,HP,VP,YF,VF,CF)
        if (EnableRuntimeStepPrint) then
            Print*,x,Time,Inter_type,ShockLu,ShockWu
        endif
        include 'output.inc'
        if (Time>=2.0d0) then
            exit
        endif
        if (Time>=Time_putout) exit

        ShockDxP=max(ShockLu/dble(NPShock),Tol1)
        call ComputeFreeHybridMetrics(NFShock,ShockLu,NFActive,ShockDxFGeom,ShockDxF)
        call UpdateFreeCelerity(NFShock,NFActive,YF,CF)
        ShockDt1=Cr*ShockDxP/max(maxval(abs(VP))+C_Waterhammer,Tol1)
        ShockDt2=Cr*ShockDxF/max(maxval(abs(VF(1:NFActive))+CF(1:NFActive)),Tol1)
        if (abs(ShockWu)>Tol1) then
            ShockDt0=0.5d0*min(ShockDxP,ShockDxF)/abs(ShockWu)
        else
            ShockDt0=1.0d30
        endif
        Deltat=min(ShockDt0,min(ShockDt1,ShockDt2))
        Deltat=min(Deltat,2.0d0-Time)
        Deltat=min(Deltat,Time_putout-Time)
        if (Deltat<=Tol2) exit
        StepAccepted=.false.
        DtTry=Deltat
        VolumeGBase=VolumeG
        PgBase=Pg
        HGasBase=H_gas
        DenGBase=DenG
        MassGVenBase=MassG_Ven

        do RetryLoc=1,16,1
            ShockT1Deltat=DtTry
            call AdvancePressurizedHybrid(NPShock,HP,VP,ShockDxP,DtTry,ShockHp2,ShockVp2, &
                HNextP,VNextP)
            call AdvanceFreeHybrid(NFShock,NFActive,YF,VF,CF,ShockDxF,DtTry,ShockYp1,ShockVp1, &
                ShockCp1,YNextF,VNextF)
            call UpdateFreeCelerity(NFShock,NFActive,YNextF,CF)

            ShockVolumeNew=max(Tol1,VolumeGBase-DtTry*Af*VNextP(NPShock))
            if (Ventilation==0) then
                ShockHGasAbs=max(Hb+Tol1,HGasBase+Hb)
                H_gas=ShockHGasAbs*(VolumeGBase/max(ShockVolumeNew,Tol1))**Air_gamma-Hb
                Pg=P_ref+Dens_ref*g*H_gas
                DenG=Pg/C_Air**2
            else
                ShockRhoGas=max(MassGVenBase/max(VolumeGBase,Tol1),1.0d-12)
                ShockDP=PgBase-P_ref
                if (abs(ShockDP)>1.0d-8 .and. VentHoleArea>1.0d-12) then
                    if (ShockDP>=0.0d0) then
                        ShockRhoUp=ShockRhoGas
                    else
                        ShockRhoUp=max(DenG0,1.0d-12)
                    endif
                    ShockMdot=sign(VentCd*VentHoleArea*sqrt(2.0d0*ShockRhoUp*abs(ShockDP)),ShockDP)
                else
                    ShockMdot=0.0d0
                endif
                MassGVenNew=max(1.0d-12,MassGVenBase-DtTry*ShockMdot)
                DenG=max(MassGVenNew/max(ShockVolumeNew,Tol1),1.0d-12)
                Pg=max(1.0d-12,Pg0*(DenG/max(DenG0,1.0d-12))**Air_gamma)
                H_gas=(Pg-P_ref)/(Dens_ref*g)
                MassG_Ven=MassGVenNew
            endif
            VolumeG=ShockVolumeNew

            IfaceIndex=max(1,min(nc,int(ShockLu/Deltax)+1))
            if (IfaceIndex>=n .and. NFActive==1) then
                ShockWuNew=0.0d0
                ShockTag='OUTLET'
                ShockResn=0.0d0
                ShockOk=.true.
                ShockHp2=max(HNextP(NPShock),Tol1)
                ShockVp2=VNextP(NPShock)
                ShockYp1=min(max(YNextF(1),Tol1),y_ref-Tol1)
                ShockVp1=VNextF(1)
                ShockCp1=max(Celerity(ShockYp1),Tol1)
            else
                call ShockTran(IfaceIndex,HNextP(NPShock),VNextP(NPShock),YNextF(1),VNextF(1),CF(1), &
                    ShockWu,H_gas,DtTry,ShockWuNew,ShockTag,ShockResn,ShockOk)
            endif
            if (ShockOk) then
                Deltat=DtTry
                StepAccepted=.true.
                exit
            endif

            VolumeG=VolumeGBase
            Pg=PgBase
            H_gas=HGasBase
            DenG=DenGBase
            MassG_Ven=MassGVenBase
            DtTry=0.5d0*DtTry
            if (DtTry<=Tol2) exit
        enddo

        HP=HNextP
        VP=VNextP
        YF=YNextF
        VF=VNextF
        if (.not.StepAccepted) exit

        Time=Time+Deltat
        ShockLu=max(Tol1,min(L-Tol1,ShockLu+ShockWu*Deltat))
        ShockWu=ShockWuNew
        x=x+1
        Write(3,'(F24.14,1X,F24.14)') Time,Pg/1000.0d0-100.0d0
    enddo

    call ProjectHybridStateToFixedGrid(NPShock,NFShock,ShockLu,HP,VP,YF,VF,CF)
    Write(99,*) 'HYB_STAGE','EXIT_RUN',Time,ShockLu,ShockWu
End Subroutine

Subroutine ComputeFreeHybridMetrics(NFShock,ShockLu,NFActive,ShockDxGeom,ShockDxComp)
    Implicit None
    Integer,Intent(in)::NFShock
    Real(8),Intent(in)::ShockLu
    Integer,Intent(out)::NFActive
    Real(8),Intent(out)::ShockDxGeom,ShockDxComp
    Real(8)::ShockLenF

    ShockLenF=max(L-ShockLu,Tol1)
    if (ShockLenF<2.0d0*Deltax) then
        NFActive=max(1,min(NFShock,ceiling(ShockLenF/max(Deltax,Tol1))))
        ShockDxGeom=max(ShockLenF/dble(NFActive),Tol1)
        ShockDxComp=Deltax
    else
        NFActive=NFShock
        ShockDxGeom=max(ShockLenF/dble(NFShock),Tol1)
        ShockDxComp=ShockDxGeom
    endif
End Subroutine

Subroutine InitializeHybridState(NPShock,NFShock,ShockLu,ShockWu,HP,VP,HNextP,VNextP, &
    YF,VF,YNextF,VNextF,CF)
    Implicit None
    Integer,Intent(out)::NPShock,NFShock
    Real(8),Intent(out)::ShockLu,ShockWu
    Real(8),Allocatable,Intent(out)::HP(:),VP(:),HNextP(:),VNextP(:)
    Real(8),Allocatable,Intent(out)::YF(:),VF(:),YNextF(:),VNextF(:),CF(:)
    Integer::K
    Real(8)::Om1Loc,Om2Loc

    NPShock=1
    do K=1,n,1
        if (H(K)<y_ref) exit
        NPShock=K
    enddo
    NPShock=max(1,min(n-1,NPShock))
    NFShock=n-NPShock

    Allocate(HP(NPShock),VP(NPShock),HNextP(NPShock),VNextP(NPShock))
    Allocate(YF(NFShock),VF(NFShock),YNextF(NFShock),VNextF(NFShock),CF(NFShock))

    do K=1,NPShock,1
        HP(K)=max(H(K),Tol1)
        VP(K)=V(K)
    enddo
    do K=1,NFShock,1
        YF(K)=min(max(H(NPShock+K),Tol1),y_ref-Tol1)
        VF(K)=V(NPShock+K)
    enddo
    call UpdateFreeCelerity(NFShock,NFShock,YF,CF)

    ShockLu=dble(NPShock)*Deltax
    ShockWu=ShockPredictW(VP(NPShock),YF(1),VF(1),0.0d0)
    Wpre=ShockPredictW(VP(NPShock),YF(1),VF(1),ShockWu)
    W=ShockWu
    if (ShockWu>0.0d0) then
        Inter_type=1
    else
        Inter_type=2
    endif

    ShockHp2=HP(NPShock)
    ShockVp2=VP(NPShock)
    ShockYp1=YF(1)
    ShockVp1=VF(1)
    ShockCp1=CF(1)

    Om1Loc=HybridPressMass(ShockHp2)
    Om2Loc=Dens_ref*Area(ShockYp1)
    int_m(1)=Om1Loc
    int_m(2)=Om1Loc*ShockVp2
    int_m(3)=Om2Loc
    int_m(4)=Om2Loc*ShockVp1
End Subroutine

Subroutine ProjectHybridStateToFixedGrid(NPShock,NFShock,ShockLu,HP,VP,YF,VF,CF)
    Implicit None
    Integer,Intent(in)::NPShock,NFShock
    Real(8),Intent(in)::ShockLu,HP(NPShock),VP(NPShock),YF(NFShock),VF(NFShock),CF(NFShock)
    Integer::K,CellID,NFActive
    Real(8)::ShockDxP,ShockDxF,ShockDxFComp,XC

    ShockDxP=max(ShockLu/dble(NPShock),Tol1)
    call ComputeFreeHybridMetrics(NFShock,ShockLu,NFActive,ShockDxF,ShockDxFComp)

    do CellID=1,n,1
        XC=(dble(CellID)-0.5d0)*Deltax
        if (XC<=ShockLu) then
            K=min(NPShock,max(1,int(XC/ShockDxP)+1))
            H(CellID)=max(HP(K),Tol1)
            V(CellID)=VP(K)
            A(CellID)=Af
            Dens(CellID)=HybridPressMass(H(CellID))/Af
            C(CellID)=C_Waterhammer
        else
            K=min(NFActive,max(1,int((XC-ShockLu)/ShockDxF)+1))
            H(CellID)=min(max(YF(K),Tol1),y_ref-Tol1)
            V(CellID)=VF(K)
            A(CellID)=Area(H(CellID))
            Dens(CellID)=Dens_ref
            C(CellID)=CF(K)
        endif
        Q(CellID)=A(CellID)*V(CellID)
    enddo

    Include 'bc.inc'
End Subroutine

Real(8) Function HybridPressMass(HLoc)
    Implicit None
    Real(8),Intent(in)::HLoc

    HybridPressMass=Dens_ref*Af*(1.0d0+g*(HLoc-D)/C_Waterhammer**2)
    HybridPressMass=max(HybridPressMass,Tol1)
End Function

Real(8) Function HybridPressHead(OmLoc)
    Implicit None
    Real(8),Intent(in)::OmLoc

    HybridPressHead=D+((OmLoc/Af-Dens_ref)*C_Waterhammer**2)/(Dens_ref*g)
End Function

Real(8) Function HybridPressAPY(OmLoc)
    Implicit None
    Real(8),Intent(in)::OmLoc

    HybridPressAPY=Dens_ref*g*Af*(HybridPressHead(OmLoc)-0.5d0*D)
End Function

Subroutine UpdateHybridGasVolume(NFShock,NFActive,ShockDxF,YF,VolumeOut)
    Implicit None
    Integer,Intent(in)::NFShock,NFActive
    Real(8),Intent(in)::ShockDxF,YF(NFShock)
    Real(8),Intent(out)::VolumeOut
    Integer::K

    VolumeOut=0.0d0
    do K=1,NFActive,1
        VolumeOut=VolumeOut+ShockDxF*max(Af-Area(YF(K)),0.0d0)
    enddo
    VolumeOut=max(VolumeOut,Tol1)
End Subroutine

Subroutine UpdateFreeCelerity(NFShock,NFActive,YF,CF)
    Implicit None
    Integer,Intent(in)::NFShock,NFActive
    Real(8),Intent(in)::YF(NFShock)
    Real(8),Intent(out)::CF(NFShock)
    Integer::K

    do K=1,NFActive,1
        CF(K)=max(Celerity(YF(K)),Tol1)
    enddo
End Subroutine

Real(8) Function ShockPredictW(VL,HR,VR,Fallback)
    Implicit None
    Real(8),Intent(in)::VL,HR,VR,Fallback
    Real(8)::ARight

    ARight=Area(min(max(HR,Tol1),y_ref-Tol1))
    if (abs(Af-ARight)<=Tol1*max(1.0d0,Af)) then
        ShockPredictW=Fallback
    else
        ShockPredictW=(Af*VL-ARight*VR)/(Af-ARight)
    endif
    if (ShockPredictW/=ShockPredictW) ShockPredictW=Fallback
End Function

Subroutine AdvancePressurizedHybrid(NPShock,HP,VP,ShockDxP,DtLoc,HpGhost,VpGhost, &
    HNext,VNext)
    Implicit None
    Integer,Intent(in)::NPShock
    Real(8),Intent(in)::HP(NPShock),VP(NPShock),ShockDxP,DtLoc,HpGhost,VpGhost
    Real(8),Intent(out)::HNext(NPShock),VNext(NPShock)
    Integer::K
    Real(8)::Om(0:NPShock+1),Qm(0:NPShock+1),Flux(2,1:NPShock+1)
    Real(8)::HLeft,OmNew,QmNew,QmHalf,Um1Loc,Um2Loc,HLoc
    Real(8)::OmFloor

    OmFloor=0.5d0*Dens_ref*Af
    do K=1,NPShock,1
        Om(K)=HybridPressMass(HP(K))
        Qm(K)=Om(K)*VP(K)
    enddo

    select case (Left_boundary)
    case (1)
        HLeft=HP(1)
        Om(0)=HybridPressMass(HLeft)
        Qm(0)=-Qm(1)
    case (3)
        HLeft=LeftPressHead
        Om(0)=HybridPressMass(HLeft)
        Qm(0)=Qm(1)
    case default
        HLeft=HP(1)
        Om(0)=HybridPressMass(HLeft)
        Qm(0)=Qm(1)
    end select

    Om(NPShock+1)=HybridPressMass(max(HpGhost,-16.0d0*D+Tol1))
    Qm(NPShock+1)=Om(NPShock+1)*VpGhost

    do K=1,NPShock+1,1
        Um1Loc=exp(0.5d0*(log(max(Om(K-1),Tol1)*max(Om(K),Tol1))+ &
            (Qm(K-1)/max(Om(K-1),Tol1)-Qm(K)/max(Om(K),Tol1))/C_Waterhammer))
        Um2Loc=Um1Loc*0.5d0*(Qm(K-1)/max(Om(K-1),Tol1)+Qm(K)/max(Om(K),Tol1)+ &
            C_Waterhammer*(log(max(Om(K-1),Tol1))-log(max(Om(K),Tol1))))
        Flux(1,K)=Um2Loc
        Flux(2,K)=Um2Loc**2/max(Um1Loc,Tol1)+HybridPressAPY(max(Um1Loc,OmFloor))
    enddo

    do K=1,NPShock,1
        OmNew=Om(K)-DtLoc/ShockDxP*(Flux(1,K+1)-Flux(1,K))
        QmNew=Qm(K)-DtLoc/ShockDxP*(Flux(2,K+1)-Flux(2,K))
        OmNew=max(OmNew,OmFloor)
        QmHalf=QmNew+0.5d0*DtLoc*(g*OmNew*(S0-Friction*QmNew*abs(QmNew)/ &
            (max(OmNew,Tol1)**2*2.0d0*D*g)))
        QmNew=QmNew+DtLoc*(g*OmNew*(S0-Friction*QmHalf*abs(QmHalf)/ &
            (max(OmNew,Tol1)**2*2.0d0*D*g)))
        HLoc=max(HybridPressHead(OmNew),Tol1)
        HNext(K)=HLoc
        VNext(K)=QmNew/max(OmNew,Tol1)
        if (HNext(K)/=HNext(K) .or. VNext(K)/=VNext(K)) then
            HNext(K)=HP(K)
            VNext(K)=VP(K)
        endif
    enddo
End Subroutine

Subroutine AdvanceFreeHybrid(NFShock,NFActive,YF,VF,CF,ShockDxF,DtLoc,YpGhost,VpGhost, &
    CpGhost,YNext,VNext)
    Implicit None
    Integer,Intent(in)::NFShock,NFActive
    Real(8),Intent(in)::YF(NFShock),VF(NFShock),CF(NFShock),ShockDxF,DtLoc, &
        YpGhost,VpGhost,CpGhost
    Real(8),Intent(out)::YNext(NFShock),VNext(NFShock)
    Integer::K
    Real(8)::U1(0:NFShock+1),U2(0:NFShock+1),CLoc(0:NFShock+1),Flux(2,1:NFShock+1)
    Real(8)::YRight,U1New,U2New,AHere,RhydLoc,QmHalf,YLoc
    Real(8)::UM1Loc,UM2Loc,VLloc,VRloc,MLloc,MRloc,SLloc,SRloc,RanSlopeLoc
    Real(8)::ULVec(2),URVec(2),FLVec(2),FRVec(2),FLLVec(2)

    U1(0)=Dens_ref*Area(min(max(YpGhost,Tol1),y_ref-Tol1))
    U2(0)=U1(0)*VpGhost
    CLoc(0)=max(CpGhost,Tol1)
    do K=1,NFActive,1
        U1(K)=Dens_ref*Area(YF(K))
        U2(K)=U1(K)*VF(K)
        CLoc(K)=max(CF(K),Tol1)
    enddo

    select case (Right_boundary)
    case (1)
        YRight=YF(NFActive)
        U1(NFActive+1)=Dens_ref*Area(YRight)
        U2(NFActive+1)=-U2(NFActive)
    case (3)
        YRight=y_ref-Tol1
        U1(NFActive+1)=Dens_ref*Area(YRight)
        U2(NFActive+1)=U2(NFActive)
    case default
        YRight=YF(NFActive)
        U1(NFActive+1)=Dens_ref*Area(YRight)
        U2(NFActive+1)=U2(NFActive)
    end select
    CLoc(NFActive+1)=max(Celerity(YRight),Tol1)

    do K=1,NFActive+1,1
        VLloc=U2(K-1)/max(U1(K-1),Tol1)
        VRloc=U2(K)/max(U1(K),Tol1)
        UM1Loc=0.5d0*(U1(K-1)+U1(K))+(U1(K-1)+U1(K))*(VLloc-VRloc)/ &
            (2.0d0*max(CLoc(K-1)+CLoc(K),Tol1))
        UM1Loc=max(UM1Loc,Tol1)
        UM2Loc=0.5d0*UM1Loc*(VLloc+VRloc)+0.5d0*UM1Loc*(F_Equation(UM1Loc,U1(K))- &
            F_Equation(UM1Loc,U1(K-1)))
        ULVec(1)=U1(K-1)
        ULVec(2)=U2(K-1)
        URVec(1)=U1(K)
        URVec(2)=U2(K)
        FLVec(1)=U2(K-1)
        FLVec(2)=U2(K-1)**2/max(U1(K-1),Tol1)+APY(max(U1(K-1),Tol1))
        FRVec(1)=U2(K)
        FRVec(2)=U2(K)**2/max(U1(K),Tol1)+APY(max(U1(K),Tol1))
        if (UM1Loc>ULVec(1)) then
            MLloc=sqrt(max((APY(UM1Loc)-APY(max(ULVec(1),Tol1)))*UM1Loc/ &
                (max(ULVec(1),Tol1)*max(UM1Loc-ULVec(1),Tol1)),Tol1))
        else
            MLloc=CLoc(K-1)
        endif
        if (UM1Loc>URVec(1)) then
            MRloc=sqrt(max((APY(UM1Loc)-APY(max(URVec(1),Tol1)))*UM1Loc/ &
                (max(URVec(1),Tol1)*max(UM1Loc-URVec(1),Tol1)),Tol1))
        else
            MRloc=CLoc(K)
        endif
        SLloc=VLloc-MLloc
        SRloc=VRloc+MRloc
        FLLVec=(SRloc*FLVec-SLloc*FRVec+SRloc*SLloc*(URVec-ULVec))/ &
            max(SRloc-SLloc,Tol1)
        RanSlopeLoc=0.0d0
        Flux(:,K)=Choice(RanSlopeLoc,SLloc,SRloc,FLVec,FRVec,FLLVec)
    enddo

    do K=1,NFActive,1
        U1New=U1(K)-DtLoc/ShockDxF*(Flux(1,K+1)-Flux(1,K))
        U2New=U2(K)-DtLoc/ShockDxF*(Flux(2,K+1)-Flux(2,K))
        U1New=max(U1New,Tol1)
        AHere=U1New/Dens_ref
        RhydLoc=max(Hydraulic(AHere),Tol1)
        QmHalf=U2New+0.5d0*DtLoc*(g*Dens_ref*AHere*(S0-Manning**2*U2New*abs(U2New)/ &
            (max(U1New,Tol1)**2*RhydLoc**(1.33333d0))))
        U2New=U2New+DtLoc*(g*Dens_ref*AHere*(S0-Manning**2*QmHalf*abs(QmHalf)/ &
            (max(U1New,Tol1)**2*RhydLoc**(1.33333d0))))
        YLoc=min(max(Depth(U1New),Tol1),y_ref-Tol1)
        YNext(K)=YLoc
        VNext(K)=U2New/max(Dens_ref*Area(YLoc),Tol1)
        if (YNext(K)/=YNext(K) .or. VNext(K)/=VNext(K)) then
            YNext(K)=YF(K)
            VNext(K)=VF(K)
        endif
    enddo
    do K=NFActive+1,NFShock,1
        YNext(K)=YF(K)
        VNext(K)=VF(K)
    enddo
End Subroutine

Subroutine ShockTran(IfaceIndex,HL,VL,HR,VR,CR,ShockWOld,HGasLoc,DtLoc,ShockWNew, &
    ShockTag,ShockResn,ShockOk)
    Implicit None
    Integer,Intent(in)::IfaceIndex
    Real(8),Intent(in)::HL,VL,HR,VR,CR,ShockWOld,HGasLoc,DtLoc
    Real(8),Intent(out)::ShockWNew,ShockResn
    Character(16),Intent(out)::ShockTag
    Logical,Intent(out)::ShockOk
    Real(8)::OmLeft,OmRight,ForceDrive

    DiagStep=x
    DiagTime=Time
    DiagIface=IfaceIndex
    Wpre=ShockPredictW(VL,HR,VR,ShockWOld)
    OmLeft=HybridPressMass(max(HL,Tol1))
    OmRight=Dens_ref*Area(min(max(HR,Tol1),y_ref-Tol1))
    ForceDrive=Dens_ref*g*Af*(HL-0.5d0*D-HGasLoc)-APY(OmRight)
    if (EnableInterfaceDebugOutput) then
        Write(99,*) 'PRE',x,Time,IfaceIndex,Inter_type,Wpre,ShockWOld,OmLeft, &
            OmLeft*VL,OmRight,OmRight*VR,HL,HR,VL,VR,CR,VR-CR,int_m(1),int_m(2), &
            int_m(3),int_m(4),0.0d0,0.0d0
    endif

    if (ShockWOld>Tol1) then
        call ShockSolveT4(HL,VL,HR,VR,HGasLoc,DtLoc,Wpre,ShockWNew,ShockResn,ShockOk)
        ShockTag='T4'
    elseif (ShockWOld<-Tol1) then
        call ShockSolveT1(HL,VL,HR,VR,CR,HGasLoc,DtLoc,Wpre,ShockWNew,ShockResn,ShockOk)
        ShockTag='T1'
        if ((.not.ShockOk) .and. Ventilation/=0 .and. HR<=max(20.0d0*Tol1,5.0d-2*D)) then
            call ShockSolveT4(HL,VL,HR,VR,HGasLoc,DtLoc,Wpre,ShockWNew,ShockResn,ShockOk)
            if (ShockOk) ShockTag='T4FB'
        endif
    elseif (ForceDrive>=0.0d0) then
        call ShockSolveT4(HL,VL,HR,VR,HGasLoc,DtLoc,Wpre,ShockWNew,ShockResn,ShockOk)
        ShockTag='T4'
    else
        call ShockSolveT1(HL,VL,HR,VR,CR,HGasLoc,DtLoc,Wpre,ShockWNew,ShockResn,ShockOk)
        ShockTag='T1'
        if ((.not.ShockOk) .and. Ventilation/=0 .and. HR<=max(20.0d0*Tol1,5.0d-2*D)) then
            call ShockSolveT4(HL,VL,HR,VR,HGasLoc,DtLoc,Wpre,ShockWNew,ShockResn,ShockOk)
            if (ShockOk) ShockTag='T4FB'
        endif
    endif

    if (ShockOk) then
        W=ShockWNew
        if (EnableInterfaceDebugOutput) then
            Write(99,*) 'CHK3',x,Time,IfaceIndex,W,Wpre,int_m(1),int_m(2),int_m(3), &
                int_m(4),ShockVp2,ShockVp1,ShockHp2,ShockYp1,0.0d0,0.0d0,0.0d0, &
                ShockResn,W-Wpre
        endif
    endif
    Write(99,*) 'HYB_TRAN','EXIT',ShockTag,ShockWNew,ShockResn,ShockOk
End Subroutine

Subroutine ShockSolveT1(HL,VL,HR,VR,CR,HGasLoc,DtLoc,WStart,ShockWNew,ShockResn,ShockOk)
    Implicit None
    Real(8),Intent(in)::HL,VL,HR,VR,CR,HGasLoc,DtLoc,WStart
    Real(8),Intent(out)::ShockWNew,ShockResn
    Logical,Intent(out)::ShockOk
    Real(8)::X(4),R(4),J(4,4),RHS(4),DX(4),XT(4),RT(4),HStep,LambdaLoc
    Real(8)::Om1Loc,Om2Loc,ForceDrive
    Real(8)::TrialNorm
    Integer::IterLoc,K,LSLoc
    Logical::StateOK,TrialOK,AcceptLoc,SolveOK

    ShockWNew=ShockPredictW(VL,HR,VR,WStart)
    ForceDrive=Dens_ref*g*Af*(HL-0.5d0*D-HGasLoc)-APY(Dens_ref*Area(min(max(HR,Tol1),y_ref-Tol1)))
    if (ShockWNew>=-Tol1 .and. ForceDrive<0.0d0) ShockWNew=-10.0d0*Tol1
    X(1)=VL
    X(2)=max(HL,y_ref+Tol1)
    X(3)=VR
    X(4)=min(max(HR,10.0d0*Tol1),y_ref-Tol1)
    ShockResn=1.0d30
    ShockOk=.false.

    call ShockResidualT1(X,HL,VL,HR,VR,CR,HGasLoc,ShockWNew,R,StateOK)
    if (.not.StateOK) return
    ShockResn=ShockNormT1(X,R,HL,VL,HR,VR,CR,HGasLoc,ShockWNew)

    do IterLoc=1,40
        if (ShockResn<1.0d-6) then
            ShockOk=.true.
            exit
        endif
        do K=1,4,1
            XT=X
            HStep=1.0d-7*max(1.0d0,abs(X(K)))
            XT(K)=XT(K)+HStep
            call ShockResidualT1(XT,HL,VL,HR,VR,CR,HGasLoc,ShockWNew,RT,TrialOK)
            if (.not.TrialOK) RT=R+1.0d6
            J(:,K)=(RT-R)/HStep
        enddo
        RHS=-R
        call SolveSmallLinear4x4(J,RHS,DX,SolveOK)
        if (.not.SolveOK) exit
        AcceptLoc=.false.
        LambdaLoc=1.0d0
        do LSLoc=1,12,1
            XT=X+LambdaLoc*DX
            call ShockResidualT1(XT,HL,VL,HR,VR,CR,HGasLoc,ShockWNew,RT,TrialOK)
            if (TrialOK) then
                TrialNorm=ShockNormT1(XT,RT,HL,VL,HR,VR,CR,HGasLoc,ShockWNew)
                if (TrialNorm<ShockResn) then
                    X=XT
                    R=RT
                    ShockResn=TrialNorm
                    AcceptLoc=.true.
                    exit
                endif
            endif
            LambdaLoc=0.5d0*LambdaLoc
        enddo
        if (.not.AcceptLoc) exit
    enddo

    if (.not.ShockOk) then
        call ShockResidualT1(X,HL,VL,HR,VR,CR,HGasLoc,ShockWNew,R,StateOK)
        if (StateOK) then
            ShockResn=ShockNormT1(X,R,HL,VL,HR,VR,CR,HGasLoc,ShockWNew)
            if (ShockResn<1.0d-6) ShockOk=.true.
        endif
    endif

    if (ShockOk) then
        ShockVp2=X(1)
        ShockHp2=max(X(2),Tol1)
        ShockVp1=X(3)
        ShockYp1=min(max(X(4),Tol1),y_ref-Tol1)
        ShockCp1=max(Celerity(ShockYp1),Tol1)
        Om1Loc=HybridPressMass(ShockHp2)
        Om2Loc=Dens_ref*Area(ShockYp1)
        int_m(1)=Om1Loc
        int_m(2)=Om1Loc*ShockVp2
        int_m(3)=Om2Loc
        int_m(4)=Om2Loc*ShockVp1
        Inter_type=3
    endif
End Subroutine

Subroutine ShockSolveT4(HL,VL,HR,VR,HGasLoc,DtLoc,WStart,ShockWNew,ShockResn,ShockOk)
    Implicit None
    Real(8),Intent(in)::HL,VL,HR,VR,HGasLoc,DtLoc,WStart
    Real(8),Intent(out)::ShockWNew,ShockResn
    Logical,Intent(out)::ShockOk
    Real(8)::X(3),R(3),J(3,3),RHS(3),DX(3),XT(3),RT(3),HStep,LambdaLoc
    Real(8)::Om1Loc,Om2Loc,ForceDrive,ARight,SF2,HpSeed
    Real(8)::TrialNorm
    Integer::IterLoc,K,LSLoc
    Logical::StateOK,TrialOK,AcceptLoc,SolveOK,ReducedOK

    ARight=Area(min(max(HR,Tol1),y_ref-Tol1))
    Om2Loc=Dens_ref*ARight
    SF2=ShockT1Friction*VL*abs(VL)/(2.0d0*D)
    HpSeed=0.5d0*D+HGasLoc+APY(Om2Loc)/(Dens_ref*g*Af)
    X(2)=max(min(HpSeed,HL),-16.0d0*D+Tol1)
    X(1)=VL-g*(X(2)-HL)/C_Waterhammer-g*(ShockT1Slope-SF2)*ShockT1Deltat
    X(3)=(Af*X(1)-ARight*VR)/max(Af-ARight,Tol1)
    HpSeed=0.5d0*D+HGasLoc+APY(Om2Loc)/(Dens_ref*g*Af)-((X(1)-VR)*(X(3)+VR))/g
    X(2)=max(min(HpSeed,HL),-16.0d0*D+Tol1)
    X(1)=VL-g*(X(2)-HL)/C_Waterhammer-g*(ShockT1Slope-SF2)*ShockT1Deltat
    X(3)=(Af*X(1)-ARight*VR)/max(Af-ARight,Tol1)
    ForceDrive=Dens_ref*g*Af*(HL-0.5d0*D-HGasLoc)-APY(Om2Loc)
    if (X(3)<=Tol1 .and. ForceDrive>=0.0d0) X(3)=10.0d0*Tol1
    ShockResn=1.0d30
    ShockOk=.false.

    call ShockSolveT4Reduced(HL,VL,HR,VR,HGasLoc,X,ShockResn,ReducedOK)
    if (ReducedOK) ShockOk=.true.

    if (.not.ShockOk) then
        do IterLoc=1,80,1
            call ShockResidualT4(X,HL,VL,HR,VR,HGasLoc,R,StateOK)
            if (.not.StateOK) then
                if (EnableInterfaceDebugOutput .and. IterLoc==1) then
                    Write(99,*) 'T4SEEDBAD',x,Time,DiagIface,X(1),X(2),X(3)
                endif
                exit
            endif
            if (EnableInterfaceDebugOutput .and. IterLoc==1) then
                Write(99,*) 'T4SEED',x,Time,DiagIface,X(1),X(2),X(3),R(1),R(2),R(3)
            endif
            ShockResn=ShockNormT4(X,R,HL,VL,HR,VR,HGasLoc)
            if (ShockResn<1.0d-6) then
                ShockOk=.true.
                exit
            endif
            do K=1,3,1
                XT=X
                HStep=1.0d-7*max(1.0d0,abs(X(K)))
                XT(K)=XT(K)+HStep
                call ShockResidualT4(XT,HL,VL,HR,VR,HGasLoc,RT,TrialOK)
                if (.not.TrialOK) RT=R+1.0d6
                J(:,K)=(RT-R)/HStep
            enddo
            RHS=-R
            call SolveSmallLinear3x3(J,RHS,DX,SolveOK)
            if (.not.SolveOK) exit
            AcceptLoc=.false.
            LambdaLoc=1.0d0
            do LSLoc=1,16,1
                XT=X+LambdaLoc*DX
                call ShockResidualT4(XT,HL,VL,HR,VR,HGasLoc,RT,TrialOK)
                if (TrialOK) then
                    TrialNorm=ShockNormT4(XT,RT,HL,VL,HR,VR,HGasLoc)
                    if (TrialNorm<ShockResn) then
                        X=XT
                        R=RT
                        ShockResn=TrialNorm
                        AcceptLoc=.true.
                        exit
                    endif
                endif
                LambdaLoc=0.5d0*LambdaLoc
            enddo
            if (.not.AcceptLoc) exit
        enddo
    endif

    if (.not.ShockOk) then
        call ShockResidualT4(X,HL,VL,HR,VR,HGasLoc,R,StateOK)
        if (StateOK) then
            ShockResn=ShockNormT4(X,R,HL,VL,HR,VR,HGasLoc)
            if (ShockResn<1.0d-6) ShockOk=.true.
            if ((.not.ShockOk) .and. Ventilation/=0 .and. ShockResn<1.05d-4) ShockOk=.true.
            if ((.not.ShockOk) .and. Ventilation/=0 .and. DiagIface>=n .and. ShockResn<1.5d-4) ShockOk=.true.
        endif
    endif

    if (ShockOk) then
        ShockVp2=X(1)
        ShockHp2=max(X(2),-16.0d0*D+Tol1)
        ShockVp1=VR
        ShockYp1=min(max(HR,Tol1),y_ref-Tol1)
        ShockCp1=max(Celerity(ShockYp1),Tol1)
        ShockWNew=X(3)
        Om1Loc=HybridPressMass(ShockHp2)
        Om2Loc=Dens_ref*Area(ShockYp1)
        int_m(1)=Om1Loc
        int_m(2)=Om1Loc*ShockVp2
        int_m(3)=Om2Loc
        int_m(4)=Om2Loc*ShockVp1
        Inter_type=1
    endif
End Subroutine

Subroutine ShockSolveT4Reduced(HL,VL,HR,VR,HGasLoc,XOut,ShockResn,ShockOk)
    Implicit None
    Real(8),Intent(in)::HL,VL,HR,VR,HGasLoc
    Real(8),Intent(out)::XOut(3),ShockResn
    Logical,Intent(out)::ShockOk
    Integer,Parameter::NSample=128,NExpand=4
    Real(8)::HMin,HMax,HLoc,HPrev,FLoc,FPrev,HpBest,VBest,WBest,BestAbs,HScale
    Real(8)::HpBestScan,VBestScan,WBestScan,BestAbsScan
    Real(8)::VLoc,WLoc,HLeft,HRight,FLeft,FRight,HMid,FMid,VMid,WMid
    Real(8)::RLoc(3)
    Integer::I,IterLoc,ExpandLoc
    Logical::StateOK,PrevOK,BracketFound,BestSet,BestSetScan

    ShockOk=.false.
    ShockResn=1.0d30
    XOut=0.0d0
    HScale=2.0d0
    HMin=-HScale*D+Tol1
    HMax=min(HL,D)
    if (HMax<=HMin+Tol2) return

    BestAbs=1.0d300
    BestSet=.false.

    do ExpandLoc=1,NExpand,1
        HMin=-HScale*D+Tol1
        PrevOK=.false.
        BracketFound=.false.
        BestAbsScan=1.0d300
        BestSetScan=.false.

        do I=0,NSample,1
            HLoc=HMin+(HMax-HMin)*dble(I)/dble(NSample)
            call ShockT4ReducedResidual(HLoc,HL,VL,HR,VR,HGasLoc,VLoc,WLoc,FLoc,StateOK)
            if (.not.StateOK) cycle
            if (abs(FLoc)<BestAbsScan) then
                BestAbsScan=abs(FLoc)
                HpBestScan=HLoc
                VBestScan=VLoc
                WBestScan=WLoc
                BestSetScan=.true.
            endif
            if (abs(FLoc)<BestAbs) then
                BestAbs=abs(FLoc)
                HpBest=HLoc
                VBest=VLoc
                WBest=WLoc
                BestSet=.true.
            endif
            if (PrevOK) then
                if (FPrev==0.0d0 .or. FLoc==0.0d0 .or. FPrev*FLoc<0.0d0) then
                    HLeft=HPrev
                    HRight=HLoc
                    FLeft=FPrev
                    FRight=FLoc
                    BracketFound=.true.
                    exit
                endif
            endif
            HPrev=HLoc
            FPrev=FLoc
            PrevOK=.true.
        enddo

        if (BracketFound) then
            do IterLoc=1,60,1
                HMid=0.5d0*(HLeft+HRight)
                call ShockT4ReducedResidual(HMid,HL,VL,HR,VR,HGasLoc,VMid,WMid,FMid,StateOK)
                if (.not.StateOK) exit
                if (abs(FMid)<Tol2 .or. abs(HRight-HLeft)<1.0d-10) then
                    HpBest=HMid
                    VBest=VMid
                    WBest=WMid
                    BestSet=.true.
                    exit
                endif
                if (FLeft*FMid<=0.0d0) then
                    HRight=HMid
                    FRight=FMid
                else
                    HLeft=HMid
                    FLeft=FMid
                endif
                HpBest=HMid
                VBest=VMid
                WBest=WMid
                BestSet=.true.
            enddo
            exit
        endif

        if (.not.BestSetScan) exit
        if (abs(HpBestScan-HMin)>1.0d-10) exit
        HScale=2.0d0*HScale
    enddo

    if (.not.BestSet) return

    call ShockResidualT4((/VBest,HpBest,WBest/),HL,VL,HR,VR,HGasLoc,RLoc,StateOK)
    if (.not.StateOK) return
    XOut(1)=VBest
    XOut(2)=HpBest
    XOut(3)=WBest
    ShockResn=ShockNormT4(XOut,RLoc,HL,VL,HR,VR,HGasLoc)
    if (ShockResn<1.0d-6) ShockOk=.true.
    if ((.not.ShockOk) .and. Ventilation/=0 .and. ShockResn<1.05d-4) ShockOk=.true.
    if ((.not.ShockOk) .and. Ventilation/=0 .and. DiagIface>=n .and. ShockResn<1.5d-4) ShockOk=.true.
End Subroutine

Subroutine ShockT4ReducedResidual(HpLoc,HL,VL,HR,VR,HGasLoc,VpLoc,WLoc,FLoc,StateOK)
    Implicit None
    Real(8),Intent(in)::HpLoc,HL,VL,HR,VR,HGasLoc
    Real(8),Intent(out)::VpLoc,WLoc,FLoc
    Logical,Intent(out)::StateOK
    Real(8)::ARight,Om2Loc,SF2,XLoc(3),RLoc(3)

    StateOK=.false.
    VpLoc=0.0d0
    WLoc=0.0d0
    FLoc=0.0d0
    if (HpLoc<-16.0d0*D+Tol1) return
    ARight=Area(min(max(HR,Tol1),y_ref-Tol1))
    if (abs(Af-ARight)<=Tol2) return
    Om2Loc=Dens_ref*ARight
    SF2=ShockT1Friction*VL*abs(VL)/(2.0d0*D)
    VpLoc=VL-g*(HpLoc-HL)/C_Waterhammer-g*(ShockT1Slope-SF2)*ShockT1Deltat
    WLoc=(Af*VpLoc-ARight*VR)/max(Af-ARight,Tol1)
    XLoc(1)=VpLoc
    XLoc(2)=HpLoc
    XLoc(3)=WLoc
    call ShockResidualT4(XLoc,HL,VL,HR,VR,HGasLoc,RLoc,StateOK)
    if (.not.StateOK) return
    FLoc=RLoc(3)
End Subroutine

Subroutine ShockResidualT1(X,HL,VL,HR,VR,CR,HGasLoc,WLoc,RLoc,StateOK)
    Implicit None
    Real(8),Intent(in)::X(4),HL,VL,HR,VR,CR,HGasLoc,WLoc
    Real(8),Intent(out)::RLoc(4)
    Logical,Intent(out)::StateOK
    Real(8)::ARight,CRight,RhydRight,SF1,SF2,Om1Loc,Om2Loc

    StateOK=.false.
    RLoc=0.0d0
    if (X(2)<=Tol1) return
    if (X(4)<=Tol1 .or. X(4)>=y_ref-Tol1) return
    if (X(1)/=X(1) .or. X(2)/=X(2) .or. X(3)/=X(3) .or. X(4)/=X(4)) return

    ARight=Area(X(4))
    CRight=max(Celerity(X(4)),Tol1)
    RhydRight=max(Hydraulic(ARight),Tol1)
    Om1Loc=HybridPressMass(X(2))
    Om2Loc=Dens_ref*ARight
    SF1=ShockT1Manning**2*VR*abs(VR)/(RhydRight**(1.33333d0))
    SF2=ShockT1Friction*VL*abs(VL)/(2.0d0*D)
    RLoc(1)=X(1)-VL+g*(X(2)-HL)/C_Waterhammer+g*(ShockT1Slope-SF2)*ShockT1Deltat
    RLoc(2)=X(3)-VR-g*X(4)/CRight+g*HR/max(CR,Tol1)+g*(ShockT1Slope-SF1)*ShockT1Deltat
    RLoc(3)=Af*X(1)-ARight*X(3)-WLoc*(Af-ARight)
    RLoc(4)=Dens_ref*g*Af*(X(2)-0.5d0*D)-APY(Om2Loc)-Dens_ref*g*Af*HGasLoc- &
        Dens_ref*Af*X(1)*(X(3)-X(1))+Dens_ref*Af*WLoc*(X(3)-X(1))
    if (RLoc(1)/=RLoc(1) .or. RLoc(2)/=RLoc(2) .or. RLoc(3)/=RLoc(3) .or. RLoc(4)/=RLoc(4)) return
    StateOK=.true.
End Subroutine

Subroutine ShockResidualT4(X,HL,VL,HR,VR,HGasLoc,RLoc,StateOK)
    Implicit None
    Real(8),Intent(in)::X(3),HL,VL,HR,VR,HGasLoc
    Real(8),Intent(out)::RLoc(3)
    Logical,Intent(out)::StateOK
    Real(8)::ARight,Om1Loc,Om2Loc,SF2

    StateOK=.false.
    RLoc=0.0d0
    if (X(2)<-16.0d0*D+Tol1) return
    if (X(1)/=X(1) .or. X(2)/=X(2) .or. X(3)/=X(3)) return

    ARight=Area(min(max(HR,Tol1),y_ref-Tol1))
    Om1Loc=HybridPressMass(X(2))
    Om2Loc=Dens_ref*ARight
    SF2=ShockT1Friction*VL*abs(VL)/(2.0d0*D)
    RLoc(1)=X(1)-VL+g*(X(2)-HL)/C_Waterhammer+g*(ShockT1Slope-SF2)*ShockT1Deltat
    RLoc(2)=Af*X(1)-ARight*VR-X(3)*(Af-ARight)
    RLoc(3)=Dens_ref*g*Af*(X(2)-0.5d0*D)-APY(Om2Loc)-Dens_ref*g*Af*HGasLoc- &
        Dens_ref*Af*X(3)*(VR-X(1))-Dens_ref*Af*VR*(VR-X(1))
    if (RLoc(1)/=RLoc(1) .or. RLoc(2)/=RLoc(2) .or. RLoc(3)/=RLoc(3)) return
    StateOK=.true.
End Subroutine

Real(8) Function ShockNormT1(X,RLoc,HL,VL,HR,VR,CR,HGasLoc,WLoc)
    Implicit None
    Real(8),Intent(in)::X(4),RLoc(4),HL,VL,HR,VR,CR,HGasLoc,WLoc
    Real(8)::ARight,CRight,Om2Loc,Scale1,Scale2,Scale3,Scale4

    ARight=Area(X(4))
    CRight=max(Celerity(X(4)),Tol1)
    Om2Loc=Dens_ref*ARight
    Scale1=max(1.0d0,abs(VL),abs(X(1)),abs(g*(X(2)-HL)/C_Waterhammer))
    Scale2=max(1.0d0,abs(VR),abs(X(3)),abs(g*X(4)/CRight),abs(g*HR/max(CR,Tol1)))
    Scale3=max(Tol1,abs(Af*X(1))+abs(ARight*X(3))+abs(WLoc*(Af-ARight)))
    Scale4=max(1.0d0,abs(Dens_ref*g*Af*(X(2)-0.5d0*D))+abs(APY(Om2Loc))+ &
        abs(Dens_ref*g*Af*HGasLoc)+abs(Dens_ref*Af*X(1)*(X(3)-X(1)))+ &
        abs(Dens_ref*Af*WLoc*(X(3)-X(1))))
    ShockNormT1=max(abs(RLoc(1))/Scale1,max(abs(RLoc(2))/Scale2,max(abs(RLoc(3))/Scale3, &
        abs(RLoc(4))/Scale4)))
End Function

Real(8) Function ShockNormT4(X,RLoc,HL,VL,HR,VR,HGasLoc)
    Implicit None
    Real(8),Intent(in)::X(3),RLoc(3),HL,VL,HR,VR,HGasLoc
    Real(8)::ARight,Om2Loc,Scale1,Scale2,Scale3

    ARight=Area(min(max(HR,Tol1),y_ref-Tol1))
    Om2Loc=Dens_ref*ARight
    Scale1=max(1.0d0,abs(VL),abs(X(1)),abs(g*(X(2)-HL)/C_Waterhammer))
    Scale2=max(Tol1,abs(Af*X(1))+abs(ARight*VR)+abs(X(3)*(Af-ARight)))
    Scale3=max(1.0d0,abs(Dens_ref*g*Af*(X(2)-0.5d0*D))+abs(APY(Om2Loc))+ &
        abs(Dens_ref*g*Af*HGasLoc)+abs(Dens_ref*Af*X(3)*(VR-X(1)))+ &
        abs(Dens_ref*Af*VR*(VR-X(1))))
    ShockNormT4=max(abs(RLoc(1))/Scale1,max(abs(RLoc(2))/Scale2,abs(RLoc(3))/Scale3))
End Function

Subroutine SolveSmallLinear4x4(AMat,BVec,XVec,SolveOK)
    Implicit None
    Real(8),Intent(in)::AMat(4,4),BVec(4)
    Real(8),Intent(out)::XVec(4)
    Logical,Intent(out)::SolveOK
    Real(8)::M(4,4),RHS(4),RowTmp(4),Fac,ValTmp
    Integer::I,K,P

    M=AMat
    RHS=BVec
    SolveOK=.true.
    do K=1,3,1
        P=K-1+maxloc(abs(M(K:4,K)),1)
        if (abs(M(P,K))<1.0d-14) then
            SolveOK=.false.
            return
        endif
        if (P/=K) then
            RowTmp=M(K,:)
            M(K,:)=M(P,:)
            M(P,:)=RowTmp
            ValTmp=RHS(K)
            RHS(K)=RHS(P)
            RHS(P)=ValTmp
        endif
        do I=K+1,4,1
            Fac=M(I,K)/M(K,K)
            M(I,K:4)=M(I,K:4)-Fac*M(K,K:4)
            RHS(I)=RHS(I)-Fac*RHS(K)
        enddo
    enddo
    if (abs(M(4,4))<1.0d-14) then
        SolveOK=.false.
        return
    endif
    XVec(4)=RHS(4)/M(4,4)
    do I=3,1,-1
        XVec(I)=(RHS(I)-sum(M(I,I+1:4)*XVec(I+1:4)))/M(I,I)
    enddo
End Subroutine

Subroutine SolveSmallLinear3x3(AMat,BVec,XVec,SolveOK)
    Implicit None
    Real(8),Intent(in)::AMat(3,3),BVec(3)
    Real(8),Intent(out)::XVec(3)
    Logical,Intent(out)::SolveOK
    Real(8)::M(3,3),RHS(3),RowTmp(3),Fac,ValTmp
    Integer::I,K,P

    M=AMat
    RHS=BVec
    SolveOK=.true.
    do K=1,2,1
        P=K-1+maxloc(abs(M(K:3,K)),1)
        if (abs(M(P,K))<1.0d-14) then
            SolveOK=.false.
            return
        endif
        if (P/=K) then
            RowTmp=M(K,:)
            M(K,:)=M(P,:)
            M(P,:)=RowTmp
            ValTmp=RHS(K)
            RHS(K)=RHS(P)
            RHS(P)=ValTmp
        endif
        do I=K+1,3,1
            Fac=M(I,K)/M(K,K)
            M(I,K:3)=M(I,K:3)-Fac*M(K,K:3)
            RHS(I)=RHS(I)-Fac*RHS(K)
        enddo
    enddo
    if (abs(M(3,3))<1.0d-14) then
        SolveOK=.false.
        return
    endif
    XVec(3)=RHS(3)/M(3,3)
    do I=2,1,-1
        XVec(I)=(RHS(I)-sum(M(I,I+1:3)*XVec(I+1:3)))/M(I,I)
    enddo
End Subroutine

Subroutine HybridUpdateNegativeTrackState()
    Implicit None

    if (int_m(1)<=0.0d0 .or. int_m(3)<=0.0d0) return
    if (.not.allocated(NegTrackPendingReadyMap)) return
    if (DiagIface<lbound(NegTrackPendingReadyMap,1) .or. DiagIface>ubound(NegTrackPendingReadyMap,1)) return
    NegTrackPendingReadyMap(DiagIface)=.true.
    NegTrackPendingStepMap(DiagIface)=DiagStep
    NegTrackPendingOm1Map(DiagIface)=int_m(1)
    NegTrackPendingQ1Map(DiagIface)=int_m(2)
    NegTrackPendingOm2Map(DiagIface)=int_m(3)
    NegTrackPendingQ2Map(DiagIface)=int_m(4)
    NegTrackPendingWMap(DiagIface)=W
End Subroutine

Subroutine HybridPrepareNegativeTrackForStep()
    Implicit None
    Integer::ITrack

    if (.not.allocated(NegTrackReadyMap)) then
        call HybridClearNegativeTrackSelection()
        return
    endif

    if (NegTrackPreparedStep/=DiagStep) then
        NegTrackReadyMap=.false.
        NegTrackStepMap=0
        do ITrack=lbound(NegTrackPendingReadyMap,1),ubound(NegTrackPendingReadyMap,1)
            if (NegTrackPendingReadyMap(ITrack)) then
                if (NegTrackPendingStepMap(ITrack)<DiagStep) then
                    NegTrackReadyMap(ITrack)=.true.
                    NegTrackStepMap(ITrack)=NegTrackPendingStepMap(ITrack)
                    NegTrackOm1Map(ITrack)=NegTrackPendingOm1Map(ITrack)
                    NegTrackQ1Map(ITrack)=NegTrackPendingQ1Map(ITrack)
                    NegTrackOm2Map(ITrack)=NegTrackPendingOm2Map(ITrack)
                    NegTrackQ2Map(ITrack)=NegTrackPendingQ2Map(ITrack)
                    NegTrackWMap(ITrack)=NegTrackPendingWMap(ITrack)
                endif
                NegTrackPendingReadyMap(ITrack)=.false.
            endif
        enddo
        NegTrackPreparedStep=DiagStep
    endif

    call HybridLoadNegativeTrackForIface(DiagIface)
End Subroutine

Subroutine HybridClearNegativeTrackSelection()
    Implicit None

    NegTrackReady=.false.
    NegTrackStep=0
    NegTrackIface=0
    NegTrackOm1=0.0d0
    NegTrackQ1=0.0d0
    NegTrackOm2=0.0d0
    NegTrackQ2=0.0d0
    NegTrackW=0.0d0
End Subroutine

Subroutine HybridLoadNegativeTrackForIface(IfaceWanted)
    Implicit None
    Integer,Intent(in)::IfaceWanted
    Integer::ITrack,BestIface,BestDelta

    call HybridClearNegativeTrackSelection()
    if (.not.allocated(NegTrackReadyMap)) return
    if (IfaceWanted<lbound(NegTrackReadyMap,1) .or. IfaceWanted>ubound(NegTrackReadyMap,1)) return
    BestIface=0
    BestDelta=huge(1)

    do ITrack=max(lbound(NegTrackReadyMap,1),IfaceWanted-1),min(ubound(NegTrackReadyMap,1),IfaceWanted+1)
        if (NegTrackReadyMap(ITrack)) then
            if (abs(ITrack-IfaceWanted)<BestDelta) then
                BestIface=ITrack
                BestDelta=abs(ITrack-IfaceWanted)
            endif
        endif
    enddo
    if (BestIface==0) return

    NegTrackReady=.true.
    NegTrackStep=NegTrackStepMap(BestIface)
    NegTrackIface=BestIface
    NegTrackOm1=NegTrackOm1Map(BestIface)
    NegTrackQ1=NegTrackQ1Map(BestIface)
    NegTrackOm2=NegTrackOm2Map(BestIface)
    NegTrackQ2=NegTrackQ2Map(BestIface)
    NegTrackW=NegTrackWMap(BestIface)
End Subroutine
