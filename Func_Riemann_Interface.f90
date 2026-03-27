Function Riemann_Interface(OmigaL,OmigaR,Qm_L,Qm_R) Result(UM)
    use variables
    Implicit none
    Real(8)::OmigaL,OmigaR,Qm_L,Qm_R,HL,HR,VL,VR !HP是我们需要求解到的值，中间态的水深
    Real(8)::UM(2)
    Real(8),External::APY,Depth,Density,Area,SOLU,Celerity,F_Equation,Phi_Circular
    Real(8)::F1,F2,F3,HP1,HP2,HP3,OmigaP1,OmigaP2,OmigaP3
    Real(8)::k1,k2,k3,kk
    Real(8)::R1,R2,R1P,R2P,H1,H2,H1T,H2T,Omiga1,Omiga2,U1,U2,Omiga1T,Omiga2T,U1T,U2T
    Real(8)::DH1,DH2,Det,Delta1,Delta2,Lambda,ResNorm,TrialNorm,Scale1,Scale2,Scale3,HgasLocal
    Real(8)::J11,J12,J21,J22
    Real(8)::ResNormInit,ScanBestNorm,ScanBestH1,ScanBestH2,ScanBestR1,ScanBestR2
    Real(8)::ScanH1Min,ScanH1Max,ScanH2Min,ScanH2Max,RootBestAbs,RootBestH1,RootBestH2,RootBestEnergy
    Real(8)::WState,RightCharState,ShadowRightCharState
    Real(8)::SavedWpre
    Real(8)::Mode2WNearestRoot,Mode2WNearestDelta
    Real(8)::S1SeedH1,S1SeedH2,S1SeedNorm,S2SeedH1,S2SeedH2,S2SeedW,S2SeedNorm
    Real(8)::RuntimeCharState,RuntimeTrackU2,RuntimeWTrack
    Integer::Iter,LSIter,IterLast,LSIterLast,ScanValidCount,RootSampleCount,RootBracketCount
    Integer::RootScanNScanH1,RootScanNScanH2,RootScanMaxBisect
    Integer::Mode2WRootHitCount
    Logical::PositiveSolved,NegativeSolved,ValidState,AcceptStep,UseRightCharEq,RootHasSignChange,EnableStrategy2,UseShadowRightChar
    Logical::HasS1Seed,HasS2Seed,RuntimeTrackedReady,RuntimeWTrackedReady
    Character(32)::NegativeExitReason,NegativeRootStatus,RuntimeTrackStatus,RuntimeWTrackStatus

    VL=Qm_L/OmigaL
    VR=Qm_R/OmigaR
    HL=Depth(OmigaL)
    HR=Depth(OmigaR)
    Inter_type=1

    goto 100

    kk=APY(OmigaR)
    HP1=max(30.0d0,HL+Tol1)
    HP2=max(Tol1,min(0.95d0*y_ref+tolerance,y_ref-Tol1))

    OmigaP2=Density(HP2)*Area(HP2)
    F2=VL+C_Waterhammer*LOG(OmigaL/OmigaP2)-VR-sqrt((APY(OmigaP2)-APY(OmigaR))*(OmigaP2-OmigaR)/(OmigaP2*OmigaR))
    if (F2 /= F2) then
        HP2=max(Tol1,min(SOLU(APY(OmigaR))+Tol1,y_ref-Tol1))
    endif

    DO Iter=1,200 !二分法
        HP3=0.5d0*(HP1+HP2)

        call EvaluatePositiveState(HP1,OmigaP1,k1)
        call EvaluatePositiveState(HP2,OmigaP2,k2)
        call EvaluatePositiveState(HP3,OmigaP3,k3)

        F1=VL+C_Waterhammer*LOG(OmigaL/OmigaP1)-VR-sqrt((k1-kk)*(OmigaP1-OmigaR)/(OmigaP1*OmigaR))
        F2=VL+C_Waterhammer*LOG(OmigaL/OmigaP2)-VR-sqrt((k2-kk)*(OmigaP2-OmigaR)/(OmigaP2*OmigaR))
        F3=VL+C_Waterhammer*LOG(OmigaL/OmigaP3)-VR-sqrt((k3-kk)*(OmigaP3-OmigaR)/(OmigaP3*OmigaR))

        if (F1 /= F1 .or. F2 /= F2 .or. F3 /= F3) exit
        if (F1*F2>0.0d0) exit

        if (abs(F3)<1.0d-5 .or. abs(HP1-HP2)<Tol1) then
            PositiveSolved=.true.
            OmigaP1=OmigaP3
            U1=VL+C_Waterhammer*LOG(OmigaL/OmigaP1)
            UM(1)=OmigaP1
            UM(2)=OmigaP1*U1
            if (abs(OmigaP1-OmigaR)>Tol1*max(1.0d0,abs(OmigaR))) then
                w=(OmigaP1*U1-OmigaR*VR)/(OmigaP1-OmigaR)
            else
                w=Wpre
                if (EnableRuntimeStepPrint) Print*,"模式1兜底"
            endif
            int_m=(/UM(1),UM(2),OmigaR,Qm_R/)
            return
        endif

        If(F1*F3<0.0d0) Then
            HP2=HP3
        Else
            HP1=HP3
        Endif
    Enddo

100 continue
    HgasLocal=0.0d0
    HasS1Seed=.false.
    HasS2Seed=.false.
    S1SeedNorm=-1.0d0
    S2SeedNorm=-1.0d0
    Mode2WRootHitCount=0
    Mode2WNearestRoot=0.0d0
    Mode2WNearestDelta=-1.0d0
    Scale1=max(1.0d0,abs(APY(OmigaL)),abs(APY(OmigaR)))
    PositiveSolved=.false.
    NegativeSolved=.false.
    UseShadowRightChar=.false.
    ShadowRightCharState=0.0d0
    call PrepareNegativeTrackForStep()
    if (EnableNegativeShadowDiagnostics) then
        call WriteNegativeRootDiag(.true.,2,'M2ROOT')
        call WriteNegativeTrackShadowDiag()
        call WriteNegativeRightCharShadowDiag()
        call WriteNegativeRightCharTryDiag()
    endif
    call GetTrackedRightChar(RuntimeCharState,RuntimeTrackU2,RuntimeTrackedReady,RuntimeTrackStatus)
    if (NegUseTrackedRightCharRuntime .and. RuntimeTrackedReady) then
        UseShadowRightChar=.true.
        ShadowRightCharState=RuntimeCharState
        call TryNegativeMode(.true.,2,'M2FAILT','M2SCANT',NegativeSolved,UM)
        if (.not.NegativeSolved) then
            call GetTrackedW(RuntimeWTrack,RuntimeWTrackedReady,RuntimeWTrackStatus)
            if (RuntimeWTrackedReady) then
                SavedWpre=Wpre
                Wpre=RuntimeWTrack
                call TryNegativeMode(.true.,2,'M2FAILTW','M2SCANTW',NegativeSolved,UM)
                Wpre=SavedWpre
            endif
            if (.not.NegativeSolved) then
                if (.not.EnableNegativeShadowDiagnostics) call WriteNegativeRootDiag(.true.,2,'M2ROOT')
                UseShadowRightChar=.false.
                ShadowRightCharState=0.0d0
                call TryNegativeMode(.true.,2,'M2FAIL','M2SCAN',NegativeSolved,UM)
            endif
        endif
    else
        if (.not.EnableNegativeShadowDiagnostics) call WriteNegativeRootDiag(.true.,2,'M2ROOT')
        call TryNegativeMode(.true.,2,'M2FAIL','M2SCAN',NegativeSolved,UM)
    endif

    if (.not.NegativeSolved) then
        UM(1)=OmigaL
        UM(2)=Qm_L
        int_m=(/OmigaL,Qm_L,OmigaR,Qm_R/)
    else
        call UpdateNegativeTrackState()
    endif

    return

contains

    subroutine SetNegativeModeControl(UseRightCharEqWanted,InterTypeWanted)
        Implicit none
        Logical,Intent(in)::UseRightCharEqWanted
        Integer,Intent(in)::InterTypeWanted

        UseRightCharEq=UseRightCharEqWanted
        Inter_type=InterTypeWanted
        WState=Wpre
        if (UseShadowRightChar) then
            RightCharState=ShadowRightCharState
        else
            RightCharState=VR-Phi_Circular(OmigaR)
        endif
        RootScanNScanH1=61
        RootScanNScanH2=81
        RootScanMaxBisect=60
        if (UseRightCharEq) then
            Scale2=max(1.0d0,abs(VL),abs(VR),abs(Phi_Circular(OmigaL)),abs(Phi_Circular(OmigaR)), &
                abs(RightCharState))
        else
            Scale2=max(1.0d0,abs(OmigaL*VL*VL+2.0d0*OmigaL*g*HL), &
                abs(OmigaR*VR*VR+2.0d0*OmigaR*g*HR),abs(HgasLocal))
            Scale2=max(Scale2,abs(APY(OmigaL)),abs(APY(OmigaR)))
        endif
    end subroutine

    subroutine GetNegativeWWindow(SpanFactor,WMin,WMax)
        Implicit none
        Real(8),Intent(in)::SpanFactor
        Real(8),Intent(out)::WMin,WMax
        Real(8)::SpeedScale,WSpan
        Real(8),External::Celerity

        SpeedScale=max(abs(VL),abs(VR),Celerity(HR),abs(Wpre),1.0d-3)
        WSpan=max(5.0d-2,SpanFactor*SpeedScale)
        WMin=Wpre-WSpan
        WMax=Wpre+WSpan
    end subroutine

    subroutine WriteNegativeRootDiag(UseRightCharEqWanted,InterTypeWanted,RootTag)
        Implicit none
        Logical,Intent(in)::UseRightCharEqWanted
        Integer,Intent(in)::InterTypeWanted
        Character(*),Intent(in)::RootTag

        call SetNegativeModeControl(UseRightCharEqWanted,InterTypeWanted)
        Mode2WRootHitCount=0
        Mode2WNearestRoot=0.0d0
        Mode2WNearestDelta=-1.0d0
        call CheckNegativeExistence(RootBestAbs,RootBestH1,RootBestH2,RootBestEnergy, &
            RootSampleCount,RootBracketCount,RootHasSignChange,NegativeRootStatus)
        HasS1Seed=.false.
        S1SeedNorm=-1.0d0
        if (RootSampleCount>0 .and. RootBestAbs>=0.0d0) then
            HasS1Seed=.true.
            S1SeedH1=RootBestH1
            S1SeedH2=RootBestH2
            S1SeedNorm=RootBestAbs
        endif
        if (EnableNegativeShadowDiagnostics) then
            Write(99,*) trim(RootTag),trim(NegativeRootStatus),RootSampleCount,RootBracketCount, &
                RootHasSignChange,RootBestAbs,RootBestH1,RootBestH2,RootBestEnergy
        endif
        if (UseRightCharEqWanted .and. trim(NegativeRootStatus)/='has_root') then
            call WriteNegativeWScanDiag('M2WSCAN')
        endif
    end subroutine

    subroutine WriteNegativeTrackShadowDiag()
        Implicit none
        Real(8)::WTrack,SavedWpre,TrackBestAbs,TrackBestH1,TrackBestH2,TrackBestEnergy,DenTrack
        Integer::TrackSampleCount,TrackBracketCount
        Logical::TrackHasSignChange,TrackReadyNow
        Character(32)::TrackRootStatus

        if (.not.EnableNegativeShadowDiagnostics) return

        WTrack=Wpre
        TrackBestAbs=-1.0d0
        TrackBestH1=0.0d0
        TrackBestH2=0.0d0
        TrackBestEnergy=0.0d0
        TrackSampleCount=0
        TrackBracketCount=0
        TrackHasSignChange=.false.
        TrackReadyNow=.false.
        TrackRootStatus='not_ready'

        if (NegTrackReady) then
            if (abs(DiagIface-NegTrackIface)>1) then
                WTrack=NegTrackW
                TrackRootStatus='stale_iface'
            else
                DenTrack=NegTrackOm2-NegTrackOm1
                if (abs(DenTrack)>Tol1*max(1.0d0,abs(NegTrackOm1),abs(NegTrackOm2))) then
                    WTrack=(NegTrackQ2-NegTrackQ1)/DenTrack
                    if (WTrack==WTrack) then
                        TrackReadyNow=.true.
                        SavedWpre=Wpre
                        Wpre=WTrack
                        call SetNegativeModeControl(.true.,2)
                        call CheckNegativeExistence(TrackBestAbs,TrackBestH1,TrackBestH2,TrackBestEnergy, &
                            TrackSampleCount,TrackBracketCount,TrackHasSignChange,TrackRootStatus)
                        Wpre=SavedWpre
                        call SetNegativeModeControl(.true.,2)
                    else
                        WTrack=NegTrackW
                        TrackRootStatus='nan_track'
                    endif
                else
                    WTrack=NegTrackW
                    TrackRootStatus='degenerate'
                endif
            endif
        endif

        Write(99,*) 'M2WTRACK',DiagStep,DiagTime,DiagIface,Wpre,WTrack,TrackReadyNow, &
            trim(NegativeRootStatus),RootBestAbs,RootSampleCount,RootBracketCount, &
            trim(TrackRootStatus),TrackBestAbs,TrackSampleCount,TrackBracketCount,NegTrackIface,NegTrackW
    end subroutine

    subroutine UpdateNegativeTrackState()
        Implicit none

        if (int_m(1)<=0.0d0 .or. int_m(3)<=0.0d0) return
        if (.not.allocated(NegTrackPendingReadyMap)) return
        if (DiagIface<lbound(NegTrackPendingReadyMap,1) .or. DiagIface>ubound(NegTrackPendingReadyMap,1)) return
        NegTrackPendingReadyMap(DiagIface)=.true.
        NegTrackPendingStepMap(DiagIface)=DiagStep
        NegTrackPendingOm1Map(DiagIface)=int_m(1)
        NegTrackPendingQ1Map(DiagIface)=int_m(2)
        NegTrackPendingOm2Map(DiagIface)=int_m(3)
        NegTrackPendingQ2Map(DiagIface)=int_m(4)
        NegTrackPendingWMap(DiagIface)=w
    end subroutine

    subroutine PrepareNegativeTrackForStep()
        Implicit none
        Integer::ITrack

        if (.not.allocated(NegTrackReadyMap)) then
            call ClearNegativeTrackSelection()
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

        call LoadNegativeTrackForIface(DiagIface)
    end subroutine

    subroutine ClearNegativeTrackSelection()
        Implicit none

        NegTrackReady=.false.
        NegTrackStep=0
        NegTrackIface=0
        NegTrackOm1=0.0d0
        NegTrackQ1=0.0d0
        NegTrackOm2=0.0d0
        NegTrackQ2=0.0d0
        NegTrackW=0.0d0
    end subroutine

    subroutine LoadNegativeTrackForIface(IfaceWanted)
        Implicit none
        Integer,Intent(in)::IfaceWanted
        Integer::ITrack,BestIface,BestDelta

        call ClearNegativeTrackSelection()
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
    end subroutine

    subroutine GetTrackedRightChar(CharTrack,TrackU2,TrackReadyNow,TrackStatus)
        Implicit none
        Real(8),Intent(out)::CharTrack,TrackU2
        Logical,Intent(out)::TrackReadyNow
        Character(32),Intent(out)::TrackStatus

        CharTrack=VR-Phi_Circular(OmigaR)
        TrackU2=0.0d0
        TrackReadyNow=.false.
        TrackStatus='not_ready'

        if (NegTrackReady) then
            if (NegTrackOm2>0.0d0) then
                TrackU2=NegTrackQ2/NegTrackOm2
                CharTrack=TrackU2-Phi_Circular(NegTrackOm2)
            endif
            if (abs(DiagIface-NegTrackIface)>1) then
                TrackStatus='stale_iface'
            elseif (NegTrackOm2<=0.0d0) then
                TrackStatus='degenerate'
            elseif (CharTrack==CharTrack) then
                TrackReadyNow=.true.
                TrackStatus='ready'
            else
                CharTrack=VR-Phi_Circular(OmigaR)
                TrackStatus='nan_track'
            endif
        endif
    end subroutine

    subroutine GetTrackedW(WTrack,TrackReadyNow,TrackStatus)
        Implicit none
        Real(8),Intent(out)::WTrack
        Logical,Intent(out)::TrackReadyNow
        Character(32),Intent(out)::TrackStatus
        Real(8)::DenTrack

        WTrack=Wpre
        TrackReadyNow=.false.
        TrackStatus='not_ready'

        if (NegTrackReady) then
            if (abs(DiagIface-NegTrackIface)>1) then
                TrackStatus='stale_iface'
            else
                DenTrack=NegTrackOm2-NegTrackOm1
                if (abs(DenTrack)>Tol1*max(1.0d0,abs(NegTrackOm1),abs(NegTrackOm2))) then
                    WTrack=(NegTrackQ2-NegTrackQ1)/DenTrack
                    if (WTrack==WTrack) then
                        TrackReadyNow=.true.
                        TrackStatus='ready'
                    else
                        WTrack=Wpre
                        TrackStatus='nan_track'
                    endif
                else
                    TrackStatus='degenerate'
                endif
            endif
        endif
    end subroutine

    subroutine WriteNegativeRightCharShadowDiag()
        Implicit none
        Real(8)::CharRaw,CharTrack,TrackBestAbs,TrackBestH1,TrackBestH2,TrackBestEnergy
        Real(8)::SavedShadowRightCharState,TrackU2
        Integer::TrackSampleCount,TrackBracketCount
        Logical::TrackHasSignChange,TrackReadyNow
        Logical::SavedUseShadowRightChar
        Character(32)::TrackRootStatus

        if (.not.EnableNegativeShadowDiagnostics) return

        call SetNegativeModeControl(.true.,2)
        CharRaw=RightCharState
        TrackBestAbs=-1.0d0
        TrackBestH1=0.0d0
        TrackBestH2=0.0d0
        TrackBestEnergy=0.0d0
        TrackSampleCount=0
        TrackBracketCount=0
        TrackHasSignChange=.false.
        TrackReadyNow=.false.
        TrackRootStatus='not_ready'
        SavedUseShadowRightChar=UseShadowRightChar
        SavedShadowRightCharState=ShadowRightCharState

        call GetTrackedRightChar(CharTrack,TrackU2,TrackReadyNow,TrackRootStatus)
        if (TrackReadyNow) then
            UseShadowRightChar=.true.
            ShadowRightCharState=CharTrack
            call SetNegativeModeControl(.true.,2)
            call CheckNegativeExistence(TrackBestAbs,TrackBestH1,TrackBestH2,TrackBestEnergy, &
                TrackSampleCount,TrackBracketCount,TrackHasSignChange,TrackRootStatus)
            UseShadowRightChar=SavedUseShadowRightChar
            ShadowRightCharState=SavedShadowRightCharState
            call SetNegativeModeControl(.true.,2)
        endif

        Write(99,*) 'M2RCHAR',DiagStep,DiagTime,DiagIface,CharRaw,CharTrack,TrackReadyNow, &
            trim(NegativeRootStatus),RootBestAbs,RootSampleCount,RootBracketCount, &
            trim(TrackRootStatus),TrackBestAbs,TrackSampleCount,TrackBracketCount,NegTrackIface,NegTrackOm2,TrackU2
    end subroutine

    subroutine WriteNegativeRightCharTryDiag()
        Implicit none
        Real(8)::CharRaw,CharTrack,TrackU2,UMShadow(2),TryW,SavedW,SavedIntM(4)
        Real(8)::SavedShadowRightCharState
        Logical::TrackReadyNow,ShadowSolved,SavedUseShadowRightChar
        Character(32)::TryStatus

        if (.not.EnableNegativeShadowDiagnostics) return

        CharRaw=VR-Phi_Circular(OmigaR)
        CharTrack=CharRaw
        TrackU2=0.0d0
        TrackReadyNow=.false.
        ShadowSolved=.false.
        TryStatus='not_ready'
        TryW=Wpre
        UMShadow(1)=OmigaL
        UMShadow(2)=Qm_L
        SavedW=w
        SavedIntM=int_m
        SavedUseShadowRightChar=UseShadowRightChar
        SavedShadowRightCharState=ShadowRightCharState

        call GetTrackedRightChar(CharTrack,TrackU2,TrackReadyNow,TryStatus)
        if (TrackReadyNow) then
            UseShadowRightChar=.true.
            ShadowRightCharState=CharTrack
            call TryNegativeMode(.true.,2,'M2RFAIL','M2RSCAN',ShadowSolved,UMShadow)
            TryW=w
            TryStatus=NegativeExitReason
            UseShadowRightChar=SavedUseShadowRightChar
            ShadowRightCharState=SavedShadowRightCharState
        endif

        w=SavedW
        int_m=SavedIntM
        call SetNegativeModeControl(.true.,2)
        Write(99,*) 'M2RTRY',DiagStep,DiagTime,DiagIface,TrackReadyNow,trim(TryStatus),CharRaw,CharTrack, &
            ShadowSolved,TryW,UMShadow(1),UMShadow(2)
    end subroutine

    subroutine WriteNegativeWScanDiag(WTag)
        Implicit none
        Character(*),Intent(in)::WTag
        Integer,Parameter::NWScan=21
        Integer::IW,RootHitCount,SampleCountLoc,BracketCountLoc
        Integer::SavedNScanH1,SavedNScanH2,SavedMaxBisect
        Real(8)::WMin,WMax,WTrial,WBest,NormBest,H1Best,H2Best,EnergyBest
        Real(8)::WNearestRoot,NearestRootDelta,NormLoc,H1Loc,H2Loc,EnergyLoc
        Logical::HasSignChangeLoc
        Character(32)::RootStatusLoc

        call GetNegativeWWindow(1.0d0,WMin,WMax)
        RootHitCount=0
        WBest=Wpre
        NormBest=huge(1.0d0)
        H1Best=0.0d0
        H2Best=0.0d0
        EnergyBest=0.0d0
        WNearestRoot=0.0d0
        NearestRootDelta=-1.0d0

        SavedNScanH1=RootScanNScanH1
        SavedNScanH2=RootScanNScanH2
        SavedMaxBisect=RootScanMaxBisect
        RootScanNScanH1=21
        RootScanNScanH2=31
        RootScanMaxBisect=30

        do IW=0,NWScan-1
            WTrial=WMin+(WMax-WMin)*dble(IW)/dble(NWScan-1)
            WState=WTrial
            call CheckNegativeExistence(NormLoc,H1Loc,H2Loc,EnergyLoc,SampleCountLoc,BracketCountLoc, &
                HasSignChangeLoc,RootStatusLoc)
            if (NormLoc>=0.0d0 .and. NormLoc<NormBest) then
                NormBest=NormLoc
                WBest=WTrial
                H1Best=H1Loc
                H2Best=H2Loc
                EnergyBest=EnergyLoc
            endif
            if (trim(RootStatusLoc)=='has_root') then
                RootHitCount=RootHitCount+1
                if (NearestRootDelta<0.0d0 .or. abs(WTrial-Wpre)<NearestRootDelta) then
                    NearestRootDelta=abs(WTrial-Wpre)
                    WNearestRoot=WTrial
                endif
            endif
        enddo

        RootScanNScanH1=SavedNScanH1
        RootScanNScanH2=SavedNScanH2
        RootScanMaxBisect=SavedMaxBisect
        WState=Wpre

        if (NormBest==huge(1.0d0)) NormBest=-1.0d0
        Mode2WRootHitCount=RootHitCount
        Mode2WNearestRoot=WNearestRoot
        Mode2WNearestDelta=NearestRootDelta
        if (EnableNegativeShadowDiagnostics) then
            Write(99,*) trim(WTag),RootHitCount,Wpre,WMin,WMax,WNearestRoot,NearestRootDelta, &
                WBest,NormBest,H1Best,H2Best,EnergyBest
        endif
    end subroutine

    subroutine FindMode2ClosestRetryW(WRetry,FoundW)
        Implicit none
        Real(8),Intent(out)::WRetry
        Logical,Intent(out)::FoundW
        Integer,Parameter::NWRefine=41
        Integer::IW,SampleCountLoc,BracketCountLoc
        Integer::SavedNScanH1,SavedNScanH2,SavedMaxBisect
        Real(8)::WTrial,NormLoc,H1Loc,H2Loc,EnergyLoc
        Logical::HasSignChangeLoc
        Character(32)::RootStatusLoc

        WRetry=Mode2WNearestRoot
        FoundW=.false.
        if (Mode2WRootHitCount<=0 .or. Mode2WNearestDelta<0.0d0) return

        SavedNScanH1=RootScanNScanH1
        SavedNScanH2=RootScanNScanH2
        SavedMaxBisect=RootScanMaxBisect
        RootScanNScanH1=31
        RootScanNScanH2=41
        RootScanMaxBisect=40

        do IW=1,NWRefine
            WTrial=Wpre+(Mode2WNearestRoot-Wpre)*dble(IW)/dble(NWRefine)
            WState=WTrial
            call CheckNegativeExistence(NormLoc,H1Loc,H2Loc,EnergyLoc,SampleCountLoc,BracketCountLoc, &
                HasSignChangeLoc,RootStatusLoc)
            if (trim(RootStatusLoc)=='has_root') then
                WRetry=WTrial
                FoundW=.true.
                exit
            endif
        enddo

        RootScanNScanH1=SavedNScanH1
        RootScanNScanH2=SavedNScanH2
        RootScanMaxBisect=SavedMaxBisect
        WState=Wpre
        if (.not.FoundW) WRetry=Mode2WNearestRoot
    end subroutine

    subroutine TryNegativeMode(UseRightCharEqWanted,InterTypeWanted,FailTag,ScanTag,ModeSolved,UMMode)
        Implicit none
        Logical,Intent(in)::UseRightCharEqWanted
        Integer,Intent(in)::InterTypeWanted
        Character(*),Intent(in)::FailTag,ScanTag
        Logical,Intent(out)::ModeSolved
        Real(8),Intent(out)::UMMode(2)
        Real(8)::SeedNorm,TrackSeedH1,TrackSeedH2
        Real(8)::CurveH1,CurveH2,CurveEnergy,WRetry
        Logical::SeedStateOK,CurveSolved,RetryWFound

        call SetNegativeModeControl(UseRightCharEqWanted,InterTypeWanted)

        ModeSolved=.false.
        UMMode(1)=OmigaL
        UMMode(2)=Qm_L
        w=Wpre
        Det=0.0d0
        Delta1=0.0d0
        Delta2=0.0d0
        Lambda=0.0d0
        TrialNorm=-1.0d0
        ResNorm=-1.0d0
        ResNormInit=-1.0d0
        IterLast=0
        LSIterLast=0
        NegativeExitReason='init1'
        H1=max(HL,y_ref+Tol1)
        H2=min(max(HR,Tol1),y_ref-Tol1)
        if (H2<=Tol1) H2=max(10.0d0*Tol1,0.5d0*y_ref)

        call EvaluateNegativeState(H1,H2,R1,R2,Omiga1,U1,Omiga2,U2,ValidState)
        if (.not.ValidState) then
            NegativeExitReason='init1_invalid'
            H1=y_ref+max(Tol1,0.1d0*D)
            H2=max(10.0d0*Tol1,min(0.5d0*y_ref,y_ref-Tol1))
            call EvaluateNegativeState(H1,H2,R1,R2,Omiga1,U1,Omiga2,U2,ValidState)
            if (.not.ValidState) NegativeExitReason='init2_invalid'
        endif
        if (ValidState) ResNormInit=max(abs(R1)/Scale1,abs(R2)/Scale2)
        if (NegTrackReady) then
            if (NegTrackOm1>0.0d0 .and. NegTrackOm2>0.0d0) then
                TrackSeedH1=max(Depth(NegTrackOm1),y_ref+Tol1)
                TrackSeedH2=min(max(Depth(NegTrackOm2),10.0d0*Tol1),y_ref-Tol1)
                call EvaluateNegativeState(TrackSeedH1,TrackSeedH2,R1P,R2P,Omiga1T,U1T,Omiga2T,U2T,SeedStateOK)
                if (SeedStateOK) then
                    SeedNorm=max(abs(R1P)/Scale1,abs(R2P)/Scale2)
                    if ((.not.ValidState) .or. ResNormInit<0.0d0 .or. SeedNorm<ResNormInit) then
                        H1=TrackSeedH1
                        H2=TrackSeedH2
                        R1=R1P
                        R2=R2P
                        Omiga1=Omiga1T
                        Omiga2=Omiga2T
                        U1=U1T
                        U2=U2T
                        ValidState=.true.
                        ResNormInit=SeedNorm
                    endif
                endif
            endif
        endif
        if (HasS1Seed) then
            call EvaluateNegativeState(S1SeedH1,S1SeedH2,R1P,R2P,Omiga1T,U1T,Omiga2T,U2T,SeedStateOK)
            if (SeedStateOK) then
                SeedNorm=max(abs(R1P)/Scale1,abs(R2P)/Scale2)
                if ((.not.ValidState) .or. ResNormInit<0.0d0 .or. SeedNorm<ResNormInit) then
                    H1=S1SeedH1
                    H2=S1SeedH2
                    R1=R1P
                    R2=R2P
                    Omiga1=Omiga1T
                    Omiga2=Omiga2T
                    U1=U1T
                    U2=U2T
                    ValidState=.true.
                    ResNormInit=SeedNorm
                endif
            endif
        endif

        if (ValidState) then
            NegativeExitReason='iterating'
            do Iter=1,60
                IterLast=Iter
                ResNorm=max(abs(R1)/Scale1,abs(R2)/Scale2)
                if (ResNorm<1.0d-8) then
                    ModeSolved=.true.
                    NegativeExitReason='conv_1e8'
                    exit
                endif

                DH1=max(1.0d-6,1.0d-6*abs(H1))
                DH2=max(1.0d-6,1.0d-6*abs(H2))
                if (H2+DH2>=y_ref-Tol1) DH2=-DH2
                if (H2+DH2<=Tol1) DH2=0.5d0*H2

                call EvaluateNegativeState(H1+DH1,H2,R1P,R2P,Omiga1T,U1T,Omiga2T,U2T,ValidState)
                if (.not.ValidState) then
                    NegativeExitReason='jac_h1_invalid'
                    exit
                endif
                J11=(R1P-R1)/DH1
                J21=(R2P-R2)/DH1
                call EvaluateNegativeState(H1,H2+DH2,R1P,R2P,Omiga1T,U1T,Omiga2T,U2T,ValidState)
                if (.not.ValidState) then
                    NegativeExitReason='jac_h2_invalid'
                    exit
                endif
                J12=(R1P-R1)/DH2
                J22=(R2P-R2)/DH2
                Det=J11*J22-J12*J21
                if (abs(Det)<=1.0d-12) then
                    NegativeExitReason='det_small'
                    exit
                endif

                Delta1=(-R1*J22+R2*J12)/Det
                Delta2=(R1*J21-J11*R2)/Det
                Lambda=1.0d0
                AcceptStep=.false.
                LSIterLast=0

                do LSIter=1,20
                    LSIterLast=LSIter
                    H1T=H1+Lambda*Delta1
                    H2T=H2+Lambda*Delta2
                    if (H1T<=y_ref+Tol1 .or. H2T<=Tol1 .or. H2T>=y_ref-Tol1) then
                        Lambda=0.5d0*Lambda
                        cycle
                    endif

                    call EvaluateNegativeState(H1T,H2T,R1P,R2P,Omiga1T,U1T,Omiga2T,U2T,ValidState)
                    if (.not.ValidState) then
                        Lambda=0.5d0*Lambda
                        cycle
                    endif

                    TrialNorm=max(abs(R1P)/Scale1,abs(R2P)/Scale2)
                    if (TrialNorm<ResNorm) then
                        H1=H1T
                        H2=H2T
                        R1=R1P
                        R2=R2P
                        Omiga1=Omiga1T
                        Omiga2=Omiga2T
                        U1=U1T
                        U2=U2T
                        AcceptStep=.true.
                        exit
                    endif
                    Lambda=0.5d0*Lambda
                enddo

                if (.not.AcceptStep) then
                    NegativeExitReason='line_search_fail'
                    exit
                endif
            enddo

            ResNorm=max(abs(R1)/Scale1,abs(R2)/Scale2)
            if (ResNorm<1.0d-5) then
                ModeSolved=.true.
                if (NegativeExitReason/='conv_1e8') NegativeExitReason='conv_1e5'
            elseif (NegativeExitReason=='iterating') then
                NegativeExitReason='iter_limit'
            endif
        endif

        if (.not.ModeSolved) then
            call SolveNegativeMode1Root(CurveH1,CurveH2,CurveEnergy,CurveSolved)
            if (CurveSolved) then
                call EvaluateNegativeState(CurveH1,CurveH2,R1,R2,Omiga1,U1,Omiga2,U2,ValidState)
                if (ValidState) then
                    H1=CurveH1
                    H2=CurveH2
                    ResNorm=max(abs(R1)/Scale1,abs(R2)/Scale2)
                    if (ResNorm<1.0d-5) then
                        ModeSolved=.true.
                        NegativeExitReason='curve_root'
                    endif
                endif
            endif
        endif

        if (.not.ModeSolved .and. Mode2WRootHitCount>0 .and. Mode2WNearestDelta>=0.0d0) then
            call FindMode2ClosestRetryW(WRetry,RetryWFound)
            WState=WRetry
            call SolveNegativeMode1Root(CurveH1,CurveH2,CurveEnergy,CurveSolved)
            if (CurveSolved) then
                call EvaluateNegativeState(CurveH1,CurveH2,R1,R2,Omiga1,U1,Omiga2,U2,ValidState)
                if (ValidState) then
                    H1=CurveH1
                    H2=CurveH2
                    ResNorm=max(abs(R1)/Scale1,abs(R2)/Scale2)
                    if (ResNorm<1.0d-5) then
                        ModeSolved=.true.
                        w=WState
                        NegativeExitReason='wscan_root'
                    endif
                endif
            endif
            if (.not.ModeSolved) then
                WState=Wpre
                w=Wpre
            endif
        endif

        if (ModeSolved) then
            int_m=(/Omiga1,Omiga1*U1,Omiga2,Omiga2*U2/)
            call SelectNegativeGodunovState(Omiga2,U2,UMMode,ValidState)
            if (.not.ValidState) then
                UMMode(1)=Omiga2
                UMMode(2)=Omiga2*U2
            endif
        else
            if (EnableInterfaceDebugOutput .or. EnableNegativeShadowDiagnostics) then
                call ScanNegativeResiduals(ScanH1Min,ScanH1Max,ScanH2Min,ScanH2Max, &
                    ScanBestNorm,ScanBestH1,ScanBestH2,ScanBestR1,ScanBestR2,ScanValidCount)
                Write(99,*) trim(FailTag),NegativeExitReason,IterLast,LSIterLast,ResNormInit,ResNorm, &
                    R1,R2,H1,H2,Det,Delta1,Delta2,Lambda,TrialNorm,Wpre,HL,HR,VL,VR
                Write(99,*) trim(ScanTag),ScanValidCount,ScanH1Min,ScanH1Max,ScanH2Min,ScanH2Max, &
                    ScanBestNorm,ScanBestH1,ScanBestH2,ScanBestR1,ScanBestR2
            endif
        endif
    end subroutine

    subroutine EvaluatePositiveState(HPVal,OmigaPVal,KVal)
        Implicit none
        Real(8),Intent(in)::HPVal
        Real(8),Intent(out)::OmigaPVal,KVal
        Real(8),External::APY,Density,Area

        OmigaPVal=Density(HPVal)*Area(HPVal)
        KVal=APY(OmigaPVal)
    end subroutine

    subroutine EvaluateNegativeState(H1P,H2P,RMom,REne,Omiga1Loc,U1Loc,Omiga2Loc,U2Loc,StateOK)
        Implicit none
        Real(8),Intent(in)::H1P,H2P
        Real(8),Intent(out)::RMom,REne,Omiga1Loc,U1Loc,Omiga2Loc,U2Loc
        Logical,Intent(out)::StateOK
        Real(8)::JRel
        Real(8),External::APY,Density,Area,Phi_Circular

        StateOK=.false.
        RMom=0.0d0
        REne=0.0d0
        Omiga1Loc=0.0d0
        Omiga2Loc=0.0d0
        U1Loc=0.0d0
        U2Loc=0.0d0

        if (H1P<=y_ref+Tol1) return
        if (H2P<=Tol1 .or. H2P>=y_ref-Tol1) return

        Omiga1Loc=Density(H1P)*Area(H1P)
        Omiga2Loc=Density(H2P)*Area(H2P)
        if (Omiga1Loc<=0.0d0 .or. Omiga2Loc<=0.0d0) return

        U1Loc=VL+C_Waterhammer*LOG(OmigaL/Omiga1Loc)
        JRel=Omiga1Loc*(U1Loc-WState)
        U2Loc=WState+JRel/Omiga2Loc

        RMom=JRel*JRel/Omiga1Loc+APY(Omiga1Loc)-JRel*JRel/Omiga2Loc-APY(Omiga2Loc)
        if (UseRightCharEq) then
            REne=U2Loc-Phi_Circular(Omiga2Loc)-RightCharState
        else
            REne=Omiga1Loc*(U1Loc-WState)**2+2.0d0*Omiga1Loc*g*H1P- &
                Omiga2Loc*(U2Loc-WState)**2-2.0d0*Omiga2Loc*g*H2P
        endif
        if (RMom/=RMom .or. REne/=REne) return

        StateOK=.true.
    end subroutine

    subroutine SelectNegativeGodunovState(Omiga2Loc,U2Loc,UMSel,StateOK)
        Implicit none
        Real(8),Intent(in)::Omiga2Loc,U2Loc
        Real(8),Intent(out)::UMSel(2)
        Logical,Intent(out)::StateOK
        Real(8)::Q2Loc,CLoc,RLoc,VRight,UStar(2)

        Q2Loc=Omiga2Loc*U2Loc
        UMSel(1)=Omiga2Loc
        UMSel(2)=Q2Loc
        StateOK=.false.
        if (Omiga2Loc<=0.0d0 .or. OmigaR<=0.0d0) return

        if (UseRightCharEq) then
            StateOK=.true.
            return
        endif

        CLoc=Celerity(Depth(Omiga2Loc))
        RLoc=Celerity(HR)
        VRight=VR

        if (VRight-RLoc>=0.0d0) then
            UMSel(1)=OmigaR
            UMSel(2)=Qm_R
            StateOK=.true.
        elseif (U2Loc-CLoc<=0.0d0) then
            StateOK=.true.
        else
            call SolveFreeSurfaceStar(Omiga2Loc,OmigaR,Q2Loc,Qm_R,UStar,StateOK)
            if (StateOK) then
                UMSel(1)=UStar(1)
                UMSel(2)=UStar(2)
            endif
        endif
    end subroutine

    subroutine SolveFreeSurfaceStar(U1L,U1R,U2L,U2R,UStar,StateOK)
        Implicit none
        Real(8),Intent(in)::U1L,U1R,U2L,U2R
        Real(8),Intent(out)::UStar(2)
        Logical,Intent(out)::StateOK
        Real(8)::VLfs,VRfs,m,n,mn,Fm,Fn,Fmn,Omiga_m,Omiga_n,Omiga_mn,OmigaM,QmM
        Integer::IterFS
        Real(8),External::Area,Density,F_Equation

        StateOK=.false.
        UStar(1)=U1L
        UStar(2)=U2L
        if (U1L<=0.0d0 .or. U1R<=0.0d0) return

        VLfs=U2L/U1L
        VRfs=U2R/U1R
        m=1.0d-5
        n=max(2.0d0*Tol1,y_ref-Tol1)
        mn=0.5d0*(m+n)

        do IterFS=1,200
            mn=0.5d0*(m+n)
            Omiga_m=Area(m)*Density(m)
            Omiga_n=Area(n)*Density(n)
            Omiga_mn=Area(mn)*Density(mn)
            Fm=F_Equation(Omiga_m,U1L)+F_Equation(Omiga_m,U1R)+VRfs-VLfs
            Fn=F_Equation(Omiga_n,U1L)+F_Equation(Omiga_n,U1R)+VRfs-VLfs
            Fmn=F_Equation(Omiga_mn,U1L)+F_Equation(Omiga_mn,U1R)+VRfs-VLfs
            if (Fm/=Fm .or. Fn/=Fn .or. Fmn/=Fmn) return
            if (Fm*Fn>0.0d0) return

            if (Fm*Fmn<0.0d0) then
                n=mn
            else
                m=mn
            endif

            if (abs(m-n)<1.0d-5) exit
        enddo

        OmigaM=Area(mn)*Density(mn)
        QmM=0.5d0*OmigaM*(VRfs+VLfs)+0.5d0*OmigaM*(F_Equation(OmigaM,U1R)-F_Equation(OmigaM,U1L))
        UStar(1)=OmigaM
        UStar(2)=QmM
        StateOK=.true.
    end subroutine

    subroutine ScanNegativeResiduals(H1Min,H1Max,H2Min,H2Max,BestNorm,BestH1,BestH2,BestR1,BestR2,ValidCount)
        Implicit none
        Real(8),Intent(out)::H1Min,H1Max,H2Min,H2Max,BestNorm,BestH1,BestH2,BestR1,BestR2
        Integer,Intent(out)::ValidCount
        Integer,Parameter::NScan1=21,NScan2=21
        Integer::I1,I2
        Real(8)::H1S,H2S,R1S,R2S,Omiga1S,Omiga2S,U1S,U2S,NormS
        Logical::StateOK

        H1Min=y_ref+max(Tol1,1.0d-3*D)
        H1Max=max(HL+abs(HL-HR)+2.0d0*D,2.0d0*max(HL,y_ref),y_ref+3.0d0*D)
        if (H1Max<=H1Min) H1Max=H1Min+max(1.0d-3,0.5d0*D)
        H2Min=max(10.0d0*Tol1,1.0d-5)
        H2Max=y_ref-Tol1
        if (H2Max<=H2Min) H2Min=0.5d0*H2Max

        BestNorm=huge(1.0d0)
        BestH1=0.0d0
        BestH2=0.0d0
        BestR1=0.0d0
        BestR2=0.0d0
        ValidCount=0

        do I1=0,NScan1-1
            H1S=H1Min+(H1Max-H1Min)*dble(I1)/dble(NScan1-1)
            do I2=0,NScan2-1
                H2S=H2Min+(H2Max-H2Min)*dble(I2)/dble(NScan2-1)
                call EvaluateNegativeState(H1S,H2S,R1S,R2S,Omiga1S,U1S,Omiga2S,U2S,StateOK)
                if (.not.StateOK) cycle
                NormS=max(abs(R1S)/Scale1,abs(R2S)/Scale2)
                ValidCount=ValidCount+1
                if (NormS<BestNorm) then
                    BestNorm=NormS
                    BestH1=H1S
                    BestH2=H2S
                    BestR1=R1S
                    BestR2=R2S
                endif
            enddo
        enddo

        if (ValidCount==0) BestNorm=-1.0d0
    end subroutine

    subroutine SolveMomentumRootForH1(H1Val,H2Best,EnergyBest,FoundRoot,BracketCount)
        Implicit none
        Real(8),Intent(in)::H1Val
        Real(8),Intent(out)::H2Best,EnergyBest
        Logical,Intent(out)::FoundRoot
        Integer,Intent(out)::BracketCount
        Integer::I2,IBis,NScanH2,MaxBisect
        Real(8)::H2Min,H2Max,H2Prev,H2Cur,H2Left,H2Right,H2Mid
        Real(8)::RMomPrev,REnePrev,RMomCur,REneCur,RMomMid,REneMid,RMomLeft,RMomRight
        Real(8)::Omiga1Tmp,Omiga2Tmp,U1Tmp,U2Tmp,BestNormE,NormE
        Logical::PrevValid,StateOK,MidOK

        NScanH2=max(3,RootScanNScanH2)
        MaxBisect=max(5,RootScanMaxBisect)
        H2Best=0.0d0
        EnergyBest=0.0d0
        FoundRoot=.false.
        BracketCount=0
        H2Min=max(10.0d0*Tol1,1.0d-5)
        H2Max=y_ref-Tol1
        if (H2Max<=H2Min) return

        BestNormE=huge(1.0d0)
        PrevValid=.false.

        do I2=0,NScanH2-1
            H2Cur=H2Min+(H2Max-H2Min)*dble(I2)/dble(NScanH2-1)
            call EvaluateNegativeState(H1Val,H2Cur,RMomCur,REneCur,Omiga1Tmp,U1Tmp,Omiga2Tmp,U2Tmp,StateOK)
            if (.not.StateOK) then
                PrevValid=.false.
                cycle
            endif

            if (abs(RMomCur)<=1.0d-10*Scale1) then
                BracketCount=BracketCount+1
                NormE=abs(REneCur)/Scale2
                if (NormE<BestNormE) then
                    BestNormE=NormE
                    H2Best=H2Cur
                    EnergyBest=REneCur
                    FoundRoot=.true.
                endif
            elseif (PrevValid .and. RMomPrev*RMomCur<0.0d0) then
                BracketCount=BracketCount+1
                H2Left=H2Prev
                H2Right=H2Cur
                RMomLeft=RMomPrev
                RMomRight=RMomCur
                H2Mid=0.5d0*(H2Left+H2Right)
                RMomMid=0.0d0
                REneMid=0.0d0
                MidOK=.true.

                do IBis=1,MaxBisect
                    H2Mid=0.5d0*(H2Left+H2Right)
                    call EvaluateNegativeState(H1Val,H2Mid,RMomMid,REneMid,Omiga1Tmp,U1Tmp,Omiga2Tmp,U2Tmp,MidOK)
                    if (.not.MidOK) exit
                    if (abs(RMomMid)<=1.0d-10*Scale1 .or. abs(H2Right-H2Left)<=1.0d-8) exit
                    if (RMomLeft*RMomMid<=0.0d0) then
                        H2Right=H2Mid
                        RMomRight=RMomMid
                    else
                        H2Left=H2Mid
                        RMomLeft=RMomMid
                    endif
                enddo

                if (MidOK) then
                    NormE=abs(REneMid)/Scale2
                    if (NormE<BestNormE) then
                        BestNormE=NormE
                        H2Best=H2Mid
                        EnergyBest=REneMid
                        FoundRoot=.true.
                    endif
                endif
            endif

            H2Prev=H2Cur
            RMomPrev=RMomCur
            REnePrev=REneCur
            PrevValid=.true.
        enddo
    end subroutine

    subroutine CheckNegativeExistence(BestNormE,BestH1,BestH2,BestEnergy,SampleCount,BracketCount,HasSignChange,RootStatus)
        Implicit none
        Real(8),Intent(out)::BestNormE,BestH1,BestH2,BestEnergy
        Integer,Intent(out)::SampleCount,BracketCount
        Logical,Intent(out)::HasSignChange
        Character(32),Intent(out)::RootStatus
        Integer::I1,BracketCountLoc
        Integer::NScanH1
        Real(8)::H1Min,H1Max,H1Val,H2Root,EnergyRoot,NormE,PrevEnergy
        Logical::FoundRoot,PrevValid,HasExactRoot

        NScanH1=max(3,RootScanNScanH1)
        H1Min=y_ref+max(Tol1,1.0d-3*D)
        H1Max=max(HL+abs(HL-HR)+2.0d0*D,2.0d0*max(HL,y_ref),y_ref+3.0d0*D)
        if (H1Max<=H1Min) H1Max=H1Min+max(1.0d-3,0.5d0*D)

        BestNormE=huge(1.0d0)
        BestH1=0.0d0
        BestH2=0.0d0
        BestEnergy=0.0d0
        SampleCount=0
        BracketCount=0
        HasSignChange=.false.
        HasExactRoot=.false.
        PrevValid=.false.
        RootStatus='inconclusive'

        do I1=0,NScanH1-1
            H1Val=H1Min+(H1Max-H1Min)*dble(I1)/dble(NScanH1-1)
            call SolveMomentumRootForH1(H1Val,H2Root,EnergyRoot,FoundRoot,BracketCountLoc)
            BracketCount=BracketCount+BracketCountLoc
            if (.not.FoundRoot) then
                PrevValid=.false.
                cycle
            endif

            SampleCount=SampleCount+1
            NormE=abs(EnergyRoot)/Scale2
            if (NormE<BestNormE) then
                BestNormE=NormE
                BestH1=H1Val
                BestH2=H2Root
                BestEnergy=EnergyRoot
            endif
            if (NormE<=1.0d-6) HasExactRoot=.true.

            if (PrevValid) then
                if (PrevEnergy==0.0d0 .or. EnergyRoot==0.0d0 .or. PrevEnergy*EnergyRoot<0.0d0) then
                    HasSignChange=.true.
                endif
            endif
            PrevEnergy=EnergyRoot
            PrevValid=.true.
        enddo

        if (SampleCount==0) then
            BestNormE=-1.0d0
            return
        endif

        if (HasExactRoot .or. HasSignChange) then
            RootStatus='has_root'
        elseif (SampleCount>=3) then
            RootStatus='no_interior_root'
        endif
    end subroutine

    subroutine SolveNegativeMode1Root(H1Root,H2Root,EnergyRoot,FoundRoot)
        Implicit none
        Real(8),Intent(out)::H1Root,H2Root,EnergyRoot
        Logical,Intent(out)::FoundRoot
        Integer::I1,IBis,NScanH1,MaxBisect,BracketCountLoc
        Real(8)::H1Min,H1Max,H1Prev,H1Cur,H1Left,H1Right,H1Mid
        Real(8)::H2Cur,H2Mid,EnergyPrev,EnergyCur,EnergyLeft,EnergyMid,BestNormE,NormE
        Logical::PrevValid,StateOK

        NScanH1=max(5,2*RootScanNScanH1-1)
        MaxBisect=max(10,RootScanMaxBisect)
        H1Min=y_ref+max(Tol1,1.0d-3*D)
        H1Max=max(HL+abs(HL-HR)+2.0d0*D,2.0d0*max(HL,y_ref),y_ref+3.0d0*D)
        if (H1Max<=H1Min) H1Max=H1Min+max(1.0d-3,0.5d0*D)

        H1Root=0.0d0
        H2Root=0.0d0
        EnergyRoot=0.0d0
        FoundRoot=.false.
        BestNormE=huge(1.0d0)
        PrevValid=.false.

        do I1=0,NScanH1-1
            H1Cur=H1Min+(H1Max-H1Min)*dble(I1)/dble(NScanH1-1)
            call SolveMomentumRootForH1(H1Cur,H2Cur,EnergyCur,StateOK,BracketCountLoc)
            if (.not.StateOK) then
                PrevValid=.false.
                cycle
            endif

            NormE=abs(EnergyCur)/Scale2
            if (NormE<BestNormE) then
                BestNormE=NormE
                H1Root=H1Cur
                H2Root=H2Cur
                EnergyRoot=EnergyCur
            endif
            if (NormE<=1.0d-8) then
                FoundRoot=.true.
                return
            endif

            if (PrevValid) then
                if (EnergyPrev==0.0d0 .or. EnergyCur==0.0d0 .or. EnergyPrev*EnergyCur<0.0d0) then
                    H1Left=H1Prev
                    H1Right=H1Cur
                    EnergyLeft=EnergyPrev
                    H1Mid=0.5d0*(H1Left+H1Right)
                    H2Mid=H2Cur
                    EnergyMid=EnergyCur
                    StateOK=.true.

                    do IBis=1,MaxBisect
                        H1Mid=0.5d0*(H1Left+H1Right)
                        call SolveMomentumRootForH1(H1Mid,H2Mid,EnergyMid,StateOK,BracketCountLoc)
                        if (.not.StateOK) exit
                        if (abs(EnergyMid)<=1.0d-8*Scale2 .or. abs(H1Right-H1Left)<=1.0d-8) exit
                        if (EnergyLeft==0.0d0 .or. EnergyLeft*EnergyMid<=0.0d0) then
                            H1Right=H1Mid
                        else
                            H1Left=H1Mid
                            EnergyLeft=EnergyMid
                        endif
                    enddo

                    if (StateOK) then
                        H1Root=H1Mid
                        H2Root=H2Mid
                        EnergyRoot=EnergyMid
                        FoundRoot=.true.
                        return
                    endif
                endif
            endif

            H1Prev=H1Cur
            EnergyPrev=EnergyCur
            PrevValid=.true.
        enddo

        if (BestNormE<=1.0d-5) FoundRoot=.true.
    end subroutine

    subroutine InitializeNegativeStrategy2()
        Implicit none

        Inter_type=2
        UseRightCharEq=.true.
        WState=Wpre
        Scale2=max(1.0d0,abs(VL),abs(VR),abs(Phi_Circular(OmigaL)),abs(Phi_Circular(OmigaR)), &
            abs(VR-Phi_Circular(OmigaR)))
        Scale3=max(1.0d0,abs(OmigaL*VL*VL+2.0d0*OmigaL*g*HL), &
            abs(OmigaR*VR*VR+2.0d0*OmigaR*g*HR),abs(HgasLocal))
        Scale3=max(Scale3,abs(APY(OmigaL)),abs(APY(OmigaR)))
    end subroutine

    subroutine EvaluateNegativeStateS2(H1P,H2P,WLoc,RMom,RChar,REne,Omiga1Loc,U1Loc,Omiga2Loc,U2Loc,StateOK)
        Implicit none
        Real(8),Intent(in)::H1P,H2P,WLoc
        Real(8),Intent(out)::RMom,RChar,REne,Omiga1Loc,U1Loc,Omiga2Loc,U2Loc
        Logical,Intent(out)::StateOK
        Real(8)::JRel
        Real(8),External::APY,Density,Area,Phi_Circular

        StateOK=.false.
        RMom=0.0d0
        RChar=0.0d0
        REne=0.0d0
        Omiga1Loc=0.0d0
        Omiga2Loc=0.0d0
        U1Loc=0.0d0
        U2Loc=0.0d0

        if (H1P<=y_ref+Tol1) return
        if (H2P<=Tol1 .or. H2P>=y_ref-Tol1) return

        Omiga1Loc=Density(H1P)*Area(H1P)
        Omiga2Loc=Density(H2P)*Area(H2P)
        if (Omiga1Loc<=0.0d0 .or. Omiga2Loc<=0.0d0) return

        U1Loc=VL+C_Waterhammer*LOG(OmigaL/Omiga1Loc)
        JRel=Omiga1Loc*(U1Loc-WLoc)
        U2Loc=WLoc+JRel/Omiga2Loc

        RMom=JRel*JRel/Omiga1Loc+APY(Omiga1Loc)-JRel*JRel/Omiga2Loc-APY(Omiga2Loc)
        RChar=U2Loc-Phi_Circular(Omiga2Loc)-VR+Phi_Circular(OmigaR)
        REne=Omiga1Loc*(U1Loc-WLoc)**2+2.0d0*Omiga1Loc*g*H1P- &
            Omiga2Loc*(U2Loc-WLoc)**2-2.0d0*Omiga2Loc*g*H2P
        if (RMom/=RMom .or. RChar/=RChar .or. REne/=REne) return

        StateOK=.true.
    end subroutine

    subroutine SolveLinearSystem3x3(A11,A12,A13,A21,A22,A23,A31,A32,A33,B1,B2,B3,X1,X2,X3,SolveOK)
        Implicit none
        Real(8),Intent(in)::A11,A12,A13,A21,A22,A23,A31,A32,A33,B1,B2,B3
        Real(8),Intent(out)::X1,X2,X3
        Logical,Intent(out)::SolveOK
        Real(8)::DetA,Det1,Det2,Det3

        X1=0.0d0
        X2=0.0d0
        X3=0.0d0
        SolveOK=.false.
        DetA=A11*(A22*A33-A23*A32)-A12*(A21*A33-A23*A31)+A13*(A21*A32-A22*A31)
        if (abs(DetA)<=1.0d-12) return
        Det1=B1*(A22*A33-A23*A32)-A12*(B2*A33-A23*B3)+A13*(B2*A32-A22*B3)
        Det2=A11*(B2*A33-A23*B3)-B1*(A21*A33-A23*A31)+A13*(A21*B3-B2*A31)
        Det3=A11*(A22*B3-B2*A32)-A12*(A21*B3-B2*A31)+B1*(A21*A32-A22*A31)
        X1=Det1/DetA
        X2=Det2/DetA
        X3=Det3/DetA
        SolveOK=.true.
    end subroutine

    subroutine ScanNegativeResidualsS2(H1Min,H1Max,H2Min,H2Max,WMin,WMax,BestNorm,BestH1,BestH2,BestW, &
        BestR1,BestR2,BestR3,ValidCount)
        Implicit none
        Real(8),Intent(out)::H1Min,H1Max,H2Min,H2Max,WMin,WMax,BestNorm,BestH1,BestH2,BestW, &
            BestR1,BestR2,BestR3
        Integer,Intent(out)::ValidCount
        Integer,Parameter::NScan1=11,NScan2=11,NScanW=21
        Integer::I1,I2,IW
        Real(8)::H1S,H2S,WS,R1S,R2S,R3S,Omiga1S,Omiga2S,U1S,U2S,NormS
        Logical::StateOK

        H1Min=y_ref+max(Tol1,1.0d-3*D)
        H1Max=max(HL+abs(HL-HR)+2.0d0*D,2.0d0*max(HL,y_ref),y_ref+3.0d0*D)
        if (H1Max<=H1Min) H1Max=H1Min+max(1.0d-3,0.5d0*D)
        H2Min=max(10.0d0*Tol1,1.0d-5)
        H2Max=y_ref-Tol1
        if (H2Max<=H2Min) H2Min=0.5d0*H2Max
        call GetNegativeWWindow(1.0d0,WMin,WMax)

        BestNorm=huge(1.0d0)
        BestH1=0.0d0
        BestH2=0.0d0
        BestW=0.0d0
        BestR1=0.0d0
        BestR2=0.0d0
        BestR3=0.0d0
        ValidCount=0

        do IW=0,NScanW-1
            WS=WMin+(WMax-WMin)*dble(IW)/dble(NScanW-1)
            do I1=0,NScan1-1
                H1S=H1Min+(H1Max-H1Min)*dble(I1)/dble(NScan1-1)
                do I2=0,NScan2-1
                    H2S=H2Min+(H2Max-H2Min)*dble(I2)/dble(NScan2-1)
                    call EvaluateNegativeStateS2(H1S,H2S,WS,R1S,R2S,R3S,Omiga1S,U1S,Omiga2S,U2S,StateOK)
                    if (.not.StateOK) cycle
                    NormS=max(abs(R1S)/Scale1,abs(R2S)/Scale2,abs(R3S)/Scale3)
                    ValidCount=ValidCount+1
                    if (NormS<BestNorm) then
                        BestNorm=NormS
                        BestH1=H1S
                        BestH2=H2S
                        BestW=WS
                        BestR1=R1S
                        BestR2=R2S
                        BestR3=R3S
                    endif
                enddo
            enddo
        enddo

        if (ValidCount==0) BestNorm=-1.0d0
    end subroutine

    subroutine WriteNegativeRootDiagS2(RootTag)
        Implicit none
        Character(*),Intent(in)::RootTag
        Real(8)::H1Min,H1Max,H2Min,H2Max,WMin,WMax,BestNorm,BestH1,BestH2,BestW,BestR1,BestR2,BestR3
        Integer::ValidCount
        Character(32)::RootStatusLoc

        call InitializeNegativeStrategy2()
        call ScanNegativeResidualsS2(H1Min,H1Max,H2Min,H2Max,WMin,WMax,BestNorm,BestH1,BestH2,BestW, &
            BestR1,BestR2,BestR3,ValidCount)
        HasS2Seed=.false.
        S2SeedNorm=-1.0d0
        if (ValidCount>0 .and. BestNorm>=0.0d0) then
            HasS2Seed=.true.
            S2SeedH1=BestH1
            S2SeedH2=BestH2
            S2SeedW=BestW
            S2SeedNorm=BestNorm
        endif
        if (ValidCount==0) then
            RootStatusLoc='inconclusive'
        elseif (BestNorm<=1.0d-6) then
            RootStatusLoc='has_root'
        else
            RootStatusLoc='no_grid_root'
        endif
        Write(99,*) trim(RootTag),trim(RootStatusLoc),ValidCount,0,.false.,BestNorm,BestH1,BestH2, &
            BestW,BestR1,BestR2,BestR3
    end subroutine

    subroutine TryNegativeStrategy2(FailTag,ScanTag,ModeSolved,UMMode)
        Implicit none
        Character(*),Intent(in)::FailTag,ScanTag
        Logical,Intent(out)::ModeSolved
        Real(8),Intent(out)::UMMode(2)
        Real(8)::R1Loc,R2Loc,R3Loc,R1PLoc,R2PLoc,R3PLoc,H1Loc,H2Loc,H1TLoc,H2TLoc,WLoc,WTrial
        Real(8)::Omiga1Loc,Omiga2Loc,U1Loc,U2Loc,Omiga1TLoc,Omiga2TLoc,U1TLoc,U2TLoc
        Real(8)::DH1Loc,DH2Loc,DWLoc,Delta1Loc,Delta2Loc,Delta3Loc,LambdaLoc,ResNormLoc,TrialNormLoc, &
            ResNormInitLoc,SeedNormLoc
        Real(8)::J11Loc,J12Loc,J13Loc,J21Loc,J22Loc,J23Loc,J31Loc,J32Loc,J33Loc
        Real(8)::ScanH1MinLoc,ScanH1MaxLoc,ScanH2MinLoc,ScanH2MaxLoc,ScanWMinLoc,ScanWMaxLoc
        Real(8)::ScanBestNormLoc,ScanBestH1Loc,ScanBestH2Loc,ScanBestWLoc,ScanBestR1Loc,ScanBestR2Loc, &
            ScanBestR3Loc
        Integer::IterLoc,LSIterLoc,IterLastLoc,LSIterLastLoc,ScanValidCountLoc
        Logical::StateOKLoc,AcceptStepLoc,SolveOKLoc,SeedStateOKLoc
        Character(32)::NegativeExitReasonLoc

        call InitializeNegativeStrategy2()
        ModeSolved=.false.
        UMMode(1)=OmigaL
        UMMode(2)=Qm_L
        w=Wpre
        WState=Wpre
        Delta1Loc=0.0d0
        Delta2Loc=0.0d0
        Delta3Loc=0.0d0
        LambdaLoc=0.0d0
        TrialNormLoc=-1.0d0
        ResNormLoc=-1.0d0
        ResNormInitLoc=-1.0d0
        IterLastLoc=0
        LSIterLastLoc=0
        NegativeExitReasonLoc='init1'
        H1Loc=max(HL,y_ref+Tol1)
        H2Loc=min(max(HR,Tol1),y_ref-Tol1)
        if (H2Loc<=Tol1) H2Loc=max(10.0d0*Tol1,0.5d0*y_ref)
        WLoc=Wpre

        call EvaluateNegativeStateS2(H1Loc,H2Loc,WLoc,R1Loc,R2Loc,R3Loc,Omiga1Loc,U1Loc,Omiga2Loc,U2Loc,StateOKLoc)
        if (.not.StateOKLoc) then
            NegativeExitReasonLoc='init1_invalid'
            H1Loc=y_ref+max(Tol1,0.1d0*D)
            H2Loc=max(10.0d0*Tol1,min(0.5d0*y_ref,y_ref-Tol1))
            WLoc=Wpre
            call EvaluateNegativeStateS2(H1Loc,H2Loc,WLoc,R1Loc,R2Loc,R3Loc,Omiga1Loc,U1Loc,Omiga2Loc,U2Loc,StateOKLoc)
            if (.not.StateOKLoc) NegativeExitReasonLoc='init2_invalid'
        endif
        if (StateOKLoc) ResNormInitLoc=max(abs(R1Loc)/Scale1,abs(R2Loc)/Scale2,abs(R3Loc)/Scale3)
        if (HasS2Seed) then
            call EvaluateNegativeStateS2(S2SeedH1,S2SeedH2,S2SeedW,R1PLoc,R2PLoc,R3PLoc,Omiga1TLoc,U1TLoc,Omiga2TLoc, &
                U2TLoc,SeedStateOKLoc)
            if (SeedStateOKLoc) then
                SeedNormLoc=max(abs(R1PLoc)/Scale1,abs(R2PLoc)/Scale2,abs(R3PLoc)/Scale3)
                if ((.not.StateOKLoc) .or. ResNormInitLoc<0.0d0 .or. SeedNormLoc<ResNormInitLoc) then
                    H1Loc=S2SeedH1
                    H2Loc=S2SeedH2
                    WLoc=S2SeedW
                    R1Loc=R1PLoc
                    R2Loc=R2PLoc
                    R3Loc=R3PLoc
                    Omiga1Loc=Omiga1TLoc
                    Omiga2Loc=Omiga2TLoc
                    U1Loc=U1TLoc
                    U2Loc=U2TLoc
                    StateOKLoc=.true.
                    ResNormInitLoc=SeedNormLoc
                endif
            endif
        endif

        if (StateOKLoc) then
            NegativeExitReasonLoc='iterating'
            do IterLoc=1,60
                IterLastLoc=IterLoc
                ResNormLoc=max(abs(R1Loc)/Scale1,abs(R2Loc)/Scale2,abs(R3Loc)/Scale3)
                if (ResNormLoc<1.0d-8) then
                    ModeSolved=.true.
                    NegativeExitReasonLoc='conv_1e8'
                    exit
                endif

                DH1Loc=max(1.0d-6,1.0d-6*abs(H1Loc))
                DH2Loc=max(1.0d-6,1.0d-6*abs(H2Loc))
                if (H2Loc+DH2Loc>=y_ref-Tol1) DH2Loc=-DH2Loc
                if (H2Loc+DH2Loc<=Tol1) DH2Loc=0.5d0*H2Loc
                DWLoc=max(1.0d-6,1.0d-6*max(abs(WLoc),abs(VL),abs(VR),1.0d0))

                call EvaluateNegativeStateS2(H1Loc+DH1Loc,H2Loc,WLoc,R1PLoc,R2PLoc,R3PLoc,Omiga1TLoc,U1TLoc,Omiga2TLoc,U2TLoc,StateOKLoc)
                if (.not.StateOKLoc) then
                    NegativeExitReasonLoc='jac_h1_invalid'
                    exit
                endif
                J11Loc=(R1PLoc-R1Loc)/DH1Loc
                J21Loc=(R2PLoc-R2Loc)/DH1Loc
                J31Loc=(R3PLoc-R3Loc)/DH1Loc

                call EvaluateNegativeStateS2(H1Loc,H2Loc+DH2Loc,WLoc,R1PLoc,R2PLoc,R3PLoc,Omiga1TLoc,U1TLoc,Omiga2TLoc,U2TLoc,StateOKLoc)
                if (.not.StateOKLoc) then
                    NegativeExitReasonLoc='jac_h2_invalid'
                    exit
                endif
                J12Loc=(R1PLoc-R1Loc)/DH2Loc
                J22Loc=(R2PLoc-R2Loc)/DH2Loc
                J32Loc=(R3PLoc-R3Loc)/DH2Loc

                call EvaluateNegativeStateS2(H1Loc,H2Loc,WLoc+DWLoc,R1PLoc,R2PLoc,R3PLoc,Omiga1TLoc,U1TLoc,Omiga2TLoc,U2TLoc,StateOKLoc)
                if (.not.StateOKLoc) then
                    NegativeExitReasonLoc='jac_w_invalid'
                    exit
                endif
                J13Loc=(R1PLoc-R1Loc)/DWLoc
                J23Loc=(R2PLoc-R2Loc)/DWLoc
                J33Loc=(R3PLoc-R3Loc)/DWLoc

                call SolveLinearSystem3x3(J11Loc,J12Loc,J13Loc,J21Loc,J22Loc,J23Loc,J31Loc,J32Loc,J33Loc, &
                    -R1Loc,-R2Loc,-R3Loc,Delta1Loc,Delta2Loc,Delta3Loc,SolveOKLoc)
                if (.not.SolveOKLoc) then
                    NegativeExitReasonLoc='det_small'
                    exit
                endif

                LambdaLoc=1.0d0
                AcceptStepLoc=.false.
                LSIterLastLoc=0
                do LSIterLoc=1,20
                    LSIterLastLoc=LSIterLoc
                    H1TLoc=H1Loc+LambdaLoc*Delta1Loc
                    H2TLoc=H2Loc+LambdaLoc*Delta2Loc
                    WTrial=WLoc+LambdaLoc*Delta3Loc
                    if (H1TLoc<=y_ref+Tol1 .or. H2TLoc<=Tol1 .or. H2TLoc>=y_ref-Tol1) then
                        LambdaLoc=0.5d0*LambdaLoc
                        cycle
                    endif

                    call EvaluateNegativeStateS2(H1TLoc,H2TLoc,WTrial,R1PLoc,R2PLoc,R3PLoc,Omiga1TLoc,U1TLoc, &
                        Omiga2TLoc,U2TLoc,StateOKLoc)
                    if (.not.StateOKLoc) then
                        LambdaLoc=0.5d0*LambdaLoc
                        cycle
                    endif

                    TrialNormLoc=max(abs(R1PLoc)/Scale1,abs(R2PLoc)/Scale2,abs(R3PLoc)/Scale3)
                    if (TrialNormLoc<ResNormLoc) then
                        H1Loc=H1TLoc
                        H2Loc=H2TLoc
                        WLoc=WTrial
                        R1Loc=R1PLoc
                        R2Loc=R2PLoc
                        R3Loc=R3PLoc
                        Omiga1Loc=Omiga1TLoc
                        Omiga2Loc=Omiga2TLoc
                        U1Loc=U1TLoc
                        U2Loc=U2TLoc
                        AcceptStepLoc=.true.
                        exit
                    endif
                    LambdaLoc=0.5d0*LambdaLoc
                enddo

                if (.not.AcceptStepLoc) then
                    NegativeExitReasonLoc='line_search_fail'
                    exit
                endif
            enddo

            ResNormLoc=max(abs(R1Loc)/Scale1,abs(R2Loc)/Scale2,abs(R3Loc)/Scale3)
            if (ResNormLoc<1.0d-5) then
                ModeSolved=.true.
                if (NegativeExitReasonLoc/='conv_1e8') NegativeExitReasonLoc='conv_1e5'
            elseif (NegativeExitReasonLoc=='iterating') then
                NegativeExitReasonLoc='iter_limit'
            endif
        endif

        if (ModeSolved) then
            w=WLoc
            WState=WLoc
            int_m=(/Omiga1Loc,Omiga1Loc*U1Loc,Omiga2Loc,Omiga2Loc*U2Loc/)
            call SelectNegativeGodunovState(Omiga2Loc,U2Loc,UMMode,StateOKLoc)
            if (.not.StateOKLoc) then
                UMMode(1)=Omiga2Loc
                UMMode(2)=Omiga2Loc*U2Loc
            endif
        else
            call ScanNegativeResidualsS2(ScanH1MinLoc,ScanH1MaxLoc,ScanH2MinLoc,ScanH2MaxLoc,ScanWMinLoc,ScanWMaxLoc, &
                ScanBestNormLoc,ScanBestH1Loc,ScanBestH2Loc,ScanBestWLoc,ScanBestR1Loc,ScanBestR2Loc,ScanBestR3Loc,ScanValidCountLoc)
            Write(99,*) trim(FailTag),NegativeExitReasonLoc,IterLastLoc,LSIterLastLoc,ResNormInitLoc,ResNormLoc, &
                R1Loc,R2Loc,R3Loc,H1Loc,H2Loc,WLoc,Delta1Loc,Delta2Loc,Delta3Loc,LambdaLoc,TrialNormLoc,Wpre,HL,HR,VL,VR
            Write(99,*) trim(ScanTag),ScanValidCountLoc,ScanH1MinLoc,ScanH1MaxLoc,ScanH2MinLoc,ScanH2MaxLoc, &
                ScanWMinLoc,ScanWMaxLoc,ScanBestNormLoc,ScanBestH1Loc,ScanBestH2Loc,ScanBestWLoc,ScanBestR1Loc,ScanBestR2Loc,ScanBestR3Loc
        endif
    end subroutine

End function
