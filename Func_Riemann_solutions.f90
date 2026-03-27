Function Riemann_solutions(U1L,U1R,U2L,U2R) Result(UM)
    use variables
	Implicit none
	Real(8)::UM(2)
	Real(8)::F,AM,AL,AR,VR,VL,F_Diff,AM2,HL,HR,CL,CR,QM,HM
	Real(8)::Phi_Circular,U1L,U1R,U2L,U2R
	Real(8)::F_Equation,Depth,Celerity,Area,Density
    Real(8)::m,n,mn,Fm,Fn,Fmn,OmigaM,QmM,Omiga_m,Omiga_n,Omiga_mn
    Integer::k
    
    VL=U2L/U1L
    VR=U2R/U1R
    HL=Depth(U1L)
    HR=Depth(U1R)
    CL=Celerity(HL)
    CR=Celerity(HR)
    k=1
    OmigaM=0.5*(U1L+U1R)+(U1L+U1R)/(2*(CL+CR))*(VL-VR)
    
    m=0.00001 !最低与最高水头
    n=y_ref
    
    Do
        mn=0.5*(m+n)
        Omiga_m=Area(m)*Density(m)
        Omiga_n=Area(n)*Density(n)
        Omiga_mn=Area(mn)*Density(mn)
        Fm=F_Equation(Omiga_m,U1L)+F_Equation(Omiga_m,U1R)+VR-VL
        Fn=F_Equation(Omiga_n,U1L)+F_Equation(Omiga_n,U1R)+VR-VL
        Fmn=F_Equation(Omiga_mn,U1L)+F_Equation(Omiga_mn,U1R)+VR-VL
        If(Fm*Fmn<0) Then
            n=mn
        Else
            m=mn        
        Endif
        
        If(abs(m-n)<0.00001) Exit
        If(Fm*Fn>0) Then
            Print*,"Riemann_solutions_错误_2"
            !Stop
        End if
    Enddo
    If(Fm*Fn>0) Then
        mn=y_ref
    Endif
    OmigaM=Area(mn)*Density(mn)
    QmM=0.5*OmigaM*(VR+VL)+0.5*OmigaM*(F_Equation(OmigaM,U1R)-F_Equation(OmigaM,U1L))
    UM(1)=OmigaM
    UM(2)=QmM
    
End Function