Function Depth(U1) Result(NEW_Result)
    use variables
	Implicit None
    
	Real(8)::U1,Theta,NEW_Result,H,U1_ref,A
	Real(8)::Dichotomization
    
    If(U1<Omiga_ref) then !
        A=U1/Dens_ref
        Theta=Dichotomization(A)
        H=0.5*(D-D*COS(0.5*Theta))
    Else
        H=((U1/A_ref-Dens_ref)*C_Waterhammer**2)/(Dens_ref*g)+y_ref
    Endif
    New_Result=H
    
End Function