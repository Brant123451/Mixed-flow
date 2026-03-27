Function Phi_Circular(Omiga) Result(Phi_Result)
    use variables     
    Implicit none
    Real(8)::Phi,K,Area,Depth
    Real(8)::H,Theta,Phi_Result,A,Omiga
    
    H=Depth(Omiga)
    K=sqrt(3.0)
    If(H<=y_ref) Then
        Theta=2*ACOS((D-2*H)/D)
        PHI_Result=sqrt(g*D/8)*(K*Theta-K/80*Theta**3+19*K/448000*Theta**5+K/10035200*Theta**7+491.0d0*K/(27.0d0*7064780800.0d0)*Theta**9)
    Else
        PHI_Result=C_Waterhammer*log(Omiga)
    Endif

End function