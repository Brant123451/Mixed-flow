Function APY(Omiga) Result(New_Result)
    use variables
    Implicit none
    
    Real(8)::H,New_Result,Area,Omiga,Depth
    Real(8)::kk
    
    kk=1.0d0
    H=Depth(Omiga)
    If(Omiga<Omiga_ref) Then !
        New_Result=Dens_ref*g/12*((3*D**2-4*D*H+4*H**2)*sqrt(H*(D-H))-3*D**2*(D-2*H)*atan(sqrt(H/(D-H)))) !+kk*A_ref*(Pg-P_ref)
    Else    !
        New_Result=(Dens_ref*g/12*((3*D**2-4*D*y_ref+4*y_ref**2)*sqrt(y_ref*(D-y_ref))-3*D**2*(D-2*y_ref)*atan(sqrt(y_ref/(D-y_ref))))/A_ref+Dens_ref*g*(H-y_ref))*A_ref
    Endif
End function APY