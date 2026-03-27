Function SOLU(APY_FreeSurface) Result(New_Result)
    use variables
    Implicit none
    
    Real(8)::H,New_Result,Area,APY_FreeSurface,Depth
    
    New_Result=(APY_FreeSurface-Dens_ref*g/12*((3*D**2-4*D*y_ref+4*y_ref**2)*sqrt(y_ref*(D-y_ref))-3*D**2*(D-2*y_ref)*atan(sqrt(y_ref/(D-y_ref)))))/(Dens_ref*g*A_ref)+y_ref

End function SOLU
