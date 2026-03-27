Function F_Equation(OmigaM,OmigaK) Result(F_Result) 
    use variables     
    Implicit none
    
    Real(8)::Phi_Circular,OmigaM,OmigaK
    Real(8)::F_Result
      
    F_Result=Phi_Circular(OmigaM)-Phi_Circular(OmigaK)
End function