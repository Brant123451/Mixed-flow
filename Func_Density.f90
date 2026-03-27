Function Density(H) Result(NEW_Result)
    use variables
	Implicit None
    
	Real(8)::A,Theta,NEW_Result,H 
    If(H<y_ref) then !自由液面流
        NEW_Result=Dens_ref
    Else !压力流
        New_Result=Dens_ref+Dens_ref*g*(H-y_ref)/C_Waterhammer**2
    Endif
End Function