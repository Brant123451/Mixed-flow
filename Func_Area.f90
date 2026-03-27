Function Area(H) Result(NEW_Result)
    use variables    
	Implicit none
    
	Real(Kind=8)::H,A,Theta,NEW_Result
    
    If(H<=y_ref) Then     !Ã÷Á÷
        Theta=2*ACOS(1-2*H/D)
        A=D**2/8*(Theta-SIN(Theta))
    Else
        Theta=2*ACOS(1-2*y_ref/D)
        A=D**2/8*(Theta-SIN(Theta))
    End if
	NEW_Result=A
End Function    