Function Celerity(H) Result(Cele_Result)
	use variables
	Implicit none
    
	Real(8)::H,Cele_Result,Theta,WID
	Real(8)::Area

    If(H<y_ref) Then !Ã÷Çþ
        Theta=2*acos((D-2*H)/D)
        WID=D*sin(Theta/2)
        Cele_Result=sqrt(g*Area(H)/WID) 
    Else  
        Cele_Result=C_Waterhammer     !ÓÐÑ¹Á÷
    End if
End function