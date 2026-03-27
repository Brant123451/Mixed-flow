Function Hydraulic(A) Result(NEW_Hydraulic)
	use variables
	Implicit none
    
	Real(8)::A,H,NEW_Hydraulic,Rhyd,Theta
	Real(8)::Dichotomization,Depth
	H=Depth(A)
    
    !If(A<0) Print*,"Function Hydraulic ÊäÈë¸ºÖµA",A
    If(H<y_ref) then
        Theta=Dichotomization(A)
        Rhyd=D/4*(1-sin(Theta)/Theta)
    Else  
        Rhyd=0.25*D      
    End if

	NEW_Hydraulic=Rhyd
End function