Function Choice(RanSlope,SL,SR,FL,FR,FLL) result(Choice_Result)
    
	Implicit none
	Real(8),Dimension(2)::FL,FR,FLL,Choice_Result
	Real(8)::RanSlope,SL,SR

	If(RanSlope<SL) Then	
		Choice_RESULT=FL
	Else	
		If(RanSlope>=SL.AND.RanSlope<SR) Then
			Choice_RESULT=FLL
		Else
			Choice_RESULT=FR
		End if
	End if

End function