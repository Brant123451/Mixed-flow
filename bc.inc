If(Left_boundary==2) then !'连续边界'
	H(0)=H(1)
	A(0)=A(1)
	Q(0)=Q(1)
	C(0)=C(1)
	V(0)=V(1)
	Dens(0)=Dens(1)
Else
	If(Left_boundary==1) then !封闭边界
		H(0)=H(1)
		A(0)=Area(H(0))
		Q(0)=-Q(1)
		V(0)=Q(0)/A(0)
		C(0)=Celerity(H(0))
		Dens(0)=Dens(1)
	Else
		If(Left_boundary==3) then !'定水深'
			H(0)=LeftPressHead
			A(0)=Area(H(0))
			Q(0)=Q(1)
			C(0)=Celerity(H(0))
			V(0)=Q(0)/A(0)
			Dens(0)=Density(H(0))
		End if
	Endif
End if

If(Right_boundary==1) then !'右侧盲端'
	H(nc)=H(n)
	A(nc)=A(n)
	Q(nc)=-Q(n)
	C(nc)=C(n)
	V(nc)=Q(nc)/A(nc)
	Dens(nc)=Dens(n)

Else 
	If(Right_boundary==2) then !'连续边界'
		H(nc)=H(n)
		A(nc)=A(n)
		Q(nc)=Q(n)
		C(nc)=C(n)
		V(nc)=V(n)
		Dens(nc)=Dens(n)
	Else
		If(Right_boundary==3) then !'定水深'
			H(nc)=11
			A(nc)=Area(H(nc))
			Q(nc)=Q(n)
			C(nc)=Celerity(H(nc))
			V(nc)=Q(nc)/A(nc)
			Dens(nc)=Density(H(nc))
		End if
	End if
	
End if