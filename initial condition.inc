Do i=0,InitIfaceCell,1
	H(i)=LeftPressHead
	A(i)=Area(H(i))
	Dens(i)=Density(H(i))

	V(i)=0.0d0
	Q(i)=A(i)*V(i)
	C(i)=Celerity(H(i))
Enddo
Do i=InitIfaceCell+1,nc,1
	H(i)=max(flow_depth,Tol1)
	A(i)=Area(H(i))
	Dens(i)=Density(H(i))
	V(i)=0.0d0
	Q(i)=A(i)*V(i)

	C(i)=Celerity(H(i))
Enddo