Do i=0,int(7.9744/deltax),1
	H(i)=0.06
	A(i)=Area(H(i))
	Dens(i)=Density(H(i))
	V(i)=0.18
	Q(i)=A(i)*V(i)
	C(i)=Celerity(H(i))
Enddo
Do i=int(7.9744/deltax)+1,nc,1
	H(i)=0.020
	A(i)=Area(H(i))
	Dens(i)=Density(H(i))
	V(i)=0.40
	Q(i)=A(i)*V(i)
	C(i)=Celerity(H(i))
Enddo