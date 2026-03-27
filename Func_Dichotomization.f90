Function Dichotomization(A) Result(Theta_Result)    !FIND THETA IN A
    use variables
    Implicit none
    
    Real(Kind=8)::Fa,Fb,Fm,A,m,n,Theta_Result

    Integer::k
    
    M=0
    N=2*Pai-2*acos((2*y_ref-D)/D)
    
    If(A<0) Then
        Print*,"Function Dichotomization 输入负值A",A
        stop
    Endif
    If(A<A_ref) then
    Do      
        Fa=M-sin(M)-8*A/D**2
        Fb=N-sin(N)-8*A/D**2
        Fm=(M+N)/2-sin((M+N)/2)-8*A/D**2
        If(abs(Fm)<0.00000001) Exit         
        If(Fa*Fm<0) Then
            N=(M+N)/2            
        Else        
            M=(M+N)/2           
        Endif
        If(k==1000) then
            Print*,"Function Dichotomization 循环错误"
            Print*,A
            !stop
        End if 
    Enddo
    Endif
            
    Theta_Result=(M+N)/2
    
End function