if(Time>(Time_step*step)) then
    
    write(2,'(A7,I6,A4,I6,A30)') str3,4*n+2,str4,2*n,str5
    write(2,*) str6
    write(2,*) str7
    
!横坐标&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&            
    write(2,*) Orix

    if(H(1)<y_ref) then
        write(2,*) Orix-H(1)*sin(Theta)
    else
        write(2,*) Orix-y_ref*sin(Theta)
    end if

    write(2,*) Orix-y_ref*sin(Theta)
            
    do i=2,n,1
            
        write(2,*) Orix+Deltax*(i-1)*cos(Theta)

        if(H(i-1)<y_ref) then
            write(2,*) Orix+Deltax*(i-1)*cos(Theta)-H(i-1)*sin(Theta)
        else
            write(2,*) Orix+Deltax*(i-1)*cos(Theta)-y_ref*sin(Theta)
        end if


        if(H(i)<y_ref) then
            write(2,*) Orix+Deltax*(i-1)*cos(Theta)-H(i)*sin(Theta)
        else
            write(2,*) Orix+Deltax*(i-1)*cos(Theta)-y_ref*sin(Theta)
        end if

        write(2,*) Orix+Deltax*(i-1)*cos(Theta)-y_ref*sin(Theta)
                
    end do

    write(2,*) Orix+L*cos(Theta)

    if(H(n)<y_ref) then
        write(2,*) Orix+L*cos(Theta)-H(n)*sin(Theta)
    else
        write(2,*) Orix+L*cos(Theta)-y_ref*sin(Theta)
    end if

    write(2,*) Orix+L*cos(Theta)-y_ref*sin(Theta)
!横坐标&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!纵坐标&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    write(2,*) Oriy
    
    if(H(1)<y_ref) then
        write(2,*) Oriy+H(1)*cos(Theta)
    else
        write(2,*) Oriy+y_ref*cos(Theta)
    end if

    write(2,*) Oriy+y_ref*cos(Theta)

    do i=2,N,1

        write(2,*) Oriy+Deltax*(i-1)*sin(Theta)

        if(H(i-1)<y_ref) then
            write(2,*) Oriy+Deltax*(i-1)*sin(Theta)+H(i-1)*cos(Theta)
        else
            write(2,*) Oriy+Deltax*(i-1)*sin(Theta)+y_ref*cos(Theta)
        end if

        if(H(i)<y_ref) then
            write(2,*) Oriy+Deltax*(i-1)*sin(Theta)+H(i)*cos(Theta)
        else
            write(2,*) Oriy+Deltax*(i-1)*sin(Theta)+y_ref*cos(Theta)
        end if

        write(2,*) Oriy+Deltax*(i-1)*sin(Theta)+y_ref*cos(Theta)

    end do

    write(2,*) Oriy+L*sin(Theta)
    
    if(H(n)<y_ref) then
        write(2,*) Oriy+L*sin(Theta)+H(n)*cos(Theta)
    else
        write(2,*) Oriy+L*sin(Theta)+y_ref*cos(Theta)
    end if

    write(2,*) Oriy+L*sin(Theta)+y_ref*cos(Theta)
!纵坐标&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!压力&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
    do i=1,n,1
        
        write(2,*) 1 !水

        if(H(i)>=y_ref) then

            write(2,*) 1 !有压流——气

        else
        
            write(2,*) 0 !明渠流——气

        end if

    end do
!压力&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!网格结点&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    write(2,'(4I4)') 1,2,5,4
    write(2,'(4I4)') 2,3,7,5

    do i=2,n-1,1

        write(2,'(4I6)') (i-1)*4,(i-1)*4+2,(i-1)*4+5,(i-1)*4+4
        write(2,'(4I6)') (i-1)*4+2,(i-1)*4+3,(i-1)*4+7,(i-1)*4+5

    end do

    write(2,'(4I6)') 4*(n-1),4*(n-1)+2,4*n+1,4*n
    write(2,'(4I6)') 4*(n-1)+2,4*(n-1)+3,4*n+2,4*n+1

    Step=Step+1
    Print*,time
end if
         
If(Time>Time_putout) exit