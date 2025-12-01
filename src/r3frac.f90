        pure subroutine r3frac(eps,v)
        implicit none
        ! arguments
        real(8), intent(in) :: eps
        real(8), intent(inout) :: v(3)
        ! local variables
        real(8) t1
        
        t1=1.d0-eps
        v(1)=v(1)-int(v(1))
        if (v(1) < 0.d0) v(1)=v(1)+1.d0
        if ((v(1) < eps).or.(v(1) > t1)) v(1)=0.d0
        v(2)=v(2)-int(v(2))
        if (v(2) < 0.d0) v(2)=v(2)+1.d0
        if ((v(2) < eps).or.(v(2) > t1)) v(2)=0.d0
        v(3)=v(3)-int(v(3))
        if (v(3) < 0.d0) v(3)=v(3)+1.d0
        if ((v(3) < eps).or.(v(3) > t1)) v(3)=0.d0
        
        end subroutine
