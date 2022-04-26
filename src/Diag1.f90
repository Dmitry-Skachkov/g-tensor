


      Program Diag1
       real(8), parameter :: PI=3.1415926535897932384626433d0 
       real(8)    :: A(3,3)            ! matrix
       real(8)    :: Q(3,3)            ! eigenvectors
       real(8)    :: W(3)              ! eigenvalues
       real(8)    :: Az(3)             ! z axis for angle calulations
       real(8)    :: Ay(3)             ! y axis
       real(8)    :: Ax(3)             ! x axis 
       real(8)    :: Ab1(3)            ! a,b,c axis from input             
       real(8)    :: Ac1(3)             
       real(8)    :: Aa1(3)             
       real(8)    :: Ab2(3)             
       real(8)    :: Ac2(3)             
       real(8)    :: Aa2(3)             
       real(8)    :: V(3,3)            ! g-tensor vectors for plot
       real(8)    :: theta,fi          ! theta vs to z, fi vs to x
       real(8)    :: theta1(3),fi1(3)  ! store theta and fi for principal axis
       real(8)    :: u1(3),v1(3),w1(3)
       real(8)    :: ux,vx,wx
       integer    :: Iz(3),Iy(3)       ! how we define z and y axis from the actual axis 
       character(1):: do
       call getarg(1,do)               ! g or A
       open(unit=2,file='g-tensor.dat')
       read(2,*) A                     ! tensor in Aa1,Ab1,Ac1 coordinates
       read(2,*)
       read(2,*) Aa1                   ! read a vector from QE input 
       read(2,*) Ab1                   !      b vector
       read(2,*) Ac1                   !      c vector
       read(2,*)
       read(2,*) Iz                    !      real crystolagrafic b vector (z) for Ga2O3, for LiGaO2 it's z (c) vector!
       read(2,*) Iy                    !      real crystolagrafic c vector (-y) for Ga2O3, for LiGaO2 it's y (b) vector
       call print_lattice(Aa1,Ab1,Ac1)
       Az(1:3) = Iz(1)*Aa1(1:3) + Iz(2)*Ab1(1:3) + Iz(3)*Ac1(1:3)      !define real crystolagrafic b(z) axis for Ga2O3, this is c(z) axis for LiGaO2
       Ay(1:3) = Iy(1)*Aa1(1:3) + Iy(2)*Ab1(1:3) + Iy(3)*Ac1(1:3)      !define real crystolagrafic c(y) axis for Ga2O3, this is b(y) axis for LiGaO2
       print *,'defined axis'
       print 111,Az
       print 112,Ay
111    format(' z=',3F10.5,' for Ga2O3 it s b axis, for LiGaO2 it s c axis')
112    format(' y=',3F10.5,' for Ga2O3 it s c axis, for LiGaO2 it s b axis')
       call calc_a(Az,Ay,Ax)           ! a* for Ga2O3, and a axis for LiGaO2

       call make_norm(Ax)              ! crystal a axis for Ga2O3 (x)
       call make_norm(Az)              ! crystal b axis (z) for Ga2O3, for LiGaO2 z axis (c)
       call make_norm(Ay)              ! crystal c axis (-y) for Ga2O3, for LiGaO2 y axis (b)
       call print_lattice2(Ax,Ay,Az)
!       print *,'check orthogonality a,b,c:'
!       call test_orthog(Ac,Ab,Aas)
!       print *,'abs(Aas)=',dsqrt(Aas(1)**2+Aas(2)**2+Aas(3)**2)

       call symm1(A)                   ! symmetrize A 
       call DSYEVJ3(A, Q, W)           ! diagonalize matrix A and calculate eigenvectors (Q) and eigenvalues (W)
       print 5,W
       print 6
!       call make_norm(Q(1:3,1))
!       call make_norm(Q(1:3,2))
!       call make_norm(Q(1:3,3))
       do j=1,3
        print 61,Q(1:3,j),dsqrt(Q(1,j)**2+Q(2,j)**2+Q(3,j)**2)
       enddo
       
       print *
       print *,'test orthogonality for eigenvectors'
       call test_orthog(Q(1:3,1),Q(1:3,2),Q(1:3,3))
       print *

!      calculate theta and fi for all principal axis
!       print 1, Aas
!       print 9, Ab
!       print 10,Ac


       Aa2 = Aa1            ! calculate Miller indexes for this axis (untransformed)
       Ab2 = Ab1
       Ac2 = Ac1

!       Aa2 = (/Aa1(1),-Aa1(2),Aa(3)/) ! calculate Miller indexes for this axis (transformed)
!       Ab2 =  Ac1
!       Ac2 = -Ab1

       print *
       do j=1,3
        print *
        print *
        print *
        print 17
        print 18,j,W(j)
        call calc_theta_fi(Q(1:3,j),Az,Ay,Ax,theta,fi)
        theta1(j) = theta
        fi1(j) = fi
        call make_range_angle(theta,fi)                  !!!   change angles !!!! -> reduce accuracy!!!!!!!! 
        print 8,j,Q(1:3,j)
        print 81,theta1(j)*180.d0/PI,theta*180.d0/PI,fi1(j)*180.d0/PI,fi*180.d0/PI 
        print 11,Q(1:3,j)*W(j)/10000.d0 
        V(1:3,j) =  Q(1:3,j)*W(j)/10000.d0
        call miller_indexes(W(j),Q(1:3,j),Aa2,Ab2,Ac2,ux,vx,wx,do)              ! calculate Miller indexes u,v,w in the initial coordinate system for VESTA
        u1(j) = ux
        v1(j) = vx
        w1(j) = wx
       enddo
       print 17
       print *
       print *
       print *
       call check_theta_fi(theta1,fi1,Q,Ax,Ay,Az)
       call check_uvw(u1,v1,w1,Aa2,Ab2,Ac2)
!       call check_directions(Aa2,Ab2,Ac2)
   5   format(/' Tensor principal values=',3F12.4)
   6   format(/' axis for principal values=')
  61   format(3F14.7,'  abs=',F25.17)
   8   format(/'axis ',I2,3F14.7)
  81   format(' theta=',F6.1,' ( ',F6.1,' )    fi=',F6.1,' ( ',F6.1,' ) ' ) !,' th_a=',F6.1,' th_c=',F6.1/)
  11   format(/' axis for plot (*g/10000) in original system x,y,z'/3F14.4) 
  17   format(80('*'))
  18   format(/' g(',I1,')=',F12.4)
      end program Diag1






!      subroutine check_directions(Aa,Ab,Ac)
!       real(8), dimension(3) :: Aa,Ab,Ac
!       real(8), dimension(3) :: V,R,VX
!       print *
!       print *,' check (1,0,0)'
!       call print_VR(1.d0,0.d0,0.d0,Aa,Aa,Ab,Ac,'a')
!       call print_VR(1.d0,0.d0,0.d0,Ab,Aa,Ab,Ac,'b')
!       call print_VR(1.d0,0.d0,0.d0,Ac,Aa,Ab,Ac,'c')
!       print *
!       print *,' check (0,1,0)'
!       call print_VR(0.d0,1.d0,0.d0,Aa,Aa,Ab,Ac,'a')
!       call print_VR(0.d0,1.d0,0.d0,Ab,Aa,Ab,Ac,'b')
!       call print_VR(0.d0,1.d0,0.d0,Ac,Aa,Ab,Ac,'c')
!       print *
!       print *,' check (0,0,1)'
!       call print_VR(0.d0,0.d0,1.d0,Aa,Aa,Ab,Ac,'a')
!       call print_VR(0.d0,0.d0,1.d0,Ab,Aa,Ab,Ac,'b')
!       call print_VR(0.d0,0.d0,1.d0,Ac,Aa,Ab,Ac,'c')
!      end subroutine check_directions



!      subroutine print_VR(u,v,w,R,Aa,Ab,Ac,Comm)
!       real(8), dimension(3) :: Aa,Ab,Ac
!       real(8), dimension(3) :: R
!       real(8)               :: u,v,w
!       character(1)          :: Comm
!       real(8), dimension(3) :: R1
!       real(8), dimension(3) :: VX
!       call calc_vector_from_uvw(u,v,w,Aa,Ab,Ac,VX)
!       call make_norm(VX)
!       R1 = R
!       call make_norm(R1)
!       print 1,u,v,w,Comm,(VX(1)*R1(1)+VX(2)*R1(2)+VX(3)*R1(3))
!1      format(' (',F2.0,',',F2.0,',',F2.0,')*',A1,'=',F15.7)
!      end subroutine print_VR





!      subroutine calc_vector_from_uvw(u,v,w,Aa,Ab,Ac,R)
!       real(8)   :: R(3)
!       real(8)   :: Aa(3),Ab(3),Ac(3)
!       real(8)   :: u,v,w
!       R(1:3) = u*Aa(1:3) + v*Ab(1:3) + w*Ac(1:3)
!      end subroutine calc_vector_from_uvw







      subroutine print_lattice2(Aas,Ac,Ab)
       real(8)    :: Aas(3),Ac(3),Ab(3)
       print 1
       print 2,Aas             
       print 3,Ac             
       print 4,Ab
1      format(/' Actual crystolagrafic axis for Ga2O3')
2      format(' x (a*) =',3F10.5)             
3      format(' y (c)  =',3F10.5)             
4      format(' z (b)  =',3F10.5)             
      end subroutine print_lattice2








      subroutine calc_a(Ab,Ac,Aa)                        ! a = [b*c]
       real(8)        :: Ab(3),Ac(3),Aa(3)
       real(8)        :: Abm,Acm
       Abm = dsqrt(Ab(1)**2+Ab(2)**2+Ab(3)**2)
       Acm = dsqrt(Ac(1)**2+Ac(2)**2+Ac(3)**2)
       Aa(1) = (Ab(2)*Ac(3) - Ab(3)*Ac(2))/(Abm*Acm)
       Aa(2) = (Ab(3)*Ac(1) - Ab(1)*Ac(3))/(Abm*Acm)
       Aa(3) = (Ab(1)*Ac(2) - Ab(2)*Ac(1))/(Abm*Acm)
      end subroutine calc_a




      subroutine symm1(A)            ! make A symmetrized
       real(8)    :: A(3,3)        
       a(1,2) = 0.5d0*(a(1,2)+a(2,1))
       a(2,1) = a(1,2)
       a(1,3) = 0.5d0*(a(1,3)+a(3,1))
       a(3,1) = a(1,3)
       a(2,3) = 0.5d0*(a(2,3)+a(3,2))
       a(3,2) = a(2,3)
      end subroutine symm1




      subroutine make_norm(Ab)
       implicit none
       real(8)         :: Ab(3)
       real(8)         :: Ar
       Ar = dsqrt(Ab(1)**2+Ab(2)**2+Ab(3)**2)
       Ab(1:3) = Ab(1:3)/Ar
      end subroutine make_norm





      subroutine calc_theta_fi(A1,Az,Ay,Ax,theta,fi)
       real(8)             :: A1(3),Ax(3),Ay(3),Az(3),A1z(3),A1xy(3),A1xym        
       real(8)             :: theta,fi
       real(8)             :: costheta,cosfi 
       costheta = (A1(1)*Az(1)+A1(2)*Az(2)+A1(3)*Az(3))                              !    abs(A1) = 1
       theta = dacos(costheta)                                                       ! angle between A1 and z       
       A1z(1:3) = costheta*Az(1:3)   
       A1xy(1:3) = A1(1:3) - A1z(1:3)
       A1xym = dsqrt(A1xy(1)**2+A1xy(2)**2+A1xy(3)**2)
       cosfi = (A1xy(1)*Ax(1)+A1xy(2)*Ax(2)+A1xy(3)*Ax(3))/A1xym
       fi    = dacos(cosfi)*dsign(1.d0,(A1xy(1)*Ay(1)+A1xy(2)*Ay(2)+A1xy(3)*Ay(3)))  ! angle between A1 projection on plane (x,y) and x with the sign of rotation around z       
      end subroutine calc_theta_fi






      subroutine print_lattice(Aa,Ab,Ac)
       real(8), parameter :: PI=3.1415926535897932384626433d0 
       real*8             :: Aa(3),Ab(3),Ac(3)
       real*8             :: Aam,Abm,Acm,sbetta
       Aam=dsqrt(Aa(1)**2+Aa(2)**2+Aa(3)**2)
       Abm=dsqrt(Ab(1)**2+Ab(2)**2+Ab(3)**2)
       Acm=dsqrt(Ac(1)**2+Ac(2)**2+Ac(3)**2)
       sbetta = dacos((Aa(1)*Ac(1)+Aa(2)*Ac(2)+Aa(3)*Ac(3))/(Aam*Acm))*180.d0/PI
       print 11
       print 12,Aa
       print 13,Ab
       print 14,Ac
       print 1
       print 2,Aam
       print 3,Abm
       print 4,Acm
       print 5,sbetta
1      format(/' Parameters of the lattice:')
2      format(' a=',F8.3)
3      format(' b=',F8.3)
4      format(' c=',F8.3)
5      format(' betta=',F12.3)
11     format(/' Vectors from QE input:')
12     format(' a=',3F10.5)
13     format(' b=',3F10.5)
14     format(' c=',3F10.5)
      end subroutine print_lattice








       subroutine make_range_angle(theta,fi)
        real(8), parameter :: PI=3.1415926535897932384626433d0 
        real(8)            :: theta,fi
        if(theta > PI*0.5d0) then
         theta = PI - theta
        endif
        if(fi > PI*0.5d0) then
         fi = -(PI - fi)
        elseif(fi < -0.5d0*PI) then
         fi =  (PI + fi)
        endif
       end subroutine make_range_angle







      subroutine test_orthog(Q1,Q2,Q3)
       real(8)  :: Q1(3),Q2(3),Q3(3)
       print 1,Q1(1)*Q2(1)+Q1(2)*Q2(2)+Q1(3)*Q2(3)
       print 2,Q2(1)*Q3(1)+Q2(2)*Q3(2)+Q2(3)*Q3(3)
       print 3,Q3(1)*Q1(1)+Q3(2)*Q1(2)+Q3(3)*Q1(3)
1      format('Q1*Q2=',F25.17)
2      format('Q2*Q3=',F25.17)
3      format('Q3*Q1=',F25.17)
      end subroutine test_orthog   
       
       




      subroutine check_theta_fi(theta,fi,Q,Ax,Ay,Az)
       real(8), parameter :: PI=3.1415926535897932384626433d0 
       real(8)            :: theta(3),fi(3)
       real(8)            :: V(3,3)
       real(8)            :: Q(3,3)
       real(8)            :: Ax(3),Ay(3),Az(3)
       print *
       print *,'test theta and fi:'
       do j= 1,3
        V(1,j) = dsin(theta(j))*dcos(fi(j)) 
        V(2,j) = dsin(theta(j))*dsin(fi(j)) 
        V(3,j) = dcos(theta(j))
       enddo
       print *,'orthonormality between three vectors built from theta and fi'
       call test_orthog(V(1:3,1),V(1:3,2),V(1:3,3))
       do j=1,3                                                                 ! move to (Ax,Ay,Az) coordinate system
        V(1:3,j) = V(1,j)*Ax(1:3) + V(2,j)*Ay(1:3) + V(3,j)*Az(1:3) 
       enddo
       print *,'check (Q*V) of three vectors built from theta and fi and initial eigenvectors'
       call test_parall(V(1:3,1),V(1:3,2),V(1:3,3),Q(1:3,1),Q(1:3,2),Q(1:3,3))
1      format('V=',3F25.17,'  abs(V)=',F25.17)
      end subroutine check_theta_fi







      subroutine test_parall(V1,V2,V3,Q1,Q2,Q3)
       real(8)            :: V1(3),V2(3),V3(3)
       real(8)            :: Q1(3),Q2(3),Q3(3)
       print 1,1,(V1(1)*Q1(1)+V1(2)*Q1(2)+V1(3)*Q1(3))
       print 1,2,(V2(1)*Q2(1)+V2(2)*Q2(2)+V2(3)*Q2(3))
       print 1,3,(V3(1)*Q3(1)+V3(2)*Q3(2)+V3(3)*Q3(3))
1      format(I3,F25.17)
      end subroutine test_parall





      subroutine miller_indexes(We,Qe,Aa,Ab,Ac,u,v,w,do)
       real(8)   :: Qe(3),We
       real(8)   :: Aa(3),Ab(3),Ac(3)
       real(8)   :: u,v,w
       real(8)   :: M(3,3),Y(3)
       character(1)    :: do
       if(do=='A') then
        c = 1.1d0
       else    ! for g-tensor
        c = We/10000.d0                                                        ! coefficient for plotting g-tensor proportional to eigenvalues
       endif
       M(1,1:3) = Aa(1:3)
       M(2,1:3) = Ab(1:3)
       M(3,1:3) = Ac(1:3)
!     test for solve3x3
!       M(1:3,1) = (/2,-1,6/)
!       M(1:3,2) = (/0,5,-3/)
!       M(1:3,3) = (/0,0,1/)
!       Qe(1:3) = (/20,-2,4/)
       call solve3x3(M,Qe,Y)                                                  ! solve 3x3 equation M*Y=Qe
       u = Y(1)
       v = Y(2)
       w = Y(3)
!!       call check_vector(u,v,w,Aa,Ab,Ac,Qe)
       print 1
       print 3
       print 2,c*u,c*v,c*w                                                    ! print g-tensor vectors for both + and - directions
       print 2,-c*u,-c*v,-c*w
       print 3
1      format(/' axis of g-tensor in lattice vector notations u,v,w (for plot in VESTA)')
2      format(3F15.7) 
3      format(80('-'))
      end subroutine miller_indexes

       




      subroutine check_vector(u,v,w,Aa,Ab,Ac,Qe)
       real(8)   :: Qe(3)
       real(8)   :: Aa(3),Ab(3),Ac(3)
       real(8)   :: u,v,w
       real(8)   :: R(3)
       R(1:3) = u*Aa(1:3) + v*Ab(1:3) + w*Ac(1:3)
       print *
       print *,'test u,v,w '
       print 1,Qe
       print 2,R
1      format('Initial vector         ',3F18.10)
2      format('Vector u*a + v*b + w*c ',3F18.10)
      end subroutine check_vector







     subroutine solve3x3(A,B,X)            ! solve 3x3 linear equation AX=B https://en.wikipedia.org/wiki/Cramer%27s_rule
      real(8)  :: A(3,3)
      real(8)  :: B(3)
      real(8)  :: X(3)
      real(8)  :: D,D1,D2,D3
      real(8)  :: Aw(3,3)
      call det3x3(A,D)
      call make_dx(A,1,B,Aw)
      call det3x3(Aw,D1)
      call make_dx(A,2,B,Aw)
      call det3x3(Aw,D2)
      call make_dx(A,3,B,Aw)
      call det3x3(Aw,D3)
      X(1) = D1/D
      X(2) = D2/D
      X(3) = D3/D
!      print *,'solution X'
!      print 1,X
!!      call check_solution(A,B,X)
1     format(3F15.8)
     end subroutine solve3x3



     subroutine make_dx(A,j,B,Aw)           ! substitute column j in A
      real(8)    :: A(3,3),B(3)
      integer    :: j
      real(8)    :: Aw(3,3)
      Aw(1:3,1:3) = A(1:3,1:3)
      do i=1,3
       Aw(j,i) = B(i)
      enddo
     end subroutine make_dx



    
     subroutine det3x3(A,D)                 ! determinant of 3x3
      real(8)  :: A(3,3)
      real(8)  :: D
      D = A(1,1)*A(2,2)*A(3,3)  &
        + A(1,2)*A(2,3)*A(3,1)  &
        + A(1,3)*A(2,1)*A(3,2)  &
        - A(1,3)*A(2,2)*A(3,1)  &
        - A(1,1)*A(2,3)*A(3,2)  &
        - A(1,2)*A(2,1)*A(3,3)  
     end subroutine det3x3




     subroutine check_solution(A,B,X)       ! A*X=B
      real(8)  :: A(3,3)
      real(8)  :: B(3),X(3)
      real(8)  :: R
      print *
      print *,'check solution A*X=B (inaccuracy)'
      do j=1,3
       R = A(1,j)*X(1)+A(2,j)*X(2)+A(3,j)*X(3)
       print 1,abs(R-B(j))
      enddo
1     format(F22.18)
     end subroutine check_solution













      subroutine check_uvw(u1,v1,w1,Aa,Ab,Ac)
       real(8)  :: u1(3),v1(3),w1(3)
       real(8)  :: Aa(3),Ab(3),Ac(3)
       real(8)  :: A1(3),A2(3),A3(3)
       print *
       print *,'test orthonormality of vectors built from u,v,w'
       A1(1:3) = u1(1)*Aa(1:3) + v1(1)*Ab(1:3) + w1(1)*Ac(1:3)
       A2(1:3) = u1(2)*Aa(1:3) + v1(2)*Ab(1:3) + w1(2)*Ac(1:3)
       A3(1:3) = u1(3)*Aa(1:3) + v1(3)*Ab(1:3) + w1(3)*Ac(1:3)
       call test_orthog(A1,A2,A3)
      end subroutine check_uvw











