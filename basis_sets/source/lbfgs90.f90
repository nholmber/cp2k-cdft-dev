!------------------------------------------------------------------------------!
  subroutine lbfgs(v,x,g,h,ldx,n,m,m0,mmax)
!------------------------------------------------------------------------------!
            !     Limited Memory - BFGS Update method
            !     J. Nocedal  ,  Math. Comp. 35, 773 (1980)                 
            !     On Input : V(*) is the vector to be multiplied by H    
            !                X(LDX,MMAX) the position vectors         
            !                G(LDX,MMAX) the gradient vectors         
            !                H(*) Diagonal approximation to inverse Hessian 
            !                S(*) an auxillary vector of length M             
            !                LDX leading dimension of X and G    
            !                N is the dimension of the problem 
            !                M is the number of updates to be performed 
            !                M0 is the position of the last position vector  
            !                   within X                      
            !                MMAX is the maximum number of vectors in the buffer
            !
            !     On Output: V(*) containes the product H(m)*V(*)            
            !------------------------------------------------------------------!

  USE basic_data_types, ONLY: dp

  IMPLICIT NONE

!in:
  integer, intent(in)   :: ldx,n,m,m0,mmax
!inout:
  real(dp)              :: v(n),x(ldx,mmax),g(ldx,mmax),h(n)
!locals:
  real(dp)              :: s(0:m-1),b
  integer               :: i,li,lip
!------------------------------------------------------------------------------!
!input test:
  if (n>ldx) stop 'error in lbfgs: n>ldx !'
  if (m>mmax) stop 'error in lbfgs: m>mmax !'

!downward recursion:
  do i=m-1,0,-1
    li=m0-(m-i)
    if (li<=0) li=li+mmax
    lip=m0-(m-i)+1
    if (lip<=0) lip=lip+mmax
    s(i)=dot_product(v,x(1:n,lip)-x(1:n,li))
    v=v-s(i)*(g(1:n,lip)-g(1:n,li))
  enddo
!multiplication with guessed inverse Hessian:
  v=h*v
!upward recursion:
  do i=0,m-1
    li=m0-(m-i)
    if (li<=0) li=li+mmax
    lip=m0-(m-i)+1
    if (lip<=0) lip=lip+mmax
    b=dot_product(v,g(1:n,lip)-g(1:n,li))
    v=v+(s(i)-b)*(x(1:n,lip)-x(1:n,li))
  enddo

end subroutine lbfgs
