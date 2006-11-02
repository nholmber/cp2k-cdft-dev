      	subroutine ppinit
	use atom,only: lmax
        use pspot
      	implicit none
	integer i,j,k,l
c
	if (pptype.eq.1) call c2a
c..BHS->general form
	if (pptype.le.2) then
  	  do i=1,ERFnr
	    PPNerf(i)=sqrt(PPNerf(i))
	  enddo
	  do l=0,PPlmax+1
	    do i=1,EXPnr(l)
	      PPNexp(i+EXPnr(l),l)=PPNexp(i,l)
	      PPRexp(i,l)=0
	      PPRexp(i+EXPnr(l),l)=2
	    enddo
	    EXPnr(l)=2*EXPnr(l)
	  enddo
	endif
c..GTH-seperable pseudopotential
	if (pptype.eq.5) then
          do i=1,EXPnr(0)
            PPNexp(i,0)=PPNexp(1,0)
            PPRexp(i,0)=2*(i-1)
            do l=1,lmax
	      PPLexp(i,l)=PPLexp(i,0)
	      PPNexp(i,l)=PPNexp(i,0)
	      PPRexp(i,l)=PPRexp(i,0)
	    enddo
	  enddo
	endif
c
      	return
      	end
