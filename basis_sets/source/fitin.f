C  From inet!cs.utexas.edu!cline Tue Oct 31 17:10:31 CST 1989
C  Received: from mojave.cs.utexas.edu by cs.utexas.edu (5.59/1.44)
C  	id AA29509; Tue, 31 Oct 89 17:11:51 CST
C  Posted-Date: Tue, 31 Oct 89 17:10:31 CST
C  Message-Id: <8910312310.AA04442@mojave.cs.utexas.edu>
C  Received: by mojave.cs.utexas.edu (14.5/1.4-Client)
C  	id AA04442; Tue, 31 Oct 89 17:10:34 cst
C  Date: Tue, 31 Oct 89 17:10:31 CST
C  X-Mailer: Mail User's Shell (6.5 4/17/89)
C  From: cline@cs.utexas.edu (Alan Cline)
C  To: ehg@research.att.com
C  Subject: New FITPACK Subset for netlib
C  
C  
C  This new version of FITPACK distributed by netlib is about 20% of 
C  the total package in terms of characters, lines of code, and num-
C  ber of subprograms. However, these 25 subprograms represent about
C  95% of usages of the package.  What has been omitted are such ca-
C  pabilities as:
C    1. Automatic tension determination,
C    2. Derivatives, arclengths, and enclosed areas for planar 
C       curves,
C    3. Three dimensional curves,
C    4. Special surface fitting using equispacing assumptions,
C    5. Surface fitting in annular, wedge, polar, toroidal, lunar,
C       and spherical geometries,
C    6. B-splines in tension generation and usage,
C    7. General surface fitting in three dimensional space.
C  
C  (The code previously circulated in netlib is less than 10% of the
C  total  package  and is more than a decade old.  Its usage is dis-
C  couraged.)
C  
C  Please note:  Two versions of the subroutine snhcsh are included.
C  Both serve the same purpose:  obtaining approximations to certain
C  hyperbolic trigonometric-like functions.  The first is less accu-
C  rate (but more efficient) than the second.  Installers should se- 
C  lect the one with the precision they desire.
C  
C  Interested parties can obtain the entire package on disk or  tape
C  from Pleasant  Valley Software, 8603 Altus Cove, Austin TX (USA),
C  78759 at a cost of $495 US. A 340 page manual  is  available  for
C  $30  US  per  copy.  The  package  includes  examples and machine
C  readable documentation.
C  
c------------------------------------------------------------
c------------------------------------------------------------
c------------------------------------------------------------
      subroutine curv1 (n,x,y,slp1,slpn,islpsw,yp,temp,
     *                  sigma,ierr)
c
      implicit real*8 (a-h,o-z)
      integer n,islpsw,ierr
      real*8 x(n),y(n),slp1,slpn,yp(n),temp(n),sigma
c
c                                 coded by alan kaylor cline
c                           from fitpack -- january 26, 1987
c                        a curve and surface fitting package
c                      a product of pleasant valley software
c                  8603 altus cove, austin, texas 78759, usa
c
c this subroutine determines the parameters necessary to
c compute an interpolatory spline under tension through
c a sequence of functional values. the slopes at the two
c ends of the curve may be specified or omitted.  for actual
c computation of points on the curve it is necessary to call
c the function curv2.
c
c on input--
c
c   n is the number of values to be interpolated (n.ge.2).
c
c   x is an array of the n increasing abscissae of the
c   functional values.
c
c   y is an array of the n ordinates of the values, (i. e.
c   y(k) is the functional value corresponding to x(k) ).
c
c   slp1 and slpn contain the desired values for the first
c   derivative of the curve at x(1) and x(n), respectively.
c   the user may omit values for either or both of these
c   parameters and signal this with islpsw.
c
c   islpsw contains a switch indicating which slope data
c   should be used and which should be estimated by this
c   subroutine,
c          = 0 if slp1 and slpn are to be used,
c          = 1 if slp1 is to be used but not slpn,
c          = 2 if slpn is to be used but not slp1,
c          = 3 if both slp1 and slpn are to be estimated
c              internally.
c
c   yp is an array of length at least n.
c
c   temp is an array of length at least n which is used for
c   scratch storage.
c
c and
c
c   sigma contains the tension factor. this value indicates
c   the curviness desired. if abs(sigma) is nearly zero
c   (e.g. .001) the resulting curve is approximately a
c   cubic spline. if abs(sigma) is large (e.g. 50.) the
c   resulting curve is nearly a polygonal line. if sigma
c   equals zero a cubic spline results.  a standard value
c   for sigma is approximately 1. in absolute value.
c
c on output--
c
c   yp contains the values of the second derivative of the
c   curve at the given nodes.
c
c   ierr contains an error flag,
c        = 0 for normal return,
c        = 1 if n is less than 2,
c        = 2 if x-values are not strictly increasing.
c
c and
c
c   n, x, y, slp1, slpn, islpsw and sigma are unaltered.
c
c this subroutine references package modules ceez, terms,
c and snhcsh.
c
c-----------------------------------------------------------
c
      nm1 = n-1
      np1 = n+1
      ierr = 0
      if (n .le. 1) go to 8
      if (x(n) .le. x(1)) go to 9
c
c denormalize tension factor
c
      sigmap = abs(sigma)*float(n-1)/(x(n)-x(1))
c
c approximate end slopes
c
      if (islpsw .ge. 2) go to 1
      slpp1 = slp1
      go to 2
    1 delx1 = x(2)-x(1)
      delx2 = delx1+delx1
      if (n .gt. 2) delx2 = x(3)-x(1)
      if (delx1 .le. 0. .or. delx2 .le. delx1) go to 9
      call ceez (delx1,delx2,sigmap,c1,c2,c3,n)
      slpp1 = c1*y(1)+c2*y(2)
      if (n .gt. 2) slpp1 = slpp1+c3*y(3)
    2 if (islpsw .eq. 1 .or. islpsw .eq. 3) go to 3
      slppn = slpn
      go to 4
    3 delxn = x(n)-x(nm1)
      delxnm = delxn+delxn
      if (n .gt. 2) delxnm = x(n)-x(n-2)
      if (delxn .le. 0. .or. delxnm .le. delxn) go to 9
      call ceez (-delxn,-delxnm,sigmap,c1,c2,c3,n)
      slppn = c1*y(n)+c2*y(nm1)
      if (n .gt. 2) slppn = slppn+c3*y(n-2)
c
c set up right hand side and tridiagonal system for yp and
c perform forward elimination
c
    4 delx1 = x(2)-x(1)
      if (delx1 .le. 0.) go to 9
      dx1 = (y(2)-y(1))/delx1
      call terms (diag1,sdiag1,sigmap,delx1)
      yp(1) = (dx1-slpp1)/diag1
      temp(1) = sdiag1/diag1
      if (n .eq. 2) go to 6
      do 5 i = 2,nm1
        delx2 = x(i+1)-x(i)
        if (delx2 .le. 0.) go to 9
        dx2 = (y(i+1)-y(i))/delx2
        call terms (diag2,sdiag2,sigmap,delx2)
        diag = diag1+diag2-sdiag1*temp(i-1)
        yp(i) = (dx2-dx1-sdiag1*yp(i-1))/diag
        temp(i) = sdiag2/diag
        dx1 = dx2
        diag1 = diag2
    5   sdiag1 = sdiag2
    6 diag = diag1-sdiag1*temp(nm1)
      yp(n) = (slppn-dx1-sdiag1*yp(nm1))/diag
c
c perform back substitution
c
      do 7 i = 2,n
        ibak = np1-i
    7   yp(ibak) = yp(ibak)-temp(ibak)*yp(ibak+1)
      return
c
c too few points
c
    8 ierr = 1
      return
c
c x-values not strictly increasing
c
    9 ierr = 2
      return
      end
c
c------------------------------------------------------------
c------------------------------------------------------------
c------------------------------------------------------------
      function curv2 (t,n,x,y,yp,sigma)
c
      implicit real*8 (a-h,o-z)
      integer n
      real*8 t,x(n),y(n),yp(n),sigma
c
c                                 coded by alan kaylor cline
c                           from fitpack -- january 26, 1987
c                        a curve and surface fitting package
c                      a product of pleasant valley software
c                  8603 altus cove, austin, texas 78759, usa
c
c this function interpolates a curve at a given point
c using a spline under tension. the subroutine curv1 should
c be called earlier to determine certain necessary
c parameters.
c
c on input--
c
c   t contains a real value to be mapped onto the interpo-
c   lating curve.
c
c   n contains the number of points which were specified to
c   determine the curve.
c
c   x and y are arrays containing the abscissae and
c   ordinates, respectively, of the specified points.
c
c   yp is an array of second derivative values of the curve
c   at the nodes.
c
c and
c
c   sigma contains the tension factor (its sign is ignored).
c
c the parameters n, x, y, yp, and sigma should be input
c unaltered from the output of curv1.
c
c on output--
c
c   curv2 contains the interpolated value.
c
c none of the input parameters are altered.
c
c this function references package modules intrvl and
c snhcsh.
c
c-----------------------------------------------------------
c
c determine interval
c
      im1 = intrvl(t,x,n)
      i = im1+1
c
c denormalize tension factor
c
      sigmap = abs(sigma)*float(n-1)/(x(n)-x(1))
c
c set up and perform interpolation
c
      del1 = t-x(im1)
      del2 = x(i)-t
      dels = x(i)-x(im1)
      sum = (y(i)*del1+y(im1)*del2)/dels
      if (sigmap .ne. 0.) go to 1
      curv2 = sum-del1*del2*(yp(i)*(del1+dels)+yp(im1)*
     *        (del2+dels))/(6.*dels)
      return
    1 sigdel = sigmap*dels
      call snhcsh (ss,dummy,sigdel,-1)
      call snhcsh (s1,dummy,sigmap*del1,-1)
      call snhcsh (s2,dummy,sigmap*del2,-1)
      curv2 = sum+(yp(i)*del1*(s1-ss)+yp(im1)*del2*(s2-ss))/
     *            (sigdel*sigmap*(1.+ss))
      return
      end
c
      function curvd (t,n,x,y,yp,sigma)
c
      implicit real*8 (a-h,o-z)
      integer n
      real*8 t,x(n),y(n),yp(n),sigma
c
c                                 coded by alan kaylor cline
c                           from fitpack -- january 26, 1987
c                        a curve and surface fitting package
c                      a product of pleasant valley software
c                  8603 altus cove, austin, texas 78759, usa
c
c this function differentiates a curve at a given point
c using a spline under tension. the subroutine curv1 should
c be called earlier to determine certain necessary
c parameters.
c
c on input--
c
c   t contains a real value at which the derivative is to be
c   determined.
c
c   n contains the number of points which were specified to
c   determine the curve.
c
c   x and y are arrays containing the abscissae and
c   ordinates, respectively, of the specified points.
c
c   yp is an array of second derivative values of the curve
c   at the nodes.
c
c and
c
c   sigma contains the tension factor (its sign is ignored).
c
c the parameters n, x, y, yp, and sigma should be input
c unaltered from the output of curv1.
c
c on output--
c
c   curvd contains the derivative value.
c
c none of the input parameters are altered.
c
c this function references package modules intrvl and
c snhcsh.
c
c-----------------------------------------------------------
c
c determine interval
c
      im1 = intrvl(t,x,n)
      i = im1+1
c
c denormalize tension factor
c
      sigmap = abs(sigma)*float(n-1)/(x(n)-x(1))
c
c set up and perform differentiation
c
      del1 = t-x(im1)
      del2 = x(i)-t
      dels = x(i)-x(im1)
      sum = (y(i)-y(im1))/dels
      if (sigmap .ne. 0.) go to 1
      curvd = sum+(yp(i)*(2.*del1*del1-del2*(del1+dels))-
     *             yp(im1)*(2.*del2*del2-del1*(del2+dels)))
     *             /(6.*dels)
      return
    1 sigdel = sigmap*dels
      call snhcsh (ss,dummy,sigdel,-1)
      call snhcsh (dummy,c1,sigmap*del1,1)
      call snhcsh (dummy,c2,sigmap*del2,1)
      curvd = sum+(yp(i)*(c1-ss)-yp(im1)*(c2-ss))/
     *        (sigdel*sigmap*(1.+ss))
      return
      end
      function curvi (xl,xu,n,x,y,yp,sigma)
c
      implicit real*8 (a-h,o-z)
      integer n
      real*8 xl,xu,x(n),y(n),yp(n),sigma
c
c                                 coded by alan kaylor cline
c                           from fitpack -- january 26, 1987
c                        a curve and surface fitting package
c                      a product of pleasant valley software
c                  8603 altus cove, austin, texas 78759, usa
c
c this function integrates a curve specified by a spline
c under tension between two given limits. the subroutine
c curv1 should be called earlier to determine necessary
c parameters.
c
c on input--
c
c   xl and xu contain the upper and lower limits of inte-
c   gration, respectively. (sl need not be less than or
c   equal to xu, curvi (xl,xu,...) .eq. -curvi (xu,xl,...) ).
c
c   n contains the number of points which were specified to
c   determine the curve.
c
c   x and y are arrays containing the abscissae and
c   ordinates, respectively, of the specified points.
c
c   yp is an array from subroutine curv1 containing
c   the values of the second derivatives at the nodes.
c
c and
c
c   sigma contains the tension factor (its sign is ignored).
c
c the parameters n, x, y, yp, and sigma should be input
c unaltered from the output of curv1.
c
c on output--
c
c   curvi contains the integral value.
c
c none of the input parameters are altered.
c
c this function references package modules intrvl and
c snhcsh.
c
c-----------------------------------------------------------
c
c denormalize tension factor
c
      sigmap = abs(sigma)*float(n-1)/(x(n)-x(1))
c
c determine actual upper and lower bounds
c
      xxl = xl
      xxu = xu
      ssign = 1.
      if (xl .lt. xu) go to 1
      xxl = xu
      xxu = xl
      ssign = -1.
      if (xl .gt. xu) go to 1
c
c return zero if xl .eq. xu
c
      curvi = 0.
      return
c
c search for proper intervals
c
    1 ilm1 = intrvl (xxl,x,n)
      il = ilm1+1
      ium1 = intrvl (xxu,x,n)
      iu = ium1+1
      if (il .eq. iu) go to 8
c
c integrate from xxl to x(il)
c
      sum = 0.
      if (xxl .eq. x(il)) go to 3
      del1 = xxl-x(ilm1)
      del2 = x(il)-xxl
      dels = x(il)-x(ilm1)
      t1 = (del1+dels)*del2/(2.*dels)
      t2 = del2*del2/(2.*dels)
      sum = t1*y(il)+t2*y(ilm1)
      if (sigma .eq. 0.) go to 2
      call snhcsh (dummy,c1,sigmap*del1,2)
      call snhcsh (dummy,c2,sigmap*del2,2)
      call snhcsh (ss,cs,sigmap*dels,3)
      sum = sum+((dels*dels*(cs-ss/2.)-del1*del1*(c1-ss/2.))
     *           *yp(il)+del2*del2*(c2-ss/2.)*yp(ilm1))/
     *          (sigmap*sigmap*dels*(1.+ss))
      go to 3
    2 sum = sum-t1*t1*dels*yp(il)/6.
     *         -t2*(del1*(del2+dels)+dels*dels)*yp(ilm1)/12.
c
c integrate over interior intervals
c
    3 if (iu-il .eq. 1) go to 6
      ilp1 = il+1
      do 5 i = ilp1,ium1
        dels = x(i)-x(i-1)
        sum = sum+(y(i)+y(i-1))*dels/2.
        if (sigma .eq. 0.) go to 4
        call snhcsh (ss,cs,sigmap*dels,3)
        sum = sum+(yp(i)+yp(i-1))*dels*(cs-ss/2.)/
     *            (sigmap*sigmap*(1.+ss))
        go to 5
    4   sum = sum-(yp(i)+yp(i-1))*dels*dels*dels/24.
    5   continue
c
c integrate from x(iu-1) to xxu
c
    6 if (xxu .eq. x(ium1)) go to 10
      del1 = xxu-x(ium1)
      del2 = x(iu)-xxu
      dels = x(iu)-x(ium1)
      t1 = del1*del1/(2.*dels)
      t2 = (del2+dels)*del1/(2.*dels)
      sum = sum+t1*y(iu)+t2*y(ium1)
      if (sigma .eq. 0.) go to 7
      call snhcsh (dummy,c1,sigmap*del1,2)
      call snhcsh (dummy,c2,sigmap*del2,2)
      call snhcsh (ss,cs,sigmap*dels,3)
      sum = sum+(yp(iu)*del1*del1*(c1-ss/2.)+yp(ium1)*
     *          (dels*dels*(cs-ss/2.)-del2*del2*(c2-ss/2.)))
     *         /(sigmap*sigmap*dels*(1.+ss))
      go to 10
    7 sum = sum-t1*(del2*(del1+dels)+dels*dels)*yp(iu)/12.
     *         -t2*t2*dels*yp(ium1)/6.
      go to 10
c
c integrate from xxl to xxu
c
    8 delu1 = xxu-x(ium1)
      delu2 = x(iu)-xxu
      dell1 = xxl-x(ium1)
      dell2 = x(iu)-xxl
      dels = x(iu)-x(ium1)
      deli = xxu-xxl
      t1 = (delu1+dell1)*deli/(2.*dels)
      t2 = (delu2+dell2)*deli/(2.*dels)
      sum = t1*y(iu)+t2*y(ium1)
      if (sigma .eq. 0.) go to 9
      call snhcsh (dummy,cu1,sigmap*delu1,2)
      call snhcsh (dummy,cu2,sigmap*delu2,2)
      call snhcsh (dummy,cl1,sigmap*dell1,2)
      call snhcsh (dummy,cl2,sigmap*dell2,2)
      call snhcsh (ss,dummy,sigmap*dels,-1)
      sum = sum+(yp(iu)*(delu1*delu1*(cu1-ss/2.)
     *            -dell1*dell1*(cl1-ss/2.))
     *          +yp(ium1)*(dell2*dell2*(cl2-ss/2.)
     *            -delu2*delu2*(cu2-ss/2.)))/
     *          (sigmap*sigmap*dels*(1.+ss))
      go to 10
    9 sum = sum-t1*(delu2*(dels+delu1)+dell2*(dels+dell1))*
     *             yp(iu)/12.
     *         -t2*(dell1*(dels+dell2)+delu1*(dels+delu2))*
     *             yp(ium1)/12.
c
c correct sign and return
c
   10 curvi = ssign*sum
      return
      end
c------------------------------------------------------------
c------------------------------------------------------------
c------------------------------------------------------------
      subroutine ceez (del1,del2,sigma,c1,c2,c3,n)
c
      implicit real*8 (a-h,o-z)
      real*8 del1,del2,sigma,c1,c2,c3
c
c                                 coded by alan kaylor cline
c                           from fitpack -- january 26, 1987
c                        a curve and surface fitting package
c                      a product of pleasant valley software
c                  8603 altus cove, austin, texas 78759, usa
c
c this subroutine determines the coefficients c1, c2, and c3
c used to determine endpoint slopes. specifically, if
c function values y1, y2, and y3 are given at points x1, x2,
c and x3, respectively, the quantity c1*y1 + c2*y2 + c3*y3
c is the value of the derivative at x1 of a spline under
c tension (with tension factor sigma) passing through the
c three points and having third derivative equal to zero at
c x1. optionally, only two values, c1 and c2 are determined.
c
c on input--
c
c   del1 is x2-x1 (.gt. 0.).
c
c   del2 is x3-x1 (.gt. 0.). if n .eq. 2, this parameter is
c   ignored.
c
c   sigma is the tension factor.
c
c and
c
c   n is a switch indicating the number of coefficients to
c   be returned. if n .eq. 2 only two coefficients are
c   returned. otherwise all three are returned.
c
c on output--
c
c   c1, c2, and c3 contain the coefficients.
c
c none of the input parameters are altered.
c
c this subroutine references package module snhcsh.
c
c-----------------------------------------------------------
c
      if (n .eq. 2) go to 2
      if (sigma .ne. 0.) go to 1
      del = del2-del1
c
c tension .eq. 0.
c
      c1 = -(del1+del2)/(del1*del2)
      c2 = del2/(del1*del)
      c3 = -del1/(del2*del)
      return
c
c tension .ne. 0.
c
    1 call snhcsh (dummy,coshm1,sigma*del1,1)
      call snhcsh (dummy,coshm2,sigma*del2,1)
      delp = sigma*(del2+del1)/2.
      delm = sigma*(del2-del1)/2.
      call snhcsh (sinhmp,dummy,delp,-1)
      call snhcsh (sinhmm,dummy,delm,-1)
      denom = coshm1*(del2-del1)-2.*del1*delp*delm*
     *        (1.+sinhmp)*(1.+sinhmm)
      c1 = 2.*delp*delm*(1.+sinhmp)*(1.+sinhmm)/denom
      c2 = -coshm2/denom
      c3 = coshm1/denom
      return
c
c two coefficients
c
    2 c1 = -1./del1
      c2 = -c1
      return
      end
c
c------------------------------------------------------------
c------------------------------------------------------------
c------------------------------------------------------------
      function intrvl (t,x,n)
c
      implicit real*8 (a-h,o-z)
      integer n
      real*8 t,x(n)
c
c                                 coded by alan kaylor cline
c                           from fitpack -- january 26, 1987
c                        a curve and surface fitting package
c                      a product of pleasant valley software
c                  8603 altus cove, austin, texas 78759, usa
c
c this function determines the index of the interval
c (determined by a given increasing sequence) in which
c a given value lies.
c
c on input--
c
c   t is the given value.
c
c   x is a vector of strictly increasing values.
c
c and
c
c   n is the length of x (n .ge. 2).
c
c on output--
c
c   intrvl returns an integer i such that
c
c          i =  1       if         e   t .lt. x(2)  ,
c          i =  n-1     if x(n-1) .le. t            ,
c          otherwise       x(i)  .le. t .le. x(i+1),
c
c none of the input parameters are altered.
c
c-----------------------------------------------------------
c
      save i
      data i /1/
c
      tt = t
c
c check for illegal i
c
      if (i .ge. n) i = n/2
c
c check old interval and extremes
c
      if (tt .lt. x(i)) then
        if (tt .le. x(2)) then
          i = 1
          intrvl = 1
          return
        else
          il = 2
          ih = i
        end if
      else if (tt .le. x(i+1)) then
        intrvl = i
        return
      else if (tt .ge. x(n-1)) then
        i = n-1
        intrvl = n-1
        return
      else
        il = i+1
        ih = n-1
      end if
c
c binary search loop
c
    1 i = (il+ih)/2
      if (tt .lt. x(i)) then
         ih = i
      else if (tt .gt. x(i+1)) then
         il = i+1
      else
         intrvl = i
         return
      end if
      go to 1
      end
c
c------------------------------------------------------------
c------------------------------------------------------------
c------------------------------------------------------------
      subroutine snhcsh (sinhm,coshm,x,isw)
c
      implicit real*8 (a-h,o-z)
      integer isw
      real*8 sinhm,coshm,x
c
c                                 coded by alan kaylor cline
c                           from fitpack -- january 26, 1987
c                        a curve and surface fitting package
c                      a product of pleasant valley software
c                  8603 altus cove, austin, texas 78759, usa
c
c this subroutine returns approximations to
c       sinhm(x) = sinh(x)/x-1
c       coshm(x) = cosh(x)-1
c and
c       coshmm(x) = (cosh(x)-1-x*x/2)/(x*x)
c with relative error less than 4.0e-14.
c
c on input--
c
c   x contains the value of the independent variable.
c
c   isw indicates the function desired
c           = -1 if only sinhm is desired,
c           =  0 if both sinhm and coshm are desired,
c           =  1 if only coshm is desired,
c           =  2 if only coshmm is desired,
c           =  3 if both sinhm and coshmm are desired.
c
c on output--
c
c   sinhm contains the value of sinhm(x) if isw .le. 0 or
c   isw .eq. 3 (sinhm is unaltered if isw .eq.1 or isw .eq.
c   2).
c
c   coshm contains the value of coshm(x) if isw .eq. 0 or
c   isw .eq. 1 and contains the value of coshmm(x) if isw
c   .ge. 2 (coshm is unaltered if isw .eq. -1).
c
c and
c
c   x and isw are unaltered.
c
c-----------------------------------------------------------
c
      data sp14/.227581660976348e-7/,
     *     sp13/.612189863171694e-5/,
     *     sp12/.715314759211209e-3/,
     *     sp11/.398088289992973e-1/,
     *     sq12/.206382701413725e-3/,
     *     sq11/-.611470260009508e-1/,
     *     sq10/.599999999999986e+1/
      data sp25/.129094158037272e-9/,
     *     sp24/.473731823101666e-7/,
     *     sp23/.849213455598455e-5/,
     *     sp22/.833264803327242e-3/,
     *     sp21/.425024142813226e-1/,
     *     sq22/.106008515744821e-3/,
     *     sq21/-.449855169512505e-1/,
     *     sq20/.600000000268619e+1/
      data sp35/.155193945864942e-9/,
     *     sp34/.511529451668737e-7/,
     *     sp33/.884775635776784e-5/,
     *     sp32/.850447617691392e-3/,
     *     sp31/.428888148791777e-1/,
     *     sq32/.933128831061610e-4/,
     *     sq31/-.426677570538507e-1/,
     *     sq30/.600000145086489e+1/
      data sp45/.188070632058331e-9/,
     *     sp44/.545792817714192e-7/,
     *     sp43/.920119535795222e-5/,
     *     sp42/.866559391672985e-3/,
     *     sp41/.432535234960858e-1/,
     *     sq42/.824891748820670e-4/,
     *     sq41/-.404938841672262e-1/,
     *     sq40/.600005006283834e+1/
      data cp5/.552200614584744e-9/,
     *     cp4/.181666923620944e-6/,
     *     cp3/.270540125846525e-4/,
     *     cp2/.206270719503934e-2/,
     *     cp1/.744437205569040e-1/,
     *     cq2/.514609638642689e-4/,
     *     cq1/-.177792255528382e-1/,
     *     cq0/.200000000000000e+1/
      data zp4/.664418805876835e-8/,
     *     zp3/.218274535686385e-5/,
     *     zp2/.324851059327161e-3/,
     *     zp1/.244515150174258e-1/,
     *     zq2/.616165782306621e-3/,
     *     zq1/-.213163639579425e0/,
     *     zq0/.240000000000000e+2/
c
      ax = abs(x)
      if (isw .ge. 0) go to 5
c
c sinhm approximation
c
      if (ax .gt. 3.9) go to 2
      xs = ax*ax
      if (ax .gt. 2.2) go to 1
c
c sinhm approximation on (0.,2.2)
c
      sinhm = xs*((((sp14*xs+sp13)*xs+sp12)*xs+sp11)*xs+1.)/
     .             ((sq12*xs+sq11)*xs+sq10)
      return
c
c sinhm approximation on (2.2,3.9)
c
    1 sinhm = xs*(((((sp25*xs+sp24)*xs+sp23)*xs+sp22)*xs+sp21)
     .        *xs+1.)/((sq22*xs+sq21)*xs+sq20)
      return
    2 if (ax .gt. 5.1) go to 3
c
c sinhm approximation on (3.9,5.1)
c
      xs = ax*ax
      sinhm = xs*(((((sp35*xs+sp34)*xs+sp33)*xs+sp32)*xs+sp31)
     .        *xs+1.)/((sq32*xs+sq31)*xs+sq30)
      return
    3 if (ax .gt. 6.1) go to 4
c
c sinhm approximation on (5.1,6.1)
c
      xs = ax*ax
      sinhm = xs*(((((sp45*xs+sp44)*xs+sp43)*xs+sp42)*xs+sp41)
     .        *xs+1.)/((sq42*xs+sq41)*xs+sq40)
      return
c
c sinhm approximation above 6.1
c
    4 expx = exp(ax)
      sinhm = (expx-1./expx)/(ax+ax)-1.
      return
c
c coshm and (possibly) sinhm approximation
c
    5 if (isw .ge. 2) go to 7
      if (ax .gt. 2.2) go to 6
      xs = ax*ax
      coshm = xs*(((((cp5*xs+cp4)*xs+cp3)*xs+cp2)*xs+cp1)
     .        *xs+1.)/((cq2*xs+cq1)*xs+cq0)
      if (isw .eq. 0) sinhm = xs*((((sp14*xs+sp13)*xs+sp12)
     .          *xs+sp11)*xs+1.)/((sq12*xs+sq11)*xs+sq10)
      return
    6 expx = exp(ax)
      coshm = (expx+1./expx)/2.-1.
      if (isw .eq. 0) sinhm = (expx-1./expx)/(ax+ax)-1.
      return
c
c coshmm and (possibly) sinhm approximation
c
    7 xs = ax*ax
      if (ax .gt. 2.2) go to 8
      coshm = xs*((((zp4*xs+zp3)*xs+zp2)*xs+zp1)*xs+1.)/
     .             ((zq2*xs+zq1)*xs+zq0)
      if (isw .eq. 3) sinhm = xs*((((sp14*xs+sp13)*xs+sp12)
     .          *xs+sp11)*xs+1.)/((sq12*xs+sq11)*xs+sq10)
      return
    8 expx = exp(ax)
      coshm = ((expx+1./expx-xs)/2.-1.)/xs
      if (isw .eq. 3) sinhm = (expx-1./expx)/(ax+ax)-1.
      return
      end
c------------------------------------------------------------
c------------------------------------------------------------
c------------------------------------------------------------
      subroutine terms (diag,sdiag,sigma,del)
c
      implicit real*8 (a-h,o-z)
      real*8 diag,sdiag,sigma,del
c
c                                 coded by alan kaylor cline
c                           from fitpack -- january 26, 1987
c                        a curve and surface fitting package
c                      a product of pleasant valley software
c                  8603 altus cove, austin, texas 78759, usa
c
c this subroutine computes the diagonal and superdiagonal
c terms of the tridiagonal linear system associated with
c spline under tension interpolation.
c
c on input--
c
c   sigma contains the tension factor.
c
c and
c
c   del contains the step size.
c
c on output--
c
c                sigma*del*cosh(sigma*del) - sinh(sigma*del)
c   diag = del*--------------------------------------------.
c                     (sigma*del)**2 * sinh(sigma*del)
c
c                   sinh(sigma*del) - sigma*del
c   sdiag = del*----------------------------------.
c                (sigma*del)**2 * sinh(sigma*del)
c
c and
c
c   sigma and del are unaltered.
c
c this subroutine references package module snhcsh.
c
c-----------------------------------------------------------
c
      if (sigma .ne. 0.) go to 1
      diag = del/3.
      sdiag = del/6.
      return
    1 sigdel = sigma*del
      call snhcsh (sinhm,coshm,sigdel,0)
      denom = sigma*sigdel*(1.+sinhm)
      diag = (coshm-sinhm)/denom
      sdiag = sinhm/denom
      return
      end

