*******************************************
*******************************************

!           GRID-1.00

*******************************************
*******************************************

      subroutine grid

      include 'param-1.00.inc'
      include 'com-1.00.inc'

!     s grid parameters: gridding for parallel dynamics

      real s(nz,nf)
      real qs(nz,nf)
      real xs(nz,nf),ys(nz,nf),zs(nz,nf)
      real brs(nz,nf)

!     p grid parameters: gridding for perpendicular dynamics

      real qp(nzp1,nfp1),pp(nzp1,nfp1)
      real xp(nzp1,nfp1,2),yp(nzp1,nfp1,2),zp(nzp1,nfp1,2)
      real brp(nzp1,nfp1),blatp(nzp1,nfp1),blonp(nzp1,nfp1)
      real grp(nzp1,nfp1),glatp(nzp1,nfp1),glonp(nzp1,nfp1)
      real altp(nzp1,nfp1)

      do j = 1,nf
        k    = nf + 1 - j
        dx   = float(j-1)/float(nf-1)
        x2   = dx * sinh ( gamp )
        x3   = alog ( x2 + sqrt ( x2**2 + 1. ) ) / gamp
        grad_inp(k) = rmin + ( rmax - rmin ) * ( 1. - x3 )
      enddo

      call igrf_sub(glon_in)

! MS: Incorporate eccentric dipole fit to IGRF.  The fit was done by 
! taken swaths in geographic longitude and fitting the best eccentric 
! dipole.  We want the simulation plane to go through the point 
! specified in the namelist, but the fit is done in terms of 
! geographic longitude at the equator.  So iterate: Take a guess at 
! an initial geographic longitude and find the eccentric dipole 
! parameters. Find the magnetic longitude (use gtob_norm) of the 
! namelist point and a point at the equator.  Are the two magnetic 
! longitudes the same?  If not, adjust the guess for the geographic 
! longitude and try again.  Rinse, repeat, wipe hands on pants. 
 
! MS: Note that since the magnetic field is tilted (and SAMI models 
! one magnetic longitude) points far away from the equator will not be 
! correct.  They are at different geographic longitudes so their IGRF 
! parameters will be different (and hence the g to b to g conversions 
! will be different). I'm not yet sure how significant a problem this 
! is. 
 
       phitemp = 0. 

       do j = 1,20 
         call igrf_sub(phitemp) 
         call gtob ( temp1,temp2,temp3,grad_in+re,glat_in,glon_in ) 
         call gtob ( brad,blon0,blat,grad_in+re,0.,phitemp ) 
         phitemp = phitemp + (temp2-blon0) + 360. 
         phitemp = mod(phitemp,360.) 
       enddo 
 
! MS: Check for convergence 

       if (abs(temp2-blon0) .gt. 1.e-2) then 
         print *,'Failure to converge to initial magnetic latitude.' 
         stop 
       endif 


      call gtob (  brad,blon0,blat
     .            ,grad_in+re
     .            ,glat_in
     .            ,glon_in          )


      pvalue0  = brad / re / cos(blat*po180) ** 2
      rgmin    = altmin + re
        
      call btog ( pvalue0*re,blon0,0.
     .           ,grad_max
     .           ,glat_max
     .           ,glon_max             )

!     s grid

      do j = 1,nf

        i    = 0
        delp = .01 / re
        delg = grad_inp(j)

        do while ( abs(delg) .gt. .01 )
          if ( grad_inp(j) + re .lt. grad_max .and. j .eq. 1 ) then
                pvalue = pvalue0 - i * delp 
          else
                pvalue = pvalue0 + i * delp
          endif
          call btog ( pvalue*re,blon0,0.
     .               ,grad
     .               ,glat
     .               ,glon                 )
          delg = grad_inp(j) + re - grad
          i = i + 1
        enddo

        pvalue0 = pvalue

!        print *, ' j,grad_inp(j),pvalue',j,grad_inp(j),pvalue

!       north minimum: find q value at northern most point

        blata   = 0.
        blatb   = 60.
        blatc   = 0.5 * ( blata + blatb )
        rba     = pvalue * re * cos ( blata * po180 ) ** 2
        rbb     = pvalue * re * cos ( blatb * po180 ) ** 2
        rbc     = pvalue * re * cos ( blatc * po180 ) ** 2
        call  btog ( rba,blon0,blata,rga,glat,glon )
        call  btog ( rbb,blon0,blatb,rgb,glat,glon )
        call  btog ( rbc,blon0,blatc,rgc,glat,glon )
        delrg   = .01
        delrc   = abs ( rgc - rgmin )
        qminn   = ( re / rbc ) ** 2 * sin ( blatc * po180 )

        iminn   = 0

        do while ( delrc .gt. delrg .and. iminn .lt. 20)
          iminn   = iminn + 1
          if ( rgc .lt. rgmin ) blatb = blatc
          if ( rgc .gt. rgmin ) blata = blatc
          blatc   = 0.5 * ( blata + blatb )
          rbc     = pvalue * re * cos ( blatc * po180 ) ** 2
          call  btog ( rbc,blon0,blatc,rgc,glat,glon )
          delrc   = abs ( rgc - rgmin )
          qminn   = ( re / rbc ) ** 2 * sin ( blatc * po180 )
        enddo

!       south minimum: find q value at southern most point

        blata   = 0.
        blatb   = -60.
        blatc   = 0.5 * ( blata + blatb )
        rba     = pvalue * re * cos ( blata * po180 ) ** 2
        rbb     = pvalue * re * cos ( blatb * po180 ) ** 2
        rbc     = pvalue * re * cos ( blatc * po180 ) ** 2
        call  btog ( rba,blon0,blata,rga,glat,glon )
        call  btog ( rbb,blon0,blatb,rgb,glat,glon )
        call  btog ( rbc,blon0,blatc,rgc,glat,glon )

        delrg   = .01
        delrc   = abs ( rgc - rgmin )
        qmins   = ( re / rbc ) ** 2 * sin ( blatc * po180 )

        imins   = 0

        do while ( delrc .gt. delrg )
          imins   = imins + 1
          if ( rgc .lt. rgmin ) blatb = blatc
          if ( rgc .gt. rgmin ) blata = blatc
          blatc   = 0.5 * ( blata + blatb )
          rbc     = pvalue * re * cos ( blatc * po180 ) ** 2
          call  btog ( rbc,blon0,blatc,rgc,glat,glon )
          delrc   = abs ( rgc - rgmin )
          qmins   = ( re / rbc ) ** 2 * sin ( blatc * po180 )
          rbms    = rbc
        enddo

!       determine distribution of points (q,s) along field line
!       see red book p. 249 (millward et al.)

        r    = rbms

        nzh   = (nz-1)/2
        delqs = qmins
        delqn = qminn

!       firs half of field line

        do i = 1,nzh+1

          x    = -1. + ( float(i) - 1 ) / float(nzh)
          xx   = x * sinh ( gams * qmins )
          qtmp = -alog ( xx + sqrt ( xx**2 + 1. ) ) / gams

!         x    = -1. + 2. * ( float(i) - 1 ) / float(nz-1)
!         xx   = x/2. * ( sinh(gams*qminn) * (x+1) + 
!     .                   sinh(gams*qmins) * (x-1) )
!         qtmp = alog ( xx + sqrt ( xx**2 + 1. ) ) / gams
 
          r = re*qp_solve(qtmp,pvalue)
         
          qs(i,j)       = qtmp
          ps(i,j)       = pvalue
          brs(i,j)      = r
          br_norm       = r / re
          blats(i,j)  = asin ( qs(i,j) * br_norm ** 2 ) * rtod
          call btog ( brs(i,j),blon0,blats(i,j),
     .                grs(i,j),glats(i,j),glons(i,j) )

!         we need to calculate some quantities which depend on the usual
!         angular definition of theta

          theta = acos ( qs(i,j) * (brs(i,j)/re) ** 2 )  ! theta is in radians

!         get the dimensionless magnetic field: bm = b/b0 
!         and sine / cosine of dip angle

          bms(i,j) = sqrt ( ( 2. * cos(theta) ) ** 2 
     .                           + sin(theta)   ** 2 ) / br_norm ** 3

!          sindips(i,j) = 2.* cos (theta) / ( bms(i,j) * br_norm **3 )
!          cosdips(i,j) =     sin (theta) / ( bms(i,j) * br_norm **3 )

          call vector( brs(i,j),blon0,blats(i,j),
     .                 arg(i,j),athg(i,j),aphig(i,j) ) 


!         we define the dimensional grid points given by s = q * re

          s(i,j)  = qs(i,j)  * re 
          alts(i,j) = grs(i,j) - re   ! altitude above earth

!         grid in cartesian coordinates

          xs(i,j)    = brs(i,j) * cos( blats(i,j)*po180 ) *
     .                            sin( blon0     *po180 )
          ys(i,j)    = brs(i,j) * sin( blats(i,j)*po180 )
          zs(i,j)    = brs(i,j) * cos( blats(i,j)*po180 ) *
     .                            cos( blon0     *po180 )

        enddo

!       second half of field line

        do i = nz,nzh+2,-1

          x    = ( float(i) - ( float(nzh) + 1.) ) / float(nzh)
          xx   = x * sinh ( gams * qminn )
          qtmp = alog ( xx + sqrt ( xx**2 + 1. ) ) / gams

!         x    = -1. + 2. * ( float(i) - 1 ) / float(nz-1)
!         xx   = x/2. * ( sinh(gams*qminn) * (x+1) + 
!     .                   sinh(gams*qmins) * (x-1) )
!         qtmp = alog ( xx + sqrt ( xx**2 + 1. ) ) / gams
 
          r = re*qp_solve(qtmp,pvalue)
         
          qs(i,j)       = qtmp
          ps(i,j)       = pvalue
          brs(i,j)      = r
          br_norm       = r / re
          blats(i,j)  = asin ( qs(i,j) * br_norm ** 2 ) * rtod
          call btog ( brs(i,j),blon0,blats(i,j),
     .                grs(i,j),glats(i,j),glons(i,j) )

!         we need to calculate some quantities which depend on the usual
!         angular definition of theta

          theta = acos ( qs(i,j) * (brs(i,j)/re) ** 2 )  ! theta is in radians

!         get the dimensionless magnetic field: bm = b/b0 
!         and sine / cosine of dip angle

          bms(i,j) = sqrt ( ( 2. * cos(theta) ) ** 2 
     .                           + sin(theta)   ** 2 ) / br_norm ** 3

!          sindips(i,j) = 2.* cos (theta) / ( bms(i,j) * br_norm **3 )
!          cosdips(i,j) =     sin (theta) / ( bms(i,j) * br_norm **3 )

          call vector( brs(i,j),blon0,blats(i,j),
     .                 arg(i,j),athg(i,j),aphig(i,j) ) 

!         we define the dimensional grid points given by s = q * re

          s(i,j)  = qs(i,j)  * re 
          alts(i,j) = grs(i,j) - re   ! altitude above earth

!         grid in cartesian coordinates

          xs(i,j)    = brs(i,j) * cos( blats(i,j)*po180 ) *
     .                            sin( blon0     *po180 )
          ys(i,j)    = brs(i,j) * sin( blats(i,j)*po180 )
          zs(i,j)    = brs(i,j) * cos( blats(i,j)*po180 ) *
     .                            cos( blon0     *po180 )

        enddo

      enddo ! j loop (number of field lines)

!     following distances are in cm
     
      do j = 1,nf
        do i = 2,nz-1
          ds(i,j)   = ( s(i,j)   - s(i-1,j) ) * 1.e5
          d2s(i,j)  = ( s(i+1,j) - s(i-1,j) ) * 1.e5
          d22s(i,j) = .5 * d2s(i,j)
        enddo
        ds(1,j)    = ds(2,j)
        ds(nz,j)   = ds(nz-1,j)
        d2s(1,j)   = d2s(2,j)
        d2s(nz,j)  = d2s(nz-1,j)
        d22s(1,j)  = d22s(2,j)
        d22s(nz,j) = d22s(nz-1,j)
      enddo

!     calculate dels: grid length along field line using cartesian geometry
!     and convert to cm from km

      do j = 1,nf
        do i = 1,nz-1
          x2  = ( xs(i+1,j) - xs(i,j) ) ** 2
          y2  = ( ys(i+1,j) - ys(i,j) ) ** 2
          zz2 = ( zs(i+1,j) - zs(i,j) ) ** 2
          dels(i,j) = sqrt ( x2 + y2 + zz2 ) * 1.e5
        enddo
        dels(nz,j) = dels(nz-1,j) 
      enddo

!     now do gridding for p mesh

!     calculate new set of q's  (qp)    

      do j = 1,nf
        do i = 1,nz-1
          qp(i+1,j) = 0.5 * ( qs(i,j) + qs(i+1,j) )
        enddo
        qp(1,j)    = qs(1,j)  + 0.5 * ( qs(1,j) - qs(2,j) )
        qp(nzp1,j) = qs(nz,j) + 0.5 * ( qs(nz,j) - qs(nz-1,j) )
      enddo
      do i = 1,nzp1
        qp(i,nfp1) = qp(i,nf) + ( qp(i,nf) - qp(i,nf-1) )
      enddo

      do i = 1,nz
        do j = 1,nf-1
          pp(i,j+1) = 0.5 * ( ps(i,j) + ps(i,j+1) )
        enddo
        pp(i,1)    = ps(i,1)  + 0.5 * ( ps(i,1)  - ps(i,2)    ) 
        pp(i,nfp1) = ps(i,nf) + 0.5 * ( ps(i,nf) - ps(i,nf-1) )
      enddo
      do j = 1,nfp1
        pp(nzp1,j) = pp(nz,j) + 0.5 * ( pp(nz,j) - pp(nz-1,j) )
      enddo

!     now calculate the grid for nfp1 to 1 

      do j = nfp1,1,-1
        do i = 1,nzp1
          pvalue   = pp(i,j)
          if ( j .eq. nfp1 .and. i .ne. nzp1 ) then
            r    = brs(i,nf)
          elseif ( j .eq. nfp1 .and. i .eq. nzp1 ) then
            r    = brs(nz,nf)
          else
            r    = brp(i,j+1)
          endif
          qtmp = qp(i,j)
          r    = re*qp_solve(qtmp,pvalue)

          brp(i,j)    = r
          br_norm     = r / re
          blatp(i,j)  = asin ( qp(i,j) * br_norm ** 2 ) * rtod
          blonp(i,j)  = blon0

          call btog ( brp(i,j),blonp(i,j),blatp(i,j),
     .                grp(i,j),glatp(i,j),glonp(i,j) )

          altp(i,j)   = grp(i,j) - re

          delphi =  5.
          blon1  =  blon0 - delphi
          if ( blon1 .ge. 360. ) blon1 = blon1 - 360.
          blon2  =  blon0 + delphi
          if ( blon2 .ge. 360. ) blon2 = blon2 - 360.

!         grid in cartesian coordinates

          xp(i,j,1)   = brp(i,j) * cos( blatp(i,j) * po180 ) *
     .                             sin( blon1      * po180 ) 
          yp(i,j,1)   = brp(i,j) * sin( blatp(i,j) * po180 )
          zp(i,j,1)   = brp(i,j) * cos( blatp(i,j) * po180 ) *
     .                             cos( blon1      * po180 ) 

          xp(i,j,2)   = brp(i,j) * cos( blatp(i,j) * po180 ) *
     .                             sin( blon2      * po180 ) 
          yp(i,j,2)   = brp(i,j) * sin( blatp(i,j) * po180 )
          zp(i,j,2)   = brp(i,j) * cos( blatp(i,j) * po180 ) *
     .                             cos( blon2      * po180 ) 
        enddo
      enddo

!     calculate geometric parameters: cell area and volume

!     vol is cell volume

      call volume ( xp,yp,zp )

!     areap is cell face area in j-direction (p)
!     areas is cell face area in i-direction (s)

      call area ( xp,yp,zp )

!     xdels is line distance in i-direction (s)
!     xdelp is line distance in j-direction (p)
!     xdelh is line distance in k-direction (phi)

      call line ( xp,yp,zp )

!     normal: calculates normal to cell face in s-direction
!     vnormal: calcuates e x b direction

      call  normal ( xp,yp,zp )
      call vnormal ( qs,blon0,brs )

      end

*******************************************
*******************************************

!            gtob

*******************************************
*******************************************

      subroutine gtob(brad,blond,blatd,grad,glatd,glond)

!     conversion from geographic to 
!     offset centered dipole coordinates
!     brad: radius in the offset dipole system
!     grad: radius in the geocentric system
!     (g joyce june 1998)

      include 'param-1.00.inc'
      include 'com-1.00.mod.inc'

!     coordinates of dipole in km relative to center of earth

!      x0 = -392. 
!      y0 =  258.
!      z0 =  179.

!     convert to radians

      glat = glatd * po180
      glon = glond * po180

!     rotate in longitude

      xg = grad * cos ( glat ) * cos ( glon - plon )
      yg = grad * cos ( glat ) * sin ( glon - plon )
      zg = grad * sin ( glat )

!     rotate in latitude

      xmm = xg * sin ( plat ) - zg * cos ( plat )
      ymm = yg
      zmm = xg * cos ( plat ) + zg * sin ( plat )

!     calculate offset in geographic polar (p. 248 redbook) 

      d = sqrt(x0**2 + y0**2 + z0**2)
      thd = acos( z0/d )
      phd = atan2( y0,x0 )

!     change to colatitude

      cotn = pie/2. - plat
      phin = plon

!     calculate offset in cartesian rotated system

      thp = acos(cos(cotn)*cos(thd) 
     .     + sin(cotn)*sin(thd)*cos(phd-phin))
      cosphp = ( cos(cotn)*cos(thp) - cos(thd) ) 
     .         / (sin(cotn) * sin(thp))
      sinphp = sin(thd) * sin(phd-phin) / sin(thp)

      x0p = d*sin(thp)*cosphp
      y0p = d*sin(thp)*sinphp
      z0p = d*cos(thp)


!     geomagnetic latitude and longitude in degrees

      xm = xmm - x0p
      ym = ymm - y0p
      zm = zmm - z0p

!     magnetic lat and long converted to degrees

      brad  = sqrt ( xm ** 2 + ym ** 2 + zm ** 2 )
      blatd = asin ( zm / brad ) * rtod
      blond = ( atan2 ( ym/brad,xm/brad ) ) * rtod
      if (blond .lt. -0.01) blond = blond+360.

      return
      end

*******************************************
*******************************************

!            btog

*******************************************
*******************************************


      subroutine btog ( brad,blond,blatd,grad,glatd,glond )

!     conversion from centered dipole to geographic coordinates 
!     plat,plon =  coords of north magnetic pole (in param2d.9.14e.inc)
!     the geomagnetic longitude is measured east from the 
!     geomagnetic prime meridian at 291 degrees east geographic.
!     brad: radius in the offset dipole frame
!     grad: radius in the geocentric frame
!     (g joyce june 1998)

      include 'param-1.00.inc'
      include 'com-1.00.mod.inc'

      real brad,blond,blatd,grad,glatd,glond

!     coordinates of dipole in km relatvie to center of earth

!      x0 = -392. 
!      y0 =  258.
!      z0 =  179.

!     convert magnetic lat and long to radians

      blonr = blond * po180 
      blatr = blatd * po180 

!     position of point in geomagnetic coords
!     first get xm ym zm in the eccentric dipole system

      xm = brad *cos ( blatr ) * cos ( blonr )
      ym = brad *cos ( blatr ) * sin ( blonr )
      zm = brad *sin ( blatr )

!     calculate offset in geographic polar (p. 248 redbook) 

      d = sqrt(x0**2 + y0**2 + z0**2)
      thd = acos( z0/d )
      phd = atan2( y0,x0 )

!     change to colatitude

      cotn = pie/2. - plat
      phin = plon

!     calculate offset in cartesian rotated system

      thp = acos(cos(cotn)*cos(thd) 
     .      + sin(cotn)*sin(thd)*cos(phd-phin))
      cosphp = ( cos(cotn)*cos(thp) - cos(thd) ) 
     .         / (sin(cotn) * sin(thp))
      sinphp = sin(thd) * sin(phd-phin) / sin(thp)

      x0p = d*sin(thp)*cosphp
      y0p = d*sin(thp)*sinphp
      z0p = d*cos(thp)

!     next shift to the tilted dipole

      xmm = xm + x0p
      ymm = ym + y0p
      zmm = zm + z0p

!     r is invariant under the rotations of the tilted dipole

      grad = sqrt ( xmm ** 2 + ymm ** 2 + zmm ** 2 )

!     rotate coords in north-south direction

      xg =  xmm * sin ( plat ) + zmm * cos ( plat )
      yg =  ymm
      zg = -xmm * cos ( plat ) + zmm * sin ( plat )

!     geographic latitude and longitude converted to degrees

      glatd = asin ( zg / grad ) * rtod
      glond = ( plon + atan2 ( yg/grad,xg/grad ) ) * rtod

      if (glond .ge. 360) glond = glond-  360.

      return
      end

*******************************************
*******************************************

!             volume

*******************************************
*******************************************

        subroutine volume(x,y,z)

!       calculate cell volume
!       break each cell into
!       twelve tetrahedrons and use the formula: 
!           V = (1/6) a . ( b x c )
!       where
!           a: vector from A to B
!           b: vector from A to C
!           c: vector from A to D
!       and node A is the 'midpoint' coordinate

        include 'param-1.00.inc'
        include 'com-1.00.mod.inc'

        real x(nzp1,nfp1,2),y(nzp1,nfp1,2),z(nzp1,nfp1,2)
        real voli(nz,nf),volj(nz,nf),volk(nz,nf)

!       volume from sidei

        do k = 1,1
          do j = 1,nf
            do i = 1,nz

              xmid  = 0.5 * ( x(i,j,k) + x(i+1,j+1,k+1) )
              ymid  = 0.5 * ( y(i,j,k) + y(i+1,j+1,k+1) )
              zmid  = 0.5 * ( z(i,j,k) + z(i+1,j+1,k+1) )
              
              ax1 = x(i,j,k) - xmid
              ay1 = y(i,j,k) - ymid
              az1 = z(i,j,k) - zmid

              bx1 = x(i,j+1,k) - xmid
              by1 = y(i,j+1,k) - ymid
              bz1 = z(i,j+1,k) - zmid

              cx1 = x(i,j,k+1) - xmid
              cy1 = y(i,j,k+1) - ymid
              cz1 = z(i,j,k+1) - zmid

              dx1 =    by1 * cz1 - bz1 * cy1
              dy1 = -( bx1 * cz1 - bz1 * cx1 )
              dz1 =    bx1 * cy1 - by1 * cx1

              v1 = 0.166667 * 
     .             abs(( ax1 * dx1 + ay1 * dy1 + az1 * dz1 ))

              ax2 = x(i,j+1,k+1) - xmid
              ay2 = y(i,j+1,k+1) - ymid
              az2 = z(i,j+1,k+1) - zmid

              v2 = 0.166667 * 
     .             abs(( ax2 * dx1 + ay2 * dy1 + az2 * dz1 ))

              ax1 = x(i+1,j,k) - xmid
              ay1 = y(i+1,j,k) - ymid
              az1 = z(i+1,j,k) - zmid

              bx1 = x(i+1,j+1,k) - xmid
              by1 = y(i+1,j+1,k) - ymid
              bz1 = z(i+1,j+1,k) - zmid

              cx1 = x(i+1,j,k+1) - xmid
              cy1 = y(i+1,j,k+1) - ymid
              cz1 = z(i+1,j,k+1) - zmid

              dx1 =    by1 * cz1 - bz1 * cy1
              dy1 = -( bx1 * cz1 - bz1 * cx1 )
              dz1 =    bx1 * cy1 - by1 * cx1

              v3 = 0.166667 * 
     .             abs(( ax1 * dx1 + ay1 * dy1 + az1 * dz1 ))
       
              ax2 = x(i+1,j+1,k+1) - xmid
              ay2 = y(i+1,j+1,k+1) - ymid
              az2 = z(i+1,j+1,k+1) - zmid

              v4 = 0.166667 * 
     .             abs(( ax2 * dx1 + ay2 * dy1 + az2 * dz1 ))

              voli(i,j) = v1 + v2 + v3 + v4
               
            enddo
          enddo
        enddo

!       volume from sidej

        do k = 1,1
          do j = 1,nf
            do i = 1,nz

              xmid = 0.5 * ( x(i,j,k) + x(i+1,j+1,k+1) )
              ymid = 0.5 * ( y(i,j,k) + y(i+1,j+1,k+1) )
              zmid = 0.5 * ( z(i,j,k) + z(i+1,j+1,k+1) )

              ax1 = x(i,j,k) - xmid
              ay1 = y(i,j,k) - ymid
              az1 = z(i,j,k) - zmid

              bx1 = x(i+1,j,k) - xmid
              by1 = y(i+1,j,k) - ymid
              bz1 = z(i+1,j,k) - zmid

              cx1 = x(i,j,k+1) - xmid
              cy1 = y(i,j,k+1) - ymid
              cz1 = z(i,j,k+1) - zmid

              dx1 =    by1 * cz1 - bz1 * cy1
              dy1 = -( bx1 * cz1 - bz1 * cx1 )
              dz1 =    bx1 * cy1 - by1 * cx1

              v1 = 0.166667 * 
     .             abs(( ax1 * dx1 + ay1 * dy1 + az1 * dz1 ))

              ax2 = x(i+1,j,k+1) - xmid
              ay2 = y(i+1,j,k+1) - ymid
              az2 = z(i+1,j,k+1) - zmid

              v2 = 0.166667 * 
     .             abs(( ax2 * dx1 + ay2 * dy1 + az2 * dz1 ))

              ax1 = x(i,j+1,k) - xmid
              ay1 = y(i,j+1,k) - ymid
              az1 = z(i,j+1,k) - zmid

              bx1 = x(i+1,j+1,k) - xmid
              by1 = y(i+1,j+1,k) - ymid
              bz1 = z(i+1,j+1,k) - zmid

              cx1 = x(i,j+1,k+1) - xmid
              cy1 = y(i,j+1,k+1) - ymid
              cz1 = z(i,j+1,k+1) - zmid

              dx1 =    by1 * cz1 - bz1 * cy1
              dy1 = -( bx1 * cz1 - bz1 * cx1 )
              dz1 =    bx1 * cy1 - by1 * cx1

              v3 = 0.166667 * 
     .             abs(( ax1 * dx1 + ay1 * dy1 + az1 * dz1 ))
       
              ax2 = x(i+1,j+1,k+1) - xmid
              ay2 = y(i+1,j+1,k+1) - ymid
              az2 = z(i+1,j+1,k+1) - zmid

              v4 = 0.166667 * 
     .             abs(( ax2 * dx1 + ay2 * dy1 + az2 * dz1 ))

              volj(i,j) = v1 + v2 + v3 + v4 
 
            enddo
          enddo
        enddo

!       volume from sidek 

        do k = 1,1
          do j = 1,nf
            do i = 1,nz

              xmid = 0.5 * ( x(i,j,k) + x(i+1,j+1,k+1) )
              ymid = 0.5 * ( y(i,j,k) + y(i+1,j+1,k+1) )
              zmid = 0.5 * ( z(i,j,k) + z(i+1,j+1,k+1) )

              ax1 = x(i,j,k) - xmid
              ay1 = y(i,j,k) - ymid
              az1 = z(i,j,k) - zmid

              bx1 = x(i+1,j,k) - xmid
              by1 = y(i+1,j,k) - ymid
              bz1 = z(i+1,j,k) - zmid

              cx1 = x(i,j+1,k) - xmid
              cy1 = y(i,j+1,k) - ymid
              cz1 = z(i,j+1,k) - zmid

              dx1 =    by1 * cz1 - bz1 * cy1
              dy1 = -( bx1 * cz1 - bz1 * cx1 )
              dz1 =    bx1 * cy1 - by1 * cx1

              v1 = 0.166667 * 
     .             abs(( ax1 * dx1 + ay1 * dy1 + az1 * dz1 ))

              ax2 = x(i+1,j+1,k) - xmid
              ay2 = y(i+1,j+1,k) - ymid
              az2 = z(i+1,j+1,k) - zmid

              v2 = 0.166667 * 
     .             abs(( ax2 * dx1 + ay2 * dy1 + az2 * dz1 ))

              ax1 = x(i,j,k+1) - xmid
              ay1 = y(i,j,k+1) - ymid
              az1 = z(i,j,k+1) - zmid

              bx1 = x(i+1,j,k+1) - xmid
              by1 = y(i+1,j,k+1) - ymid
              bz1 = z(i+1,j,k+1) - zmid

              cx1 = x(i,j+1,k+1) - xmid
              cy1 = y(i,j+1,k+1) - ymid
              cz1 = z(i,j+1,k+1) - zmid

              dx1 =    by1 * cz1 - bz1 * cy1
              dy1 = -( bx1 * cz1 - bz1 * cx1 )
              dz1 =    bx1 * cy1 - by1 * cx1

              v3 = 0.166667 * 
     .             abs(( ax1 * dx1 + ay1 * dy1 + az1 * dz1 ))
       
              ax2 = x(i+1,j+1,k+1) - xmid
              ay2 = y(i+1,j+1,k+1) - ymid
              az2 = z(i+1,j+1,k+1) - zmid

              v4 = 0.166667 * 
     .             abs(( ax2 * dx1 + ay2 * dy1 + az2 * dz1 ))

              volk(i,j) = v1 + v2 + v3 + v4 

            enddo
          enddo
        enddo

        do j = 1,nf
          do i = 1,nz
            vol(i,j) = 1.e15 *
     .                 ( voli(i,j) + volj(i,j) + volk(i,j) )

          enddo
        enddo

        return
        end

*******************************************
*******************************************

!            area

*******************************************
*******************************************


        subroutine area(x,y,z)

!       calculate areas of cell sides
!       break each quadrilateral side into
!       two triangles and use the formula: 
!           A = (1/2)|a x b|
!       where
!           a: vector from A to B
!           b: vector from A to C

        include 'param-1.00.inc'
        include 'com-1.00.mod.inc'

        real x(nzp1,nfp1,2),y(nzp1,nfp1,2),z(nzp1,nfp1,2)

!       sidei (s-direction)

        do k = 1,1
          do j = 1,nf
            do i = 1,nzp1

              ax1 = x(i,j+1,k) - x(i,j,k)
              ay1 = y(i,j+1,k) - y(i,j,k)
              az1 = z(i,j+1,k) - z(i,j,k)

              bx1 = x(i,j,k+1) - x(i,j,k)
              by1 = y(i,j,k+1) - y(i,j,k)
              bz1 = z(i,j,k+1) - z(i,j,k)

              cx1 =    ay1 * bz1 - az1 * by1
              cy1 = -( ax1 * bz1 - az1 * bx1 )
              cz1 =    ax1 * by1 - ay1 * bx1

              a1 = 0.5 * sqrt ( cx1*cx1 + cy1*cy1 + cz1*cz1 )

              ax2 = x(i,j+1,k) - x(i,j+1,k+1)
              ay2 = y(i,j+1,k) - y(i,j+1,k+1)
              az2 = z(i,j+1,k) - z(i,j+1,k+1)

              bx2 = x(i,j,k+1) - x(i,j+1,k+1)
              by2 = y(i,j,k+1) - y(i,j+1,k+1)
              bz2 = z(i,j,k+1) - z(i,j+1,k+1)

              cx2 =    ay2 * bz2 - az2 * by2
              cy2 = -( ax2 * bz2 - az2 * bx2 )
              cz2 =    ax2 * by2 - ay2 * bx2

              a2 = 0.5 * sqrt ( cx2*cx2 + cy2*cy2 + cz2*cz2 )

              areas(i,j) = ( a1 + a2 ) * 1.e10

            enddo
          enddo
        enddo

!       sidej (p-direction)

        do k = 1,1
          do j = 1,nfp1
            do i = 1,nz

              ax1 = x(i+1,j,k) - x(i,j,k)
              ay1 = y(i+1,j,k) - y(i,j,k)
              az1 = z(i+1,j,k) - z(i,j,k)

              bx1 = x(i,j,k+1) - x(i,j,k)
              by1 = y(i,j,k+1) - y(i,j,k)
              bz1 = z(i,j,k+1) - z(i,j,k)

              cx1 =    ay1 * bz1 - az1 * by1
              cy1 = -( ax1 * bz1 - az1 * bx1 )
              cz1 =    ax1 * by1 - ay1 * bx1

              a1 = 0.5 * sqrt ( cx1*cx1 + cy1*cy1 + cz1*cz1 )
       
              ax2 = x(i+1,j,k) - x(i+1,j,k+1)
              ay2 = y(i+1,j,k) - y(i+1,j,k+1)
              az2 = z(i+1,j,k) - z(i+1,j,k+1)

              bx2 = x(i,j,k+1) - x(i+1,j,k+1)
              by2 = y(i,j,k+1) - y(i+1,j,k+1)
              bz2 = z(i,j,k+1) - z(i+1,j,k+1)

              cx2 =    ay2 * bz2 - az2 * by2
              cy2 = -( ax2 * bz2 - az2 * bx2 )
              cz2 =    ax2 * by2 - ay2 * bx2

              a2 = 0.5 * sqrt ( cx2*cx2 + cy2*cy2 + cz2*cz2 )
       
              areap(i,j) = ( a1 + a2 ) * 1.e10

            enddo
          enddo
        enddo

        return
        end

*******************************************
*******************************************

!            line

*******************************************
*******************************************


        subroutine line(x,y,z)

!       calculate length of cell sides

        include 'param-1.00.inc'
        include 'com-1.00.mod.inc'

        real x(nzp1,nfp1,2),y(nzp1,nfp1,2),z(nzp1,nfp1,2)
        real xdelh(nzp1,nfp1)

!       xdels (s-direction)

        do k = 1,2
          do j = 1,nfp1
            do i = 1,nz

              ax1 = x(i+1,j,k) - x(i,j,k)
              ay1 = y(i+1,j,k) - y(i,j,k)
              az1 = z(i+1,j,k) - z(i,j,k)

              xdels(i,j,k) = sqrt ( ax1*ax1 + ay1*ay1 + az1*az1 )
     .                       * 1.e5

            enddo
          enddo
        enddo

!       xdelp (p-direction)

        do k = 1,2
          do j = 1,nf
            do i = 1,nzp1

              ax1 = x(i,j+1,k) - x(i,j,k)
              ay1 = y(i,j+1,k) - y(i,j,k)
              az1 = z(i,j+1,k) - z(i,j,k)

              xdelp(i,j,k) = sqrt ( ax1*ax1 + ay1*ay1 + az1*az1 )
     .                       * 1.e5

            enddo
          enddo
        enddo

!       xdelh (phi-direction)

          do j = 1,nfp1
            do i = 1,nzp1

              ax1 = x(i,j,2) - x(i,j,1)
              ay1 = y(i,j,2) - y(i,j,1)
              az1 = z(i,j,2) - z(i,j,1)

              xdelh(i,j) = sqrt ( ax1*ax1 + ay1*ay1 + az1*az1 )
     .                     * 1.e5

            enddo
          enddo

        return
        end

*******************************************
*******************************************

!            normal

*******************************************
*******************************************


        subroutine normal ( x,y,z )

!       calculate unit normal direction to cell face
!       normal: c = a x b / |a x b|

        include 'param-1.00.inc'
        include 'com-1.00.mod.inc'

        real x(nzp1,nfp1,2),y(nzp1,nfp1,2),z(nzp1,nfp1,2)

!       norms (normal to cell face in s-direction)

        do j = 1,nf
          do i = 1,nzp1
            
            ax1 = x(i,j+1,2) - x(i,j,1)
            ay1 = y(i,j+1,2) - y(i,j,1)
            az1 = z(i,j+1,2) - z(i,j,1)

            bx1 = x(i,j,2) - x(i,j+1,1)
            by1 = y(i,j,2) - y(i,j+1,1)
            bz1 = z(i,j,2) - z(i,j+1,1)

            cx1 =   ay1 * bz1 - az1 * by1
            cy1 = -(ax1 * bz1 - az1 * bx1)
            cz1 =   ax1 * by1 - ay1 * bx1

            ca1 = sqrt( cx1*cx1 + cy1*cy1 + cz1*cz1 )

            xnorms(i,j) = cx1/ca1
            ynorms(i,j) = cy1/ca1
            znorms(i,j) = cz1/ca1

          enddo
        enddo

!       normp (normal to cell face in p-direction)

        do j = 1,nfp1
          do i = 1,nz
            
            ax1 = x(i,j,2) - x(i+1,j,1)
            ay1 = y(i,j,2) - y(i+1,j,1)
            az1 = z(i,j,2) - z(i+1,j,1)

            bx1 = x(i+1,j,2) - x(i,j,1)
            by1 = y(i+1,j,2) - y(i,j,1)
            bz1 = z(i+1,j,2) - z(i,j,1)

            cx1 =   ay1 * bz1 - az1 * by1
            cy1 = -(ax1 * bz1 - az1 * bx1)
            cz1 =   ax1 * by1 - ay1 * bx1

            ca1 = sqrt( cx1*cx1 + cy1*cy1 + cz1*cz1 )

            xnormp(i,j) = cx1/ca1
            ynormp(i,j) = cy1/ca1
            znormp(i,j) = cz1/ca1
          enddo
        enddo

!       norml (normal to cell face in phi-direction: longitude)

!        do j = 1,nf
!          do i = 1,nzp1
            
!            ax1 = x(i+1,j+1,1) - x(i,j,1)
!            ay1 = y(i+1,j+1,1) - y(i,j,1)
!            az1 = z(i+1,j+1,1) - z(i,j,1)

!            bx1 = x(i,j+1,1) - x(i+1,j,1)
!            by1 = y(i,j+1,1) - y(i+1,j,1)
!            bz1 = z(i,j+1,1) - z(i+1,j,1)

!            cx1 =   ay1 * bz1 - az1 * by1
!            cy1 = -(ax1 * bz1 - az1 * bx1)
!            cz1 =   ax1 * by1 - ay1 * bx1

!            ca1 = sqrt( cx1*cx1 + cy1*cy1 + cz1*cz1 )

!            xnorml(i,j) = cx1/ca1
!            ynorml(i,j) = cy1/ca1
!            znorml(i,j) = cz1/ca1

!          enddo
!        enddo

        return
        end

*******************************************
*******************************************

!            vnormal

*******************************************
*******************************************


      subroutine vnormal ( qs,blon0,brs )

!     calculate e x b velocity direction
!     change p but keep q constant

      include 'param-1.00.inc'
      include 'com-1.00.mod.inc'

      real qs(nz,nf),qss(nzp1,nf),pss(nzp1,nf)
      real brs1(nzp1,nf),blats1(nzp1,nf)
      real brs2(nzp1,nf),blats2(nzp1,nf)
      real xs1(nzp1,nf),ys1(nzp1,nf),zs1(nzp1,nf)
      real xs2(nzp1,nf),ys2(nzp1,nf),zs2(nzp1,nf)
      real brs(nz,nf)
 
!      open (unit=12, file='xs1.dat')
!      open (unit=13, file='ys1.dat')
!      open (unit=15, file='xs2.dat')
!      open (unit=16, file='ys2.dat')

      delps = .001

      do j = 1,nf
        do i = 1,nz-1
          qss(i+1,j) = 0.5 * ( qs(i,j) + qs(i+1,j) )
        enddo
        qss(1,j)    = qs(1,j)  + 0.5 * ( qs(1,j) - qs(2,j) )
        qss(nzp1,j) = qs(nz,j) + 0.5 * ( qs(nz,j) - qs(nz-1,j) )
      enddo

      do j = 1,nf
        do i = 1,nz
          pss(i+1,j) = ps(i,j)
        enddo
      enddo
      do j = 1,nf
        pss(1,j) = pss(2,j)
      enddo

      do j = 1,nf
        do i = 1,nz+1

!         grid 1 (at pss(i,j))

          qtmp   = qss(i,j)
          pvalue = pss(i,j)
          r = re*qp_solve(qtmp,pvalue)
         
          brs1(i,j)      = r
          br_norm        = r / re
          blats1(i,j)    = asin ( qss(i,j) * br_norm ** 2 ) * rtod

          xs1(i,j)    = brs1(i,j) * cos( blats1(i,j)*po180 ) *
     .                              sin( blon0      * po180 ) 
          ys1(i,j)    = brs1(i,j) * sin( blats1(i,j)*po180 )
          zs1(i,j)    = brs1(i,j) * cos( blats1(i,j) * po180 ) *
     .                              cos( blon0      * po180 ) 

!         grid 2 (at pss(i,j)+delps)

          qtmp   = qss(i,j)
          pvalue = pss(i,j) + delps
          r = re*qp_solve(qtmp,pvalue)
         
          brs2(i,j)      = r
          br_norm        = r / re
          blats2(i,j)    = asin ( qss(i,j) * br_norm ** 2 ) * rtod

!         grid in cartesian coordinates

          xs2(i,j)    = brs2(i,j) * cos( blats2(i,j)*po180 ) *
     .                              sin( blon0      * po180 ) 
          ys2(i,j)    = brs2(i,j) * sin( blats2(i,j)*po180 )
          zs2(i,j)    = brs2(i,j) * cos( blats2(i,j) * po180 ) *
     .                              cos( blon0       * po180 ) 

        enddo
      enddo

!     direction of e x b velocity

      do j = 1,nf
        do i = 1,nzp1
          ax1 = xs2(i,j) - xs1(i,j)
          ay1 = ys2(i,j) - ys1(i,j)
          az1 = zs2(i,j) - zs1(i,j)
          a1  = sqrt ( ax1*ax1 + ay1*ay1 + az1*az1 )
          vnx(i,j) = ax1 / a1
          vny(i,j) = ay1 / a1
          vnz(i,j) = az1 / a1
        enddo
      enddo

      do i = 1,nzp1
        vnx(i,nfp1) = vnx(i,nf)
        vny(i,nfp1) = vny(i,nf)
        vnz(i,nfp1) = vnz(i,nf)
      enddo
      

      return
      end


********************************************************

        function qp_solve(q,p)

        real qp_solve,q,p
        real term1,term2

!       MS: Functionally unnecessary, but kind of neat.  To get r from
!       (q,p) one needs to find the root of a fourth-order polynomial.
!       The code did this with a standard root-finding method,
!       but, of course, there is a closed-form solution.  Since the
!       polynomial has no second or third order terms, the answer
!       isn't completely ugly, and the code below gives the exact
!       result.  Note the special case for q = 0, i.e. points on
!       the magnetic equator.
  
!        if (p .eq. 0.) then
!          print *,'p = 0 in qp_solve'
!          stop
!        endif

!        if (q .eq. 0.) then
!          qp_solve = p
!        else
!        term1 = ( (sqrt(27.) + sqrt(27.+256.*q**2*p**4)) / 
!     .             (16.*abs(q)*p**2) ) ** (1./3.)
!        term2 = 2./(abs(q)*sqrt(3.)) * (term1 - 1./term1)
!        qp_solve = 0.5*sqrt(term2)*(sqrt(2./(p*q**2*sqrt(term2**3)) - 
!     .                       1.) - 1.)
!      endif


! MS: The formula is actually my third attempt.  The first had huge
! cancellation problems when q was small, the second had smaller
! problems when term0 was large. This should be well-behaved
! everywhere.  Massive algebra was involved along the way.

        term0 = 256./27.*q**2*p**4
        term1 = ( 1. + sqrt(1. + term0) )**(2./3.)
        term2 = term0**(1./3.)
        term3 = 0.5 * ( (term1**2 + term1*term2 + 
     .                term2**2)/term1 )**(3./2.)
        qp_solve = p * (4.*term3) / (1.+term3) / 
     .                 ( 1.+sqrt(2.*term3-1.) )


        end 

!    IGRF_SUB.F
!    contributed by Paul Bernhardt (NRL)

       subroutine igrf_sub(phi0)

       include 'param-1.00.inc'
       include 'com-1.00.mod.inc'

       phi = phi0 * pie / 180.

       x0  = -353.8674146847993 -
     .        80.0112901196596*Cos(phi) +
     .        165.9453553374371*Cos(2*phi) + 
     .        96.90976911496519*Cos(3*phi) - 
     .        41.234237919645906*Cos(4*phi) - 
     .        5.647625148538237*Cos(5*phi) + 
     .        29.81898130515779*Cos(6*phi) + 
     .        27.310156089558348*Cos(7*phi) + 
     .        7.017930973481185*Cos(8*phi) - 
     .        7.511889133487608*Cos(9*phi) - 
     .        6.720946815524082*Cos(10*phi) - 
     .        0.5967604884638268*Cos(11*phi) + 
     .        3.9340250001598633*Cos(12*phi) + 
     .        3.232375199053169*Cos(13*phi) + 
     .        0.1418342157838855*Cos(14*phi) - 
     .        1.7399223106946722*Cos(15*phi) - 
     .        1.2396197007828484*Cos(16*phi) + 
     .        0.24319862206273882*Cos(17*phi) + 
     .        1.0043588925588962*Cos(18*phi) + 
     .        0.6050359806604826*Cos(19*phi) - 
     .        0.2175679881607454*Cos(20*phi) - 
     .        0.5970638247030868*Cos(21*phi) - 
     .        0.3095886741276035*Cos(22*phi) + 
     .        0.24106351981877516*Cos(23*phi) - 
     .        39.488705071918595*Sin(phi) - 
     .        35.12750391544765*Sin(2*phi) - 
     .        52.61217403191678*Sin(3*phi) + 
     .        35.372208218215505*Sin(4*phi) + 
     .        62.517949955768565*Sin(5*phi) + 
     .        33.70515829546757*Sin(6*phi) - 
     .        6.367067330825964*Sin(7*phi) - 
     .        16.84168815021257*Sin(8*phi) - 
     .        8.77073458124951*Sin(9*phi) + 
     .        2.4509417362535326*Sin(10*phi) + 
     .        5.964643505586716*Sin(11*phi) + 
     .        2.7430884455100206*Sin(12*phi) - 
     .        1.469044512466587*Sin(13*phi) - 
     .        2.58647152886297*Sin(14*phi) - 
     .        0.9057773601920639*Sin(15*phi) + 
     .        0.9500051070510113*Sin(16*phi) + 
     .        1.2595308671288294*Sin(17*phi) + 
     .        0.30670958561561906*Sin(18*phi) - 
     .        0.6041955011153372*Sin(19*phi) - 
     .        0.669980882605414*Sin(20*phi) - 
     .        0.08600912796198656*Sin(21*phi) + 
     .        0.4475307568863845*Sin(22*phi) + 
     .        0.4538479001719651*Sin(23*phi)



       y0 = 357.71572803711445 + 
     .      179.10000315892776*Cos(phi) + 
     .      0.48625012094252057*Cos(2*phi) - 
     .      18.222451520290683*Cos(3*phi) - 
     .      74.19071886155683*Cos(4*phi) - 
     .      36.523066435387676*Cos(5*phi) - 
     .      1.3391741086911595*Cos(6*phi) + 
     .      21.128815605008093*Cos(7*phi) + 
     .      12.30309152643703*Cos(8*phi) + 
     .      0.7072007210810681*Cos(9*phi) - 
     .      4.253716697410791*Cos(10*phi) - 
     .      1.9043150223068182*Cos(11*phi) + 
     .      1.5369582665526593*Cos(12*phi) + 
     .      2.322273352128386*Cos(13*phi) + 
     .      0.8660859526717453*Cos(14*phi) - 
     .      0.6517850740970705*Cos(15*phi) - 
     .      0.8233191201701169*Cos(16*phi) - 
     .      0.06587490166894897*Cos(17*phi) + 
     .      0.5590511561227545*Cos(18*phi) + 
     .      0.4857427492770479*Cos(19*phi) + 
     .      0.003220351289183923*Cos(20*phi) - 
     .      0.31210096865413983*Cos(21*phi) - 
     .      0.19887427994494658*Cos(22*phi) + 
     .      0.13672910587356485*Cos(23*phi) + 
     .      416.65605532707093*Sin(phi) + 
     .      97.64092937463448*Sin(2*phi) + 
     .      89.90576641057676*Sin(3*phi) + 
     .      19.563789031651595*Sin(4*phi) + 
     .      46.61385007838707*Sin(5*phi) + 
     .      46.23624804545088*Sin(6*phi) + 
     .      19.328984217463248*Sin(7*phi) + 
     .      2.4598324280936414*Sin(8*phi) - 
     .      3.4043388007500432*Sin(9*phi) + 
     .      2.111663037458501*Sin(10*phi) + 
     .      4.723429973505645*Sin(11*phi) + 
     .      3.7939791153184026*Sin(12*phi) + 
     .      0.783837405015646*Sin(13*phi) - 
     .      0.8951808216876126*Sin(14*phi) - 
     .      0.5853291438883146*Sin(15*phi) + 
     .      0.49374732010654354*Sin(16*phi) + 
     .      0.9299824387332063*Sin(17*phi) + 
     .      0.48650275322307124*Sin(18*phi) - 
     .      0.15694935678709498*Sin(19*phi) - 
     .      0.36063573138052973*Sin(20*phi) - 
     .      0.08389592287106021*Sin(21*phi) + 
     .      0.25573941044924925*Sin(22*phi) + 
     .      0.2852013162673946*Sin(23*phi)

       z0 = 250.43165707492653 + 
     .      340.66122483264905*Cos(phi) + 
     .      69.95634836232071*Cos(2*phi) + 
     .      44.71168620310048*Cos(3*phi) + 
     .      5.572752558766106*Cos(4*phi) - 
     .      1.4339461502279616*Cos(5*phi) - 
     .      3.217930812956232*Cos(6*phi) - 
     .      1.7751137903566045*Cos(7*phi) - 
     .      1.5857495900135135*Cos(8*phi) - 
     .      0.7461261912547706*Cos(9*phi) + 
     .      0.14563859321432496*Cos(10*phi) + 
     .      0.3315895855881013*Cos(11*phi) + 
     .      0.11379367536280918*Cos(12*phi) - 
     .      0.15395572288202422*Cos(13*phi) - 
     .      0.19461671510187845*Cos(14*phi) - 
     .      0.06361599694133094*Cos(15*phi) + 
     .      0.055766203925038796*Cos(16*phi) + 
     .      0.06413644582150604*Cos(17*phi) + 
     .      0.0016750432363388439*Cos(18*phi) - 
     .      0.04182349014199574*Cos(19*phi) - 
     .      0.031182905113287293*Cos(20*phi) + 
     .      0.001753974707015248*Cos(21*phi) + 
     .      0.013973913699992777*Cos(22*phi) + 
     .      0.001281301856688131*Cos(23*phi) - 
     .      30.176529635127274*Sin(phi) + 
     .      3.4340205951414946*Sin(2*phi) - 
     .      6.8123984277392795*Sin(3*phi) - 
     .      11.668750059562694*Sin(4*phi) - 
     .      0.5810580595917775*Sin(5*phi) - 
     .      1.5870410791195877*Sin(6*phi) - 
     .      1.5244529990716116*Sin(7*phi) - 
     .      0.559954677080775*Sin(8*phi) + 
     .      0.5194385178951049*Sin(9*phi) + 
     .      0.415495066076052*Sin(10*phi) + 
     .      0.04521386570333395*Sin(11*phi) - 
     .      0.24535417687477926*Sin(12*phi) - 
     .      0.18203325048713365*Sin(13*phi) + 
     .      0.018878146825439902*Sin(14*phi) + 
     .      0.12360222004213386*Sin(15*phi) + 
     .      0.07257334462042853*Sin(16*phi) - 
     .      0.022737330338338836*Sin(17*phi) - 
     .      0.05776054926047047*Sin(18*phi) - 
     .      0.023788034168565193*Sin(19*phi) + 
     .      0.018325492914075404*Sin(20*phi) + 
     .      0.02376235757596904*Sin(21*phi) + 
     .      0.002628561346701197*Sin(22*phi) - 
     .      0.00943276960175865*Sin(23*phi)

       bb0 = 29229.88363541812 - 
     .      1115.0281463875144*Cos(phi) + 
     .      211.9965588481871*Cos(2*phi) - 
     .      145.87863848381699*Cos(3*phi) - 
     .      185.94991131071188*Cos(4*phi) - 
     .      94.76388090038132*Cos(5*phi) - 
     .      21.781240606511467*Cos(6*phi) + 
     .      10.608481026851262*Cos(7*phi) + 
     .      26.445003680471622*Cos(8*phi) + 
     .      9.033803335738316*Cos(9*phi) - 
     .      5.121597505263254*Cos(10*phi) - 
     .      8.699077903900074*Cos(11*phi) - 
     .      2.3758078601113084*Cos(12*phi) + 
     .      3.003659043878612*Cos(13*phi) + 
     .      3.3447365817081867*Cos(14*phi) + 
     .      0.46321746934859276*Cos(15*phi) - 
     .      1.6640994428565439*Cos(16*phi) - 
     .      1.3856408382560013*Cos(17*phi) + 
     .      0.06865975473204816*Cos(18*phi) + 
     .      0.8650543900984624*Cos(19*phi) + 
     .      0.4841643209453366*Cos(20*phi) - 
     .      0.22968004188392457*Cos(21*phi) - 
     .      0.38869426768468657*Cos(22*phi) - 
     .      0.00489072283424378*Cos(23*phi) - 
     .      343.1161166643772*Sin(phi) - 
     .      214.46870730780248*Sin(2*phi) - 
     .      218.4373355152967*Sin(3*phi) - 
     .      139.32266779124643*Sin(4*phi) - 
     .      40.479773535930704*Sin(5*phi) + 
     .      21.102185585291537*Sin(6*phi) + 
     .      27.763390345914974*Sin(7*phi) + 
     .      1.6622501124443485*Sin(8*phi) - 
     .      13.49400563169258*Sin(9*phi) - 
     .      10.573966878144494*Sin(10*phi) + 
     .      0.0010493967978287462*Sin(11*phi) + 
     .      5.211690149433707*Sin(12*phi) + 
     .      3.1952678133195556*Sin(13*phi) - 
     .      1.1454324071564117*Sin(14*phi) - 
     .      2.728164501523932*Sin(15*phi) - 
     .      1.2850993411366254*Sin(16*phi) + 
     .      0.7013933219589532*Sin(17*phi) + 
     .      1.1699781168734327*Sin(18*phi) + 
     .      0.301352221634545*Sin(19*phi) - 
     .      0.5291265885714342*Sin(20*phi) - 
     .      0.49114080833339335*Sin(21*phi) + 
     .      0.050614256681022596*Sin(22*phi) + 
     .      0.28881769341899144*Sin(23*phi)

       thn = 0.19622243171731607 - 
     .      0.006649098360752184*Cos(phi) - 
     .      0.01956517122663376*Cos(2*phi) - 
     .      0.008660216931195966*Cos(3*phi) - 
     .      0.006137819025589696*Cos(4*phi) - 
     .      0.0075964729910884795*Cos(5*phi) - 
     .      0.0016461860924103853*Cos(6*phi) + 
     .      0.0007294454409291066*Cos(7*phi) + 
     .      0.0023685663829071343*Cos(8*phi) + 
     .      0.0003889280477075606*Cos(9*phi) - 
     .      0.0006979450055063492*Cos(10*phi) - 
     .      0.000885537681495029*Cos(11*phi) - 
     .      0.00014878907992479928*Cos(12*phi) + 
     .      0.0003679941251949079*Cos(13*phi) + 
     .      0.0003130319474168302*Cos(14*phi) - 
     .      0.000016633051122379423*Cos(15*phi) - 
     .      0.00020931570552896087*Cos(16*phi) - 
     .      0.0001303595920774504*Cos(17*phi) + 
     .      0.0000388948168387536*Cos(18*phi) + 
     .      0.00010454833096758408*Cos(19*phi) + 
     .      0.00003792479727719437*Cos(20*phi) - 
     .      0.00004613565179395065*Cos(21*phi) - 
     .      0.00005090440921535995*Cos(22*phi) + 
     .      7.923475578463114e-6*Cos(23*phi) - 
     .      0.05678320125406135*Sin(phi) - 
     .      0.023022365710781457*Sin(2*phi) - 
     .      0.018030195256608213*Sin(3*phi) - 
     .      0.013827655583181043*Sin(4*phi) - 
     .      0.0061990077703202836*Sin(5*phi) + 
     .      0.0017061925607212282*Sin(6*phi) + 
     .      0.0023352555360747727*Sin(7*phi) - 
     .      0.0002812394123411336*Sin(8*phi) - 
     .      0.0014632759741499833*Sin(9*phi) - 
     .      0.0009892215014679607*Sin(10*phi) + 
     .      0.00012383435618247*Sin(11*phi) + 
     .      0.0005411721277811578*Sin(12*phi) + 
     .      0.00024600375235788244*Sin(13*phi) - 
     .      0.0002043069395354885*Sin(14*phi) - 
     .      0.0003033173283262694*Sin(15*phi) - 
     .      0.00010053028847115513*Sin(16*phi) + 
     .      0.0001067032704224899*Sin(17*phi) + 
     .      0.00012248942798267554*Sin(18*phi) + 
     .      5.7013723087397106e-6*Sin(19*phi) - 
     .      0.00007782884382598787*Sin(20*phi) - 
     .      0.0000523747335689761*Sin(21*phi) + 
     .      0.000021462140784394788*Sin(22*phi) + 
     .      0.00004552716184377109*Sin(23*phi)

       phn = -1.1840115989900575 - 
     .      0.2501052377367187*Cos(phi) - 
     .      0.0870108497349657*Cos(2*phi) + 
     .      0.05524474825558274*Cos(3*phi) + 
     .      0.06571634574967282*Cos(4*phi) + 
     .      0.054497177647939085*Cos(5*phi) - 
     .      0.010914652319742967*Cos(6*phi) - 
     .      0.017494764756306118*Cos(7*phi) - 
     .      0.001677474136905055*Cos(8*phi) + 
     .      0.009094152785717784*Cos(9*phi) + 
     .      0.00778609301072576*Cos(10*phi) - 
     .      0.000017569943641948104*Cos(11*phi) - 
     .      0.00427665106371633*Cos(12*phi) - 
     .      0.0031248120532568428*Cos(13*phi) + 
     .      0.0004018612433600332*Cos(14*phi) + 
     .      0.002025672163296362*Cos(15*phi) + 
     .      0.001115968407376618*Cos(16*phi) - 
     .      0.0005287555269570307*Cos(17*phi) - 
     .      0.0011232678915695546*Cos(18*phi) - 
     .      0.0005000380314326315*Cos(19*phi) + 
     .      0.00036696679782018067*Cos(20*phi) + 
     .      0.000649640645505565*Cos(21*phi) + 
     .      0.0002778176080716778*Cos(22*phi) - 
     .      0.0002750904248179207*Cos(23*phi) + 
     .      0.03756894419604881*Sin(phi) - 
     .      0.005007302851653333*Sin(2*phi) + 
     .      0.09805567331135825*Sin(3*phi) - 
     .      0.015507188585201165*Sin(4*phi) - 
     .      0.05239031186985234*Sin(5*phi) - 
     .      0.01826942178430312*Sin(6*phi) + 
     .      0.0022905823113473396*Sin(7*phi) + 
     .      0.018062815520885563*Sin(8*phi) + 
     .      0.00816576532607775*Sin(9*phi) - 
     .      0.0018070882685123845*Sin(10*phi) - 
     .      0.005789120556186378*Sin(11*phi) - 
     .      0.0018598667086740701*Sin(12*phi) + 
     .      0.0023061252186550446*Sin(13*phi) + 
     .      0.002787011388680176*Sin(14*phi) + 
     .      0.0006208190253855384*Sin(15*phi) - 
     .      0.001298615463680044*Sin(16*phi) - 
     .      0.0012840198090914765*Sin(17*phi) - 
     .      0.00011250361770319753*Sin(18*phi) + 
     .      0.00076532973526117*Sin(19*phi) + 
     .      0.0006628511134355898*Sin(20*phi) - 
     .      6.212956128694389e-6*Sin(21*phi) - 
     .      0.000505326773031062*Sin(22*phi) - 
     .      0.0004597641346084901*Sin(23*phi)

         bb0 = 1.e-5 * bb0 ! convert to gauss

!         print *,'phi0,x0,y0,z0,thn,phn,b0',phi0,x0,y0,z0,thn,phn,b0


       plat = 0.5 * pie - thn
       plon = 2. * pie + phn

       return
       end

*************************************************

!            vector 

! MS 9/15/05 
! This subroutine finds the coefficients needed to calculate the   
! parallel components of vectors.  The math is shamelessly copied from   
! the CTIP paper in the red book.  Inputs are the magnetic coordinates 
! of a point, outputs are the coefficients at that point. 
 
      subroutine vector(brad,blond,blatd,varg,vathg,vaphig) 
 
      include 'param-1.00.inc'
      include 'com-1.00.mod.inc'
 
      real brad,blond,blatd,blonr,grad,glat,glon 
      real thetaprime,cosphiprime,sinphiprime 
      real sinomega,cosomega,coplat,theta 
      real arm,athm,aphim,ax,ay,az,arp,athp,aphip,varg,vathg,vaphig 
 
!     Convert blatd and plat to colatitude, blon to radians 

      theta  = pie/2. - blatd*po180 
      coplat = pie/2. - plat 
      blonr  = blond*po180 
 
      arm = 2.* cos (theta) /  sqrt ( ( 2. * cos(theta) ) ** 2 +  
     .        sin(theta)   ** 2 ) 
      athm =     sin (theta) / sqrt ( ( 2. * cos(theta) ) ** 2 +
     .        sin(theta)   ** 2 ) 
       aphim = 0. 
 
      ax = arm * sin(theta)*cos(blonr) + athm*cos(theta)*cos(blonr) 
      ay = arm * sin(theta)*sin(blonr) + athm*cos(theta)*sin(blonr) 
      az = arm * cos(theta) - athm*sin(theta) 
 
!     Find geographic coordinates of point 

      call btog(brad,blond,blatd,grad,glat,glon) 
      glat = pie/2. - glat*po180 
      glon = glon * po180 
 
!     See p.243 of red book 

      thetaprime = acos(cos(coplat)*cos(glat) + 
     .             sin(coplat)*sin(glat)*cos(glon-plon)) 

!     Account for possibility of no tilt. 

      if (coplat .eq. 0.) then 
        cosphiprime = 1. 
        sinphiprime = 0. 
      else 
        cosphiprime = ( cos(coplat)*cos(thetaprime) - cos(glat) ) / 
     .              (sin(coplat) * sin(thetaprime)) 
        sinphiprime = sin(glat) * sin(glon-plon) / sin(thetaprime) 
      endif 
 
      arp = ax * sin(thetaprime)*cosphiprime + 
     .      ay*sin(thetaprime)*sinphiprime + 
     .      az*cos(thetaprime) 
      athp = ax * cos(thetaprime)*cosphiprime + 
     .        ay*cos(thetaprime)*sinphiprime - 
     .        az*sin(thetaprime) 
      aphip = -ax*sinphiprime + ay*cosphiprime 
 
      sinomega = -sin(coplat)*sin(glon - plon) / sin(thetaprime) 

!     MS: This formula is wrong in the red book.  Corrected here. 

      cosomega = (cos(coplat) - cos(glat)*cos(thetaprime)) / 
     .           (sin(thetaprime) * sin(glat)) 
 
      varg  = arp 
      vathg = athp*cosomega + aphip*sinomega 

!     Because of a typographical error this formula does not appear in 
!     the red book. 

      vaphig = -athp*sinomega + aphip*cosomega 
 
      return 
      end
