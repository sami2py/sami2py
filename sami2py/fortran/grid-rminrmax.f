*******************************************
*******************************************

!         GRID_RMIN_RMAX

*******************************************
*******************************************

      program main   

      parameter ( po180  = 1.745329e-02 )
      parameter ( re     = 6370.        )

      print *,'input: grad_in_min,grad_in_max,glat_in,glon_in'
      read *,grad_in_min,grad_in_max,glat_in,glon_in

!     print out input data
 
      print *,''
      print *,'Input data:'
      print *,''
      print *,'grad_in_min:',grad_in_min
      print *,'grad_in_max:',grad_in_max
      print *,'glat_in:',glat_in
      print *,'glon_in:',glon_in
      print *,''

!     rmin
      
      call gtob (  brad,blon,blat
     .            ,grad_in_min+re
     .            ,glat_in
     .            ,glon_in          )

      pvalue  = brad / re / cos(blat*po180) ** 2

      call btog ( pvalue*re,blon,0.
     .           ,grad_max
     .           ,glat_max
     .           ,glon_max             )

      print *,''
      print *,'pvalue: ',pvalue
      print *,'rmin: ',grad_max-re

!     rmax
      
      call gtob (  brad,blon,blat
     .            ,grad_in_max+re
     .            ,glat_in
     .            ,glon_in          )

      pvalue  = brad / re / cos(blat*po180) ** 2

      call btog ( pvalue*re,blon,0.
     .           ,grad_max
     .           ,glat_max
     .           ,glon_max             )

      print *,''
      print *,'pvalue: ',pvalue
      print *,'rmax: ',grad_max-re
      print *,''

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

      parameter ( plat = 1.375, plon = 5.079  )
      parameter ( po180  = 1.745329e-02 )
      parameter ( rtod   = 57.295780    )

!     coordinates of dipole in km relative to center of earth

      x0 = -392. 
      y0 =  258.
      z0 =  179.

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

!     geomagnetic latitude and longitude in degrees

      xm = xmm - x0
      ym = ymm - y0
      zm = zmm - z0

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

      parameter ( plat = 1.375, plon = 5.079  )
      parameter ( po180  = 1.745329e-02 )
      parameter ( rtod   = 57.295780    )

      real brad,blond,blatd,grad,glatd,glond

!     coordinates of dipole in km relatvie to center of earth

      x0 = -392. 
      y0 =  258.
      z0 =  179.

!     convert magnetic lat and long to radians

      blonr = blond * po180 
      blatr = blatd * po180 

!     position of point in geomagnetic coords
!     first get xm ym zm in the eccentric dipole system

      xm = brad *cos ( blatr ) * cos ( blonr )
      ym = brad *cos ( blatr ) * sin ( blonr )
      zm = brad *sin ( blatr )

!     next shift to the tilted dipole

      xmm = xm + x0
      ymm = ym + y0
      zmm = zm + z0

!     r is invariant under the rotations of the tilted dipole

      grad = sqrt ( xmm ** 2 + ymm ** 2 + zmm ** 2 )

!     rotate coords in north-south direction

      xg =  xmm * sin ( plat ) + zmm * cos ( plat )
      yg =  ymm
      zg = -xmm * cos ( plat ) + zmm * sin ( plat )

!     geographic latitude and longitude converted to degrees

      glatd = asin ( zg / grad ) * rtod
      glond = ( plon + atan2 ( yg/grad,xg/grad ) ) * rtod

      if (glond .ge. 360) glond = glond - 360.

      return
      end

