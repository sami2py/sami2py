!     *******************************************
!     *******************************************

!                  SAMI2-1.00

!     *******************************************
!     *******************************************


!                   LICENSE

!     *******************************************
!     *******************************************

!     I hereby agree to the following terms governing the use and
!     redistribution of the SAMI2 software release written and
!     developed by J.D. Huba, G. Joyce and M. Swisdak.

!     Redistribution and use in source and binary forms, with or
!     without modification, are permitted provided that (1) source code
!     distributions retain this paragraph in its entirety, (2) distributions
!     including binary code include this paragraph in its entirety in
!     the documentation or other materials provided with the distribution,
!     (3) improvements, additions and upgrades to the software will be
!     provided to NRL Authors in computer readable form, with an unlimited,
!     royalty-free license to use these improvements, additions and upgrades,
!     and the authority to grant unlimited royalty-free sublicenses to these
!     improvements and (4) all published research
!     using this software display the following acknowledgment ``This
!     work uses the SAMI2 ionosphere model written and developed
!     by the Naval Research Laboratory.''

!     Neither the name of NRL or its contributors, nor any entity of the
!     United States Government may be used to endorse or promote products
!     derived from this software, nor does the inclusion of the NRL written
!     and developed software directly or indirectly suggest NRL's or the
!     United States Government's endorsement of this product.


!     THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR IMPLIED
!     WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
!     MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.


!     *******************************************
!     *******************************************
!
      program main

      include 'param-1.00.inc'
      include 'com-1.00.inc'

!     open input files

      open ( unit=10, file='sami2py-1.00.namelist'  )
      open ( unit=11, file='exb.inp')
      open ( unit=20, file='deni-init.inp')
      open ( unit=30, file='ichem.inp')
      open ( unit=50, file='phabsdt.inp')
      open ( unit=60, file='phiondt.inp')
      open ( unit=61, file='phionnt.inp')
      open ( unit=65, file='euvflux.inp')
      open ( unit=66, file='thetant.inp')
      open ( unit=67, file='zaltnt.inp')

      call initial

!     open output files

      if ( fmtout ) then
        call open_f
      else
        call open_u
      endif

      ntm   = 0
      istep = 0
!      call output ( hrinit,ntm,istep )

      close (10)
      close (11)
      close (20)
      close (30)
      close (50)
      close (60)
      close (61)
      close (65)
      close (66)
      close (67)
      close (68)

!     time loop


      hrut    = hrinit
      timemax = hrmax * sphr
      istep   = 0
      tprnt   = 0.
      tneut   = 0.
      time    = 0.

      call neutambt (hrinit)



      do while (      istep .le. maxstep
     .          .and. time  .lt. timemax  )

!       parallel transport

        do nfl = nf,1,-1
          call zenith (hrut,nfl)
          call transprt (nfl)
        enddo

!       perpendicular transport

        call exb(hrut)

!       time/step advancement

        istep  = istep + 1
        time   = time  + dt
        hrut   = time / sphr + hrinit
        tprnt  = tprnt + dt / sphr
        tneut  = tneut + dt / sphr

        call courant ( hrut )

!       update neutrals

        if ( tneut .ge. 0.25 ) then
          call neutambt (hrut)
          tneut = 0.
        endif

!       output data

        if ( tprnt .ge. dthr .and. hrut .gt. hrpr+hrinit) then
          ntm      = ntm + 1
          call output ( hrut,ntm,istep )
          tprnt   = 0.
        elseif ( tprnt .ge. dthr ) then
          print *,'no output yet -- hour = ',hrut
          tprnt   = 0.
        endif

      enddo    ! end time loop

!     close files

      close (70)
      close (71)
      close (72)
      close (73)
      close (75)
      close (78)
      close (81)
      close (82)
      close (83)
      close (84)
      close (85)
      close (86)
      close (87)
      close (88)

      close (90)
      close (91)
      close (92)
      close (93)

      stop
      end


*******************************************
*******************************************

!            initial

*******************************************
*******************************************

      subroutine initial

      include 'param-1.00.inc'
      include 'com-1.00.inc'

      real f1026(nz,nf,91),f584(nz,nf,91),
     .     f304 (nz,nf,91),f1216(nz,nf,91)
      real zi(29),denii(29,7),d(9),temp(2),app(2),whm07(2)
      real phionr(linesuv,5),fluxdat(linesuv,2)

      namelist / go / fmtout,maxstep,hrmax,dt0,dthr,hrpr,
     .                grad_in,glat_in,glon_in,
     .                fejer,
     .                rmin,rmax,
     .                altmin,
     .                fbar,f10p7,ap,
     .                year,day,mmass,
     .                nion1,nion2,hrinit,tvn0,tvexb0,ve01,
     .                gams,gamp,snn,stn,denmin,alt_crit,cqe,
     .                Tinf_scl,euv_scl,hwm_mod,
     .                outn
!     outn added to include the neutral parameters (JMS)

!     read in parameters and initial ion density data

      read(10,go)

      dt = dt0

      ami(pthp)  = 1.
      ami(pthep) = 4.
      ami(ptnp)  = 14.
      ami(ptop)  = 16.
      ami(ptn2p) = 28.
      ami(ptnop) = 30.
      ami(pto2p) = 32.

      amn(pth)  = 1.
      amn(pthe) = 4.
      amn(ptn)  = 14.
      amn(pto)  = 16.
      amn(ptn2) = 28.
      amn(ptno) = 30.
      amn(pto2) = 32.

      alpha0(pth)  = 0.67
      alpha0(pthe) = 0.21
      alpha0(ptn)  = 1.10
      alpha0(pto)  = 0.79
      alpha0(ptn2) = 1.76
      alpha0(ptno) = 1.74
      alpha0(pto2) = 1.59

      do i = 1,7
        aap(i) = ap
      enddo

!     read in alternate ExB data

      do i = 1,10
        read(11,*) fourierA(i),fourierB(i)
      enddo

!     read in initial density data

      do i = 1,29
        read(20,102) zi(i),(denii(i,j),j=1,7)
 102    format(1x,f7.1,1p7e8.1)
      enddo


!     read in chemistry data
!     in format statement 104 need to 'hardwire' nneut (= 7)

      do k = 1,nchem
        read(30,103) (ichem(k,j),j=1,3)
 103    format(3i3)
      enddo

!     build reaction matrix

      do k = 1,nchem
        do j = 1,nneut
          do i = nion1,nion2
            ireact(i,j,k) = 0
            if (      i .eq. ichem(k,1)
     .          .and. j .eq. ichem(k,2) ) ireact(i,j,k) = 1.
          enddo
        enddo
      enddo

!     generate the mesh data

      call grid

!     output grid data

      if ( fmtout ) then
        open ( unit=69, file='zaltf.dat'   ,form='formatted' )
        open ( unit=76, file='glatf.dat'   ,form='formatted' )
        open ( unit=77, file='glonf.dat'   ,form='formatted' )
        write(69,100) alts
        write(76,100) glats
        write(77,100) glons
        close(69)
        close(76)
        close(77)
      else
        open ( unit=69, file='zaltu.dat'   ,form='unformatted' )
        open ( unit=76, file='glatu.dat'   ,form='unformatted' )
        open ( unit=77, file='glonu.dat'   ,form='unformatted' )
        write(69) alts
        write(76) glats
        write(77) glons
        close(69)
        close(76)
        close(77)
      endif

 100  format (1x,1p1e16.6)

! MS: chicrit is the zenith angle below which the Sun is visible.
! For points on the surface this is just pi/2, but at higher
! altitudes it is bigger.

      do j = 1,nf
        do i = 1,nz
          coschicrit(i,j) = cos(pie -
     .                asin( 1./ (1. + alts(i,j)/re) ))
        enddo
      enddo


!     put deni on mesh via linear interpolation
!     and put on lower limit

!     initialize all ion densities
      j0 = 1
      do n = 1,nion
        do l = 1,nf
          do i = 1,nz
            j = 1
            do while (  alts(i,l) .ge. zi(j) .and. j .le. 28 )
              j0 = j
              j  = j + 1
            enddo
            if ( n .eq. 1 ) k = pthp
            if ( n .eq. 2 ) k = pthep
            if ( n .eq. 3 ) k = ptnp
            if ( n .eq. 4 ) k = ptop
            if ( n .eq. 5 ) k = ptn2p
            if ( n .eq. 6 ) k = ptnop
            if ( n .eq. 7 ) k = pto2p
            slope   = ( denii(j0+1,n) - denii(j0,n) )
     .                / ( zi   (j0+1)   - zi   (j0) )
            deni(i,l,k) = denii(j0,n) + ( alts(i,l) - zi(j0) ) * slope
            deni(i,l,k) = amax1 ( deni(i,l,k) , denmin )
            if ( alts(i,l) .gt. zi(29) ) then
              if ( n .eq. 1 )  then
                nn = pthp
                deni(i,l,k) = denii(29,n)
              else
                deni(i,l,k) = denmin
              endif
            endif
          enddo
        enddo
      enddo


!     initialize neutrals

!     neutral density and temperature


      do j = 1,nf
        do i = 1,nz
          hrl = hrinit + glons(i,j) / 15.
          if ( hrl .ge. 24. ) hrl = hrl - 24.
          call msistim ( int(year),int(day),hrl,glons(i,j),iyd,sec )
          call gtd7 ( iyd,sec,alts(i,j),glats(i,j),glons(i,j),
     .                hrl,fbar,f10p7,aap,mmass,Tinf_scl,d,temp )

          denn(i,j,pth )  = snn(pth)  * d(7)
          denn(i,j,pthe)  = snn(pthe) * d(1)
          denn(i,j,ptn )  = snn(ptn)  * d(8)
          denn(i,j,pto )  = snn(pto)  * d(2)
          denn(i,j,ptn2)  = snn(ptn2) * d(3) + 1.e-30
          denn(i,j,pto2)  = snn(pto2) * d(4) + 1.e-30
          tn(i,j)         = stn * temp(2)
          denn(i,j,ptno)  = 0.4 * exp( -3700. / tn(i,j) )
     .                      * denn(i,j,pto2)
     .                      + 5.0e-7 * denn(i,j,pto)
        enddo
      enddo

!     electron and ion temperature initialization

      do k = nion1,nion2
        do j = 1,nf
          do i = 1,nz
            te(i,j)   = tn(i,j)
            ti(i,j,k) = tn(i,j)
          enddo
        enddo
      enddo

!     initialize ion velocity to zero

      do k = nion1,nion2
        do j = 1,nf
          do i = 1,nz
            vsi(i,j,k)     = 0.
            sumvsi(i,j,k)  = 0.
          enddo
        enddo
      enddo

!     neutral winds (convert from m/s to cm/s)

      do j = 1,nf
        do i = 1,nz
          app(1)   = ap
          app(2)   = ap
          hrl = hrinit + glons(i,j) / 15.
          if ( hrl .ge. 24. ) hrl = hrl - 24.
          call msistim ( int(year),int(day),hrl,glons(i,j),iyd,sec )
          if (hwm_mod .eq. 7) then
            call hwm07 ( iyd,sec,alts(i,j),glats(i,j),glons(i,j),
     .                hrl,fbar,f10p7,app,whm07                )
          else if (hwm_mod .eq. 14) then
            call hwm14 ( iyd,sec,alts(i,j),glats(i,j),glons(i,j),
     .                hrl,fbar,f10p7,app,whm07                )
          else
            call gws5 ( iyd,sec,alts(i,j),glats(i,j),glons(i,j),
     .                hrl,fbar,f10p7,app,whm07                )
          endif
          v(i,j)   = 100. * whm07(1) * tvn0
          u(i,j)   = 100. * whm07(2) * tvn0
        enddo
      enddo

!     read in photoabsorption rates

      do i = 1,linesuv
        read (50,105) (sigabsdt(i,j), j=1,3)
 105    format (3f7.2)
      enddo

      do j = 1,3
        do i = 1,linesuv
          sigabsdt(i,j) = tm18 * sigabsdt(i,j)
        enddo
      enddo

!     initialize photoionization rates to zero

      do j = 1,nneut
        do i = 1,linesuv
          sigidt(i,j)  = 0.
        enddo
        do i = 1,linesnt
          sigint(i,j)  = 0.
        enddo
      enddo

!     read in daytime photoionization line data
!     (only n, o, he, n_2, o_2)

      do i = 1,linesuv
        read(60,106) (phionr(i,j), j=1,5)
        sigidt(i,ptn ) = phionr(i,1)
        sigidt(i,pto ) = phionr(i,2)
        sigidt(i,pthe) = phionr(i,3)
        sigidt(i,ptn2) = phionr(i,4)
        sigidt(i,pto2) = phionr(i,5)
      enddo
 106  format(5f7.2)

      do j = 1,nion
        do i = 1,linesuv
          sigidt(i,j) = tm18 * sigidt(i,j)
        enddo
      enddo

!     read in nighttime photoionization line data
!     (only o, n_2, n0, o_2)

      do i = 1,linesnt
        read(61,106) (phionr(i,j), j=1,4)
        sigint(i,pto ) = phionr(i,1)
        sigint(i,ptn2) = phionr(i,2)
        sigint(i,ptno) = phionr(i,3)
        sigint(i,pto2) = phionr(i,4)
      enddo

      do j = 1,nion
        do i = 1,linesnt
          sigint(i,j) = tm18 * sigint(i,j)
        enddo
      enddo

!     read in f74113, ai data and set euv flux
!     (from richards et al., jgr 99, 8981, 1994)

      p  = 0.5 * ( f10p7 + fbar )

      do i = 1,linesuv
        read (65,107) (fluxdat(i,j),j=1,2)
        f74   = fluxdat(i,1)
        ai    = fluxdat(i,2)
        xflux = 1. + ai * ( p - 80.)
        if ( xflux .lt. 0.8 ) xflux = 0.8
        flux(i) = f74 * xflux * 1.e9 * euv_scl
!        if ( flux(i) .lt. 0 ) flux(i) = 0.
!        print *,'i,flux',i,flux(i)
      enddo
 107  format (f6.3,1pe11.4)

!      stop

!     read in angles for nighttime deposition fluxes

      do i = 1,linesnt
        read(66,108) (thetant(i,j), j=1,4)
      enddo
 108  format (4f7.1)

!     read in min/max altitude for nighttime deposition fluxes
!       zaltnt(i,1): zmin(i)
!       zaltnt(i,2): zmax(i)

      do i = 1,linesnt
        read(67,108) (zaltnt(i,j), j=1,2)
      enddo
 109  format (2f7.1)

!     call nighttime euv flux subroutines
!     (lyman beta 1026, he i 584, he ii 304, lyman alpha 1216)

      do j = 1,nf
        call sf1026 ( f1026,1,j )
        call sf584  ( f584 ,2,j )
        call sf304  ( f304 ,3,j )
        call sf1216 ( f1216,4,j )
        do k = 1,91
          do i = 1,nz
            fluxnt(i,j,k,1) = f1026(i,j,k)
            fluxnt(i,j,k,2) = f584 (i,j,k)
            fluxnt(i,j,k,3) = f304 (i,j,k)
            fluxnt(i,j,k,4) = f1216(i,j,k)
          enddo
        enddo
      enddo

!     intialize diagnostic variables to 0

      do j = 1,nf
        do i = 1,nz
          u1(i,j) = 0.
          u2(i,j) = 0.
          u3(i,j) = 0.
          u4(i,j) = 0.
          u5(i,j) = 0.
        enddo
      enddo

      do k = 1,nion
        do j = 1,nf
          do i = 1,nz
            t1(i,j,k) = 0.
            t2(i,j,k) = 0.
            t3(i,j,k) = 0.
          enddo
        enddo
      enddo

      print *,' finished initialization'

      return
      end

*******************************************
*******************************************

!            neutambt

*******************************************
*******************************************


!     calculate neutral densities and temperature
!     from nrlmsise00

      subroutine neutambt (hrut)


      include 'param-1.00.inc'
      include 'com-1.00.inc'

      real d(9),temp(2)
      real whm07(2),app(2)



!     no obtained from eq. (128) - bailey and balan (red book)

!     neutral density and temperature

!     input:
!        iyd - year and day as yyddd
!        sec - ut(sec)
!        alt - altitude(km) (greater than 85 km)
!        glat - geodetic latitude(deg)
!        glong - geodetic longitude(deg)
!        stl - local apparent solar time(hrs)
!        f107a - 3 month average of f10.7 flux
!        f107 - daily f10.7 flux for previous day
!        ap - magnetic index(daily) or when sw(9)=-1. :
!           - array containing:
!             (1) daily ap
!             (2) 3 hr ap index for current time
!             (3) 3 hr ap index for 3 hrs before current time
!             (4) 3 hr ap index for 6 hrs before current time
!             (5) 3 hr ap index for 9 hrs before current time
!             (6) average of eight 3 hr ap indicies from 12 to 33 hrs prior
!                    to current time
!             (7) average of eight 3 hr ap indicies from 36 to 59 hrs prior
!                    to current time
!        mass - mass number (only density for selected gas is
!                 calculated.  mass 0 is temperature.  mass 48 for all.
!     output:
!        d(1) - he number density(cm-3)
!        d(2) - o number density(cm-3)
!        d(3) - n2 number density(cm-3)
!        d(4) - o2 number density(cm-3)
!        d(5) - ar number density(cm-3)
!        d(6) - total mass density(gm/cm3)
!        d(7) - h number density(cm-3)
!        d(8) - n number density(cm-3)
!        d(9) - anomalous O (see msis)
!        t(1) - exospheric temperature
!        t(2) - temperature at alt


      do j = 1,nf
        do i = 1,nz
          hrl = hrut + glons(i,j) / 15.
          if ( hrl .ge. 24. ) hrl = hrl - 24.
          call msistim ( int(year),int(day),hrl,glons(i,j),iyd,sec )
          call gtd7 ( iyd,sec,alts(i,j),glats(i,j),glons(i,j),
     .                hrl,fbar,f10p7,aap,mmass,Tinf_scl,d,temp )
          denn(i,j,pth )  = snn(pth)  * d(7)
          denn(i,j,pthe)  = snn(pthe) * d(1)
          denn(i,j,ptn )  = snn(ptn)  * d(8)
          denn(i,j,pto )  = snn(pto)  * d(2)
          denn(i,j,ptn2)  = snn(ptn2) * d(3) + 1.e-30
          denn(i,j,pto2)  = snn(pto2) * d(4) + 1.e-30
          tn(i,j)         = stn * temp(2)
          denn(i,j,ptno)  = 0.4 * exp( -3700. / tn(i,j) )
     .                      * denn(i,j,pto2)
     .                      + 5.0e-7 * denn(i,j,pto)
        enddo
      enddo

!     neutral winds

!        iyd - year and day as yyddd
!        sec - ut(sec)  (not important in lower atmosphere)
!        alt - altitude(km)
!        glat - geodetic latitude(deg)
!        glong - geodetic longitude(deg)
!        stl - local apparent solar time(hrs)
!        f107a - 3 month average of f10.7 flux (use 150 in lower atmos.)
!        f107 - daily f10.7 flux for previous day ( " )
!        ap - two element array with
!             ap(1) = magnetic index(daily) (use 4 in lower atmos.)
!             ap(2)=current 3hr ap index (used only when sw(9)=-1.)
!     note:  ut, local time, and longitude are used independently in the
!            model and are not of equal importance for every situation.
!            for the most physically realistic calculation these three
!            variables should be consistent.
!      output
!        w(1) = meridional (m/sec + northward)
!        w(2) = zonal (m/sec + eastward)

      do j = 1,nf
        do i = 1,nz
          app(1)   = ap
          app(2)   = ap
          hrl = hrut + glons(i,j) / 15.
          if ( hrl .ge. 24. ) hrl = hrl - 24.
          call msistim ( int(year),int(day),hrl,glons(i,j),iyd,sec )
          if (hwm_mod .eq.7) then
            call hwm07 ( iyd,sec,alts(i,j),glats(i,j),glons(i,j),
     .                hrl,fbar,f10p7,app,whm07                )
          else if (hwm_mod .eq.14) then
            call hwm14 ( iyd,sec,alts(i,j),glats(i,j),glons(i,j),
     .                hrl,fbar,f10p7,app,whm07                )
          else
            call gws5 ( iyd,sec,alts(i,j),glats(i,j),glons(i,j),
     .                hrl,fbar,f10p7,app,whm07                )
          endif
          v(i,j)   = 100. * whm07(1) * tvn0  ! convert to cm/sec
          u(i,j)   = 100. * whm07(2) * tvn0   ! convert to cm/sec
        enddo
      enddo

      return
      end

*******************************************
*******************************************

!            transprt

*******************************************
*******************************************

      subroutine transprt (nfl)

      include 'param-1.00.inc'
      include 'com-1.00.inc'

      real prod(nz,nion),loss(nz,nion),lossr,
     .     phprodr(nz,nion),chrate(nz,nchem),
     .     chloss(nz,nion),chprod(nz,nion),relossr(nz,nion)
      real deni_old(nz,nion),te_old(nz),ti_old(nz,nion),vsi_old(nz,nion)
      real tvn(nz)
      real nuin(nz,nion,nneut),
     .     nuij(nz,nion,nion),sumnuj(nz,nion)
      real vsin(nz,nion),vsidn(nz,nion),denin(nz,nion)
      real ten(nz),tin(nz,nion)

!     calculation of production and loss
!       phprodr: photo production rates
!       chrate:  chemical rates (ichem)
!       chloss:  chemical loss term
!       chprod:  chemical production term
!       relossr: recombination loss rates

!     initialize tvn and gs

      do i = 1,nz
        tvn(i) = 0.
        gs(i)  = 0.
      enddo

      do i = 1,nz

        ne(i,nfl)   = 1.
        te_old(i)   = te(i,nfl)
        do j = nion1,nion2
          deni_old(i,j) = deni(i,nfl,j)
          ne(i,nfl)     = ne(i,nfl) + deni(i,nfl,j)
          ti_old(i,j)   = ti(i,nfl,j)
          vsi_old(i,j)  = vsi(i,nfl,j)
        enddo

       enddo

       call photprod ( cx,phprodr,nfl   )         ! calculates phprodr
       call chemrate ( chrate,nfl               ) ! calculates chrate
       call chempl   ( chrate,chloss,chprod,nfl ) ! calcualtes chloss,chprod
       call recorate ( relossr,nfl              ) ! calculates relossr

       do i = 1,nz

        do j = nion1,nion2
          prod  (i,j) =  phprodr(i,j) * denn(i,nfl,j)
     .                   + chprod(i,j)
          lossr       =  relossr(i,j) * deni(i,nfl,j) * ne(i,nfl)
     .                   + chloss(i,j)
          loss (i,j)  =  lossr / deni(i,nfl,j)
        enddo

!       gravity and neutral wind
!       modified 9/19/05 (MS)

        gs(i)   =  gzero * arg(i,nfl)
     .             * ( re / (re + alts(i,nfl)) ) ** 2

        tvn(i)  = (  v(i,nfl) * athg(i,nfl)
     .             - u(i,nfl) * aphig(i,nfl) )

        tvn(i)    = tvn0 * tvn(i) ! tvn0 used to modify tvn

        u3(i,nfl) = u(i,nfl)
        u4(i,nfl) = v(i,nfl)
        u5(i,nfl) = tvn(i)

      enddo

      call update ( tvn,nuin,sumnuj,nuij,nfl )

      do i = 1,nz
        do nni = nion1,nion2
          sumvsi(i,nfl,nni) = 0.
          do nj = nion1,nion2
          sumvsi(i,nfl,nni) =   sumvsi(i,nfl,nni)
     .                     + nuij(i,nni,nj)*vsi(i,nfl,nj)
          enddo
        enddo
      enddo

!     define new arrays for velocity and density

      do ni = nion1,nion2
        do i = 1,nz
          vsin (i,ni) = vsi(i,nfl,ni)
          vsidn(i,ni) = vsid(i,nfl,ni)
          denin(i,ni) = deni(i,nfl,ni)
        enddo
      enddo

!     update variables

      do ni = nion1,nion2

        call vsisolv ( vsin(1,ni),vsidn(1,ni),vsi_old(1,ni)
     .                ,sumnuj(1,ni),nfl )

! compensating filter

        call smoothz ( vsin(1,ni), 1 )

!       put stuff back into velocity array

        do i = 1,nz
          vsi(i,nfl,ni)  = vsin(i,ni)
          vsid(i,nfl,ni) = vsidn(i,ni)
        enddo

        if ( nfl .eq. 1) then
          do i = 1,nz
            vsi(i,1,ni) = vsi(i,2,ni)
          enddo
        endif

        if ( nfl .eq. nf ) then
          do i = 1,nz
            vsi(i,nf,ni) = vsi(i,nf-1,ni)
          enddo
        endif

        call densolv2 ( ni,denin(1,ni)
     .       ,prod(1,ni),loss(1,ni),deni_old(1,ni),nfl )

!       put stuff back into density array

        do i = 1,nz
          deni(i,nfl,ni) = denin(i,ni)
        enddo

!       put floor on density

        do i = 1,nz
          deni(i,nfl,ni) = amax1 ( deni(i,nfl,ni), denmin )
        enddo

      enddo

!     define new arrays for temperature

      do ni = nion1,nion2
        do i = 1,nz
          tin(i,ni)  = ti(i,nfl,ni)
        enddo
      enddo

      do i = 1,nz
        ten(i)  = te(i,nfl)
      enddo

!     temperatures (with floors and warnings)

      call etemp  (ten,te_old,phprodr,nfl)
      do i = 1,nz
        te(i,nfl)  = amax1(tn(i,nfl),ten(i))
        te(i,nfl)  = amin1(te(i,nfl),1.e4)
        if ( te(i,nfl) .lt. 0 ) then
          print *,' T(e) negative: i,nfl',i,nfl
          stop
        endif
      enddo

      call htemp  (tin(1,pthp) ,ti_old(1,pthp) ,tvn,nuin,nfl)
      do i = 1,nz
        ti(i,nfl,pthp)  = amax1(tn(i,nfl),tin(i,pthp))
        ti(i,nfl,pthp)  = amin1(ti(i,nfl,pthp),1.e4)
        if ( ti(i,nfl,pthp) .lt. 0 ) then
          print *,' T(H) negative: i,nfl',i,nfl
          stop
        endif
      enddo

      call hetemp (tin(1,pthep),ti_old(1,pthep),tvn,nuin,nfl)
      do i = 1,nz
        ti(i,nfl,pthep)  = amax1(tn(i,nfl),tin(i,pthep))
        ti(i,nfl,pthep)  = amin1(ti(i,nfl,pthep),1.e4)
        if ( ti(i,nfl,pthep) .lt. 0 ) then
          print *,' T(He) negative: i,nfl',i,nfl
          stop
        endif
      enddo

      call otemp  (tin(1,ptop) ,ti_old(1,ptop) ,tvn,nuin,nfl)
      do i = 1,nz
        ti(i,nfl,ptop)  = amax1(tn(i,nfl),tin(i,ptop))
        ti(i,nfl,ptop)  = amin1(ti(i,nfl,ptop),1.e4)
        if ( ti(i,nfl,ptop) .lt. 0 ) then
          print *,' T(O) negative: i,nfl',i,nfl
          stop
        endif
      enddo

      do i = 1,nz
        ti(i,nfl,ptnp )    = ti(i,nfl,ptop)
        ti(i,nfl,ptn2p)    = ti(i,nfl,ptop)
        ti(i,nfl,ptnop)    = ti(i,nfl,ptop)
        ti(i,nfl,pto2p)    = ti(i,nfl,ptop)
      enddo

      return
      end


*******************************************
*******************************************

!            photprod

*******************************************
*******************************************

!     photoproduction rates

      subroutine photprod ( cxl,phprodr,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'

      real cxl(nz,nf)
      real phprodr(nz,nion),xmass(3)
      integer idx(3)

!     scale height of neutral atmosphere

      hcof = 1.e-5 * bolt / ( gzero * amu * re ** 2 )

      do iz = 1,nz

         coschi = cxl(iz,nfl)

         do j = nion1,nion2
            phprodr ( iz,j ) = 0.
         enddo

!     only consider o, n2, o2 for absorption

         idx(1) = pto
         idx(2) = ptn2
         idx(3) = pto2

         rp    = alts(iz,nfl) + re
         rp2   = rp * rp

!         if ( coschi .ge. 0. ) then ! sun is up
         if ( coschi .ge. coschicrit(iz,nfl) ) then ! sun is up

!     daytime deposition

            do i = 1,3
               hscale   = hcof * tn(iz,nfl) * rp2 / amn(idx(i))
               xscale   = rp / hscale
               y1       = sqrt ( .5 * xscale ) * abs(coschi)
               ch1      = atm_chapman(xscale,rtod*acos(coschi))
               if (ch1 .gt. 1.e22) ch1 = 1.e22
               xmass(i) = denn(iz,nfl,idx(i)) * hscale * ch1 * 1.e5
            enddo

            do l=1,linesuv
               exa =   xmass(1) * sigabsdt(l,1)
     .              + xmass(2) * sigabsdt(l,2)
     .              + xmass(3) * sigabsdt(l,3)
               if(exa .gt. 35.) exa = 35.
               flx = flux(l) * exp(-exa)
               do j=nion1,nion2
                  phprodr(iz,j) = phprodr(iz,j) + sigidt(l,j) * flx
               enddo
            enddo

            ang    = acos ( coschi )
            itheta = nint ( ang / po180 ) - 90
            itheta = int ( amax1 ( float(itheta), 1. ) )
            do l = 1,linesnt
               do j=nion1,nion2
                  phprodr(iz,j) =   phprodr(iz,j)
     .                 + sigint(l,j) * fluxnt(iz,nfl,itheta,l)
               enddo
            enddo
!     nighttime deposition

         else                   ! sun is dowm

            ang    = acos ( coschi )
            itheta = nint ( ang / po180 ) - 90
            itheta = int ( amax1 ( float(itheta), 1. ) )
            do l = 1,linesnt
               do j=nion1,nion2
                  phprodr(iz,j) =   phprodr(iz,j)
     .                 + sigint(l,j) * fluxnt(iz,nfl,itheta,l)
               enddo
            enddo

         endif
      enddo

      return
      end

*******************************************
*******************************************

!            chemrate

*******************************************
*******************************************

!     chemical producation and loss rates
!     bb: bailley and balan (red book, 1996)

      subroutine chemrate ( chrate,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'

      real chrate(nz,nchem)

      do iz = 1,nz

      ti300o = ti(iz,nfl,ptop) / 300.

      chrate (iz,1) = 2.2e-11
     .            * sqrt( ti(iz,nfl,pthp) )       ! h+ + o --> o+ + h (bb)

      chrate (iz,2) = 3.5e-10                        ! he+ + n2 --> n2+ + he (bb)

      chrate (iz,3) = 8.5e-10                        ! he+ + n2 --> n+ + n + he (schunk)

      chrate (iz,4) = 8.0e-10                        ! he+ + o2 --> o+ + o + he (bb)

      chrate (iz,5) = 2.0e-10                        ! he+ + o2 --> o2+ + he

      chrate (iz,6) = 2.0e-10                        ! n+ + o2 --> no+ + o  (schunk)

      chrate (iz,7) = 4.0e-10                        ! n+ + o2 --> o2+ + n(2d) (schunk)

      chrate (iz,8) = 1.0e-12                        ! n+ + 0 --> o+ + n

      chrate (iz,9) = 2.0e-11                        ! n+ + no --> no+ + o (schunk)

      chrate(iz,10) = 2.5e-11
     .             * sqrt( tn(iz,nfl) )           ! o+ + h --> h+ + o   (bb)

      chrate(iz,11) = 1.533e-12 -                    ! o+ + n2 --> no+ + n (bb)
     .             5.920e-13 * ti300o +
     .             8.600e-14 * ti300o ** 2
      if ( ti(iz,nfl,ptop) .gt. 1700 )
     .  chrate(iz,11) = 2.730e-12 -
     .                1.155e-12 * ti300o +
     .                1.483e-13 * ti300o ** 2

      chrate(iz,12) = 2.820e-11 -                    ! o+ + o2 --> o2+ + o
     .             7.740e-12 * ti300o +
     .             1.073e-12 * ti300o ** 2 -
     .             5.170e-14 * ti300o ** 3 +
     .             9.650e-16 * ti300o ** 4

      chrate(iz,13) = 1.0e-12                        ! o+ + no --> no+ + o

      chrate(iz,14) = 1.4e-10 / ti300o ** .44        ! n2+ + o --> no+ + n(2d) (bb)

      chrate(iz,15) = 5.0e-11 / sqrt( ti300o )       ! n2+ + o2 --> o2+ + n2 (schunk)

      chrate(iz,16) = 1.0e-14                        ! n2+ + o2 --> no+ + no

      chrate(iz,17) = 3.3e-10                        ! n2+ + no --> no+ + n2 (schunk)

      chrate(iz,18) = 1.2e-10                        ! o2+ + n --> no+ + o (schunk)

      chrate(iz,19) = 2.5e-10                        ! o2+ + n(2d) --> n+ + o2

      chrate(iz,20) = 4.4e-10                        ! o2+ + no --> no+ + o2 (bb)

      chrate(iz,21) = 5.0e-16                        ! o2+ + n2 --> no+ + no (schunk)


      enddo

      return
      end

*******************************************
*******************************************

!            recorate

*******************************************
*******************************************

!     recombination rates
!     bb: bailley and balan (red book, 1996)

      subroutine recorate ( relossr,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'

      real relossr(nz,nion)

      do iz = 1,nz

        te300 = te(iz,nfl) / 300.

        relossr(iz,pthp)  = 4.43e-12 / te300 ** .7
        relossr(iz,pthep) = relossr(iz,pthp)
        relossr(iz,ptnp)  = relossr(iz,pthp)
        relossr(iz,ptop)  = relossr(iz,pthp)
        relossr(iz,ptn2p) = 1.8e-7 / te300 ** 0.39     !   (schunk)
        relossr(iz,ptnop) = 4.2e-7 / te300 ** 0.85     !   (bb)
        relossr(iz,pto2p) = 1.6e-7 / te300 ** 0.55     !   (schunk)

      enddo

      return
      end

*******************************************
*******************************************

!          chempl

*******************************************
*******************************************

!     chemical loss (chloss) and production (chprod)

!     chrate: chemical reaction rates calculated in chemrate
!     ichem: input data file showing loss, neutral, production
!            species for each reaction

      subroutine chempl ( chrate,chloss,chprod,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'

      real chrate(nz,nchem),chloss(nz,nion),chprod(nz,nion)

      do i = nion1,nion2
        do iz = 1,nz
          chloss(iz,i)   = 0.
          chprod(iz,i)   = 0.
        enddo
      enddo

      do k = 1,nchem
        il = ichem(k,1) ! ion species (reacting) loss
        in = ichem(k,2) ! neutral species reacting
        ip = ichem(k,3) ! ion species produced
        do iz = 1,nz
           chem  = denn(iz,nfl,in) * chrate(iz,k)
           tdeni = deni(iz,nfl,il) * chem
           chloss(iz,il) = tdeni + chloss(iz,il)
           chprod(iz,ip) = tdeni + chprod(iz,ip)
        enddo
      enddo

      return
      end

*******************************************
*******************************************

!            update

*******************************************
*******************************************

      subroutine update ( tvn,nuin,sumnuj,nuij,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'

      real nuin(nz,nion,nneut),nuij(nz,nion,nion)
      real nuint(nz,nion)
      real sumnuj(nz,nion),nufacij,nufacin
      real tvn(nz)
      real k0,mi

!     ion-neutral collision frequency

!     initialize everything to 0

      do nn = 1,nneut
        do ni = nion1,nion2
          do iz = 1,nz
            nuin (iz,ni,nn) = 0.
            nuint(iz,ni)    = 0.
          enddo
        enddo
      enddo

!     collision frequencies/factors

!     hydrogen (H)

      ni = pthp
      do nn = 1,nneut
         do i = 1,nz
          if ( nn .eq. pto ) then

!     MS: According to both the SAMI2 paper and the red book the
!     temperature used here should be the H+ temperature, and not the
!     hybrid used in the other terms. I've changed that.

!            teff    = 0.5 * ( ti(i,nfl,ni) + tn(i,nfl) ) ! original
            teff    = ti(i,nfl,ni)
            fac     = ( 1.00 - .047 * alog10(teff) ) ** 2
            tfactor = sqrt(teff) * fac
            nuin(i,ni,nn)  = 6.61e-11 * denn(i,nfl,nn) * tfactor
          else
            amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
            amimn   = amn(nn) / ( ami(ni) + amn(nn) )
            nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
            nuin(i,ni,nn) = nufacin * denn(i,nfl,nn)
          endif
          nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
      enddo

!     helium (He)

      ni = pthep
      do nn = 1,nneut
         do i = 1,nz
          amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
          amimn   = amn(nn) / ( ami(ni) + amn(nn) )
          nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
          nuin(i,ni,nn) = nufacin * denn(i,nfl,nn)
          nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
      enddo

!     nitrogen (N)

      ni = ptnp
      do nn = 1,nneut
         do i = 1,nz
          amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
          amimn   = amn(nn) / ( ami(ni) + amn(nn) )
          nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
          nuin(i,ni,nn) = nufacin * denn(i,nfl,nn)
          nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
      enddo

!     oxygen (O)

      ni = ptop
      do nn = 1,nneut
         do i = 1,nz
          if ( nn .eq. pto ) then
            teff    = 0.5 * ( ti(i,nfl,ni) + tn(i,nfl) )
            fac     = ( 1.04 - .067 * alog10(teff) ) ** 2
            tfactor = sqrt(teff) * fac
            nuin(i,ni,nn)  = 4.45e-11 * denn(i,nfl,nn) * tfactor
          else
            amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
            amimn   = amn(nn) / ( ami(ni) + amn(nn) )
            nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
            nuin(i,ni,nn) = nufacin * denn(i,nfl,nn)
          endif
          nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
      enddo

!     nitrogen 2(N2)

      ni = ptn2p
      do nn = 1,nneut
         do i = 1,nz
          if ( nn .eq. ptn2 ) then
            teff    = 0.5 * ( ti(i,nfl,ni) + tn(i,nfl) )
!       MS: According to the SAMI2 paper the first coefficient in
!       fac should be 1.04, not 1.00.
!            fac     = ( 1.00 - .069 * alog10(teff) ) ** 2 ! original
            fac     = ( 1.04 - .069 * alog10(teff) ) ** 2
            tfactor = sqrt(teff) * fac
            nuin(i,ni,nn) = 5.14e-11 * denn(i,nfl,nn) * tfactor
          else
            amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
            amimn   = amn(nn) / ( ami(ni) + amn(nn) )
            nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
            nuin(i,ni,nn) = nufacin * denn(i,nfl,nn)
          endif
          nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
      enddo

!     nitrous oxide (N0)

      ni = ptnop
      do nn = 1,nneut
         do i = 1,nz
          amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
          amimn   = amn(nn) / ( ami(ni) + amn(nn) )
          nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
          nuin(i,ni,nn) = nufacin * denn(i,nfl,nn)
          nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
      enddo

!     oxygen 2(O2)

      ni = pto2p
      do nn = 1,nneut
         do i = 1,nz
          if ( nn .eq. pto2 ) then
            teff    = 0.5 * ( ti(i,nfl,ni) + tn(i,nfl) )
            fac     = ( 1.00 - .073 * alog10(teff) ) ** 2
            tfactor = sqrt(teff) * fac
            nuin(i,ni,nn) = 2.59e-11 * denn(i,nfl,nn) * tfactor
          else
            amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
            amimn   = amn(nn) / ( ami(ni) + amn(nn) )
            nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
            nuin(i,ni,nn) = nufacin * denn(i,nfl,nn)
          endif
          nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
      enddo

!     ion-ion collision frequency

      do nj = nion1,nion2
        do ni = nion1,nion2
          do i = 1,nz
            nuij(i,ni,nj) = 0.
          enddo
        enddo
      enddo

      do nj = nion1,nion2
        do ni = nion1,nion2
          if ( ni .ne. nj ) then
             do i = 1,nz
              alame1  = ( ami(ni) + ami(nj) ) * evtok /
     .                ( ami(ni)*ti(i,nfl,nj) + ami(nj)*ti(i,nfl,ni) )
              alame2  = deni(i,nfl,ni) * evtok / ti(i,nfl,ni) +
     .                  deni(i,nfl,nj) * evtok / ti(i,nfl,nj)
              if ( alame2 .lt. 0 ) then
                print *,'ni,i,nj,nfl,tii,tij,alame1,alame2,nii,nij',
     .                   ni,i,nj,nfl,ti(i,nfl,ni),ti(i,nfl,nj),
     .                   alame1,alame2,
     .                   deni(i,nfl,ni),deni(i,nfl,nj)
                print *,' code has gone bad ...'
                stop
              endif
              alame   = alame1 * sqrt(alame2)
              alam    = 23. - alog(alame)
              amufac  = (ami(nj)/ami(ni))/(ami(ni) +ami(nj))
              nufacij = 9.2e-2*alam*sqrt(amufac)
              nuij(i,ni,nj) =  nufacij * deni(i,nfl,nj)
     .                        / sqrt( ti(i,nfl,ni)**3 )
          enddo
          endif
        enddo
      enddo


 100  format(1x,2e12.2)

!     sumnuj: sum of ion-ion coll freq and nuin

      do ni = nion1,nion2
        do i = 1,nz
          sumnuj(i,ni) = 0.
          do nj = nion1,nion2
            sumnuj(i,ni) = sumnuj(i,ni) + nuij(i,ni,nj)
          enddo
          sumnuj(i,ni) = sumnuj(i,ni) + nuint(i,ni)
        enddo
      enddo

!     update ne

!     modified do loop to speed up

!      do i = 1,nz
!      ne(i,nfl) = 1.
!        do ni = nion1,nion2
!          ne(i,nfl) = ne(i,nfl) + deni(i,nfl,ni)
!        enddo
!      enddo

      do ni = nion1,nion2
        do i = 1,nz
          ne(i,nfl) = 1.
        enddo
      enddo

      do ni = nion1,nion2
        do i = 1,nz
          ne(i,nfl) = ne(i,nfl) + deni(i,nfl,ni)
        enddo
      enddo

!     get a new value for vsid

!     reversed order of do loops

!        do i = 2,nz-1              ! original
!          do ni = nion1,nion2      ! original

      do ni = nion1,nion2
        do i = 2,nz-1
          mi    = amu * ami(ni)
          k0    = bolt / mi
          term1 = nuint(i,ni) * tvn(i) + sumvsi(i,nfl,ni) + gs(i)
          pip   = 0.5 * (   deni(i+1,nfl,ni) * ti(i+1,nfl,ni)
     .                    + deni(i,nfl,ni)   * ti(i,nfl,ni)   )
          pim   = 0.5 * (   deni(i,nfl,ni)   * ti(i,nfl,ni)
     .                    + deni(i-1,nfl,ni) * ti(i-1,nfl,ni) )
          denid =
     .             (        deni(i-1,nfl,ni)
     .               + 4. * deni(i,nfl,ni)
     .               +      deni(i+1,nfl,ni)  ) / 6.
          term2 =  - bms(i,nfl) * k0 /  denid
     .             * ( pip - pim ) / d22s(i,nfl)
          pep   = 0.5 * (   ne(i+1,nfl) * te(i+1,nfl)
     .                    + ne(i,nfl)   * te(i,nfl)   )
          pem   = 0.5 * (   ne(i,nfl)   * te(i,nfl)
     .                    + ne(i-1,nfl) * te(i-1,nfl) )
          dened =
     .  ( ne(i-1,nfl) + 4. * ne(i,nfl) + ne(i+1,nfl) ) / 6.
          term3 =  - bms(i,nfl) * k0 /  dened
     .             * ( pep - pem ) / d22s(i,nfl)

          vsid(i,nfl,ni)  =  term1 + term2 + term3

          if ( deni(i,nfl,ni) .le. .0001*ne(i,nfl) )
     .      vsid(i,nfl,ni) =   vsid(i,nfl,ni)
     .                       * exp ( -.0001*ne(i,nfl)/deni(i,nfl,ni) )

        enddo
      enddo

!     fix up end points for vsid

      do ni = nion1,nion2
        vsid (1,nfl,ni)    = vsid (2,nfl,ni)
        vsid (nz,nfl,ni)   = vsid (nz-1,nfl,ni)
      enddo

!     calculate collisional ion velocity
!     not used; simply a diagnostic

!      do i = 1,nz
!        do ni = nion1,nion2
!          vsic(i,nfl,ni) = vsid(i,nfl,ni) / sumnuj(i,ni)
!        enddo
!      enddo

      return
      end



*******************************************
*******************************************

!            htemp

*******************************************
*******************************************

      subroutine htemp ( tti,tiold,tvn,nuin,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'

      real tiold(nz),kapi(nz),s1i(nz),s2i(nz),s3i(nz),s4i(nz),s5i(nz)
      real tvn(nz),nuin(nz,nion,nneut),s6i(nz),s7i(nz),tti(nz)
      real divvexb(nz)

      convfac = amu / bolt / 3.

      do i = 1,nz
        s1i(i)  = 0.
        s2i(i)  = 0.
        s3i(i)  = 0.
        s4i(i)  = 0.
        s5i(i)  = 0.
        s6i(i)  = 0.
        s7i(i)  = 0.
        kapi(i) = 0.
      enddo

      do i = 1,nz

        kapi(i) = 4.6e+4 * sqrt ( ti(i,nfl,pthp)**5 ) / ne(i,nfl) *
     .            deni(i,nfl,pthp) / sqrt(ami(pthp))

        kapi(i)  = 0.6667 * kapi(i)  * evtok

!       neutrals

        do nn = 1,nneut
          redmass =
     .     ami(pthp) * amn(nn) / ( ami(pthp) + amn(nn) ) ** 2
          s2i(i) = s2i(i) + 2. * nuin(i,pthp,nn) * redmass
          s3i(i) = s3i(i)
     .      + convfac * amn(nn)
     .                * abs ( vsi(i,nfl,pthp) - tvn(i) ) ** 2
     .      * 2. * nuin(i,pthp,nn) * redmass
        enddo

        s1i(i) = s2i(i) * tn(i,nfl)

!       electrons

        s4i(i) = 7.7e-6 * ne(i,nfl) / ami(pthp)
     .                   / te(i,nfl) / sqrt(te(i,nfl))
     .                   * .66667 * evtok
        s5i(i) = s4i(i) * te(i,nfl)

!       other ions

        do ni = nion1,nion2
!          if ( ni .ne. ptop ) then
          if ( ni .ne. pthp ) then
            tfac    =    ti(i,nfl,pthp) / ami(pthp)
     .                +  ti(i,nfl,ni) / ami(ni)
            xs6i    = 3.3e-4 * deni(i,nfl,ni) / ami(pthp) / ami(ni)
     .                / tfac / sqrt(tfac) * .66667 * evtok
            xs7i    = xs6i * ti(i,nfl,ni)
            s6i(i) = s6i(i) + xs6i
            s7i(i) = s7i(i) + xs7i
          endif
        enddo

      enddo

! MS: Neglected term, divergence of ExB drift
! Divergence of the ExB drift; requires equatorial drift

      nzh = (nz+1)/2
      vexbeq = vexb(nzh,nfl)
      do i = 1,nz
        divvexb(i) = 6.*vexbeq /
     .               (ps(i,nfl)*re*1.e5) *
     .               cos(blats(i,nfl)*po180)**2 *
     .               (1.+sin(blats(i,nfl)*po180)**2) /
     .               (1.+3.*sin(blats(i,nfl)*po180)**2)**2
        s2i(i) = s2i(i) - 0.6667 * divvexb(i)
      enddo

      call tisolv(tti,tiold,kapi,s1i,s2i,s3i,s4i,s5i,s6i,s7i,pthp,nfl)

      return
      end

*******************************************
*******************************************

!            hetemp

*******************************************
*******************************************

      subroutine hetemp ( tti,tiold,tvn,nuin,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'

      real tiold(nz),kapi(nz),s1i(nz),s2i(nz),s3i(nz),s4i(nz),s5i(nz)
      real tvn(nz),nuin(nz,nion,nneut),s6i(nz),s7i(nz),tti(nz)
      real divvexb(nz)

      convfac = amu / bolt / 3.

      do i = 1,nz
        s1i(i)  = 0.
        s2i(i)  = 0.
        s3i(i)  = 0.
        s4i(i)  = 0.
        s5i(i)  = 0.
        s6i(i)  = 0.
        s7i(i)  = 0.
        kapi(i) = 0.
      enddo

      do i = 1,nz

        kapi(i) = 4.6e+4 * sqrt ( ti(i,nfl,pthep)**5 ) / ne(i,nfl) *
     .            deni(i,nfl,pthep) / sqrt(ami(pthep))
        kapi(i)  = 0.6667 * kapi(i) * evtok

!       neutrals

        do nn = 1,nneut
          redmass =
     .     ami(pthep) * amn(nn) / ( ami(pthep) + amn(nn) ) ** 2
          s2i(i) = s2i(i) + 2. * nuin(i,pthep,nn) * redmass
          s3i(i) = s3i(i)
     .      + convfac * amn(nn)
     .                * abs ( vsi(i,nfl,pthep) - tvn(i) ) ** 2
     .      * 2. * nuin(i,pthep,nn) * redmass
        enddo

        s1i(i) = s2i(i) * tn(i,nfl)

!       electrons

        s4i(i) = 7.7e-6 * ne(i,nfl) / ami(pthep)
     .                   / te(i,nfl) / sqrt(te(i,nfl))
     .                   * .66667 * evtok
        s5i(i) = s4i(i) * te(i,nfl)

!       other ions

        do ni = nion1,nion2
!          if ( ni .ne. ptop ) then
          if ( ni .ne. pthep ) then
            tfac    =   ti(i,nfl,pthep) / ami(pthep)
     .                + ti(i,nfl,ni) / ami(ni)
            xs6i    = 3.3e-4 * deni(i,nfl,ni) / ami(pthep) / ami(ni)
     .                / tfac / sqrt(tfac) * .66667 * evtok
            xs7i    = xs6i * ti(i,nfl,ni)
            s6i(i) = s6i(i) + xs6i
            s7i(i) = s7i(i) + xs7i
          endif
        enddo

      enddo

! MS: Neglected term, divergence of ExB drift
! Divergence of the ExB drift; requires equatorial drift

      nzh = (nz+1)/2
      vexbeq = vexb(nzh,nfl)
      do i = 1,nz
        divvexb(i) = 6.*vexbeq /
     .               (ps(i,nfl)*re*1.e5) *
     .               cos(blats(i,nfl)*po180)**2 *
     .               (1.+sin(blats(i,nfl)*po180)**2) /
     .               (1.+3.*sin(blats(i,nfl)*po180)**2)**2
        s2i(i) = s2i(i) - 0.6667 * divvexb(i)
      enddo

      call tisolv(tti,tiold,kapi,s1i,s2i,s3i,s4i,s5i,s6i,s7i,pthep,nfl)

      return
      end

*******************************************
*******************************************

!            otemp

*******************************************
*******************************************

      subroutine otemp ( tti,tiold,tvn,nuin,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'

      real tiold(nz),kapi(nz),s1i(nz),s2i(nz),s3i(nz),s4i(nz),s5i(nz)
      real tvn(nz),nuin(nz,nion,nneut),s6i(nz),s7i(nz),tti(nz)
      real divvexb(nz)

      convfac = amu / bolt / 3.

      do i = 1,nz
        s1i(i)  = 0.
        s2i(i)  = 0.
        s3i(i)  = 0.
        s4i(i)  = 0.
        s5i(i)  = 0.
        s6i(i)  = 0.
        s7i(i)  = 0.
        kapi(i) = 0.
      enddo

      do i = 1,nz

        kapi(i) = 4.6e+4 * sqrt ( ti(i,nfl,ptop)**5 ) / ne(i,nfl) *
     .            deni(i,nfl,ptop) / sqrt(ami(ptop))
        kapi(i)  = 0.6667 * kapi(i) * evtok

!       neutrals

        do nn = 1,nneut
          redmass =
     .     ami(ptop) * amn(nn) / ( ami(ptop) + amn(nn) ) ** 2
          s2i(i) = s2i(i) + 2. * nuin(i,ptop,nn) * redmass
          s3i(i) = s3i(i)
     .      + convfac * amn(nn)
     .                * abs ( vsi(i,nfl,ptop) - tvn(i) ) ** 2
     .      * 2. * nuin(i,ptop,nn) * redmass
        enddo

        s1i(i) = s2i(i) * tn(i,nfl)

!       electrons

        s4i(i) = 7.7e-6 * ne(i,nfl) / ami(ptop)
     .                   / te(i,nfl) / sqrt(te(i,nfl))
     .                   * .66667 * evtok
        s5i(i) = s4i(i) * te(i,nfl)

!       other ions

        do ni = nion1,nion2
          if ( ni .ne. ptop ) then
            tfac    =    ti(i,nfl,ptop) / ami(ptop)
     .                 + ti(i,nfl,ni) / ami(ni)
            xs6i    = 3.3e-4 * deni(i,nfl,ni) / ami(ptop) / ami(ni)
     .                / tfac / sqrt(tfac) * .66667 * evtok
            xs7i    = xs6i * ti(i,nfl,ni)
            s6i(i) = s6i(i) + xs6i
            s7i(i) = s7i(i) + xs7i
          endif
        enddo

      enddo

! MS: Neglected term, divergence of ExB drift
! Divergence of the ExB drift; requires equatorial drift

      nzh = (nz+1)/2
      vexbeq = vexb(nzh,nfl)
      do i = 1,nz
        divvexb(i) = 6.*vexbeq /
     .               (ps(i,nfl)*re*1.e5) *
     .               cos(blats(i,nfl)*po180)**2 *
     .               (1.+sin(blats(i,nfl)*po180)**2) /
     .               (1.+3.*sin(blats(i,nfl)*po180)**2)**2
        s2i(i) = s2i(i) - 0.6667 * divvexb(i)
      enddo

      call tisolv(tti,tiold,kapi,s1i,s2i,s3i,s4i,s5i,s6i,s7i,ptop,nfl)

      return
      end



*******************************************
*******************************************

!            etemp

*******************************************
*******************************************

      subroutine etemp ( tte,te_old,phprodr,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'


      real tte(nz),te_old(nz),kape(nz)
      real s1e(nz),s2e(nz),s3e(nz),s4e(nz),phprodr(nz,nion)
      real s5e(nz),qphe(nz),phprod(nz)
      real qen(nz,nneut)
      real ratio(nz)
      real divvexb(nz)
      integer iz300s(nf),iz300n(nf)

      do i = 1,nz
        s1e(i)  = 0.
        s2e(i)  = 0.
        s3e(i)  = 0.
        s4e(i)  = 0.
        kape(i) = 0.
        do ni = 1,nneut
          qen(i,ni) = 0.
        enddo
      enddo

      do i = 1,nz

        fac1 = denn(i,nfl,pto)  * 1.1e-16
     .          * ( 1. + 5.7e-4 * te(i,nfl) )
        fac2 = denn(i,nfl,ptn2) * 2.82e-17 * sqrt(te(i,nfl))
     .          * ( 1  - 1.2e-4 * te(i,nfl) )
        fac3 = denn(i,nfl,pto2) * 2.2e-16
     .         * ( 1. + 3.6e-2  * sqrt(te(i,nfl)) )
        akpefac = fac1 + fac2 + fac3

        kape(i) = 7.7e5 * sqrt ( te(i,nfl)**5 ) * 0.6667 * evtok
     .      / ( 1. + 3.22e4 * ( te(i,nfl)**2 / ne(i,nfl) * akpefac) )


!       neutrals (Tn - Te) term

!       N2

!       vibrational state from red book (p. 269) milward et al.
!       removed (2/16/01)

        qen(i,ptn2) = .6667 *  evtok * denn(i,nfl,ptn2) *
     .                  ( 1.2e-19 * ( 1. - 1.2e-4 * te(i,nfl) )
     .                            * te(i,nfl) +
     .                    2.e-14 / sqrt(te(i,nfl))
     .                    + 6.5e-22 * ( tn(i,nfl) - 310 ) ** 2 *
     .                      exp(.0023*(te(i,nfl) - tn(i,nfl))) )

!       O2

        qen(i,pto2) = .6667 * evtok * denn(i,nfl,pto2) *
     .                 ( 7.9e-19 * ( 1. + 3.6e-2 * sqrt(te(i,nfl)) ) *
     .                     sqrt(te(i,nfl)) +
     .                   7.e-14 / sqrt(te(i,nfl)) )

!       O

        qen(i,pto) = .6667 * 7.2e-18 * evtok * denn(i,nfl,pto) *
     .                  sqrt(te(i,nfl))

!       H

        qen(i,pth) = .6667 * 6.3e-16 * evtok * denn(i,nfl,pth) *
     .                  ( 1. - 1.35e-4 * te(i,nfl) ) *
     .                  sqrt(te(i,nfl))

        do nn = 1,nneut
          s2e(i) = s2e(i) + qen(i,nn)
        enddo

        s1e(i) = s2e(i) * tn(i,nfl)

!       ions (Ti - Te) term

        do ni = nion1,nion2
          xs3e    = 7.7e-6 * deni(i,nfl,ni) / ami(ni)
     .                     / te(i,nfl) / sqrt(te(i,nfl))
     .                     * .66667 * evtok
          xs4e    = xs3e * ti(i,nfl,ni)
          s3e(i) = s3e(i) + xs3e
          s4e(i) = s4e(i) + xs4e
        enddo

      enddo

!     photoelectron heating
!     red book (millward et al. p. 269)

!     calculate total ion photoproduction (= photoelectron)

      do i = 1,nz
        phprod(i)   = 0.
        do ni = nion1,nion2
          phprod(i) = phprodr(i,ni) * denn(i,nfl,ni) + phprod(i)
        enddo
      enddo

! comment out: use iz300s/n calculated in grid

! iz300s/iz300n are redefined here

      do i = 1,nz
        ratio(i) = ne(i,nfl) /
     .             (0.1*denn(i,nfl,pto)+
     .              denn(i,nfl,pto2)+denn(i,nfl,ptn2))
      enddo

      i = 1
      do while ( ratio(i) .le. 3.e-3 .and. i .lt. nz )
         iz300s(nfl) = i
         i         = i + 1
      enddo

      i = nz
      do while ( ratio(i) .le. 3.e-3 .and. i .gt. 1 )
         iz300n(nfl) = i
         i         = i - 1
      enddo

      if ( iz300s(nfl) .gt. iz300n(nfl) ) then

        do i = 1,nz
            xarg =   ne(i,nfl)
     .             / (        denn(i,nfl,pto2)
     .                 +      denn(i,nfl,ptn2)
     .                 + .1 * denn(i,nfl,pto)   )
            x    = alog ( xarg )
            earg =     12.75
     .               + 6.941 * x
     .               + 1.166 * x ** 2
     .               + 0.08034 * x ** 3
     .               + 0.001996 * x ** 4
            epsi = exp ( -earg )
            qphe(i) = epsi * phprod(i)
          enddo
      else
          do i = 1,iz300s(nfl)
            xarg =   ne(i,nfl)
     .             / (        denn(i,nfl,pto2)
     .                 +      denn(i,nfl,ptn2)
     .                 + .1 * denn(i,nfl,pto)   )
            x    = alog ( xarg )
            earg =     12.75
     .               + 6.941 * x
     .               + 1.166 * x ** 2
     .               + 0.08034 * x ** 3
     .               + 0.001996 * x ** 4
            epsi = exp ( -earg )
            qphe(i) = epsi * phprod(i)
          enddo

!       smooth things at 300 km

        izs       = iz300s(nfl)
!        facts = (250.-alts(izs,nfl)) /
!     .             (alts(izs+1,nfl)-alts(izs,nfl))
        facts = (3.e-3-ratio(izs)) /
     .          (ratio(izs+1)-ratio(izs))
        ne300s = ne(izs,nfl) + (ne(izs+1,nfl)-ne(izs,nfl)) * facts
        o2300 = denn(izs,nfl,pto2) +
     .         (denn(izs+1,nfl,pto2)-denn(izs,nfl,pto2)) * facts
        n2300 = denn(izs,nfl,ptn2) +
     .         (denn(izs+1,nfl,ptn2)-denn(izs,nfl,ptn2)) * facts
        o300 = denn(izs,nfl,pto) +
     .         (denn(izs+1,nfl,pto)-denn(izs,nfl,pto)) * facts
        phprod300 = phprod(izs) +
     .        (phprod(izs+1)-phprod(izs)) * facts
        xarg300 = ne300s / ( o2300 + n2300 + 0.1*o300 )
        x300 = alog( xarg300)
        earg300 =     12.75 +
     .        6.941 * x300 +
     .        1.166 * x300 ** 2 +
     .        0.08034 * x300 ** 3 +
     .        0.001996 * x300 ** 4
        epsi300 = exp ( -earg300 )
        q0s = epsi300 * phprod300 / ne300s

        do i = iz300n(nfl),nz
          xarg =   ne(i,nfl)
     .           / (       denn(i,nfl,pto2)
     .              +      denn(i,nfl,ptn2)
     .              + .1 * denn(i,nfl,pto) )
          x    = alog ( xarg )
          earg =     12.75
     .             + 6.941 * x
     .             + 1.166 * x ** 2
     .             + 0.08034 * x ** 3
     .             + 0.001996 * x ** 4
          epsi = exp ( -earg )
          qphe(i) = epsi * phprod(i)
        enddo

        izn      = iz300n(nfl)
!        factn = (250.-alts(izn,nfl)) /
!     .             (alts(izn-1,nfl)-alts(izn,nfl))
        factn = (3.e-3-ratio(izn)) /
     .           (ratio(izn-1)-ratio(izn))
        ne300n = ne(izn,nfl) +
     .        (ne(izn-1,nfl)-ne(izn,nfl)) * factn
        o2300 = denn(izn,nfl,pto2) +
     .        (denn(izn-1,nfl,pto2)-denn(izn,nfl,pto2)) * factn
        n2300 = denn(izn,nfl,ptn2) +
     .        (denn(izn-1,nfl,ptn2)-denn(izn,nfl,ptn2)) * factn
        o300 = denn(izn,nfl,pto) +
     .        (denn(izn-1,nfl,pto)-denn(izn,nfl,pto)) * factn
        phprod300 = phprod(izn) +
     .        (phprod(izn-1)-phprod(izn)) * factn
        xarg300 = ne300n / ( o2300 + n2300 + 0.1*o300 )
        x300 = alog( xarg300)
        earg300 =     12.75 +
     .        6.941 * x300 +
     .        1.166 * x300 ** 2 +
     .        0.08034 * x300 ** 3 +
     .        0.001996 * x300 ** 4
        epsi300 = exp ( -earg300 )
        q0n = epsi300 * phprod300 / ne300n

        xbms = bms(izs,nfl) + (bms(izs+1,nfl)-bms(izs,nfl)) * facts
        xbmn = bms(izn,nfl) + (bms(izn-1,nfl)-bms(izn,nfl)) * factn


        dels300s = dels(iz300s(nfl),nfl) * facts
        dels300n = dels(iz300n(nfl)-1,nfl) * factn

        ! MS: Old code used a wasteful way to calculate xn.
        ! Cleaner version here.
        xn = 0.
        ! Set bottom integration bound to 300 km.
        xn =   xn + 0.5 * ( ne(iz300n(nfl)-1,nfl) + ne300n ) *
     .        (dels(iz300n(nfl)-1,nfl) - dels300n )
        do i =iz300n(nfl)-2,iz300s(nfl)+1,-1
           xn = xn + 0.5 * ( ne(i,nfl) + ne(i+1,nfl) ) * dels(i,nfl)
        enddo

!        cqe   = 6.e-14                     ! constant (now in namelist)

        if ( q0s .lt. 0 .or. q0n .lt. 0 ) then
          print *,' q0s = ',q0s,' q0n = ',q0n,' nfl = ',nfl
        endif

!       1/22/00

!       put in dels (arc length along field line)

        xs    = 0.

        do i = iz300s(nfl)+1,iz300n(nfl)-1
           if (i .eq. iz300s(nfl)+1) then
             xs = xs + 0.5*( ne300s + ne(i,nfl) ) *
     .           (dels(iz300s(nfl),nfl) - dels300s)
           else
              xs = xs + 0.5 * ( ne(i,nfl) + ne(i-1,nfl) )
     .                         * dels(i-1,nfl)
              xn = xn - 0.5 * ( ne(i,nfl) + ne(i-1,nfl) )
     .                         * dels(i-1,nfl)
           endif

           xints = cqe*xs
           xintn = cqe*xn
           xqs    = ne(i,nfl) * q0s * bms(i,nfl) / xbms * exp(-xints)
           xqn    = ne(i,nfl) * q0n * bms(i,nfl) / xbmn * exp(-xintn)
           qphe(i) = xqs + xqn
        enddo

      endif

      do i = 1,nz
        s5e(i) = 0.66667 * evtok * qphe(i) / ne(i,nfl) ! * .15
      enddo

! MS: Neglected term, divergence of ExB drift
! Divergence of the ExB drift; requires equatorial drift

      nzh    = (nz+1)/2
      vexbeq = vexb(nzh,nfl)
      do i = 1,nz
        divvexb(i) = 6.*vexbeq /
     .                 (ps(i,nfl)*re*1.e5) *
     .                 cos(blats(i,nfl)*po180)**2 *
     .                 (1.+sin(blats(i,nfl)*po180)**2) /
     .                 (1.+3.*sin(blats(i,nfl)*po180)**2)**2
        s2e(i) = s2e(i) - 0.6667 * divvexb(i)
      enddo

      call tesolv(tte,te_old,kape,s1e,s2e,s3e,s4e,s5e,nfl)

      return
      end


*******************************************
*******************************************

!            densolv2

*******************************************
*******************************************

      subroutine densolv2( ni,tdeni,prod,loss,oldion,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'

      real tdeni(nz)
      real oldion(nz), prod(nz), loss(nz)
      real a(nz), b(nz), c(nz), d(nz)

!     initialize

      do j = 1,nz
        a(j) = 0.
        b(j) = 0.
        c(j) = 0.
        d(j) = 0.
      enddo


      do j = 2,nz-1

      ujm1  = vsi(j-1,nfl,ni)/bms(j-1,nfl)
      uj    = vsi(j,nfl,ni)  /bms(j,nfl)
      ujp1  = vsi(j+1,nfl,ni)/bms(j+1,nfl)
      ur = .5*( uj +ujp1)
      ul = .5*( uj +ujm1)

      if (ur .ge. 0. .and. ul .ge. 0.) then
        a0 = -ul
        b0 =  ur
        c0 =  0.
      endif
      if (ur .le. 0. .and. ul .le. 0.) then
        a0 = 0.
        b0 = -ul
        c0 = ur
      endif
      if (ur .ge. 0. .and. ul .le. 0.) then
        a0 = 0.
        b0 = ur - ul
        c0 = 0.
      endif
      if (ur .le. 0. .and. ul .ge. 0.) then
        a0 = -ul
        b0 = 0.
        c0 = ur
      endif

      a(j) =  a0 * bms(j,nfl) ** 2 / d22s(j,nfl)

      b(j) = 1. / dt + loss(j) + b0 * bms(j,nfl) ** 2 / d22s(j,nfl)

      c(j) = c0 * bms(j,nfl) ** 2 / d22s(j,nfl)

      d(j) = oldion(j) / dt + prod(j)

      enddo

!     we will assume that they are determined by the production and loss
!     at both ends of the field line

!     lower bc

      a(1) = 0.
      b(1) = 1.
      c(1) = 0.
      d(1) =
     .  sqrt ( tdeni(1) * prod(1) / loss(1) ) + denmin

!     upper bc

      a(nz) = 0.
      b(nz) = 1.
      c(nz) = 0.
      d(nz) =
     .  sqrt ( tdeni(nz) * prod(nz) / loss(nz) ) + denmin


      call rtryds ( a,b,c,d,tdeni,nz )

      return
      end



*******************************************
*******************************************

!            vsisolv

*******************************************
*******************************************

      subroutine vsisolv ( vi,vid,viold,snuj,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'

      dimension a(nz), b(nz), c(nz), d(nz)
      real vi(nz),vid(nz),viold(nz),snuj(nz)

!     initialize

      do j = 1,nz
        a(j) = 0.
        b(j) = 0.
        c(j) = 0.
        d(j) = 0.
      enddo


      do j = 2,nz-1

        ujm1 = vi(j-1)
        uj   = vi(j)
        ujp1 = vi(j+1)
        ur = .25*( uj +ujp1)
        ul = .25*( uj +ujm1)

        if (ur .ge. 0. .and. ul .ge. 0.) then
          a0 = -ul
          b0 =  ur
          c0 =  0.
        endif
        if (ur .le. 0. .and. ul .le. 0.) then
          a0 = 0.
          b0 = -ul
          c0 = ur
        endif
        if (ur .ge. 0. .and. ul .le. 0.) then
          a0 = 0.
          b0 = ur - ul
          c0 = 0.
        endif
        if (ur .le. 0. .and. ul .ge. 0.) then
          a0 = -ul
          b0 = 0.
          c0 = ur
        endif

        a(j) = a0 / d22s(j,nfl) * bms(j,nfl)

        b(j) = 1/dt + snuj(j) + b0 / d22s(j,nfl) * bms(j,nfl)

        c(j) = c0 / d22s(j,nfl) * bms(j,nfl)

        d(j) = viold(j)/dt + vid(j)

      enddo

!     we will assume that the bc's are the neutral temperature
!     at both ends of the field line

!     lower bc

      a(1) = 0.
      b(1) = 1.
      c(1) = 0.
      d(1) = 0.

!     upper bc

      a(nz) = 0.
      b(nz) = 1.
      c(nz) = 0.
      d(nz) = 0.

      call rtryds(a,b,c,d,vi,nz)

      return
      end



*******************************************
*******************************************

!            tisolv

*******************************************
*******************************************

      subroutine tisolv(tti,tio,kap,s1,s2,s3,s4,s5,s6,s7,npt,nfl)

      include 'param-1.00.inc'
      include 'com-1.00.inc'

      real a(nz),b(nz),c(nz),d(nz)
      real s1(nz),s2(nz),s3(nz),tti(nz),tio(nz),kap(nz)
      real s4(nz),s5(nz),s6(nz),s7(nz)

!     initialize

      do j = 1,nz
        a(j) = 0.
        b(j) = 0.
        c(j) = 0.
        d(j) = 0.
      enddo


      do j = 2,nz-1
        ujm1 = bms(j-1,nfl)*vsi(j-1,nfl,npt)
        uj   = bms(j,nfl)  *vsi(j,nfl,npt)
        ujp1 = bms(j+1,nfl)*vsi(j+1,nfl,npt)
        ur = .5*( uj +ujp1)
        ul = .5*( uj +ujm1)

        if (ur .ge. 0. .and. ul .ge. 0.) then
          a0 = -ul
          b0 =  ur
          c0 =  0.
        endif
        if (ur .le. 0. .and. ul .le. 0.) then
          a0 = 0.
          b0 = -ul
          c0 = ur
        endif
        if (ur .ge. 0. .and. ul .le. 0.) then
          a0 = 0.
          b0 = ur - ul
          c0 = 0.
        endif
        if (ur .le. 0. .and. ul .ge. 0.) then
          a0 = -ul
          b0 = 0.
          c0 = ur
        endif

        a(j) =     a0 / d22s(j,nfl)
     .         - ( bms(j,nfl)**2 / deni(j,nfl,npt) ) / d22s(j,nfl)
     .           *.5 * ( kap(j) + kap(j-1) ) / ds(j,nfl)

        b(j) = 1. / dt + b0 / d22s(j,nfl)
     .         -.333333 * ( bms(j,nfl)
     .                     * (vsi(j+1,nfl,npt) - vsi(j-1,nfl,npt) )
     .                     + 5. * vsi(j,nfl,npt)
     .                          * (bms(j+1,nfl) - bms(j-1,nfl) ) )
     .         / d2s(j,nfl)
     .         +  ( bms(j,nfl)**2 / deni(j,nfl,npt) ) / d22s(j,nfl)
     .           *(.5* (kap(j+1) + kap(j) ) / ds(j+1,nfl)
     .         +.5 * (kap(j) + kap(j-1) ) / ds(j,nfl))
     .         + s2(j) + s4(j) + s6(j)

        c(j) =     c0 / d22s(j,nfl)
     .         - ( bms(j,nfl)**2 / deni(j,nfl,npt) ) /d22s(j,nfl)
     .           *.5 * (kap(j+1) + kap(j) ) / ds(j+1,nfl)

        d(j) = tio(j)/dt + s1(j) + s3(j) + s5(j) + s7(j)

      enddo

!     we will assume that the bc's are the neutral temperature
!     at both ends of the field line

!     lower bc

      a(1) = 0.
      b(1) = 1.
      c(1) = 0.
      d(1) = tn(1,nfl)

!     upper bc

      a(nz) = 0.
      b(nz) = 1.
      c(nz) = 0.
      d(nz) = tn(nz,nfl)

      call rtryds ( a,b,c,d,tti,nz )


      return
      end

*******************************************
*******************************************

!            tesolv

*******************************************
*******************************************

      subroutine tesolv(tte,te_old,kap,s1,s2,s3,s4,s5,nfl)

      include 'param-1.00.inc'
      include 'com-1.00.inc'

      dimension a(nz),b(nz),c(nz),d(nz)
      dimension s1(nz),s2(nz),s3(nz),s4(nz),s5(nz)
      real kap(nz),te_old(nz),tte(nz)

!     initialize

      do j = 1,nz
        a(j) = 0.
        b(j) = 0.
        c(j) = 0.
        d(j) = 0.
      enddo

!     note: ne used here is in a common block

      do j = 2,nz-1

        a(j) = - bms(j,nfl)**2 / ne(j,nfl) / d22s(j,nfl)
     .         *.5 * ( kap(j) + kap(j-1) ) / ds(j,nfl)

        b(j) = 1. / dt + bms(j,nfl)**2 / ne(j,nfl) / d22s(j,nfl)
     .        *(  .5 * (kap(j+1) + kap(j)   ) /ds(j+1,nfl)
     .           +.5 * (kap(j)   + kap(j-1) ) /ds(j,nfl)   )
     .        + s2(j) + s3(j)

        c(j) = - bms(j,nfl)**2 / ne(j,nfl) /d22s(j,nfl)
     .         *.5 * ( kap(j+1) + kap(j) )/ ds(j+1,nfl)

        d(j) = te_old(j)/dt + s1(j) + s4(j) + s5(j)

       enddo

!     we will assume that the bc's are the neutral temperature
!     at both ends of the field line

!     lower bc

      a(1) = 0.
      b(1) = 1.
      c(1) = 0.
      d(1) = tn(1,nfl)

!     upper bc

      a(nz) = 0.
      b(nz) = 1.
      c(nz) = 0.
      d(nz) = tn(nz,nfl)

      call rtryds(a,b,c,d,tte,nz)

      return
      end

*******************************************
*******************************************

!            rtryds

*******************************************
*******************************************

      subroutine rtryds(a,b,c,d,x,n)

      include 'param-1.00.inc'

!     arrays a,b,c, and d may be used for storage of alfa, beta and x
!     in the actual call of this routine, but remember, whatever you
!     use will be lost by the definition of of alfa and beta here.
!     form,  a(k)*x(k-1) + b(k)*x(k) + c(k)*x(k+1) = d(k)

!     i have modified the input sequence to the routine, but have left it
!     otherwise intact.  we may  want to eventually change this (gj)

      dimension a(nz),b(nz),c(nz),d(nz),x(nz)
      dimension alfa(nz),beta(nz)

      nm1=n-1

!     apply the boundary condition at x(1)
!     alfa(1) and beta(1) determined from b(1),c(1),d(1),a(1)=0.

      dst     = d(1)
      rb      = 1. / b(1)
      alfa(1) = -c(1) * rb
      beta(1) =   dst * rb

!     calculate the alfas and betas of k on forward sweep

      do k=2,nm1
        ast     =  a(k)
        z       =  1. / ( b(k) + ast * alfa(k-1) )
        dst     =  d(k)
        alfa(k) = -c(k) * z
        beta(k) =  ( dst - ast * beta(k-1) ) * z
      enddo

!     apply the boundary condition at x(n)
!     x(n) determined from a(n),b(n),d(n),c(n)=0.

      x(n) = ( d(n) - a(n) *beta(nm1) ) / ( b(n) + a(n) * alfa(nm1) )

!     calculate x of k from the alfas and betas on backward sweep

      do i=2,n
        k    = n + 1 - i
        x(k) = x(k+1) * alfa(k) + beta(k)
      enddo

      return
      end

*******************************************
*******************************************

!            msistim

*******************************************
*******************************************


      subroutine msistim ( iyr,iday,hrl,glong,iyd,secut )

!     msistim calculates time parameters for the
!     nrlmsise00 neutral atmosphere model.

!     the arguments are defined as follows:

!       iyr    the julian year
!       iday   the day of the year
!       hr     the local time in hours
!       glong  the geocentric longitude in degrees east
!       iyd    the year and day in the form yydd
!       secut  the universal time in seconds

      iyd    = 1000 * mod(iyr,100) + iday
      hrut   = hrl - glong /15.

      do while ( hrut .lt. 0.  )
        hrut = hrut + 24.
      enddo

      do while ( hrut. ge. 24. )
        hrut = hrut - 24.
      enddo

      secut  = hrut * 3600.

      return
      end

*******************************************
*******************************************

!            zenith

*******************************************
*******************************************

       subroutine zenith (hrut,nfl)

       include 'param-1.00.inc'
       include 'com-1.00.inc'

!      geometric variables

!      sdec: solar zenith angle
!      cx:  cos of the zenith angle

       do i = 1,nz
         hrl = hrut + glons(i,nfl) / 15.
         if ( hrl .ge. 24. ) hrl = hrl - 24.
         sdec         = rtod * asin (  sin (2.*pie*(day-dayve)/sidyr)
     .                               * sin (solinc/rtod)             )
         cossdec      = cos ( po180 * sdec )
         sinsdec      = sin ( po180 * sdec )
         clat         = cos ( po180 * glats(i,nfl) )
         slat         = sin ( po180 * glats(i,nfl) )
         cx(i,nfl)    =   slat * sinsdec
     .                  - clat * cossdec * cos ( 15.0*po180*hrl )
!         u3(i,nfl)    = cx(i,nfl)
! MS: Since we will be taking acos of this value in photprod, make
! sure that the absolute value does not minutely exceed 1 because of
! round-off error.
        if (abs(abs(cx(i,nfl))-1.) .lt. 1.e-6)
     .     cx(i,nfl) = sign(1.,cx(i,nfl))
       enddo

       return
       end


*******************************************
*******************************************

!             xerfcexp

*******************************************
*******************************************

        function xerfcexp(x)

        include 'param-1.00.inc'

          t          = 1. / (1 + pas * x)
          xerfcexp   = (   z1 * t
     .                   + z2 * t ** 2
     .                   + z3 * t ** 3
     .                   + z4 * t ** 4
     .                   + z5 * t ** 5  )

        return
        end

*******************************************
*******************************************

!            f1026

*******************************************
*******************************************

!     subroutine to calculate the nighttime flux of
!     lyman beta (1026) (note: line = 1)

      subroutine sf1026 ( f,line,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'

      real f(nz,nf,91)

      imax = 1

!     determine f for the 4 known values of theta

      do i = 1,nz
        if ( alts(i,nfl) .lt. zaltnt(line,1) ) then
          do k = 1,4
            f( i,nfl,int(thetant(line,k))+1-90 ) = 1.
          enddo
        elseif ( zaltnt(line,1) .le. alts(i,nfl) .and.
     .           alts(i,nfl) .le. zaltnt(line,2)       ) then
          f( i,nfl,int(thetant(line,1))+1-90 ) =
     .       1.4e8 * tanh ( (alts(i,nfl) - 90.) / 50. )
          f( i,nfl,int(thetant(line,2))+1-90 ) =
     .       3.8e7 * tanh ( (alts(i,nfl) - 90.) / 50. )
          f( i,nfl,int(thetant(line,3))+1-90 ) =
     .       1.4e7 * tanh ( (alts(i,nfl) - 93.) / 55. )
          f( i,nfl,int(thetant(line,4))+1-90 ) =
     .       9.2e6 * tanh ( (alts(i,nfl) - 94.) / 55. )
          imax = i
        else
          do k = 1,4
            f( i,nfl,   int(thetant(line,k))+1-90 ) =
     .      f( imax,nfl,int(thetant(line,k))+1-90 )
          enddo
        endif
      enddo

      do k = 1,4
        do i = 1,nz
          f( i,nfl,int(thetant(line,k))+1-90 ) =
     .    amax1 ( 1., f( i,nfl,int(thetant(line,k))+1-90 ) )
        enddo
      enddo

!     now interpolate to all valuse of theta (90 - 180)

      do k = 1,91
        k90 = 90 + k - 1
        ji  = 1
        ki  = int(thetant(line,1))
        do j = 1,4
          if ( k90 .gt. int(thetant(line,j)) ) then
            ji = j
            ki = int(thetant(line,ji))
          endif
        enddo
        jip1 = ji + 1
        kip1 = int(thetant(line,jip1))
        delk = float (   int(thetant(line,jip1))
     .                 - int(thetant(line,ji  )) )
        do i = 1,nz
          flog =   alog10(f(i,nfl,ki+1-90))
     .           + (k90 - ki) / delk
     .                        * (  alog10(f(i,nfl,kip1+1-90))
     .                           - alog10(f(i,nfl,ki  +1-90)) )
          f(i,nfl,k) = 10 ** flog
        enddo
      enddo

      return
      end

*******************************************
*******************************************

!            f584

*******************************************
*******************************************

!     subroutine to calculate the nighttime flux of
!     he i (584) (note: line = 2)

      subroutine sf584 ( f,line,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'

      real f(nz,nf,91)

      imax = 1

!     determine f for the 4 known values of theta

      do i = 1,nz
        if ( alts(i,nfl) .lt. zaltnt(line,1) ) then
          do k = 1,4
            f( i,nfl,int(thetant(line,k))+1-90 ) = 1.
          enddo
        elseif ( zaltnt(line,1) .le. alts(i,nfl) .and.
     .           alts(i,nfl) .le. zaltnt(line,2)       ) then
          f( i,nfl,int(thetant(line,1))+1-90 ) =
     .       1.85e5 * ( alts(i,nfl) - 170. ) ** 1.20
          f( i,nfl,int(thetant(line,2))+1-90 ) =
     .       2.60e4 * ( alts(i,nfl) - 170. ) ** 1.25
          f( i,nfl,int(thetant(line,3))+1-90 ) =
     .       2.60e3 * ( alts(i,nfl) - 170. ) ** 1.20
          f( i,nfl,int(thetant(line,4))+1-90 ) =
     .       2.60e2 * ( alts(i,nfl) - 170. ) ** 1.20
          imax = i
        else
          do k = 1,4
            f( i   ,nfl,int(thetant(line,k))+1-90 ) =
     .      f( imax,nfl,int(thetant(line,k))+1-90 )
          enddo
        endif
      enddo

      do k = 1,4
        do i = 1,nz
          f( i,nfl,int(thetant(line,k))+1-90 ) =
     .    amax1 ( 1., f( i,nfl,int(thetant(line,k))+1-90 ) )
        enddo
      enddo

!     now interpolate to all valuse of theta (90 - 180)
!     set f(i,nfl,theta=180) = 1.

      do k = 1,91
        k90 = 90 + k - 1
        ji  = 1
        ki  = int(thetant(line,1))
        do j = 1,4
          if ( k90 .gt. int(thetant(line,j)) ) then
            ji = j
            ki = int(thetant(line,ji))
          endif
        enddo
        if ( ji .ne. 4 ) then
          jip1 = ji + 1
          kip1 = int(thetant(line,jip1))
          delk = float (   int(thetant(line,jip1))
     .                   - int(thetant(line,ji  )) )
          do i = 1,nz
            flog =   alog10(f(i,nfl,ki+1-90))
     .             + (k90 - ki) / delk
     .                          * (  alog10(f(i,nfl,kip1+1-90))
     .                             - alog10(f(i,nfl,ki  +1-90)) )
            f(i,nfl,k) = 10 ** flog
          enddo
        else
          delk = float (   180
     .                   - int(thetant(line,ji  )) )
          do i = 1,nz
            flog =   alog10(f(i,nfl,ki+1-90))
     .             + (k90 - ki) / delk
     .                          * (  alog10(1.)
     .                             - alog10(f(i,nfl,ki  +1-90)) )
            f(i,nfl,k) = 10 ** flog
          enddo
        endif
      enddo

      return
      end

*******************************************
*******************************************

!            f304

*******************************************
*******************************************

!     subroutine to calculate the nighttime flux of
!     he ii (304) (note: line = 3)

      subroutine sf304 ( f,line,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'

      real f(nz,nf,91)

      imax = 1

!     determine f for the 4 known values of theta

      do i = 1,nz
        if ( alts(i,nfl) .lt. zaltnt(line,1) ) then
          do k = 1,4
            f( i,nfl,int(thetant(line,k))+1-90 ) = 1.
          enddo
        elseif ( zaltnt(line,1) .le. alts(i,nfl) .and.
     .           alts(i,nfl) .le. zaltnt(line,2)       ) then
          f( i,nfl,int(thetant(line,1))+1-90 ) =
     .       3.8e6 * tanh ( (alts(i,nfl) - 138.) / 80. )
          f( i,nfl,int(thetant(line,2))+1-90 ) =
     .       3.0e6 * tanh ( (alts(i,nfl) - 138.) / 80. )
          f( i,nfl,int(thetant(line,3))+1-90 ) =
     .       2.5e6 * tanh ( (alts(i,nfl) - 138.) / 80. )
          f( i,nfl,int(thetant(line,4))+1-90 ) =
     .       2.5e6 * tanh ( (alts(i,nfl) - 138.) / 80. )
          imax = i
        else
          do k = 1,4
            f( i,   nfl,int(thetant(line,k))+1-90 ) =
     .      f( imax,nfl,int(thetant(line,k))+1-90 )
          enddo
        endif
      enddo

      do k = 1,4
        do i = 1,nz
          f( i,nfl,int(thetant(line,k))+1-90 ) =
     .    amax1 ( 1., f( i,nfl,int(thetant(line,k))+1-90 ) )
        enddo
      enddo

!     now interpolate to all valuse of theta (90 - 180)
!     set f(i,nfl,theta=180) = 1.

      do k = 1,91
        k90 = 90 + k - 1
        ji  = 1
        ki  = int(thetant(line,1))
        do j = 1,4
          if ( k90 .gt. int(thetant(line,j)) ) then
            ji = j
            ki = int(thetant(line,ji))
          endif
        enddo
        if ( ji .ne. 4 ) then
          jip1 = ji + 1
          kip1 = int(thetant(line,jip1))
          delk = float (   int(thetant(line,jip1))
     .                   - int(thetant(line,ji  )) )
          do i = 1,nz
            flog =   alog10(f(i,nfl,ki+1-90))
     .             + (k90 - ki) / delk
     .                          * (  alog10(f(i,nfl,kip1+1-90))
     .                             - alog10(f(i,nfl,ki  +1-90)) )
            f(i,nfl,k) = 10 ** flog
          enddo
        else
          delk = float (   180
     .                   - int(thetant(line,ji  )) )
          do i = 1,nz
            flog =   alog10(f(i,nfl,ki+1-90))
     .             + (k90 - ki) / delk
     .                          * (  alog10(1.)
     .                             - alog10(f(i,nfl,ki  +1-90)) )
            f(i,nfl,k) = 10 ** flog
          enddo
        endif
      enddo

      return
      end

*******************************************
*******************************************

!            f1216

*******************************************
*******************************************

!     subroutine to calculate the nighttime flux of
!     lyman alpha (1216) (note: line = 4)

      subroutine sf1216 ( f,line,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'

      real f(nz,nf,91)

      imax = 1

!     determine f for the 4 known values of theta

      do i = 1,nz
        if ( alts(i,nfl) .lt. zaltnt(line,1) ) then
          do k = 1,4
            f( i,nfl,int(thetant(line,k))+1-90 ) = 1.
          enddo
        elseif ( zaltnt(line,1) .le. alts(i,nfl) .and.
     .           alts(i,nfl) .le. zaltnt(line,2)       ) then
          f( i,nfl,int(thetant(line,1))+1-90 ) =
     .       1.2e10 * tanh ( (alts(i,nfl) - 80.) / 50. ) + 3.e9
          f( i,nfl,int(thetant(line,2))+1-90 ) =
     .       4.0e9  * tanh ( (alts(i,nfl) - 80.) / 50. ) + 1.e9
          f( i,nfl,int(thetant(line,3))+1-90 ) =
     .       2.0e9  * tanh ( (alts(i,nfl) - 65.) / 50. ) + 1.e8
          f( i,nfl,int(thetant(line,4))+1-90 ) =
     .       1.5e9  * tanh ( (alts(i,nfl) - 75.) / 50. ) + 1.e8
          imax = i
        else
          do k = 1,4
            f( i,   nfl,int(thetant(line,k))+1-90 ) =
     .      f( imax,nfl,int(thetant(line,k))+1-90 )
          enddo
        endif
      enddo

      do k = 1,4
        do i = 1,nz
          f( i,nfl,int(thetant(line,k))+1-90 ) =
     .    amax1 ( 1., f( i,nfl,int(thetant(line,k))+1-90 ) )
        enddo
      enddo

!     now interpolate to all valuse of theta (90 - 180)

      do k = 1,91
        k90 = 90 + k - 1
        ji  = 1
        ki  = int(thetant(line,1))
        do j = 1,4
          if ( k90 .gt. int(thetant(line,j)) ) then
            ji = j
            ki = int(thetant(line,ji))
          endif
        enddo
        jip1 = ji + 1
        kip1 = int(thetant(line,jip1))
        delk = float (   int(thetant(line,jip1))
     .                 - int(thetant(line,ji  )) )
        do i = 1,nz
          flog =   alog10(f(i,nfl,ki+1-90))
     .           + (k90 - ki) / delk
     .                        * (  alog10(f(i,nfl,kip1+1-90))
     .                           - alog10(f(i,nfl,ki  +1-90)) )
          f(i,nfl,k) = 10 ** flog
        enddo
      enddo

      return
      end

*******************************************
*******************************************

!             open_u

*******************************************
*******************************************

      subroutine open_u

      include 'param-1.00.inc'
      include 'com-1.00.inc'
!     open output files (unformatted, except time.dat)

      open ( unit=70, file='time.dat'      ,form='formatted'   )
      open ( unit=71, file='deniu.dat'     ,form='unformatted' )
      open ( unit=72, file='tiu.dat'       ,form='unformatted' )
      open ( unit=73, file='vsiu.dat'      ,form='unformatted' )
      open ( unit=75, file='teu.dat'       ,form='unformatted' )
!      open ( unit=78, file='vnu.dat'       ,form='unformatted' )
!      open ( unit=90, file='vtu.dat'       ,form='unformatted' )
!      open ( unit=91, file='vru.dat'       ,form='unformatted' )
      if (outn) then
          open ( unit=92, file='dennu.dat'     ,form='unformatted' )
      endif
!      open ( unit=93, file='vexbu.dat'     ,form='unformatted' )

!     diagnostic files (unformatted)

!      open ( unit=81, file='t1u.dat'  ,form='unformatted' )
!      open ( unit=82, file='t2u.dat'  ,form='unformatted' )
!      open ( unit=83, file='t3u.dat'  ,form='unformatted' )
!      open ( unit=84, file='u1u.dat'  ,form='unformatted' )
!      open ( unit=85, file='u2u.dat'  ,form='unformatted' )
!      open ( unit=86, file='u3u.dat'  ,form='unformatted' )
      if (outn) then
          open ( unit=87, file='u4u.dat'  ,form='unformatted' )
      endif
!      open ( unit=88, file='u5u.dat'  ,form='unformatted' )

      return
      end

*******************************************
*******************************************

!             open_f

*******************************************
*******************************************

      subroutine open_f

      include 'param-1.00.inc'
      include 'com-1.00.inc'
!     open output files (formatted)

      open ( unit=70, file='time.dat'      ,form='formatted' )
      open ( unit=71, file='denif.dat'     ,form='formatted' )
      open ( unit=72, file='tif.dat'       ,form='formatted' )
      open ( unit=73, file='vsif.dat'      ,form='formatted' )
      open ( unit=75, file='tef.dat'       ,form='formatted' )
!      open ( unit=78, file='vnf.dat'       ,form='formatted' )
!      open ( unit=90, file='vtf.dat'       ,form='formatted' )
!      open ( unit=91, file='vrf.dat'       ,form='formatted' )
      if ( outn ) then
        open ( unit=92, file='dennf.dat'     ,form='formatted' )
      endif
!      open ( unit=93, file='vexbf.dat'     ,form='formatted' )

!     diagnostic files (formatted)

!      open ( unit=81, file='t1f.dat'  ,form='formatted' )
!      open ( unit=82, file='t2f.dat'  ,form='formatted' )
!      open ( unit=83, file='t3f.dat'  ,form='formatted' )
!      open ( unit=84, file='u1f.dat'  ,form='formatted' )
!      open ( unit=85, file='u2f.dat'  ,form='formatted' )
!      open ( unit=86, file='u3f.dat'  ,form='formatted' )
      if ( outn ) then
        open ( unit=87, file='u4f.dat'  ,form='formatted' )
      endif
!      open ( unit=88, file='u5f.dat'  ,form='formatted' )

      return
      end

*******************************************
*******************************************

!             output

*******************************************
*******************************************

       subroutine output ( hrut,ntm,istep )

       include 'param-1.00.inc'
       include 'com-1.00.inc'

       hr24   = mod (hrut,24.)
       totsec = hr24 * 3600.
       thr    = totsec / 3600.
       nthr   = int(thr)
       tmin   = ( thr - nthr ) * 60.
       ntmin  = int(mod(tmin,60.))
       tsec   = ( tmin - ntmin ) * 60.
       ntsec  = int(tsec)

!      put data into output variables

       write (*,*) 'istep = ',istep,'ntm = ',ntm,
     .             'time step = ',dt,' hour = ',hrut
       write (70,100) ntm,nthr,ntmin,ntsec

       if ( fmtout ) then
         write(71,101) deni
         write(72,101) ti
         write(73,101) vsi
         write(75,101) te
!         write(78,101) vn
!         write(81,101) t1
!         write(82,101) t2
!         write(83,101) t3
!         write(84,101) u1
!         write(85,101) u2
!         write(86,101) u3
!         write(87,101) u4
!         write(88,101) u5
!         write(90,101) vot
!         write(91,101) vor
!         write(92,101) denn
       endif
       if ( fmtout .and. outn ) then !added for RTGR (JMS)
          write(87,101) u4
          write(92,101) denn
       endif

       if ( .not. fmtout ) then
         write(71) deni
         write(72) ti
         write(73) vsi
         write(75) te
!         write(78) vn
!         write(81) t1
!         write(82) t2
!         write(83) t3
!         write(84) u1
!         write(85) u2
!         write(86) u3
         if (outn) then
           write(87) u4
!         write(88) u5
!         write(90) vot
!         write(91) vor
           write(92) denn
         endif
!         write(93) vexbp
       endif

 100   format(1x,4i6)
 101   format(1x,1p1e16.6)

       return
       end


*******************************************
*******************************************

!             EXB

*******************************************
*******************************************

       subroutine exb(hrut)

       include 'param-1.00.inc'
       include 'com-1.00.inc'

       real denic(nz,nf,nion)
       real tic(nz,nf,nion)
       real tec(nz,nf)
       real fluxnp(nz,nf,nion),fluxtp(nz,nf,nion)
       real fluxtep(nz,nf)
       real fluxns(nz,nf,nion),fluxts(nz,nf,nion)
       real fluxtes(nz,nf)
       real param(2)

!     define the e x b drift

      param(1) = day
      param(2) = f10p7
      nzh      = ( nz - 1 ) / 2

!     note: modification of vexb because of field line variation
!           uses cos^3/sqrt(1.+3sin^2) instead of
!           uses sin^3/sqrt(1.+3cos^2) because
!           blats = 0 at the magnetic equator
!           (as opposed to pi/2 in spherical coordinates)


      hrl = hrut + glons(nzh,1) / 15.
      if ( hrl .ge. 24. ) hrl = hrl - 24.
      call vdrift_model(hrl,glon_in,param,vd)

      do j = 1,nf
!        altfac = ( alts(nzh,j) + re ) / re ! L^2 dependence on E x B drift
        altfac = 1.
        vexb0 = 100. * vd * tvexb0 * altfac * altfac ! convert to cm/s
        do i = 1,nz
          vexb(i,j) = 0.
          blat = blats(i,j) * po180
          vexb(i,j) = vexb0 *
     .                cos(blat) * cos(blat) * cos(blat)/
     .                sqrt( 1. + 3. * sin(blat)* sin(blat) )
        enddo
      enddo

      do j = 1,nf
        do i = 1,nz
          if ( alts(i,j) .lt. alt_crit ) then
            dela      = 20.
            farg      = ( alt_crit - alts(i,j) ) / dela
            vexb(i,j) = vexb(i,j) * exp (-farg*farg)
          endif
        enddo
      enddo

      do j = 1,nf
        vexb(nzp1,j) = vexb(nz,j)
      enddo
      do i = 1,nzp1
        vexb(i,nfp1) = vexb(i,nf)
      enddo

!     sweep in p-direction

      do j = 1,nfp1
        do i = 1,nz
          vexbp(i,j) = vexb(i,j) * ( vnx(i,j) * xnormp(i,j) +
     .                               vny(i,j) * ynormp(i,j) +
     .                               vnz(i,j) * znormp(i,j)   )
        enddo
      enddo

!     sweep in s-direction

      do j = 1,nf
        do i = 1,nzp1
          vexbs(i,j) = vexb(i,j) * ( vnx(i,j) * xnorms(i,j) +
     .                               vny(i,j) * ynorms(i,j) +
     .                               vnz(i,j) * znorms(i,j)   )
        enddo
      enddo

!     output e x b velocities

      do j = 1,nf
        do i = 1,nz
          u1(i,j) = vexbp(i,j)
          u2(i,j) = vexbs(i,j)
        enddo
      enddo

!      dump velocity
!        vot: velocity in theta direction
!        vor: velocity in radial direction
!      needs to be redone ....


!       do ni = nion1,nion2
!         do j = 1,nf
!           do i = 1,nz
!             vot(i,j,ni)   =  vsi(i,j,ni) * cosdips(i,j) +
!     .                        vexbp(i,j)  * sindips(i,j)
!             vor(i,j,ni)   = -vsi(i,j,ni) * sindips(i,j) +
!     .                        vexbp(i,j)  * cosdips(i,j)
!           enddo
!         enddo
!       enddo

!      calculate conserved particle number: denic
!      and 'conserved' temperature: tic,tec

       do ni = nion1,nion2
         do j = 1,nf
           do i = 1,nz
             denic(i,j,ni) = deni(i,j,ni) * vol(i,j)
             tic(i,j,ni)   = ti(i,j,ni) * vol(i,j)
           enddo
         enddo
       enddo

       do j = 1,nf
         do i = 1,nz
           tec(i,j)   = te(i,j) * vol(i,j)
         enddo
       enddo

!      calculate flux in p-direction at interface

       do ni = nion1,nion2
         do j = 2,nf
           do i = 1,nz
             if ( vexbp(i,j) .ge. 0 ) then
               fluxnp(i,j,ni) = deni(i,j-1,ni) * vexbp(i,j)
               fluxtp(i,j,ni) = ti(i,j-1,ni)   * vexbp(i,j)
             else
               fluxnp(i,j,ni) = deni(i,j,ni) * vexbp(i,j)
               fluxtp(i,j,ni) = ti(i,j,ni)   * vexbp(i,j)
               if ( j .eq. nf ) then
                 fluxnp(i,j,ni) = fluxnp(i,j-1,ni)
     .                            * areap(i,j-1)/areap(i,j)
                 fluxtp(i,j,ni) = fluxtp(i,j-1,ni)
     .                            * areap(i,j-1)/areap(i,j)

                 if ( ni .eq. 1 ) then
                 fluxnp(i,j,ni) = (deni(i,j-1,ni) -
     .                            (areas(i,j-1) / areas(i,j-2)) *
     .                            (deni(i,j-2,ni) - deni(i,j-1,ni)))*
     .                            vexbp(i,j)

                 fluxtp(i,j,ni) = fluxtp(i,j-1,ni) -
     .                            (areas(i,j-1) / areas(i,j-2)) *
     .                            (fluxtp(i,j-2,ni) - fluxtp(i,j-1,ni))

                  endif

!                 fluxnp(i,j,ni) = fluxnp(i,j-1,ni)
!                 fluxtp(i,j,ni) = fluxtp(i,j-1,ni)
               endif
             endif
           enddo
         enddo
       enddo


       do j = 2,nf
         do i = 1,nz
           if ( vexbp(i,j) .ge. 0 ) then
             fluxtep(i,j) = te(i,j-1)   * vexbp(i,j)
           else
             fluxtep(i,j) = te(i,j)   * vexbp(i,j)
             if ( j .eq. nf ) then
                fluxtep(i,j) = fluxtep(i,j-1)
     .                        * areap(i,j-1)/areap(i,j)

                if (ni .eq. 1 )
     .          fluxtep(i,j) = fluxtep(i,j-1) -
     .                         (areas(i,j-1) / areas(i,j-2)) *
     .                         (fluxtep(i,j-2) - fluxtep(i,j-1))

!              fluxtep(i,j) = fluxtep(i,j-1)
             endif
           endif
         enddo
       enddo

!      calculate flux in s-direction at interface

       do ni = nion1,nion2
         do j = 1,nf
           do i = 2,nz
             if ( vexbs(i,j) .ge. 0 ) then
               fluxns(i,j,ni) = deni(i-1,j,ni) * vexbs(i,j)
               fluxts(i,j,ni) = ti(i-1,j,ni)   * vexbs(i,j)
             else
               fluxns(i,j,ni) = deni(i,j,ni) * vexbs(i,j)
               fluxts(i,j,ni) = ti(i,j,ni)   * vexbs(i,j)
             endif
           enddo
         enddo
       enddo

       do j = 1,nf
         do i = 2,nz
           if ( vexbs(i,j) .ge. 0 ) then
             fluxtes(i,j) = te(i-1,j)   * vexbs(i,j)
           else
             fluxtes(i,j) = te(i,j)   * vexbs(i,j)
           endif
         enddo
       enddo

!      update total particle number and density
!      and temperatures

       do ni = nion1,nion2
         do j = 2,nfm1
           do i = 2,nzm1
             denic(i,j,ni) = denic(i,j,ni)
     .                       + dt * ( areap(i,j)   * fluxnp(i,j,ni) -
     .                                areap(i,j+1) * fluxnp(i,j+1,ni) )
     .                       + dt * ( areas(i,j)   * fluxns(i,j,ni) -
     .                                areas(i+1,j) * fluxns(i+1,j,ni) )
             deni(i,j,ni)  = denic(i,j,ni) / vol(i,j)
             tic(i,j,ni) = tic(i,j,ni)
     .                       + dt * ( areap(i,j)   * fluxtp(i,j,ni) -
     .                                areap(i,j+1) * fluxtp(i,j+1,ni) )
     .                       + dt * ( areas(i,j)   * fluxts(i,j,ni) -
     .                                areas(i+1,j) * fluxts(i+1,j,ni) )
             ti(i,j,ni)  = tic(i,j,ni) / vol(i,j)
           enddo
         enddo
       enddo

       do j = 2,nfm1
         do i = 2,nzm1
           tec(i,j) = tec(i,j)
     .                     + dt * ( areap(i,j)   * fluxtep(i,j) -
     .                              areap(i,j+1) * fluxtep(i,j+1) )
     .                     + dt * ( areas(i,j)   * fluxtes(i,j) -
     .                              areas(i+1,j) * fluxtes(i+1,j) )
           te(i,j)  = tec(i,j) / vol(i,j)
         enddo
       enddo

!      fill cells at j = 1 and nf with j = 2 and nfm1

       do ni = nion1,nion2
         do i = 2,nzm1
           deni(i,1,ni)  = deni(i,2,ni)
           deni(i,nf,ni) = deni(i,nfm1,ni)
           ti(i,1,ni)    = ti(i,2,ni)
           ti(i,nf,ni)   = ti(i,nfm1,ni)
         enddo
       enddo

       do i = 2,nzm1
         te(i,1)    = te(i,2)
         te(i,nf)   = te(i,nfm1)
       enddo

!      fill cells at i = 1 and nz with i = 2 and nzm1

       do ni = nion1,nion2
         do j = 1,nf
           deni(1,j,ni)  = deni(2,j,ni)
           deni(nz,j,ni) = deni(nzm1,j,ni)
           ti(1,j,ni)    = ti(2,j,ni)
           ti(nz,j,ni)   = ti(nzm1,j,ni)
         enddo
       enddo

       do j = 1,nf
         te(1,j)    = te(2,j)
         te(nz,j)   = te(nzm1,j)
       enddo

!      E x B drift can cause excessive density
!      if the first field line peak is below alt_crit
!      and the second field line peak is above alt_crit
!      using zero gradient boundary condition -
!      this is a fix

       if ( alts(nz/2,2) .ge. alt_crit .and.
     .      alts(nz/2,1) .lt. alt_crit       ) then
         do ni = nion1,nion2
           do i = 1,nz
             deni(i,1,ni)  = denmin
             ti(i,1,ni)    = tn(i,1)
           enddo
         enddo
         do ie = 1,nz
           te(ie,1)  = tn(ie,1)
         enddo
       endif


       return
       end


*******************************************
*******************************************

!             courant

*******************************************
*******************************************

       subroutine courant ( hrut )

       include 'param-1.00.inc'
       include 'com-1.00.inc'

       hr24ut = mod(hrut,24.)
       dt00   = dt0
       hr24l = hr24ut + glon_in / 15.
       if ( hr24l .ge. 24. ) hr24l = hr24l - 24.
       if ( hr24l .gt. 18. .and. hr24l .le. 24. ) dt00 = .4 * dt0
       if ( hr24l .gt.  0. .and. hr24l .le.  6. ) dt00 = .4 * dt0

!      parallel motion

       dtnew = 1.e6
       do k = nion1,nion2
         do j = 1,nf
           do i = 1,nz
             dt1 = dels(i,j) / amax1(1.,abs(vsi(i,j,k)))
             if ( dt1 .le. dtnew ) dtnew = dt1
           enddo
         enddo
       enddo

!      perpendicular motion

       do j = 1,nf
         do i = 1,nz
           dts = xdels(i,j,1) / amax1(1.,abs(vexbs(i,j)))
           dtp = xdelp(i,j,1) / amax1(1.,abs(vexbp(i,j)))
           dt1 = amin1 ( dts,dtp )
           if ( dt1 .le. dtnew ) dtnew = dt1
         enddo
       enddo

       if ( dtnew .le. .01 ) then
         print *,' Time step too small: dtnew',dtnew
         stop
       elseif ( dtnew .ge. 5e4 ) then
         print *,' Time step too big: dtnew',dtnew
         stop
       endif
       dt = .25 * dtnew
       if ( dtnew/dt .le. 1.0  ) dt = amin1(dt00,dtnew   )
       if ( dtnew/dt .ge. 1.2  ) dt = amin1(dt00,dt * 1.2)

       return
       end


*******************************************
*******************************************

!             vdrift_model

*******************************************
*******************************************

C       ************************************************************
C       ************************************************************

	subroutine vdrift_model(xt,xl,param,y)

C       ************************************************************

C       ************************************************************
C       SUBROUTINE CALCULATES EQUATORIAL VERTICAL DRIFT AS DESCRIBED
C       IN SCHERLIESS AND FEJER, JGR, 104, 6829-6842, 1999
C       ************************************************************

C       INPUT:   XT: SOLAR LOCAL TIME
C                XL: GEOGRAPHIC LONGITUDE (+ EAST)
C
C             PARAM: 2-DIM ARRAY (DOY,F10.7CM)
C                    DOY     :Day of Year has to run from 1 to 365 (366)
C                    F10.7cm : F10.7cm solar flux
C
C       OUTPUT:   Y: EQUATORIAL VERTICAL DRIFT
C
C       JK => modified to include common data => 21 Mar 2012

C       ************************************************************

        include 'param-1.00.inc'
        include 'com-1.00.inc'

!	implicit none

        real dmlt
!        logical fejer   removed, JK

        real param(2),coeff(624),funct(6)
        real coeff1(312),coeff2(312)
	real xt,xl,y
	real bspl4,bspl4_time,bspl4_long

	integer i,j,ind,il,kk
	integer index_t/13/,dim_t/78/
	integer index_l/8/,dim_l/48/
	integer index/104/,dim/624/
	integer nfunc/6/

        data coeff1/
     *  -10.80592, -9.63722,-11.52666,  -0.05716,-0.06288,  0.03564,
     *   -5.80962, -7.86988, -8.50888, -0.05194, -0.05798, -0.00138,
     *    2.09876,-19.99896, -5.11393, -0.05370, -0.06585,  0.03171,
     *  -10.22653, -3.62499,-14.85924, -0.04023, -0.01190, -0.09656,
     *   -4.85180,-26.26264, -6.20501, -0.05342, -0.05174,  0.02419,
     *  -13.98936,-18.10416, -9.30503, -0.01969, -0.03132, -0.01984,
     *  -18.36633,-24.44898,-16.69001,  0.02033, -0.03414, -0.02062,
     *  -20.27621,-16.95623,-36.58234,  0.01445, -0.02044, -0.08297,
     *    1.44450,  5.53004,  4.55166, -0.02356, -0.04267,  0.05023,
     *    5.50589,  7.05381,  1.94387, -0.03147, -0.03548,  0.01166,
     *    3.24165, 10.05002,  4.26218, -0.03419, -0.02651,  0.07456,
     *    7.02218,  0.06708,-11.31012, -0.03252, -0.01021, -0.09008,
     *   -3.47588, -2.82534, -4.17668, -0.03719, -0.01519,  0.06507,
     *   -4.02607,-11.19563,-10.52923, -0.00592, -0.01286, -0.00477,
     *  -11.47478, -9.57758,-10.36887,  0.04555, -0.02249,  0.00528,
     *  -14.19283,  7.86422, -8.76821,  0.05758, -0.02398, -0.04075,
     *   14.58890, 36.63322, 27.57497,  0.01358, -0.02316,  0.04723,
     *   12.53122, 29.38367, 21.40356, -0.00071, -0.00553,  0.01484,
     *   18.64421, 26.27327, 18.32704,  0.00578,  0.03349,  0.11249,
     *    4.53014,  6.15099,  7.41935, -0.02860, -0.00395, -0.08394,
     *   14.29422,  9.77569,  2.85689, -0.00107,  0.04263,  0.10739,
     *    7.17246,  4.40242, -1.00794,  0.00089,  0.01436,  0.00626,
     *    7.75487,  5.01928,  4.36908,  0.03952, -0.00614,  0.03039,
     *   10.25556,  8.82631, 24.21745,  0.05492, -0.02968,  0.00177,
     *   21.86648, 24.03218, 39.82008,  0.00490, -0.01281, -0.01715,
     *   19.18547, 23.97403, 34.44242,  0.01978,  0.01564, -0.02434,
     *   26.30614, 14.22662, 31.16844,  0.06495,  0.19590,  0.05631,
     *   21.09354, 25.56253, 29.91629, -0.04397, -0.08079, -0.07903,
     *   28.30202, 16.80567, 38.63945,  0.05864,  0.16407,  0.07622,
     *   22.68528, 25.91119, 40.45979, -0.03185, -0.01039, -0.01206,
     *   31.98703, 24.46271, 38.13028, -0.08738, -0.00280,  0.01322,
     *   46.67387, 16.80171, 22.77190, -0.13643, -0.05277, -0.01982,
     *   13.87476, 20.52521,  5.22899,  0.00485, -0.04357,  0.09970,
     *   21.46928, 13.55871, 10.23772, -0.04457,  0.01307,  0.06589,
     *   16.18181, 16.02960,  9.28661, -0.01225,  0.14623, -0.01570,
     *   18.16289, -1.58230, 14.54986, -0.00375, -0.00087,  0.04991,
     *   10.00292, 11.82653,  0.44417, -0.00768,  0.15940, -0.01775,
     *   12.15362,  5.65843, -1.94855, -0.00689,  0.03851,  0.04851,
     *   -1.25167,  9.05439,  0.74164,  0.01065,  0.03153,  0.02433,
     *  -15.46799, 18.23132, 27.45320,  0.00899, -0.00017,  0.03385,
     *    2.70396, -0.87077,  6.11476, -0.00081,  0.05167, -0.08932,
     *    3.21321, -1.06622,  5.43623,  0.01942,  0.05449, -0.03084,
     *   17.79267, -3.44694,  7.10702,  0.04734, -0.00945,  0.11516,
     *    0.46435,  6.78467,  4.27231, -0.02122,  0.10922, -0.03331,
     *   15.31708,  1.70927,  7.99584,  0.07462,  0.07515,  0.08934,
     *    4.19893,  6.01231,  8.04861,  0.04023,  0.14767, -0.04308,
     *    9.97541,  5.99412,  5.93588,  0.06611,  0.12144, -0.02124,
     *   13.02837, 10.29950, -4.86200,  0.04521,  0.10715, -0.05465,
     *    5.26779,  7.09019,  1.76617,  0.09339,  0.22256,  0.09222,
     *    9.17810,  5.27558,  5.45022,  0.14749,  0.11616,  0.10418,
     *    9.26391,  4.19982, 12.66250,  0.11334,  0.02532,  0.18919,
     *   13.18695,  6.06564, 11.87835,  0.26347,  0.02858,  0.14801/

        data coeff2/
     *   10.08476,  6.14899, 17.62618,  0.09331,  0.08832,  0.28208,
     *   10.75302,  7.09244, 13.90643,  0.09556,  0.16652,  0.22751,
     *    6.70338, 11.97698, 18.51413,  0.15873,  0.18936,  0.15705,
     *    5.68102, 23.81606, 20.65174,  0.19930,  0.15645,  0.08151,
     *   29.61644,  5.49433, 48.90934,  0.70710,  0.40791,  0.26325,
     *   17.11994, 19.65380, 44.88810,  0.45510,  0.41689,  0.22398,
     *    8.45700, 34.54442, 27.25364,  0.40867,  0.37223,  0.22374,
     *   -2.30305, 32.00660, 47.75799,  0.02178,  0.43626,  0.30187,
     *    8.98134, 33.01820, 33.09674,  0.33703,  0.33242,  0.41156,
     *   14.27619, 20.70858, 50.10005,  0.30115,  0.32570,  0.45061,
     *   14.44685, 16.14272, 45.40065,  0.37552,  0.31419,  0.30129,
     *    6.19718, 18.89559, 28.24927,  0.08864,  0.41627,  0.19993,
     *    7.70847, -2.36281,-21.41381,  0.13766,  0.05113, -0.11631,
     *   -9.07236,  3.76797,-20.49962,  0.03343,  0.08630,  0.00188,
     *   -8.58113,  5.06009, -6.23262,  0.04967,  0.03334,  0.24214,
     *  -27.85742,  8.34615,-27.72532, -0.08935,  0.15905, -0.03655,
     *    2.77234,  0.14626, -4.01786,  0.22338, -0.04478,  0.18650,
     *    5.61364, -3.82235,-16.72282,  0.26456, -0.03119, -0.08376,
     *   13.35847, -6.11518,-16.50327,  0.28957, -0.01345, -0.19223,
     *   -5.37290, -0.09562,-27.27889,  0.00266,  0.22823, -0.35585,
     *  -15.29676,-18.36622,-24.62948, -0.31299, -0.23832, -0.08463,
     *  -23.37099,-13.69954,-26.71177, -0.19654, -0.18522, -0.20679,
     *  -26.33762,-15.96657,-42.51953, -0.13575, -0.00329, -0.28355,
     *  -25.42140,-14.14291,-21.91748, -0.20960, -0.19176, -0.32593,
     *  -23.36042,-23.89895,-46.05270, -0.10336,  0.03030, -0.21839,
     *  -19.46259,-21.27918,-32.38143, -0.17673, -0.15484, -0.11226,
     *  -19.06169,-21.13240,-34.01677, -0.25497, -0.16878, -0.11004,
     *  -18.39463,-16.11516,-19.55804, -0.19834, -0.23271, -0.25699,
     *  -19.93482,-17.56433,-18.58818,  0.06508, -0.18075,  0.02796,
     *  -23.64078,-18.77269,-22.77715, -0.02456, -0.12238,  0.02959,
     *  -12.44508,-21.06941,-19.36011,  0.02746, -0.16329,  0.19792,
     *  -26.34187,-19.78854,-24.06651, -0.07299, -0.03082, -0.03535,
     *  -10.71667,-26.04401,-16.59048,  0.02850, -0.09680,  0.15143,
     *  -18.40481,-23.37770,-16.31450, -0.03989, -0.00729, -0.01688,
     *   -9.68886,-20.59304,-18.46657,  0.01092, -0.07901,  0.03422,
     *   -0.06685,-19.24590,-29.35494,  0.12265, -0.24792,  0.05978,
     *  -15.32341, -9.07320,-13.76101, -0.17018, -0.15122, -0.06144,
     *  -14.68939,-14.82251,-13.65846, -0.11173, -0.14410, -0.07133,
     *  -18.38628,-18.94631,-19.00893, -0.08062, -0.14481, -0.12949,
     *  -16.15328,-17.40999,-14.08705, -0.08485, -0.06896, -0.11583,
     *  -14.50295,-16.91671,-25.25793, -0.06814, -0.13727, -0.12213,
     *  -10.92188,-14.10852,-24.43877, -0.09375, -0.11638, -0.09053,
     *  -11.64716,-14.92020,-19.99063, -0.14792, -0.08681, -0.12085,
     *  -24.09766,-16.14519, -8.05683, -0.24065, -0.05877, -0.23726,
     *  -25.18396,-15.02034,-15.50531, -0.12236, -0.09610, -0.00529,
     *  -15.27905,-19.36708,-12.94046, -0.08571, -0.09560, -0.03544,
     *   -7.48927,-16.00753,-13.02842, -0.07862, -0.10110, -0.05807,
     *  -13.06383,-27.98698,-18.80004, -0.05875, -0.03737, -0.11214,
     *  -13.67370,-16.44925,-16.12632, -0.07228, -0.09322, -0.05652,
     *  -22.61245,-21.24717,-18.09933, -0.05197, -0.07477, -0.05235,
     *  -27.09189,-21.85181,-20.34676, -0.05123, -0.05683, -0.07214,
     *  -27.09561,-22.76383,-25.41151, -0.10272, -0.02058, -0.16720/

        do i = 1,312
          coeff(i)     = coeff1(i)
          coeff(i+312) = coeff2(i)
        enddo

        xt = mod(xt,24.)

!       sinusoid e x b model

!        vpre = 100.0
!        dpre = 0.25

        if ( .not. fejer ) then
          y = ve01
          do i = 1,10
                y = y + fourierA(i)*cos(i*xt*pie/12) + fourierB(i)*sin(i*xt*pie/12)
            enddo
!          y = ve01 * sin ( 2 * pie * ( xt - 7. ) / 24. )
!          y = y + vpre*exp(-((xt-22)/dpre)**2.)
          return
        endif


!       fejer-scherliess e x b model

	call g(param,funct,xl)

C       **********************************
	y=0.
C       **********************************
	do i=1,index_t
	  do il=1,index_l
	    kk=index_l*(i-1)+il
	    do j=1,nfunc
	       ind=nfunc*(kk-1)+j
	       bspl4=bspl4_time(i,xt)*bspl4_long(il,xl)
               y=y+bspl4*funct(j)*coeff(ind)
	    end do
          end do
         end do

	end
C------------------------------------------------------------------



c       *************************************************
c       *************************************************
        real function bspl4_time(i,x1)
c       *************************************************
	implicit none

	integer i,order/4/,j,k
	real t_t(0:39)
	real x,b(20,20),x1

        data t_t/
     *          0.00,2.75,4.75,5.50,6.25,
     *          7.25,10.00,14.00,17.25,18.00,
     *          18.75,19.75,21.00,24.00,26.75,
     *          28.75,29.50,30.25,31.25,34.00,
     *          38.00,41.25,42.00,42.75,43.75,
     *          45.00,48.00,50.75,52.75,53.50,
     *          54.25,55.25,58.00,62.00,65.25,
     *          66.00,66.75,67.75,69.00,72.00/

	x=x1
        if(i.ge.0) then
          if (x.lt.t_t(i-0)) then
	     x=x+24
	  end if
	end if
	do j=i,i+order-1
	   if(x.ge.t_t(j).and.x.lt.t_t(j+1)) then
	       b(j,1)=1
	   else
	       b(j,1)=0
	   end if
	end do

	do j=2,order
	     do k=i,i+order-j
		b(k,j) = ( x - t_t(k) ) / ( t_t(k+j-1) - t_t(k) ) *
     .                   b(k,j-1)
		b(k,j) = b(k,j) +
     .                   ( t_t(k+j)-x ) / ( t_t(k+j) - t_t(k+1) ) *
     .                    b(k+1,j-1)
	     end do
	end do

	bspl4_time=b(i,order)
	end
C------------------------------------------------------------------



c       *************************************************
c       *************************************************
        real function bspl4_long(i,x1)
c       *************************************************
	implicit none

	integer i,order/4/,j,k
	real t_l(0:24)
	real x,b(20,20),x1

        data t_l/
     *          0,10,100,190,200,250,280,310,
     *          360,370,460,550,560,610,640,670,
     *          720,730,820,910,920,970,1000,1030,1080/

	x=x1
        if(i.ge.0) then
          if (x.lt.t_l(i-0)) then
	     x=x+360
	  end if
	end if
	do j=i,i+order-1
	   if(x.ge.t_l(j).and.x.lt.t_l(j+1)) then
	       b(j,1)=1
	   else
	       b(j,1)=0
	   end if
	end do

	do j=2,order
	     do k=i,i+order-j
		b(k,j)=(x-t_l(k))/(t_l(k+j-1)-t_l(k))*b(k,j-1)
		b(k,j)=b(k,j)+(t_l(k+j)-x)/(t_l(k+j)-t_l(k+1))*
     .                 b(k+1,j-1)
	     end do
	end do

	bspl4_long=b(i,order)
	end
C------------------------------------------------------------------



c       *************************************************
c       *************************************************
        subroutine g(param,funct,x)
c       *************************************************
	implicit none

        integer i
	real param(2),funct(6)
	real x,a,sigma,gauss,flux,cflux

c       *************************************************
	flux=param(2)
        if(param(2).le.75)  flux=75.
        if(param(2).ge.230) flux=230.
	cflux=flux

	a=0.
        if((param(1).ge.120).and.(param(1).le.240)) a=170.
        if((param(1).ge.120).and.(param(1).le.240)) sigma=60
        if((param(1).le.60).or.(param(1).ge.300)) a=170.
        if((param(1).le.60).or.(param(1).ge.300)) sigma=40

	if((flux.le.95).and.(a.ne.0)) then
	 gauss=exp(-0.5*((x-a)**2)/sigma**2)
         cflux=gauss*95.+(1-gauss)*flux
        end if
c       *************************************************

c       *************************************************
        do i=1,6
         funct(i)=0.
        end do
c       *************************************************

c       *************************************************
        if((param(1).ge.135).and.(param(1).le.230)) funct(1)=1
        if((param(1).le.45).or.(param(1).ge.320)) funct(2)=1
        if((param(1).gt.75).and.(param(1).lt.105)) funct(3)=1
        if((param(1).gt.260).and.(param(1).lt.290)) funct(3)=1
c       *************************************************

        if((param(1).ge.45).and.(param(1).le.75)) then  ! W-E
	 funct(2)=1.-(param(1)-45.)/30.
	 funct(3)=1-funct(2)
        end if
        if((param(1).ge.105).and.(param(1).le.135)) then  ! E-S
	 funct(3)=1.-(param(1)-105.)/30.
	 funct(1)=1-funct(3)
        end if
        if((param(1).ge.230).and.(param(1).le.260)) then  ! S-E
	 funct(1)=1.-(param(1)-230.)/30.
	 funct(3)=1-funct(1)
        end if
        if((param(1).ge.290).and.(param(1).le.320)) then  ! E-W
	 funct(3)=1.-(param(1)-290.)/30.
	 funct(2)=1-funct(3)
        end if

c       *************************************************
        funct(4)=(cflux-140)*funct(1)
        funct(5)=(cflux-140)*funct(2)
        funct(6)=(flux-140)*funct(3)
c       *************************************************

	end
C------------------------------------------------------------------


!  NOTE: all variables changed to lower case (JH) 4/27/05

C*********************************************************************C
C*                                                                   *C
C*  chapman.for                                                      *C
C*                                                                   *C
C*  Written by:  David L. Huestis, Molecular Physics Laboratory      *C
C*                                                                   *C
C*  Copyright (c) 2000  SRI International                            *C
C*  All Rights Reserved                                              *C
C*                                                                   *C
C*  This software is provided on an as is basis; without any         *C
C*  warranty; without the implied warranty of merchantability or     *C
C*  fitness for a particular purpose.                                *C
C*                                                                   *C
C*********************************************************************C
C*
C*	To calculate the Chapman Function, Ch(X,chi0), the column
C*	depth of an exponential atmosphere integrated along a line
C*	from a given point to the sun, divided by the column depth for
C*	a vertical sun.
C*
C*  USAGE:
C*
C*	  z = altitude above the surface
C*	  R = radius of the planet
C*	  H = atmospheric scale height
C*
C*	  X = (R+z)/H
C*	  chi0 = solar zenith angle (in degrees)
C*
C*	  implicit real*4(a-h,o-z)
C*	  depth = atm_chapman(X,chi0)	! analytical
C*	  depth = atm_chap_num(X,chi0)	! numerical (chi0 .le. 90)
C*
C*	  implicit real*8(a-h,o-z)
C*	  depth = atm8_chapman(X,chi0)	! analytical
C*	  depth = atm8_chap_num(X,chi0)	! numerical (chi0 .le. 90)
C*
C*  PERFORMANCE:
C*
C*	Compiled and linked using Microsoft FORTRAN 5.1, and executed
C*	in MS-DOS mode under Windows 95 on a 160 MHz PC.
C*
C*    TIMING (in microseconds, typical)
C*
C*	  120	atm_chapman and atm8_chapman for X .lt. 36
C*	   25	atm_chapman and atm8_chapman for X .ge. 36
C*	  500	atm_chap_num
C*	 5000	atm8_chap_num
C*
C*    ACCURACY (maximum relative error, 0.le.chi0.le.90, 1.le.X.le.820)
C*
C*	6.0E-7	atm_chapman and atm8_chapman for X .lt. 60
C*	1.5E-7	atm_chapman and atm8_chapman for X .ge. 60
C*	6.0E-8	atm_chap_num
C*	1.E-15	atm8_chap_num (convergence test)
C*
C*    CODING
C*
C*	No claims are made that the code is optimized for speed,
C*	accuracy, or compactness.  The principal objectives were
C*
C*	  (1) Robustness with respect to argument values
C*	  (2) Rigorous mathematical derivation and error control
C*	  (3) Maximal use of "well known" mathematical functions
C*	  (4) Ease of readability and mapping of theory to coding
C*
C*	The real*8 accuracy could be improved with more accurate
C*	representations of E1(), erfc(), I0(), I1(), K0(), K1().
C*
C*	In the course of development, many representations and
C*	approximations of the Chapman Function were attempted that
C*	failed to be robustly extendable to machine-precision.
C*
C*  INTERNET ACCESS:
C*
C*	Source: http://www-mpl.sri.com/software/chapman/chapman.html
C*	Author: mailto:david.huestis@sri.com
C*	        http://www-mpl.sri.com/bios/Huestis-DL.html
C*
C*  EDIT HISTORY:
C*
C*	01/22/2000 DLH	First complete documentation
C*
C*	01/15/2000 DLH	First complete version of chapman.for
C*
C**********************************************************************
C*
C*  THEORY:
C*
C*    INTRODUCTION
C*
C*	    This computer code models the absorption of solar radiation
C*	by an atmosphere that depends exponentionally on altitude.  In
C*	specific we calculate the effective column depth of a species
C*	of local density, n(z), from a point at a given altitude, z0,
C*	to the sun at a given solar zenith angle, chi0.  Following Rees
C*	[Re89, Section 2.2] we write the column depth for chi0 .le. 90
C*	degrees as
C*
C*   (A)  N(z0,chi0) = int{z=z0,infinity}
C*	     [ n(z)/sqrt( 1 - ( sin(chi0) * (R+z0) / (R+z) ) **2 ) dz ]
C*
C*	where R is the radius of the solid planet (e.g. Earth).  For
C*	chi0 .gt. 90 degrees we write
C*
C*	  N(z0,chi0) = 2*N(zs,90) - N(z0,180-chi0)
C*
C*	where zs = (R+z0)*sin(chi0)-R is the tangent height.
C*
C*	    For an exponential atmosphere, with
C*
C*	  n(z) = n(z0) * exp(-(z-z0)/H)
C*
C*	with a constant scale height, H, the column depth can be
C*	represented by the Chapman function, Ch(X,chi0), named after
C*	the author of the first quantitative mathematical investigation
C*	[Ch31b] trough the relation
C*
C*	  N(z0,chi0) = H * n(z0) * Ch(X,chi0)
C*
C*	where X = (R+z0)/H is a dimensionless measure of the radius
C*	of curvature, with values from about 300 to 1300 on Earth.
C*
C*
C*    APPROACH
C*
C*	    We provide function entry points for very stable and
C*	reasonably efficient evaluation of Ch(X,chi0) with full
C*	single-precision accuracy (.le. 6.0E-7 relative) for a wide
C*	range of parameters.  A 15-digit-accurate double precision
C*	numerical integration routine is also provided.
C*
C*	    Below we will develop (1) a compact asymptotic expansion of
C*	good accuracy for moderately large values of X (.gt. 36) and all
C*	values of chi0, (2) an efficient numerical integral for
C*	all values of X and chi0, and (3) an explicit analytical
C*	representation, valid for all values of X and chi0, based
C*	the differential equation satisfied by Ch(X,chi0).
C*
C*	    All three of these represent new research results as well
C*	as significant computational improvements over the previous
C*	literature, much of which is cited below.
C*
C*
C*    CHANGES OF THE VARIABLE OF INTEGRATION
C*
C*	Substituting y = (R+z)/(R+z0) - 1 we find
C*
C*   (B)  Ch(X,chi0) = X * int{y=0,infinity}
C*	     [ exp(-X*y) / sqrt( 1 - ( sin(chi0) / (1+y) )**2 ) dy ]
C*
C*	The futher substitutions s = (1+y)/sin(chi0), s0 = 1/sin(chi0)
C*	give
C*
C*   (C)  Ch(X,chi0) = X*sin(chi0) * int{s=s0,infinity}
C*	     [ exp(X*(1-sin(chi0)*s)) * s / sqrt(s**2-1) ds ]
C*
C*	From this equation we can establish that
C*
C*	  Ch(X,90) = X*exp(X)*K1(X)
C*
C*	[AS64, Equations 9.6.23 and 9.6.27].  If we now substitute
C*	s = 1/sin(lambda) we obtain
C*
C*   (D)  Ch(X,chi0) = X*sin(chi0) * int{lambda=0,chi0}
C*	    [ exp(X*(1-sin(chi0)*csc(lambda))) * csc(lambda)**2 dlambda]
C*
C*	which is the same as Chapman's original formulation [Ch31b, p486,
C*	eqn (10)].  If we first expand the square root in (B)
C*
C*	  1/sqrt(1-q) = 1 + q/( sqrt(1-q)*(1+sqrt(1-q)) )
C*
C*	with q = ( sin(chi0) / (1+y) )**2 = sin(lambda)**2, we obtain
C*	a new form of (D) without numerical sigularities and simple
C*	convergence to Ch(0,chi0) = Ch(X,0) = 1
C*
C*   (E)  Ch(X,chi0) = 1 + X*sin(chi0) * int{lambda=0,chi0}
C*	    [ exp(X*(1-sin(chi0)*csc(lambda)))
C*		/ (1 + cos(lambda) ) dlambda ]
C*
C*	Alternatively, we may substitute t**2 = y + t0**2,
C*	into Equation (B), with t0**2 = 1-sin(chi0), finding
C*
C*   (F)  Ch(X,chi0) = X * int{s=t0,infinity}
C*	    [ exp(-X*(t**2-t0**2)) * f(t,chi0) dt ]
C*
C*	where
C*
C*	  f(t,chi0) = (t**2 + sin(chi0)) / sqrt(t**2+2*sin(chi0))
C*
C*	  f(t,chi0) = (t**2-t0**2+1)/sqrt(t**2-t0**2+1+sin(chi0))
C*
C*	    Below we will use Equation (F) above to develop a
C*	compact asymptotic expansion of good accuracy for moderately
C*	large values of X (.gt. 36) and all values of chi0, Equation (E)
C*	to develop an efficient numerical integral for Ch(X,chi0) for
C*	all values of X and chi0, and Equation (C) to derive an explicit
C*	analytical representation, valid for all values of X and chi0,
C*	based on the differential equation satisfied by Ch(X,chi0).
C*
C*    atm_chapman(X,chi0) and atm8_chapman(X,chi0)
C*
C*	These routines return real*4 and real*8 values of Ch(X,chi0)
C*	selecting the asymptotic expansion or differential equation
C*	approaches, depending on the value of X.  These routines also
C*	handle the case of chi0 .gt. 90 degrees.
C*
C*    atm_chap_num(X,chi0) and atm8_chap_num(X,chi0)
C*
C*	These routines return real*4 and real*8 values of Ch(X,chi0)
C*	evaluated numerically.  They are both more accurate than the
C*	corresponding atm*_chapman() functions, but take significantly
C*	more CPU time.
C*
C*
C*    ASYMPTOTIC EXPANSION
C*
C*	From Equation (F) we expand, with t0**2 = 1-sin(chi0),
C*
C*	  f(t,chi0) = sum{n=0,3} [ C(n,chi0) * (t**2-t0**2)**n ]
C*
C*	The function atm8_chap_asy(X,chi0) evaluates integrals of the
C*	form
C*
C*	  int{t=t0,infinity} [exp(-X*(t**2-t0**2))*(t**2-t0**2)**n dt]
C*
C*	in terms of incomplete gamma functions, and sums them to
C*	compute Ch(X,chi0).  For large values of X, this results in an
C*	asymptotic expansion in negative powers of X, with coefficients
C*	that are stable for all values of chi0.
C*
C*	In contrast, the asymptotic expansions of Chapman [Ch31b,
C*	p488, Equation (22) and p490, Equation (38)], Hulburt [He39],
C*	and Swider [Sw64, p777, Equation (43)] use negative powers of
C*	X*cos(chi0)**2 or X*sin(chi0), and are accurate only for
C*	small values or large values of chi0, respectively.
C*
C*	Taking only the first term in the present expansion gives the
C*	simple formula
C*
C*	  Ch(X,chi0) = sqrt(pi*X/(1+sin(chi0))) * exp(X*(1-sin(chi0)))
C*		* erfc( sqrt(X*(1-sin(chi0))) )
C*
C*	This is slightly more accurate than the semiempirical
C*	formula of Fitzmaurice [Fi64, Equation (3)], and sightly less
C*	accurate than that of Swider [Sw64, p780, Equation (52),
C*	corrected in SG69].
C*
C*
C*    NUMERICAL INTEGRATION
C*
C*	We are integrating
C*
C*   (E)  Ch(X,chi0) = 1 + X*sin(chi0) * int{lambda=0,chi0}
C*	    [ exp(X*(1-sin(chi0)*csc(lambda)))
C*		/ ( 1 + cos(lambda) ) dlambda ]
C*
C*	The integrand is numerically very smooth, and rapidly varying
C*	only near lambda = 0.  For X .ne. 0 we choose the lower limit
C*	of numerical integration such that the integrand is
C*	exponentially small, 7.0E-13 (3.0E-20 for real*8).  The domain
C*	of integration is divided into 64 equal intervals (6000 for
C*	real*8), and integrated numerically using the 9-point closed
C*	Newton-Cotes formula from Hildebrand [Hi56a, page 75, Equation
C*	(3.5.17)].
C*
C*
C*    INHOMOGENOUS DIFFERENTIAL EQUATION
C*
C*	    The function atm8_chap_deq(X,chi0) calculates Ch(X,chi0),
C*	based on Equation (C) above, using the inhomogeneous
C*	Bessel's equation as described below.  Consider the function
C*
C*	  Z(Q) = int{s=s0,infinity} [ exp(-Q*s) / sqrt(s**2-1) ds ]
C*
C*	Differentiating with respect to Q we find that
C*
C*	  Ch(X,chi0) = - Q * exp(X) * d/dQ [ Z(Q) ]
C*
C*	with Q = X*sin(chi0), s0 = 1/sin(chi0).  Differentiating
C*	inside the integral, we find that
C*
C*	  Z"(Q) + Z'(Q)/Q - Z(Q) = sqrt(s0**2-1) * exp(-Q*s0) / Q
C*
C*	giving us an inhomogeneous modified Bessel's equation of order
C*	zero.  Following Rabenstein [Ra66, pp43-45,149] the solution
C*	of this equation can be written as
C*
C*	  Z(Q) = A*I0(Q) + B*K0(Q) - sqrt(s0**2-1)
C*	         * int{t=Q,infinity} [ exp(-t*s0)
C*		   * ( I0(Q)*K0(t) - I0(t)*K0(Q) ) dt ]
C*
C*	with coefficients A and B to be determined by matching
C*	boundary conditions.
C*
C*	    Differentiating with respect to Q we obtain
C*
C*	  Ch(X,chi0) = X*sin(chi0)*exp(X)*(
C*		- A*I1(X*sin(chi0)) + B*K1(X*sin(chi0))
C*		+ cos(chi0) * int{y=X,infinity} [ exp(-y)
C*		  * ( I1(X*sin(chi0))*K0(y*sin(chi0))
C*		    + K1(X*sin(chi0))*I0(y*sin(chi0)) ) dy ] )
C*
C*	Applying the boundary condition Ch(X,0) = 1 requires that
C*	B = 0.  Similarly, the requirement that Ch(X,chi0) approach
C*	the finite value of sec(chi0) as X approaches infinity [Ch31b,
C*	p486, Equation (12)] implies A = 0.  Thus we have
C*
C*	  Ch(X,chi0) = X*sin(chi0)*cos(chi0)*exp(X)*
C*		int{y=X,infinity} [ exp(-y)
C*		  * ( I1(X*sin(chi0))*K0(y*sin(chi0))
C*		    + K1(X*sin(chi0))*I0(y*sin(chi0)) ) dy ]
C*
C*	The function atm8_chap_deq(X,chi0) evaluates this expression.
C*	Since explicit approximations are available for I1(z) and K1(z),
C*	the remaining challenge is evaluation of the integrals
C*
C*	  int{y=X,infinity} [ exp(-y) I0(y*sin(chi0)) dy ]
C*
C*	and
C*
C*	  int{y=X,infinity} [ exp(-y) K0(y*sin(chi0)) dy ]
C*
C*	which are accomplished by term-by-term integration of ascending
C*	and descending power series expansions of I0(z) and K0(z).
C*
C*  REFERENCES:
C*
C*	AS64	M. Abramowitz and I. A. Stegun, "Handbook of
C*		Mathematical Functions," NBS AMS 55 (USGPO,
C*		Washington, DC, June 1964, 9th printing, November 1970).
C*
C*	Ch31b	S. Chapman, "The Absorption and Dissociative or
C*		Ionizing Effect of Monochromatic Radiation in an
C*		Atmosphere on a Rotating Earth: Part II. Grazing
C*		Incidence," Proc. Phys. Soc. (London), _43_, 483-501
C*		(1931).
C*
C*	Fi64	J. A. Fitzmaurice, "Simplfication of the Chapman
C*		Function for Atmospheric Attenuation," Appl. Opt. _3_,
C*		640 (1964).
C*
C*	Hi56a	F. B. Hildebrand, "Introduction to Numerical
C*		Analysis," (McGraw-Hill, New York, 1956).
C*
C*	Hu39	E. O. Hulburt, "The E Region of the Ionosphere,"
C*		Phys. Rev. _55_, 639-645 (1939).
C*
C*	PFT86	W. H. Press, B. P. Flannery, S. A. Teukolsky, and
C*		W. T. Vetterling, "Numerical Recipes," (Cambridge,
C*		1986).
C*
C*	Ra66	A. L. Rabenstein, "Introduction to Ordinary
C*		Differential Equations," (Academic, NY, 1966).
C*
C*	Re89	M. H. Rees, "Physics and Chemistry of the Upper
C*		Atmosphere," (Cambridge, 1989).
C*
C*	SG69	W. Swider, Jr., and M. E. Gardner, "On the Accuracy
C*		of Chapman Function Approximations," Appl. Opt. _8_,
C*		725 (1969).
C*
C*	Sw64	W. Swider, Jr., "The Determination of the Optical
C*		Depth at Large Solar Zenith Angles," Planet. Space
C*		Sci. _12_, 761-782 (1964).
C
C  ####################################################################
C
C	Chapman function calculated by various methods
C
C	  Ch(X,chi0) = atm_chapman(X,chi0)   : real*4 entry
C	  Ch(X,chi0) = atm8_chapman(X,chi0)  : real*8 entry
C
C	Internal service routines - user should not call, except for
C	testing.
C
C	  Ch(X,chi0) = atm8_chap_asy(X,chi0) : asymptotic expansion
C	  Ch(X,chi0) = atm8_chap_deq(X,chi0) : differential equation
C	  Ch(X,chi0) = atm_chap_num(X,chi0)  : real*4 numerical integral
C	  Ch(X,chi0) = atm8_chap_num(X,chi0) : real*8 numerical integral
C
C  ####################################################################

C  ====================================================================
C
C	These are the entries for the user to call.
C
C	chi0 can range from 0 to 180 in degrees.  For chi0 .gt. 90, the
C	product X*(1-sin(chi0)) must not be too large, otherwise we
C	will get an exponential overflow.
C
C	For chi0 .le. 90 degrees, X can range from 0 to thousands
C	without overflow.
C
C  ====================================================================

	real*4 function atm_chapman( x, chi0 )
	real*8 atm8_chapman
	atm_chapman = atm8_chapman( dble(x), dble(chi0) )
	return
	end

c  ====================================================================

	real*8 function atm8_chapman( x, chi0 )
	implicit real*8(a-h,o-z)
	parameter (rad=57.2957795130823208768d0)

	if( (x .le. 0) .or. (chi0 .le. 0) .or. (chi0 .ge. 180) ) then
	  atm8_chapman = 1
	  return
	end if

	if( chi0 .gt. 90 ) then
	  chi = 180 - chi0
	else
	  chi = chi0
	end if

	if( x .lt. 36 ) then
	  atm8_chapman = atm8_chap_deq(x,chi)
	else
	  atm8_chapman = atm8_chap_asy(x,chi)
	end if

	if( chi0 .gt. 90 ) then
	  atm8_chapman = 2*exp(x*2*sin((90-chi)/(2*rad))**2)
     *		* atm8_chap_xk1(x*sin(chi/rad)) - atm8_chapman
	end if

	return
	end

c  ====================================================================
c
c	this chapman function routine calculates
c
c	  ch(x,chi0) = atm8_chap_asy(x,chi0)
c		     = sum{n=0,3} [c(n) * int{t=t0,infinity}
c			[ exp(-x*(t**2-t0**2) * (t**2-t0**2)**n dy ] ]
c
c	with t0**2 = 1 - sin(chi0)
c
c  ====================================================================

	real*8 function atm8_chap_asy( x, chi0 )
	implicit real*8(a-h,o-z)
	parameter (rad=57.2957795130823208768d0)
	dimension c(0:3), xi(0:3), dn(0:3)
	common/atm8_chap_cm/fn(0:3)

	if( (x .le. 0) .or. (chi0 .le. 0) ) then
	  do i=0,3
	    fn(i) = 1
	  end do
	  go to 900
	end if

	sinchi = sin(chi0/rad)
	s1 = 1 + sinchi
	rx = sqrt(x)
	y0 = rx * sqrt( 2*sin( (90-chi0)/(2*rad) )**2 )

	c(0) = 1/sqrt(s1)
	fact = c(0)/s1
	c(1) = fact * (0.5d0+sinchi)
	fact = fact/s1
	c(2) = - fact * (0.125d0+0.5d0*sinchi)
	fact = fact/s1
	c(3) = fact * (0.0625d0+0.375d0*sinchi)

	call atm8_chap_gd3( y0, dn )
	fact = 2*rx
	do n=0,3
	  xi(n) = fact * dn(n)
	  fact = fact/x
	end do

	fn(0) = c(0) * xi(0)
	do i=1,3
	  fn(i) = fn(i-1) + c(i)*xi(i)
	end do

900	atm8_chap_asy = fn(3)
	return
	end

c  ====================================================================
c
c	this chapman function routine calculates
c
c	  ch(x,chi0) = atm8_chap_deq(x,chi0)
c		     = x * sin(chi0) * cos(chi0) * exp(x*sin(chi0))
c		       * int{y=x,infinity} [ exp(-y)*(
c			 i1(x*sin(chi0))*k0(y*sin(chi0))
c			 + k1(x*sin(chi0))*i0(y*sin(chi0)) ) dy ]
c
c  ====================================================================

	real*8 function atm8_chap_deq( x, chi0 )
	implicit real*8(a-h,o-z)
	parameter (rad=57.2957795130823208768d0)
	common/atm8_chap_cm/xi1,xk1,yi0,yk0

	if( (x .le. 0) .or. (chi0 .le. 0) ) go to 800
	alpha = x * sin(chi0/rad)

c  --------------------------------------------------------------------
c
c	this code fragment calculates
c
c	  yi0 = exp(x*(1-sin(chi0))) * cos(chi0) *
c		int{y=x,infinity} [ exp(-y) * i0(y*sin(chi0)) dy ]
c
c  --------------------------------------------------------------------

	yi0 = atm8_chap_yi0( x, chi0 )

c  --------------------------------------------------------------------
c
c	this code fragment calculates
c
c	  yk0 = exp(x*(1+sin(chi0))) x * sin(chi0) * cos(chi0) *
c		int{y=x,infinity} [ exp(-y) * k0(y*sin(chi0)) dy ]
c
c  --------------------------------------------------------------------

	yk0 = atm8_chap_yk0( x, chi0 )

c  --------------------------------------------------------------------
c
c	this code fragment calculates
c
c	  xi1 = exp(-x*sin(chi0)) * i1(x*sin(chi0))
c
c  --------------------------------------------------------------------

	xi1 = atm8_chap_xi1( alpha )

c  --------------------------------------------------------------------
c
c	this code fragment calculates
c
c	  xk1 = x*sin(chi0) * exp(x*sin(chi0)) * k1(x*sin(chi0))
c
c  --------------------------------------------------------------------

	xk1 = atm8_chap_xk1( alpha )

c  --------------------------------------------------------------------
c
c	combine the terms
c
c  --------------------------------------------------------------------

	atm8_chap_deq = xi1*yk0 + xk1*yi0
	go to 900

800	atm8_chap_deq = 1
900	return
	end

c  ====================================================================
c
c	this chapman function routine calculates
c
c	  ch(x,chi0) = atm_chap_num(x,chi0) = numerical integral
c
c  ====================================================================

	real*4 function atm_chap_num(x,chi0)
	implicit real*8(a-h,o-z)
	real*4 x, chi0
	parameter (rad=57.2957795130823208768d0)
	parameter (n=65,nfact=8)
	dimension factor(0:nfact)
	data factor/14175.0d0, 23552.0d0, -3712.0d0, 41984.0d0,
     *	  -18160.0d0, 41984.0d0, -3712.0d0, 23552.0d0, 7912.0d0/

	if( (chi0 .le. 0) .or. (chi0 .gt. 90) .or. (x .le. 0) ) then
	  atm_chap_num = 1
	  return
	end if

	x8 = x
	chi0rad = chi0/rad
	sinchi = sin(chi0rad)

	alpha0 = asin( (x8/(x8+28)) * sinchi )
	delta = (chi0rad - alpha0)/(n-1)

	sum = 0

	do i=1,n
	  alpha = -(i-1)*delta + chi0rad

	  if( (i .eq. 1) .or. (x .le. 0) ) then
	    f = 1/(1+cos(alpha))
	  else if( alpha .le. 0 ) then
	    f = 0
	  else
	    f = exp(-x8*(sinchi/sin(alpha)-1) ) /(1+cos(alpha))
	  end if

	  if( (i.eq.1) .or. (i.eq.n) ) then
	    fact = factor(nfact)/2
	  else
	    fact = factor( mod(i-2,nfact)+1 )
	  end if

	  sum = sum + fact*f
	end do

	atm_chap_num = 1 + x8*sinchi*sum*delta/factor(0)
	return
	end

c  ====================================================================
c
c	this chapman function routine calculates
c
c	  ch(x,chi0) = atm8_chap_num(x,chi0) = numerical integral
c
c  ====================================================================

	real*8 function atm8_chap_num(x,chi0)
	implicit real*8(a-h,o-z)
	parameter (rad=57.2957795130823208768d0)
	parameter (n=601,nfact=8)
	dimension factor(0:nfact)
	data factor/14175.0d0, 23552.0d0, -3712.0d0, 41984.0d0,
     *	  -18160.0d0, 41984.0d0, -3712.0d0, 23552.0d0, 7912.0d0/

	if( (chi0 .le. 0) .or. (chi0 .gt. 90) .or. (x .le. 0) ) then
	  atm8_chap_num = 1
	  return
	end if

	chi0rad = chi0/rad
	sinchi = sin(chi0rad)

	alpha0 = asin( (x/(x+45)) * sinchi )
	delta = (chi0rad - alpha0)/(n-1)

	sum = 0

	do i=1,n
	  alpha = -(i-1)*delta + chi0rad

	  if( (i .eq. 1) .or. (x .le. 0) ) then
	    f = 1/(1+cos(alpha))
	  else if( alpha .le. 0 ) then
	    f = 0
	  else
	    f = exp(-x*(sinchi/sin(alpha)-1) ) /(1+cos(alpha))
	  end if

	  if( (i.eq.1) .or. (i.eq.n) ) then
	    fact = factor(nfact)/2
	  else
	    fact = factor( mod(i-2,nfact)+1 )
	  end if

	  sum = sum + fact*f
	end do

	atm8_chap_num = 1 + x*sinchi*sum*delta/factor(0)
	return
	end

c  ####################################################################
c
c	the following "bessel integral" routines return various
c	combinations of integrals of bessel functions, powers,
c	and exponentials, involving trigonometric functions of chi0.
c
c	for small values of z = x*sin(chi0) we expand
c
c	  i0(z) = sum{n=0,6} [ ai0(n) * z**(2*n) ]
c	  k0(z) = -log(z)*i0(z) + sum{n=0,6} [ ak0(n) * z**(2*n) ]
c
c	for large values of z we expand in reciprocal powers
c
c	  i0(z) = exp(z) * sum{n=0,8} [ bi0(n) * z**(-n-0.5) ]
c	  k0(z) = exp(-z) * sum{n=0,6} [ bk0(n) * z**(-n-0.5) ]
c
c	the expansion coefficients are calculated from those given
c	by abramowitz and stegun [as64, pp378-9, section 9.8] and
c	press et al. [pft86, pp177-8, bessi0.for, bessk0.for].
c
c	for small values of x*sin(chi0) we break the integral
c	into two parts (with f(z) = i0(z) or k0(z)):
c
c	  int{y=x,infinity} [ exp(-y) * f(y*sin(chi0)) dy ]
c
c	    = int{y=x,x1} [ exp(-y) * f(y*sin(chi0)) dy ]
c	      + int{y=x1,infinity} [ exp(-y) * f(y*sin(chi0)) dy ]
c
c	where x1 = 3.75/sin(chi0) for i0 and 2/sin(chi0) for k0.
c
c	in the range y=x,x1 we integrate the term-by-term using
c
c	  int{z=a,b} [ exp(-z) * z**(2*n) dz ]
c	    = gamma(2*n+1,a) - gamma(2*n+1,b)
c
c	and a similar but more complicated formula for
c
c	  int{z=a,b} [ log(z) * exp(-z) * z**(2*n) dz ]
c
c	in the range y=x1,infinity we use
c
c	  int{z=b,infinity} [ exp(-z) * z**(-n-0.5) dz]
c	    = gamma(-n+0.5,b)
c
c  ####################################################################

c  ====================================================================
c
c	this bessel integral routine calculates
c
c	  yi0 = exp(x*(1-sin(chi0))) * cos(chi0) *
c		int{y=x,infinity} [ exp(-y) * i0(y*sin(chi0)) dy ]
c
c  ====================================================================

	real*8 function atm8_chap_yi0( x, chi0 )
	implicit real*8(a-h,o-z)
	parameter (rad=57.2957795130823208768d0)
	dimension qbeta(0:8), gg(0:6)
	dimension ai0(0:6), bi0(0:8)

        data ai0/ 1.0000000d+00, 2.4999985d-01, 1.5625190d-02,
     *      4.3393973d-04, 6.8012343d-06, 6.5601736d-08,
     *      5.9239791d-10/
        data bi0/ 3.9894228d-01, 4.9822200d-02, 3.1685484d-02,
     *     -8.3090918d-02, 1.8119815d+00,-1.5259477d+01,
     *      7.3292025d+01,-1.7182223d+02, 1.5344533d+02/

	theta = (90-chi0)/(2*rad)
	sint = sin(theta)
	cost = cos(theta)
	sinchi = sin(chi0/rad)
	coschi = cos(chi0/rad)
	sc1m = 2*sint**2	! = (1-sinchi)

	alpha = x * sinchi

	if( alpha .le. 0 ) then
	  atm8_chap_yi0 = 1
	else if( alpha .lt. 3.75d0 ) then
	  x1 = 3.75d0/sinchi
	  call atm8_chap_gg06( x, x1, gg )
	  if( x .le. 1 ) then
	    rho = 1
	  else
	    rho = 1/x
	  end if
	  f = (sinchi/rho)**2
	  sum = ai0(6)*gg(6)
	  do i=5,0,-1
	    sum = sum*f + ai0(i)*gg(i)
c	    write(*,1900)i,sum,gg(i)
c1900	format(i5,1p5d14.6)
	  end do
	  call atm8_chap_gq85( x1*sc1m, qbeta )
	  sum2 = bi0(8) * qbeta(8)
	  do n=7,0,-1
	    sum2 = sum2/3.75d0 + bi0(n)*qbeta(n)
	  end do
	  atm8_chap_yi0 = exp(-alpha)*coschi*sum
     *		+ exp((x-x1)*sc1m)*sum2*cost*sqrt(2/sinchi)
	else
	  call atm8_chap_gq85( x*sc1m, qbeta )
	  sum = bi0(8) * qbeta(8)
	  do n=7,0,-1
	    sum = sum/alpha + bi0(n)*qbeta(n)
	  end do
	  atm8_chap_yi0 = sum * cost * sqrt( 2 / sinchi )
	end if
	return
	end

c  ====================================================================
c
c	this bessel integral routine calculates
c
c	  yk0 = exp(x*(1+sin(chi0))) x * sin(chi0) * cos(chi0) *
c		int{y=x,infinity} [ exp(-y) * k0(y*sin(chi0)) dy ]
c
c  ====================================================================

	real*8 function atm8_chap_yk0( x, chi0 )
	implicit real*8(a-h,o-z)
	parameter (rad=57.2957795130823208768d0)
	dimension ai0(0:6), ak0(0:6), bk0(0:6)
	dimension gf(0:6), gg(0:6), qgamma(0:8)

        data ai0/ 1.0000000d+00, 2.4999985d-01, 1.5625190d-02,
     *      4.3393973d-04, 6.8012343d-06, 6.5601736d-08,
     *      5.9239791d-10/
        data ak0/ 1.1593152d-01, 2.7898274d-01, 2.5249154d-02,
     *      8.4587629d-04, 1.4975897d-05, 1.5045213d-07,
     *      2.2172596d-09/
        data bk0/ 1.2533141d+00,-1.5664716d-01, 8.7582720d-02,
     *     -8.4995680d-02, 9.4059520d-02,-8.0492800d-02,
     *      3.4053120d-02/

	theta = (90-chi0)/(2*rad)
	sint = sin(theta)
	cost = cos(theta)
	sinchi = sin(chi0/rad)
	sc1 = 1+sinchi
	coschi = sin(2*theta)

	alpha = x * sinchi
	gamma = x * sc1

	if( alpha .le. 0 ) then
	  atm8_chap_yk0 = 0
	else if( alpha .lt. 2 ) then
	  x1 = 2/sinchi
	  call atm8_chap_gfg06( x, x1, gf, gg )
	  if( x .le. 1 ) then
	    rho = 1
	  else
	    rho = 1/x
	  end if
	  sl = log(sinchi)
	  f = (sinchi/rho)**2
	  sum = -ai0(6)*gf(6) + (-sl*ai0(6)+ak0(6))*gg(6)
	  do i=5,0,-1
	    sum = sum*f - ai0(i)*gf(i) + (-sl*ai0(i)+ak0(i))*gg(i)
c	    write(*,1900)i,sum,gf(i),gg(i)
c1900	format(i5,1p5d14.6)
	  end do
	  call atm8_chap_gq85( x1*sc1, qgamma )
	  sum2 = bk0(6)*qgamma(6)
	  do i=5,0,-1
	    sum2 = sum2*0.5d0 + bk0(i)*qgamma(i)
c	    write(*,1900)i,sum2,bk0(i),qgamma(i)
	  end do
	  sum = sum + exp(x-x1-2)*sum2/sqrt(sinchi*sc1)
	  atm8_chap_yk0 = sum * exp(alpha) * alpha * coschi
	else
	  call atm8_chap_gq85( gamma, qgamma )
	  sum = bk0(6) * qgamma(6)
	  do i=5,0,-1
	    sum = sum/alpha + bk0(i)*qgamma(i)
	  end do
	  atm8_chap_yk0 = sum * sint * sqrt( 2 * sinchi ) * x
	end if

	return
	end

c  ####################################################################
c
c	the following "pure math" routines return various combinations
c	of bessel functions, powers, and exponentials.
c
c  ####################################################################

c  ====================================================================
c
c	this bessel function math routine returns
c
c	  xi1 = exp(-|z|) * i1(z)
c
c	following press et al [pft86, page 178, bessi1.for] and
c	abrahamson and stegun [as64, page 378, 9.8.3, 9.8.4].
c
c  ====================================================================

	real*8 function atm8_chap_xi1( z )
	implicit real*8(a-h,o-z)
        dimension ai1(0:6), bi1(0:8)

        data ai1/ 5.00000000d-01, 6.2499978d-02, 2.6041897d-03,
     *      5.4244512d-05, 6.7986797d-07, 5.4830314d-09,
     *      4.1909957d-11/
        data bi1/ 3.98942280d-01,-1.4955090d-01,-5.0908781d-02,
     *      8.6379434d-02,-2.0399403d+00, 1.6929962d+01,
     *     -8.0516146d+01, 1.8642422d+02,-1.6427082d+02/

	if( z .lt. 0 ) then
	  az = -z
	else if( z .eq. 0 ) then
	  atm8_chap_xi1 = 0
	  return
	else
	  az = z
	end if
	if( az .lt. 3.75d0 ) then
	  z2 = z*z
	  sum = ai1(6)
	  do i=5,0,-1
	    sum = sum*z2 + ai1(i)
	  end do
	  atm8_chap_xi1 = z*exp(-az) * sum
	else
	  sum = bi1(8)
	  do i=7,0,-1
	    sum = sum/az + bi1(i)
	  end do
	  atm8_chap_xi1 = sum*sqrt(az)/z
	end if
	return
	end

c  ====================================================================
c
c	this bessel function math routine returns
c
c	  xk1 = z * exp(+z) * k1(z)
c
c	following press et al [pft86, page 179, bessk1.for] and
c	abrahamson and stegun [as64, page 379, 9.8.7, 9.8.8].
c
c  ====================================================================

	real*8 function atm8_chap_xk1( z )
	implicit real*8(a-h,o-z)
        dimension ak1(0:6), bk1(0:6)

        data ak1/ 1.00000000d+00, 3.8607860d-02,-4.2049112d-02,
     *     -2.8370152d-03,-7.4976641d-05,-1.0781641d-06,
     *     -1.1440430d-08/
        data bk1/ 1.25331414d+00, 4.6997238d-01,-1.4622480d-01,
     *      1.2034144d-01,-1.2485648d-01, 1.0419648d-01,
     *     -4.3676800d-02/

	if( z .le. 0 ) then
	  atm8_chap_xk1 = 1
	else if( z .lt. 2 ) then
	  xz = exp(z)
	  z2 = z*z
	  sum = ak1(6)
	  do i=5,0,-1
	    sum = sum*z2 + ak1(i)
	  end do
	  atm8_chap_xk1 = xz * ( sum
     *		+ z*log(z/2)*atm8_chap_xi1(z)*xz )
	else
	  sum = bk1(6)
	  do i=5,0,-1
	    sum = sum/z + bk1(i)
	  end do
	  atm8_chap_xk1 = sum*sqrt(z)
	end if

	return
	end

c  ####################################################################
c
c	the following "pure math" routines return various combinations
c	of the error function, powers, and exponentials.
c
c  ####################################################################

c  ====================================================================
c
c	this error function math routine returns
c
c	  xerfc(x) = exp(x**2)*erfc(x)
c
c	following press et al. [pft86, p164, erfcc.for]
c
c  ====================================================================

	real*8 function atm8_chap_xerfc(x)
	implicit real*8(a-h,o-z)
        t=1.0d0/(1.0d0+0.5d0*x)
	atm8_chap_xerfc =
     *	  t*exp( -1.26551223d0 +t*(1.00002368d0 +t*( .37409196d0
     *       +t*(  .09678418d0 +t*(-.18628806d0 +t*( .27886807d0
     *	     +t*(-1.13520398d0 +t*(1.48851587d0 +t*(-.82215223d0
     *	     +t*   .17087277d0) ))))))))
        return
        end

c  ####################################################################
c
c	the following "pure math" routines return various combinations
c	of exponential integrals, powers, and exponentials.
c
c  ####################################################################

c  ====================================================================
c
c	this exponential math routine evaluates
c
c	  zxe1(x) = x*exp(x) int{y=1,infinity} [ exp(-x*y)/y dy ]
c
c	following abramowitz and stegun [as64, p229;231, equations
c	5.1.11 and 5.1.56]
c
c  ====================================================================

	real*8 function atm8_chap_zxe1(x)
	implicit real*8(a-h,o-z)
	parameter (gamma = 0.5772156649015328606d0)
	dimension ae1(0:4), be1(0:4), cein(1:10)

	data ae1/1.0d0, 8.5733287401d0, 18.0590169730d0,
     *	    8.6347608925d0, 0.2677737343d0 /
	data be1/1.0d0, 9.5733223454d0, 25.6329561486d0,
     *	    21.0996530827d0, 3.9584969228d0/
        data cein/ 1.00000000d+00,-2.50000000d-01, 5.55555556d-02,
     *    -1.0416666667d-02, 1.6666666667d-03,-2.3148148148d-04,
     *     2.8344671202d-05,-3.1001984127d-06, 3.0619243582d-07,
     *    -2.7557319224d-08/

	if( x .le. 0 ) then
	  atm8_chap_zxe1 = 0
	else if( x .le. 1 ) then
	  sum = cein(10)
	  do i=9,1,-1
	    sum = sum*x + cein(i)
	  end do
	  atm8_chap_zxe1 = x*exp(x)*( x * sum - log(x) - gamma )
	else
	  top = ae1(4)
	  bot = be1(4)
	  do i=3,0,-1
	    top = top/x + ae1(i)
	    bot = bot/x + be1(i)
	  end do
	  atm8_chap_zxe1 = top/bot
	end if
	return
	end

c  ####################################################################
c
c	the following "pure math" routines return various combinations
c	of incomplete gamma functions, powers, and exponentials.
c
c  ####################################################################

c  ====================================================================
c
c	this gamma function math routine calculates
c
c	dn(n) = int{t=z,infinity}
c		[ exp( -(t**2-z**2) ) * (t**2-z**2)**n dt ]
c
c  ====================================================================

	subroutine atm8_chap_gd3( z, dn )
	implicit real*8(a-h,o-z)
	parameter (rpi=1.7724538509055160273d0)
	dimension dn(0:3), xg(0:3)

	if( z .le. 0 ) then
	  dn(0) = rpi/2
	  do i=1,3
	    dn(i) = (i-0.5d0)*dn(i-1)
	  end do
	  return
	end if

	z2 = z*z
	if( z .ge. 7 ) r = 1/z2

	if( z .lt. 14 ) then
	  z4 = z2*z2
	  xg(0) = rpi * atm8_chap_xerfc(z)
	  xg(1) = 0.5d0*xg(0) + z
	  xg(2) = 1.5d0*xg(1) + z*z2
	  dn(0) = 0.5d0*xg(0)
	  dn(1) = 0.5d0*(xg(1)-z2*xg(0))
	  dn(2) = 0.5d0*(xg(2)-2*z2*xg(1)+z4*xg(0))
	else
	  dn(0) = ( 1 + r*(-0.5d0 +r*(0.75d0 +r*(-1.875d0
     *		+r*6.5625d0) ) ) )/(2*z)
	  dn(1) = ( 1 + r*(-1.0d0 +r*(2.25d0 +r*(-7.5d0
     *		+r*32.8125d0) ) ) )/(2*z)
	  dn(2) = ( 2 + r*(-3.0d0 +r*(9.00d0 +r*(-37.5d0
     *		+r*196.875d0) ) ) )/(2*z)
	end if

	if( z .lt. 7 ) then
	  z6 = z4*z2
	  xg(3) = 2.5d0*xg(2) + z*z4
	  dn(3) = 0.5d0*(xg(3)-3*z2*xg(2)+3*z4*xg(1)-z6*xg(0))
	else
	  dn(3) = ( 6 + r*(-12.0d0 +r*(45.0d0 +r*(-225.0d0
     *		+r*1378.125d0) ) ) )/(2*z)
	end if

	return
	end

c  ====================================================================
c
c	this gamma function math routine calculates
c
c	  gf06(n) = g(n,x) * int{y=x,z} [log(y) * exp(-y) * y**(2*n) dy]
c
c	and
c
c	  gg06(n) = g(n,x) * int{y=x,z} [ exp(-y) * y**(2*n) dy ]
c	          = g(n,x) * ( gamma(2*n+1,x) - gamma(2*n+1,z) )
c
c	for n=0,6, with g(n,x) = exp(x) * max(1,x)**(-2*n)
c
c  ====================================================================

	subroutine atm8_chap_gfg06( x, z, gf06, gg06 )
	implicit real*8 (a-h,o-z)
	parameter (gamma = 0.5772156649015328606d0)
	dimension gf06(0:6), gg06(0:6)
	dimension gh13x(13), gh13z(13), rgn(13), delta(13)
	call atm8_chap_gh13( x, x, gh13x )
	call atm8_chap_gh13( x, z, gh13z )
	if( x .le. 1 ) then
	  rho = 1
	else
	  rho = 1/x
	end if

	delta(1) = 0
	delta(2) = ( gh13x(1) - gh13z(1) ) * rho
	rgn(1) = 1
	rgn(2) = rho
	do n=2,12
	  delta(n+1) = rho*( n*delta(n) + gh13x(n) - gh13z(n) )
	  rgn(n+1) = (n*rho)*rgn(n)

	end do

	if( x .gt. 0 ) then
	  xe1_x = atm8_chap_zxe1(x)/x
	  xlog = log(x)
	end if
	if( z .gt. 0 ) then
	  xe1_z = exp(x-z)*atm8_chap_zxe1(z)/z
	  zlog = log(z)
	end if

	do k=0,6
	  n = 2*k+1
	  if( x .le. 0 ) then
	    gf06(k) = -gamma*rgn(n) + delta(n)
	  else
	    gf06(k) = xlog*gh13x(n) + rgn(n)*xe1_x + delta(n)
	  end if
	  if( z .le. 0 ) then
	    gf06(k) = gf06(k) + gamma*rgn(n)
	  else
	    gf06(k) = gf06(k) - (zlog*gh13z(n) + rgn(n)*xe1_z)
	  end if
	  gg06(k) = gh13x(n) - gh13z(n)
	end do

	return
	end

c  ====================================================================
c
c	this gamma function math routine calculates
c
c	  gg06(n) = g(n,x) * int{y=x,z} [ exp(-y) * y**(2*n) dy ]
c	          = g(n,x) * ( gamma(2*n+1,x) - gamma(2*n+1,z) )
c
c	for n=0,6, with g(n,x) = exp(x) * max(1,x)**(-2*n)
c
c  ====================================================================

	subroutine atm8_chap_gg06( x, z, gg06 )
	implicit real*8 (a-h,o-z)
	dimension gg06(0:6), gh13x(13), gh13z(13)
	call atm8_chap_gh13( x, x, gh13x )
	call atm8_chap_gh13( x, z, gh13z )
	do n=0,6
	  gg06(n) = gh13x(2*n+1) - gh13z(2*n+1)
	end do
	return
	end

c  ====================================================================
c
c	this gamma function math routine calculates
c
c	  gh13(n) = f(n,x) * int{y=z,infinity} [exp(-y) * y**(n-1) dy]
c	          = f(n,x) * gamma(n,z)
c
c	for n=1,13, with f(n,x) = exp(x) * max(1,x)**(-n+1)
c
c  ====================================================================

	subroutine atm8_chap_gh13( x, z, gh13 )
	implicit real*8 (a-h,o-z)
	dimension gh13(13), tab(12)

	if( z .le. 0 ) then
	  gh13(1) = 1
	  do n=1,12
	    gh13(n+1) = n*gh13(n)
	  end do
	  return
	end if

	if( x .le. 1 ) then
	  rho = 1
	else
	  rho = 1/x
	end if
	rhoz = rho * z
	exz = exp(x-z)
	tab(12) = exp( (x-z) + 12*log(rhoz) )
	do n=11,1,-1
	  tab(n) = tab(n+1)/rhoz
	end do
	gh13(1) = exz
	do n=1,12
	  gh13(n+1) = rho*n*gh13(n) + tab(n)
	end do
	return
	end

c  ====================================================================
c
c	this gamma function math subroutine calculates
c
c	  qn(x) = x**n * exp(x) * gamma(-n+0.5,x), n=0,8
c	    = x**n * exp(x) * int{y=x,infinity} [exp(-y)*y**(-n-0.5)dy]
c
c	for x .lt. 2 we first calculate
c
c	  q0(x) = sqrt(pi)*exp(x)*erfc(sqrt(x)) = exp(x)*gamma(0.5,x)
c
c	and use upward recursion.  else, we first calculate
c
c	  q8(x) = x**8 * exp(x) * gamma(-7.5,x)
c
c	following press et al. [pft86, pp162-63, gcf.for] and then
c	recur downward.  also see abramowitz and stegun [as64, 6.5].
c
c  ====================================================================

	subroutine atm8_chap_gq85( x, qn )
	implicit real*8(a-h,o-z)
	parameter (rpi=1.7724538509055160273d0)
	parameter (itmax=100,eps=3.0d-9)
	dimension qn(0:8)

	if( x .le. 0 ) then
	  qn(0) = rpi
	  do i=1,8
	    qn(i) = 0
	  end do
	  return
	end if

	rx = sqrt(x)

	if( x .lt. 2 ) then
	  qn(0) = rpi * atm8_chap_xerfc( rx )
	  do n=1,8
	    qn(n) = ( -rx*qn(n-1) + 1 ) * rx / ( n - 0.5d0 )
	  end do
	else
          gold=0.0d0
	  a0=1.0d0
	  a1=x
	  b0=0.0d0
	  b1=1.0d0
	  fac=1.0d0
	  do 11 n=1,itmax
	    an= (n)
	    ana=an + 7.5d0
	    a0=(a1+a0*ana)*fac
	    b0=(b1+b0*ana)*fac
	    anf=an*fac
	    a1=x*a0+anf*a1
	    b1=x*b0+anf*b1
	    fac=1./a1
            g=b1*fac
	    test = g*eps
	    del = g - gold
	    if( test .lt. 0 ) test = - test
	    if( (del .ge. -test) .and. (del .le. test) ) go to 12
	    gold=g
11        continue
12	  qn(8) = g * rx
	  do n=8,1,-1
	    qn(n-1) = ( (-n+0.5d0)*qn(n)/rx + 1 ) / rx
	  end do
	end if

	return
	end
! *********************
!
!     smoothz
!
! *********************

      subroutine smoothz(finout,ncomp)

      include 'param-1.00.inc'

      dimension finout(nz), tempz(nz)

c
c This is the binomial filter (in x space) as described in
c Birdsall appendix C.
c We have the choice of a compensating filter or not.
c if ncomp=0, no compensation, else compensation
c

c do smoothz in the z direction

       do i = 1,nz
          ip1 = i +1
          if(i .eq. nz) ip1 = 1
          im1 = i -1
          if(i .eq. 1) im1 = nz
          tempz(i) = .25*(finout(im1) +2.*finout(i)
     &                   +finout(ip1))
       enddo
       do i = 1,nz
          finout(i) = tempz(i)
       enddo

       if ( ncomp .ne. 0 ) then

c put in compensator
c the alogrithm below is equivalent to
c fftmp(i)=(1./16.)*(-ff0(i-2)+4.*ff0(i-1)+10.*ff0(i)+4.*ff0(i+1)-ff0(i+2))

c do compensation in the z direction

       const = sqrt(1.4571072)
       do i = 1,nz
          ip1 = i +1
          if(i .eq. nz) ip1 = 1
          finout(i) = const*(finout(i) -.171573*finout(ip1))
       enddo
       do i = nz,1,-1
          im1 = i -1
          if(i .eq. 1) im1 = nz
          finout(i) = const*(finout(i) -.171573*finout(im1))
       enddo

      endif

      return
      end
