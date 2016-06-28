*******************************************
*******************************************

!            PARAM-1.00.INC

*******************************************
*******************************************

!      number of altitudes

       parameter ( nf   = 98,
     .             nfp1 = nf + 1,
     .             nfm1 = nf - 1,
     .             nfm2 = nf - 2  )

!      number of grid cells
     
       parameter ( nz = 101,
     .             nzp1 = nz + 1,
     .             nzm1 = nz - 1  )

!      ion densities

       integer nion,pthp,pthep,ptnp,ptop,ptn2p,ptnop,pto2p

       parameter ( nion  = 7 )   ! number of ions
       parameter ( pthp  = 1 )   ! h+
       parameter ( pthep = 5 )   ! he+
       parameter ( ptnp  = 7 )   ! n+
       parameter ( ptop  = 2 )   ! o+
       parameter ( ptn2p = 6 )   ! n2+
       parameter ( ptnop = 3 )   ! no+
       parameter ( pto2p = 4 )   ! o2+

!      neutrals 

       integer nneut,pth,pthe,ptn,pto,ptn2,ptno,pto2

       parameter ( nneut = 7 )   ! number of neutrals
       parameter ( pth   = 1 )   ! h
       parameter ( pthe  = 5 )   ! he
       parameter ( ptn   = 7 )   ! n
       parameter ( pto   = 2 )   ! o
       parameter ( ptn2  = 6 )   ! n2
       parameter ( ptno  = 3 )   ! no
       parameter ( pto2  = 4 )   ! o2

!      number of chemical reactions

       parameter ( nchem = 21 )

!      various constants

!      ftnchek is giving some meaningless errors about precision, 
!      but i am going to lower the precision of some of these
!      variables to keep down the error messages


       parameter ( pie    = 3.1415927   )
       parameter ( po180  = 1.745329e-02 )
       parameter ( rtod   = 57.295780    )
       parameter ( tm18   = 1.e-18       ) 

       parameter ( spday  = 86400., sphr = 3600.   )

       parameter ( gzero  = 980.665, re = 6370., bmag = 0.25 )  

       parameter ( bolt   = 1.38044e-16 )
       parameter ( amu    = 1.67252e-24 )
       parameter ( evtok  = 1.1604e+04  ) 

       parameter ( linesuv = 37, linesnt = 4 )

       parameter ( dayve = 80., sidyr = 365.4, solinc = 23.5 )

!      these are for the error function

       parameter ( pas =   .3275911, z1 =  .2548295,
     .             z2  = - .284496 , z3 = 1.421413,
     .             z4  = -1.453152 , z5 = 1.0614054  )


