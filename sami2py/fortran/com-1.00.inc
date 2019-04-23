*******************************************
*******************************************

!            COM-1.00.INC

*******************************************
*******************************************

!     namelist data

      logical fejer,fmtout,outn
      real snn(nneut)

      common / nmlst / fmtout,maxstep,hrmax,dt0,dthr,hrpr,
     .                 grad_in,glat_in,glon_in,
     .                 fejer,
     .                 rmin,rmax,
     .                 altmin,
     .                 fbar,f10p7,ap,
     .                 year,day,mmass,
     .                 nion1,nion2,hrinit,tvn0,tvexb0,ve01,
     .                 gams,gamp,snn,stn,denmin,alt_crit,cqe,
     .                 Tinf_scl,euv_scl,hwm_mod,  ! Added input variables - JK
     .		       outn                       ! Added variable -JMS


!     s grid data

      real alts(nz,nf),grs(nz,nf),glats(nz,nf),glons(nz,nf)
      real bms(nz,nf),gs(nz),ps(nz,nf),blats(nz,nf)
      real coschicrit(nz,nf)
      real ds(nz,nf),d2s(nz,nf),d22s(nz,nf)
      real dels(nz,nf)
      real grad_inp (nf)
      real xnorms(nzp1,nf),ynorms(nzp1,nf),znorms(nzp1,nf)
      real xnormp(nz,nfp1),ynormp(nz,nfp1),znormp(nz,nfp1)
      real arg(nz,nf),athg(nz,nf),aphig(nz,nf)

      common / sgrid / alts,grs,glats,glons,
     .                 bms,gs,ps,blats,
     .                 ds,d2s,d22s,
     .                 dels,
     .                 grad_inp,
     .                 xnorms,ynorms,znorms,
     .                 xnormp,ynormp,znormp,coschicrit,
     .                 arg,athg,aphig

!       alts     altitude  (in km) on s mesh
!       grs      radial geographic distance to field line on s mesh
!       glats    geographic latitude on s mesh
!       glons    geographic longitude on s mesh
!       bms      normalized magnetic field on s mesh (b/b0)
!       ds,d2,
!       d22s     differential `distances' used in diff eqs
!       dels     actual arc length of grid in s direction

!     p grid data

      real delsp(nz,nfp1),vol(nz,nf)
      real areap(nz,nfp1),areas(nzp1,nf)
      real vnx(nzp1,nfp1),vny(nzp1,nfp1),vnz(nzp1,nfp1)
      real xdels(nz,nfp1,2),xdelp(nzp1,nf,2)
      real vexbs(nzp1,nf),vexbp(nz,nfp1),vexb(nzp1,nfp1)

      common / pgrid / delsp,vol,
     .                 areap,areas,
     .                 vnx,vny,vnz,
     .                 xdels,xdelp,
     .                 vexbs,vexbp,vexb

!     delsp      actual arc length of grid in s direction on p grid
!     vol        volume (i.e., area) of cell

!     alternate ExB data - JK

      real fourierA(10), fourierB(10)

      common  / altexb / fourierA, fourierB

!     chemical reaction data

      integer ichem(nchem,3)
      real ireact(nion,nneut,nchem)

      common / chem / ichem,ireact

!     variables

      real deni(nz,nf,nion),denn(nz,nf,nneut),ne(nz,nf)
      real vsi(nz,nf,nion),vsid(nz,nf,nion)
      real sumvsi(nz,nf,nion),vsic(nz,nf,nion)
      real te(nz,nf),ti(nz,nf,nion),tn(nz,nf)
      real u(nz,nf),v(nz,nf),vpi(nz,nf)

      common / var / deni,denn,ne,te,ti,tn,u,v,vsi,vsid
     .              ,sumvsi,vpi,vsic,ftc,dt

!     velocity in radial (vor) and theta (vot) directions

      real vot(nz,nf,nion),vor(nz,nf,nion)

      common / dvel / vot,vor

!     atomic masses

      real ami(nion),amn(nneut),alpha0(nneut),aap(7)

      common / atm / ami,amn,alpha0,aap

!     zenith datt

      real cx(nz,nf)

      common / geom / cx

!     photodeposition rates
!     used 3 (absorption) and 7 (nion) explicitly
!     used 4 (number of angles in nighttime deposition)

      real sigabsdt(linesuv,3),flux(linesuv),sigidt(linesuv,7)
      real sigint(linesnt,7),fluxnt(nz,nf,91,linesnt)
      real thetant(linesnt,4),zaltnt(linesnt,2)

      common / photdep / sigabsdt,flux,sigidt,
     .                   sigint,fluxnt,thetant,zaltnt

!     diagnostic variables

      real t1(nz,nf,nion),t2(nz,nf,nion),t3(nz,nf,nion)
      real u1(nz,nf),u2(nz,nf),u3(nz,nf),u4(nz,nf),u5(nz,nf)

      common / diag / t1,t2,t3,u1,u2,u3,u4,u5

      real x0,y0,z0,plat,plon,bb0
      common / igrfstuff / x0,y0,z0,plat,plon,bb0
