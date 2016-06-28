
README-1.00 for SAMI2
J.D. Huba, 3/11


****************************************************************

02/10/2011:

Changes from sami2-0.98 to sami2-1.00:

1. The vnormal subroutine has been rewritten
   (i.e., corrected) in grid-1.00.f.

2. The factor ch1 in the subroutine
   photprod has been changed from 1.e38
   to 1.e22 (line 933). The larger value 
   caused problems at low F10.7 (below 75).

3. A subroutine smoothz has been added to
   smooth the ion velocity along the geomagnetic
   field (line 780).

4. The factor of .01*ne(i,nfl) to limit the
   velocity vsid in the update subroutine
   has been changed to .0001*ne(i,nfl) (line 1398).
   This factor reduces the ion velocity to zero
   for ions that are less than .01% of the
   electron density.

5. The upper boundary condition has been
   changed to a constant flux when the E x B
   velocity is downward (i.e., negative)
   (lines 3149 - 3204).

6. The variable b0 in grid-1.00.f and
   com-1.00.inc has been changed to bb0
   (this variable is not used).

****************************************************************

This file contains the source code description,
tutorial, and graphics sections of the SAMI2 web site.


Source Code Description

  sami2-1.00.f:

    The main code that calculates the evolution
    of the low- to mid-latitude ionosphere. SAMI2
    treats the dynamic plasma and chemical evolution of 
    seven ion species (H+, He+, N+, O+, N2+, NO+, and O2+)
    in the altitude range 85 km to 20,000 km. This corresponds
    to a latitudinal extent of (+/-)62.5 degrees about
    the magnetic equator. The ion continuity and momentum 
    equations can be solved for all 7 species; the temperature 
    equation is solved for H+, He+, O+, and the electrons. 
    SAMI2  models the plasma along the earth's geomagnetic 
    field from hemisphere to hemisphere; an offset, tilted 
    dipole field is used. The code includes a modeled E x B 
    drift of the plasma, as well as ion inertia in the ion 
    momentum equation for motion along the dipole field line. 
    The code uses a fixed, nonorthogonal grid in which one 
    coordinate axis is aligned with the geomagnetic field. 
    The neutral species are specified using the empirical 
    models NRLMSISE00 and HWM93.  

    The primary output variables are the following:

      deni(nz,nf,nion,nt):    ion density
      ti(nz,nf,nion,nt):      ion temperature
      vsi(nz,nf,nion,nt):     ion velocity (along B)
      te(nz,nf,nt):           electron temperature
      time(4,nt):             time step, local time
      glat(nz,nf):            geographic latitude
      glon(nz,nf):            geographic longitude
      zalt(nz,nf):            altitude

    The array indices denote the following:

      nz:          number of mesh points along the geomagnetic field line
      nf:          number of mesh points transverse to the geomagnetic field
                   line (i.e., number of magnetic field lines)
      nion:        number of ion species (default: 7)
      time(4,nt):  time(1,nt): number of time steps
                   time(2,nt): hour
                   time(3,nt): minute
                   time(4,nt): second

    The indexing of the ion specie number 

      1: hydrogen (H)
      2: oxygen (O)
      3: nitrous oxide (NO) 
      4: molecular oxygen (O2)
      5: helium (He)
      6: molecular nitrogen (N2)
      7: nitrogen (N)

  sami2-1.00.namelist:

    The input parameters for SAMI2.

      fmtout:       True: output formatted data files
                    False: output unformatted data files
      maxstep:      The maximum number of time steps allowed.
                    Typically, this is a large number, e.g., 20000000.
                    It is set to smaller number, e.g., 10, just for
                    testing purposes.
      hrmax:        The number of hours for the run (hr). A typical run
                    is for 48 hrs; the first 24 hrs allows transients
                    to clear the system.
      dt0:          The maximum time step allowed (sec). The default
                    is 30 s. This shouldn't be changed.
      dthr:         Defines how often the data is output (hr). The
                    default value is 0.25 (i.e., the data is dumped
                    every 15 min).
      hrpr:         The time period that elapses before the 
                    data is output (hr). This should typically be
                    24 hr. 
      grad_in:      The input altitude (km).
      glat_in:      The input latitude (geographic).
      glon_in:      The input longitude (geographic).
      fejer:        True: use the Fejer/Scherliess empirical E x B drift model
                    False: use the sinusoidal E x B drift model; if this is
                           used then the magnitude of the E x B drift is given
                           by ve01.
      rmin:         Maximum altitude of the lowest field line (km).
                    A typical value is 150 km.
      rmax:         Maximum altitude of the highest field line (km).
                    This has to be less than 20,000 km. 
      altmin:       Altitude of the base of a field line (km).
                    The default is 85 km. 
      fbar:         Value of F10.7A (3 month average of F10.7).
      f10p7:        Value of F10.7.              
      ap:           Value of Ap.
      year:         Year
      day:          Day
      mmass:        Average neutral mass density. The default is 48.
      nion1:        Minimum ion specie index. The default is 1.
      nion2:        Maximum ion specie index. The default is 7.
                    However, one can use 4 and consider only the
                    dominant ions in the ionosphere (H, O, NO, O2).
                    This will speed up the run time of the code
                    by about 30%.
      hrinit:       Local time at the start of the run (hr).
                    The default is 0000 UT.
      gams:         Determines grid spacing along the geomagnetic field.
                    The default is 3. As this parameter is increased,
                    the spacing between grid points along the field
                    line increases at high altitudes. As it is decreased,
                    the spacing becomes more uniform.
      gamp:         Determines grid spacing orthogonal to the geomagnetic
                    field. The default is 3. As this parameter is increased,
                    the spacing between field lines increases at high 
                    altitudes. As it is decreased, the spacing becomes 
                    more uniform.
      tvn0:         Multiplicative factor for the neutral wind speed.
                    The default value is 1. For example, if 
                    tvn0 = 0,  the neutral wind is turned off. 
      tvexb0:       Multiplicative factor for the E x B drift velocity.
                    The default value is 1. For example, if 
                    tvexb0 = 0,  the E x B drift velocity is turned off. 
      ve01:         Maximum E x B drift velocity for the sinusoidal
                    drift model (m/sec). Typical values are
                    5 m/s - 30 m/s.
      snn:          Multiplicative factors for the neutral density.
                    The default values are 1. For example, the neutral
                    O density can be decreased by 75% by setting this
                    parameter to 1.,1.,1.,.75,1.,1.,1.
                     where oxygen is the fourth specie.
      stn:          Multiplicative factor for the neutral temperature.
                    The default value is 1. 
      denmin:       Miniumum ion density allowed. The default value
                    is 1.e-6.
      alt_crit:     The E x B drift is exponentially decreased below
                    this altitude with a scale length 20 km. The
                    default value is 150 km. [This is done to allow
                    rmin to be less than 150 km without using an
                    extremely small time step.]
      cqe:          Constant used in the subroutine 'etemp' associated
                    with photoelectron heating. The typical range is
                    3e-14 -- 8e-14. The higher this value, the lower
                    the electron temperature above 300 km.

  grid-1.00.f:
        
    Sets up the nonorthogonal grid.

  grid-rminrmax.f:

    An auxiliary program that determines the values of 
    rmin and rmax used in sami2-1.00.namelist
    for a given geographic latitude, longitude, and range
    of altitudes at this position.

  param-1.00.inc:

    Parameters used in SAMI2. In general, the only
    parameters that need to be changed are 

    nz (the number of grid points along the geomagnetic field line), and 
    nf (the number of field lines).

    If these parameters are changed, the code must be recompiled.

  com-1.00.inc:

    The SAMI2 variables passed through common statements.

  README-1.00:
   
    An ASCII file that covers the Source Code Description,
    Tutorial, and Graphics sections.

  makesami2-1.00:
    
    A makefile for SAMI2-1.00 that includes the Intel,
    Absoft, Lahey, and Portland Group fortran compilers.

  deni-init.inp:

    Input file that provides initial values of the
    ion densities that are interpolated to the 
    SAMI2 grid.

  euvflux.inp:
    
    Input file that provides the EUV flux. This is
    the EUVAC model developed by Phil Richards
    (Richards et al., J. Geophys. Res., 99, 8981, 1994).

  ichem.inp:

    Input file that identifies the chemically reacting 
    ion and neutral species.

  phabsdt.inp:
  
    Input file that provides the photoabsorption 
    cross-sections associated with O, N2, and O2.
    The data for O are from Bailey and Balan,
    STEP Handbook of Ionospheric Models,
    ed. R. Schunk, p. 184, 1996. The data for N2
    and O2 are from Richards et al., J. Geophys. Res., 
    99, 8981, 1994.

  phiondt.inp:
        
    Input file that provides the daytime photoionization
    cross-sections associated with He, N, O, N2, and O2.
    The data for He are from Bailey and Balan,
    STEP Handbook of Ionospheric Models,
    ed. R. Schunk, p. 184, 1996. The data for N, O, N2,
    and O2 are from Richards et al., J. Geophys. Res.,
    99, 8981, 1994.

  phionnt.inp:
        
    Input file that provides the nighttime photoionization
    cross-sections associated with O, N2, NO, and O2.
    This data is obatained from Oran et al.,
    NRL Memo Report 3984, Naval Research
    Laboratory, Washington, DC, 1979.

  thetant.inp:

    Input file that provides angular information
    regarding nighttime photoionization processes.
    Based on Strobel et al., The nighttime
    ionosphere: E region and lower F
    region, . Geophys. Res., 79, 3171, 1974.

  zaltnt.inp:
        
    Input file that provides altitude information regarding 
    nighttime photoionization processes.
    Based on Strobel et al., The nighttime
    ionosphere: E region and lower F
    region, J. Geophys. Res., 79, 3171, 1974.

  nrlmsise00.f:

    An empirical model that provides the neutral densities
    and temperature. It has been developed by A. Hedin
    and M. Picone. It is an improvement over the MSIS-86
    model [Hedin, J. Geophys. Res., 92d, 4649, 1987].

  hwm93.f:

    An empirical model that provides the neutral wind
    velocity [Hedin et al., <i>J. Geophys. Res., 96</i>,
    7657, 1991].

SAMI2-1.00 Tutorial

  Setting up a code run:

    Only two files need to be edited to set up a code run:
    param-1.00.inc and sami2-1.00.namelist.

    1. param-1.00.inc:
  
       The include file param-1.00.inc contains the parameters 
       nz and nf. These should be the only parameters
       that need to be changed.

         nz: Number of grid points along the geomagnetic field; 
             it must be an odd integer. 
         nf: Number of field lines; 
             it must be greater than or equal to 3 (the 
             first and last field lines are boundary cells).

       A typical value of nz is 201 although the code
       will run with a value 101 as long as rmax < 10,000 km. 
       Increasing the number of grid points along the field 
       decreases the time step, and will lead to longer run times. 
       If any of these parameters is changed the code must be 
       recompiled to effect the change.
    
    2. sami2-1.00.namelist

       The namelist file sami2-1.00.namelist contains the geophysical 
       input conditions. The inputs can be changed without recompiling
       the code. The inputs are described above. For basic runs, the user
       should only change the following parameters:

         fmtout    
         hrmax
         dthr
         hrpr
         grad_in   
         glat_in   
         glon_in   
         fejer     
         rmin      
         rmax      
         fbar      
         f10p7     
         ap        
         year      
         day 
         nion2      
         ve01
         cqe

  Compiling the code:

    SAMI2 is written in Fortran 77. It has been compiled and 
    run using the following Fortran compilers: Absoft, Lahey,
    the Portland Group, and the `free Fortran' compiler g77. 
    Typical compilation lines for
    each compiler are as follows where sami2.x is the 
    executable file. 

    1. Absoft (v7.5)

       f77 hwm93.f nrlmsise00.f grid-1.00.f sami2-1.00.f -s -O -N3 -o sami2.x 

       Absoft (> v9.0)

       f77 hwm93.f nrlmsise00.f grid-1.00.f sami2-1.00.f -s -O -o sami2.x 

    2. Lahey

       lf95 hwm93.f nrlmsise00.f grid-1.00.f sami2-1.00.f --sav -O -o sami2.x

    3. Portland Group

       pgf77 -Msav -fast  hwm93.f nrlmsise00.f grid-1.00.f sami2-1.00.f -o sami2.x

    4. Intel

       ifort hwm93.f nrlmsise00.f grid-1.00.f sami2-1.00.f -O2 -save -o sami2.x

    In the above, the options -s, -Msav, --sav, and -fno-automatic
    save the values of the local variables in subroutines, and the options 
    -O, -fast optimize the program. For Absoft Fortran 77 (v7.5)
    it is also necessary to use the option -N3 so that 
    unformatted data is written in a form that can be read by 
    the IDL procedure read-u.pro discussed in the Graphics Section.

    Note: It is important that SAMI2 be compiled with
          the save or static option.

  Running the code:

    If the code has compiled successfully, simply type
    sami2.x and the code should run. To test the code run the
    example cases described below. The only data files generated 
    by the default version of sami2-1.00.f
    are denif.dat, glatf.dat, glonf.dat, time.dat, zaltf.dat.
    Other data files can be generated by uncommenting 
    the appropriate write statements in the `output' subroutine in 
    sami2-1.00.f.
       
  Test cases:

    Example 1

    This is a short run for a single field line to see if the code
    runs 10 or so time steps and outputs the data in the
    appropriate format. 

    1.  Copy the file example1.namelist to sami2-1.00.namelist.
    2.  Copy the file example1.param to param-1.00.inc.
    3.  Compile the code.
    4.  Run the code. 
    5.  The output to the screen should look something like this

    finished initialization
   istep =   1 ntm =   1 time step =   12.0000  hour =   8.333334E-03
   istep =   2 ntm =   2 time step =   12.0000  hour =   1.166667E-02
   istep =   3 ntm =   3 time step =   12.0000  hour =   1.500000E-02
   istep =   4 ntm =   4 time step =   12.0000  hour =   1.833333E-02
   istep =   5 ntm =   5 time step =   12.0000  hour =   2.166667E-02
   istep =   6 ntm =   6 time step =   12.0000  hour =   2.500000E-02
   istep =   7 ntm =   7 time step =   12.0000  hour =   2.833333E-02
   istep =   8 ntm =   8 time step =   12.0000  hour =   3.166667E-02
   istep =   9 ntm =   9 time step =   12.0000  hour =   3.500000E-02
   istep =   10 ntm =   10 time step =   12.0000  hour =   3.833333E-02
   istep =   11 ntm =   11 time step =   12.0000  hour =   4.166667E-02

        There can be small differences depending on the compiler used.
        The data files (e.g., denif.dat) should be readable using a 
        text editor. 

    6.  Edit the file sami2-1.00.namelist so that
        fmtout = .false. and rerun the code. The output to 
        the screen should be the same but now the output files 
        are in binary format (e.g., deniu.dat).

    Example 2

    This is a longer run for the single field line case.
    It is for the same geophysical conditions
    as Example 1 but is run for 48 hrs with the data dumped
    every 15 minutes. The default is to output formatted data
    files.

    1.  Copy the file example2.namelist to sami2-1.00.namelist.
    2.  Copy the file example2.param to param-1.00.inc.
    3.  Compile the code.
    4.  Run the code. 
        (This takes approximately 16 s on a 2.83 GHZ Intel chip
         using Intel fortran.)
    5.  Examples of the output are presented in the
        Graphics Section of the SAMI2 web site.

    Example 3

    This is a much longer run: the number of field lines is
    nf = 60. The default is to output formatted data files. 
    This run will generate a file denif.dat that is
    approximately 126 MB. This can be changed by editing the
    sami2-1.00.namelist and setting fmtout = .false.
    after Step 1 below. This will generate a file
    deniu.dat that is approximately 31 MB. 
    (Note: The data files must be formatted to use the `free' 
    version of IDL for graphics. See `Some IDL issues' in the
    Graphics Section of the SAMI2 web site.)

    1.  Copy the file example3.namelist to sami2-1.00.namelist.
    2.  Copy the file example3.param to param-1.00.inc.
    3.  Compile the code.
    4.  Run the code. 
        (This takes approximately 5m22s on a 2.83 GHZ Intel chip
         using Intel fortran.)
    5.  Examples of the output are presented in the
        Graphics Section of the SAMI2 web site.


    Example 4

    This example shows how to use the auxiliary code grid-rminrmax.f.
    SAMI2 sets up the grid using the parameters rmin (the altitude 
    of the lowest field line at its magnetic equator) and rmax
    (the altitude of the highest field line at its magnetic equator). 
    The purpose of grid-rminrmax.f is to determine these altitudes 
    by specifying the input parameters glat_in, glon_in, grad_in_min, and
    grad_in_max where grad_in_min and grad_in_max are the minimum 
    and maximum altitudes at a given latitude and longitude of interest.
    
    For instance, say you are interested in doing a simulation
    relevant to Arecibo observations in the altitude range
    100 - 2000 km. The geographic latitude and longitude for
    Arecibo are 18.3 and 293.25, respectively.

    1. Compile the code grid-rminrmax.f
       (e.g., f77 grid-rminrmax.f -o grid-rminrmax.x)
    2. Run the code by typing 
         grid-rminrmax.x
    3. The following should appear: 
         input: grad_in_min,grad_in_max,glat_in,glon_in
    4. Enter the following:
         100. 2000. 18.3 293.25
    5. The following should appear: 
         Input data:
         grad_in_min:  100.000 
         grad_in_max:  2000.00 
         glat_in:  18.3000 
         glon_in:  293.250
         pvalue:   1.31998 
         rmin:   1651.64 
         pvalue:   1.71185 
         rmax:   4149.03
    6. The values of rmin and rmax are now given. These values 
       can be used in the namelist file. However, it is usually best
       to set rmin to a much lower value (e.g., 150 km) and rmax
       to a somewhat large value (e.g., 5000 km) to account for the 
       E x B drift of the plasma.
  
SAMI2 Graphics

  IDL graphics procedures

    We use IDL for graphical output of SAMI2 data. See `IDL Issues'
    at the end of ths file for some information about IDL.
    The following IDL procedures can be used for SAMI2 graphics. 
    Read the comments at the beginning of each procedure for needed 
    inputs at the idl prompt. 

      read-f.pro
 
          This procedure reads in formatted data files.

      read-u.pro

          This procedure reads in unformatted data files.

      contour2d.pro

          This procedure makes a color contour plot of the function f 
          (defined by the user at the idl prompt) 
          as a function of geographic latitude and 
          altitude at a given local time. It 
          plots to the screen and to a postscript file
          contour.ps.

      zalt1d.pro

          This procedure makes a line plot of the function f 
          (defined by the user at the idl prompt) 
          as a function of altitude at
          a given local time and geographic latitude.
          It plots to the screen and to a postscript file 
          zalt1d.ps.

      glat1d.pro

          This procedure makes a line plot of the function f 
          (defined by the user at the idl prompt) 
          as a function of latitude at
          a given local time and altitude.
          It plots to the screen and to a postscript file 
          glat1d.ps.

  Test cases:

    The various IDL generated plots are on the SAMI2 web site.
  
    Example 2
    
    This is an example of graphical output for the data
    generated from Example 2 in the Tutorial
    Section. This is essentially a single field line
    case so the procedures contour.pro,
    zalt1d.pro, glat1d.pro are not useful.

    1. Edit the file read-f.pro so that nf = 3 and nt = 95.
       (Note: The value of nt can be found from the file
       time.dat; it is the first number on the last line.)
    2. Start IDL - enter the commands that are after the IDL.

       IDL> .run read-f
       IDL> ntm=75
       IDL> plot,dene(50:100,1,ntm),zalt(50:100,1),yrange=[0,2000],/xlog,charsize=1.6,xr=[1.e4,1.e7]
            (This is a graph of the electron density 
            as a function of altitude along
            a geomagnetic field that passes over Millstone Hill at
            an altitude 300 km at 1900 UT.)
       IDL> ntm=15
       IDL> plot,dene(50:100,1,ntm),zalt(50:100,1),yrange=[0,2000],/xlog,charsize=1.6,xr=[1.e3,1.e6]
       IDL> oplot,deni(50:100,1,0,ntm),zalt(50:100,1),line=2
       IDL> oplot,deni(50:100,1,1,ntm),zalt(50:100,1),line=3
            (This is a graph of the electron density (solid), 
            O+ density (dash-dot), and H+ density (dash)
            as a function of altitude along
            a geomagnetic field that passes over Millstone Hill at
            an altitude 300 km at about 0400 UT.)

    Example 3
    
    This is an example of graphical output for the data
    generated from Example 3 in the Tutorial Section.

    1. Edit the file read-f.pro or read-u.pro
         so that nf = 60 depending on whether the data 
         file generated is formatted (denif.dat) or
         unformatted (deniu.dat), and nt = 190.
    2. Start IDL - enter the commands that are after the IDL prompt.

       IDL> run read-f or .run read-u
       IDL> f = alog10(dene)
       IDL> ntm = 99
       IDL> .run contour2d
            (This is a color contour plot of the electron density
            as a function of latitude and altitude.)
       IDL> glat0=18.3
       IDL> .run zalt1d
            (This is  a plot of the electron density as a function of
            altitude over Arecibo.)
       IDL> zalt0=350
       IDL> .run glat1d
            (This is a plot of the electron density as a function of
            latitude at a fixed altitude 350 km.)

  Some IDL issues:

     IDL is an excellent (and expensive) graphics software package from
     ITT (http://www.ittvis.com/language/en-US/ProductsServices/IDL.aspx). 
     If you do not have access to a licensed copy, you can download IDL or
     obtain a trial CD from ITT. 
     This trial version only runs for a limited time (7 minutes) 
     and has several limitations (e.g., you cannot write JPEG or GIF 
     files).

     The IDL procedures provided can be edited to obtain 
     plots with titles, axes labels, different scales, etc.
     Also, the procedure contour2d.pro
     can be edited to show grid points, magnetic field lines,
     and contour levels.
  



