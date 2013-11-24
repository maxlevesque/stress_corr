program stress_corr_ext

!************************************************************
!************************************************************

    implicit none

    integer, parameter :: nummax=648
    integer, parameter :: nfilemax=10
    integer, parameter :: ncorrtimemax=500000

    double precision xx,xy,xz,yy,yz,zz
    double precision dtime,dtimerec,dtime2
    double precision du

    integer norm
    integer nfilecf
    integer ncorrtime,ncorr,ncorrcall
    integer nrun,npereng,nstep,npervel,nperfri
    integer num,nanion,ncation,nionunit
    integer nspec
    integer iwrite, nfiles, ntot, ncnt
    integer i,j
    integer nfilecnt
    integer nskip(nfilemax)
 
!================================================
! Double precision arrays in common.
  
    common/coords/xx,xy,xz,yy,yz,zz
    common/cfnorm/norm(0:ncorrtimemax)
 
!================================================
! Double precisions in common.
    common/timest/dtime,dtimerec,dtime2

    common/tttttt/nfilecf(nfilemax)

!================================================
! Integers in common.
    common/corrf1/ncorrtime,ncorr,ncorrcall
    common/runtim/nrun,npereng,nstep,npervel,nperfri
    common/sizeof/num,nanion,ncation,nionunit
    common/specie/nspec

    character(len=80) posfile(nfilemax),posfile2(nfilemax)

!=====================================================================
! Input correlation function matrix.
    open(10,file='stress_corr.inpt',status='old')

       read(10,*) dtime
       read(10,*) iwrite
       dtime=dtime*float(iwrite)*2.418d-05
       read(10,*)nfiles
       write(6,*)nfiles
       do i=1,nfiles

          read(10,'(a)')posfile(i)
          write(6,*) posfile(i)
          read(10,'(a)')posfile2(i)
          write(6,*) posfile2(i)
          read(10,*)nfilecf(i)
          read(10,*)nskip(i)

       enddo
       read(10,*)ncorr

    close(10)

!=====================================================================
    do j=0,ncorr
        norm(j)=0
    enddo

! Auxillary variables.
    ncorrtime=1
    ncorrcall=0

    ntot=0
    ncnt=0
    nfilecnt=1
    do i=1,nfiles
       ntot=ntot+nfilecf(i)

       open(11,file=posfile(i),status='old',form='formatted')
       open(15,file=posfile2(i),status='old',form='formatted')

       do j=1,nfilecf(i)

          write(6,*)j,ncnt,nfilecnt,ncorrcall

          if ((mod(j,nskip(i))).eq.0) then

             read(11,*)du,xy,xz,yz
             read(15,*)du,xx,yy,zz

             call ircfcalc
             ncorrtime=mod(ncorrtime,ncorr)
             ncorrcall=ncorrcall+1
             ncorrtime=ncorrtime+1

             ncnt=ncnt+1

          else

             read(11,*)du,du,du
             read(15,*)du,du,du
             ncnt=ncnt+1

          endif

       enddo

       close(11)
       close(15)

       nfilecnt=nfilecnt+1

    enddo

    call ircfdump

    end program stress_corr_ext


!************************************************************

    subroutine ircfcalc

    implicit none

    integer, parameter :: nummax=648
    integer, parameter :: ncorrtimemax=500000

    double precision xx,xy,xz,yy,yz,zz

    double precision cf1, cf2, cf3, cf4, cf5

    double precision store1, store2, store3, store4, store5

    integer norm
    integer ncorrtime, ncorr, ncorrcall
    integer num, nanion, ncation, nionunit
    integer m, nt

! Double precision arrays in common.
 
    common/coords/xx,xy,xz,yy,yz,zz
 

    common/ircorr/cf1(0:ncorrtimemax)
    common/ircor2/cf2(0:ncorrtimemax)
    common/ircor3/cf3(0:ncorrtimemax)
    common/ircor4/cf4(0:ncorrtimemax)
    common/ircor5/cf5(0:ncorrtimemax)
 
    common/irstor/store1(ncorrtimemax)
    common/irsto2/store2(ncorrtimemax)
    common/irsto3/store3(ncorrtimemax)
    common/irsto3/store4(ncorrtimemax)
    common/irsto3/store5(ncorrtimemax)
 
! Integer arrays in common.
    common/cfnorm/norm(0:ncorrtimemax)

! Integers in common.
    common/corrf1/ncorrtime,ncorr,ncorrcall
    common/sizeof/num,nanion,ncation,nionunit

! Store anion and cation correlation functions.

    store1(ncorrtime)=xy 
    store2(ncorrtime)=xz 
    store3(ncorrtime)=yz 
    store4(ncorrtime)=(xx-yy)/2.0d0 
    store5(ncorrtime)=(2.0d0*zz-xx-yy)/sqrt(12.0d0) 

    do m=1,ncorrtime

       nt=ncorrtime-m

! Normalisation for the correlation functions.
       norm(nt)=norm(nt)+1

! Auto-correlation terms.
          cf1(nt)=cf1(nt)+store1(ncorrtime)*store1(m)
          cf2(nt)=cf2(nt)+store2(ncorrtime)*store2(m) 
          cf3(nt)=cf3(nt)+store3(ncorrtime)*store3(m) 
          cf4(nt)=cf4(nt)+store4(ncorrtime)*store4(m) 
          cf5(nt)=cf5(nt)+store5(ncorrtime)*store5(m) 

    enddo

! Secondary loop.
    if (ncorrcall.gt.ncorr) then

       do m=1+ncorrtime,ncorr

          nt=ncorrtime-m+ncorr

! Normalisation for the correlation functions.
          norm(nt)=norm(nt)+1

! Auto-correlation terms.

             cf1(nt)=cf1(nt)+store1(ncorrtime)*store1(m)
             cf2(nt)=cf2(nt)+store2(ncorrtime)*store2(m)
             cf3(nt)=cf3(nt)+store3(ncorrtime)*store3(m)
             cf4(nt)=cf4(nt)+store4(ncorrtime)*store4(m)
             cf5(nt)=cf5(nt)+store5(ncorrtime)*store5(m) 

       enddo

    endif
!
! array pointers are updated in main routine
!
    return

    end subroutine ircfcalc

!************************************************************

    subroutine ircfdump

!************************************************************

    implicit none

    integer, parameter :: ncorrtimemax=500000

    double precision cf1, cf2, cf3, cf4, cf5
    double precision dtime, dtimerec, dtime2
    double precision t, xnormrec

    integer norm

    integer ncorrtime, ncorr, ncorrcall
    integer j
!================================================
! Double precision arrays.
    common/ircorr/cf1(0:ncorrtimemax)
    common/ircor2/cf2(0:ncorrtimemax)
    common/ircor3/cf3(0:ncorrtimemax)
    common/ircor4/cf4(0:ncorrtimemax)
    common/ircor5/cf5(0:ncorrtimemax)
  
!================================================
! Integer arrays.
    common/cfnorm/norm(0:ncorrtimemax)

!================================================
! Integers in common.
    common/corrf1/ncorrtime,ncorr,ncorrcall
    common/timest/dtime,dtimerec,dtime2

    open(60,file='scorr.out_xy')
    open(61,file='scorr.out_xz')
    open(62,file='scorr.out_yz')
    open(63,file='scorr.out_xx-yy')
    open(64,file='scorr.out_2zz-xx-yy')

    do j=0,ncorr

       if (norm(j).ne.0) then
          xnormrec=1.0d0/float(norm(j))
       else
          xnormrec=0.0d0
       endif
          t=j*dtime
          write(60,*)t,cf1(j)*xnormrec
          write(61,*)t,cf2(j)*xnormrec
          write(62,*)t,cf3(j)*xnormrec
          write(63,*)t,cf4(j)*xnormrec
          write(64,*)t,cf5(j)*xnormrec
    enddo

    close(60)
    close(61)
    close(62)
    close(63)
    close(64)

    write(6,*)
    write(6,*)'**** Finished correlation functions written out ****'
    write(6,*)

    return

end subroutine ircfdump
 
