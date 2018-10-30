!*********************************************************************
!x=windowed data
!isize=window length
!zlyap=output lmax values
!C invocation syntax
! lamd_(double*,uunsigned int*,double*,int*,int*,float*,double*,double*,int*,double*,int*)
!SUBROUTINE LAMD(x,isize,zlyap,DIMM_,tau_,dt_,bb_,cc_,mult_,angmx1_,fiduc1_)

!        SUBROUTINE LAMD(x,isize,zlyap)
        SUBROUTINE LAMD(x,isize,zlyap,DIMMi,taui,dti,bbi,cci,multi,angmx1i,fiduc1i)
          
        IMPLICIT REAL*8(A-H,O-Z)

        !INTEGER DIMM,tau,EVOLV,zmult,fiduc1,fiduc2,difind
	!	real dt

        !INTEGER DIMMi,taui,fiduc1i,multi
        !        real dti

        INTEGER DIMM,tau,EVOLV,zmult,fiduc1,fiduc2,difind,DIMMi,taui,fiduc1i,multi
        real dt,dti

        DIMENSION x(500048),pt1(20),pt2(20),z(500000,20)

!C       DIMENSION y(500048),z1(500000,20)

	common/abc/ DIMM,tau,EVOLV,dt,IDIST

        !bb=0.2
        bb=bbi

        !cc=0.03
        cc=cci

        !mult=5
        mult=multi

        !angmx1=0.2
        angmx1=angmx1i

	!fiduc1=1
        fiduc1=fiduc1i

        !dt=0.005 !sampling time
        dt=dti

	!DIMM=7 !embedding dimension
        DIMM=DIMMi

        !tau=4 !time delay
        tau=taui

        EVOLV=(DIMM-1)*TAU/2 !evolution distance

        IDIST=(DIMM-1)*tau

!*

!*Construction of the Z array of the phase space

!*

        imax=isize-(dimm-1)*tau

        do  i=1,imax

        do  j=1,dimm

        z(i,j)=x(i+(j-1)*tau)

!C       z1(i,j)=y(i+(j-1)*tau)

	enddo

	enddo







!*

!*Calculate useful size of datafile

!*

       npt=isize-dimm*tau-evolv



!*

!*First fiducial point

!*Find nearest neighbor to first data point

!*



       ind=fiduc1



!******* dont take a point too close in time to the first fiducial point

       di=1.d38

       do 30 i=(dimm-1)*tau,npt

       d=0.0d0

       do 20 j=1,dimm

20     d=d+(z(ind,j)-z(i,j))**2

       d=dsqrt(d)

       if (d.ge.di) goto 30

       di=d

       ind2=i

30     continue

!c-----------------------------



       zmult=1

       anglmx=angmx1

       sum=0.0d0

       its=0



!*

!*Get coordinates of evolved points and compute the new

!*separation length between them

!*

40     continue



        do 50 j=1,dimm

        pt1(j)=z(ind+evolv,j)

        pt2(j)=z(ind2+evolv,j)

50      continue



        df=0.0d0

        do 60 j=1,dimm

        df=df+(pt1(j)-pt2(j))**2

60      continue

        df=dsqrt(df)

!*

!*Compute (w.r. to power of 2 ) and update the Lyapunov exponent

!*

       its=its+1



       if (di.eq.0.0d0) di=1.0d-7

       if (df.eq.0.0d0) df=1.0d-7



       dl=df/di



       sump=dlog(dl)/(dfloat(evolv)*dt)

       sum=sum+sump





       zlyap=sum/dfloat(its)



       fiduc2=ind



!*Find out the replacement point for the evolved secondary one

!*Search over longer distances which must always be less than

!*the product ZMULT*SCALMX and greater than noise level SCALMN



       indold=ind2

       zmult=1

       anglmx=angmx1



!*****Compute distance between fiducial point and candidate

!*****   and the magnitude of the fiducial point vector

       dnewmx=0.0d0

       do 69  i=1,npt

       difind=iabs(i-ind-evolv)

!*------------------------------------------------------------------

       if ((difind.le.tau).or.(difind.gt.IDIST)) goto 69

!ccc

!CCC      Having the above command as a comment and the following one in,

!CCC   then normalisation over the whole segment - Great effect of

!CCC   prominent spikes that may be present - Use ANGMX=0.7 (larger L

!CCC   in value, save money with no difference at the spikes POSITION)

!*-------------------------------------------------------------------



          dnew=0.0d0

          do 65 j=1,dimm

          dnew=dnew+(pt1(j)-z(i,j))**2

65        continue

          dnew=dsqrt(dnew)



          if (dnew.gt.dnewmx) dnewmx=dnew

69        continue



!********* Values for SCALMN,SCALMX

70      thmin=3.1415d0

        scalmn=dnewmx*cc

        scalmx=dnewmx*bb*zmult

!*---------------------------------------------------------------------

        do 100  i=1,npt

        difind=iabs(i-ind-evolv)

        if (difind.LE.(dimm-1)*tau) goto 100

!*****Compute distance between fiducial point and candidate and compare

          dnew=0.0d0

          do 80 j=1,dimm

          dnew=dnew+(pt1(j)-z(i,j))**2

80        continue

          dnew=dsqrt(dnew)

          if ((dnew.gt.scalmx).or.(dnew.lt.scalmn)) goto 100



!*****Find angular change old to new vector

         dot=0.0d0

         do 90 j=1,dimm

          dot=dot+(pt1(j)-z(i,j))*(pt1(j)-pt2(j))

90       continue



        if (df.eq.0.0d0) df=1.0d-5

        if (dnew.eq.0.0d0) dnew=1.0d-5

        cth=(dot/(dnew*df))

        if (cth.gt.1.0 .or. cth.lt.-1.0 ) cth=1.0d0

        th=dacos(cth)

!*Save point with smallest angular change AND distance so far

        if (th.gt.thmin) goto 100

          thmin=th

          dii=dnew

          ind2=I

100    continue



        if (thmin.lt.anglmx) goto 110

!*

!* I cant find a replacement with a proper angle deviation

!* so I look at longer distances

!*

       zmult=zmult+1

       if (zmult.le.mult) goto 70

!*

!* I cant find a replacement with a proper angle deviation

!* with 5*scale distance so I double the search angle and

!* start from the beginning w.r.t. the distance from the

!* fudicial point

!*

       zmult=1

       anglmx=2*anglmx

       if (anglmx.le.1.0d0) goto 70

       dii=df

       ind2=indold+evolv

110    ind=ind+evolv



!*

!* Program ends when fiducial trajectory hits end of file

!*

        if (ind.ge.npt) goto 1000

        di=dii

        goto 40



1000     return

        END





      

!*********************************************************

!*********************************************************
