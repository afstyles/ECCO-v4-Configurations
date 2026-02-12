C $Header: /u/gcmpack/MITgcm/pkg/ecco/ecco_phys.F,v 1.18 2017/04/03 23:16:38 ou.wang Exp $
C $Name:  $

#include "ECCO_OPTIONS.h"

      SUBROUTINE ECCO_PHYS( myTime, myIter, myThid )

c     ==================================================================
c     SUBROUTINE ecco_phys
c     ==================================================================
c
c     ==================================================================
c     SUBROUTINE ecco_phys
c     ==================================================================

      implicit none

c     == global variables ==

#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "GMREDI_OPTIONS.h"
#include "GMREDI.h"
#include "GRID.h"
#ifdef ALLOW_ECCO
# include "ecco.h"
#endif
#ifdef ALLOW_PTRACERS
# include "PTRACERS_SIZE.h"
# include "PTRACERS_FIELDS.h"
#endif

c     == routine arguments ==

      _RL     myTime
      INTEGER myIter, myThid

c     == local variables ==

      INTEGER bi,bj
      INTEGER i,j,k
      INTEGER jmin,jmax
      INTEGER imin,imax
#ifdef ALLOW_GENCOST_CONTRIBUTION
      INTEGER kgen, kgen3d, itr, kp1, maskp1
      _RL areavolTile(nSx,nSy), areavolGlob
      _RL tmpfld, tmpvol, tmpmsk, tmpmsk2, tmpmskW, tmpmskS
      _RL tmp_sigmsk, tmpsig, tmpsig_lower, tmpsig_upper
      _RL delPsi, uRes, vRes, tmpfldS, tmpfldW
#endif

c- note defined with overlap here, not needed, but more efficient
      _RL trVolW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL trVolS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL trHeatW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL trHeatS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL trSaltW(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL trSaltS(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)

#ifdef ATMOSPHERIC_LOADING
#ifdef ALLOW_IB_CORR
      _RL AREAsumTile(nSx,nSy),PLOADsumTile(nSx,nSy)
      _RL tmpfac
      CHARACTER*(MAX_LEN_MBUF) msgBuf
#endif
      _RL sIceLoadFac
#endif
#ifdef ALLOW_PSBAR_STERIC
      _RL RHOInSituLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL VOLsumTile(nSx,nSy),RHOsumTile(nSx,nSy)
#endif

c need to include halos for find_rho_2d
      iMin = 1-OLx
      iMax = sNx+OLx
      jMin = 1-OLy
      jMax = sNy+OLy

#ifdef ALLOW_PSBAR_STERIC

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
          do k = 1,nr
            CALL FIND_RHO_2D(
     I                iMin, iMax, jMin, jMax, k,
     I                theta(1-OLx,1-OLy,k,bi,bj),
     I                salt (1-OLx,1-OLy,k,bi,bj),
     O                RHOInSituLoc(1-OLx,1-OLy,k,bi,bj),
     I                k, bi, bj, myThid )
          enddo
       enddo
      enddo

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
          RHOsumTile(bi,bj)=0. _d 0
          VOLsumTile(bi,bj)=0. _d 0
          VOLsumGlob=0. _d 0
          RHOsumGlob=0. _d 0
          do k = 1,nr
            do j = 1,sNy
              do i =  1,sNx
                RHOsumTile(bi,bj)=RHOsumTile(bi,bj)+
     &            (rhoConst+RHOInSituLoc(i,j,k,bi,bj))*
     &            hFacC(i,j,k,bi,bj)*drF(k)*rA(i,j,bi,bj)
                VOLsumTile(bi,bj)=VOLsumTile(bi,bj)+
     &            hFacC(i,j,k,bi,bj)*drF(k)*rA(i,j,bi,bj)
              enddo
            enddo
          enddo
        enddo
      enddo
      CALL GLOBAL_SUM_TILE_RL( VOLsumTile, VOLsumGlob, myThid )
      CALL GLOBAL_SUM_TILE_RL( RHOsumTile, RHOsumGlob, myThid )
      RHOsumGlob=RHOsumGlob/VOLsumGlob

      if (RHOsumGlob_0.GT.0. _d 0) then
        sterGloH=VOLsumGlob_0/globalArea
     &        *(1. _d 0 - RHOsumGlob/RHOsumGlob_0)
      else
        sterGloH=0. _d 0
      endif

      WRITE(msgBuf,'(A,1PE21.14)') ' sterGloH= ', sterGloH
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )

#endif

#ifdef ATMOSPHERIC_LOADING
#ifdef ALLOW_IB_CORR
      tmpfac = recip_rhoConst*recip_gravity
      ploadbar = 0. _d 0
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
          PLOADsumTile(bi,bj)=0. _d 0
          AREAsumTile(bi,bj)=0. _d 0
          AREAsumGlob=0. _d 0
          PLOADsumGlob=0. _d 0
          ploadbar = 0. _d 0
          do j = 1,sNy
            do i =  1,sNx
              PLOADsumTile(bi,bj)=PLOADsumTile(bi,bj)+
     &          pload(i,j,bi,bj)*
     &          maskC(i,j,1,bi,bj)*rA(i,j,bi,bj)
              AREAsumTile(bi,bj)=AREAsumTile(bi,bj)+
     &          maskC(i,j,1,bi,bj)*rA(i,j,bi,bj)
            enddo
          enddo
        enddo
      enddo
      CALL GLOBAL_SUM_TILE_RL( AREAsumTile, AREAsumGlob, myThid )
      CALL GLOBAL_SUM_TILE_RL( PLOADsumTile, PLOADsumGlob, myThid )
      ploadbar=PLOADsumGlob/AREAsumGlob

      WRITE(msgBuf,'(A,1PE21.14)') ' ploadbar= ', ploadbar
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
#endif

      sIceLoadFac=zeroRL
      IF ( useRealFreshWaterFlux ) sIceLoadFac=recip_rhoConst
#endif

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
            do j = jmin,jmax
              do i =  imin,imax
                m_eta(i,j,bi,bj)=
     &                (etan(i,j,bi,bj)
#ifdef ATMOSPHERIC_LOADING
     &                +sIceLoad(i,j,bi,bj)*sIceLoadFac
#endif
#ifdef ALLOW_PSBAR_STERIC
     &                +sterGloH
#endif
     &                ) * maskC(i,j,1,bi,bj)

                m_bp(i,j,bi,bj)=
     &                (phiHydLow(i,j,bi,bj)
#ifdef ALLOW_PSBAR_STERIC
C add back the correction due to the global mean steric ssh change,
C     i.e. sterGloH computed in ecco_phys.F (units converted from m to m2/s2)
     &                 +sterGloH * gravity
#endif
     &                ) * maskC(i,j,1,bi,bj)

#ifdef ATMOSPHERIC_LOADING
#ifdef ALLOW_IB_CORR
                m_eta_ib(i,j,bi,bj)=
#ifdef INCLUDE_PLOADBAR
     &           (-pload(i,j,bi,bj))*tmpfac
#else
     &           (ploadbar-pload(i,j,bi,bj))*tmpfac
#endif
     &           * maskC(i,j,1,bi,bj)
                m_eta_dyn(i,j,bi,bj)=
     &           (m_eta(i,j,bi,bj) - m_eta_ib(i,j,bi,bj))
     &           * maskC(i,j,1,bi,bj)

                m_bp_nopabar(i,j,bi,bj)=
     &                (phiHydLow(i,j,bi,bj)
#ifdef ALLOW_PSBAR_STERIC
C add back the correction due to the global mean steric ssh change,
C     i.e. sterGloH computed in ecco_phys.F (units converted from m to m2/s2)
     &                 +sterGloH * gravity
#ifdef INCLUDE_PLOADBAR
#else
     &                 - ploadbar * recip_rhoConst
#endif
#endif
     &                ) * maskC(i,j,1,bi,bj)
#endif
#endif
              enddo
            enddo
        enddo
      enddo

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics .AND. myIter.GE.0 ) THEN
        CALL DIAGNOSTICS_FILL( m_eta, 'SSHNOIBC', 0,1, 0,1,1, myThid )
        CALL DIAGNOSTICS_SCALE_FILL( m_bp, recip_gravity, 1,
     &                         'OBPGMAP ', 0,1, 0,1,1, myThid )
#ifdef ATMOSPHERIC_LOADING
#ifdef ALLOW_IB_CORR
        CALL DIAGNOSTICS_FILL( m_eta_ib,
     &                         'SSHIBC  ', 0,1, 0,1,1, myThid )
        CALL DIAGNOSTICS_FILL( m_eta_dyn,
     &                         'SSH     ', 0,1, 0,1,1, myThid )
        CALL DIAGNOSTICS_SCALE_FILL( m_bp_nopabar, recip_gravity, 1,
     &                         'OBP     ', 0,1, 0,1,1, myThid )
#endif /* ALLOW_IB_CORR */
#endif /* ATMOSPHERIC_LOADING */
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
          do k = 1,nr
            do j = 1,sNy
              do i =  1,sNx
                m_UE(i,j,k,bi,bj)=0. _d 0
                m_VN(i,j,k,bi,bj)=0. _d 0
              enddo
            enddo
          enddo
        enddo
      enddo

      CALL ROTATE_UV2EN_RL(
     U          uVel, vVel, m_UE, m_VN,
     I          .TRUE., .TRUE., .FALSE., Nr, mythid )


c--   trVol : volume flux    --- [m^3/sec] (order of 10^6 = 1 Sv)
c--   trHeat: heat transport --- [Watt] (order of 1.E15 = PW)
c--   trSalt: salt transport --- [kg/sec] (order 1.E9 equiv. 1 Sv in vol.)
c--       convert from [ppt*m^3/sec] via rhoConst/1000.
c--       ( 1ppt = 1000*[mass(salt)]/[mass(seawater)] )

c-- init
      call ecco_zero(trVol,Nr,zeroRL,myThid)
      call ecco_zero(trHeat,Nr,zeroRL,myThid)
      call ecco_zero(trSalt,Nr,zeroRL,myThid)

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
c-- init: if done with overlap, more efficient. But since overwritten
c immediately afterward, init is probably not needed.
          do k = 1,nr
            do j = 1-OLy,sNy+Oly
              do i = 1-OLx,sNx+OLx
                trVolW(i,j,k)=0. _d 0
                trVolS(i,j,k)=0. _d 0
                trHeatW(i,j,k)=0. _d 0
                trHeatS(i,j,k)=0. _d 0
                trSaltW(i,j,k)=0. _d 0
                trSaltS(i,j,k)=0. _d 0
              enddo
            enddo
          enddo
          do k = 1,nr
            do j = 1,sNy
              do i =  1,sNx
                trVolW(i,j,k) =
     &                 uVel(i,j,k,bi,bj)*hFacW(i,j,k,bi,bj)
     &                *dyG(i,j,bi,bj)*drF(k)*msktrVolW(i,j,bi,bj)
     &                *maskInW(i,j,bi,bj)
                trVolS(i,j,k) =
     &                 vVel(i,j,k,bi,bj)*hFacS(i,j,k,bi,bj)
     &                *dxG(i,j,bi,bj)*drF(k)*msktrVolS(i,j,bi,bj)
     &                *maskInS(i,j,bi,bj)

                trHeatW(i,j,k) = trVolW(i,j,k)
     &                *(theta(i,j,k,bi,bj)+theta(i-1,j,k,bi,bj))*halfRL
     &                *HeatCapacity_Cp*rhoConst
                trHeatS(i,j,k) = trVolS(i,j,k)
     &                *(theta(i,j,k,bi,bj)+theta(i,j-1,k,bi,bj))*halfRL
     &                *HeatCapacity_Cp*rhoConst

                trSaltW(i,j,k) = trVolW(i,j,k)
     &                *(salt(i,j,k,bi,bj)+salt(i-1,j,k,bi,bj))*halfRL
     &                *rhoConst/1000.
                trSaltS(i,j,k) = trVolS(i,j,k)
     &                *(salt(i,j,k,bi,bj)+salt(i,j-1,k,bi,bj))*halfRL
     &                *rhoConst/1000.
c now summing
                trVol(i,j,k,bi,bj)=trVolW(i,j,k)+trVolS(i,j,k)
                trHeat(i,j,k,bi,bj)=trHeatW(i,j,k)+trHeatS(i,j,k)
                trSalt(i,j,k,bi,bj)=trSaltW(i,j,k)+trSaltS(i,j,k)
                
              enddo
            enddo
          enddo
        enddo
      enddo

#ifdef ALLOW_GENCOST_CONTRIBUTION

      do kgen=1,NGENCOST

      kp1 = MIN(k+1,nr)
      maskp1 = 1.
      IF (k.GE.nr) maskp1 = 0.
      
      itr = gencost_itracer(kgen)

      call ecco_zero(gencost_storefld(1-OLx,1-OLy,1,1,kgen),
     &     1,zeroRL,myThid)

      do bj=myByLo(myThid),myByHi(myThid)
       do bi=myBxLo(myThid),myBxHi(myThid)
         areavolTile(bi,bj)=0. _d 0
       enddo
      enddo
      areavolGlob=0. _d 0

      do bj=myByLo(myThid),myByHi(myThid)
       do bi=myBxLo(myThid),myBxHi(myThid)
        do j = 1,sNy
         do i =  1,sNx
c---------
          do k = 1,nr
            tmpvol=hFacC(i,j,k,bi,bj)*drF(k)*rA(i,j,bi,bj)
c
            tmpmsk=0. _d 0
            if (.NOT.gencost_msk_is3d(kgen)) then
              tmpmsk=gencost_mskCsurf(i,j,bi,bj,kgen)*
     &               gencost_mskVertical(k,kgen)
#ifdef ALLOW_GENCOST3D
            else
              kgen3d=gencost_msk_pointer3d(kgen)
              tmpmsk=gencost_mskC(i,j,k,bi,bj,kgen3d)
#endif
            endif

C ---- If density mask is enabled, use it here ----
            IF ( maskC(i,j,k,bi,bj).EQ.oneRS .AND.
     &           gencost_useDensityMask(kgen) ) THEN
C            - first, calculate the scalar density
             CALL FIND_RHO_SCALAR(
     I              theta(i,j,k,bi,bj),
     I              salt(i,j,k,bi,bj),
     I              gencost_refPressure(kgen),
     O              tmpsig,
     I              myThid )
C            - subtract 1000 to get sigma
             tmpsig = tmpsig - 1000. _d 0
C            - now, tmpmsk is sigmoid times this value
             tmpsig_lower = 0.5 + 0.5*tanh(gencost_tanhScale(kgen)
     &           *(tmpsig-gencost_sigmaLow(kgen)))
             tmpsig_upper = 0.5 - 0.5*tanh(gencost_tanhScale(kgen)
     &           *(tmpsig-gencost_sigmaHigh(kgen)))
C             - update mask value based on the sigmoid function
             tmp_sigmsk = tmpsig_lower*tmpsig_upper
             tmpmsk = tmpmsk*tmp_sigmsk
            ENDIF
C ---- end of density mask (but tmpmsk is used below)

c
            tmpfld=0. _d 0
            tmpmsk2=0. _d 0
            if (gencost_barfile(kgen)(1:15).EQ.'m_boxmean_theta') then
              tmpfld=theta(i,j,k,bi,bj)
              if (tmpmsk.NE.0. _d 0) tmpmsk2=1. _d 0
            elseif (gencost_barfile(kgen)(1:14).EQ.'m_boxmean_salt')
     &        then
              tmpfld=salt(i,j,k,bi,bj)
              if (tmpmsk.NE.0. _d 0) tmpmsk2=1. _d 0
#ifdef ALLOW_PTRACERS
            elseif (gencost_barfile(kgen)(1:17).EQ.'m_boxmean_ptracer')
     &        then
              tmpfld=pTracer(i,j,k,bi,bj,itr)
              if (tmpmsk.NE.0. _d 0) tmpmsk2=1. _d 0
#endif
            endif
c
            gencost_storefld(i,j,bi,bj,kgen) =
     &          gencost_storefld(i,j,bi,bj,kgen)
     &          +tmpmsk*tmpfld
            areavolTile(bi,bj)=areavolTile(bi,bj)
     &          +tmpmsk2*eccoVol_0(i,j,k,bi,bj)
c
          enddo
c---------
          tmpmsk=maskC(i,j,1,bi,bj)*gencost_mskCsurf(i,j,bi,bj,kgen)
          tmpfld=0. _d 0
          tmpmsk2=0. _d 0
          if (gencost_barfile(kgen)(1:13).EQ.'m_boxmean_eta' .OR.
     &        gencost_barfile(kgen)(1:13).EQ.'m_boxmean_obp') then
           if (gencost_barfile(kgen)(1:17).EQ.'m_boxmean_eta_dyn') then
            tmpfld=m_eta_dyn(i,j,bi,bj)
           else if (gencost_barfile(kgen)(1:13).EQ.'m_boxmean_eta')
     &      then
            tmpfld=m_eta(i,j,bi,bj)
           else if (gencost_barfile(kgen)(1:17).EQ.'m_boxmean_obpgmap')
     &      then
c OBPGMAP (in m; converted from m^2/s^2 to m)
            tmpfld=m_bp(i,j,bi,bj)*recip_gravity
           else if (gencost_barfile(kgen)(1:13).EQ.'m_boxmean_obp')
     &      then
c OBP (in m; converted from m^2/s^2 to m)
            tmpfld=m_bp_nopabar(i,j,bi,bj)*recip_gravity
           endif
           if (tmpmsk.NE.0. _d 0) tmpmsk2=1. _d 0
          endif
c
          gencost_storefld(i,j,bi,bj,kgen) =
     &        gencost_storefld(i,j,bi,bj,kgen)
     &        +tmpmsk*tmpfld
c---------
          do k = 1,nr
c
            tmpmskW=0. _d 0
            tmpmskS=0. _d 0
            if (.NOT.gencost_msk_is3d(kgen)) then
              tmpmskW=gencost_mskWsurf(i,j,bi,bj,kgen)
     &          *gencost_mskVertical(k,kgen)
              tmpmskS=gencost_mskSsurf(i,j,bi,bj,kgen)
     &          *gencost_mskVertical(k,kgen)
#ifdef ALLOW_GENCOST3D
            else
              kgen3d=gencost_msk_pointer3d(kgen)
              tmpmskW=gencost_mskW(i,j,k,bi,bj,kgen3d)
              tmpmskS=gencost_mskS(i,j,k,bi,bj,kgen3d)
#endif
            endif
C ---- If density mask is enabled, use it here ----
            IF (( maskW(i,j,k,bi,bj).EQ.oneRS .OR.
     &            maskS(i,j,k,bi,bj).EQ.oneRS) .AND.
     &           gencost_useDensityMask(kgen) ) THEN
C            - first, calculate the scalar density
             CALL FIND_RHO_SCALAR(
     I              theta(i,j,k,bi,bj),
     I              salt(i,j,k,bi,bj),
     I              gencost_refPressure(kgen),
     O              tmpsig,
     I              myThid )
C            - subtract 1000 to get sigma
             tmpsig = tmpsig - 1000. _d 0
C            - now, tmpmsk is sigmoid times this value
             tmpsig_lower = 0.5 + 0.5*tanh(gencost_tanhScale(kgen)
     &           *(tmpsig-gencost_sigmaLow(kgen)))
             tmpsig_upper = 0.5 - 0.5*tanh(gencost_tanhScale(kgen)
     &           *(tmpsig-gencost_sigmaHigh(kgen)))
C             - update mask value based on the sigmoid function
             tmp_sigmsk = tmpsig_lower*tmpsig_upper
             tmpmskW = tmpmskW*tmp_sigmsk
             tmpmskS = tmpmskS*tmp_sigmsk
            ENDIF
c
            if (gencost_barfile(kgen)(1:13).EQ.'m_horflux_vol') then
              gencost_storefld(i,j,bi,bj,kgen) =
     &          gencost_storefld(i,j,bi,bj,kgen)
     &          +uVel(i,j,k,bi,bj)*tmpmskW
     &          +vVel(i,j,k,bi,bj)*tmpmskS
            endif

            tmpfldW=0. _d 0
            tmpfldS=0. _d 0
            tmpmsk2=0. _d 0
            if (gencost_barfile(kgen)(1:15).EQ.'m_horflux_theta')
     &        then
              tmpfldW=(theta(i,j,k,bi,bj)+theta(i-1,j,k,bi,bj))*halfRL
     &           *HeatCapacity_Cp*rhoConst
              tmpfldS=(theta(i,j,k,bi,bj)+theta(i,j-1,k,bi,bj))*halfRL
     &           *HeatCapacity_Cp*rhoConst
              tmpmsk2=1. _d 0
            elseif (gencost_barfile(kgen)(1:14).EQ.'m_horflux_salt')
     &         then
              tmpfldW=(salt(i,j,k,bi,bj)+salt(i-1,j,k,bi,bj))*halfRL
              tmpfldS=(salt(i,j,k,bi,bj)+salt(i,j-1,k,bi,bj))*halfRL
              tmpmsk2=1. _d 0
            elseif ((gencost_barfile(kgen)(1:12).EQ.'m_horflux_fw')
     &               .AND.(gencost_sRef(kgen).NE.0. _d 0))
     &         then
              tmpfldW=(gencost_sRef(kgen)-
     &        (salt(i,j,k,bi,bj)+salt(i-1,j,k,bi,bj))*halfRL)
     &        /gencost_sRef(kgen)
              tmpfldS=(gencost_sRef(kgen)-
     &        (salt(i,j,k,bi,bj)+salt(i,j-1,k,bi,bj))*halfRL)
     &        /gencost_sRef(kgen)
              tmpmsk2=1. _d 0
            endif

            if (tmpmsk2.EQ.1.) then
               kp1 = MIN(k+1,nr)
               maskp1 = 1.
               IF (k.GE.nr) maskp1 = 0.
            
               delPsi = GM_PsiX(i,j,kp1,bi,bj)*maskp1
     &            - GM_PsiX(i,j, k, bi,bj)
               uRes = uVel(i,j,k,bi,bj)
     &            + delPsi*recip_drF(k)*_recip_hFacW(i,j,k,bi,bj) 
               delPsi = GM_PsiY(i,j,kp1,bi,bj)*maskp1
     &            - GM_PsiY(i,j, k, bi,bj)
               vRes = vVel(i,j,k,bi,bj)
     &            + delPsi*recip_drF(k)*_recip_hFacS(i,j,k,bi,bj)

               gencost_storefld(i,j,bi,bj,kgen) =
     &         gencost_storefld(i,j,bi,bj,kgen)
     &          +uRes*tmpmskW*tmpfldW
     &          +vRes*tmpmskS*tmpfldS  

            endif
          enddo
c---------
         enddo
        enddo
       enddo
      enddo

      enddo

#endif /* ALLOW_GENCOST_CONTRIBUTION */


      return
      end
