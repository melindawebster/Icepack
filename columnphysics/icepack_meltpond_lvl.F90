!=======================================================================

! Level-ice meltpond parameterization
!
! This meltpond parameterization was developed for use with the delta-
! Eddington radiation scheme, and only affects the radiation budget in
! the model.  That is, although the pond volume is tracked, that liquid
! water is not used elsewhere in the model for mass budgets or other
! physical processes.
!
! authors Elizabeth Hunke (LANL)
!         David Hebert (NRL Stennis)
!         Olivier Lecomte (Univ. Louvain)

      module icepack_meltpond_lvl

      use icepack_kinds
      use icepack_parameters, only: c0, c1, c2, c10, p01, p5, puny
      use icepack_parameters, only: viscosity_dyn, rhoi, rhos, rhow, Timelt, Tffresh, Lfresh
      use icepack_parameters, only: gravit, depressT, rhofresh, kice, pndaspect, use_smliq_pnd
      use icepack_parameters, only: ktherm, frzpnd, dpscale, hi_min
      use icepack_parameters, only: pndhyps
      use icepack_tracers,    only: nilyr
      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      implicit none

      private
      public :: compute_ponds_lvl, pond_hypsometry

!=======================================================================

      contains

!=======================================================================

      subroutine compute_ponds_lvl(dt,                   &
                                   rfrac,  meltt, melts, &
                                   frain,  Tair,  fsurfn,&
                                   dhs,    ffrac,        &
                                   aicen,  vicen, vsnon, &
                                   qicen,  sicen,        &
                                   Tsfcn,  alvl,         &
                                   apnd,   hpnd,  ipnd,  &
                                   meltsliqn, frpndn,    &
                                   rfpndn, ilpndn)

      real (kind=dbl_kind), intent(in) :: &
         dt          ! time step (s)

      real (kind=dbl_kind), intent(in) :: &
         Tsfcn, &    ! surface temperature (C)
         alvl,  &    ! fraction of level ice
         rfrac, &    ! water fraction retained for melt ponds
         meltt, &    ! top melt rate (m/s)
         melts, &    ! snow melt rate (m/s)
         frain, &    ! rainfall rate (kg/m2/s)
         Tair,  &    ! air temperature (K)
         fsurfn,&    ! atm-ice surface heat flux  (W/m2)
         aicen, &    ! ice area fraction
         vicen, &    ! ice volume (m)
         vsnon, &    ! snow volume (m)
         meltsliqn   ! liquid contribution to meltponds in dt (kg/m^2)

      real (kind=dbl_kind), intent(inout) :: &
         apnd, hpnd, ipnd, &
         frpndn, &   ! pond drainage rate due to freeboard constraint (m/step)
         rfpndn, &   ! runoff rate due to rfrac (m/step)
         ilpndn      ! pond loss/gain due to ice lid (m/step)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         qicen, &  ! ice layer enthalpy (J m-3)
         sicen     ! salinity (ppt)

      real (kind=dbl_kind), intent(in) :: &
         dhs       ! depth difference for snow on sea ice and pond ice

      real (kind=dbl_kind), intent(out) :: &
         ffrac     ! fraction of fsurfn over pond used to melt ipond

      ! local temporary variables

      real (kind=dbl_kind) :: &
         dhpond, &    ! change in hpond (m)
         dvn_temp     ! local variable for change in volume due to rfrac

      real (kind=dbl_kind), dimension (nilyr) :: &
         Tmlt      ! melting temperature (C)

      real (kind=dbl_kind) :: &
         hi                     , & ! ice thickness (m)
         hs                     , & ! snow depth (m)
         dTs                    , & ! surface temperature diff for freeze-up (C)
         Tp                     , & ! pond freezing temperature (C)
         Ts                     , & ! surface air temperature (C)
         apondn                 , & ! local pond area
         hpondn                 , & ! local pond depth (m)
         vpondn                 , & ! pond volume per category area (m)
         dvpondn                , & ! change in pond volume per category area (m)
         hlid                   , & ! refrozen lid thickness
         dhlid                  , & ! change in refrozen lid thickness
         bdt                    , & ! 2 kice dT dt / (rhoi Lfresh)
         alvl_tmp               , & ! level ice fraction of ice area
         draft, deltah, pressure_head, perm, drain ! for permeability

      real (kind=dbl_kind), parameter :: &
         Td       = c2          , & ! temperature difference for freeze-up (C)
         rexp     = p01             ! pond contraction scaling

      character(len=*),parameter :: subname='(compute_ponds_lvl)'

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      vpondn = hpnd * alvl * apnd
      ffrac = c0

      !-----------------------------------------------------------------
      ! Identify grid cells where ponds can be
      !-----------------------------------------------------------------

      if (aicen*alvl > puny**2) then

         hi = vicen/aicen
         hs = vsnon/aicen
         alvl_tmp = alvl

         if (hi < hi_min) then

            !--------------------------------------------------------------
            ! Remove ponds on thin ice
            !--------------------------------------------------------------
            apondn = c0
            hpondn = c0
            vpondn = c0
            hlid = c0

         else

            !-----------------------------------------------------------
            ! initialize pond area as fraction of ice
            !-----------------------------------------------------------
            apondn = apnd*alvl_tmp

            !-----------------------------------------------------------
            ! update pond volume
            !-----------------------------------------------------------
            ! add melt water
            ! note melt from deformed area is routed onto level area
            if (use_smliq_pnd) then
               dvpondn = rfrac/rhofresh*(meltt*rhoi + meltsliqn)
            else
               dvpondn = rfrac/rhofresh*(meltt*rhoi + melts*rhos + frain*dt)
            endif
            ! Track lost meltwater dvn is volume of meltwater (m3/m2) captured
            ! over entire grid cell area. Multiply by (1-rfrac)/rfrac to get
            ! loss over entire area. And divide by aicen to get loss per unit
            ! category area (for consistency with melttn, frpndn, etc)
            rfpndn = dvpondn * (c1-rfrac) / rfrac
            dvn_temp = dvpondn

            ! shrink pond volume under freezing conditions
            if (trim(frzpnd) == 'cesm') then
               Tp = Timelt - Td
               dTs = max(Tp - Tsfcn,c0)
               dvpondn = dvpondn - vpondn * (c1 - exp(rexp*dTs/Tp))

            else
               ! trim(frzpnd) == 'hlid' Stefan approximation
               ! assumes pond is fresh (freezing temperature = 0 C)
               ! and ice grows from existing pond ice
               hlid = ipnd
               if (dvpondn == c0) then ! freeze pond
                  Ts = Tair - Tffresh
                  if (Ts < c0) then
                     ! if (Ts < -c2) then ! as in meltpond_cesm
                     bdt = -c2*Ts*kice*dt/(rhoi*Lfresh)
                     dhlid = p5*sqrt(bdt)                  ! open water freezing
                     if (hlid > dhlid) dhlid = p5*bdt/hlid ! existing ice
                     dhlid = min(dhlid, hpnd*rhofresh/rhoi)
                     hlid = hlid + dhlid
                  else
                     dhlid = c0 ! to account for surface inversions
                  endif
               else ! convert refrozen pond ice back to water
                  dhlid = max(fsurfn*dt / (rhoi*Lfresh), c0) ! > 0
                  dhlid = -min(dhlid, hlid) ! < 0
                  hlid = max(hlid + dhlid, c0)
                  if (hs - dhs < puny) then ! pond ice is snow-free
                     ffrac = c1 ! fraction of fsurfn over pond used to melt ipond
                     if (fsurfn > puny) &
                          ffrac = min(-dhlid*rhoi*Lfresh/(dt*fsurfn), c1)
                  endif
               endif
               dvpondn = dvpondn - dhlid*apondn*rhoi/rhofresh
            endif

            ! Track lost/gained meltwater per unit category area from pond 
            ! lid freezing/melting. Note sign flip relative to dvn convention
            ilpndn = dvn_temp - dvpondn
            
            !-----------------------------------------------------------
            ! update pond area and depth
            !-----------------------------------------------------------
            call pond_hypsometry(hpond=hpondn, apond=apondn, vpond=vpondn, &
                                 dvpond=dvpondn, aicen=aicen, alvl=alvl_tmp)

            ! limit pond depth to maintain nonnegative freeboard
            dhpond = min(((rhow-rhoi)*hi - rhos*hs)/rhofresh - hpondn, c0)
            call pond_hypsometry(hpond=hpondn, apond=apondn, dhpond=dhpond, &
                                 alvl=alvl_tmp)
            ! at this point apondn is the fraction of the entire category 
            ! (level + deformed) with ponds on it
            frpndn = - dhpond * apondn
            
            vpondn = hpondn*apondn
            ! note, this implies that if ponds fully drain or freeze their
            ! depressions cease to exist and the lid ice also ceases to exist
            if (vpondn <= c0) then
               vpondn = c0
               apondn = c0
               hpondn = c0
               hlid = c0
            endif

            !-----------------------------------------------------------
            ! drainage due to permeability (flushing)
            ! setting dpscale = 0 turns this off
            ! NOTE this uses the initial salinity and melting T profiles
            !-----------------------------------------------------------

            if (ktherm /= 2 .and. hpondn > c0 .and. dpscale > puny) then
               draft = (rhos*hs + rhoi*hi)/rhow + hpondn
               deltah = hpondn + hi - draft
               pressure_head = gravit * rhow * max(deltah, c0)
               Tmlt(:) = -sicen(:) * depressT
               call brine_permeability(qicen, &
                    sicen, Tmlt, perm)
               if (icepack_warnings_aborted(subname)) return
               drain = perm*pressure_head*dt / (viscosity_dyn*hi) * dpscale
               deltah = min(drain, hpondn)
               dvpondn = -deltah*apondn
               vpondn = vpondn + dvpondn
               apondn = max(c0, min(apondn &
                    + 0.5*dvpondn/(pndaspect*apondn), alvl_tmp))
               hpondn = c0
               if (apondn > puny) hpondn = vpondn/apondn
            endif

         endif

         !-----------------------------------------------------------
         ! Reload tracer array
         !-----------------------------------------------------------

         hpnd = hpondn
         apnd = apondn / alvl_tmp
         if (trim(frzpnd) == 'hlid') ipnd = hlid

      endif

      end subroutine compute_ponds_lvl

!=======================================================================

! determine the liquid fraction of brine in the ice and the permeability

      subroutine brine_permeability(qicen, salin, Tmlt, perm)

      use icepack_therm_shared, only: calculate_Tin_from_qin

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         qicen, &  ! enthalpy for each ice layer (J m-3)
         salin, &  ! salinity (ppt)
         Tmlt      ! melting temperature (C)

      real (kind=dbl_kind), intent(out) :: &
         perm      ! permeability (m^2)

      ! local variables

      real (kind=dbl_kind) ::   &
         Sbr       ! brine salinity

      real (kind=dbl_kind), dimension(nilyr) ::   &
         Tin, &    ! ice temperature (C)
         phi       ! liquid fraction

      integer (kind=int_kind) :: k

      character(len=*),parameter :: subname='(brine_permeability)'

      !-----------------------------------------------------------------
      ! Compute ice temperatures from enthalpies using quadratic formula
      !-----------------------------------------------------------------

      do k = 1,nilyr
         Tin(k) = calculate_Tin_from_qin(qicen(k),Tmlt(k))
      enddo

      !-----------------------------------------------------------------
      ! brine salinity and liquid fraction
      !-----------------------------------------------------------------

      do k = 1,nilyr
         Sbr = c1/(1.e-3_dbl_kind - depressT/Tin(k)) ! Notz thesis eq 3.6
         phi(k) = salin(k)/Sbr ! liquid fraction
         if (phi(k) < 0.05) phi(k) = c0 ! impermeable
      enddo

      !-----------------------------------------------------------------
      ! permeability
      !-----------------------------------------------------------------

      perm = 3.0e-8_dbl_kind * (minval(phi))**3

      end subroutine brine_permeability

!=======================================================================

! compute the changes in pond area and depth

      subroutine pond_hypsometry(hpond, apond, vpond, dhpond, dvpond, aicen, alvl)

      real (kind=dbl_kind), intent(inout) :: &
         hpond     ! pond depth tracer
      
      real (kind=dbl_kind), intent(inout), optional :: &
         apond, &  ! pond fractional area of category (apnd*alvl for lvl ponds)
         vpond     ! pond volume per unit category area (apnd*alvl*hpnd)  

      real (kind=dbl_kind), intent(in), optional :: &
         dhpond, & ! incoming change in pond depth (may be converted to dv)
         dvpond, & ! incoming change in pond volume per unit category area
         aicen, &  ! category fractional area
         alvl      ! category fraction level ice
      
      ! local variables
      
      real (kind=dbl_kind) :: &
         dv, &     ! local variable for change in pond volume
         vp        ! local variable for pond volume per unit category area
      
      character(len=*),parameter :: subname='(pond_hypsometry)'

      ! Behavior is undefined if dhpond and dvpond are both present
      if (present(dhpond) .and. present(dvpond)) then
         call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         call icepack_warnings_add(subname//" dhpond and dvpond cannot both be input" )
         return
      endif
      ! or both absent
      if ((.not. present(dhpond)) .and. (.not. present(dvpond))) then
         call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
         call icepack_warnings_add(subname//" dhpond and dvpond cannot both be absent" )
         return
      endif

      if (trim(pndhyps) == 'none') then
         if (present(dhpond)) then ! simply change hpond by dhpond
            hpond = hpond + dhpond
         else ! dvpond must be present (due to input checking above)
            ! Check that apond is present
            if (.not. present(apond)) then
               call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
               call icepack_warnings_add(subname//" apond must be present if we are modifying apond" )
               return
            endif
            if ((.not. present(aicen)) .or. (.not. present(alvl))) then
               call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
               call icepack_warnings_add(subname//" missing aicen and/or alvl")
               return
            endif
            vpond = vpond + dvpond
            if (vpond <= c0) then
               vpond = c0
               apond = c0
            endif

            if (apond*aicen > puny) then ! existing ponds
               apond = max(c0, min(alvl, &
                    apond + 0.5*dvpond/(pndaspect*apond)))
               hpond = c0
               if (apond > puny) &
                  hpond = vpond/apond

            elseif (alvl*aicen > c10*puny) then ! new ponds
               apond = min (sqrt(vpond/pndaspect), alvl)
               hpond = pndaspect * apond ! Possible loss of meltwater if apondn == alvl_tmp

            else           ! melt water runs off deformed ice
               apond = c0
               hpond = c0 ! Loss of meltwater for very deformed ice
            endif
            apond = max(apond, c0)
         endif
      elseif (trim(pndhyps) == 'fixed') then
         if (.not. present(apond)) then
            call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            call icepack_warnings_add(subname//" apond must be present if we are modifying apond" )
            return
         endif
         if (.not. present(alvl)) then
            call icepack_warnings_setabort(.true.,__FILE__,__LINE__)
            call icepack_warnings_add(subname//" missing alvl")
            return
         endif
         ! Get the change in volume
         if (present(dvpond)) then
            dv = dvpond
         else ! compute change in volume from change in depth
            dv = dhpond * apond
         endif
         if (present(vpond)) then
            vp = vpond
            vpond = max(vpond + dv, c0)
         else ! compute initial pond volume from area and depth
            vp = apond * hpond
         endif
         vp = vp + dv ! update pond volume
         ! Compute pond area assuming that apond*pndaspect = hpond
         if (vp <= puny) then
            apond = c0
            hpond = c0
         else
            apond = min(sqrt(vp/pndaspect), alvl)
            ! preserve pond volume if pond fills all available level area
            hpond = c0
            if (apond > puny) hpond = vp/apond
         endif
      endif
      
      end subroutine pond_hypsometry

!=======================================================================

      end module icepack_meltpond_lvl

!=======================================================================
