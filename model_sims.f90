module models_mod
  use global_params_mod
  private
  public :: etasi_sim

contains

 subroutine etasi_sim(lt_bg, ln_bg, nbg, pm_set, t_sim, lat_sim, lon_sim, mag_sim, nsim, seed)

   implicit none

   real(8), intent(in) :: lt_bg(:), ln_bg(:)
   integer, intent(in) :: nbg

   real(8), intent(in) :: pm_set(:)
   real(8), allocatable, intent(out) :: t_sim(:), lat_sim(:), lon_sim(:), mag_sim(:)
   integer, intent(out) :: nsim
   integer, intent(inout) :: seed
   
   real(8) :: pp, cc, alpha, K, dd0, gamma, qdec, tau1, dr1, bval, stdv, bg_rate
   real(8), allocatable :: t_all(:), lat_all(:), lon_all(:), mag_all(:)
   real(8), allocatable :: t_tg(:), lat_tg(:), lon_tg(:), mag_tg(:)
   integer, allocatable :: keep_idx(:)
   integer :: nevent, n_tg, nbg_tot
   integer :: i, j, ii, jj, jk, nn
   real(8) :: rr, qq, qmin, qsup
   real(8) :: x, ll, deltar, theta
   real(8) :: muaft, tt0, xx0, yy0, time
   real(8) :: dt, dr, prob, ln, lt
   logical :: keep_event
   
   pp = pm_set(1) ! 1.104D0
   cc = pm_set(2) * t_day_to_sec ! 0.04D0 * t_day_to_sec
   alpha = pm_set(3) ! 2.418D0
   K = pm_set(4) ! 0.074D0
   dd0 = 10**pm_set(5) ! 9.35D-6
   gamma = pm_set(6) ! 1.452D0
   qdec = pm_set(7) ! 1.715D0
   tau1 = pm_set(8) ! 58.0D0
   dr1 = pm_set(9) ! 50.0D0
   bg_rate = pm_set(10) 
   bval = pm_set(11) ! 1.06D0
   
   stdv = 0.4D0

   qmin = mc
   qsup = msup

   allocate(t_all(max_events), lat_all(max_events), lon_all(max_events), mag_all(max_events))
   allocate(t_tg(max_events), lat_tg(max_events), lon_tg(max_events), mag_tg(max_events))

   do i = 1 , max_events
    mag_all(i) = -100
   enddo
   
   nevent = 0
   n_tg = 0
   ! Generate background events
   nbg_tot = int(bg_rate*tlast*(lat_max-lat_min)*(lon_max-lon_min))  ! Number of expected main shocks
   do while (nevent < nbg_tot .and. nevent < max_events)
     rr=ran2(seed)
     qq = qmin - log10(1.0D0 - rr) / bval

     jk = 1 + int(ran2(seed) * nbg)
     if(jk.lt.1) jk=1
     if(jk.gt.nbg) jk=nbg
     lt = lt_bg(jk)+2*0.01*ran2(seed)-0.01
     ln = ln_bg(jk)+2*0.01*ran2(seed)-0.01
     if (lt <= lat_min .or. lt >= lat_max) cycle
     if (ln <= lon_min .or. ln >= lon_max) cycle
     nevent = nevent + 1
     lat_all(nevent) = lt
     lon_all(nevent) = ln
     mag_all(nevent) = min(qq,qsup)
     t_all(nevent) = ran2(seed) * tlast
     
     if(t_all(nevent).ge.tc)then
       n_tg = n_tg + 1  ! Number of mothers in target region
       t_tg(n_tg) = t_all(nevent)
       lat_tg(n_tg) = lat_all(nevent)
       lon_tg(n_tg) = lon_all(nevent)
       mag_tg(n_tg) = mag_all(nevent)
     endif
   enddo  ! Loop on main

   !!!! Generate aftershocks !!!!
   do j = 1, max_events
     if(mag_all(j) < 0 .or. nevent >= max_events) exit

     muaft = K* exp(alpha * (mag_all(j) - qmin))
     nn = zbqlpoi(muaft, seed)
     tt0 = t_all(j)
     xx0 = lat_all(j)
     yy0 = lon_all(j)
     
     do i = 1 , nn
      time = cc * (ran2(seed)**(1.D0 / (1.0D0-pp))) - cc
      if (tt0 + time > tlast) then
       goto 77
      endif
      
      ll= dd0 * exp(gamma * (mag_all(j) - qmin))
      x = ran2(seed)
      deltar = sqrt(ll * (x**(1.0D0 / (1.0D0-qdec))) - ll)!/100.
      theta = 2.0D0 * ran2(seed) * 3.14159265358979D0
      
      if (xx0 + deltar * sin(theta) > lat_min .and. xx0 + deltar * sin(theta) < lat_max) then
       if (yy0 + deltar * cos(theta) > lon_min .and. yy0 + deltar * cos(theta) < lon_max) then
        
        nevent = nevent + 1    
        rr = ran2(seed)
        qq = qmin - log10(1.0D0 - rr) / bval
        mag_all(nevent) = min(qq, qsup)
        t_all(nevent) = tt0 + time
        lat_all(nevent) = xx0 + deltar * sin(theta)
        lon_all(nevent) = yy0 + deltar * cos(theta)
        
        if(t_all(nevent) >= tc)then
          n_tg = n_tg + 1  ! Number of aftershocks in target region
          t_tg(n_tg) = t_all(nevent)
          lat_tg(n_tg) = lat_all(nevent)
          lon_tg(n_tg) = lon_all(nevent)
          mag_tg(n_tg) = mag_all(nevent)
        endif
       
       endif
      endif
    
      77 continue
    
     enddo

     if (nevent >= max_events) exit       ! hard safety break    
    enddo ! end aftershock generation
     
    !! Subroutine that sorts events based in t_sim array
    call hpsort4(n_tg, t_tg, lon_tg, lat_tg, mag_tg)
    
    allocate(keep_idx(n_tg))
    keep_idx(1) = 1
    nsim = 1
    do jj = 2 , n_tg
     keep_event = .true.
     do ii = jj-1 , 1, -1
       dt = t_tg(jj) - t_tg(ii)
       if(dt > tau1) exit

       call distkm_sm(lat_tg(jj), lon_tg(jj), lat_tg(ii), lon_tg(ii), dr)
       if(dr > dr1) cycle
       
       prob = 0.5D0 + 0.5D0 * erf((mag_tg(jj) - mag_tg(ii)) / stdv)
       if (ran2(seed) > prob) then
         keep_event = .false.
         exit
       endif
     enddo

     if (keep_event) then
       nsim = nsim + 1
       keep_idx(nsim) = jj
     endif
   enddo
   
   allocate(t_sim(nsim), lat_sim(nsim), lon_sim(nsim), mag_sim(nsim))
   do i = 1, nsim
     jj = keep_idx(i)
     t_sim(i) = t_tg(jj)
     lat_sim(i)  = lat_tg(jj)
     lon_sim(i)  = lon_tg(jj)
     mag_sim(i)  = mag_tg(jj)
   enddo

   deallocate(t_all, lat_all, lon_all, mag_all, t_tg, lat_tg, lon_tg, mag_tg, keep_idx)

  end subroutine etasi_sim

  ! ************************************************************************
  !    FULL ETAS-BC
  !    Background -> connector chain -> multi-generation aftershocks
  !
  !    Definitions:
  !    - Background         : initial event of the sequence
  !    - Connectors         : muBC chain
  !    - Mainshock          : last connector of the sequence
  !    - Aftershocks        : ETAS multi-generational cascade
  !
  !    Final output flags:
  !      0 = background
  !      1 = pre-mainshock connectors
  !      2 = aftershocks with t < tM
  !      3 = mainshock
  !      4 = aftershocks with t > tM
  !
  !    Output:
  !      fort.89  : t, m, x, y, kflag
  !      fort.191 : n_fore_tot, n_conn_pre, n_af_pre, mM, iseq
  !
  !    NOTES:
  !    - final catalog sorted temporally
  !    - here I also keep post-mainshock aftershocks
  !************************************************************************
  subroutine etasbc_sim(lt_bg, ln_bg, nbg, pm_set, time_sim, lat_sim, lon_sim, mag_sim, nsim, seed)
    real(8), intent(in) :: lt_bg(nbg),ln_bg(nbg), pm_set(:)
    integer, intent(in) :: nbg
    integer, intent(inout) :: seed
    real(8), allocatable, intent(out) :: time_sim(:), lat_sim(:), lon_sim(:), mag_sim(:)
    integer, intent(out) :: nsim
  
    !-------------------- global catalog -----------------------
    real(8) :: itime(max_events), q(max_events), xx(max_events), yy(max_events)
    integer :: kflag(max_events), seqid(max_events)
   
    !-------------------- connector parameters -----------------
    real(8) :: phi0_conn, pp_conn, cc_conn
    real(8) :: DD_conn, gamma_conn, qdec_conn, bb_conn
    real(8) :: qmin_conn, qsup, bg_rate

    !-------------------- aftershock parameters ----------------
    real(8) :: AAAf, alpha_af, pp_af, cc_af
    real(8) :: DD_af, gamma_af, qdec_af, bb_af, qmin_af

    !-------------------- variables --------------
    integer :: nevent, nevc, nmain
    integer :: i, jk, istart, iend, i_main, ip, n_af, iaf
    integer :: iseq, j_conn, nn, kconn
    
    real(8) :: final_time
    real(8) :: time, tt0, xx0d, yy0d, qq0d
    real(8) :: theta, deltar, rr, x, ll, r2
    real(8) :: xxxx, yyyy, mu_af, u_af, dt_af, tchild_af
    real(8) :: tmain, mmain, pi
    parameter (pi=3.141592653589793d0)

    real(8) :: dt, dr, prob, tau1, dr1, stdv
    integer :: ii, jj
    logical :: keep_event
    integer, allocatable :: keep_idx(:)
    integer :: n_tg
    real(8), allocatable :: t_tg(:), lat_tg(:), lon_tg(:), mag_tg(:)

    !***********************************************************

    nmain = int(bg_rate*tlast*(lat_max-lat_min)*(lon_max-lon_min))  ! Number of expected main shocks

    !--- connector/background magnitudes
    qmin_conn = mc
    qsup      = msup
    bb_conn   = pm_set(12) ! /dlog(10.d0)
    bg_rate   = pm_set(10)

    !--- connectors
    phi0_conn = pm_set(18) ! 0.018d0
    pp_conn   = pm_set(13) ! 1.22d0
    cc_conn   = pm_set(14) ! 0.01d0 * t_day_to_sec
    DD_conn   = 10**pm_set(15) ! 0.018d0
    gamma_conn= pm_set(16) ! 1.440d0
    qdec_conn = pm_set(17) ! 1.550d0

    !--- aftershocks:
    qmin_af   = mc
    pp_af     = pm_set(1) !1.104d0 ! 1.2d0
    alpha_af  = pm_set(3) ! 1.05d0  ! 0.97d0
    AAAf      = pm_set(4) ! 0.074d0 ! 0.10d0
    cc_af     = pm_set(2) ! 0.04d0 * t_day_to_sec ! 0.01d0 * t_day_to_sec
    DD_af     = 10**pm_set(5) ! 9.35D-6 * 1E4 ! 0.006d0
    gamma_af  = pm_set(6) ! 1.452d0 ! 0.85d0
    qdec_af   = pm_set(7) ! 1.715d0 ! 1.506d0
    bb_af     = pm_set(11) ! 1.06d0 ! /dlog(10.d0)

    !--- aftershock incompleteness:
    tau1 = pm_set(8) ! 58.0d0
    dr1  = pm_set(9) ! 50.0d0
    stdv = sqrt(2.0d0)*0.28284d0 ! 0.4d0

    final_time = tlast

    nevent    = 0
    nevc      = nbg

    !======================================================================
    !  1) GENERATE THE SEQUENCES
    !======================================================================
    do iseq=1,nmain

      !-------- SEQUENCE BACKGROUND
      istart = nevent + 1
      nevent = nevent + 1
      if(nevent.gt.max_events) goto 888

      time = ran2(seed)*final_time

      jk = 1 + int(ran2(seed)*nevc)
      if(jk.lt.1) jk=1
      if(jk.gt.nevc) jk=nevc
      xx(nevent) = lt_bg(jk)
      yy(nevent) = ln_bg(jk)

      rr = ran2(seed)
      qq0d = qmin_conn - (1.d0/bb_af)*dlog10(1.d0-rr)
      if(qq0d.gt.qsup) qq0d = qsup

      itime(nevent) = time
      q(nevent)     = qq0d
      kflag(nevent) = 0
      seqid(nevent) = iseq
      if(time.ge.tc) then
        n_tg = n_tg + 1  ! Number of mothers in target region
        t_tg(n_tg) = time
        lat_tg(n_tg) = xx(nevent)
        lon_tg(n_tg) = yy(nevent)
        mag_tg(n_tg) = qq0d
      endif
      !======================================================================
      !  2a) CONNECTOR CHAIN
      !======================================================================
      j_conn = istart
      do kconn=1,1000

        qq0d = dble(q(j_conn))
        tt0  = itime(j_conn)
        xx0d = dble(xx(j_conn))
        yy0d = dble(yy(j_conn))

        !----------- decide whether the chain continues
        nn=0    
        if(phi0_conn.ge.ran2(seed)) nn=1
        nn=min(1,nn)
        if(nn.eq.0) goto 205

        !----------- connector time
        time = cc_conn*((ran2(seed))**(1.d0/(1.d0-pp_conn))) &
                       - cc_conn
        if(tt0+time.gt.final_time) goto 205

        !----------- connector space
        ll = DD_conn*dexp(gamma_conn*(qq0d-qmin_conn))
        x  = ran2(seed)
        r2 = ll*(x**(1.d0/(1.d0-qdec_conn))) - ll
        if(r2.lt.0.d0) r2 = 0.d0
        deltar = dsqrt(r2)/100.d0

        theta = ran2(seed)*2.d0*pi
        xxxx  = xx0d + deltar*dsin(theta)
        yyyy  = yy0d + deltar*dcos(theta)

        !----------- connector magnitude: child >= parent
        rr   = ran2(seed)
        qq0d = qq0d - (1.d0/bb_conn)*dlog10(1.d0-rr)
        if(qq0d.gt.qsup) qq0d = qsup

        !----------- add connector
        nevent = nevent + 1
        if(nevent.gt.max_events) goto 888

        itime(nevent) = tt0 + time
        q(nevent)     = qq0d
        xx(nevent)    = xxxx
        yy(nevent)    = yyyy
        kflag(nevent) = 1
        seqid(nevent) = iseq

        j_conn    = nevent

        if(itime(nevent).ge.tc) then
          n_tg = n_tg + 1  ! Number of connectors in target region
          t_tg(n_tg) = itime(nevent)
          lat_tg(n_tg) = xx(nevent)
          lon_tg(n_tg) = yy(nevent)
          mag_tg(n_tg) = q(nevent)
        endif

      enddo

      205 continue

      !-------- last connector = mainshock
      i_main = j_conn
      tmain  = itime(i_main)
      mmain  = dble(q(i_main))

      !======================================================================
      !  2b) MULTI-GENERATION AFTERSHOCKS
      !======================================================================
      iend = nevent
      ip   = istart

      700     continue
      do while (ip.le.iend)

        mu_af = AAAf * exp((alpha_af*(dble(q(ip))-qmin_af)))
        if(mu_af.le.0.d0) then
          ip = ip + 1
          goto 700
        endif

        n_af = ZBQLPOI(mu_af, seed)
        if(n_af.le.0) then
          ip = ip + 1
          goto 700
        endif

        tt0  = itime(ip)
        xx0d = dble(xx(ip))
        yy0d = dble(yy(ip))
        qq0d = dble(q(ip))

        do iaf=1,n_af

          if(nevent.gt.max_events) goto 888

          !-------------- Omori-Utsu time
          u_af  = ran2(seed)
          dt_af = cc_af * &
              ( (1.d0-u_af)**(-1.d0/(pp_af-1.d0)) - 1.d0 )

          tchild_af = tt0 + dt_af
          if(tchild_af.gt.final_time) goto 771

          !-------------- aftershock space
          ll = DD_af*dexp(gamma_af*(qq0d-qmin_af))
          x  = ran2(seed)
          r2 = ll*(x**(1.d0/(1.d0-qdec_af))) - ll
          if(r2.lt.0.d0) r2 = 0.d0
          deltar = dsqrt(r2)/100.d0
          theta = ran2(seed)*2.d0*pi
          xxxx  = xx0d + deltar*dsin(theta)
          yyyy  = yy0d + deltar*dcos(theta)

          !-------------- aftershock magnitude
          rr   = ran2(seed)
          qq0d = qmin_af - (1.d0/bb_af)*dlog10(1.d0-rr)
          if(qq0d.gt.qsup) qq0d = qsup

          !-------------- add aftershock
          nevent = nevent + 1
          if(nevent.gt.max_events) goto 888

          itime(nevent) = tchild_af
          q(nevent)     = qq0d
          xx(nevent)    = xxxx
          yy(nevent)    = yyyy
          kflag(nevent) = 2
          seqid(nevent) = iseq

          iend    = nevent

          if(itime(nevent).ge.tc) then
            n_tg = n_tg + 1  ! Number of aftershocks in target region
            t_tg(n_tg) = itime(nevent)
            lat_tg(n_tg) = xx(nevent)
            lon_tg(n_tg) = yy(nevent)
            mag_tg(n_tg) = q(nevent)
          endif

          771     continue
        enddo

        ip = ip + 1
      enddo

    enddo
    
    !! Subroutine that sorts events based in t_sim array
    call hpsort4(n_tg, t_tg, lat_tg, lon_tg, mag_tg)
    
    allocate(keep_idx(n_tg))
    keep_idx(1) = 1
    nsim = 1
    do jj = 2 , n_tg
     keep_event = .true.
     do ii = jj-1 , 1, -1
       dt = t_tg(jj) - t_tg(ii)
       if(dt > tau1) exit

       call distkm_sm(lat_tg(jj), lon_tg(jj), lat_tg(ii), lon_tg(ii), dr)
       if(dr > dr1) cycle
       
       prob = 0.5D0 + 0.5D0 * erf((mag_tg(jj) - mag_tg(ii)) / stdv)
       if (ran2(seed) > prob) then
         keep_event = .false.
         exit
       endif
     enddo

     if (keep_event) then
       nsim = nsim + 1
       keep_idx(nsim) = jj
     endif
   enddo
   
   allocate(time_sim(nsim), lat_sim(nsim), lon_sim(nsim), mag_sim(nsim))
   do i = 1, nsim
     jj = keep_idx(i)
     time_sim(i) = t_tg(jj)
     lat_sim(i)  = lat_tg(jj)
     lon_sim(i)  = lon_tg(jj)
     mag_sim(i)  = mag_tg(jj)
   enddo
 
   return
    
    888  print*,'STOP: array limit exceeded.'
    stop
  end subroutine etasbc_sim

  subroutine distkm_sm(x1, y1, x3, y3, dr)
    implicit none
    real(8), intent(in)  :: x1, y1, x3, y3
    real(8), intent(out) :: dr
    real(8), parameter :: prad = 3.14159265358979D0 / 180.0D0, dr0 = 0.01D0
    real(8) :: phi1, phi2, dphi, dlambda, a
    phi1 = x1 * prad;  phi2 = x3 * prad
    dphi = phi2 - phi1;  dlambda = (y3 - y1) * prad
    a = sin(dphi/2.0D0)**2 + cos(phi1)*cos(phi2)*sin(dlambda/2.0D0)**2
    dr = 2.0D0 * asin(sqrt(a)) * 6370.0D0
    if (dr == 0.0D0) dr = dr0
  end subroutine

  !======================================================================
  !  Poisson (Knuth)
  !======================================================================
  integer function zbqlpoi(mu, seed1) result(npoi)
    implicit none
    real(8), intent(in) :: mu
    integer, intent(inout) :: seed1
    real(8) :: L, p
    integer :: k

    if (mu <= 0.0D0) then
      npoi = 0
      return
    end if

    L = exp(-mu)
    k = 0
    p = 1.0D0

    do
      k = k + 1
      p = p * ran2(seed1)
      if (p <= L) exit
    end do

    npoi = k - 1
  end function zbqlpoi

  real(8) function ran2(idum) result(r)
    implicit none
    integer, intent(inout) :: idum
    integer, parameter :: IM1=2147483563, IM2=2147483399
    integer, parameter :: IMM1=IM1-1, IA1=40014, IA2=40692
    integer, parameter :: IQ1=53668, IQ2=52774, IR1=12211, IR2=3791
    integer, parameter :: NTAB=32, NDIV=1+IMM1/NTAB
    real(8), parameter :: AM=1.0D0/IM1, EPS=1.2D-7, RNMX=1.0D0-EPS
    integer, save :: idum2 = 123456789
    integer, save :: iv(NTAB) = 0, iy = 0
    integer :: j, k

    if (idum <= 0) then
      idum = max(-idum, 1)
      idum2 = idum
      do j = NTAB + 8, 1, -1
        k = idum / IQ1
        idum = IA1 * (idum - k * IQ1) - k * IR1
        if (idum < 0) idum = idum + IM1
        if (j <= NTAB) iv(j) = idum
      end do
      iy = iv(1)
    end if

    k = idum / IQ1
    idum = IA1 * (idum - k * IQ1) - k * IR1
    if (idum < 0) idum = idum + IM1

    k = idum2 / IQ2
    idum2 = IA2 * (idum2 - k * IQ2) - k * IR2
    if (idum2 < 0) idum2 = idum2 + IM2

    j = 1 + iy / NDIV
    iy = iv(j) - idum2
    iv(j) = idum
    if (iy < 1) iy = iy + IMM1

    r = min(AM * real(iy,8), RNMX)
  end function ran2

  subroutine hpsort4(n, time, lat, lon, mag)
   implicit none
   integer, intent(in) :: n
   real(8), intent(inout) :: time(n), lat(n), lon(n), mag(n)

   integer :: ir, l, i, j
   real(8) :: time_ra, lat_ra, lon_ra, mag_ra

   if (n <= 1) return

   l = n / 2 + 1
   ir = n

   10  continue
   if (l > 1) then
     l = l - 1
     time_ra = time(l)
     lat_ra = lat(l)
     lon_ra = lon(l)
     mag_ra = mag(l)
   else
     time_ra = time(ir)
     lat_ra = lat(ir)
     lon_ra = lon(ir)
     mag_ra = mag(ir)

     time(ir) = time(1)
     lat(ir) = lat(1)
     lon(ir) = lon(1)
     mag(ir) = mag(1)

     ir = ir - 1
     if (ir == 1) then
       time(1) = time_ra
       lat(1) = lat_ra
       lon(1) = lon_ra
       mag(1) = mag_ra
       return
     end if
   end if

    i = l
    j = l + l

    20  if (j <= ir) then
    if (j < ir) then
      if (time(j) < time(j + 1)) j = j + 1
    end if
    if (time_ra < time(j)) then
      time(i) = time(j)
      lat(i) = lat(j)
      lon(i) = lon(j)
      mag(i) = mag(j)
      i = j
      j = j + j
    else
      j = ir + 1
    end if
    goto 20
   end if

   time(i) = time_ra
   lat(i) = lat_ra
   lon(i) = lon_ra
   mag(i) = mag_ra
   goto 10
   
  end subroutine hpsort4

end module models_mod
