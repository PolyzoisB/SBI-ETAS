module models_mod
  use global_params_mod
  use rng_mod, only: rng_uniform, rng_poisson
  private
  public :: etasi_sim

contains

 subroutine etasi_sim(lt_bg, ln_bg, pm_set, t_sim, lat_sim, lon_sim, mag_sim, nsim)

   implicit none

   real(8), intent(in) :: lt_bg(:), ln_bg(:)
   
   real(8), intent(in) :: pm_set(:)
   real(8), allocatable, intent(out) :: t_sim(:), lat_sim(:), lon_sim(:), mag_sim(:)
   integer, intent(out) :: nsim
      
   real(8) :: pp, cc, alpha, K, dd0, gamma, qdec, tau1, dr1, bval, stdv, bg_rate
   real(8), allocatable :: t_all(:), lat_all(:), lon_all(:), mag_all(:)
   real(8), allocatable :: t_tg(:), lat_tg(:), lon_tg(:), mag_tg(:)
   integer, allocatable :: keep_idx(:)
   integer :: nevent, n_tg, nbg_tot, nbg
   integer :: i, j, ii, jj, jk, nn
   real(8) :: rr, q0, qq, qc, qsup
   real(8) :: x, ll, deltar, theta
   real(8) :: muaft, tt0, xx0, yy0, time
   real(8) :: dt, dr, prob, ln, lt
   logical :: keep_event
   
   pp = pm_set(1) ! 1.104D0
   cc = pm_set(2) * t_day_to_sec ! 0.04D0 * t_day_to_sec
   alpha = pm_set(3) ! 2.418D0
   K = pm_set(4) ! 0.074D0
   dd0 = pm_set(5) ! 9.35D-6
   gamma = pm_set(6) ! 1.452D0
   qdec = pm_set(7) ! 1.715D0
   tau1 = pm_set(8) ! 58.0D0
   dr1 = pm_set(9) ! 50.0D0
   bg_rate = pm_set(10) 
   bval = pm_set(11) ! 1.06D0
   stdv = pm_set(12)

   q0 = m0
   qc = mc
   qsup = msup

   allocate(t_all(max_events), lat_all(max_events), lon_all(max_events), mag_all(max_events))
   allocate(t_tg(max_events), lat_tg(max_events), lon_tg(max_events), mag_tg(max_events))
  
   nevent = 0
   n_tg = 0
   
   ! Generate background events
   nbg = size(lt_bg)
   nbg_tot = int(bg_rate*tlast*(lat_max_0-lat_min_0)*(lon_max_0-lon_min_0)) ! Number of expected main shocks

   do while (nevent < nbg_tot .and. nevent < max_events)
     rr=rng_uniform()
     qq = q0 - log10(1.0D0 - rr) / bval

     jk = 1 + int(rng_uniform() * nbg)
     if(jk.lt.1) jk=1
     if(jk.gt.nbg) jk=nbg

     lt = lt_bg(jk)+2*0.01*rng_uniform()-0.01
     ln = ln_bg(jk)+2*0.01*rng_uniform()-0.01
     
     if (lt <= lat_min_0 .or. lt >= lat_max_0) cycle
     if (ln <= lon_min_0 .or. ln >= lon_max_0) cycle
     
     nevent = nevent + 1
     lat_all(nevent) = lt
     lon_all(nevent) = ln
     mag_all(nevent) = min(qq,qsup)
     t_all(nevent) = rng_uniform() * tlast
     
     if (t_all(nevent) >= tc.and.mag_all(nevent) >= qc)then
      if(lat_all(nevent) > lat_min.and.lat_all(nevent) < lat_max)then
       if(lon_all(nevent) > lon_min.and.lon_all(nevent) < lon_max)then
        n_tg = n_tg + 1               ! Number of mothers in target region
        t_tg(n_tg) = t_all(nevent)
        lat_tg(n_tg) = lat_all(nevent)
        lon_tg(n_tg) = lon_all(nevent)
        mag_tg(n_tg) = mag_all(nevent)
       endif
      endif
     endif

   enddo  ! Loop on main

   !!!! Generate aftershocks !!!!
   j = 1

   do while (j <= nevent .and. nevent < max_events)
   
     muaft = K* exp(alpha * (mag_all(j) - q0))
     nn = rng_poisson(muaft)
     
     tt0 = t_all(j)
     xx0 = lat_all(j)
     yy0 = lon_all(j)
     
     do i = 1 , nn
      time = cc * (rng_uniform()**(1.D0 / (1.0D0-pp))) - cc
      
      if (tt0 + time > tlast) then
       goto 77
      endif
      
      ll= dd0 * exp(gamma * (mag_all(j) - qc))
      x = rng_uniform()
      deltar = sqrt(ll * (x**(1.0D0 / (1.0D0-qdec))) - ll)!/100.
      theta = 2.0D0 * rng_uniform() * 3.14159265358979D0
      
      if (xx0 + deltar * sin(theta) > lat_min_0 .and. xx0 + deltar * sin(theta) < lat_max_0) then
       if (yy0 + deltar * cos(theta) > lon_min_0 .and. yy0 + deltar * cos(theta) < lon_max_0) then
        
        if (nevent >= max_events) exit
        
        nevent = nevent + 1    
        rr = rng_uniform()
        qq = q0 - log10(1.0D0 - rr) / bval
        
        mag_all(nevent) = min(qq, qsup)
        t_all(nevent) = tt0 + time
        lat_all(nevent) = xx0 + deltar * sin(theta)
        lon_all(nevent) = yy0 + deltar * cos(theta)
        
        if(t_all(nevent) >= tc.and.qq >= qc)then
         if(lat_all(nevent) > lat_min.and.lat_all(nevent) < lat_max)then
          if(lon_all(nevent) > lon_min.and.lon_all(nevent) < lon_max)then
            n_tg = n_tg + 1  ! Number of aftershocks in target region
            t_tg(n_tg) = t_all(nevent)
            lat_tg(n_tg) = lat_all(nevent)
            lon_tg(n_tg) = lon_all(nevent)
            mag_tg(n_tg) = mag_all(nevent)
          endif
         endif
        endif
       
       endif
      endif
    
      77 continue
     enddo
     
     j = j + 1
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
       if (rng_uniform() > prob) then
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
