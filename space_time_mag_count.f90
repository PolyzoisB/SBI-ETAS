! summary_stats_mod.f90
module space_time_mag_count_mod
  use global_params_mod, only: ncmagn, nctime, ncspace, ncmain, nctimef, ncspacef, &
                               ncmagnf, step_mgn, step_mgnf, thtime, thml, thtimef, &
                               thspacef, thmlf, mmain
  implicit none
  contains

  subroutine count_fore_aft(time, lat, lon, mag, nflag, n, &
                            thspace_in, nfore_norm, naft_norm, nmain, flag_main)
    integer, intent(in)    :: n
    real(8), intent(in)    :: time(n), lat(n), lon(n), mag(n)
    integer, intent(inout) :: nflag(n), flag_main(ncmain)
    real(8), intent(in)    :: thspace_in(ncspace)
    real(8), intent(out)   :: nfore_norm(ncmagnf,nctimef,ncspacef,ncmain)
    real(8), intent(out)   :: naft_norm(ncmagn,nctime,ncspace,ncmain)
    integer, intent(out)   :: nmain(ncmain)

    integer :: nfore(ncmagnf,nctimef,ncspacef,ncmain), naft(ncmagn,nctime,ncspace,ncmain)
    real(8) :: thspace_loc(ncspace), ddt, ddr, tmax, rmax, mmain_cut
    integer :: i, j, iqmax, indmain(ncmain) 

    nmain = 0;  nfore = 0;  naft = 0;  indmain = 0

    mmain_cut = mmain
    do i = 1, n
      if (nflag(i) /= 0 .or. mag(i) < mmain_cut) cycle
      if (int(mag(i)) == 4) then;      iqmax = 1
      else if (int(mag(i)) == 5) then;  iqmax = 2
      else if (int(mag(i)) >= 6) then;  iqmax = 3
      else; cycle
      end if
      nmain(iqmax) = nmain(iqmax) + 1

      !! Aftershocks (forward in time) !!
      thspace_loc = thspace_in
      thspace_loc(ncspace) = 0.01D0 * (10.0D0**(0.5D0 * mag(i)))
      rmax = maxval(thspace_loc)
      tmax = maxval(thtime)
      do j = i + 1, n
        if (nflag(j) == 0 .or. nflag(j) == 2) cycle
        ddt = (time(j) - time(i)) 
        if (ddt > tmax) exit
        call distkm_sm(lat(j), lon(j), lat(i), lon(i), ddr)
        if (ddr > rmax) cycle
        nflag(j) = 2
        call bin_event(mag(j), ddt, ddr, thtime, thspace_loc, thml,&
                        nctime, ncspace, ncmagn, step_mgn, naft(:,:,:,iqmax))
      end do

      !! Foreshocks (backward in time) !!
      rmax = maxval(thspacef)
      tmax = maxval(thtimef)
      do j = i - 1, 1, -1
        if (nflag(j) == 0 .or. nflag(j) == 2) cycle
        ddt = (time(i) - time(j)) 
        if (ddt > tmax) exit
        call distkm_sm(lat(j), lon(j), lat(i), lon(i), ddr)
        if (ddr > rmax) cycle
        nflag(j) = 2
        call bin_event(mag(j), ddt, ddr, thtimef, thspacef, thmlf, &
                        nctimef, ncspacef, ncmagnf, step_mgnf, nfore(:,:,:,iqmax))
      end do
    end do

    ! Normalize
    do iqmax = 1, ncmain
      if (nmain(iqmax) == 0) cycle
      indmain(iqmax) = 1
      nfore_norm(:,:,:,iqmax) = real(nfore(:,:,:,iqmax), 8) / real(nmain(iqmax), 8)
      naft_norm(:,:,:,iqmax)  = real(naft(:,:,:,iqmax), 8)  / real(nmain(iqmax), 8)
    end do
    flag_main = indmain
  end subroutine

  subroutine bin_event(m, dt, dr, tht, ths, thm, nt, ns, nm, step, counts)
    real(8), intent(in)    :: m, dt, dr, tht(nt), ths(ns), thm(nm), step
    integer, intent(in)    :: nt, ns, nm
    integer, intent(inout) :: counts(nm, nt, ns)
    integer :: it, is, im
    do im = 1, nm
      if (m < thm(im) .or. m >= thm(im) + step) cycle
      do it = 1, nt
        if (dt > tht(it)) cycle
        do is = 1, ns
          if (dr > ths(is)) cycle
          counts(im, it, is) = counts(im, it, is) + 1
        end do
      end do
    end do
  end subroutine

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

end module space_time_mag_count_mod
