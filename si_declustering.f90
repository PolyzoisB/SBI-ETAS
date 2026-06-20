module susceptibility_index_mod
  use global_params_mod, only: pthmin, pthmax, bin0, max_thr, max_si_points, df,&
                               max_events, susc_index_file_name, output_cat, final_results
  implicit none
  private
  public :: si_threshold

contains

  !==================================================================
  ! Main entry point: compute NN threshold and background rate
  !==================================================================
  subroutine si_threshold(time, xxe, yye, mag, n_event, &
                          nn_thr, sim_thr, nbg, bval, lt_bg, ln_bg)
    integer, intent(in)  :: n_event
    real(8), intent(in)  :: time(n_event), mag(n_event), xxe(n_event), yye(n_event)
    real(8), intent(in)  :: bval
    real(8), intent(out) :: lt_bg(n_event), ln_bg(n_event), nn_thr, sim_thr
    integer, intent(out) :: nbg

    real(8) :: ti(n_event)
    integer :: kr(-max_thr:max_thr), nlink(-max_thr:max_thr)
    integer :: id_pr(n_event)
    integer :: ipmin, ipmax, ip_thr, kmin_shft
    real(8) :: sindex(max_si_points), si_thresh(n_event)
    integer :: k_si(max_si_points), num_points, i, id
    integer :: best_idx, nn(-max_thr:max_thr), ip_si(max_si_points)
    
    ! Relative times
    ti = time - time(1)

    ! Step 1: compute pairwise similarities → nlink, kr, id_pr
    call compute_similarities(ti, xxe, yye, mag, n_event, nlink, kr,&
                              id_pr, ipmin, ipmax, nn, bval)
    
    ! Step 2: build susceptibility index curve
    call compute_si_curve(nlink, kr, n_event, ipmin, ipmax,&
                           sindex, k_si, ip_si, si_thresh, num_points)
    
    ! Step 3: find optimal threshold via smoothed minimum
    call find_si_minimum(sindex, k_si, num_points, best_idx)
    
    if (best_idx < 1 .or. best_idx > num_points) then
       print *, 'ERROR: invalid optimal threshold index.'
       stop 1
    end if
    kmin_shft = k_si(best_idx)
    ip_thr = ip_si(best_idx)                              
    
    nbg = 0
    do i = 1, n_event
     if (i == 1) then
       id = 0
       lt_bg(i) = xxe(1);  ln_bg(i) = yye(1)
     else if (id_pr(i) >= ip_thr) then
       id = 1
     else
       id = 0
       lt_bg(nbg + 1) = xxe(i);  ln_bg(nbg + 1) = yye(i)
     end if

     if (id == 0) nbg = nbg + 1
    end do

    sim_thr = si_thresh(best_idx)
    nn_thr  = 1.0D0 / sim_thr
    ! bg_rate = real(nbg, 8) / ti(n_event)

  end subroutine

  subroutine compute_similarities(ti, xxe, yye, mag, n, &
                                  nlink, kr, id_pr, ipmin, ipmax, nn, bval)
   implicit none
   integer, intent(in)  :: n
   real(8), intent(in)  :: ti(n), xxe(n), yye(n), mag(n)
   integer, intent(out) :: nlink(-max_thr:max_thr), kr(-max_thr:max_thr)
   integer, intent(out) :: id_pr(n), ipmin, ipmax
   integer, intent(out) :: nn(-max_thr:max_thr)
   real(8), intent(in)  :: bval
   
   integer :: i, j, ip, ipf, ipmaxt
   integer :: nlink_inc(-max_thr:max_thr)  ! for deferred prefix sum
   real(8) :: sim, dt_resc, dr_resc
   ! integer :: unit_number
   real(8) :: maxdt(n), maxdr(n)

   nn = 0;  nlink = 0;  kr = 0;  nlink_inc = 0
   ipmin = int(dlog(pthmin) / bin0)
   ipmax = int(dlog(pthmax) / bin0)
   id_pr = ipmin
   ! open (unit_number,file=rescaled_distances,status="replace",action="write") 
   do i = 2, n
     ipmaxt = ipmin
     do j = i - 1, 1, -1
       call bp_similarity(xxe(i), yye(i), ti(i), xxe(j), yye(j), ti(j), &
                          mag(j), bval, sim, dt_resc, dr_resc)
       ipf = int(dlog(sim) / bin0)
       if (ipf < ipmin) ipf = ipmin
       if (ipf > ipmax) ipf = ipmax
       if(ipf.gt.ipmaxt)then
          maxdt(i)=dt_resc
          maxdr(i)=dr_resc
       endif
       ipmaxt = max(ipf, ipmaxt)
        
       ! Only increment at ipf; prefix-sum later
       nlink_inc(ipf) = nlink_inc(ipf) + 1
     end do
     ! write(unit_number,*) maxdt(i)/(365.25d0*24.d0*3600.d0), maxdr(i) ! in (years,km)
     
     nn(ipmaxt) = nn(ipmaxt)+1
     id_pr(i) = ipmaxt
     do ip = ipmaxt + 1, ipmax
       kr(ip) = kr(ip) + 1
     end do
   
   end do
   ! close(unit_number)

   ! Suffix sum: nlink(ip) = sum of nlink_inc(ip..ipmax)
   nlink(ipmax) = nlink_inc(ipmax)
   do ip = ipmax - 1, ipmin, -1
      nlink(ip) = nlink(ip + 1) + nlink_inc(ip)
   end do

  end subroutine compute_similarities

  !==================================================================
  ! BP similarity (returns 1/eta = similarity)
  !==================================================================
  subroutine bp_similarity(x1, y1, t1, x2, y2, t2, q2, bval, sim, dt_resc, dr_resc)
   implicit none
   real(8), intent(in)  :: x1, y1, t1, x2, y2, t2, q2, bval
   real(8), intent(out) :: sim, dt_resc, dr_resc

   real(8) :: dr, dt, phi1, phi2, dphi, dlambda, a
   real(8), parameter :: dt0 = 1.0D0, dr0 = 0.01D0
   real(8), parameter :: prad = 3.14159d0 / 180.0D0

   ! Haversine distance (km)
   phi1 = x1 * prad;  phi2 = x2 * prad
   dphi = phi2 - phi1;  dlambda = (y2 - y1) * prad
   a = sin(dphi/2.0D0)**2 + cos(phi1)*cos(phi2)*sin(dlambda/2.0D0)**2
   a = min(1.0D0, max(0.0D0, a))
   dr = 2.0D0 * asin(sqrt(a)) * 6370.0D0

   dt = t1 - t2
   ! Regularize independently
   if (dt == 0.0D0) dt = dt0 ! 1 sec
   if (dr == 0.0D0) dr = dr0 ! 10 mt
   ! if (dt == 0.0D0.or.dr == 0.0D0)then;  dt = dt + dt0;  dr = dr + dr0;  end if

   sim = 1.0D0 / (dt * (dr**df) * 10.0D0**(-bval * q2))
   dt_resc = dt*(10.d0**(-bval*q2*0.5d0))
   dr_resc = (dr**df)*(10.d0**(-bval*q2*0.5d0))
  end subroutine bp_similarity

  !==================================================================
  ! Build SI curve from nlink/kr histograms
  !==================================================================
  subroutine compute_si_curve(nlink, kr, n_event, ipmin, ipmax, &
                              sindex, k_si, ip_si, si_thresh, num_points)
    integer, intent(in)  :: nlink(-max_thr:max_thr), kr(-max_thr:max_thr)
    integer, intent(in)  :: n_event, ipmin, ipmax
    real(8), intent(out) :: sindex(max_si_points), si_thresh(max_si_points)
    integer, intent(out) :: k_si(max_si_points), ip_si(max_si_points), num_points
    
    integer, parameter :: smooth_step = 10
    integer :: ip, dkr, dnl
    real(8) :: zmexp

    num_points = 0
    do ip = ipmax, ipmin + smooth_step, -1
      if (kr(ip) == 1) exit
      dnl = nlink(ip - smooth_step) - nlink(ip)
      if (dnl <= 0) cycle
      dkr = kr(ip) - kr(ip - smooth_step)
      num_points = num_points + 1
      if (num_points > max_si_points) then
        print *, 'ERROR: max_points is too small for the SI curve.'
        stop 1
      end if
      k_si(num_points) = kr(ip)
      ip_si(num_points) = ip
      zmexp = real(dkr, 8) / real(dnl, 8)
      sindex(num_points) = zmexp * real(nlink(ip), 8) / real(n_event, 8)
      si_thresh(num_points) = exp(ip * bin0)
    end do
  end subroutine compute_si_curve

  !==================================================================
  ! Find SI minimum via progressive smoothing
  !==================================================================
  subroutine find_si_minimum(sindex, k_si, num_points, best_idx)
   implicit none
   real(8), intent(in) :: sindex(max_si_points)
   integer, intent(in) :: k_si(max_si_points), num_points
   integer, intent(out) :: best_idx

   real(8) :: smoothed(max_si_points), zbest, zmin, zmax1
   integer :: maxima(max_si_points), minima(max_si_points)
   integer :: nmx, nmn, sw, i, j, hw, nsum, tot_max
   integer :: imin, imax_left, ibest, kmin_shft
   real(8), parameter :: smooth_limit = 0.10D0, min_peak_pct = 0.05D0, max_peak_pct = 0.95D0

   tot_max = k_si(1)

   do sw = 2, int(num_points * smooth_limit)
     nmx = 0;  nmn = 0;  hw = sw / 2

     ! Moving average (correct denominator)
     do i = 1 + hw, num_points - hw
       smoothed(i) = 0.0D0;  nsum = 0
       do j = i - hw, i + hw
         smoothed(i) = smoothed(i) + sindex(j);  nsum = nsum + 1
       end do
       smoothed(i) = smoothed(i) / real(nsum, 8)
     end do

     ! Detect extrema
     do i = 1 + hw, num_points - hw
       if (smoothed(i) > smoothed(i-1) .and. smoothed(i) > smoothed(i+1) .and. &
           k_si(i) > int(tot_max*min_peak_pct) .and. k_si(i) < int(tot_max*max_peak_pct)) then
         nmx = nmx + 1;  maxima(nmx) = i
       else if (smoothed(i) < smoothed(i-1) .and. smoothed(i) < smoothed(i+1) .and. &
                k_si(i) > int(tot_max*min_peak_pct) .and. k_si(i) < int(tot_max*max_peak_pct)) then
         nmn = nmn + 1;  minima(nmn) = i
       end if
     end do

     if (nmn == 1 .and. nmx == 2) then
        imin = minima(1);  imax_left = maxima(2)
        zmin = smoothed(imin);  zmax1 = smoothed(imax_left)
        ! Shift 30% toward left peak in log-space
        zbest = exp(log(zmin) + (log(zmax1) - log(zmin)) * 0.3D0)
        do i = imin, num_points - hw
          if (smoothed(i) >= zbest) then;  ibest = i;  exit;  end if
        end do
        kmin_shft = k_si(ibest)
        best_idx = ibest
        return
     end if

     ! Fallback at smoothing limit
     if (sw == int(num_points * smooth_limit)) then
       call fallback_minimum(smoothed, k_si, num_points, hw,&
                             minima, nmn, best_idx)
       return
     end if
   end do
  end subroutine find_si_minimum

  !==================================================================
  ! Fallback when ideal 2-max/1-min not found
  !==================================================================
  subroutine fallback_minimum(smoothed, k_si, np, hw, minima, nmn, best_idx)
    implicit none
    real(8), intent(in) :: smoothed(max_si_points)
    integer, intent(in) :: k_si(max_si_points), np, hw, minima(max_si_points), nmn
    integer, intent(out) :: best_idx
    real(8) :: zmin, zmax1, zbest
    integer :: imin, imax_left, ibest, i, kmin_shft

    ! Global minimum
    zmin = 1.0D10;  imin = 1
    do i = 1, nmn
      if (smoothed(minima(i)) < zmin) then
        zmin = smoothed(minima(i));  imin = minima(i)
      end if
    end do
    ! Left peak (to the right of minimum in index space)
    zmax1 = -1.0D0;  imax_left = imin
    do i = imin + 1, np - hw
      if (smoothed(i) > zmax1) then
        zmax1 = smoothed(i);  imax_left = i
      end if
    end do
    ! Shift 30%
    zbest = exp(log(zmin) + (log(zmax1) - log(zmin)) * 0.3D0)
    ibest = imin
    do i = imin, np - hw
      if (smoothed(i) >= zbest) then;  ibest = i;  exit;  end if
    end do
    kmin_shft = k_si(ibest)
    best_idx = ibest
  end subroutine fallback_minimum

  !==================================================================
  ! Write SI curve to file
  ! Note: This subroutine writes the SI curve to a file,
  !       including the number of pairs in each similarity bin (nn) for reference.
  !==================================================================
  subroutine write_si_curve(sindex, k_si, ip_si, si_thresh, num_points, nn)
   implicit none
   real(8), intent(in) :: sindex(max_si_points), si_thresh(max_si_points)
   integer, intent(in) :: k_si(max_si_points), ip_si(max_si_points), num_points
   integer, intent(in) :: nn(-max_thr:max_thr)

   integer :: unit_number, ios, i

   open(newunit=unit_number, file=susc_index_file_name, status='replace', action='write', iostat=ios)
   if (ios /= 0) then
     print *, 'ERROR: cannot open SI output file: ', trim(susc_index_file_name)
     stop 1
   end if

   do i = 1, num_points
     write(unit_number, *) k_si(i), sindex(i), nn(ip_si(i)), si_thresh(i)
   end do

   close(unit_number)
  end subroutine write_si_curve

  subroutine write_output_catalog(ti, xxe, yye, mag, id_pr, n_event, ip_thr, nbg)
   implicit none
   real(8), intent(in) :: ti(max_events), xxe(max_events), yye(max_events), mag(max_events)
   integer, intent(in) :: id_pr(max_events), n_event, ip_thr
   integer, intent(out) :: nbg

   integer :: unit_number, ios, i, id

   open(newunit=unit_number, file=output_cat, status='replace', action='write', iostat=ios)
   if (ios /= 0) then
     print *, 'ERROR: cannot open catalog output file: ', trim(output_cat)
     stop 1
   end if

   nbg = 0

   do i = 1, n_event
     if (i == 1) then
       id = 0
     else if (id_pr(i) >= ip_thr) then
       id = 1
     else
       id = 0
     end if

     if (id == 0) nbg = nbg + 1
     write(unit_number, *) ti(i), xxe(i), yye(i), mag(i), id
   end do

   close(unit_number)
  end subroutine write_output_catalog

  subroutine write_summary(kmin_shft, similarity_threshold, ip_thr)
   implicit none
   integer, intent(in) :: kmin_shft, ip_thr
   real(8), intent(in) :: similarity_threshold

   integer :: unit_number, ios

   open(newunit=unit_number, file=final_results, status='replace', action='write', iostat=ios)
   if (ios /= 0) then
     print *, 'ERROR: cannot open summary output file: ', trim(final_results)
     stop 1
   end if
   write(unit_number, *) kmin_shft, ip_thr, similarity_threshold
    ! write(unit_number, '(A,I12)') 'Number of background events: ', kmin_shft
    ! write(unit_number, '(A,I12)') 'Threshold index: ', ip_thr
    ! write(unit_number, '(A,ES24.16)') 'Similarity threshold (1/sec*km): ', similarity_threshold
    close(unit_number)
  end subroutine write_summary

end module susceptibility_index_mod
