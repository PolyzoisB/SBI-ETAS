! ************************************************************************************************!
!    Program with steps                                                                           !
!    1. Give as input the earthquake catalog format [elapsed time, lat, lon, dep, mag]            !
!    2. Compute Susceptibility index, bg-rate and b-value                                         !
!    4. Compute normalized foreshock and aftershock stats                                         !
!    4. Export i) Input catalog stats ii) Input summary statistics                                ! 
!*************************************************************************************************!
program summary_stats
 use susceptibility_index_mod  
 use global_params_mod         
 use nn_cluster_mod
 use space_time_mag_count_mod
 implicit none
  
 real(8) :: thspace(ncspace)
 real(8) :: nn_thr, sim_thr
 integer, parameter :: bg_unit = 100

 ! aftershock/foreshock counts
 real(8) :: nfore_true(ncmagnf,nctimef,ncspacef,ncmain)
 real(8) :: naft_true(ncmagn,nctime,ncspace,ncmain)
 integer :: nmain_true(ncmain)

 ! Catalog variables
 real(8), allocatable :: time_true(:), lat_true(:), lon_true(:), mag_true(:)
 real(8), allocatable :: lt_bg_est(:), ln_bg_est(:)
 integer :: ii,nbg,ntrue
 real(8) :: bg_rate,bval, tr_to_tot, br_sup, br_inf, n_sup, n_inf, total_time
 
 ! Initialize the spatial intervals
 thspace(1:5) = (/3.0D0, 10.0D0, 20.0D0, 40.0D0, 0.0D0/)

 call load_input_catalog(input_catalog, time_true, lat_true, lon_true, mag_true, ntrue, bval)
 allocate(lt_bg_est(ntrue), ln_bg_est(ntrue))
 
 print*,'Number of events in the catalog: ',ntrue,' b-value: ',bval

 ! Compute Susceptibility Index and background rate
 call si_threshold(time_true, lat_true, lon_true, mag_true,&
                    ntrue, nn_thr, sim_thr, nbg, bval, lt_bg_est, ln_bg_est)
 print*, ' Estimated background events: ', nbg
 ! Export the background coordinates
 open(bg_unit,file=bg_coords,status='replace')
 do ii=1,nbg
   write(bg_unit,*)lt_bg_est(ii),ln_bg_est(ii)
 enddo
 close(bg_unit)
 
 tr_to_tot = (ntrue-nbg)*1. / ntrue*1.
 
 br_sup = 0.95
 br_inf = min(tr_to_tot, 0.99)
 
 n_sup = int(ntrue*(1.+0.4)) ! Upper limit for the number of events in the simulations
 n_inf = int(ntrue*(1.-0.4)) ! Lower limit for the number of events in the simulations
 
 total_time = time_true(ntrue) - time_true(1)
 bg_rate = nbg / ((lat_max-lat_min)*(lon_max-lon_min)*total_time) ! in events/sec/deg^2
 
 call compute_norm_stats(time_true, lat_true, lon_true, mag_true, ntrue, &
                          bval, nn_thr, thspace, &
                          nfore_true, naft_true, nmain_true)
  
 call write_summary_stats(input_aft_summary_stats, input_fore_summary_stats, &
                           naft_true, nfore_true, nmain_true)

 ! Export the catalog statistics
 call export_input_catalog_stats(input_catalog_stats, true_model_stats, bval, nn_thr, bg_rate, tr_to_tot, &
                                br_sup, br_inf, n_sup, n_inf, nbg, ntrue)

 contains
 
 subroutine load_input_catalog(path, time, lat, lon, mag, ncat, bval)
    character(len=*), intent(in) :: path
    real(8), allocatable, intent(out) :: time(:), lat(:), lon(:), mag(:)
    real(8), intent(out) :: bval
    integer, intent(out) :: ncat

    real(8) :: t, lt, ln, dep, mg, mn_mag
    integer :: i

    ncat = 0
    open(12, file=path, status='old')
    do
      read(12, *, end=100) t, lt, ln, dep, mg
      if (mg >= mc) ncat = ncat + 1
    end do
    100 close(12)
    if (ncat <= 0) stop 'No events loaded from input catalog.'
    allocate(time(ncat), lat(ncat), lon(ncat), mag(ncat))

    open(11, file=path, status='old')
    i = 0
    do
      read(11, *, end=101) t, lt, ln, dep, mg
      if (mg < mc) cycle
      i = i + 1
      time(i) = t
      lat(i)  = lt
      lon(i)  = ln
      mag(i)  = mg
      mn_mag = mn_mag + mag(i)
    end do
    101 close(11)
    mn_mag = mn_mag / (ncat*1.)
    ! MLE b-value estimation
    bval = log10(exp(1.))/(mn_mag-mc)
  end subroutine load_input_catalog

  subroutine compute_norm_stats(time, lat, lon, mag, ncat, bval_fixed, nc_fixed, thspace_in, &
                                nfore_norm, naft_norm, nmain)
    real(8), intent(in) :: time(ncat), lat(ncat), lon(ncat), mag(ncat)
    integer, intent(in) :: ncat
    real(8), intent(in) :: bval_fixed, nc_fixed, thspace_in(ncspace)
    real(8), intent(out) :: nfore_norm(ncmagnf,nctimef,ncspacef,ncmain)
    real(8), intent(out) :: naft_norm(ncmagn,nctime,ncspace,ncmain)
    integer, intent(out) :: nmain(ncmain)

    real(8) :: ti(ncat), lati(ncat), loni(ncat), magi(ncat)
    integer :: idi(ncat), parent(ncat), nflag(ncat)
    integer :: i, nev, flag_main(ncmain)

    nev = 0
    do i = 1, ncat
      if (mag(i) >= mcl) then
        nev = nev + 1
        ti(nev) = time(i)
        lati(nev) = lat(i)
        loni(nev) = lon(i)
        magi(nev) = mag(i)
        idi(nev) = i
      end if
    end do

    nfore_norm = 0.0D0
    naft_norm = 0.0D0

    if (nev <= 0) return

    !call cluster_analysis(ti, lati, loni, magi, nev, bval_fixed, nc_fixed, &
    !                      thspace_in, nfore_norm, naft_norm, indmain, nmain)  
    
    call nn_assign(ti, lati, loni, magi, nev, &
                  bval_fixed, nc_fixed, parent)

    call find_mainshocks(magi, parent, nev, idi, ncat, nflag)

    call count_fore_aft(time, lat, lon, mag, nflag, ncat, &
                       thspace_in, nfore_norm, naft_norm, nmain, flag_main)

  end subroutine compute_norm_stats

  subroutine write_summary_stats(path_a, path_f, naft, nfore, nmain)
    character(len=*), intent(in) :: path_a, path_f
    real(8), intent(in) :: naft(ncmagn,nctime,ncspace,ncmain)
    real(8), intent(in) :: nfore(ncmagnf,nctimef,ncspacef,ncmain)
    integer, intent(in) :: nmain(ncmain)

    integer :: iqmax, imfor, itime, ispace

    open(102, file=path_a, status='replace')
    do iqmax = 1, ncmain
      do imfor = 1, ncmagn
        do itime = 1, nctime
          do ispace = 1, ncspace
            write(102,*) naft(imfor,itime,ispace,iqmax), nmain(iqmax), iqmax + 3, &
                         thtime(itime)/t_day_to_sec, thspace(ispace), thml(imfor)
          end do
        end do
      end do
    end do
    close(102)
    pause
    
    open(103, file=path_f, status='replace')
    do iqmax = 1, ncmain
      do imfor = 1, ncmagnf
        do itime = 1, nctimef
          do ispace = 1, ncspacef
            write(103,*) nfore(imfor,itime,ispace,iqmax), nmain(iqmax), iqmax + 3, &
                         thtimef(itime)/t_day_to_sec, thspacef(ispace), thmlf(imfor)
          end do
        end do
      end do
    end do
    close(103)
  end subroutine write_summary_stats

  subroutine export_input_catalog_stats(cat_path, model_path, bval, nc, bg_rate, tr_to_tot, &
                                      br_sup, br_inf, n_sup, n_inf, nbg, ncat)
    character(len=*), intent(in) :: cat_path, model_path
    real(8), intent(in) :: bval, nc, bg_rate, tr_to_tot, br_sup, br_inf, n_sup, n_inf
    integer, intent(in) :: nbg, ncat

    open(201, file=cat_path, status='replace')
    write(201,*) bval
    write(201,*) bg_rate
    write(201,*) nc
    write(201,*) br_sup
    write(201,*) br_inf
    write(201,*) n_sup
    write(201,*) n_inf
    close(201)

    open(202, file=model_path, status='replace')
    write(202,*) 'Learning rate: ', lr
    write(202,*) 'Total time (years): ', tlast/t_year_to_sec, 'Auxilary window: ', tc/t_year_to_sec
    write(202,*) 'Branching Ratio: ', tr_to_tot, ' b-value (MLE): ', bval
    write(202,*) 'Decision Boundary: ', nc, ' Bg_rate (1/sec/deg^2): ', bg_rate
    write(202,*) 'Nbg: ', nbg, ' Total events: ', ncat
    write(202,*) 'Branching upper limit: ', br_sup, ' Branching lower limit: ', br_inf
    write(202,*) 'Cat size upper limit: ', n_sup, ' Cat size lower limit: ', n_inf
    close(202)

  end subroutine export_input_catalog_stats
  
end program summary_stats
