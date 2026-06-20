! ************************************************************************************************!
!    Program with steps                                                                           !
!    1. Import data with format [elapsed time (sec), longitude, latitude, depth, magnitude]       !
!    2. Compute Susceptibility index, bg-rate and b-value                                         !
!    3. Compute normalized foreshock and aftershock statistics from input catalog                 !
!    4. Compute normalized foreshock and aftershock statistics from simulated catalog             !
!    5. Repeat step 4. K1 times and compute cost function                                         !
!    6. Update parameter set of ETAS/ETAMS and repeat step 5.                                     !
!    7. Repeat step 6. K0 times                                                                   !
!    8. Export i) Estimated parameter set of ETAS/ETAMS ii) Cost function                         !
!*************************************************************************************************!
program sbi_estimation
 use global_params_mod      ! Import global parameters
 use models_mod             ! Import ETASI/ETABC simulation subroutine
 use nn_cluster_mod         ! Import clustering statistics subroutine
 use space_time_mag_count_mod ! Import space-time-magnitude count statistics subroutine
 implicit none
 !! COMPUTATIONAL TIME !!
 real :: elapsed_time
 integer :: clock_rate,start_time, end_time
 ! Export the Monte Carlo statistics
 character(len=200) :: monte_carlo_stats 
 ! Export initial parameters
 character(len=200) :: params_conv 
 ! Export the Monte Carlo sampling
 character(len=200) :: mc_sampling
 !! CLUSTER ANALYSIS !!
 ! Spatial intervals in km aft/for statistics
 real(8) :: thspace(ncspace)
 real(8) :: nc
 ! aftershock/foreshock counts
 real(8), dimension(ncmagnf,nctimef,ncspacef,ncmain) :: nfore_true
 real(8), dimension(ncmagn,nctime,ncspace,ncmain) :: naft_true
 integer :: im, imf, it, is
 real(8) :: f_nfore, f_nmainf, f_mmainf, f_tf, f_sf, f_mf
 real(8) :: f_naft, f_nmain, f_mmain, f_t, f_s, f_m
 ! Catalog variables
 real(8) :: lt_bg(max_events), ln_bg(max_events)
 real(8) :: bg_rate,bval, x
  ! ETAS Simulations
 integer :: ireal, seed1, i
 ! Monte carlo 
 real(8) :: alim(num_param,num_param),pm_init(num_param)
 real(8) :: cost_fn1, cost_new, param_set(num_param), param_out(num_param)
 real(8) :: br_sup , br_inf
 integer :: itc, stable_count
 integer :: flag_upd
 real(8) :: n_sup2, n_inf2 
 character(len=20) :: name(num_param)
 character(len=32) :: arg
 integer :: ios,seed,init_seed,nbg

 call get_command_argument(1, arg)   ! read first command-line argument
 read(arg, *, iostat=ios) seed
 if (ios /= 0) then
    seed = 4045  
    print *, "Error: could not read seed from command line."
    ! stop
 end if
 
 seed1 = -seed
 init_seed = seed
 !! Export files for inference step
 !write(monte_carlo_stats, '(A,I0,A)') 'results/cost_fn_', init_seed, 'lr015.txt'
 write(params_conv, '(A,I0,A)') 'results/sc_mc2.5/params_', init_seed,'.txt'
 !write(mc_sampling, '(A,I0,A)') 'results/mc_sampling_', init_seed, 'lr015.txt'

 !! Initialize the spatial intervals !!
 thspace(1:5) = (/3.0D0, 10.0D0, 20.0D0, 40.0D0, 0.0d0/) ! in km

 call import_bg_catalog(bg_coords, lt_bg, ln_bg, nbg)
 
 !! Import true foreshock and aftershock statistics
 open(101,file=input_aft_summary_stats,status='old')
   do im = 1, ncmain
    do imf = 1, ncmagn
      do it = 1, nctime
        do is = 1, ncspace
          read(101, *) f_naft, f_nmain, f_mmain, f_t, f_s, f_m
          naft_true(imf, it, is, im) = f_naft
        end do
      end do
    end do
  end do
  close(101)

  open(102,file=input_fore_summary_stats,status='old')
   do im = 1, ncmain
    do imf = 1, ncmagnf
      do it = 1, nctimef
        do is = 1, ncspacef
          read(102, *) f_nfore, f_nmainf, f_mmainf, f_tf, f_sf, f_mf
          nfore_true(imf, it, is, im) = f_nfore
        end do
      end do
    end do
  end do

 !! Import catalog statistics
 open(102,file=input_catalog_stats,status='old')
 read(102,*) bval
 read(102,*) bg_rate
 read(102,*) nc
 read(102,*) br_sup
 read(102,*) br_inf
 read(102,*) n_sup2
 read(102,*) n_inf2
 close(102)
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !       MONTE CARLO ESTIMATION              !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 !! INITIALIZE PARAMETERS !!
 106 call param_gen(param_set,alim,name,bg_rate,bval)
 pm_init = param_set
 
 call costfn(lt_bg,ln_bg,param_set,naft_true,nfore_true,cost_fn1,flag_upd)

 ! if many simulations out of bounds restart
 if(flag_upd == 1)then 
    print*,"Too many trials with too few or too many events, restart.."
    flush(6)
    go to 106
 endif

 !! MONTE CARLO Sampling !!
 !  open(60,file=monte_carlo_stats,status='replace') 
 !  open(61,file=mc_sampling,status='replace')
 !  do i=1,num_param
 !      write(60,*)0,0,0,name(i),param_set(i),0,0
 !  enddo
 

 ! Start the timer
 call system_clock(start_time)

 itc = 0
 stable_count = 0
 do ireal=1,K0 
  print*,"Monte carlo iteration: ", ireal
  flush(6)
  cost_new = cost_fn1
  itc = itc + 1
  call update_block([1,2], alim, param_set, name, itc, cost_fn1, param_out)
  param_set = param_out
  itc = itc + 1
  call update_block([3,4], alim, param_set, name, itc, cost_fn1, param_out)
  param_set = param_out
  itc = itc + 1
  call update_block([5,6,7], alim, param_set, name, itc, cost_fn1, param_out)
  param_set = param_out
       
  if (abs(cost_new - cost_fn1) < epsilon) then
        stable_count = stable_count + 1
  else
        stable_count = 0  ! Reset if condition not met
  end if
  print*,"Stable count: ", stable_count, " after iteration: ", ireal
  print*,"Current cost: ", cost_fn1
  flush(6)
  ! if (stable_count > int(conv_thr/2))flag_check = 1 ! switch to different learning rate
  if (stable_count >= conv_thr) then
    print*, "Stopping early at iteration ", ireal, " due to convergence."
    exit
  endif
 enddo 
 !close(60)
 !close(61)
 
 ! Stop the timer
 call system_clock(end_time, clock_rate)
 elapsed_time = real(end_time - start_time) / real(clock_rate)

 ! write best parameters to file
 open(103,file=params_conv,status='replace')
 ! Write all parameter values, initial values, cost, and elapsed time in a single row (no string labels)
 write(103,*) (param_set(i), i=1,num_param), (pm_init(i), i=1,num_param), cost_fn1, elapsed_time, ireal
 close(103)

 contains
 
 subroutine import_bg_catalog(path, lt_bg, ln_bg, nbg)
  character(len=*), intent(in) :: path
  real(8), intent(out) :: lt_bg(:), ln_bg(:)
  integer, intent(out) :: nbg

  real(8) :: x, y
  integer :: ii, jj

  jj = 0
  open(35,file=path,status='old')
  do ii=1,max_events
    read(35,*,end=199)x,y
    if(x.gt.lat_min.and.x.lt.lat_max)then
      if(y.gt.lon_min.and.y.lt.lon_max)then
        jj = jj + 1
        lt_bg(jj) = x
        ln_bg(jj) = y
      endif
    endif
  enddo
  199 close(35)
  nbg = jj
 end subroutine import_bg_catalog
 
 subroutine param_gen(param_set,param_bounds,name,bg_rate,bval)
  implicit none
  real(8), intent(out) :: param_set(num_param),param_bounds(num_param,num_param)
  character (LEN=20):: name(num_param)
  real(8), intent(in) :: bg_rate, bval
  real(8) :: n_branch,pmmax,pmmin
  integer :: i  
   
  !! p-value !!
  param_bounds(1,1) = 1.01
  param_bounds(1,2) = 1.35
  !107 param_set(1) = 1.15 
  107 param_set(1)= param_bounds(1,1)+ran2(seed1)*(param_bounds(1,2)-param_bounds(1,1))

  !! c-value (days) !!
  param_bounds(2,1) = 0.001
  param_bounds(2,2) = 0.06
  !param_set(2)=0.02
  param_set(2) = param_bounds(2,1)+ran2(seed1)*(param_bounds(2,2)-param_bounds(2,1))   
   
  !! alpha-value (exp) !!
  param_bounds(3,1) = 1.30
  param_bounds(3,2) = 2.80
  !param_set(3)=2.1 
  param_set(3) = param_bounds(3,1)+ran2(seed1)*(param_bounds(3,2)-param_bounds(3,1))
   
  !! K-value !!
  param_bounds(4,1) = 0.001
  param_bounds(4,2) = 0.37
  !param_set(4) = 0.09 
  param_set(4) = param_bounds(4,1)+ran2(seed1)*(param_bounds(4,2)-param_bounds(4,1))
   
  !! D-value (deg) !!
  param_bounds(5,1) = -6 ! 0.0000001
  param_bounds(5,2) = -4 ! 0.0001
  param_set(5)= (-6.+ran2(seed1)*(-4.+6.)) ! param_bounds(5,1)+ran2(seed1)*(param_bounds(5,2)-param_bounds(5,1))
  ! param_set(5)=log10(5E-5) 
   
  !! gamma-value (exp) !!
  param_bounds(6,1) = 0.5
  param_bounds(6,2) = 2.4
  ! param_set(6) = 1.5 
  param_set(6) = param_bounds(6,1)+ran2(seed1)*(param_bounds(6,2)-param_bounds(6,1))
   
  !! q-exponent !!
  param_bounds(7,1) = 1.05
  param_bounds(7,2) = 2.05
  !param_set(7) = 1.55 
  param_set(7) = param_bounds(7,1)+ran2(seed1)*(param_bounds(7,2)-param_bounds(7,1))
   
  !!! Aftershock incompleteness !!!
  !! tau-value (sec) !!
  param_bounds(8,1) = 0
  param_bounds(8,2) = 300
  !param_set(8) = 200. 
  param_set(8) = 0.0d0 ! param_bounds(8,1)+ran2(seed1)*(param_bounds(8,2)-param_bounds(8,1))

  !! dr-value (km) !!
  param_bounds(9,1) = 30
  param_bounds(9,2) = 70
  param_set(9)=50. ! param_bounds(9,1)+ran2(seed1)*(param_bounds(9,2)-param_bounds(9,1))
 
  !! Bg-rate (1/sec/deg^2) !!
  param_bounds(10,1) = 1E-15
  param_bounds(10,2) = 1E-1
  param_set(10)=bg_rate   
   
  !! b-value !!
  param_bounds(11,1) = 0.5
  param_bounds(11,2) = 1.5
  param_set(11) = bval 

  ! !! bf-value !!
  ! ! param_bounds(12,1) = 0.5
  ! ! param_bounds(12,2) = 1.5
  ! param_set(12) = 1.0 ! param_bounds(12,1)+ran2(seed1)*(param_bounds(12,2)-param_bounds(12,1))

  ! !! pf-value !!
  ! param_bounds(13,1) = 0.5
  ! param_bounds(13,2) = 1.5
  ! param_set(13) = 1.5 ! param_bounds(13,1)+ran2(seed1)*(param_bounds(13,2)-param_bounds(13,1))
    
  ! !! cf-value !!
  ! param_bounds(14,1) = 0.5
  ! param_bounds(14,2) = 1.5
  ! param_set(14) = 0.01 ! param_bounds(14,1)+ran2(seed1)*(param_bounds(14,2)-param_bounds(14,1))
    
  ! !! ddf-value !!
  ! param_bounds(15,1) = 0.5
  ! param_bounds(15,2) = 1.5
  ! param_set(15) = log10(5E-5) ! param_bounds(15,1)+ran2(seed1)*(param_bounds(15,2)-param_bounds(15,1))
    
  ! !! gammaf-value !!
  ! param_bounds(16,1) = 0.5
  ! param_bounds(16,2) = 1.5
  ! param_set(16) = 1.0 ! param_bounds(16,1)+ran2(seed1)*(param_bounds(16,2)-param_bounds(16,1))
    
  ! !! qf-value !!
  ! param_bounds(17,1) = 0.5
  ! param_bounds(17,2) = 1.5
  ! param_set(17) = 1.8 !param_bounds(17,1)+ran2(seed1)*(param_bounds(17,2)-param_bounds(17,1))
    
  ! !! phi-value !!
  ! param_bounds(18,1) = 0.
  ! param_bounds(18,2) = 0.4
  ! param_set(18) = param_bounds(18,1)+ran2(seed1)*(param_bounds(18,2)-param_bounds(18,1))
    
  ! Safety check
  do i=1,num_param
    pmmin = param_bounds(i,1)
    pmmax = param_bounds(i,2)
    if(param_set(i).lt.pmmin.or.param_set(i).gt.pmmax)then
      print*,"Error: parameter ", i, " out of bounds during initialization."
      stop
    endif
    if(pmmax <= pmmin)then  
      print*, 'Error: pmmax must be greater than pmmin.'
      stop
    endif
  enddo

  !! Parameter names !!
  name(1) = 'p'
  name(2) = 'c'
  name(3) = 'alpha'
  name(4) = 'K'
  name(5) = 'd'
  name(6) = 'gamma'
  name(7) = 'q'
  name(8) = 'tau-ETASI'
  name(9) = 'dr-ETASI'
  name(10) = 'bg-rate'
  name(11) = 'b-value'
  ! name(12) = 'bf-value'
  ! name(13) = 'pf-value'
  ! name(14) = 'cf-value'
  ! name(15) = 'df-value'
  ! name(16) = 'gammaf-value'
  ! name(17) = 'qf-value'
  ! name(18) = 'phi-value'
   
  ! Branching ratio condition
  call branch_rt(param_set(4),param_set(3),bval,n_branch)   
  if(n_branch.ge.br_sup.or.n_branch.lt.br_inf)go to 107
  print*,"Initial branching ratio: ",n_branch
 end subroutine param_gen
 
 subroutine update_parameter(pm, pmmin, pmmax, pm_new)
  implicit none
  ! Inputs
  REAL(8), INTENT(IN)  :: pm       ! Current parameter value
  REAL(8), INTENT(IN)  :: pmmin    ! Minimum allowed value
  REAL(8), INTENT(IN)  :: pmmax    ! Maximum allowed value
  ! Output
  REAL(8), INTENT(OUT) :: pm_new   ! Updated parameter

  ! Local variables
  REAL(8) :: pm_scaled, delta_scaled, u, pm_scaled_new
  
  ! Normalize current value to [0, 1]
  pm_scaled = (pm - pmmin) / (pmmax - pmmin)

  ! Generate a uniform random number u ~ U(0, 1)
  234 x = ran2(seed1) 
  u = ran2(seed1)
  !! Reflecting boundaries !!
  delta_scaled = lr*u*sign(1d0,x-0.5)
  pm_scaled_new = pm_scaled + delta_scaled
  if (pm_scaled_new < 0.0d0.or.pm_scaled_new > 1.0d0) go to 234
  
  ! Update in scaled space and clamp to [0,1]
  !! Truncated distribution !!
  ! lower-bound max( pm_scaled-lr, 0.0d0 ) and upper-bound min( pm_scaled+lr, 1.0d0)
  !pm_scaled_new = max(pm_scaled - lr(1), 0.0d0) + & 
  !    u * (min(pm_scaled + lr(1), 1.0d0) - max(pm_scaled - lr(1), 0.0d0))
  !! Cencored distribution !!
  ! delta_scaled = lr(1)*u*sign(1d0,x-0.5)
  !pm_scaled_new = pm_scaled + delta_scaled
  !pm_scaled_new = max(0.0, min(1.0, pm_scaled_new))
    
  
  !! logit transform !!
  !z = log(pm_scaled / (1.0d0 - pm_scaled))
  ! Additive update
  !z = z + delta_scaled
  ! back to [0,1]
  !pm_scaled_new = 1.0d0 / (1.0d0 + exp(-z))
  
  ! Convert back to original space
  pm_new = pmmin + pm_scaled_new * (pmmax - pmmin)

 end subroutine update_parameter

 subroutine costfn(lat_bg,lon_bg,param_set,true_stats_aft,true_stats_fore,cost_fn,flag_upd)
   implicit none
   real(8), intent(in) :: lat_bg(:),lon_bg(:)
   real(8), intent(in) :: param_set(num_param)
   real(8), intent(in) :: true_stats_aft(ncmagn,nctime,ncspace,ncmain), true_stats_fore(ncmagnf,nctimef,ncspacef,ncmain)
   real(8), intent(out) :: cost_fn
   
   integer :: kloop
   real(8), allocatable :: t_sim(:),lat_sim(:),lon_sim(:),mag_sim(:)
   real(8) :: nfore_sim(ncmagnf,nctimef,ncspacef,ncmain),naft_sim(ncmagn,nctime,ncspace,ncmain)
   real(8) :: naft_sum(ncmagn,nctime,ncspace,ncmain),nfore_sum(ncmagnf,nctimef,ncspacef,ncmain)
   real(8) :: nevent_avg, bval
   integer :: iqmax,imfor,itime,ispace,nev_sim,nevent_sum, nout, flag_upd
   integer :: nmain_sims(ncmain), nmain(ncmain), flag_main(ncmain)
   
   !!!!!!!! MODEL SIMULATIONS !!!!!!!!!!!!!!!!!!
   flag_upd = 0
   cost_fn = 0.
   naft_sum = 0.0
   nfore_sum = 0.0
   nevent_avg = 0.0
   nevent_sum = 0
   nmain_sims = 0
   nout = 0
   bval = param_set(11) ! b-value
   
   !! Compute average summary statistics !!
   do kloop=1,K1
     !15 call etams_sim_cat(lat_bg,lon_bg,param_set,t_sim,lat_sim,lon_sim,mag_sim,nev_sim)
     15 call etasi_sim(lat_bg,lon_bg, nbg, param_set,t_sim,lat_sim,lon_sim,mag_sim,nev_sim,seed1)
     if(nev_sim.gt.n_sup2.or.nev_sim.lt.n_inf2)then
      nout = nout + 1 
      if(nout.gt.10*K1)then
        print*,"Too many trials with too few or too many events, new update.."
        flush(6)
        flag_upd = 1
        goto 16
      endif
      goto 15
     endif

     call compute_norm_stats(t_sim, lat_sim, lon_sim, mag_sim, nev_sim, &
                          bval, nc, thspace, &
                          nfore_sim, naft_sim, nmain, flag_main)
   
     ! Accumulate naft_norm across iterations
     ! grid cells associated with zero mainshocks in this iteration are not counted
     ! to avoid biasing the average lower
     naft_sum = naft_sum + naft_sim
     nfore_sum = nfore_sum + nfore_sim
     nmain_sims = nmain_sims + flag_main ! count how many times each mainshock class had at least one event
     nevent_sum = nevent_sum + nev_sim
   enddo
   
   nevent_avg = nevent_sum / (K1*1.)
   !! Print the number of trials !!
   print*, "Number out of bounds: ", nout, "Events per sim: ", int(nevent_avg)
   flush(6)
   
   !! Cost function !!
   do iqmax=1,ncmain
      do imfor=1,ncmagn
        do itime=1,nctime
          do ispace=1,ncspace
           naft_sim(imfor,itime,ispace,iqmax) = naft_sum(imfor,itime,ispace,iqmax) / max(1,nmain_sims(iqmax))*1.0D0
           if(true_stats_aft(imfor,itime,ispace,iqmax).ne.0.)then
                cost_fn = cost_fn + ((naft_sim(imfor,itime,ispace,iqmax) / &
                        true_stats_aft(imfor,itime,ispace,iqmax)) - 1 )**2
                !print*,cost_fn, imfor, itime, ispace, iqmax, naft_sim(imfor,itime,ispace,iqmax), &
                !        true_stats_aft(imfor,itime,ispace,iqmax)
           endif
          enddo
        enddo
      enddo

      ! do imfor=1,ncmagnf
      !   do itime=1,nctimef
      !     do ispace=1,ncspacef
      !        nfore_sim(imfor,itime,ispace,iqmax) = nfore_sum(imfor,itime,ispace,iqmax) / max(1,nmain_sims(iqmax))*1.0D0
      !       if(true_stats_fore(imfor,itime,ispace,iqmax).ne.0.)then
      !            cost_fn = cost_fn + ((nfore_norm(imfor,itime,ispace,iqmax) / &
      !                      true_stats_fore(imfor,itime,ispace,iqmax)) - 1 )**2
      !            !print*,cost_fn, imfor, itime, ispace, iqmax, nfore_norm(imfor,itime,ispace,iqmax), &
      !            !        true_stats_fore(imfor,itime,ispace,iqmax)
      !       endif
      !     enddo
      !   enddo
      ! enddo

    enddo

    16 return 

 end subroutine costfn

 subroutine compute_norm_stats(time, lat, lon, mag, ncat, bval_fixed, nc_fixed, thspace_in, &
                                nfore_norm, naft_norm, nmain, flag_main)
    implicit none
    real(8), intent(in) :: time(ncat), lat(ncat), lon(ncat), mag(ncat)
    integer, intent(in) :: ncat
    real(8), intent(in) :: bval_fixed, nc_fixed, thspace_in(ncspace)
    real(8), intent(out) :: nfore_norm(ncmagnf,nctimef,ncspacef,ncmain)
    real(8), intent(out) :: naft_norm(ncmagn,nctime,ncspace,ncmain)
    integer, intent(out) :: nmain(ncmain), flag_main(ncmain)

    real(8) :: ti(ncat), lati(ncat), loni(ncat), magi(ncat)
    integer :: idi(ncat), parent(ncat), nflag(ncat)
    integer :: i, nev

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

  subroutine branch_rt(K,alpha,bval,n_branch)
  implicit none
  real(8), intent(in) :: K,alpha,bval
  real(8), intent(out) :: n_branch
  real(8) :: beta,max_mag,min_mag

  beta = bval * log(10.0)  ! Convert b-value to beta
  max_mag = msup
  min_mag = mc
  ! Calculate branching ratio Seif et al
  if(alpha.eq.beta)then
    n_branch = (K * beta * (max_mag - min_mag)) / &
          (1 - exp(-beta * (max_mag - min_mag)))
  else
    n_branch = (K * beta * (1 - exp(-(max_mag - min_mag) * (beta - alpha)))) / &
                ((beta - alpha) * (1 - exp(-beta * (max_mag - min_mag))))
  endif
  
 end subroutine branch_rt

 subroutine distkm(x1, y1, x3, y3, dr)
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
 
  !  subroutine etas(dt,dr,magn,prob)
  !    implicit none
  !    real*8 :: dt,dr,magn,prob,espo,dm,fdr
  !    real*8, parameter :: alpha = 1.06,K = 0.038,p = 1.27
  !    real*8, parameter :: gamma = 0.85,c = 0.024,qdec = 1.506
  !    real*8, parameter :: dd = 0.006 
  !      espo=10**(alpha*(magn-mc))
  !      dm=dd*10**(gamma*(magn-mc))
  !      fdr=((dr**2+dm)**(-qdec))*(dm**(-(1-qdec)))
  !      prob=espo*(dt+c)**(-p)*fdr
  !      prob=prob*K*(c**(p-1))*(p-1)*(-(1-qdec))*(1./pr)   
  !  end subroutine etas 

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
 
 subroutine update_block(block_indices, alim, paramin, name, itc, cost_fn1, paramout)
   implicit none
   integer, intent(in) :: block_indices(:), itc
   real(8), intent(inout) :: cost_fn1
   real(8), intent(in) :: alim(num_param,num_param), paramin(num_param)
   character(len=*), intent(in) :: name(num_param)
   real(8), intent(out) :: paramout(num_param)
   real(8) :: cost_fn, adum(num_param), param_set(num_param)
   integer :: i, ii, flag_upd, cnt
   real(8) :: n_branch, cost_fn_old

   param_set = paramin
   cnt = 0
   do
    do i=1, size(block_indices)
      ii = block_indices(i)
      adum(ii) = param_set(ii)
      call update_parameter(adum(ii), alim(ii,1), alim(ii,2), param_set(ii))
    enddo

    ! Check branching ratio once after updating the whole block
    call branch_rt(param_set(4),param_set(3),param_set(11),n_branch)
    if (n_branch.ge.br_sup .or. n_branch.lt.br_inf) then
      ! revert and retry another proposal for this block
      do i=1, size(block_indices)
        ii = block_indices(i)
        param_set(ii) = adum(ii)
      end do
      cycle
    end if

    call costfn(lt_bg,ln_bg,param_set,naft_true,nfore_true,cost_fn,flag_upd)
    if(flag_upd == 0)exit  ! success
    cnt = cnt + 1
    if (cnt > 10) then
      print*, "Too many re-tries, exit the optimization procedure !!!"
      flush(6)
      stop
    end if
    ! reset params on failure
    do i=1, size(block_indices)
      ii = block_indices(i)
      param_set(ii) = adum(ii)
    end do
   
   end do
   
   !  do i=1, size(block_indices)
   !       ii = block_indices(i)
   !       write(61,*)ireal,itc,cost_fn,name(ii),param_set(ii),adum(ii),cost_fn1
   !       call flush(61)  
   !  end do

   ! accept/reject based on cost
   !ratio = cost_fn1 / cost_fn
   !u = ran2(seed1)
   !if(u.lt.ratio) then    ! Metropolis-Hastings acceptance
   if (cost_fn < cost_fn1) then
      cost_fn_old = cost_fn1
      cost_fn1 = cost_fn
      !do i=1, size(block_indices)
        !ii = block_indices(i)
        !write(60,*)ireal,itc,cost_fn1,cost_fn_old,name(ii),param_set(ii),adum(ii)
        !call flush(60)
      !end do
   else
      do i=1, size(block_indices)
         ii = block_indices(i)
         param_set(ii) = adum(ii)
      end do
   end if
   
   paramout = param_set

 end subroutine update_block

end program sbi_estimation
