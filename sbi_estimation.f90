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
 use iso_fortran_env, only: real64
 use global_params_mod        ! Import global parameters
 use rng_mod, only: rng_init, rng_uniform
 use models_mod               ! Import ETASI simulation subroutine
 use nn_cluster_mod           ! Import clustering statistics subroutine
 use space_time_mag_count_mod ! Import space-time-magnitude count statistics subroutine
 implicit none
 !! COMPUTATIONAL TIME !!
 real :: elapsed_time
 integer :: clock_rate,start_time, end_time
 ! Export parameters
 character(len=200) :: params_conv 
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
 real(8), allocatable :: lt_bg(:), ln_bg(:)
 real(8) :: bg_rate , b_val
  ! ETAS Simulations
 integer :: ireal, i
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
 
 !! Initialize the spatial intervals !!
 thspace(1:5) = (/3.0D0, 10.0D0, 20.0D0, 40.0D0, 0.0d0/) ! in km

 call get_command_argument(1, arg)   ! read first command-line argument
 read(arg, *, iostat=ios) seed
 if (ios /= 0) then
    seed = 4000  
    print *, "Set seed for random generator 4000: could not read seed from command line."
 end if
 
 call rng_init(seed)
 init_seed = seed
 !! Export files for inference step
 write(params_conv, '(A,I0,A)') 'results/params_', init_seed,'.txt'
 
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
  close(102)

  !! Import catalog statistics
  open(103,file=input_catalog_stats,status='old')
  read(103,*) b_val
  read(103,*) bg_rate
  read(103,*) nc
  read(103,*) br_sup
  read(103,*) br_inf
  read(103,*) n_sup2
  read(103,*) n_inf2
  close(103)
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !       MONTE CARLO ESTIMATION              !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 !! INITIALIZE PARAMETERS !!
 106 call param_gen(param_set,alim,name,bg_rate,b_val)
 pm_init = param_set
 
 call costfn(lt_bg,ln_bg,param_set,naft_true,nfore_true,cost_fn1,flag_upd)

 ! if many simulations out of bounds restart
 if(flag_upd == 1)then 
    print*,"Too many trials with too few or too many events, restart.."
    flush(6)
    go to 106
 endif

 ! Start the timer
 call system_clock(start_time)

 itc = 0
 stable_count = 0
 do ireal = 1,K0 
  print*,"Monte carlo iteration: ", ireal
  flush(6)
  cost_new = cost_fn1
  itc = itc + 1
  call update_block([1,2,8], alim, param_set, cost_fn1, param_out)
  param_set = param_out
  itc = itc + 1
  call update_block([3,4,8], alim, param_set, cost_fn1, param_out)
  param_set = param_out
  itc = itc + 1
  call update_block([5,6,7], alim, param_set, cost_fn1, param_out)
  param_set = param_out
       
  if (abs(cost_new - cost_fn1) < epsilon) then
        stable_count = stable_count + 1
  else
        stable_count = 0  ! Reset if condition not met
  end if
  print*,"Stable count: ", stable_count, " after iteration: ", ireal
  print*,"Current cost: ", cost_fn1
  flush(6)
  if (stable_count >= conv_thr) then
    print*, "Stopping early at iteration ", ireal, " due to convergence."
    exit
  endif
 enddo 
 
 ! Stop the timer
 call system_clock(end_time, clock_rate)
 elapsed_time = real(end_time - start_time) / real(clock_rate)

 ! write best parameters to file
 open(103,file=params_conv,status='replace')
 ! Write all parameter values, initial values, cost, and elapsed time in a single row (no string labels)
 write(103,*) (param_set(i), i=1,num_param), (pm_init(i), i=1,num_param), cost_fn1, elapsed_time, ireal
 close(103)

 contains
 
 subroutine import_bg_catalog(path, lat_bg, lon_bg, nbg)
   implicit none
   character(len=*), intent(in) :: path
   real(8), allocatable, intent(out) :: lat_bg(:), lon_bg(:)
   integer, intent(out) :: nbg

   integer :: unit_number
   integer :: ios
   integer :: i
   real(8) :: lat_value, lon_value

   ! First pass: count the number of coordinate pairs.
   open(newunit=unit_number, file=path, status='old', action='read', iostat=ios)

   if (ios /= 0) then
      print *, 'ERROR: cannot open background catalog: ', trim(path)
      error stop
   end if

   nbg = 0
   do
      read(unit_number, *, iostat=ios) lat_value, lon_value

      if (ios < 0) exit

      if (ios > 0) then
        close(unit_number)
        print *, 'ERROR: malformed row in background catalog: ', &
               trim(path)
        error stop
      endif
      
      nbg = nbg + 1
   end do

   if (nbg == 0) then
      close(unit_number)
      print *, 'ERROR: background catalog contains no coordinates: ', &
               trim(path)
      error stop
   end if

   ! Allocate only the required amount of memory.
   allocate(lat_bg(nbg), lon_bg(nbg), stat=ios)

   ! Second pass: read the coordinates.
   rewind(unit_number)

   do i = 1, nbg
      read(unit_number, *, iostat=ios) lat_bg(i), lon_bg(i)

      if (ios /= 0) then
         close(unit_number)
         print *, 'ERROR: failed while reading background row ', i
         error stop
      end if
   end do

   close(unit_number)

 end subroutine import_bg_catalog

 subroutine param_gen(param_set,param_bounds,name,bg_rate,b_val)
  implicit none
  real(8), intent(out) :: param_set(num_param),param_bounds(num_param,num_param)
  character (LEN=20):: name(num_param)
  real(8), intent(in) :: bg_rate, b_val
  real(8) :: n_branch,pmmax,pmmin
  integer :: i  
   
  !! p-value !!
  param_bounds(1,1) = 1.01
  param_bounds(1,2) = 1.40
  !107 param_set(1) = 1.15 
  107 param_set(1)= param_bounds(1,1)+rng_uniform()*(param_bounds(1,2)-param_bounds(1,1))

  print*, n_branch

  !! c-value (days) !!
  param_bounds(2,1) = 0.001
  param_bounds(2,2) = 0.1
  !param_set(2)=0.02
  param_set(2) = param_bounds(2,1)+rng_uniform()*(param_bounds(2,2)-param_bounds(2,1))   
   
  !! alpha-value (exp) !!
  param_bounds(3,1) = 1.30
  param_bounds(3,2) = 2.80
  !param_set(3)=2.1 
  param_set(3) = param_bounds(3,1)+rng_uniform()*(param_bounds(3,2)-param_bounds(3,1))
   
  !! K-value !!
  param_bounds(4,1) = 0.001
  param_bounds(4,2) = 0.50
  !param_set(4) = 0.09 
  param_set(4) = param_bounds(4,1)+rng_uniform()*(param_bounds(4,2)-param_bounds(4,1))
   
  !! D-value (deg) !!
  param_bounds(5,1) = 1E-7
  param_bounds(5,2) = 1E-4
  param_set(5)= (param_bounds(5,1)+rng_uniform()*(param_bounds(5,2)-param_bounds(5,1)))
   
  !! gamma-value (exp) !!
  param_bounds(6,1) = 0.5
  param_bounds(6,2) = 2.4
  ! param_set(6) = 1.5 
  param_set(6) = param_bounds(6,1)+rng_uniform()*(param_bounds(6,2)-param_bounds(6,1))
   
  !! q-exponent !!
  param_bounds(7,1) = 1.05
  param_bounds(7,2) = 2.5
  !param_set(7) = 1.55 
  param_set(7) = param_bounds(7,1)+rng_uniform()*(param_bounds(7,2)-param_bounds(7,1))
   
  !!! Aftershock incompleteness !!!
  !! tau-value (sec) !!
  param_bounds(8,1) = 0
  param_bounds(8,2) = 300
  !param_set(8) = 200. 
  param_set(8) = param_bounds(8,1)+rng_uniform()*(param_bounds(8,2)-param_bounds(8,1))

  !! dr-value (km) !!
  param_bounds(9,1) = 30
  param_bounds(9,2) = 70
  param_set(9) = 50. ! param_bounds(9,1)+rng_uniform()*(param_bounds(9,2)-param_bounds(9,1))
 
  !! stdv-value !!
  param_bounds(12,1) = 0.1 
  param_bounds(12,2) = 0.6
  param_set(12) = 0.4

  !! Bg-rate (1/sec/deg^2) !!
  param_bounds(10,1) = 1E-15
  param_bounds(10,2) = 1E-1
  param_set(10) = bg_rate   
   
  !! b-value !!
  param_bounds(11,1) = 0.5
  param_bounds(11,2) = 1.5
  param_set(11) = b_val 
  
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
  name(12) = 'stdv'

  ! Branching ratio condition
  call branch_rt(param_set(4),param_set(3),b_val,n_branch)   
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
  REAL(8) :: pm_scaled, delta_scaled, u, x, pm_scaled_new
  
  ! Normalize current value to [0, 1]
  pm_scaled = (pm - pmmin) / (pmmax - pmmin)

  ! Generate a uniform random number u ~ U(0, 1)
  234 x = rng_uniform() 
  u = rng_uniform()
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
   integer, intent(out) :: flag_upd
   
   integer :: kloop
   real(8), allocatable :: t_sim(:),lat_sim(:),lon_sim(:),mag_sim(:)
   real(8) :: nfore_sim(ncmagnf,nctimef,ncspacef,ncmain),naft_sim(ncmagn,nctime,ncspace,ncmain)
   real(8) :: naft_sum(ncmagn,nctime,ncspace,ncmain),nfore_sum(ncmagnf,nctimef,ncspacef,ncmain)
   real(8) :: nevent_avg, b_val
   integer :: iqmax,imfor,itime,ispace,nev_sim,nevent_sum, nout
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
   b_val = param_set(11) ! b-value

   !! Compute average summary statistics !!
   do kloop=1,K1
     15 call etasi_sim(lat_bg, lon_bg, param_set, t_sim, lat_sim, lon_sim, mag_sim, nev_sim)
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
                          b_val, nc, thspace, &
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
   !print*, "Number of simulations out of bounds: ", nout, "Events per sim: ", int(nevent_avg)
   !flush(6)
   
   !! Cost function !!
   do iqmax=1,ncmain

    if (nmain_sims(iqmax) == 0) then
      flag_upd = 1
      cost_fn = huge(cost_fn)
      return
    endif
    
    do imfor=1,ncmagn
      do itime=1,nctime
        do ispace=1,ncspace
           naft_sim(imfor,itime,ispace,iqmax) = naft_sum(imfor,itime,ispace,iqmax) / real(nmain_sims(iqmax), kind=real64)
           if(true_stats_aft(imfor,itime,ispace,iqmax).ne.0.)then
                cost_fn = cost_fn + ((naft_sim(imfor,itime,ispace,iqmax) / &
                        true_stats_aft(imfor,itime,ispace,iqmax)) - 1 )**2
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
    nmain = 0
    flag_main = 0

    if (nev <= 0) return

    !call cluster_analysis(ti, lati, loni, magi, nev, bval_fixed, nc_fixed, &
    !                      thspace_in, nfore_norm, naft_norm, indmain, nmain)  
    
    call nn_assign(ti, lati, loni, magi, nev, &
                  bval_fixed, nc_fixed, parent)

    call find_mainshocks(magi, parent, nev, idi, ncat, nflag)

    call count_fore_aft(time, lat, lon, mag, nflag, ncat, &
                       thspace_in, nfore_norm, naft_norm, nmain, flag_main)

  end subroutine compute_norm_stats

  subroutine branch_rt(K,alpha,b_val,n_branch)
  implicit none
  real(8), intent(in) :: K,alpha,b_val
  real(8), intent(out) :: n_branch
  real(8) :: beta,max_mag,min_mag

  beta = b_val * log(10.0)  ! Convert b-value to beta
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
 
 subroutine update_block(block_indices, alim, paramin, cost_fn1, paramout)
   implicit none
   integer, intent(in) :: block_indices(:)
   real(8), intent(inout) :: cost_fn1
   real(8), intent(in) :: alim(num_param,num_param), paramin(num_param)
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
   
   ! accept/reject based on cost
   !ratio = cost_fn1 / cost_fn
   !u = rng_uniform()
   !if(u.lt.ratio) then    ! Metropolis-Hastings acceptance
   if (cost_fn < cost_fn1) then
      cost_fn_old = cost_fn1
      cost_fn1 = cost_fn
   else
      do i=1, size(block_indices)
         ii = block_indices(i)
         param_set(ii) = adum(ii)
      end do
   end if
   
   paramout = param_set

 end subroutine update_block

end program sbi_estimation
