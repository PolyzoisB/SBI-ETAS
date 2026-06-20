module global_params_mod
  !! FREE PARAMETERS !!
  implicit none
  
  !! PRE-PROCESS STEP !!
  
  !==============================
  ! File paths
  !==============================
  character(len=*), parameter :: input_catalog = 'datasets/sc_1981_31_8_2023_all_conv.txt'
  character(len=*), parameter :: susc_index_file_name = 'results/susc_index_pavl_v4.txt'
  character(len=*), parameter :: output_cat = 'results/output_catalog_pavl_v4.txt'
  character(len=*), parameter :: final_results = 'results/final_results_pavl_v4.txt'
  character(len=*), parameter :: rescaled_distances = 'results/rescaled_distances_pavl_v4.txt'

  ! Export the foreshock and aftershock summary statistics
  character(len=*), parameter :: input_aft_summary_stats = 'true_aft_sum_stats.txt'
  character(len=*), parameter :: input_fore_summary_stats = 'results/sc_mc2.5/true_fore_sum_stats.txt'

  character(len=*), parameter :: input_catalog_stats = 'results/sc_mc2.5/true_catalog_stats.txt'
  character(len=*), parameter :: true_model_stats = 'results/sc_mc2.5/true_model_stats.txt'
  character(len=*), parameter :: bg_coords = 'datasets/bg_coords_sc_mc3.0.txt'
  character(len=*), parameter :: bg_coords_1 = '/home/eugenio/Polyzois/etas/etas_simulation/datasets/cat_main_m2.dat'
  character(len=*), parameter :: bg_coords_2 = '/home/eugenio/Polyzois/etas/etas_simulation/datasets/cat_main_m4.dat'

  integer, parameter :: max_events = 1500000
    
  !==============================
  ! Susceptibility index parameters
  !==============================
  integer, parameter :: max_thr = 6000
  real(8), parameter :: pthmin = 1.0D-18
  real(8), parameter :: pthmax = 1.0D15
  real(8), parameter :: bin0 = 0.05D0

  !==============================
  ! Metric parameters (Baiesi–Paczuski)
  !==============================
  real(8), parameter :: df   = 1.6d0     ! fractal dimension
  ! real(8), parameter :: bval = 1.0d0     ! Gutenberg–Richter b-value

  !==============================
  ! Smoothing parameters
  !==============================
  integer, parameter :: max_si_points = 10000
    
  !! INFERENCE STEP !!
  !! Number of parameters to estimate
  integer, parameter :: num_param = 11 
  !! Temporal scales !!
  real(8), parameter :: t_day_to_sec = 24*3600.0, t_year_to_sec = 365.25*24.*3600.0
  
  !! Number of classes for aftershock statistics
  integer, parameter :: nctime=6, ncspace=5, ncmagn=2, ncmain=3
  !! Temporal intervals for aftershock statistics (days) !!
  real(8), parameter :: thtime(nctime) = (/0.25*t_day_to_sec, 0.5*t_day_to_sec, &
                        0.75*t_day_to_sec, 1.0*t_day_to_sec, 3.0*t_day_to_sec, 10.0*t_day_to_sec/) 
  !! Magnitude aftershock classes for aftershock statistics (left boundary) !!
  real(8), parameter :: thml(ncmagn) = (/3.0, 3.5 /)           
  !! Magnitude interval for aftershock statistics !!
  real(8), parameter :: step_mgn=0.5
  
  !! Number of classes for foreshock statistics !!
  integer, parameter :: nctimef=3, ncspacef=3, ncmagnf=1
  !! Temporal intervals for foreshock statistics (days) !!
  real(8), parameter :: thtimef(nctimef) = (/1*t_day_to_sec, 5.0*t_day_to_sec, 10.0*t_day_to_sec/)
  !! Spatial intervals for foreshock statistics (km) !!
  real(8), parameter :: thspacef(ncspacef) = (/10.0, 20.0, 40.0/)
  !! Magnitude foreshock classes for foreshock statistics (left boundary) !!
  real(8), parameter :: thmlf(ncmagnf) = (/3.0/)
  !! Magnitude interval for foreshock statistics !!
  real(8), parameter :: step_mgnf=1.0

  !!!!!! SIMULATION BOUNDARIES !!!!!!
  !! lat/lon bounds for simulations (degrees) - for California, USA !!
  !real(8), parameter :: lat_min = 32.0, lat_max = 37.0, lon_min = 114.0, lon_max = 121.0 
  real(8), parameter :: lat_min = 32.0, lat_max = 37.0, lon_min = -121.0, lon_max = -114.0 
  !! Temporal boundaries for simulations (years) !!
  !real(8), parameter :: tc = 5.0*365.25*24.*3600., tlast = 25.*365.25*24.*3600. 
  real(8), parameter :: tc = 5.0*t_year_to_sec, tlast = 42.67*t_year_to_sec+tc
  !! Magnitude boundaries for simulations !!
  real(8), parameter :: mc = 3.0, msup = 7.5

  !!!!!!! OTHER PARAMETERS !!!!!!
  !! Magnitude boundaries for cluster analysis (approximation)
  real(8), parameter :: mmain = 4.0, mcl = 3.5
  
  !! Other parameters
  real(8), parameter :: pr=3.14159/180.
  real(8), parameter :: epsilon = 1.0E-2        !! Convergence criterion
  integer, parameter :: K0 = 200                !! Monte carlo iterations
  integer, parameter :: K1 = 30                 !! Average statistics
  integer, parameter :: conv_thr = 60           !! Consecutive stable iterations for convergence
  ! integer, parameter :: seed = 5099           ! Seed for random number generator
  real(8) :: lr = 0.1d0                        !! Learning rate for parameters update

end module global_params_mod
