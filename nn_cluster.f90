! nn_cluster.f90
module nn_cluster_mod
  use global_params_mod
  implicit none
contains

  subroutine nn_assign(time, lat, lon, mag, n, bval, nc_thr, parent)
    !! For each event, find its nearest neighbor (BP metric) using temporal filter.
    !! parent(i) = index of NN father in sub-catalog (0 = root/background)
    integer, intent(in)  :: n
    real(8), intent(in)  :: time(n), lat(n), lon(n), mag(n), bval, nc_thr
    integer, intent(out) :: parent(n)
    integer :: i, j, cnt, id, jmin
    integer :: nn_list(n)
    real(8) :: tmax_list(n), tmax, nij, nij_min, ddt, ddr

    parent = 0;  cnt = 0
    do i = 2, n
      nij_min = 1.0D12;  jmin = 0
      ! Add i-1 to active list
      tmax = time(i-1) + 50.0D0 * (10.0D0**(1.4D0 * mag(i-1)))
      if (tmax > time(i)) then
        cnt = cnt + 1
        nn_list(cnt) = i - 1
        tmax_list(cnt) = tmax
      end if
      ! Prune expired entries
      j = 1
      do while (j <= cnt)
        if (tmax_list(j) < time(i)) then
          nn_list(j:cnt-1)   = nn_list(j+1:cnt)
          tmax_list(j:cnt-1) = tmax_list(j+1:cnt)
          cnt = cnt - 1
        else
          j = j + 1
        end if
      end do
      ! Search active list for nearest neighbor
      do j = 1, cnt
        id = nn_list(j)
        ddt = time(i) - time(id) + 1.0D0
        call distkm_nn(lat(i), lon(i), lat(id), lon(id), ddr)
        nij = ddt * (ddr**df) * 10.0D0**(-bval * mag(id))
        if (nij < nij_min) then
          nij_min = nij;  jmin = id
        end if
      end do
      if (nij_min < nc_thr) parent(i) = jmin
    end do
  end subroutine nn_assign

  subroutine find_mainshocks(mag_sub, parent, nev, idi, nevent, nflag)
    !! BFS cluster reconstruction on sub-catalog.
    !! Sets nflag in full-catalog space: 0 = mainshock, 1 = other.
    integer, intent(in)  :: nev, nevent, parent(nev), idi(nev)
    real(8), intent(in)  :: mag_sub(nev)
    integer, intent(out) :: nflag(nevent)
    ! Sparse children lists
    integer :: nchildren(nev), offset(nev), children_flat(nev)
    integer :: queue(nev)
    logical :: visited(nev)
    integer :: i, j, k, nevc, jmax, kkk
    real(8) :: qmax

    nflag = 1           ! default: everything is "not a mainshock"
    nchildren = 0
    ! Build children from parent array — O(n) instead of index_matr(n,n)
    ! Pass 1: count children per parent
    do i = 2, nev
      if (parent(i) > 0) then
        nchildren(parent(i)) = nchildren(parent(i)) + 1
      end if
    end do

    ! Build offsets (exclusive prefix sum)
    offset(1) = 0
    do i = 2, nev
      offset(i) = offset(i-1) + nchildren(i-1)
    end do

    ! Pass 2: fill flat children array
    nchildren = 0  ! reset to use as running counter
    do i = 2, nev
      if (parent(i) > 0) then
        j = parent(i)
        nchildren(j) = nchildren(j) + 1
        children_flat(offset(j) + nchildren(j)) = i
      end if
    end do

    ! do i = 2, nev
    !   if (parent(i) > 0) then
    !     j = parent(i)
    !     nchildren(j) = nchildren(j) + 1
    !     children(j, nchildren(j)) = i
    !   end if
    ! end do
    
    ! BFS over connected components
    visited = .false.
    do i = 1, nev
      if (visited(i)) cycle
      nevc = 1;  queue(1) = i;  qmax = -1.0D0;  jmax = i
      k = 1
      do while (k <= nevc)
        j = queue(k);  visited(j) = .true.
        if (mag_sub(j) >= qmax) then
          qmax = mag_sub(j);  jmax = j
        end if
        do kkk = 1, nchildren(j)
          nevc = nevc + 1
          queue(nevc) = children_flat(offset(j) + kkk)
        end do
        k = k + 1
      end do
      ! Mark the mainshock in full-catalog indexing
      nflag(idi(jmax)) = 0
    end do

  end subroutine find_mainshocks

  subroutine distkm_nn(x1, y1, x3, y3, dr)
    implicit none
    real(8), intent(in)  :: x1, y1, x3, y3
    real(8), intent(out) :: dr
    real(8), parameter :: pr = 3.14159265358979D0 / 180.0D0, dr0 = 0.01D0
    real(8) :: phi1, phi2, dphi, dlambda, a
    phi1 = x1 * pr;  phi2 = x3 * pr
    dphi = phi2 - phi1;  dlambda = (y3 - y1) * pr
    a = sin(dphi/2.0D0)**2 + cos(phi1)*cos(phi2)*sin(dlambda/2.0D0)**2
    dr = 2.0D0 * asin(sqrt(a)) * 6370.0D0
    if (dr == 0.0D0) dr = dr0
  end subroutine

  subroutine cluster_analysis(time, lat, lon, mag, nevent, bval, nc_thr, thspace_in, &
                            nfore_norm, naft_norm, indmain, nmain)
    implicit none

    integer, intent(in) :: nevent
    real(8), intent(in) :: time(nevent), lat(nevent), lon(nevent), mag(nevent)
    real(8), intent(in) :: bval, nc_thr, thspace_in(ncspace)
    real(8), intent(out) :: nfore_norm(ncmagnf,nctimef,ncspacef,ncmain)
    real(8), intent(out) :: naft_norm(ncmagn,nctime,ncspace,ncmain)
    integer, intent(out) :: indmain(ncmain), nmain(ncmain)

    integer :: parent(nevent)
    integer :: child_count(nevent), fill_count(nevent), offset(nevent), children_flat(nevent)
    integer :: queue(nevent)
    logical :: visited(nevent)

    integer :: nfore(ncmagnf,nctimef,ncspacef,ncmain)
    integer :: naft(ncmagn,nctime,ncspace,ncmain)

    real(8) :: thspace_loc(ncspace), ddt, ddr, qmax
    integer :: i, j, k, kkk, nevc, jmax, iqmax

    integer :: cnt, id, jmin, nn_list(nevent)
    real(8) :: tmax, nij, tmax_list(nevent), nij_min
    
    nfore_norm = 0.0D0
    naft_norm = 0.0D0
    indmain = 0
    nmain = 0
    nfore = 0
    naft = 0

    if (nevent <= 0) return

    parent = 0;  cnt = 0
    do i = 2, nevent
      nij_min = 1.0D12;  jmin = 0
      ! Add i-1 to active list
      tmax = time(i-1) + 50.0D0 * (10.0D0**(2.2D0 * mag(i-1)))
      if (tmax > time(i)) then
        cnt = cnt + 1
        nn_list(cnt) = i - 1
        tmax_list(cnt) = tmax
      end if
      ! Prune expired entries
      j = 1
      do while (j <= cnt)
        if (tmax_list(j) < time(i)) then
          nn_list(j:cnt-1)   = nn_list(j+1:cnt)
          tmax_list(j:cnt-1) = tmax_list(j+1:cnt)
          cnt = cnt - 1
        else
          j = j + 1
        end if
      end do
      ! Search active list for nearest neighbor
      do j = 1, cnt
        id = nn_list(j)
        ddt = time(i) - time(id) + 1.0D0
        call distkm_nn(lat(i), lon(i), lat(id), lon(id), ddr)
        nij = ddt * (ddr**df) * 10.0D0**(-bval * mag(id))
        if (nij < nij_min) then
          nij_min = nij;  jmin = id
        end if
      end do
      if (nij_min < nc_thr) parent(i) = jmin
    end do

    child_count = 0
    do i = 2, nevent
      if (parent(i) > 0) then
        child_count(parent(i)) = child_count(parent(i)) + 1
      end if
    end do

    offset(1) = 0
    do i = 2, nevent
      offset(i) = offset(i - 1) + child_count(i - 1)
    end do

    fill_count = 0
    children_flat = 0
    do i = 2, nevent
      if (parent(i) > 0) then
        j = parent(i)
        fill_count(j) = fill_count(j) + 1
        children_flat(offset(j) + fill_count(j)) = i
      end if
    end do

    visited = .false.

    do i = 1, nevent
      if (visited(i)) cycle

      nevc = 1
      queue(1) = i
      jmax = i
      qmax = -1.0D0
      k = 1

      do while (k <= nevc)
        j = queue(k)
        visited(j) = .true.

        if (mag(j) >= qmax) then
          qmax = mag(j)
          jmax = j
        end if

        do kkk = 1, fill_count(j)
          nevc = nevc + 1
          queue(nevc) = children_flat(offset(j) + kkk)
        end do

        k = k + 1
      end do

      if (qmax < mmain) cycle

      if (int(qmax) == 4) then
        iqmax = 1
      else if (int(qmax) == 5) then
        iqmax = 2
      else
        iqmax = 3
      end if

      nmain(iqmax) = nmain(iqmax) + 1

      thspace_loc = thspace_in
      thspace_loc(ncspace) = 0.01D0 * (10.0D0**(0.5D0 * mag(jmax)))

      do k = 1, nevc
        j = queue(k)
        if (j == jmax) cycle

        ddt = time(j) - time(jmax)
        call distkm_nn(lat(j), lon(j), lat(jmax), lon(jmax), ddr)

        if (ddt > 0.0D0) then
          call bin_event(mag(j), ddt, ddr, thtime, thspace_loc, thml, &
                        nctime, ncspace, ncmagn, step_mgn, naft(:,:,:,iqmax))
        else if (ddt < 0.0D0) then
          call bin_event(mag(j), -ddt, ddr, thtimef, thspacef, thmlf, &
                        nctimef, ncspacef, ncmagnf, step_mgnf, nfore(:,:,:,iqmax))
        end if
      end do
    end do

    do iqmax = 1, ncmain
      if (nmain(iqmax) == 0) cycle
      indmain(iqmax) = 1
      nfore_norm(:,:,:,iqmax) = real(nfore(:,:,:,iqmax), 8) / real(nmain(iqmax), 8)
      naft_norm(:,:,:,iqmax)  = real(naft(:,:,:,iqmax), 8)  / real(nmain(iqmax), 8)
    end do
  
  end subroutine cluster_analysis

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

  ! subroutine fb_bpzb(time,lat,lon,mag,nevent,ti,lati,loni,magi,idi,nev,bval,nfore_norm,naft_norm,indmain,nc,thspace_in)
  !  implicit none
  !  integer, intent(in) :: nevent,nev,idi(nev)
  !  real(8), intent(in) :: time(nevent),lat(nevent),lon(nevent),mag(nevent),bval,thspace_in(ncspace)
  !  real(8), intent(in) :: ti(nev), lati(nev), loni(nev), magi(nev), nc
  !  real(8), intent(out) :: nfore_norm(ncmagn,nctime,ncspace,ncmain),naft_norm(ncmagn,nctime,ncspace,ncmain)
  !  integer, intent(out) :: indmain(ncmain)
  !  real(8) :: nij_min,ddr,ddt,nij,qmax,thspace(ncspace)
  !  integer :: ntrigger,i,j,ntr,jmin,k,nevc,jmax,id,cnt,kkk
  !  integer :: iflag(nev),icluster(nev),invtrig(nev),nnafter(nev),itrigger(nev),nflag(nevent)
  !  integer :: n_main(ncmain),index_matr(nev,nev)
  !  integer :: nfore(ncmagn,nctime,ncspace,ncmain),naft(ncmagn,nctime,ncspace,ncmain)
  !  real(8) :: tmax,rmax
  !  integer :: nn_list(nev)
  !  real(8) :: tmax_list(nev)
  !  integer :: iqmax, itime, imfor, ispace, jjj
  !   thspace = thspace_in
  !   !! Initialize arrays !!
  !   n_main = 0
  !   nfore = 0
  !   naft = 0
  !   nfore_norm = 0
  !   naft_norm = 0
  !   index_matr = 0
  !   iflag=0
  !   nflag=1
  !   nnafter=0
  !   itrigger=0
  !   icluster=0
  !   invtrig=0
  !   cnt = 0 
  !   ntrigger=0
  !   do i=2,nev
  !     nij_min=1E12
  !     tmax=ti(i-1)+50*(10**(1.4*magi(i-1)))
  !     ! find in the list the events with tmax(j)<time(i) and take out
  !     if(tmax>ti(i))then
  !        cnt = cnt + 1
  !        nn_list(cnt) = i-1    ! index of the event
  !        tmax_list(cnt) = tmax 
  !     endif
  !     do jjj = cnt, 1, -1
  !        if (tmax_list(jjj) < ti(i)) then
  !           ! Shift remaining elements in tmax_list to remove the current event
  !           nn_list(jjj:cnt - 1) = nn_list(jjj + 1:cnt)
  !           tmax_list(jjj:cnt - 1) = tmax_list(jjj + 1:cnt)
  !           cnt = cnt - 1
  !        end if
  !     end do
  !     !do j=i-1,1,-1
  !     do j=1,cnt
  !       id = nn_list(j)  
  !       ddt=ti(i)-ti(id)+1.       ! seconds
  !       call distkm_nn(lati(i),loni(i),lati(id),loni(id),ddr)
  !       ! bval should change???
  !       nij=ddt*(ddr**df)*10**(-bval*magi(id)) ! BP metric in sec*km
  !       if(nij.lt.nij_min)then
  !         nij_min=nij ! minimum value of nij
  !         jmin=id     ! index of the nearest neighbor
  !       endif
  !     enddo
     
  !     if(nij_min.lt.nc)then
  !       ! father(i)=jmin ! Father of i
  !       if(iflag(jmin).eq.0)then       ! First time event jmin is a father
  !         iflag(jmin)=1                ! Flag the father event jmin
  !         ntrigger=ntrigger+1          ! Count the number of fathers
  !         itrigger(ntrigger)=jmin      ! Index of the father
  !         invtrig(jmin)=ntrigger       ! Number of the father jmin
  !       endif
  !       ntr=invtrig(jmin)
  !       nnafter(ntr)=nnafter(ntr)+1    ! Count offsprings of father jmin
  !       index_matr(ntr,nnafter(ntr))=i ! The index-i of each offspring of jmin father
  !     endif
  !   enddo
   
  !   do i=1,nev
  !     iflag(i)=0
  !   enddo
    
  !   do i=1,nev
  !     if(iflag(i).eq.0)then                     ! i is a father
  !       nevc=1
  !       icluster(nevc)=i
  !       k=1
  !       qmax=-1.
  !       do while (k.le.nevc+1)   
  !         j=icluster(k)
  !         iflag(j)=1
  !         if(magi(j).ge.qmax)then
  !           jmax=j      ! index of the mainshock
  !           qmax=magi(j) ! maximum magnitude
  !         endif
  !         if(invtrig(j).ne.0)then                ! Event j has offsprings
  !          if(nnafter(invtrig(j)).ge.1)then       ! number of offsprings >= 1
  !           do kkk=1,nnafter(invtrig(j))         ! loop over offsprings of father j
  !             nevc=nevc+1
  !             icluster(nevc)=index_matr(invtrig(j),kkk) ! get the offspring index
  !           enddo
  !          endif
  !         endif
  !         k=k+1
  !         if(k.gt.nevc)exit                  ! if no more offspring exit while loop
  !       enddo                                ! close while loop
  !       ! Flag mainshocks if it exceeds threshold mmain
  !       ! if(qmax.ge.mcl)nflag(idi(jmax))=0
  !       nflag(idi(jmax))=0
  !      endif  
  !   enddo
    
  !   do i=1,nevent
  !    if(nflag(i) /= 0.or.mag(i).lt.mmain)cycle
  !    if(int(mag(i)).eq.4)iqmax=1
  !    if(int(mag(i)).eq.5)iqmax=2
  !    if(int(mag(i)).ge.6)iqmax=3
  !    n_main(iqmax)=n_main(iqmax)+1        ! count mainshocks with M4, M5 and M6+
  !    thspace(5)=0.01*(10**(0.5*mag(i)))
     
  !    rmax=maxval(thspace) ! km
  !    tmax=maxval(thtime)  ! days
  !    ! count aftershocks
  !    do j=i+1,nevent
  !      ! Avoid counting a mainshock or aftershock already counted 
  !      if(nflag(j)==0.or.nflag(j)==2)cycle 
  !      ddt=(time(j)-time(i))/(3600.*24.) ! days
  !      if(ddt.gt.tmax)exit
  !      call distkm_nn(lat(j),lon(j),lat(i),lon(i),ddr)
  !      if(ddr>rmax)cycle
  !      nflag(j) = 2 ! flag as aftershock
  !      do itime=1,nctime
  !        do ispace=1,ncspace
  !          do imfor=1,ncmagn
  !            if(ddt.le.thtime(itime))then
  !              if(ddr.le.thspace(ispace))then
  !                if(mag(j).ge.thml(imfor).and.mag(j).lt.thml(imfor)+0.5)then
  !                  naft(imfor,itime,ispace,iqmax)=naft(imfor,itime,ispace,iqmax)+1
  !                endif
  !              endif
  !            endif
  !          enddo
  !        enddo
  !      enddo 
  !    enddo
     
  !    ! count foreshocks
  !    do j=i-1,1,-1
  !     ! Avoid counting a mainshock or aftershock already counted 
  !     if(nflag(j)==0.or.nflag(j)==2)cycle 
  !     ddt=(time(i)-time(j))/(3600.*24.) ! days
  !     if(ddt.gt.tmax)exit
  !     call distkm_nn(lat(j),lon(j),lat(i),lon(i),ddr)
  !     if(ddr.gt.rmax)cycle
  !     nflag(j) = 2 ! flag as foreshock
  !     do itime=1,nctime
  !       do ispace=1,ncspace
  !         do imfor=1,ncmagn
  !           if(ddt.le.thtime(itime))then
  !             if(ddr.le.thspace(ispace))then
  !               if(mag(j).ge.thml(imfor).and.mag(j).lt.thml(imfor)+0.5)then
  !                 nfore(imfor,itime,ispace,iqmax)=nfore(imfor,itime,ispace,iqmax)+1
  !               endif
  !             endif
  !           endif
  !         enddo
  !       enddo
  !     enddo 
  !    enddo

  !  enddo

  !  indmain = 0
  !  !! Normalize the foreshock and aftershock counts (events per mainshock) !!
  !  do iqmax=1,ncmain
  !   if(n_main(iqmax).ne.0)then
  !    indmain(iqmax) = 1 ! flag that indicates that mainshock class iqmax has mainshocks
  !    do imfor=1,ncmagn
  !      do itime=1,nctime
  !        do ispace=1,ncspace
  !        nfore_norm(imfor,itime,ispace,iqmax)=(nfore(imfor,itime,ispace,iqmax)*1.)/(n_main(iqmax)*1.)
  !        naft_norm(imfor,itime,ispace,iqmax)=(naft(imfor,itime,ispace,iqmax)*1.)/(n_main(iqmax)*1.)
  !        enddo
  !      enddo
  !    enddo
  !   endif
  !  enddo

  ! end subroutine fb_bpzb

end module nn_cluster_mod
