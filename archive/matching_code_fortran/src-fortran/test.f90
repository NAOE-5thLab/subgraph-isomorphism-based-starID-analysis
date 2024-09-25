
program src_fortran
    use subgraph_isomorphism_starid_detector, only: matching_set_for_analysis
    !
    implicit none
    integer :: FILEID_I = 1
    integer :: FILEID_HR = 2
    integer :: FILEID_RA = 3
    integer :: FILEID_DE = 4
    integer :: FILEID_Vmag = 5
    integer :: FILEID_Ipair = 6
    integer :: FILEID_Theta = 7
    !
    integer :: n_obs = 4
    !integer :: obs_setID(10) = (/918, 641, 638, 1024, 731, 960, 752, 1098, 775, 941/)
    integer :: obs_setID(4) = (/842, 939, 922, 886/)
    double precision :: s_hat_set(4, 3)
    double precision :: epsilon = 1.0d0 * 3.14d0/180d0
    !
    integer :: N_starset, N_pairset, N_candi
    integer, allocatable :: I_star(:), HR(:)
    double precision, allocatable :: RA(:), DE(:), Vmag(:)
    integer, allocatable :: pairidset(:, :)
    double precision, allocatable :: thetaset(:)
    !
    integer, allocatable :: N_candi_setid(:), candi_setid(:, :, :)
    double precision, allocatable :: time(:)
    !
    integer :: i
    double precision :: x, y
    double precision :: alpha, delta

    ! LOAD CSV !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! catalog !!!
    open(FILEID_I, file='starDB_I.csv')
    N_starset = 0
    do
        read (FILEID_I, *, end=100) x
        N_starset = N_starset + 1
    end do
100 continue
    !
    allocate(I_star(N_starset))
    allocate(HR(N_starset))
    allocate(RA(N_starset))
    allocate(DE(N_starset))
    allocate(Vmag(N_starset))
    ! ID
    rewind (FILEID_I)
    do i = 1, N_starset
        read (FILEID_I, *) I_star(i)
    end do
    close (FILEID_I)
    ! HR
    open(FILEID_HR, file='starDB_HR.csv')
    do i = 1, N_starset
        read (FILEID_HR, *) HR(i)
    end do
    close (FILEID_HR)
    ! RA
    open(FILEID_RA, file='starDB_RA.csv')
    do i = 1, N_starset
        read (FILEID_RA, *) RA(i)
    end do
    close (FILEID_RA)
    ! DE
    open(FILEID_DE, file='starDB_DE.csv')
    do i = 1, N_starset
        read (FILEID_DE, *) DE(i)
    end do
    close (FILEID_DE)
    ! Vmag
    open(FILEID_Vmag, file='starDB_Vmag.csv')
    do i = 1, N_starset
        read (FILEID_Vmag, *) Vmag(i)
    end do
    close (FILEID_Vmag)
    
    !!! db !!!
    open(FILEID_Ipair, file='I_pair_FOV.csv')
    N_pairset = 0
    do
        read (FILEID_Ipair, *, end=200) x, y
        N_pairset = N_pairset + 1
    end do
200 continue
    !
    allocate(pairidset(N_pairset, 2))
    allocate(thetaset(N_pairset))
    ! Theta_pair_FOV
    rewind (FILEID_Ipair)
    do i = 1, N_pairset
        read (FILEID_Ipair, *) pairidset(i, 1), pairidset(i, 2)
    end do
    close (FILEID_Ipair)
    ! Theta_pair_FOV
    open(FILEID_Theta, file='Theta_pair_FOV.csv')
    do i = 1, N_pairset
        read (FILEID_Theta, *) thetaset(i)
    end do
    close (FILEID_Theta)
    
    ! obs stars !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, n_obs
        alpha = RA(obs_setID(i)+1)
        delta = DE(obs_setID(i)+1)
        s_hat_set(i, 1) = cos(alpha) * cos(delta)
        s_hat_set(i, 2) = sin(alpha) * cos(delta)
        s_hat_set(i, 3) = sin(delta)
    end do
    print *, obs_setID(:n_obs)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! obs_setID(i) -> obs_setID(i)+1 !!!
    !!! python�Ƃ̔z��v�f�ԍ��̈Ⴂ       !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! matching !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(N_candi_setid(n_obs-1))
    allocate(candi_setid(n_obs-1, N_pairset, n_obs))
    allocate(time(n_obs))
    N_candi = N_pairset
    call matching_set_for_analysis(&
            N_candi_setid, candi_setid, time, &
            N_candi, n_obs, s_hat_set(:n_obs, :), epsilon, &
            N_starset, RA, DE, N_pairset, thetaset, pairidset)
    
    do i = 1, n_obs-1
        print *, i+1, N_candi_setid(i), time(i+1) - time(1)
    end do
    print *, 'end'

end program src_fortran
    
    

