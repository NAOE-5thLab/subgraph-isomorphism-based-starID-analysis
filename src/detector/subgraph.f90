module subgraph_isomorphism_starid_detector
    implicit none;

contains

    subroutine matching_set_for_analysis(&
            N_candi_setid, candi_setid, time, &
            n_obs, s_hat_set, epsilon, &
            N_set, RA_set, DE_set, N_pairset, thetaset, pairidset)
        implicit none;
        integer, intent(in) :: n_obs
        double precision, intent(in) :: s_hat_set(n_obs, 3)
        double precision, intent(in) :: epsilon
        integer, intent(in) :: N_set, N_pairset
        double precision, intent(in) :: RA_set(N_set), DE_set(N_set)
        double precision, intent(in) :: thetaset(N_pairset)
        integer, intent(in) :: pairidset(N_pairset, 2)
        integer, intent(out) :: N_candi_setid(n_obs-1)
        integer, intent(out) :: candi_setid(n_obs-1, N_pairset, n_obs)
        integer, intent(out) :: time(n_obs)
        ! 
        integer :: n
        ! 
        do n = 2, n_obs
            call cpu_time(time(n-1))
            if (n == 2) then
                call matching_pair(&
                    N_candi_setid(n - 1), candi_setid(n - 1, :, :n), &
                    s_hat_set, epsilon, N_pairset, thetaset, pairidset)
            else if (n == 3) then
                call matching_triangle(&
                    N_candi_setid(n - 1), candi_setid(n - 1, :, :n), &
                    s_hat_set(:n, :), epsilon, &
                    N_set, RA_set, DE_set, N_pairset, thetaset, pairidset)
            else
                call matching_multiangle(&
                    N_candi_setid(n - 1), candi_setid(n - 1, :, :n), &
                    n, N_candi_setid(n - 2), candi_setid(n - 2, :, :n - 1), &
                    s_hat_set(:n, :), epsilon, N_pairset, thetaset, pairidset)
            end if
        end do
        call cpu_time(time(n_obs))
    end subroutine matching_set_for_analysis

    subroutine matching_set(&
            N_candi_setid, candi_setid, &
            n_obs, s_hat_set, epsilon, &
            N_set, RA_set, DE_set, N_pairset, thetaset, pairidset)
        implicit none;
        integer, intent(in) :: n_obs
        double precision, intent(in) :: s_hat_set(n_obs, 3)
        double precision, intent(in) :: epsilon
        integer, intent(in) :: N_set, N_pairset
        double precision, intent(in) :: RA_set(N_set), DE_set(N_set)
        double precision, intent(in) :: thetaset(N_pairset)
        integer, intent(in) :: pairidset(N_pairset, 2)
        integer, intent(out) :: N_candi_setid
        integer, intent(out) :: candi_setid(N_pairset, n_obs)
        !
        integer :: pre_N_candi_setid
        integer, allocatable :: pre_candi_setid(:, :)
        integer :: N_candi_setid_temp
        integer, allocatable :: candi_setid_temp(:, :)
        ! 
        integer :: n
        !
        if (n_obs == 2) then
            call matching_pair(&
                N_candi_setid, candi_setid, &
                s_hat_set, epsilon, N_pairset, thetaset, pairidset)
        else if (n_obs > 2) then
            do n = 3, n_obs
                if (n == 3) then
                    allocate(candi_setid_temp(N_pairset, n))
                    call matching_triangle(&
                        N_candi_setid_temp, candi_setid_temp, &
                        s_hat_set(:n, :), epsilon, &
                        N_set, RA_set, DE_set, N_pairset, thetaset, pairidset)
                else
                    deallocate(candi_setid_temp)
                    allocate(candi_setid_temp(N_pairset, n))
                    call matching_multiangle(&
                        N_candi_setid_temp, candi_setid_temp, &
                        n, pre_N_candi_setid, pre_candi_setid, &
                        s_hat_set(:n, :), epsilon, N_pairset, thetaset, pairidset)
                    deallocate(pre_candi_setid)
                end if
                allocate(pre_candi_setid(N_pairset, n))
                pre_N_candi_setid = N_candi_setid_temp
                pre_candi_setid = candi_setid_temp
            end do
            N_candi_setid = N_candi_setid_temp
            candi_setid = candi_setid_temp
        else
            print*, 'error : n_obs'
        end if

    end subroutine matching_set

    subroutine matching_pair(&
            N_candi_pairid_ij, candi_pairid_ij, &
            s_hat_set, epsilon, N_pairset, thetaset, pairidset)
        implicit none;
        double precision, intent(in) :: s_hat_set(2, 3)
        double precision, intent(in) :: epsilon
        integer, intent(in) :: N_pairset
        double precision, intent(in) :: thetaset(N_pairset)
        integer, intent(in) :: pairidset(N_pairset, 2)
        integer, intent(out) :: N_candi_pairid_ij
        integer, intent(out) :: candi_pairid_ij(N_pairset, 2)
        !
        double precision :: s_hat_i(3), s_hat_j(3), theta_hat_ij
        !
        s_hat_i = s_hat_set(1, :)
        s_hat_j = s_hat_set(2, :)
        call inter_angle_3dvec(theta_hat_ij, s_hat_i, s_hat_j)
        call get_pairids_of_candi_theta(&
            N_candi_pairid_ij, candi_pairid_ij, &
            theta_hat_ij, epsilon, &
            N_pairset, thetaset, pairidset)
    end subroutine matching_pair

    subroutine matching_triangle(&
            N_candi_triangleid_ijk, candi_triangleid_ijk, &
            s_hat_set, epsilon, &
            N_set, RA_set, DE_set, N_pairset, thetaset, pairidset)
        implicit none;
        double precision, intent(in) :: s_hat_set(3, 3)
        double precision, intent(in) :: epsilon
        integer, intent(in) :: N_set, N_pairset
        double precision, intent(in) :: RA_set(N_set), DE_set(N_set)
        double precision, intent(in) :: thetaset(N_pairset)
        integer, intent(in) :: pairidset(N_pairset, 2)
        integer, intent(out) :: N_candi_triangleid_ijk
        integer, intent(out) :: candi_triangleid_ijk(N_pairset, 3)
        !
        double precision :: s_hat_i(3), s_hat_j(3), s_hat_k(3)
        double precision :: theta_hat_ij, theta_hat_ik, theta_hat_jk
        integer :: N_candi_triangleid_ijk_temp, candi_triangleid_ijk_temp(N_pairset, 3)
        integer :: N_candi_pairid_ij, candi_pairid_ij(N_pairset, 2)
        integer :: N_candi_pairid_ik, candi_pairid_ik(N_pairset, 2)
        integer :: N_candi_pairid_jk, candi_pairid_jk(N_pairset, 2)
        integer :: candi_pairid_ij_1d(N_pairset), candi_pairid_ik_1d(N_pairset)
        integer :: N_candi_pairid_ij_1d_uni, candi_pairid_ij_1d_uni(N_pairset)
        integer :: N_candi_pairid_ik_1d_uni, candi_pairid_ik_1d_uni(N_pairset)
        integer :: N_candi_starid_i, candi_starid_i(N_pairset)
        integer :: N_candi_starid_j, candi_starid_j(N_pairset)
        integer :: N_candi_starid_k, candi_starid_k(N_pairset)
        double precision :: obs_sign, candi_sign
        double precision :: vec_i(3), vec_j(3), vec_k(3)
        integer :: i, i1, i2, i3, i4, i5, i6, i7, i8
        integer :: count, count_tiangleid
        !
        s_hat_i = s_hat_set(1, :)
        s_hat_j = s_hat_set(2, :)
        s_hat_k = s_hat_set(3, :)
        call inter_angle_3dvec(theta_hat_ij, s_hat_i, s_hat_j)
        call inter_angle_3dvec(theta_hat_ik, s_hat_i, s_hat_k)
        call inter_angle_3dvec(theta_hat_jk, s_hat_j, s_hat_k)
        call get_pairids_of_candi_theta(&
            N_candi_pairid_ij, candi_pairid_ij, &
            theta_hat_ij, epsilon, N_pairset, thetaset, pairidset)
        call get_pairids_of_candi_theta(&
            N_candi_pairid_ik, candi_pairid_ik, &
            theta_hat_ik, epsilon, N_pairset, thetaset, pairidset)
        call get_pairids_of_candi_theta(&
            N_candi_pairid_jk, candi_pairid_jk, &
            theta_hat_jk, epsilon, N_pairset, thetaset, pairidset)
        ! select i
        count = 0
        ! convert 1d ij
        candi_pairid_ij_1d(:N_candi_pairid_ij) = candi_pairid_ij(:N_candi_pairid_ij, 1)
        candi_pairid_ij_1d(N_candi_pairid_ij+1:2*N_candi_pairid_ij) = candi_pairid_ij(:N_candi_pairid_ij, 2)
        call unique_intarray1d(N_candi_pairid_ij_1d_uni, candi_pairid_ij_1d_uni, &
                2*N_candi_pairid_ij, candi_pairid_ij_1d(:2*N_candi_pairid_ij))
        ! convert 1d ik
        candi_pairid_ik_1d(:N_candi_pairid_ik) = candi_pairid_ik(:N_candi_pairid_ik, 1)
        candi_pairid_ik_1d(N_candi_pairid_ik+1:2*N_candi_pairid_ik) = candi_pairid_ik(:N_candi_pairid_ik, 2)
        call unique_intarray1d(N_candi_pairid_ik_1d_uni, candi_pairid_ik_1d_uni, &
                2*N_candi_pairid_ik, candi_pairid_ik_1d(:2*N_candi_pairid_ij))
        !
        do i1 = 1, N_candi_pairid_ij_1d_uni
            if (any(candi_pairid_ij_1d_uni(i1) == candi_pairid_ik_1d_uni(1:N_candi_pairid_ik_1d_uni))) then
                count = count + 1
                candi_starid_i(count) = candi_pairid_ij_1d_uni(i1)
            end if
        end do
        N_candi_starid_i = count
        !
        count_tiangleid = 0
        do i2 = 1, N_candi_starid_i
            ! select j-th candicate
            count = 0
            do i3 = 1, N_candi_pairid_ij
                if (candi_pairid_ij(i3, 1) == candi_starid_i(i2)) then
                    count = count + 1
                    candi_starid_j(count) = candi_pairid_ij(i3, 2)
                else if (candi_pairid_ij(i3, 2) == candi_starid_i(i2)) then
                    count = count + 1
                    candi_starid_j(count) = candi_pairid_ij(i3, 1)
                end if
            end do
            N_candi_starid_j = count
            ! select k-th candicate
            count = 0
            do i4 = 1, N_candi_pairid_ik
                if (candi_pairid_ik(i4, 1) == candi_starid_i(i2)) then
                    count = count + 1
                    candi_starid_k(count) = candi_pairid_ik(i4, 2)
                else if (candi_pairid_ik(i4, 2) == candi_starid_i(i2)) then
                    count = count + 1
                    candi_starid_k(count) = candi_pairid_ik(i4, 1)
                end if
            end do
            N_candi_starid_k = count
            ! select j&k-th pair stars candidate
            do i5 = 1, N_candi_starid_j
                do i6 = 1, N_candi_starid_k    
                    if (candi_starid_j(i5) < candi_starid_k(i6)) then
                        do i7 = 1, N_candi_pairid_jk
                            if (candi_pairid_jk(i7, 1) == candi_starid_j(i5) &
                                    .and. candi_pairid_jk(i7, 2) == candi_starid_k(i6)) then
                                count_tiangleid = count_tiangleid + 1
                                candi_triangleid_ijk_temp(count_tiangleid, 1) = candi_starid_i(i2)
                                candi_triangleid_ijk_temp(count_tiangleid, 2) = candi_starid_j(i5)
                                candi_triangleid_ijk_temp(count_tiangleid, 3) = candi_starid_k(i6)
                                exit;
                            end if
                        end do        
                    else
                        do i7 = 1, N_candi_pairid_jk
                            if (candi_pairid_jk(i7, 1) == candi_starid_k(i6) &
                                    .and. candi_pairid_jk(i7, 2) == candi_starid_j(i5)) then
                                count_tiangleid = count_tiangleid + 1
                                candi_triangleid_ijk_temp(count_tiangleid, 1) = candi_starid_i(i2)
                                candi_triangleid_ijk_temp(count_tiangleid, 2) = candi_starid_j(i5)
                                candi_triangleid_ijk_temp(count_tiangleid, 3) = candi_starid_k(i6)
                                exit;
                            end if
                        end do  
                    end if                    
                end do
            end do
        end do
        N_candi_triangleid_ijk_temp = count_tiangleid
        ! check specular sign
        count_tiangleid = 0
        call triangle_specular_sign(obs_sign, s_hat_i, s_hat_j, s_hat_k)
        do i8 = 1, N_candi_triangleid_ijk_temp
            call equatorial2vec(vec_i, &
                    RA_set(candi_triangleid_ijk_temp(i8, 1) + 1), &
                    DE_set(candi_triangleid_ijk_temp(i8, 1) + 1))
            call equatorial2vec(vec_j, &
                    RA_set(candi_triangleid_ijk_temp(i8, 2) + 1), &
                    DE_set(candi_triangleid_ijk_temp(i8, 2) + 1))
            call equatorial2vec(vec_k, &
                    RA_set(candi_triangleid_ijk_temp(i8, 3) + 1), &
                    DE_set(candi_triangleid_ijk_temp(i8, 3) + 1))
            ! 
            call triangle_specular_sign(candi_sign, vec_i, vec_j, vec_k)
            if (obs_sign*candi_sign > 0) then
                count_tiangleid = count_tiangleid + 1
                candi_triangleid_ijk(count_tiangleid, :) = candi_triangleid_ijk_temp(i8, :)
            end if
        end do
        N_candi_triangleid_ijk = count_tiangleid
    end subroutine matching_triangle

    subroutine matching_multiangle(&
            N_candi_setid, candi_setid, &
            n_obs, pre_N_candi_setid, pre_candi_setid, &
            s_hat_set, epsilon, N_pairset, thetaset, pairidset)
        implicit none;
        integer, intent(in) :: n_obs
        integer, intent(in) :: pre_N_candi_setid
        integer, intent(in) :: pre_candi_setid(N_pairset, n_obs-1)
        double precision, intent(in) :: s_hat_set(n_obs, 3)
        double precision, intent(in) :: epsilon
        integer, intent(in) :: N_pairset
        double precision, intent(in) :: thetaset(N_pairset)
        integer, intent(in) :: pairidset(N_pairset, 2)
        integer, intent(out) :: N_candi_setid
        integer, intent(out) :: candi_setid(N_pairset, n_obs)
        !
        integer :: N_candi_pairid_in(n_obs-1), candi_pairid_in(n_obs-1, N_pairset, 2)
        double precision :: s_hat_i(3), s_hat_n(3), theta_hat_in
        integer :: candi_starid_i
        integer :: N_candi_starid_n_from_i, candi_starid_n_from_i(N_pairset)
        integer :: N_candi_starid_n_temp, candi_starid_n_temp(N_pairset)
        integer :: N_candi_starid_n, candi_starid_n(N_pairset)
        !
        integer :: i1, i2, i3, i4, i5, i6
        integer :: count, count_setid
        !
        count_setid = 0
        ! matching pair stars including last star
        do i1 = 1, n_obs-1
            s_hat_i = s_hat_set(i1, :)
            s_hat_n = s_hat_set(n_obs, :)
            call inter_angle_3dvec(theta_hat_in, s_hat_i, s_hat_n)
            call get_pairids_of_candi_theta(&
                N_candi_pairid_in(i1), candi_pairid_in(i1, :, :), &
                theta_hat_in, epsilon, N_pairset, thetaset, pairidset)
        end do
        ! determine candidate of last star in each set of matched stars except last
        do i2 = 1, pre_N_candi_setid
            do i3 = 1, n_obs-1
                candi_starid_i = pre_candi_setid(i2, i3)
                count = 0
                do i4 = 1, N_candi_pairid_in(i3)
                    if (candi_pairid_in(i3, i4, 1) == candi_starid_i) then
                        count = count + 1
                        candi_starid_n_from_i(count) = candi_pairid_in(i3, i4, 2)
                    else if (candi_pairid_in(i3, i4, 2) == candi_starid_i) then
                        count = count + 1
                        candi_starid_n_from_i(count) = candi_pairid_in(i3, i4, 1)
                    end if
                end do
                N_candi_starid_n_from_i = count
                !
                if (i3 == 1) then
                    N_candi_starid_n = N_candi_starid_n_from_i
                    candi_starid_n = candi_starid_n_from_i
                else
                    count = 0
                    do i5 = 1, N_candi_starid_n
                        if (any(candi_starid_n(i5) == candi_starid_n_from_i(1:N_candi_starid_n_from_i))) then
                            count = count + 1
                            candi_starid_n_temp(count) = candi_starid_n(i5)
                        end if
                    end do
                    N_candi_starid_n = count
                    candi_starid_n = candi_starid_n_temp
                end if
                !
                if (N_candi_starid_n == 0) exit
            end do
            !
            do i6 = 1, N_candi_starid_n
                count_setid = count_setid + 1
                candi_setid(count_setid, :) = [pre_candi_setid(i2, :), candi_starid_n(i6)]
            end do
        end do
        N_candi_setid = count_setid
    end subroutine matching_multiangle

    subroutine get_pairids_of_candi_theta(&
            N_candiset, candi_pairidset, &
            theta_hat, epsilon, N_pairset, thetaset, pairidset)
        implicit none;
        double precision, intent(in) :: theta_hat
        double precision, intent(in) :: epsilon
        integer, intent(in) :: N_pairset
        double precision, intent(in) :: thetaset(N_pairset)
        integer, intent(in) :: pairidset(N_pairset, 2)
        integer, intent(out) :: N_candiset
        integer, intent(out) :: candi_pairidset(N_pairset, 2)
        !
        integer :: ii
        integer :: count
        double precision :: y_lower_bound, y_upper_bound
        !
        y_lower_bound = theta_hat - epsilon
        y_upper_bound = theta_hat + epsilon
        count = 0
        do ii = 1, N_pairset
            if (y_lower_bound < thetaset(ii) .and. thetaset(ii) < y_upper_bound) then
                count = count + 1
                candi_pairidset(count, :) = pairidset(ii, :)
            end if
        end do
        N_candiset = count
    end subroutine get_pairids_of_candi_theta

    subroutine inter_angle_3dvec(interangle, vec1, vec2)
        implicit none;
        double precision, intent(in) :: vec1(3), vec2(3)
        double precision, intent(out) :: interangle
        double precision :: inner, norm_vec1, norm_vec2
        !
        inner = dot_product(vec1, vec2)
        norm_vec1 = norm2(vec1)
        norm_vec2 = norm2(vec2)
        interangle = acos(inner/(norm_vec1*norm_vec2))
    end subroutine inter_angle_3dvec

    subroutine triangle_specular_sign(sign_value, vec1, vec2, vec3)
        implicit none;
        double precision, intent(in) :: vec1(3), vec2(3), vec3(3)
        double precision, intent(out) :: sign_value
        double precision :: cross(3), inner
        ! vec2 \times vec3
        cross(1) = vec2(2)*vec3(3) - vec2(3)*vec3(2)
        cross(2) = vec2(3)*vec3(1) - vec2(1)*vec3(3)
        cross(3) = vec2(1)*vec3(2) - vec2(2)*vec3(1)
        ! vec1 \cdot cross
        inner = vec1(1)*cross(1) + vec1(2)*cross(2) + vec1(3)*cross(3)
        ! sign(inner)
        sign_value = sign(1.0, inner)
    end subroutine triangle_specular_sign

    subroutine equatorial2vec(vec, alpha, delta)
        double precision, intent(in) :: alpha, delta
        double precision, intent(out) :: vec(3)
        !
        vec(1) = cos(alpha) * cos(delta)
        vec(2) = sin(alpha) * cos(delta)
        vec(3) = sin(delta)
    end subroutine equatorial2vec

    subroutine unique_intarray1d(N_array1d_uni, array1d_uni, N_array1d, array1d)
        integer, intent(in) :: N_array1d
        integer, intent(in) :: array1d(N_array1d)
        integer, intent(out) :: N_array1d_uni
        integer, intent(out) :: array1d_uni(N_array1d)
        !
        integer :: array1d_sort(N_array1d)
        integer :: i, count
        !!! sort !!!
        array1d_sort = array1d
        call quicksort(array1d_sort, N_array1d, 1, N_array1d)
        !!! erase duplicate !!!
        ! init
        count = 1
        array1d_uni(count) = array1d_sort(count)
        ! loop
        do i = 2, N_array1d
            if (array1d_sort(i) /= array1d_uni(count)) then
                count = count + 1
                array1d_uni(count) = array1d_sort(i)
            end if
        end do
        N_array1d_uni = count

    contains

        subroutine partition(pivot, array, N, left, right)
            implicit none;
            integer, intent(in) :: N, left, right
            integer, intent(inout) :: array(N)
            integer, intent(out) :: pivot
            ! 
            integer :: i, temp, last
            !
            last = left
            do i = left+1, right
                if (array(i) < array(left)) then
                    last = last + 1
                    ! swap
                    temp = array(last)
                    array(last) = array(i)
                    array(i) = temp
                end if
            end do
            ! swap
            temp = array(left)
            array(left) = array(last)
            array(last) = temp
            pivot = last
        end subroutine partition

        recursive subroutine quicksort(array, N, left, right)
            implicit none;
            integer, intent(in) :: N, left, right
            integer, intent(inout) :: array(N)
            integer :: pivot
            !
            call partition(pivot, array, N, left, right)
            if (left < pivot-1) call quicksort(array, N, left, pivot-1)
            if (pivot+1 < right) call quicksort(array, N, pivot+1, right)
        end subroutine quicksort

    end subroutine unique_intarray1d

end module subgraph_isomorphism_starid_detector



        ! do n = 2, n_obs
        !     if (n == 2) then
        !         allocate(candi_setid_temp(N_pairset, n))
        !         call matching_pair(&
        !             N_candi_setid_temp, candi_setid_temp, &
        !             s_hat_set(1:2, :), epsilon, N_pairset, thetaset, pairidset)
        !     else
        !         deallocate(candi_setid_temp)
        !         allocate(candi_setid_temp(N_pairset, n))
        !         count_setid = 0
        !         ! matching pair stars including last star
        !         allocate(N_candi_pairid_in(n-1))
        !         allocate(candi_pairid_in(n-1, N_pairset, 2))
        !         do i1 = 1, n-1
        !             s_hat_i = s_hat_set(i1, :)
        !             s_hat_n = s_hat_set(n, :)
        !             call inter_angle_3dvec(theta_hat_in, s_hat_i, s_hat_n)
        !             call get_pairids_of_candi_theta(&
        !                 N_candi_pairid_in(i1), candi_pairid_in(i1, :, :), &
        !                 theta_hat_in, epsilon, N_pairset, thetaset, pairidset)
        !         end do
        !         ! determine candidate of last star in each set of matched stars except last
        !         do i2 = 1, pre_N_candi_setid
        !             do i3 = 1, n-1
        !                 candi_starid_i = pre_candi_setid(i2, i3)
        !                 count = 0
        !                 do i4 = 1, N_candi_pairid_in(i3)
        !                     if (candi_pairid_in(i3, i4, 1) == candi_starid_i) then
        !                         count = count + 1
        !                         candi_starid_n_from_i(count) = candi_pairid_in(i3, i4, 2)
        !                     else if (candi_pairid_in(i3, i4, 2) == candi_starid_i) then
        !                          count = count + 1
        !                         candi_starid_n_from_i(count) = candi_pairid_in(i3, i4, 1)
        !                     end if
        !                 end do
        !                 N_candi_starid_n_from_i = count
        !                 !
        !                 if (i3 == 1) then
        !                     N_candi_starid_n = N_candi_starid_n_from_i
        !                     candi_starid_n = candi_starid_n_from_i
        !                 else
        !                     count = 0
        !                     do i5 = 1, N_candi_starid_n
        !                         if (any(candi_starid_n(i5) == candi_starid_n_from_i(1:N_candi_starid_n_from_i))) then
        !                             count = count + 1
        !                             candi_starid_n_temp(count) = candi_starid_n(i5)
        !                         end if
        !                     end do
        !                     N_candi_starid_n = count
        !                     candi_starid_n = candi_starid_n_temp
        !                 end if
        !                 !
        !                 if (N_candi_starid_n == 0) exit
        !             end do
        !             ! 
        !             do i6 = 1, N_candi_starid_n
        !                 count_setid = count_setid + 1
        !                 candi_setid_temp(count_setid, :) = [pre_candi_setid(i6, :), candi_starid_n(i6)]
        !             end do
        !         end do
        !         N_candi_setid_temp = count_setid
        !         !
        !         deallocate(N_candi_pairid_in)
        !         deallocate(candi_pairid_in)
        !         deallocate(pre_candi_setid)
        !     end if
        !     allocate(pre_candi_setid(N_pairset, n))
        !     pre_N_candi_setid = N_candi_setid_temp
        !     pre_candi_setid = candi_setid_temp
        ! end do
        ! N_candi_setid = N_candi_setid_temp
        ! candi_setid = candi_setid_temp




                    ! count_setid = 0
                    ! ! matching pair stars including last star
                    ! allocate(N_candi_pairid_in(n-1))
                    ! allocate(candi_pairid_in(n-1, N_pairset, 2))
                    ! do i1 = 1, n-1
                    !     s_hat_i = s_hat_set(i1, :)
                    !     s_hat_n = s_hat_set(n, :)
                    !     call inter_angle_3dvec(theta_hat_in, s_hat_i, s_hat_n)
                    !     call get_pairids_of_candi_theta(&
                    !         N_candi_pairid_in(i1), candi_pairid_in(i1, :, :), &
                    !         theta_hat_in, epsilon, N_pairset, thetaset, pairidset)
                    ! end do
                    ! ! determine candidate of last star in each set of matched stars except last
                    ! do i2 = 1, pre_N_candi_setid
                    !     do i3 = 1, n-1
                    !         candi_starid_i = pre_candi_setid(i2, i3)
                    !         count = 0
                    !         do i4 = 1, N_candi_pairid_in(i3)
                    !             if (candi_pairid_in(i3, i4, 1) == candi_starid_i) then
                    !                 count = count + 1
                    !                 candi_starid_n_from_i(count) = candi_pairid_in(i3, i4, 2)
                    !             else if (candi_pairid_in(i3, i4, 2) == candi_starid_i) then
                    !                 count = count + 1
                    !                 candi_starid_n_from_i(count) = candi_pairid_in(i3, i4, 1)
                    !             end if
                    !         end do
                    !         N_candi_starid_n_from_i = count
                    !         !
                    !         if (i3 == 1) then
                    !             N_candi_starid_n = N_candi_starid_n_from_i
                    !             candi_starid_n = candi_starid_n_from_i
                    !         else
                    !             count = 0
                    !             do i5 = 1, N_candi_starid_n
                    !                 if (any(candi_starid_n(i5) == candi_starid_n_from_i(1:N_candi_starid_n_from_i))) then
                    !                     count = count + 1
                    !                     candi_starid_n_temp(count) = candi_starid_n(i5)
                    !                 end if
                    !             end do
                    !             N_candi_starid_n = count
                    !             candi_starid_n = candi_starid_n_temp
                    !         end if
                    !         !
                    !         if (N_candi_starid_n == 0) exit
                    !     end do
                    !     ! 
                    !     do i6 = 1, N_candi_starid_n
                    !         count_setid = count_setid + 1
                    !         candi_setid_temp(count_setid, :) = [pre_candi_setid(i6, :), candi_starid_n(i6)]
                    !     end do
                    ! end do
                    ! N_candi_setid_temp = count_setid
                    !
                    ! deallocate(N_candi_pairid_in)
                    ! deallocate(candi_pairid_in)