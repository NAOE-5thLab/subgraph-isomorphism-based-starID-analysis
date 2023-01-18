module subgraph_isomorphism_starid_detector
    implicit none;

contains

    subroutine matching_set_for_analysis(&
            N_candi_setid, candi_setid, time, &
            N_candi, n_obs, s_hat_set, epsilon, &
            N_set, RA_set, DE_set, N_pairset, thetaset, pairidset)
        implicit none;
        integer, intent(in) :: n_obs
        double precision, intent(in) :: s_hat_set(n_obs, 3)
        double precision, intent(in) :: epsilon
        integer, intent(in) :: N_set, N_pairset, N_candi
        double precision, intent(in) :: RA_set(N_set), DE_set(N_set)
        double precision, intent(in) :: thetaset(N_pairset)
        integer, intent(in) :: pairidset(N_pairset, 2)
        integer, intent(out) :: N_candi_setid(n_obs-1)
        integer, intent(out) :: candi_setid(n_obs-1, N_candi, n_obs)
        double precision, intent(out) :: time(n_obs)
        ! 
        integer :: n
        ! 
        do n = 2, n_obs
            call cpu_time(time(n-1))
            if (n == 2) then
                call matching_pair(&
                    N_candi_setid(n - 1), candi_setid(n - 1, :, :n), &
                    N_candi, s_hat_set(:n, :), epsilon, N_pairset, thetaset, pairidset)
            else if (n == 3) then
                call matching_triangle(&
                    N_candi_setid(n - 1), candi_setid(n - 1, :, :n), &
                    N_candi, s_hat_set(:n, :), epsilon, &
                    N_set, RA_set, DE_set, N_pairset, thetaset, pairidset)
            else
                call matching_multiangle(&
                    N_candi_setid(n - 1), candi_setid(n - 1, :, :n), &
                    N_candi, n, N_candi_setid(n - 2), candi_setid(n - 2, :, :n - 1), &
                    s_hat_set(:n, :), epsilon, N_pairset, thetaset, pairidset)
            end if
        end do
        call cpu_time(time(n_obs))
    end subroutine matching_set_for_analysis

    subroutine matching_pair(&
            N_candi_pairid_ij, candi_pairid_ij, &
            N_candi, s_hat_set, epsilon, N_pairset, thetaset, pairidset)
        implicit none;
        double precision, intent(in) :: s_hat_set(2, 3)
        double precision, intent(in) :: epsilon
        integer, intent(in) :: N_pairset, N_candi
        double precision, intent(in) :: thetaset(N_pairset)
        integer, intent(in) :: pairidset(N_pairset, 2)
        integer, intent(out) :: N_candi_pairid_ij
        integer, intent(out) :: candi_pairid_ij(N_candi, 2)
        !
        double precision :: s_hat_i(3), s_hat_j(3), theta_hat_ij
        !
        s_hat_i = s_hat_set(1, :)
        s_hat_j = s_hat_set(2, :)
        call inter_angle_3dvec(theta_hat_ij, s_hat_i, s_hat_j)
        call get_pairids_of_candi_theta(&
            N_candi_pairid_ij, candi_pairid_ij, N_candi, &
            theta_hat_ij, epsilon, &
            N_pairset, thetaset, pairidset)
    end subroutine matching_pair

    subroutine matching_triangle(&
            N_candi_triangleid_ijk, candi_triangleid_ijk, &
            N_candi, s_hat_set, epsilon, &
            N_set, RA_set, DE_set, N_pairset, thetaset, pairidset)
        implicit none;
        double precision, intent(in) :: s_hat_set(3, 3)
        double precision, intent(in) :: epsilon
        integer, intent(in) :: N_set, N_pairset, N_candi
        double precision, intent(in) :: RA_set(N_set), DE_set(N_set)
        double precision, intent(in) :: thetaset(N_pairset)
        integer, intent(in) :: pairidset(N_pairset, 2)
        integer, intent(out) :: N_candi_triangleid_ijk
        integer, intent(out) :: candi_triangleid_ijk(N_candi, 3)
        !
        double precision :: s_hat_i(3), s_hat_j(3), s_hat_k(3)
        double precision :: theta_hat_ij, theta_hat_ik, theta_hat_jk
        integer :: N_candi_triangleid_ijk_temp, candi_triangleid_ijk_temp(N_candi, 3)
        integer :: N_candi_pairid_ij, candi_pairid_ij(N_candi, 2)
        integer :: N_candi_pairid_ik, candi_pairid_ik(N_candi, 2)
        integer :: N_candi_pairid_jk, candi_pairid_jk(N_candi, 2)
        double precision :: obs_sign, candi_sign
        double precision :: vec_i(3), vec_j(3), vec_k(3)
        integer :: i1, i2, i3, i4
        integer :: count_tiangleid
        integer :: candi_i
        integer :: candi_j(N_set)
        integer :: candi_k(N_set)
        integer :: stage(N_set)
        !
        s_hat_i = s_hat_set(1, :)
        s_hat_j = s_hat_set(2, :)
        s_hat_k = s_hat_set(3, :)
        call inter_angle_3dvec(theta_hat_ij, s_hat_i, s_hat_j)
        call inter_angle_3dvec(theta_hat_ik, s_hat_i, s_hat_k)
        call inter_angle_3dvec(theta_hat_jk, s_hat_j, s_hat_k)
        call get_pairids_of_candi_theta(&
            N_candi_pairid_ij, candi_pairid_ij, N_candi, &
            theta_hat_ij, epsilon, N_pairset, thetaset, pairidset)
        call get_pairids_of_candi_theta(&
            N_candi_pairid_ik, candi_pairid_ik, N_candi, &
            theta_hat_ik, epsilon, N_pairset, thetaset, pairidset)
        call get_pairids_of_candi_theta(&
            N_candi_pairid_jk, candi_pairid_jk, N_candi, &
            theta_hat_jk, epsilon, N_pairset, thetaset, pairidset)
        !
        candi_j = 0
        candi_k = 0
        stage = 0
        do i1 = 1, N_candi_pairid_ij
            candi_i = candi_pairid_ij(i1, 1)
            if (stage(candi_i) == 0) then
                stage(candi_i) = 1
                candi_j(candi_i) = candi_pairid_ij(i1, 2)
            end if
            candi_i = candi_pairid_ij(i1, 2)
            if (stage(candi_i) == 0) then
                stage(candi_i) = 1
                candi_j(candi_i) = candi_pairid_ij(i1, 1)
            end if
        end do
        do i2 = 1, N_candi_pairid_ik
            candi_i = candi_pairid_ik(i2, 1)
            if (stage(candi_i) == 1) then
                stage(candi_i) = 2
                candi_k(candi_i) = candi_pairid_ik(i2, 2)
            end if
            candi_i = candi_pairid_ik(i2, 2)
            if (stage(candi_i) == 1) then
                stage(candi_i) = 2
                candi_k(candi_i) = candi_pairid_ik(i2, 1)
            end if
        end do
        count_tiangleid = 0
        do candi_i = 1, N_set
            if (stage(candi_i) == 2) then
                do i3 = 1, N_candi_pairid_jk
                    if (candi_j(candi_i) == candi_pairid_jk(i3, 1)) then
                        if (candi_k(candi_i) == candi_pairid_jk(i3, 2)) then
                            stage(candi_i) = 3
                            count_tiangleid = count_tiangleid + 1
                            candi_triangleid_ijk_temp(count_tiangleid, 1) = candi_i
                            candi_triangleid_ijk_temp(count_tiangleid, 2) = candi_j(candi_i)
                            candi_triangleid_ijk_temp(count_tiangleid, 3) = candi_k(candi_i)
                        end if
                    end if
                    if (candi_j(candi_i) == candi_pairid_jk(i3, 2)) then
                        if (candi_k(candi_i) == candi_pairid_jk(i3, 1)) then
                            stage(candi_i) = 3
                            count_tiangleid = count_tiangleid + 1
                            candi_triangleid_ijk_temp(count_tiangleid, 1) = candi_i
                            candi_triangleid_ijk_temp(count_tiangleid, 2) = candi_j(candi_i)
                            candi_triangleid_ijk_temp(count_tiangleid, 3) = candi_k(candi_i)
                        end if
                    end if
                end do
            end if
        end do
        N_candi_triangleid_ijk_temp = count_tiangleid
        ! check specular sign
        count_tiangleid = 0
        call triangle_specular_sign(obs_sign, s_hat_i, s_hat_j, s_hat_k)
        do i4 = 1, N_candi_triangleid_ijk_temp
            call equatorial2vec(vec_i, &
                    RA_set(candi_triangleid_ijk_temp(i4, 1) + 1), &
                    DE_set(candi_triangleid_ijk_temp(i4, 1) + 1))
            call equatorial2vec(vec_j, &
                    RA_set(candi_triangleid_ijk_temp(i4, 2) + 1), &
                    DE_set(candi_triangleid_ijk_temp(i4, 2) + 1))
            call equatorial2vec(vec_k, &
                    RA_set(candi_triangleid_ijk_temp(i4, 3) + 1), &
                    DE_set(candi_triangleid_ijk_temp(i4, 3) + 1))
            ! 
            call triangle_specular_sign(candi_sign, vec_i, vec_j, vec_k)
            if (obs_sign*candi_sign > 0) then
                count_tiangleid = count_tiangleid + 1
                candi_triangleid_ijk(count_tiangleid, :) = candi_triangleid_ijk_temp(i4, :)
            end if
        end do
        N_candi_triangleid_ijk = count_tiangleid
    end subroutine matching_triangle

    subroutine matching_multiangle(&
            N_candi_setid, candi_setid, &
            N_candi, n_obs, pre_N_candi_setid, pre_candi_setid, &
            s_hat_set, epsilon, N_pairset, thetaset, pairidset)
        implicit none;
        integer, intent(in) :: n_obs
        integer, intent(in) :: pre_N_candi_setid
        integer, intent(in) :: pre_candi_setid(N_candi, n_obs-1)
        double precision, intent(in) :: s_hat_set(n_obs, 3)
        double precision, intent(in) :: epsilon
        integer, intent(in) :: N_pairset, N_candi
        double precision, intent(in) :: thetaset(N_pairset)
        integer, intent(in) :: pairidset(N_pairset, 2)
        integer, intent(out) :: N_candi_setid
        integer, intent(out) :: candi_setid(N_candi, n_obs)
        !
        integer :: N_candi_pairid_in(n_obs-1), candi_pairid_in(n_obs-1, N_candi, 2)
        double precision :: s_hat_i(3), s_hat_n(3), theta_hat_in
        integer :: candi_starid_i
        integer :: N_candi_starid_n_from_i, candi_starid_n_from_i(N_candi)
        integer :: N_candi_starid_n_temp, candi_starid_n_temp(N_candi)
        integer :: N_candi_starid_n, candi_starid_n(N_candi)
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
                N_candi_pairid_in(i1), candi_pairid_in(i1, :, :), N_candi, &
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
            N_candi, theta_hat, epsilon, N_pairset, thetaset, pairidset)
        implicit none;
        double precision, intent(in) :: theta_hat
        double precision, intent(in) :: epsilon
        integer, intent(in) :: N_pairset
        integer, intent(in) :: N_candi
        double precision, intent(in) :: thetaset(N_pairset)
        integer, intent(in) :: pairidset(N_pairset, 2)
        integer, intent(out) :: N_candiset
        integer, intent(out) :: candi_pairidset(N_candi, 2)
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
        sign_value = sign(1.0d0, inner)
    end subroutine triangle_specular_sign

    subroutine equatorial2vec(vec, alpha, delta)
        implicit none;
        double precision, intent(in) :: alpha, delta
        double precision, intent(out) :: vec(3)
        !
        vec(1) = cos(alpha) * cos(delta)
        vec(2) = sin(alpha) * cos(delta)
        vec(3) = sin(delta)
    end subroutine equatorial2vec

    subroutine unique_intarray1d(N_array1d_uni, array1d_uni, N_array1d, array1d)
        implicit none;
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
