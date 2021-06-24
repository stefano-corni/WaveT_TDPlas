module dictionary_entry
    use constants

    implicit none

    type dict_entry
         integer(i4b)   :: n_values
         character(flg), allocatable :: keys(:)
         character(flg), allocatable  :: user_values(:)
         character(flg), allocatable  :: internal_values(:)
    end type


    public dict_entry, dict_entry_init, d_entry_convert_to_internal, d_entry_check_allowed_val


    contains


    subroutine dict_entry_init(d_entry, n_values, keys, user_values, internal_values)
        integer :: n_values
        character(flg), dimension(2), intent(in) :: keys
        character(flg), dimension(n_values), intent(in) :: user_values
        character(flg), dimension(n_values), intent(in) :: internal_values
        type(dict_entry), intent(inout) :: d_entry
        allocate(d_entry%keys(2))
        allocate(d_entry%user_values(n_values))
        allocate(d_entry%internal_values(n_values))
        d_entry%n_values = n_values
        d_entry%keys = keys
        d_entry%user_values = user_values
        d_entry%internal_values = internal_values
    end subroutine


    subroutine d_entry_check_allowed_val(d_entry, chosen_value)
        type(dict_entry),intent(in) :: d_entry
        character(flg) :: chosen_value
        if(any(chosen_value.eq.d_entry%user_values)) then
        else
            write(*,*), "ERROR: value for keyword:", d_entry%keys(1), "is wrong."
            stop
        end if
    end subroutine


    character(flg) function d_entry_convert_to_internal(d_entry, user_value) result(internal_value)
        type(dict_entry), intent(in) :: d_entry
        character(flg), intent(in)  :: user_value
        integer index_value
        call my_findloc(d_entry%user_values, d_entry%n_values, user_value, index_value)
        internal_value = d_entry%internal_values(index_value)
    end function


    subroutine my_findloc(array, array_size, array_value, array_index)
        integer, intent(in)  :: array_size
        character(*), dimension(array_size), intent(in) :: array
        character(*), intent(in)  :: array_value
        integer, intent(out) :: array_index
        integer ::  i
        do i=1,array_size
           if(array(i).eq.array_value) then
               array_index = i
           endif
        enddo
    end subroutine


end module
