program gau2gau_out

    !-------------------------------------------------
    ! This source is part of gau_Subs_gau interface
    !
    ! Description:
    ! ------------
    ! Reads input data file from Gaussian caller and
    ! generates a Gaussian input. 
    !
    ! The program expects to be called as:
    !
    !  gau2gau_out <job1.fchk> <job2.fchk> <output_data> <message_file>
    !
    ! The output data filie must contain 
    ! (see https://gaussian.com/external/):
    ! Items                       Pseudo Code                           Line Format
    ! energy, dipole-moment (xyz) E, Dip(I), I=1,3                      4D20.12
    ! gradient on atom (xyz)      FX(J,I), J=1,3; I=1,NAtoms            3D20.12
    ! polarizability              Polar(I), I=1,6                       3D20.12
    ! dipole derivatives          DDip(I), I=1,9*NAtoms                 3D20.12
    ! force constants             FFX(I), I=1,(3*NAtoms*(3*NAtoms+1))/2 3D20.12
    !--------------------------------------------------------------------------------

    implicit none

    stop

    contains

    subroutine read_fchk(unt,section,data_type,N,A,I,error_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS 
        !==============================================================
        !Description
        ! Generic SR to read any section of the checkpoint
        ! Enter deallocated arrays
        !Arguments
        ! unt (int;in): unit number of the log file
        ! section(char,in): name of the section to be read
        ! data_type(char,out): Integer (I) or Real (R) data read
        ! N(int,in): Number of elements to be read
        ! A(real,dimension(:)): Real array to store real data
        ! I(integer,dimension(:)): Int array to store int data
        ! error_flag(integer,out): 0: success
        !                          1: section not found
        !==============================================================

        integer,intent(in)                                       :: unt
        character(len=*),intent(in)                              :: section
        character(len=1),intent(out)                             :: data_type
        integer,intent(out)                                      :: N
        double precision, dimension(:), allocatable, intent(out) :: A
        integer,dimension(:), allocatable, intent(out)           :: I
        integer,intent(out),optional                             :: error_flag

        !Local stuff
        !=============
        character(len=240) :: line=""
        character(len=42)  :: section_full
        character(len=1)   :: is_array
        character(len=40)  :: cdata
        !I/O
        integer :: IOstatus
        
        
        ! Search section
        if (present(error_flag)) error_flag = 0
        do
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) then
                    if (present(error_flag)) error_flag=1
                    rewind(unt)
                    return
                endif
                ! 2) Found what looked for!      
                if ( INDEX(line,trim(adjustl(section))) /= 0 ) then
                    read(line,'(A42)') section_full
                    if (adjustl(section_full) == adjustl(section)) exit
                endif
        enddo

        !Get info about section from line just read
        read(line,'(A42,X,A1,3X,A1,X,A)') section_full, data_type, is_array, cdata
        if (is_array /= "N") then
            !Is not an array
            N=1
            if ( data_type == "R" ) then 
                allocate( A(1:1) )
                read(cdata,*) A(1)
            elseif ( data_type == "I" ) then
                allocate( I(1:1) )
                read(cdata,*) I(1) 
            endif
        else
            read(cdata,*) N
            if ( data_type == "R" ) then
                allocate( A(1:N) )
                read(unt,*) A(1:N)
            elseif ( data_type == "I" ) then
                allocate( I(1:N) )
                read(unt,*) I(1:N)
            endif
        endif 

        rewind(unt)
        return

    end subroutine read_fchk

end program gau2gau_out


