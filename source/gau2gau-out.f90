program gau2gau_out

    !-------------------------------------------------
    ! This source is part of gau_Subs_gau interface
    !
    ! Description:
    ! ------------
    ! Reads input data file from Gaussian caller and
    ! generates a Gaussian input. 
    !
    ! The progra expects to be called as:
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

    ! System vars
    real(8) :: energy
    integer :: Nat
    ! OUTPUT DATA MATRICES
    real(8),dimension(1:3)             :: dipole
    real(8),dimension(:,:),allocatable :: grad    !(3,Nat)
    real(8),dimension(1:6)             :: polarizability
    real(8),dimension(:,:),allocatable :: Ddipole !(9,Nat)
    real(8),dimension(:),allocatable   :: hessian !((3*Nat*(3*Nat+1))/2)

    ! IO
    character(len=150) :: job1_fchk, job2_fchk, outfile, msgfile
    integer            :: I_FC1=10, &
                          I_FC2=11, &
                          O_OUT=20, &
                          O_MSG=21

    ! Auxiliars
    ! Matrices to read data from fchk
    integer             :: N
    real(8),dimension(:),allocatable :: A
    integer,dimension(:),allocatable :: IA
    integer             :: err
    character           :: dtype
    ! Couters
    integer :: i, j, k


    ! Input data from argument list
    call getarg(1,job1_fchk)
    call getarg(2,job2_fchk)
    call getarg(3,outfile)
    call getarg(4,msgfile)

    ! Open files
    open(I_FC1,file=adjustl(job1_fchk),status='old')
    open(I_FC2,file=adjustl(job2_fchk),status='old')
    open(O_OUT,file=adjustl(outfile)  ,status='replace')
    open(O_MSG,file=adjustl(msgfile)  ,status='replace')
 

    ! READ DATA AND COMPUTE SUBSTRACTIONS
    !0. Get Nat
    call read_fchk(I_FC1,'Number of atoms',dtype,N,A,IA,err)
    if (err==0) then
        Nat = IA(1)
        deallocate(IA)
    else
        write(O_MSG,*) 'ERROR: Number of atoms not available on fchk'
        stop
    endif
    allocate(grad(3,Nat),Ddipole(9,Nat),hessian((3*Nat*(3*Nat+1))/2))

    !1a. Energy
    call read_fchk(I_FC1,'Total Energy',dtype,N,A,IA,err)
    if (err==0) then
        energy = A(1)
        deallocate(A)
    else
        write(O_MSG,*) 'WARNING: Total energy not available for job1. Set to 0.0'
        energy = 0.d0
    endif
    call read_fchk(I_FC2,'Total Energy',dtype,N,A,IA,err)
    if (err==0) then
        energy = energy - A(1)
        deallocate(A)
    else
        write(O_MSG,*) 'WARNING: Total energy not available for job2. Set to 0.0'
        energy = energy - 0.d0
    endif
    !1b. Dipole moment (take the one for job1)
    call read_fchk(I_FC1,'Dipole Moment',dtype,N,A,IA,err)
    if (err==0) then
        dipole(1:3) = A(1:3)
        deallocate(A)
    else
        write(O_MSG,*) 'WARNING: Dipole Moment not available for job1. Set to 0.0'
        dipole(1:3) = 0.d0
    endif
    
    !2. Forces
    call read_fchk(I_FC1,'Cartesian Gradient',dtype,N,A,IA,err)
    if (err==0) then
        do i=1,Nat
            grad(1,i) = A(3*i-2)
            grad(2,i) = A(3*i-1)
            grad(3,i) = A(3*i)
        enddo
        deallocate(A)
    else
        write(O_MSG,*) 'WARNING: Cartesian Gradient not available for job1. Set to 0.0'
        grad = 0.d0
    endif
    call read_fchk(I_FC2,'Cartesian Gradient',dtype,N,A,IA,err)
    if (err==0) then
        do i=1,Nat
            grad(1,i) = grad(1,i) - A(3*i-2)
            grad(2,i) = grad(2,i) - A(3*i-1)
            grad(3,i) = grad(3,i) - A(3*i)
        enddo
        deallocate(A)
    else
        write(O_MSG,*) 'WARNING: Cartesian Gradient not available for job2. Set to 0.0'
        do i=1,Nat
            grad(1,i) = grad(1,i) - 0.d0 
            grad(2,i) = grad(2,i) - 0.d0 
            grad(3,i) = grad(3,i) - 0.d0 
        enddo
    endif
    


    !3. Polizabilities and dipole derivatives (taken from job1)
    call read_fchk(I_FC1,'Polarizability',dtype,N,A,IA,err)
    if (err==0) then
        polarizability(1:6) = A(1:6)
        deallocate(A)
    else
        write(O_MSG,*) 'WARNING: Polarizability not available for job1. Set to 0.0'
        polarizability(1:6) = 0.d0
    endif
    call read_fchk(I_FC1,'Dipole Derivatives',dtype,N,A,IA,err)
    if (err==0) then
        k=0
        do i=1,9
        do j=1,Nat
            k=k+1
            Ddipole(i,j) = A(k)
        enddo
        enddo
        deallocate(A)
    else
        write(O_MSG,*) 'WARNING: Dipole Derivatives not available for job1. Set to 0.0'
        Ddipole(1:6,1:Nat) = 0.d0
    endif


    !4. Hessian 
    call read_fchk(I_FC1,'Cartesian Force Constants',dtype,N,A,IA,err)
    if (err==0) then
        hessian(:) = A(:)
        deallocate(A)
    else
        write(O_MSG,*) 'WARNING: Cartesian Force Constants not available for job1. Set to 0.0'
        hessian(:) = 0.d0
    endif
    call read_fchk(I_FC2,'Cartesian Force Constants',dtype,N,A,IA,err)
    if (err==0) then
        do i=1,N
            hessian(i) = hessian(i) - A(i)
        enddo
        deallocate(A)
    else
        write(O_MSG,*) 'WARNING: Cartesian Force Constants not available for job2. Set to 0.0'
        do i=1,3*Nat*(3*Nat+1)/2
            hessian(i) = hessian(i) - 0.d0
        enddo
    endif

    ! Close input files
    close(I_FC1)
    close(I_FC2)



    ! WRITE DATA
    !1. Energy and dipole moment
    write(O_OUT,'(4D20.12)') energy, dipole(1:3)

    !2. Forces
    do i=1,Nat
        write(O_OUT,'(3D20.12)') grad(1:3,i)
    enddo

    !3. Polizabilities and dipole derivatives
    write(O_OUT,'(3D20.12)') polarizability(1:3)
    write(O_OUT,'(3D20.12)') polarizability(4:6)
    do i=1,Nat
        write(O_OUT,'(3D20.12)') Ddipole(1:3,i)
        write(O_OUT,'(3D20.12)') Ddipole(4:6,i)
        write(O_OUT,'(3D20.12)') Ddipole(7:9,i)
    enddo        

    !4. Hessian
    do i=1,3*Nat*(3*Nat+1)/2,3
        write(O_OUT,'(3D20.12)') hessian(i:i+2)
    enddo

    close(O_OUT)


    !WRITE MESSAGE
    write(O_msg,'(A,/)') 'gauSgau  The output file for Gaussian has been updated'

    close(O_msg)

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


