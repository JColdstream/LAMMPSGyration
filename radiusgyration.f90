program radiusgyration
use omp_lib

implicit none
include 'radiusgyration.inc'
integer(8):: iter

frame = 1
call start
open(unit = 110, file='test.xyz', status='unknown')

do while ( trajstatus .eq. 0 .and. (lastframe .eq. 0 .or. frame .le. lastframe) )
  if (mod(frame, 50) .eq. 0) write(*,*) frame
  ! write(*, *) frame
  call readheader
  call readcoordinates
  if ( trajstatus .ne. 0 ) cycle
  call initialcom
  call iteratecom
  call comrg
  write (110, *) nanalyse
  write (110, *)
  do iter = 1, nanalyse
    write(110, *) 'C', atom(:, iter)
  enddo
  frame = frame+1
enddo

! test .xyz file to check molecule rebuilding
! open(unit = 110, file='test.xyz', status='unknown')
! do iter = 1, nanalyse
!   write(110, *) 'C', atom(:, iter)
! enddo

write(*, *) 'Final frame: ', frame-1

contains

subroutine start

call get_command_argument(1, inputfile)
if (len_trim(inputfile) == 0) then
  write(*,*) 'Provide an input file.'
  call exit
endif

call readinput

! open lammps trajectory file
open(11, file = trajfile, status='old')
! get ntotal
call readheader

! allocate arrays
call readmasses
call nanalyse_setup
allocate(atom(3, nanalyse))
allocate(atomtype(nanalyse))
rewind(11)

call skipframe
frame = frame+nskip

open(unit = 102, file='rg.out', status='unknown')
open(unit = 103, file='anisotropy.out', status='unknown')
open(unit = 104, file='abs_anisotropy.out', status='unknown')

end subroutine start

! reads analysis.input input file
subroutine readinput
open(10, file = inputfile, status='old')
read(10, *)
read(10, *)
read(10, '(a)') trajfile
read(10, *) massfile
read(10, *) nskip
read(10, *) lastframe
read(10, *) ntypes
read(10, '(a)') strtypesanalyse
typesanalyse = string_to_integers(strtypesanalyse, ",")
write(*,*) trajfile

end subroutine readinput

subroutine readmasses
  integer(sp):: i, j

  allocate(mass(ntypes))

  open(massf, file = massfile, status='old')
  read(massf, *)
  read(massf, *)
  do i = 1, ntypes
    read(massf, *) j, mass(i)
    ! write(*,*) j, mass(i)
  enddo
  close(massf)
  write(*,*) 'Masses read.'
end subroutine readmasses

! skips frames that don't need to be analysed
subroutine skipframe
integer:: i
do i = 1, nskip*(ntotal+9)
  read(11, *)
enddo
write(*,*) "SKIPPED FRAMES : ", nskip
write(*,*)
end subroutine skipframe

! reads in box length, timestep and natoms from the LAMMPS trajectory headers
subroutine readheader
read(11, *, iostat = trajstatus) headertext(1)   
read(11, *, iostat = trajstatus) timestep
read(11, *, iostat = trajstatus) headertext(2)
read(11, *, iostat = trajstatus) ntotal
read(11, *, iostat = trajstatus) headertext(3)
read(11, *, iostat = trajstatus) xmin, xmax
read(11, *, iostat = trajstatus) ymin, ymax
read(11, *, iostat = trajstatus) zmin, zmax
read(11, *, iostat = trajstatus) headertext(4)

!if(natom  .ne. ntotal) stop 'mismatched number of atoms'

lx = xmax-xmin
ly = ymax-ymin
lz = zmax-zmin

! Edit pbc conditions if it needs to be used for a non-cubic box.
if ((lx .ne. ly) .or. (ly .ne. lz)) then
  write(*, *) 'Script is not suitable for non-cubic boxes.'
  call exit
endif

!write(*,*) timestep
!write(*,*) xmin, xmax
!write(*,*) ymin, ymax
!write(*,*) zmin, zmax

end subroutine readheader

subroutine nanalyse_setup
  integer(sp):: i, atom_type, nindex
  nanalyse = 0
  masstotal = 0.0_dp
  do i = 1, ntotal
    read(11, *) nindex, atom_type
    if (atom_type .gt. ntypes) then
      write(*,*) "Number of atom types mismatch." 
      call exit
    endif
    if ( any(typesanalyse .eq. atom_type) ) then
      masstotal = masstotal+mass(atom_type)
      nanalyse = nanalyse+1
    endif
  enddo
  ! Array to keep track of which atoms we need to calculate.
  allocate(molanalyse(nanalyse))
end subroutine nanalyse_setup

subroutine readcoordinates
  integer(sp):: i, j, nindex, temp_atom_type
  real(dp):: tempx, tempy, tempz
  j = 1
  do i = 1, ntotal
    read(11, *, iostat = trajstatus)  nindex, temp_atom_type, tempx, tempy, tempz

    if (any(typesanalyse .eq. temp_atom_type)) then
      atom(:, j) = (/tempx, tempy, tempz/)
      atomtype(j) = temp_atom_type
      j = j+1
    endif
  enddo

 ! zero coordinates
 atom(1, :) = atom(1, :)-xmin
 atom(2, :) = atom(2, :)-ymin
 atom(3, :) = atom(3, :)-zmin
end subroutine readcoordinates

subroutine initialcom
  integer(8):: mol, i

  com = 0.0_dp
  do mol = 1, nanalyse
    do i = 1, 3
      com(i) = com(i) + atom(i, mol)*mass(atomtype(mol))
    enddo
  enddo
  com = com/masstotal
end subroutine initialcom

subroutine iteratecom
  integer(8):: i, mol
  real(dp), dimension(3):: temp_com, dxyz, com_diff
  real(dp):: com_diff_mag, drcom 
  logical:: final_com

  final_com = .False.

  do while (final_com .eqv. .False.)
    temp_com = com
    com = 0.0_dp
    do mol = 1, nanalyse
      do i = 1, 3
        ! dxyz(i) = dxyz(i) + (lx*anint(dxyz(i)/lx))
        com(i) = com(i) + (atom(i, mol) - lx*anint((atom(i, mol)-temp_com(i))/lx)) * mass(atomtype(mol))
        ! com(i) = com(i) + atom(i, mol) - lx*anint((atom(i, mol)-temp_com(i))/lx)
      enddo
    enddo
    com = com/masstotal
    !com = com/nanalyse
    do i = 1, 3
      if (com(i) .gt. lx) then
        com(i) = com(i) - lx
      elseif (com(i) .lt. 0.0) then
        com(i) = com(i) + lx
      endif
    enddo
    ! write(*,*) com
    com_diff = com-temp_com
    ! check how close the new centre of geometry is to the previous iteration
    com_diff_mag = 0.0_dp
    do i = 1, 3
      com_diff_mag = com_diff_mag+com_diff(i)**2
    enddo
    if (com_diff_mag .eq. 0.0_dp) then
      final_com = .True.
    endif
  enddo
  ! write(*,*) com
end subroutine iteratecom

subroutine comrg

  integer(8):: i, j, mol, moli, molj, minaxis, maxaxis
  real(dp):: tempq, drsq, dr, drcom, dotprod
  ! ellipsoid semiaxes
  real(dp):: ellipseaxes(3)
  real(dp):: a, b, c
  real(dp):: testrg

  double precision, dimension(3):: dxyz, ev

  !!! LAPACK CODE FOR DSYEV ADAPTED FROM !!!
  !!! https://www.intel.com/content/www/us/en/docs/onemkl/code-samples-lapack/2022-1/dsyev-example-fortran.html !!!

  ! LAPACK ARGUMENTS

  integer:: LWORK, INFO
  ! 1000 is max length of work array
  double precision:: W(3), WORK(1000) 
  ! END OF LAPACK ARGUMENTS

  rg = 0.0
  do mol = 1, nanalyse

    ! converts atom coordinates to the minimum image vector to centre of mass
    do i = 1, 3
      drcom = atom(i, mol) - com(i)
      if ( abs(drcom) .gt. lx/2 ) then
        atom(i, mol) = atom(i, mol) - lx*anint(drcom/lx)
      endif
      ! set origin to the centre of mass then calculate the rg
      atom(i, mol) = atom(i, mol) - com(i)
      rg = rg + ( mass(atomtype(mol))*( atom(i, mol) )**2 )
    enddo

    ! calculates the inertia tensor using the dot product, and the outer product
    dotprod = dot_product( atom(:,mol), atom(:,mol) )

    do i = 1, 3
      do j = 1, 3
        if (i .eq. j) then
          itensor(j, i) = itensor(j, i) + mass(atomtype(mol)) * ( dotprod-atom(i, mol)*atom(j, mol) ) 
        else
          itensor(j, i) = itensor(j, i) - mass(atomtype(mol)) * ( atom(i, mol)*atom(j, mol) ) 
        endif
      enddo
    enddo

  enddo

  !write(*,*) itensor

  ! find optimal workspace
  LWORK = -1
  CALL DSYEV('Vectors', 'Upper', 3, itensor, 3, W, WORK, LWORK, INFO)
  LWORK = MIN( 1000, INT( WORK( 1 ) ) )

  ! solve eigenvalue problem
  CALL DSYEV('Vectors', 'Upper', 3, itensor, 3, W, WORK, LWORK, INFO)
  
  !write(*,*) "Eigenvalues done."
  !write(*,*) W

  a = sqrt( 5*( W(2) + W(3) - W(1) ) / ( 2*masstotal ) )
  b = sqrt( 5*( W(3) + W(1) - W(2) ) / ( 2*masstotal ) )
  c = sqrt( 5*( W(1) + W(2) - W(3) ) / ( 2*masstotal ) )

  ellipseaxes(1) = a 
  ellipseaxes(2) = b
  ellipseaxes(3) = c

  testrg = sqrt( ( a**2+b**2+c**2 ) / 5 )

  ! convert to ratios instead of absolute values
  minaxis = minloc(ellipseaxes, 1)
  maxaxis = maxloc(ellipseaxes, 1)
  write(*,*) minaxis, maxaxis
  ellipseaxes = ellipseaxes/ellipseaxes(minaxis)

  if ( minaxis .eq. 1 ) then
    if ( maxaxis .eq. 2 ) then
      write(103, *) timestep, ellipseaxes(2), ellipseaxes(3)
    elseif ( maxaxis .eq. 3 ) then
      write(103, *) timestep, ellipseaxes(3), ellipseaxes(2)
    endif

  elseif ( minaxis .eq. 2 ) then
    if ( maxaxis .eq. 1 ) then
      write(103, *) timestep, ellipseaxes(1), ellipseaxes(3)
    elseif ( maxaxis .eq. 3 ) then
      write(103, *) timestep, ellipseaxes(3), ellipseaxes(1)
    endif

  elseif ( minaxis .eq. 3 ) then
    if ( maxaxis .eq. 1 ) then
      write(103, *) timestep, ellipseaxes(1), ellipseaxes(2)
    elseif ( maxaxis .eq. 2 ) then
      write(103, *) timestep, ellipseaxes(2), ellipseaxes(1)
    endif
  endif
      

  write(102, *) timestep, testrg
  write(104, *) timestep, a, b, c

  rgnrm = rgnrm+1.0_dp
          
end subroutine comrg

function string_to_integers(str, sep) result(a)
    integer, allocatable:: a(:)
    character(*):: str
    character:: sep
    integer:: i, n_sep

    n_sep = 0
    do i = 2, len_trim(str)
      if (str(i:i)==sep .and. str(i-1:i-1)/=sep) then
        n_sep = n_sep+1
        str(i:i) = ','
       end if
    end do
    allocate(a(n_sep+1))
    read(str, *) a
end function

end program radiusgyration
