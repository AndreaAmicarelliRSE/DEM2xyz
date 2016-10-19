!-------------------------------------------------------------------------------
! "DEM2xyz v.1.0" 
! Copyright 2016 (RSE SpA)
! "DEM2xyz v.1.0" authors and email contact are provided on 
! the documentation file.
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Description. “DEM2xyz v.1.0” (RSE SpA) reads a “DEM” file and writes the 
!              associated DEM (or DTM) in a corresponding “xyz” file, possibly 
!              changing the spatial resolution (as requested by the user).
!-------------------------------------------------------------------------------
PROGRAM DEM2xyz
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
double precision,dimension(:,:),allocatable :: mat_z
double precision :: dx   
integer :: i,j,n_col,n_raw,res_fact,n_points
character(len=100) :: char_aux
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
! Subroutine info
print *,"!!!DEM2xyz is running:"
print *,"This subroutine reads a DEM file"
print *,"    and writes a corresponding xyz file." 
print *,"Algorithm:"
print *,"   1) Reading DEM file "
print *,"   2) Writing xyz file "
! 1) Reading DEM (start)
print *,"1)  Reading DEM file "
open (1,file='DEM.txt') 
read (1,'(a14,i15)') char_aux,n_col
read (1,'(a14,i15)') char_aux,n_raw
read (1,'')
read (1,'')
read (1,'(a14,f15.7,i15)') char_aux,dx,res_fact
read (1,'')
allocate(mat_z(n_raw,n_col))
mat_z = 0.d0
do i=1,n_raw
   read (1,*) mat_z(i,:)
enddo
close(1)
print *,"End Reading DEM file "
! End 1) Reading DEM file (end)
! 2) Writing xyz file (start)
print *,"2)  Writing xyz file "
open (2,file="xyz.txt")
write (2,'(a)') '        X(m)       Y(m)        z(m)'
n_points = n_raw*n_col/res_fact/res_fact
write (2,'(i12)') n_points
do i=1,n_raw,res_fact
   do j=1,n_col,res_fact
      write (2,'(3(F12.4))') (j-1)*dx+dx/2,(n_raw+1-i)*dx-dx/2,mat_z(i,j) 
   end do   
end do
close (2)
print *,"End Writing xyz file "
! End 2) Writing xyz file (end)
!------------------------
! Deallocations
!------------------------
deallocate(mat_z)
print *,"!!!  DEM2xyz has terminated "
end program DEM2xyz
