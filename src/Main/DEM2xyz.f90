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
!              changing the spatial resolution (as requested by the user). In 
!              case the absolute value of the mean latitude is provided with
!              a non-negative value, the following conversion takes place 
!              "(lon,lat) in (°) to (X,Y) in (m)". In this case, an 
!              interpolation (weighted on the square of the distance) is 
!              carried out to provide a uniform Cartesian output grid in (X,Y).
!-------------------------------------------------------------------------------
PROGRAM DEM2xyz
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
integer :: i_in,j_in,i_out,j_out,n_col_in,n_col_out,n_raw,res_fact,n_points_in
integer :: n_points_out
double precision :: dx,dy,abs_mean_latitude,denom,distance,x_in,x_out,y_in,y_out
character(len=100) :: char_aux
double precision,dimension(:,:),allocatable :: mat_z_in,mat_z_out
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
write(*,*) "***DEM2xyz is running:"
write(*,*) "This subroutine reads a DEM file"
write(*,*) "    and writes a corresponding xyz file." 
write(*,*) "Algorithm:"
write(*,*) "   1) Reading DEM file and pre-processing "
write(*,*) "   2) Eventual grid interpolation and writing xyz file "
write(*,*) "1)  Reading DEM file and pre-processing "
open(1,file='DEM.txt') 
read(1,'(a14,i15)') char_aux,n_col_in
read(1,'(a14,i15)') char_aux,n_raw
read(1,'(a)')
read(1,'(a)')
read(1,'(a14,f15.7,i15,f15.7)') char_aux,dy,res_fact,abs_mean_latitude
if (abs_mean_latitude>=0.d0) then
abs_mean_latitude = abs_mean_latitude / 180.d0 * 3.1415926
! Conversion (lon,lat) in (°) to (X,Y) in (m) for dx and dy
! Linear unit discretization along the same parallel/latitude due to changing 
! in longitude according to the DEM input discretization 
   dx = dy * (111412.84d0 * dcos(abs_mean_latitude) - 93.5d0 * dcos(3.d0 *     &
        abs_mean_latitude) + 0.118d0 * dcos(5.d0 * abs_mean_latitude))
! Linear unit discretization along the same meridian/longitude due to changing 
! in latitude according to the DEM input discretization 
   dy = dy * (111132.92d0 - 559.82d0 * dcos(2.d0 * abs_mean_latitude) +        &
        1.175d0 * dcos(4.d0 * abs_mean_latitude) - 0.0023d0 * dcos(6.d0 *      &
        abs_mean_latitude))
   n_col_out = floor(n_col_in * dx / dy)
   else
      dx = dy
      n_col_out = n_col_in
endif
read(1,'(a)')
allocate(mat_z_in(n_raw,n_col_in))
allocate(mat_z_out(n_raw,n_col_out))
mat_z_in = 0.d0
mat_z_out = 0.d0
do i_in=1,n_raw
   read(1,*) mat_z_in(i_in,:)
enddo
close(1)
write(*,*) "End Reading DEM file and pre-processing "
write(*,*) "2)  Eventual grid interpolation and writing xyz file "
open(2,file="xyz.txt")
write(2,'(a)') '        x(m)       y(m)        z(m)        z   '
n_points_in = n_raw * n_col_in / res_fact / res_fact
n_points_out = n_raw * n_col_out / res_fact / res_fact
write(*,'(a,i15)') 'Number of vertices in the input "DEM" file: ',n_points_in
write(*,'(a,i15)') 'Number of vertices in the output "xyz" file: ',n_points_out
do j_out=1,n_col_out,res_fact
   do i_out=1,n_raw,res_fact
      x_out = (j_out - 1) * dy + dy / 2.d0
      y_out = (n_raw + 1 - i_out) * dy - dy / 2.d0
      if (abs_mean_latitude>=0.d0) then      
! Interpolation: inverse of the distance**2
         denom = 0.d0
         do j_in=1,n_col_in
            do i_in=1,n_raw
               x_in = (j_in - 1) * dx + dx / 2.d0
               y_in = (n_raw + 1 - i_in) * dy - dy / 2.d0
               distance = dsqrt((x_in - x_out) ** 2 + (y_in - y_out) ** 2)
               if (distance<=(dx*10.d0)) then
                  mat_z_out(i_out,j_out) = mat_z_out(i_out,j_out) +            &
                                           mat_z_in(i_in,j_in) / distance ** 2
                  denom = denom + 1.d0 / distance ** 2
               endif
            enddo
         enddo
         if (denom/=0.d0) mat_z_out(i_out,j_out) = mat_z_out(i_out,j_out) /    &
                                                   denom
         else
! No interpolation (dx=dy)
            mat_z_out(i_out,j_out) = mat_z_in(i_out,j_out)
      endif
      write(2,'(4(F12.4))') x_out,y_out,mat_z_out(i_out,j_out),                &
         mat_z_out(i_out,j_out)
   enddo
enddo
close(2)
write(*,*) "End Eventual grid interpolation and writing xyz file "
!------------------------
! Deallocations
!------------------------
deallocate(mat_z_in)
deallocate(mat_z_out)
write(*,*) "***  DEM2xyz has terminated "
end program DEM2xyz
