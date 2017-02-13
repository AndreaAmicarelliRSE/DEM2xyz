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
!              The height of the DEM points which belong to the digging/filling 
!              regions (provided in input) is modified. After this treatment, 
!              each digging/filling region has null slope. Variables:
!              n_col_in: number of columns in the input DEM
!              n_col_out: number of columns in the output DEM
!              n_raw: number of raws in the input/output DEM
!              res_fact: resolution factor (ratio between the output and the 
!                 input spatial resolution)
!              n_points_in: number of input vertices
!              n_points_out: number of output vertices
!              n_digging_regions: number of digging/filling regions
!              dx,dy: spatial resolution -final values in (m)-
!              abs_mean_latitude: absolute value of the latitude of the DEM 
!                 barycentre
!              x_in,y_in: horizontal coordinates in input -(m) or (°)- 
!              x_out,y_out: horizontal coordinates in output (m)
!              n_digging_vertices(n_digging_regions): number of vertices of the 
!                 polygon representing a digging/filling region (3-6)
!              z_digging_regions(n_digging_regions): height of the 
!                 digging/filling region
!              digging_certices(n_digging_regions,6,2): output X/Y-coordinates 
!                 (m) of the vertices of the digging/filling regions
!              mat_z_in(n_raw_n_col_in): input DEM
!              mat_z_out(n_raw_n_col_out): output DEM
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
integer :: n_points_out,i_aux,j_aux,n_digging_regions,i_reg,test
double precision :: dx,dy,abs_mean_latitude,denom,distance,x_in,x_out,y_in,y_out
double precision :: point(2)
character(len=100) :: char_aux
integer,dimension(:),allocatable :: n_digging_vertices
double precision,dimension(:),allocatable :: z_digging_regions
double precision,dimension(:,:),allocatable :: mat_z_in,mat_z_out
double precision,dimension(:,:,:),allocatable :: digging_vertices
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine point_inout_convex_non_degenerate_polygon(point,n_sides,         &
                                                        point_pol_1,           &
                                                        point_pol_2,           &
                                                        point_pol_3,           &
                                                        point_pol_4,           &
                                                        point_pol_5,           &
                                                        point_pol_6,test)
      implicit none
      integer(4),intent(in) :: n_sides
      double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
      double precision,intent(in) :: point_pol_3(2),point_pol_4(2)
      double precision,intent(in) :: point_pol_5(2),point_pol_6(2)
      integer(4),intent(inout) :: test
      double precision :: dis1,dis2
      double precision :: normal(2)
   end subroutine point_inout_convex_non_degenerate_polygon
end interface
interface
   subroutine point_inout_hexagon(point,point_pol_1,point_pol_2,point_pol_3,   &
                                  point_pol_4,point_pol_5,point_pol_6,test)
      implicit none
      double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
      double precision,intent(in) :: point_pol_3(2),point_pol_4(2)
      double precision,intent(in) :: point_pol_5(2),point_pol_6(2)
      integer(4),intent(inout) :: test
      integer(4) :: test1,test2,test3,test4
   end subroutine point_inout_hexagon
end interface
interface
   subroutine point_inout_pentagon(point,point_pol_1,point_pol_2,point_pol_3,  &
                                   point_pol_4,point_pol_5,test)
      implicit none
      double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
      double precision,intent(in) :: point_pol_3(2),point_pol_4(2)
      double precision,intent(in) :: point_pol_5(2)
      integer(4),intent(inout) :: test
      integer(4) :: test1,test2,test3
   end subroutine point_inout_pentagon
end interface
interface
   subroutine point_inout_quadrilateral(point,point_pol_1,point_pol_2,         &
                                        point_pol_3,point_pol_4,test)
      implicit none
      double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
      double precision,intent(in) :: point_pol_3(2),point_pol_4(2)
      integer(4),intent(inout) :: test
      integer(4) :: test1,test2
   end subroutine point_inout_quadrilateral
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
n_digging_regions = 0
!------------------------
! Statements
!------------------------
! Subroutine info
write(*,*) "***DEM2xyz is running:"
write(*,*) "This subroutine reads a DEM file"
write(*,*) "    and writes a corresponding xyz file." 
write(*,*) "Algorithm:"
write(*,*) "   1) Reading DEM file, DEM2xyz main input file and pre-processing."
write(*,*) "   2) Eventual grid interpolation, eventual digging/filling DEM ", &
   "regions and writing xyz file. "
write(*,*) "1)  Reading DEM file, DEM2xyz main input file and pre-processing. "
open(1,file='DEM.dem')
read(1,'(a14,i15)') char_aux,n_col_in
read(1,'(a14,i15)') char_aux,n_raw
read(1,'(a)')
read(1,'(a)')
read(1,'(a14,f15.7)') char_aux,dy
allocate(mat_z_in(n_raw,n_col_in))
mat_z_in = 0.d0
read(1,'(a)')
do i_in=1,n_raw
   read(1,*) mat_z_in(i_in,:)
enddo
close(1)
open(1,file='DEM2xyz.inp')
read(1,*) res_fact,abs_mean_latitude,n_digging_regions
if (n_digging_regions>0) then
   allocate(n_digging_vertices(n_digging_regions))
   allocate(z_digging_regions(n_digging_regions))
   allocate(digging_vertices(n_digging_regions,6,2))
   n_digging_vertices(:) = 0
   z_digging_regions(:) = 0.d0
   digging_vertices(:,:,:) = 0.d0
   do i_reg=1,n_digging_regions
      read(1,*) z_digging_regions(i_reg),n_digging_vertices(i_reg)
      do j_aux=1,n_digging_vertices(i_reg)
         read (1,*) digging_vertices(i_reg,j_aux,1:2)
      enddo
   enddo
endif
close(1)
if (abs_mean_latitude>=0.d0) then
! Conversion degrees to radians 
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
allocate(mat_z_out(n_raw,n_col_out))
mat_z_out = 0.d0
write(*,*) "End Reading DEM file, DEM2xyz main input file  and pre-processing. "
write(*,*) "2) Eventual grid interpolation, eventual digging/filling DEM ",    &
   "regions and writing xyz file. "
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
         j_aux = nint((j_out - 0.5) * real(dy / dx) + 0.5)
         i_aux = i_out
         do j_in=(j_aux-1),(j_aux+1)
            do i_in=(i_aux-1),(i_aux+1)
               if ((i_in<1).or.(i_in>n_raw).or.(j_in<1).or.(j_in>n_col_in))    &
                  cycle
               x_in = (j_in - 1) * dx + dx / 2.d0
               y_in = (n_raw + 1 - i_in) * dy - dy / 2.d0
               distance = dsqrt((x_in - x_out) ** 2 + (y_in - y_out) ** 2)
               if (distance<=dy) then
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
      do i_reg=1,n_digging_regions
         test = 0
         point(1) = x_out
         point(2) = y_out
         select case(n_digging_vertices(i_reg))
            case(3)
               call point_inout_convex_non_degenerate_polygon(point,3,         &
                  digging_vertices(i_reg,1,1:2),digging_vertices(i_reg,2,1:2), &
                  digging_vertices(i_reg,3,1:2),point,point,point,test)
            case(4)
               call point_inout_quadrilateral(point,                           &
                  digging_vertices(i_reg,1,1:2),digging_vertices(i_reg,2,1:2), &
                  digging_vertices(i_reg,3,1:2),digging_vertices(i_reg,4,1:2), &
                  test)
            case(5)
               call point_inout_pentagon(point,                                &
                  digging_vertices(i_reg,1,1:2),digging_vertices(i_reg,2,1:2), &
                  digging_vertices(i_reg,3,1:2),digging_vertices(i_reg,4,1:2), &
                  digging_vertices(i_reg,5,1:2),test)
            case(6)
               call point_inout_hexagon(point,digging_vertices(i_reg,1,1:2),   &
                  digging_vertices(i_reg,2,1:2),digging_vertices(i_reg,3,1:2), &
                  digging_vertices(i_reg,4,1:2),digging_vertices(i_reg,5,1:2), &
                  digging_vertices(i_reg,6,1:2),test)
         endselect
         if (test>0) then
            mat_z_out(i_out,j_out) = z_digging_regions(i_reg)
         endif
      enddo
      write(2,'(4(F12.4))') x_out,y_out,mat_z_out(i_out,j_out),                &
         mat_z_out(i_out,j_out)
   enddo
enddo
close(2)
write(*,*) "End Eventual grid interpolation, eventual digging/filling DEM ",   &
   "regions and writing xyz file. "
!------------------------
! Deallocations
!------------------------
deallocate(mat_z_in)
deallocate(mat_z_out)
deallocate(n_digging_vertices)
deallocate(z_digging_regions)
deallocate(digging_vertices)
write(*,*) "***  DEM2xyz has terminated. "
end program DEM2xyz
