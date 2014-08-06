!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Generates a few simple meshes and the associated conntectivity relations
! This would be expected to be replaced with a preprocessor step/read from file
!-------------------------------------------------------------------------------
module mesh_generator_mod

use reference_element_mod, only : nfaces,   nedges,   nverts, &
                                 nfaces_h, nedges_h, nverts_h
use constants_mod,         only : r_def                                 

implicit none

! global numbers of entities in a single 2D layer
integer :: nedge_h_g, nvert_h_g
! global numbers of entities in for full 3D domains
integer :: nface_g, nedge_g, nvert_g

! In terminology of Logg 08 these are:
! cell_next    => MeshConnectivity(3,3) ( cells incident to cells )
! vert_on_cell => MeshConnectivity(3,0) ( vertices incident to cells )
! mesh_vertex  => MeshGeometry

! This is the minimal set of information, from which all other connectivity can be computed

integer,          allocatable :: cell_next(:,:)
integer,          allocatable :: vert_on_cell(:,:)
real(kind=r_def), allocatable :: mesh_vertex(:,:)

! Extra connectivity for easy dofmap computation
! In terminology of Logg 08 these are:
! face_on_cell => MeshConnectivity(3,2) ( faces incident to cells )
! edge_on_cell => MeshConnectivity(3,1) ( edges incident to cells )
integer, allocatable ::  face_on_cell(:,:), edge_on_cell(:,:)

! together this gives all the cell -> d connectivity (3,d), d=0,1,2

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> Subroutine to allocate mesh arrays for connectivity
!> @param[in] ncells the number of cells in a horizontal layer
!> @param[in] nlayers the number of vertical layers
subroutine mesh_generator_init(ncells,nlayers)
!-----------------------------------------------------------------------------
! Subroutine to allocate connectivity
!-----------------------------------------------------------------------------

! number of cells in a horizontal layer
  integer, intent(in) :: ncells 
! number of vertical layers
  integer, intent(in) :: nlayers   

  allocate ( cell_next(ncells*nlayers,nfaces) )
  allocate ( vert_on_cell(ncells*nlayers,nverts) )

end subroutine mesh_generator_init

!> Generate a biperiodic domain of size 'nx * ny' where 'nx * ny = ncells'
!>
!> @param ncells Number of cells in a layer.
!> @param nx, ny Number of cells in X and Y.
!> @param nlayers Number of vertical layers
!> @param dx, dy, dz Cell width in X, Y and Z direction.
!>
subroutine mesh_generator_biperiodic( ncells, nx, ny, nlayers, dx, dy, dz )

  use log_mod, only : log_event, log_scratch_space, &
                      LOG_LEVEL_DEBUG, LOG_LEVEL_ERROR

  integer,              intent( in ) :: ncells
  integer,              intent( in ) :: nx, ny
  integer,              intent( in ) :: nlayers
  real( kind = r_def ), intent( in ) :: dx, dy, dz

  ! Loop indices
  integer         :: i, j, k, id, jd

  if ( nx * ny /= ncells ) then
    write( log_scratch_space, '( A, I0, A, I0)' ) &
         'Incorrect number of elements in mesh_generator_biperiodic', &
         nx * ny, ' cf. ', ncells
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    stop
  end if

! topologically a torus
  nedge_h_g = 2*nx*ny
  nvert_h_g = nx*ny
  
  nface_g = nedge_h_g*nlayers + ncells*(nlayers + 1)
  nedge_g = nedge_h_g*(nlayers + 1) + nlayers*nvert_h_g
  nvert_g = nvert_h_g*(nlayers + 1)
  
! allocate coordinate array
  allocate ( mesh_vertex(nvert_g,3) )

  id = 1
  do j=1,ny
    do i=1,nx
! j-1 cell (South face)
      cell_next(id,1) = id - nx
! i+1 cell (East face)
      cell_next(id,2) = id + 1
! j+1 cell (North face)
      cell_next(id,3) = id + nx      
! i-1 cell (West face)
      cell_next(id,4) = id - 1
! k-1 cell (bottom face)
      cell_next(id,5) = id - nx*ny
! k+1 cell (top face)
      cell_next(id,6) = id + nx*ny
      
      id = id + 1      
    end do
  end do

! Now do periodicity/connectivity along edges 
! South
  id = 1
  do i=1,nx
    cell_next(id,1) = nx*ny-nx+i
    id = id + 1
  end do
  
! North  
  id = nx*ny
  do i=nx,1,-1
    cell_next(id,3) = i
    id = id - 1
  end do
  
  id = 1
  do j=1,ny
! West     
    cell_next(id,4) = id + nx -1
    id = id + nx -1
! East
    cell_next(id,2) = id - nx + 1
    id = id + 1
  end do
  
  ! perform vertical extrusion for connectivity 
  do k=1,nlayers-1
    do i=1,nx*ny
      id = i + k*nx*ny      
      jd = id - nx*ny        
      do j=1,nfaces
        cell_next(id,j) = cell_next(jd,j) + nx*ny
      end do
    end do
  end do
  
  ! Set connectivity at lower/upper boundary to some dummy cell
  do i=1,nx*ny
    cell_next(i,5) = 0
    j =  i+(nlayers-1)*nx*ny
    cell_next(j,6) = 0
  end do
      
! compute vertices on cell
  id = 1
  do j=1,ny
    do i=1,nx
      vert_on_cell(id,1) = id
      vert_on_cell(id,2) = cell_next(id,2)
      vert_on_cell(id,3) = cell_next(cell_next(id,3),2)
      vert_on_cell(id,4) = cell_next(id,3)
      
      jd = id + nx*ny
      vert_on_cell(id,5) = jd
      vert_on_cell(id,6) = cell_next(jd,2)
      vert_on_cell(id,7) = cell_next(cell_next(jd,3),2)
      vert_on_cell(id,8) = cell_next(jd,3)
      
      id = id + 1
    end do
  end do
  
! perform vertical extrusion for connectivity 
  do k=1,nlayers-1
    do i=1,nx*ny
      id = i + k*nx*ny      
      jd = id - nx*ny        
      do j=1,nverts
        vert_on_cell(id,j) =  vert_on_cell(jd,j) + nx*ny
      end do
    end do
  end do
  
! Compute vertices
  id = 1
  do k=1,nlayers+1
    do j=1,ny
      do i=1,nx
        mesh_vertex(id,1) = real(i-1 - nx/2)*dx 
        mesh_vertex(id,2) = real(j-1 - ny/2)*dy 
        mesh_vertex(id,3) = real(k-1)*dz
        id = id + 1
      end do
    end do
  end do
   
! Diagnostic information  
  call log_event( 'grid connectivity', LOG_LEVEL_DEBUG )
  do i = 1, nx * ny * nlayers
    write( log_scratch_space,'(7i6)' ) i, &
                            cell_next(i,1), cell_next(i,2), cell_next(i,3), &
                            cell_next(i,4), cell_next(i,5), cell_next(i,6)
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
  end do
  call log_event( 'verts on cells', LOG_LEVEL_DEBUG )
  do i = 1, nx * ny * nlayers
    write( log_scratch_space, '(9i6)' ) i, &
                             vert_on_cell(i,1), vert_on_cell(i,2), &
                             vert_on_cell(i,3), vert_on_cell(i,4), &
                             vert_on_cell(i,5), vert_on_cell(i,6), &
                             vert_on_cell(i,7), vert_on_cell(i,8)
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
  end do

  call log_event( 'vert coords', LOG_LEVEL_DEBUG )
  do i = 1, nvert_g
    write( log_scratch_space, '(i6,4f8.4)' ) &
         i, mesh_vertex(i,1), mesh_vertex(i,2), mesh_vertex(i,3)
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
  end do

end subroutine mesh_generator_biperiodic

!> Generate a biperiodic domain of size nx * ny where nx * ny = ncells.
!>
!> <pre>
!> Location of panels:
!>             .....
!>            :  3  |
!>            :     |
!>             -----
!>      -----  .....  .....  -----
!>     |  5  :|  1  :|  2  :|  4  :
!>     |     :|     :|     :|     :
!>      .....  -----  -----  .....
!>             .....
!>            |  6  :
!>            |     :
!>             -----
!>
!>     Solid lines: left and bottom edges of panel
!>     Dotted lines: right and top edges of panel
!>     Currently reads in data created from John Thuburns gengrid_cube program
!> </pre>
!>
!> @param filename Data file for cubed sphere data.
!> @param ncells Total number of cells in horizontal layer.
!> @param nlayers Number of vertical layers.
!> @param dz Vertical grid spacing.
!>

subroutine mesh_generator_cubedsphere( filename, ncells, nlayers, dz )

  use log_mod,              only: log_event, log_scratch_space, &
                                  LOG_LEVEL_INFO, LOG_LEVEL_ERROR
  use ugrid_2d_mod,         only: ugrid_2d_type
  use ugrid_file_mod,       only: ugrid_file_type 
  use ncdf_quad_mod,        only: ncdf_quad_type 
  use coord_algorithms_mod, only: llr2xyz
  
  implicit none

  character(*),     intent(in) :: filename
  integer,          intent(in) :: ncells
  integer,          intent(in) :: nlayers
  real(kind=r_def), intent(in) :: dz

! file unit for mesh data file
  integer, parameter :: mesh_data_unit = 44
! Loop indices
  integer :: i, j, k, vert, id, jd
! data from file
  integer :: nvert_in, nface_in, nedge_in
  integer :: num_nodes_per_face, num_nodes_per_edge, num_edges_per_face
! lat/long coordinates
  !2D ugrid and generator strategy
  type(ugrid_2d_type)                 :: ugrid_2d
  class(ugrid_file_type), allocatable :: file_handler
  real(kind=r_def) :: long, lat, r

! topologically a cube
  nedge_h_g = 2*ncells
  nvert_h_g = (ncells + 2)
  
  nface_g = nedge_h_g*nlayers + ncells*(nlayers + 1)
  nedge_g = nedge_h_g*(nlayers + 1) + nvert_h_g*nlayers
  nvert_g = nvert_h_g*(nlayers + 1)

! allocate coordinate array
  allocate ( mesh_vertex(nvert_g,3) )

  allocate(ncdf_quad_type :: file_handler)
  call ugrid_2d%set_file_handler(file_handler)
  call ugrid_2d%read_from_file(trim(filename))

  call ugrid_2d%get_dimensions(                &
     num_nodes          = nvert_in,            &
     num_edges          = nedge_in,            &
     num_faces          = nface_in,            &
     num_nodes_per_face = num_nodes_per_face,  &
     num_edges_per_face = num_edges_per_face,  &
     num_nodes_per_edge = num_nodes_per_edge)

  if ( nface_in /= ncells ) then
    write( log_scratch_space, '(A, I0, A, I0)' ) &
         'Number of cells in ugrid file does not match: ', &
         nface_in, ' vs. ', ncells
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    stop
  end if
  if ( nvert_in /= nvert_h_g ) then
    write( log_scratch_space, '(A, I0, A, I0)' ) &
         'Number of vertices in ugrid file does not match: ', &
         nvert_in, ' vs. ', nvert_h_g
    stop
  end if  

  !Get coordinates of vertices
  call ugrid_2d%get_node_coords_transpose( mesh_vertex(1:nvert_h_g,1:3) )
  call ugrid_2d%get_face_node_connectivity_transpose(vert_on_cell(1:nface_in,1:num_nodes_per_face))
  call ugrid_2d%get_face_face_connectivity_transpose(cell_next   (1:nface_in,1:num_edges_per_face))
  
   !mesh_vertex(1:nvert_h_g,3) = earth_radius

! add connectivity for up/down
  do j=1,ncells
    cell_next(j,5) = j - ncells
    cell_next(j,6) = j + ncells
  end do
       
! perform vertical extrusion for connectivity 
  do k=1,nlayers-1
    do i=1,ncells
      id = i + k*ncells      
      jd = id - ncells        
      do j=1,nfaces
        cell_next(id,j) = cell_next(jd,j) + ncells
      end do
    end do
  end do
    
! Set connectivity at lower/upper boundary to some dummy cell
  do i=1,ncells
    cell_next(i,5) = 0
    j =  i+(nlayers-1)*ncells
    cell_next(j,6) = 0
  end do
    
! perform vertical extrusion for vertices
  do j=1,nvert_in
    do k=1,nlayers                       
      mesh_vertex(j+k*nvert_in,1) = mesh_vertex(j,1)
      mesh_vertex(j+k*nvert_in,2) = mesh_vertex(j,2) 
      mesh_vertex(j+k*nvert_in,3) = mesh_vertex(j,3) + dz*real(k)
    end do
  end do

! Convert (long,lat,r) -> (x,y,z)
  do j=1,nvert_in
    do k=0,nlayers
      long = mesh_vertex(j+k*nvert_in,1)
      lat  = mesh_vertex(j+k*nvert_in,2)
      r    = mesh_vertex(j+k*nvert_in,3)
      call llr2xyz(long,lat,r,mesh_vertex(j+k*nvert_in,1),                  &
                              mesh_vertex(j+k*nvert_in,2),                  &
                              mesh_vertex(j+k*nvert_in,3))
    end do
  end do
  
! assign vertices to cells
  do j=1,ncells
    vert_on_cell(j,5) = vert_on_cell(j,1) + nvert_in
    vert_on_cell(j,6) = vert_on_cell(j,2) + nvert_in
    vert_on_cell(j,7) = vert_on_cell(j,3) + nvert_in
    vert_on_cell(j,8) = vert_on_cell(j,4) + nvert_in
! do vertical extrusion
    do k=1,nlayers-1
      do vert=1,8
        vert_on_cell(j+k*ncells,vert) = vert_on_cell(j,vert) + k*nvert_in
      end do
    end do
  end do  
 
! Diagnostic information  
  call log_event( 'grid connectivity', LOG_LEVEL_INFO )
  do i = 1, ncells * nlayers
    write( log_scratch_space, '(7i6)' ) i, &
                              cell_next(i,1), cell_next(i,2), cell_next(i,3), &
                              cell_next(i,4), cell_next(i,5), cell_next(i,6)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  end do
  call log_event( 'verts on cells', LOG_LEVEL_INFO )
  do i = 1, ncells * nlayers
    write( log_scratch_space, '(9i6)' ) i, &
                              vert_on_cell(i,1), vert_on_cell(i,2), &
                              vert_on_cell(i,3), vert_on_cell(i,4), &
                              vert_on_cell(i,5), vert_on_cell(i,6), &
                              vert_on_cell(i,7), vert_on_cell(i,8)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  end do

  call log_event( 'vert coords', LOG_LEVEL_INFO )
  do i = 1, nvert_g
    write( log_scratch_space, '(i6,4e16.8)' ) &
         i, mesh_vertex(i,1), mesh_vertex(i,2), mesh_vertex(i,3)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  end do 

end subroutine mesh_generator_cubedsphere

!> Compute mesh connectivity.
!>
!> <pre>
!> cells->faces (3,2)
!> cells->edges (3,1)
!> </pre>
!>
!> @param ncells Total number of cells in horizontal layer.
!>
subroutine mesh_connectivity( ncells )

  use log_mod, only : log_event, log_scratch_space, &
                      LOG_LEVEL_INFO, LOG_LEVEL_ERROR

  integer, intent( in ) :: ncells
 
  integer :: i, j, k, l, m, n, ij, inxt, inxtnxt, jnxt, vert
! Number of entities for a single layer  
  integer :: nedge_layer, nface_layer

  nedge_layer = 2*nedge_h_g + nvert_h_g
  nface_layer = nedge_h_g + 2*ncells

  allocate( face_on_cell(ncells,nfaces) )
  allocate( edge_on_cell(ncells,nedges) )
  face_on_cell(:,:) = 0
  edge_on_cell(:,:) = 0
  
  ij = 1
  do i=1,ncells
    do j=1,nfaces_h
      if ( face_on_cell(i,j) == 0 ) then
        face_on_cell(i,j) = ij
! find matching face
        inxt = cell_next(i,j)
        do k=1,nfaces_h
          if ( cell_next(inxt,k) == i ) then
            jnxt = k
          end if
        end do
        face_on_cell(inxt,jnxt) = ij
        ij = ij + 1
      end if                 
    end do
    do j=nfaces_h+1,nfaces
      face_on_cell(i,j) = ij
      ij = ij + 1
    end do
  end do
  if ( maxval(face_on_cell) /= nface_layer ) then
    call log_event( 'Error computing face on cell connectivity', &
                    LOG_LEVEL_ERROR )
    stop
  end if
  
  ij = 1
  do i=1,ncells
! horizontal edges  
    do j=1,nedges_h
      if ( edge_on_cell(i,j) == 0 ) then  
        edge_on_cell(i,j) = ij
        edge_on_cell(i,j+nedges_h+nverts_h) = ij+1
! find matching edge in neighbouring cell
        inxt = cell_next(i,j)
        do k=1,nedges_h
          if ( cell_next(inxt,k) == i ) then
            jnxt = k
          end if
        end do
        edge_on_cell(inxt,jnxt) = ij
        edge_on_cell(inxt,jnxt+nedges_h+nverts_h) = ij+1
        ij = ij + 2
      end if
    end do
! vertical edges 
    do j=nedges_h+1,nedges_h+nverts_h
      if ( edge_on_cell(i,j) == 0 ) then
        edge_on_cell(i,j) = ij
! find matching edge on two neighbouring cells 
! this edge is an extrusion of the corresponding vertex
        vert = vert_on_cell(i,j-nedges_h)
        do k=1,nedges_h
          inxt = cell_next(i,k)
          do l=1,nverts_h
            if ( vert_on_cell(inxt,l) == vert ) then
              edge_on_cell(inxt,l+nedges_h) = ij
              do m=1,nedges_h
                inxtnxt = cell_next(inxt,m)
                do n=1,nverts_h
                  if ( vert_on_cell(inxtnxt,n) == vert ) then
                    edge_on_cell(inxtnxt,n+nedges_h) = ij 
                  end if
                end do
              end do
            end if
          end do
        end do
        ij = ij + 1
      end if
    end do
  end do

  if ( maxval(edge_on_cell) /= nedge_layer ) then
    write( log_scratch_space, '(A, I0, A, I0)' ) &
         'Error computing edge on cell connectivity: ', &
         maxval( edge_on_cell ), ' vs. ', nedge_layer
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    stop
  end if

end subroutine mesh_connectivity

subroutine llr2xyz(long,lat,r,x,y,z)
!-------------------------------------------------------------------------------
!  Subroutine to convert longitude and latitude to cartesian coordinates
!-------------------------------------------------------------------------------
      
  real(kind=r_def), intent(in)  :: long,lat,r
  real(kind=r_def), intent(out) :: x,y,z
  real(kind=r_def)              :: cln,sln,clt,slt

  sln=sin(long)
  cln=cos(long)
  slt=sin(lat)
  clt=cos(lat)

  x=r*cln*clt
  y=r*sln*clt
  z=r*slt

  return
      
end  subroutine llr2xyz
!-------------------------------------------------------------------------------
      
subroutine xyz2llr(x,y,z,long,lat,r)
  
!-------------------------------------------------------------------------------
!  Subroutine to convert cartesian coordinates to longitude and latitude
!------------------------------------------------------------------------------- 

  use constants_mod, only: pi


  real(kind=r_def), intent(in)  :: x, y, z
  real(kind=r_def), intent(out) :: long, lat, r
  real(kind=r_def)              :: tln, tlt
  real(kind=r_def)              :: tol = 10e-8_r_def

  if ( abs(x) < tol ) then
    if ( y >= 0.0_r_def ) then
      long = 0.5_r_def*pi
    else
      long = 1.5_r_def*pi
    end if
  else
    tln = y/x
    long = atan(tln)
    if ( x < 0.0_r_def ) then
      long = long + pi
    end if
    if ( long < 0.0_r_def ) then
      long = long + 2.0_r_def*pi
    end if
  end if

  r = sqrt(x*x + y*y)
  if ( abs(r) < tol ) then
    if (z > 0.0_r_def ) then
      lat =  0.5_r_def*pi
    else
      lat = -0.5_r_def*pi
    end if
  else
    tlt = z/r
    lat = atan(tlt)
  end if
  
  r = sqrt(x*x + y*y + z*z)  

end subroutine xyz2llr

!-------------------------------------------------------------------------------
!> @brief Subroutine to return the vertex coordinates for a column of cells
!> @param[in] cell the horizontal cell index
!> @param[in] ncells the number of horizontal cells
!> @param[in] nlayers the number of vertical layers
!> @param[out] vert_coords the coordinates of the vertices
subroutine get_cell_coords(cell, ncells, nlayers, vert_coords)

!-------------------------------------------------------------------------------
!  Subroutine to return coordinates of vertices on cell 
!------------------------------------------------------------------------------- 

  integer, intent(in) :: cell, ncells, nlayers
  real(kind=r_def), intent(out) :: vert_coords(3,nverts,nlayers)
  
  integer :: i, j, k, cell_idx, v
  
  do k = 1,nlayers
    cell_idx = cell+(k-1)*ncells
    do j = 1,nverts
      v = vert_on_cell(cell_idx,j)
      do i = 1,3      
        vert_coords(i,j,k) = mesh_vertex(v,i)
      end do
    end do
  end do

end subroutine get_cell_coords

end module mesh_generator_mod
