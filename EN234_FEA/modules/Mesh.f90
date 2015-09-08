module Mesh
   use Types

!  Data type for nodes   
   type node
      sequence
      integer :: flag                          ! Integer identifier
      integer :: coord_index                   ! Index of first coordinate in coordinate array
      integer :: n_coords                      ! Total no. coordinates for the node
      integer :: dof_index                     ! Index of first DOF in dof array
      integer :: n_dof                         ! Total no. of DOF for node
      integer :: displacement_map_index        ! Index of displacement node map
      integer :: n_displacements               ! No. displacement DOFs
   end type node
!  Data type for elements   
   type element
      sequence
      integer :: flag                          ! Integer identifier for element
      integer :: connect_index                 ! Index of first node on element in connectivity(:)
      integer :: n_nodes                       ! No. nodes on the element
      integer :: state_index                   ! Index of first state variable in element_state_variables(:)
      integer :: n_states                      ! No. state variables
      integer :: element_property_index        ! Index of first element property in element_properties(:)
      integer :: n_element_properties          ! No. element properties
      integer :: density_index                 ! Index of density value
   end type element
   
   type zone
      sequence
      integer :: start_element                 ! First element in a zone
      integer :: end_element                   ! Last element in a zone
   end type zone
   
   integer, save :: n_zones                          ! Number of zones
   integer, save :: n_nodes                          ! Total number of nodes
   integer, save :: n_elements                       ! Total number of elements
   integer, save :: length_coords                    ! Length of coordinate array
   integer, save :: length_dofs                      ! Length of nodal DOF array
   integer, save :: length_connectivity              ! Length of connectivity array
   integer, save :: length_element_properties        ! Length of element property array
   integer, save :: length_densities                 ! Length of density array
   integer, save :: length_state_variables           ! Length of state variable array
   integer, save :: length_displacement_map          ! Length of array mapping nodal DOF to displacements
   integer, save :: n_mesh_parameters                ! Parameeters controlling a user-subrouine generated mesh
   
   integer, save, allocatable :: displacement_map(:) ! Array storing mapping of nodal DOF to displacements
   integer, save, allocatable :: connectivity(:)     ! Array storing element connectivity

   real (prec), save :: nodal_force_norm             ! Norm of nodal generalized force
   real (prec), save :: unbalanced_force_norm        ! Norm of out-of-balance forces
   real (prec), save :: correction_norm              ! Norm of solution correction

   real (prec), save, allocatable :: element_properties(:)            ! List of element properties
   real (prec), save, allocatable :: densities(:)                     ! List of density values for zones
   real (prec), save, allocatable :: initial_state_variables(:)       ! Element state variables at the start of a time increment
   real (prec), save, allocatable :: updated_state_variables(:)       ! Element state variables at the end of a time increment
   real (prec), save, allocatable :: coords(:)                        ! List of nodal coordinates
   real (prec), save, allocatable :: dof_total(:)                     ! List of accumulated DOF
   real (prec), save, allocatable :: dof_increment(:)                 ! List of increment in DOF
   real (prec), save, allocatable :: velocity(:)                      ! Velocity (for explicit dynamics)
   real (prec), save, allocatable :: acceleration(:)                  ! Acceleration (for explicit dynamics)
   real (prec), save, allocatable :: lumped_mass(:)                   ! Lumped mass matrix (for explicit dynamics)
   real (prec), save, allocatable :: rforce(:)
   real (prec), save, allocatable :: mesh_subroutine_parameters(:)    ! List of parameters controlling a user-subroutine generated mesh
   
   type (node), save, allocatable :: node_list(:)
   type (element), save, allocatable :: element_list(:)
   type (zone), save, allocatable :: zone_list(:)
   
   character (len=100), allocatable :: zone_namelist(:)         ! Names of zones in mesh

   logical, save, allocatable :: element_deleted(:)             ! Flag listing deleted elements in an explicit dynamic simulation



contains

   subroutine extract_element_data(lmn,flag,n_nodes,node_list,n_properties,properties,n_state_variables,initial_svars,updated_svars)

      use Types
      use ParamIO

      implicit none

      integer, intent( in )  ::   lmn

      integer, intent( out ) :: node_list(*)
      integer, intent( out ) :: flag
      integer, intent( out ) :: n_nodes
      integer, intent( out ) :: n_properties
      integer, intent( out ) :: n_state_variables

      real (prec), intent( out ) ::  properties(*)
      real (prec), intent( out ) ::  initial_svars(*)
      real (prec), intent( out ) ::  updated_svars(*)

      !  Function to extract data for element number lmn
      !
      if (lmn>n_elements) then
          write(IOW,*) ' Error in subroutine extract_element_data '
          write(IOW,*) ' Element number ',lmn,' exceeds the number of elements in the mesh '
          stop
      endif

      flag = element_list(lmn)%flag
      n_nodes = element_list(lmn)%n_nodes
      n_state_variables = element_list(lmn)%n_states
      n_properties = element_list(lmn)%n_element_properties

!      if (size(properties)< n_properties ) then
!         write(IOW,*) ' Error in subroutine extract_element_data '
!         write(IOW,*) ' The element has ',n_properties,' properties '
!         write(IOW,*) ' but the length of the property array is ',size(properties)
!         stop
!      endif
!
!      if (size(node_list)< n_nodes ) then
!         write(IOW,*) ' Error in subroutine extract_element_data '
!         write(IOW,*) ' The element has ',n_nodes,' nodes '
!         write(IOW,*) ' but the length of the node list array is ',size(properties)
!         stop
!      endif
!
!      if (size(initial_svars)< n_state_variables ) then
!         write(IOW,*) ' Error in subroutine extract_element_data '
!         write(IOW,*) ' The element has ',n_state_variables,' state variables '
!         write(IOW,*) ' but the length of the initial state variable array is ',size(initial_svars)
!         stop
!      endif
!
!
!      if (size(updated_svars)< n_state_variables ) then
!         write(IOW,*) ' Error in subroutine extract_element_data '
!         write(IOW,*) ' The element has ',n_state_variables,' state variables '
!         write(IOW,*) ' but the length of the initial state variable array is ',size(updated_svars)
!         stop
!      endif

      node_list(1:n_nodes) = connectivity(element_list(lmn)%connect_index:element_list(lmn)%connect_index+n_nodes-1)
      if (n_properties>0) then
          properties(1:n_properties) = element_properties(element_list(lmn)%element_property_index: &
              element_list(lmn)%element_property_index+n_properties-1)
      endif
      if (n_state_variables>0) then
          initial_svars(1:n_state_variables) = initial_state_variables(element_list(lmn)%state_index: &
              element_list(lmn)%state_index+n_state_variables-1)
          updated_svars(1:n_state_variables) = updated_state_variables(element_list(lmn)%state_index: &
              element_list(lmn)%state_index+n_state_variables-1)
      endif


   end subroutine extract_element_data

   subroutine extract_node_data(nn,flag,n_coords,nodal_coords,n_dof,nodal_dof_increment,nodal_dof_total)

      use Types
      use ParamIO

      implicit none

      integer, intent( in )  ::  nn

      integer, intent( out ) :: n_coords
      integer, intent( out ) :: flag
      integer, intent( out ) :: n_dof

      real (prec), intent( out ) ::  nodal_coords(*)
      real (prec), intent( out ) ::  nodal_dof_increment(*)
      real (prec), intent( out ) ::  nodal_dof_total(*)

      !  Function to extract data for node number nn
      !
      if (nn>n_nodes) then
          write(IOW,*) ' Error in subroutine exatract_node_data '
          write(IOW,*) ' Node number ',nn,' exceeds the number of nodes in the mesh '
          stop
      endif

      flag = node_list(nn)%flag
      n_coords = node_list(nn)%n_coords
      n_dof = node_list(nn)%n_dof

!      if (size(nodal_coords)< n_coords ) then
!         write(IOW,*) ' Error in subroutine extract_node_data '
!         write(IOW,*) ' The node has ',n_coords,' coordinates '
!         write(IOW,*) ' but the length of the coordinate array is ',size(nodal_coords)
!         stop
!      endif
!
!      if (size(nodal_dof_increment)< n_dof ) then
!         write(IOW,*) ' Error in subroutine extract_node_data '
!         write(IOW,*) ' The node has ',n_dof,' DOF '
!         write(IOW,*) ' but the length of the total DOF array is ',size(nodal_dof_total)
!         stop
!      endif
!
!      if (size(nodal_dof_total)< n_dof ) then
!         write(IOW,*) ' Error in subroutine extract_node_data '
!         write(IOW,*) ' The node has ',n_dof,' DOF '
!         write(IOW,*) ' but the length of the DOF increment array is ',size(nodal_dof_increment)
!         stop
!      endif

      nodal_coords(1:n_coords) = coords(node_list(nn)%coord_index:node_list(nn)%coord_index+n_coords-1)
      nodal_dof_increment(1:n_dof) = dof_increment(node_list(nn)%dof_index: &
                                                  node_list(nn)%dof_index + n_dof-1)
      nodal_dof_total(1:n_dof) = dof_total(node_list(nn)%dof_index: &
                                                  node_list(nn)%dof_index + n_dof-1)


    end subroutine extract_node_data


end module
