!     Subroutines for basic 3D linear elastic elements 



!==========================SUBROUTINE el_timoshenko_3d ==============================
subroutine el_timoshenko_3d(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn ([u1,v1,w1,phi1,psi1,theta1,...])
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
  
    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          

    ! Local Variables
     integer      :: kint

    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  intForce(6)                       ! [N, Mz, My, Qy, Qz, Mx]
    real (prec)  ::  strain(6), dstrain(6)             ! generalized strains [e, ky, kz, gy, gz, kx]
    real (prec)  ::  D(6,6)                            ! Material matrix
    real (prec)  ::  B(6,12)                           ! matrix to convert dof to generalized strains
    real (prec)  ::  N(18), dNdx(18)                   ! array of shape functions and derivatives [N1, N2, Hv1, Hv2, Hw1, Hw2, Ht1, Ht2, Hp1, Hp2, Gv1, Gv2, Gw1, Gw2, Gt1, Gt2, Gp1, Gp2]
    real (prec)  ::  dxidx, det, w, xis(2), xi         ! values for gaussian quadrature
    real (prec)  ::  E, G, A, k, Iz, Iy                ! Material properties
    real (prec)  ::  L                                 ! Calculated Length
    real (prec)  :: ay, by, az, bz                     ! constants in shape functions
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Shear Modulus
    !     element_properties(3)         Area
    !     element_properties(4)         Timoshenko beam coefficient
    !     element_properties(5,6)       Moments of area

    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))


    element_residual = 0.d0
    element_stiffness = 0.d0
	
	! pull out properties and set up D matrix and shape function constatns
    E = element_properties(1)
    G = element_properties(2)
    A = element_properties(3)
    k = element_properties(4)
    Iz = element_properties(5)
    Iy = element_properties(6)

    D = 0.d0
    D(1,1) = E*A
    D(2,2) = E*Iy
    D(3,3) = E*Iz
    D(4,4) = k*G*A
    D(5,5) = k*G*A
    D(6,6) = k*G*(Iy+Iz)

    L = sqrt((x(1,1)-x(1,2))**2+(x(2,1)-x(2,2))**2+(x(3,1)-x(3,2))**2)

    ay = 12.d0*E*Iy/(k*G*A*L**2)
    by = 1.d0/(1.d0-ay)

    az = 12.d0*E*Iz/(k*G*A*L**2)
    bz = 1.d0/(1.d0-az)


    ! set up interation points & etc
    xis(1) = -sqrt(1.d0/3.d0)
    xis(2) = sqrt(1.d0/3.d0)
    w = 1.d0
    dxidx = 1.d0/L !note difference from typical because of luo's convention of xi = 0->1
    det = L/2.d0

    do kint = 1,2
    ! calculate shape functions and derivatives
    N = 0.d0
    dNdx = 0.d0
    xi = (xis(kint)+1.d0)/2.d0 !convert to luo's convention

    N(1) = 1-xi
    N(2) = xi
    N(3) = by*(2.d0*xi**3-3.d0*xi**2+ay*xi+1.d0-ay)
    N(4) = by*(-2.d0*xi**3+3.d0*xi**2-ay*xi)
    N(5) = bz*(2.d0*xi**3-3.d0*xi**2+az*xi+1.d0-az)
    N(6) = bz*(-2.d0*xi**3+3.d0*xi**2-az*xi)
    N(7) = L*by*(xi**3+(ay/2.d0-2.d0)*xi**2+(1.d0-ay/2.d0)*xi)
    N(8) = L*by*(xi**3-(1.d0+ay/2.d0)*xi**2+ay*xi/2.d0)
    N(9) = L*bz*(xi**3+(az/2.d0-2.d0)*xi**2+(1.d0-az/2.d0)*xi)
    N(10) = L*bz*(xi**3-(1.d0+az/2.d0)*xi**2+az*xi/2.d0)
    N(11) = 6.d0*by/L*(xi**2-xi)
    N(12) = 6.d0*by/L*(-(xi**2)+xi)
    N(13) = 6.d0*bz/L*(xi**2-xi)
    N(14) = 6.d0*bz/L*(-(xi**2)+xi)
    N(15) = by*(3.d0*xi**2+(ay-4.d0)*xi+1.d0-ay)
    N(16) = by*(3.d0*xi**2-(ay+2.d0)*xi)
    N(17) = bz*(3.d0*xi**2+(az-4.d0)*xi+1.d0-az)
    N(18) = bz*(3.d0*xi**2-(az+2.d0)*xi)

    dNdx(1) = -1
    dNdx(2) = 1
    dNdx(3) = by*(ay - 6.d0*xi + 6.d0*xi**2)
    dNdx(4) = by*(-ay + 6.d0*xi - 6.d0*xi**2)
    dNdx(5) = bz*(az - 6.d0*xi + 6.d0*xi**2)
    dNdx(6) = bz*(-az + 6.d0*xi - 6.d0*xi**2)
    dNdx(7) = by*L*(1.d0 - ay/2.d0 + 2.d0*(-2.d0 + ay/2.d0)*xi + 3.d0*xi**2)
    dNdx(8) = by*L*(ay/2.d0 - 2.d0*(1.d0 + ay/2.d0)*xi + 3.d0*xi**2)
    dNdx(9) = bz*L*(1.d0 - az/2.d0 + 2.d0*(-2.d0 + az/2.d0)*xi + 3.d0*xi**2)
    dNdx(10) = bz*L*(az/2.d0 - 2.d0*(1.d0 + az/2.d0)*xi + 3.d0*xi**2)
    dNdx(11) = (6.d0*by*(-1.d0 + 2.d0*xi))/L
    dNdx(12) = (6.d0*by*(1.d0 - 2.d0*xi))/L
    dNdx(13) = (6.d0*bz*(-1.d0 + 2.d0*xi))/L
    dNdx(14) = (6.d0*bz*(1.d0 - 2.d0*xi))/L
    dNdx(15) = by*(-4.d0 + ay + 6.d0*xi)
    dNdx(16) = by*(-2.d0 - ay + 6.d0*xi)
    dNdx(17) = bz*(-4.d0 + az + 6.d0*xi)
    dNdx(18) = bz*(-2.d0 - az + 6.d0*xi)

    dNdx = dNdx*dxidx


    ! Populate B matrix
    B = 0.d0
    B(1,1:7:6) = dNdx(1:2)
    B(2,2:8:6) = -dNdx(11:12)
    B(2,6:12:6) = -dNdx(15:16)
    B(3,3:9:6) = dNdx(13:14)
    B(3,5:11:6) = dNdx(17:18)
    B(4,2:8:6) = dNdx(3:4)-N(11:12)
    B(4,6:12:6) = dNdx(7:8)-N(15:16)
    B(5,3:9:6) = dNdx(5:6)+N(13:14)
    B(5,5:11:6) = dNdx(9:10)+N(17:18)
    B(6,4:10:6) = dNdx(1:2)

    element_stiffness = element_stiffness + w*det*matmul(transpose(B),matmul(D,B))

    !guessing here
    strain = matmul(B,dof_total)
    dstrain = matmul(B,dof_increment)
    intForce = matmul(D,strain+dstrain)
    element_residual = element_residual - w*det*matmul(transpose(B),intForce)
    end do
  
    return
end subroutine el_timoshenko_3d


!==========================SUBROUTINE fieldvars_timoshenko_3d ==============================
subroutine fieldvars_timoshenko_3d(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only: dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only: dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! # field variables

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step
             
    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables
  
    ! Local Variables
    logical      :: strcmp
  
    integer      :: n_points,kint,k

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  B_aug(6,length_dof_array)         ! array to augment B in the case of B-bar elements
    real (prec)  ::  vol                               ! volume of element
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    real (prec)  :: p, smises                          ! Pressure and Mises stress
    !
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio
  
    return
end subroutine fieldvars_timoshenko_3d

