!     Subroutines for basic 3D linear elastic elements 



!==========================SUBROUTINE el_dynfrac ==============================
subroutine el_dynfrac(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
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
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
  
    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          

    ! Local Variables
    integer      :: n_points,kint

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  B_aug(6,length_dof_array)         ! array to augment B in the case of B-bar elements
    real (prec)  ::  vol                               ! volume of element
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0
	
    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1+xnu) 
    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    D(1:3,1:3) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44

    ! need to define volume averages in the case of b-bar -- cannot do inline with B calculation :(
    if (element_identifier == 1002) then
        dNbardx = 0.d0
        vol = 0.d0
        do kint = 1, n_points
            call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
            dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
            call invert_small(dxdxi,dxidx,determinant)
            dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

            dNbardx = dNbardx + dNdx*w(kint)*determinant
            vol = vol + w(kint)*determinant
        end do

        dNbardx = (1/vol)*dNbardx
    end if

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        if (element_identifier == 1002) then
            !B-Bar element
            B_aug = 0.d0
            B_aug(1,1:3*n_nodes-2:3) = dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1)
            B_aug(1,2:3*n_nodes-1:3) = dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2)
            B_aug(1,3:3*n_nodes:3) = dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3)
            B_aug(2,1:3*n_nodes-2:3) = dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1)
            B_aug(2,2:3*n_nodes-1:3) = dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2)
            B_aug(2,3:3*n_nodes:3) = dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3)
            B_aug(3,1:3*n_nodes-2:3) = dNbardx(1:n_nodes,1)-dNdx(1:n_nodes,1)
            B_aug(3,2:3*n_nodes-1:3) = dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2)
            B_aug(3,3:3*n_nodes:3) = dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3)
            B = B + (1.d0/3.d0)*B_aug
        end if
        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
        stress = matmul(D,strain+dstrain)
        write(IOW,*) 'stress', stress
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            + matmul(transpose(B(1:6,1:3*n_nodes)),matmul(D,B(1:6,1:3*n_nodes)))*w(kint)*determinant

    end do
  
    return
end subroutine el_dynfrac


!==========================SUBROUTINE el_dynfrac_dynamic ==============================
subroutine el_dynfrac_dynamic(lmn, element_identifier, n_nodes, node_property_list, &                    ! Input variables
    n_properties, element_properties,element_coords, length_coord_array, &                               ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                       ! Input variables
    n_state_variables, initial_state_variables, &                                                        ! Input variables
    updated_state_variables,element_residual,element_deleted)                                            ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : dNbardy => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Globals, only : TIME, DTIME
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
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  stress_0, emat, Vf for each element, as in hw
               
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          
    logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete element

    ! Local Variables
    integer      :: n_points,kint

    real (prec)  ::  B(6,length_dof_array)
    real (prec)  ::  B_aug(6,length_dof_array)         ! augmentation matrix for B-bar
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  E, xnu, Y, e0, m, q1, q2, q3      ! Material properties
    real (prec)  ::  fN, eN, sN, fc, fF                ! Material properties (cont)
    real (prec)  ::  eta, deleta                       ! volume averaged things
    real (prec)  ::  delF(3,3), Fmid(3,3)              ! deformation gradient tensors (and determinant)
    real (prec)  ::  Fmid_inv(3,3), J                  ! deformation gradient tensors (and determinant) cont
    real (prec)  ::  vol                               ! element volume
    real (prec)  ::  delL(3,3), delLkk                 ! velocity gradient and it's trace
    real (prec)  ::  dNdy(n_nodes,3)                   ! Shape function derivatives w.r.t deformed coords
    real (prec)  ::  I(3,3)                            ! Identity
    real (prec)  ::  dispmid(length_dof_array)         ! midpoint displacement
    real (prec)  ::  delLbar(3,3)                      ! corrected velocity gradient
    real (prec)  ::  delFlkFklinv                      ! looks like I had a seizure, but I swear I'm ok.
    real (prec)  ::  dele(3,3), delW(3,3), delR(3,3)   ! strain, spin, and rotation increments
    real (prec)  ::  delR_var(3,3), delR_var_inv(3,3)  ! stuff for computing delR
    real (prec)  ::  delR_var_det                      ! stuff for computing delR (cont)
    real (prec)  ::  stress_0(6), stress_1(6)          ! stress before and after step
    real (prec)  ::  ebar_0, ebar_1, Vf_0, Vf_1
    !
    !     Subroutine to compute element force vector for a linear elastodynamic problem
    !     El props are:

    E = element_properties(1)
    xnu = element_properties(2)
    Y = element_properties(3)
    e0 = element_properties(4)
    m = element_properties(5)
    q1 = element_properties(6)
    q2 = element_properties(7)
    q3 = element_properties(8)
    fN = element_properties(9)
    eN = element_properties(10)
    sN = element_properties(11)
    fc = element_properties(12)
    fF = element_properties(13)
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)
  
    I = 0.d0
    I(1,1) = 1.d0
    I(2,2) = 1.d0
    I(3,3) = 1.d0

    dispmid = dof_total + .5d0*dof_increment
    vol = 0.d0
    eta = 0.d0
    deleta = 0.d0
    dNbardy = 0.d0
    !     --  Loop over integration points to calculate vol averaged quantities
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        ! calculate delF
        delF = 0.d0
        delF(1,:) = matmul(dof_increment(1:3*n_nodes-2:3),dNdx(1:n_nodes,1:3))
        delF(2,:) = matmul(dof_increment(2:3*n_nodes-1:3),dNdx(1:n_nodes,1:3))
        delF(3,:) = matmul(dof_increment(3:3*n_nodes:3),dNdx(1:n_nodes,1:3))

        ! calculate Fmid
        Fmid = 0.d0
        Fmid_inv = 0.d0
        J = 0.d0
        Fmid(1,:) = matmul(dispmid(1:3*n_nodes-2:3),dNdx(1:n_nodes,1:3))
        Fmid(2,:) = matmul(dispmid(2:3*n_nodes-1:3),dNdx(1:n_nodes,1:3))
        Fmid(3,:) = matmul(dispmid(3:3*n_nodes:3),dNdx(1:n_nodes,1:3))
        Fmid = Fmid + I
        call invert_small(Fmid, Fmid_inv,J)

        delL = matmul(delF,Fmid_inv)
        delLkk = delL(1,1)+delL(2,2)+delL(3,3)

        !element volume contribution
        vol = vol + w(kint)*determinant

        !contributions to volume averaged jacobian and average volumetric strain
        eta = eta + J*w(kint)*determinant
        deleta = deleta + J*delLkk*w(kint)*determinant

        ! calculate dNdy
        dNdy = 0.d0
        dNdy(:,1) = matmul(dNdx(1:n_nodes,:),Fmid_inv(:,1))
        dNdy(:,2) = matmul(dNdx(1:n_nodes,:),Fmid_inv(:,2))
        dNdy(:,3) = matmul(dNdx(1:n_nodes,:),Fmid_inv(:,3))

        !volume averaged spatial shape function derivatives contribution
        dNbardy = dNbardy + dNdy*w(kint)*determinant
    end do

    !calculate volume averages
    eta = eta/vol
    deleta = deleta/(eta*vol)
    dNbardy = dNbardy/(eta*vol)

    element_deleted = .true.
    !     --  Loop over integration points to calc velocity and rotation and stress
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        ! calculate delF
        delF = 0.d0
        delF(1,:) = matmul(dof_increment(1:3*n_nodes-2:3),dNdx(1:n_nodes,1:3))
        delF(2,:) = matmul(dof_increment(2:3*n_nodes-1:3),dNdx(1:n_nodes,1:3))
        delF(3,:) = matmul(dof_increment(3:3*n_nodes:3),dNdx(1:n_nodes,1:3))

        ! calculate Fmid
        Fmid = 0.d0
        Fmid_inv = 0.d0
        J = 0.d0
        Fmid(1,:) = matmul(dispmid(1:3*n_nodes-2:3),dNdx(1:n_nodes,1:3))
        Fmid(2,:) = matmul(dispmid(2:3*n_nodes-1:3),dNdx(1:n_nodes,1:3))
        Fmid(3,:) = matmul(dispmid(3:3*n_nodes:3),dNdx(1:n_nodes,1:3))
        Fmid(1,1) = Fmid(1,1)+1.d0
        Fmid(2,2) = Fmid(2,2)+1.d0
        Fmid(3,3) = Fmid(3,3)+1.d0
        call invert_small(Fmid, Fmid_inv,J)

        !corrected velocity gradient increment
        delLbar = matmul(delF,Fmid_inv)
        delFlkFklinv = delLbar(1,1) + delLbar(2,2) + delLbar(3,3)
        delLbar = delLbar + I*(deleta-delFlkFklinv)/3.d0

        !strain and spin increments
        dele = (delLbar + transpose(delLbar))/2.d0
        delW = (delLbar - transpose(delLbar))/2.d0

        !rotation increment
        delR_var = I-delW/2.d0
        delR_var_inv = 0.d0
        delR_var_det = 0.d0
        call invert_small(delR_var,delR_var_inv,delR_var_det)
        delR = delR_var_inv*(I+delW/2.d0)

        ! pull out state variables for this point
        stress_0(1:6) = initial_state_variables((kint-1)*8+1:(kint-1)*8+6)
        ebar_0 = initial_state_variables((kint-1)*8+7)
        Vf_0 = initial_state_variables(kint*8)
        stress_1 = 0.d0
        ebar_1 = 0.d0
        Vf_1 = 0.d0
        !call stress update subroutine
        call stressupdate(n_properties, element_properties, dele, delR, stress_0, ebar_0, Vf_0, stress_1, ebar_1, Vf_1)

        ! update state vars
        updated_state_variables((kint-1)*8+1:(kint-1)*8+6) = stress_1(1:6)
        updated_state_variables((kint-1)*8+7) = ebar_1
        updated_state_variables(kint*8) = Vf_1

        !check if elt hasn't failed
        if (Vf_1 < fF) then
            element_deleted = .false.
        endif

        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)

        !B-Bar element
        B_aug = 0.d0
        B_aug(1,1:3*n_nodes-2:3) = dNbardy(1:n_nodes,1)-dNdy(1:n_nodes,1)
        B_aug(1,2:3*n_nodes-1:3) = dNbardy(1:n_nodes,2)-dNdy(1:n_nodes,2)
        B_aug(1,3:3*n_nodes:3) = dNbardy(1:n_nodes,3)-dNdy(1:n_nodes,3)
        B_aug(2,1:3*n_nodes-2:3) = dNbardy(1:n_nodes,1)-dNdy(1:n_nodes,1)
        B_aug(2,2:3*n_nodes-1:3) = dNbardy(1:n_nodes,2)-dNdy(1:n_nodes,2)
        B_aug(2,3:3*n_nodes:3) = dNbardy(1:n_nodes,3)-dNdy(1:n_nodes,3)
        B_aug(3,1:3*n_nodes-2:3) = dNbardy(1:n_nodes,1)-dNdy(1:n_nodes,1)
        B_aug(3,2:3*n_nodes-1:3) = dNbardy(1:n_nodes,2)-dNdy(1:n_nodes,2)
        B_aug(3,3:3*n_nodes:3) = dNbardy(1:n_nodes,3)-dNdy(1:n_nodes,3)
        B = B + (1.d0/3.d0)*B_aug
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress_1)*w(kint)*determinant

    end do
  
    return
end subroutine el_dynfrac_dynamic


!==========================SUBROUTINE fieldvars_dynfrac ==============================
subroutine fieldvars_dynfrac(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
    use Element_Utilities, only : dNbardy => vol_avg_shape_function_derivatives_3D
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

    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  p, smises, sdev(6)                ! Pressure and Mises stress
    real (prec)  ::  B(6,length_dof_array)
    real (prec)  ::  B_aug(6,length_dof_array)         ! augmentation matrix for B-bar
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  E, xnu, Y, e0, m, q1, q2, q3      ! Material properties
    real (prec)  ::  fN, eN, sN, fc, fF                ! Material properties (cont)
    real (prec)  ::  eta, deleta                       ! volume averaged things
    real (prec)  ::  delF(3,3), Fmid(3,3)              ! deformation gradient tensors (and determinant)
    real (prec)  ::  Fmid_inv(3,3), J                  ! deformation gradient tensors (and determinant) cont
    real (prec)  ::  vol                               ! element volume
    real (prec)  ::  delL(3,3), delLkk                 ! velocity gradient and it's trace
    real (prec)  ::  dNdy(n_nodes,3)                   ! Shape function derivatives w.r.t deformed coords
    real (prec)  ::  I(3,3)                            ! Identity
    real (prec)  ::  dispmid(length_dof_array)         ! midpoint displacement
    real (prec)  ::  delLbar(3,3)                      ! corrected velocity gradient
    real (prec)  ::  delFlkFklinv                      ! looks like I had a seizure, but I swear I'm ok.
    real (prec)  ::  dele(3,3), delW(3,3), delR(3,3)   ! strain, spin, and rotation increments
    real (prec)  ::  delR_var(3,3), delR_var_inv(3,3)  ! stuff for computing delR
    real (prec)  ::  delR_var_det                      ! stuff for computing delR (cont)
    real (prec)  ::  stress_0(6), stress_1(6)          ! stress before and after step
    real (prec)  ::  ebar_0, ebar_1, Vf_0, Vf_1
    !
    !     Subroutine to compute element force vector for a linear elastodynamic problem
    !     El props are:

    E = element_properties(1)
    xnu = element_properties(2)
    Y = element_properties(3)
    e0 = element_properties(4)
    m = element_properties(5)
    q1 = element_properties(6)
    q2 = element_properties(7)
    q3 = element_properties(8)
    fN = element_properties(9)
    eN = element_properties(10)
    sN = element_properties(11)
    fc = element_properties(12)
    fF = element_properties(13)

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    I = 0.d0
    I(1,1) = 1.d0
    I(2,2) = 1.d0
    I(3,3) = 1.d0

    dispmid = dof_total + .5d0*dof_increment
    vol = 0.d0
    eta = 0.d0
    deleta = 0.d0
    dNbardy = 0.d0
    !     --  Loop over integration points to calculate vol averaged quantities
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        ! calculate delF
        delF = 0.d0
        delF(1,:) = matmul(dof_increment(1:3*n_nodes-2:3),dNdx(1:n_nodes,1:3))
        delF(2,:) = matmul(dof_increment(2:3*n_nodes-1:3),dNdx(1:n_nodes,1:3))
        delF(3,:) = matmul(dof_increment(3:3*n_nodes:3),dNdx(1:n_nodes,1:3))

        ! calculate Fmid
        Fmid = 0.d0
        Fmid_inv = 0.d0
        J = 0.d0
        Fmid(1,:) = matmul(dispmid(1:3*n_nodes-2:3),dNdx(1:n_nodes,1:3))
        Fmid(2,:) = matmul(dispmid(2:3*n_nodes-1:3),dNdx(1:n_nodes,1:3))
        Fmid(3,:) = matmul(dispmid(3:3*n_nodes:3),dNdx(1:n_nodes,1:3))
        Fmid = Fmid + I
        call invert_small(Fmid, Fmid_inv,J)

        delL = matmul(delF,Fmid_inv)
        delLkk = delL(1,1)+delL(2,2)+delL(3,3)

        !element volume contribution
        vol = vol + w(kint)*determinant

        !contributions to volume averaged jacobian and average volumetric strain
        eta = eta + J*w(kint)*determinant
        deleta = deleta + J*delLkk*w(kint)*determinant

        ! calculate dNdy
        dNdy = 0.d0
        dNdy(:,1) = matmul(dNdx(1:n_nodes,:),Fmid_inv(:,1))
        dNdy(:,2) = matmul(dNdx(1:n_nodes,:),Fmid_inv(:,2))
        dNdy(:,3) = matmul(dNdx(1:n_nodes,:),Fmid_inv(:,3))

        !volume averaged spatial shape function derivatives contribution
        dNbardy = dNbardy + dNdy*w(kint)*determinant
    end do

    !calculate volume averages
    eta = eta/vol
    deleta = deleta/(eta*vol)
    dNbardy = dNbardy/(eta*vol)

    !     --  Loop over integration points to calc velocity and rotation and stress
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        ! calculate delF
        delF = 0.d0
        delF(1,:) = matmul(dof_increment(1:3*n_nodes-2:3),dNdx(1:n_nodes,1:3))
        delF(2,:) = matmul(dof_increment(2:3*n_nodes-1:3),dNdx(1:n_nodes,1:3))
        delF(3,:) = matmul(dof_increment(3:3*n_nodes:3),dNdx(1:n_nodes,1:3))

        ! calculate Fmid
        Fmid = 0.d0
        Fmid_inv = 0.d0
        J = 0.d0
        Fmid(1,:) = matmul(dispmid(1:3*n_nodes-2:3),dNdx(1:n_nodes,1:3))
        Fmid(2,:) = matmul(dispmid(2:3*n_nodes-1:3),dNdx(1:n_nodes,1:3))
        Fmid(3,:) = matmul(dispmid(3:3*n_nodes:3),dNdx(1:n_nodes,1:3))
        Fmid(1,1) = Fmid(1,1)+1.d0
        Fmid(2,2) = Fmid(2,2)+1.d0
        Fmid(3,3) = Fmid(3,3)+1.d0
        call invert_small(Fmid, Fmid_inv,J)

        !corrected velocity gradient increment
        delLbar = matmul(delF,Fmid_inv)
        delFlkFklinv = delLbar(1,1) + delLbar(2,2) + delLbar(3,3)
        delLbar = delLbar + I*(deleta-delFlkFklinv)/3.d0

        !strain and spin increments
        dele = (delLbar + transpose(delLbar))/2.d0
        delW = (delLbar - transpose(delLbar))/2.d0

        !rotation increment
        delR_var = I-delW/2.d0
        delR_var_inv = 0.d0
        delR_var_det = 0.d0
        call invert_small(delR_var,delR_var_inv,delR_var_det)
        delR = delR_var_inv*(I+delW/2.d0)

        ! pull out state variables for this point
        stress_0(1:6) = initial_state_variables((kint-1)*8+1:(kint-1)*8+6)
        ebar_0 = initial_state_variables((kint-1)*8+7)
        Vf_0 = initial_state_variables(kint*8)
        stress_1 = 0.d0
        ebar_1 = 0.d0
        Vf_1 = 0.d0
        !call stress update subroutine
        call stressupdate(n_properties, element_properties, dele, delR, stress_0, ebar_0, Vf_0, stress_1, ebar_1, Vf_1)

        ! update state vars

        stress = stress_1/J
        p = sum(stress(1:3))/3.d0
        sdev = stress
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(6)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'Vf',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + Vf_1*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
 
    end do
  
    return
end subroutine fieldvars_dynfrac

subroutine stressupdate(n_properties, element_properties, dstrain, delR, stress_0, ebar_0, Vf_0, stress_1, ebar_1, Vf_1)

    use Types
    use Element_Utilities, only : invert_small
    use Globals, only : TIME, DTIME
    implicit none

    integer, intent( in ) :: n_properties
    real (prec), intent( in ) :: element_properties(n_properties)
    real (prec), intent( in ) :: dstrain(3,3)
    real (prec), intent( in ) :: delR(3,3)
    real (prec), intent( in ) :: stress_0(6) !initial stress [S11, S22, S33, S12, S13, S23]
    real (prec), intent( in ) :: ebar_0
    real (prec), intent( in ) :: Vf_0
    real (prec), intent( out ) :: stress_1(6)
    real (prec), intent( out ) :: ebar_1
    real (prec), intent( out ) :: Vf_1

    !Local variables

    integer      :: iter                                ! counter for N-R loop

    real (prec)  :: deleij(3,3)                         ! deviatoric strain increment
    real (prec)  :: tij(3,3)                            ! deviatoric stress
    real (prec)  :: Smat(3,3)                           ! stress tensor to make eqns simpler
    real (prec)  :: S_star(3,3), p_star, sigma_e_star   ! elastic predictor stuff
    real (prec)  :: E, xnu, Y, e0, m, q1, q2, q3        ! Material properties
    real (prec)  :: fN, eN, sN, fc, fF                  ! Material properties (cont)
    real (prec)  :: delebar                             ! change in ebar
    real (prec)  :: corr, tol                           ! N-R params
    real (prec)  :: delee, delev, ddelee, ddelev        ! solution results for N-R and their increments
    real (prec)  :: se, Dse, p, Dp                      ! expressions and their derivatives for N-R
    real (prec)  :: phi, DphiDse, DphiDp                ! phi and its first derivatives
    real (prec)  :: D2phi, D2phiDse, D2phiDp            ! phi's second derivatives
    real (prec)  :: F1, DF1Ddelee, DF1Ddelev            ! F1 and its derivatives
    real (prec)  :: F2, DF2Ddelee, DF2Ddelev            ! F2 and its derivatives
    real (prec)  :: f_star                              ! stuff for phi
    real (prec)  :: K(2,2), K_inv(2,2), K_det           ! K in N-R
    real (prec)  :: R(2)                                ! R in N-R
    real (prec)  :: dw(2)                               ! N-R increment

    ! begin
    E = element_properties(1)
    xnu = element_properties(2)
    Y = element_properties(3)
    e0 = element_properties(4)
    m = element_properties(5)
    q1 = element_properties(6)
    q2 = element_properties(7)
    q3 = element_properties(8)
    fN = element_properties(9)
    eN = element_properties(10)
    sN = element_properties(11)
    fc = element_properties(12)
    fF = element_properties(13)

    stress_1 = 0.d0
    ! elastic predictors
    deleij = dstrain
    deleij(1,1) = deleij(1,1) - (dstrain(1,1) + dstrain(2,2) + dstrain(3,3))/3.d0
    deleij(2,2) = deleij(2,2) - (dstrain(1,1) + dstrain(2,2) + dstrain(3,3))/3.d0
    deleij(3,3) = deleij(3,3) - (dstrain(1,1) + dstrain(2,2) + dstrain(3,3))/3.d0

    !set up tij
    tij = 0.d0
    tij(1,1) = stress_0(1)
    tij(2,2) = stress_0(2)
    tij(3,3) = stress_0(3)
    tij(1,2) = stress_0(4)
    tij(2,1) = stress_0(4)
    tij(1,3) = stress_0(5)
    tij(3,1) = stress_0(5)
    tij(2,3) = stress_0(6)
    tij(3,2) = stress_0(6)
    tij = tij-sum(stress_0)/3.d0

    S_star(:,:) = E/(1.d0-xnu)*deleij + matmul(matmul(delR,tij),transpose(delR))
    p_star = sum(stress_0)/3.d0 + E/(3.d0*(1.d0-2.d0*xnu))*(dstrain(1,1) + dstrain(2,2) + dstrain(3,3))

    !compute phi
    sigma_e_star = sqrt(3.d0/2.d0)*norm2(S_star)
    if (Vf_0 < fc) then
        f_star = Vf_0
    else if (Vf_0 > fF) then
        f_star = fc + (((q1+sqrt(q1**2-q3))/q3)-fc)/(fF-fc)*(Vf_0-fc)
    else
        f_star = (q1+sqrt(q1**2-q3))/q3
    endif

    phi = sqrt(sigma_e_star**2/Y**2+2.d0*q1*f_star*cosh(3.d0*q2*p_star/(2.d0*Y))-(1.d0+q3*f_star**2))

    !compute stress
    if (phi>1.d-8) then !plastic -- N-R iteration
        !write(*,*) 'plastic increment'
        !initial guesses
        delev = 0.d0
        delee = 0.d0
        corr = 1.d0 !prior to entering loop
        tol = 1.d-7
        iter = 0 !counter

        do while(corr>tol) !N-R loop
            iter = iter + 1

            se = sigma_e_star-(3.d0/2.d0)*E/(1.d0-xnu)*delee !sigma_e and its derivative
            Dse = -(3.d0/2.d0)*E/(1.d0-xnu)

            p = p_star-E/(3.d0*(1.d0-2.d0*xnu))*delev !p and its derivative
            Dp = -E/(3.d0*(1.d0-2.d0*xnu))

            phi = sqrt(se**2/Y**2+2.d0*q1*f_star*cosh(3.d0/2.d0*q2*p/Y)-(1.d0+q3*f_star**2)) !phi and its derivatives
            DphiDse = se/(Y**2*sqrt(-1.d0-f_star**2*q3+se**2/Y**2+2.d0*f_star*q1*cosh((3.d0*p*q2)/(2.d0*Y))))
            DphiDp = (3.d0*f_star*q1*q2*sinh((3.d0*p*q2)/(2.d0*Y)))/ &
                (2.d0*Y*sqrt(-1.d0 - f_star**2*q3 + se**2/Y**2 + 2.d0*f_star*q1*cosh((3.d0*p*q2)/(2.d0*Y))))
            D2phi = (-3.d0*f_star*q1*q2*se*sinh((3.d0*p*q2)/(2.d0*Y)))/ &
                (2.d0*Y**3*(-1.d0 - f_star**2*q3 + se**2/Y**2 + 2.d0*f_star*q1*cosh((3.d0*p*q2)/(2.d0*Y)))**(3.d0/2.d0))
            D2phiDse = -((1.d0 + f_star**2*q3 - 2.d0*f_star*q1*cosh((3.d0*p*q2)/(2.d0*Y)))/ &
                (Y**2*(-1.d0 - f_star**2*q3 + se**2/Y**2 + 2.d0*f_star*q1*cosh((3.d0*p*q2)/(2.d0*Y)))**(3.d0/2.d0)))
            D2phiDp = (9*f_star*q1*q2**2*(2*(se**2 - (1 + f_star**2*q3)*Y**2)*cosh((3*p*q2)/(2*Y)) + &
                f_star*q1*Y**2*(3.d0 + cosh((3.d0*p*q2)/Y))))/ &
                (8.d0*Y**4.d0*(-1.d0 - f_star**2*q3 + se**2/Y**2 + 2.d0*f_star*q1*cosh((3.d0*p*q2)/(2.d0*Y)))**(3.d0/2.d0))

            F1 = sqrt(DphiDse**2+(2.d0/9.d0)*DphiDp**2)*(delee/(DTIME*e0))-DphiDse*phi**m !F1 and its derivatives
            DF1Ddelee = sqrt(2.d0/9.d0*DphiDp**2+DphiDse**2)/(DTIME*e0)- &
                phi**m*Dse*D2phiDse+(delee*Dse*(2.d0*DphiDp*D2phi+9.d0*DphiDse*D2phiDse))/ &
                (3.d0*DTIME*e0*sqrt(2.d0*DphiDp**2+9.d0*DphiDse**2))
            DF1Ddelev = (1.d0/3.d0)*Dp*(-3.d0*phi**m*D2phi+(delee*(2.d0*DphiDp*D2phiDp+ &
                9.d0*DphiDse*D2phi))/(DTIME*e0*sqrt(2.d0*DphiDp**2+9.d0*DphiDse**2)))

            F2 = sqrt(DphiDse**2+(2.d0/9.d0)*DphiDp**2)*(delee/(DTIME*e0))-DphiDp*phi**m !F2 and its derivatives
            DF2Ddelee = sqrt(2.d0/9.d0*DphiDp**2+DphiDse**2)/(DTIME*e0)- &
                phi**m*Dse*D2phi+(delee*Dse*(2.d0*DphiDp*D2phi+9.d0*DphiDse*D2phiDse))/ &
                (3.d0*DTIME*e0*sqrt(2.d0*DphiDp**2+9.d0*DphiDse**2))
            DF2Ddelev = (1.d0/3.d0)*Dp*(-3.d0*phi**m*D2phiDp+(delee*(2.d0*DphiDp*D2phiDp+ &
                9.d0*DphiDse*D2phi))/(DTIME*e0*sqrt(2.d0*DphiDp**2+9.d0*DphiDse**2)))

            K = 0.d0
            K(1,1) = DF1Ddelee
            K(1,2) = DF1Ddelev
            K(2,1) = DF2Ddelee
            K(2,2) = DF2Ddelev
            R(1) = -F1
            R(2) = -F2

            call invert_small(K,K_inv,K_det)

            dw = matmul(K_inv,R)
            ddelee = dw(1)
            ddelev = dw(2)

            corr = sqrt(ddelee**2+ddelev**2)
            delee = delee + ddelee
            delev = delev + ddelev
        end do

        Smat = S_star-delee*(E/(1.d0+xnu))*(3.d0/2.d0)*(S_star/sigma_e_star)
        stress_1 = 0.d0
        stress_1(1) = Smat(1,1)
        stress_1(2) = Smat(2,2)
        stress_1(3) = Smat(3,3)
        stress_1(4) = Smat(1,2)
        stress_1(5) = Smat(1,3)
        stress_1(6) = Smat(2,3)
        stress_1(1:3) = stress_1(1:3) + (p_star-E/(3.d0*(1.d0-2.d0*xnu))*delev)

        delebar = e0*DTIME/(1-Vf_0)*phi**m*sqrt(DphiDse**2+(2.d0/9.d0)*DphiDp**2)*(se*DphiDse+(1.d0/3.d0)*DphiDp*p)
        ebar_1 = ebar_0+delebar
        Vf_1 = 1.d0+(Vf_0-1.d0)*exp(-delev)+fN*delebar/(sN*sqrt(2.d0*pi))*exp(-(1.d0/2.d0)*((ebar_0-eN)/sN)**2)


    else !elastic-simple
        !write(*,*) 'elastic increment'
        stress_1 = 0.d0
        stress_1(1) = S_star(1,1)
        stress_1(2) = S_star(2,2)
        stress_1(3) = S_star(3,3)
        stress_1(4) = S_star(1,2)
        stress_1(5) = S_star(1,3)
        stress_1(6) = S_star(2,3)
        stress_1(1:3) = stress_1(1:3) + p_star
        delebar = 0.d0
        Vf_1 = 0.d0
        ebar_1 = 0.d0

    endif

    return
end subroutine stressupdate
