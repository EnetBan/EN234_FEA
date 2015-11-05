!     Subroutines for basic 3D linear elastic elements 



!==========================SUBROUTINE el_hyperelast_3d ==============================
subroutine el_hyperelast_3d(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
    integer      :: n_points,kint, ind

    real (prec)  ::  disp(3*n_nodes)                   ! displacement vector contains [u11, u12, u13, u21, etc]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  stressk(6)                        ! kirchoff stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  vol                               ! volume of element
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  mu1, K1                           ! Material properties
    real (prec)  ::  F(3,3), Finv(3,3), J              ! Deformation gradient matrix, its inverse, and its determinant
    real (prec)  ::  BC(3,3), BCinv(3,3), BCdet, Bkk   ! Left C-G Tensor and its inverse, determinant, and trace
    real (prec)  ::  BCv(6), BCinvv(6), Iv(6)          ! vector forms of BC and its inverse and the identity matrix [B11,B22,B33,B12,B13,B23]
    real (prec)  ::  dNdy(n_nodes,3)                   ! Shape function derivatives w.r.t deformed coords
    real (prec)  ::  M1(6,6)                           ! Matrix for computing D
    real (prec)  ::  G(6,9)                            ! G matrix for stiffness
    real (prec)  ::  Bstar(9, 3*n_nodes)               ! B* matrix
    real (prec)  ::  sig(3*n_nodes,3*n_nodes)          ! sigma matrix
    real (prec)  ::  S(3,n_nodes), Pvec(3*n_nodes)     ! various matrices for computing sigma
    real (prec)  ::  Pmat(3*n_nodes,3*n_nodes), Svec(3*n_nodes)
    real (prec)  ::  Smat(3*n_nodes,3*n_nodes)
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D hyperelastic elements
    !     El props are:

    !     element_properties(1)         Shear Modulus
    !     element_properties(2)         Bulk Modulus

    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0
	
    mu1 = element_properties(1)
    K1 = element_properties(2)

    disp = dof_total+dof_increment
    dNdx = 0.d0

    !     --  Loop over integration points
    do kint = 1, n_points
        ! calculate shapefunction derivatives
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        ! calculate Fij
        F = 0.d0
        Finv = 0.d0
        J = 0.d0
        F(1,:) = matmul(disp(1:3*n_nodes-2:3),dNdx(1:n_nodes,1:3))
        F(2,:) = matmul(disp(2:3*n_nodes-1:3),dNdx(1:n_nodes,1:3))
        F(3,:) = matmul(disp(3:3*n_nodes:3),dNdx(1:n_nodes,1:3))
        F(1,1) = F(1,1)+1.d0
        F(2,2) = F(2,2)+1.d0
        F(3,3) = F(3,3)+1.d0
        call invert_small(F, Finv, J)                         ! determinant and inverse
        ! calculate left C-G tensor
        BC = 0.d0
        BCinv = 0.d0
        BCdet = 0.d0
        BC = matmul(F,transpose(F))
        call invert_small(BC, BCinv, BCdet)
        Bkk = BC(1,1) + BC(2,2) + BC(3,3)
        BCv(1) = BC(1,1)
        BCv(2) = BC(2,2)
        BCv(3) = BC(3,3)
        BCv(4) = BC(1,2)
        BCv(5) = BC(1,3)
        BCv(6) = BC(2,3)
        BCinvv(1) = BCinv(1,1)
        BCinvv(2) = BCinv(2,2)
        BCinvv(3) = BCinv(3,3)
        BCinvv(4) = BCinv(1,2)
        BCinvv(5) = BCinv(1,3)
        BCinvv(6) = BCinv(2,3)
        Iv = 0.d0
        Iv(1:3) = 1

        ! calculate dNdy
        dNdy = 0.d0
        dNdy(:,1) = matmul(dNdx(1:n_nodes,:),Finv(:,1))
        dNdy(:,2) = matmul(dNdx(1:n_nodes,:),Finv(:,2))
        dNdy(:,3) = matmul(dNdx(1:n_nodes,:),Finv(:,3))

        !calculate kirchoff stress
        stressk = 0.d0
        stressk(1:3) = (mu1/(J**2)**(1.d0/3.d0))*(BCv(1:3)-Bkk/3.d0)+J*K1*(J-1.d0)
        stressk(4:6) = (mu1/(J**2)**(1.d0/3.d0))*BCv(4:6)

        ! B is in relation to y
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

        ! construct D
        D = 0.d0
        M1 = 0.d0
        M1(1,1) = 1.d0
        M1(2,2) = 1.d0
        M1(3,3) = 1.d0
        M1(4,4) = .5d0
        M1(5,5) = .5d0
        M1(6,6) = .5d0

        D = (mu1/(J**2)**(1.d0/3.d0))*M1+(mu1/(3*(J**2)**(1.d0/3.d0)))*(Bkk/3 &
            * spread(Iv,dim=2,ncopies=6)*spread(BCinvv,dim=1,ncopies=6) &
            - spread(Iv,dim=2,ncopies=6)*spread(Iv,dim=1,ncopies=6) &
            - spread(BCv,dim=2,ncopies=6)*spread(BCinvv,dim=1,ncopies=6)) &
            + K1*J*(J-.5d0)*(spread(Iv,dim=2,ncopies=6)*spread(Bcinvv,dim=1,ncopies=6))

        ! fill G
        G = 0.d0
        G(1,1) = 2*BC(1,1)
        G(1,4) = 2*BC(1,2)
        G(1,6) = 2*BC(1,3)
        G(2,2) = 2*BC(2,2)
        G(2,5) = 2*BC(1,3)
        G(2,8) = 2*BC(2,3)
        G(3,3) = 2*BC(3,3)
        G(3,7) = 2*BC(1,3)
        G(3,9) = 2*BC(1,3)
        G(4,1:2) = 2*BC(1,2)
        G(4,4) = 2*BC(2,2)
        G(4,5) = 2*BC(1,1)
        G(4,6) = 2*BC(2,3)
        G(4,8) = 2*BC(1,3)
        G(5,1) = 2*BC(1,3)
        G(5,3) = 2*BC(1,3)
        G(5,4) = 2*BC(2,3)
        G(5,6) = 2*BC(3,3)
        G(5,7) = 2*BC(1,1)
        G(5,9) = 2*BC(1,2)
        G(6,2:3) = 2*BC(2,3)
        G(6,5) = 2*BC(1,3)
        G(6,7) = 2*BC(1,2)
        G(6,8) = 2*BC(3,3)
        G(6,9) = 2*BC(2,2)

        ! construct B*
        Bstar = 0.d0
        Bstar(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        Bstar(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        Bstar(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        Bstar(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        Bstar(5,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        Bstar(6,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        Bstar(7,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        Bstar(8,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        Bstar(9,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)

        ! construct sigma
        sig = 0.d0
        Pmat = 0.d0
        Smat = 0.d0
        S = reshape(matmul(transpose(B),stressk),(/3,length_dof_array/3/))
        ind = 0
        do ind = 1,n_nodes
            Pvec = reshape(spread(transpose(dNdy(ind:ind,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
            Pmat(3*ind-2:3*ind,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
            Svec = reshape(spread(S(1:3,ind),dim=2,ncopies=n_nodes),(/3*n_nodes/))
            Smat(3*ind-2:3*ind,1:3*n_nodes) = spread(Svec,dim=1,ncopies=3)
        end do
        sig = Pmat*transpose(Smat)

        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stressk)*w(kint)*determinant
        !write(IOW,*) ' '
        !write(IOW,*) 'disp', dof_total
        !write(IOW,*) 'dispinc', dof_increment
        !write(IOW,*) 'F', F(1,:)
        !write(IOW,*) '  ', F(2,:)
        !write(IOW,*) '  ', F(3,:)
        !write(IOW,*) 'BC', BC(1,:)
        !write(IOW,*) '  ', BC(2,:)
        !write(IOW,*) '  ', BC(3,:)
        !write(IOW,*) 'dNdy', dNdy
        !write(IOW,*) 'J', J
        !write(IOW,*) 'stressk', stressk
        !write(IOW,*) 'term1', (mu1/(J**2)**(1.d0/3.d0))*(BCv(1:3)-Bkk/3.d0)
        !write(IOW,*) 'term2', J*K1*(J-1.d0)
        !write(IOW,*) 'B', B
        !write(IOW,*) 'D', D
        !write(IOW,*) 'G', G
        !write(IOW,*) 'B*', Bstar
        !write(IOW,*) 'sigma', sig(1,:)
        !write(IOW,*) 'resid', element_residual(1:8)
        !write(IOW,*) '    ', element_residual(9:16)
        !write(IOW,*) '    ', element_residual(17:24)
        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            + (matmul(transpose(B),matmul(D,matmul(G,Bstar)))-sig) &
            *w(kint)*determinant
    end do

    return
end subroutine el_hyperelast_3d


!==========================SUBROUTINE el_hyperelast_3d_dynamic ==============================
subroutine el_hyperelast_3d_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
               
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          
    logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete element

    ! Local Variables
    integer      :: n_points,kint

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    !
    !     Subroutine to compute element force vector for a linear elastodynamic problem
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
	
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

        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
      
        stress = matmul(D,strain+dstrain)
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

    end do
  
    return
end subroutine el_hyperelast_3d_dynamic


!==========================SUBROUTINE fieldvars_hyperelast_3d ==============================
subroutine fieldvars_hyperelast_3d(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
  
    integer      :: n_points,kint,k, ind
    real (prec)  ::  disp(3*n_nodes)                   ! displacement vector contains [u11, u12, u13, u21, etc]
    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  stressk(6)                        ! kirchoff stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  vol                               ! volume of element
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  mu1, K1                           ! Material properties
    real (prec)  ::  F(3,3), Finv(3,3), J              ! Deformation gradient matrix, its inverse, and its determinant
    real (prec)  ::  BC(3,3), BCinv(3,3), BCdet, Bkk   ! Left C-G Tensor and its inverse, determinant, and trace
    real (prec)  ::  BCv(6), BCinvv(6), Iv(6)          ! vector forms of BC and its inverse and the identity matrix [B11,B22,B33,B12,B13,B23]
    real (prec)  ::  dNdy(n_nodes,3)                   ! Shape function derivatives w.r.t deformed coords
    real (prec)  ::  M1(6,6)                           ! Matrix for computing D
    real (prec)  ::  G(6,9)                            ! G matrix for stiffness
    real (prec)  ::  Bstar(9, 3*n_nodes)               ! B* matrix
    real (prec)  ::  sig(3*n_nodes,3*n_nodes)        ! sigma matrix
    real (prec)  ::  S(3,n_nodes), Pvec(3*n_nodes) ! various matrices for computing sigma
    real (prec)  ::  Pmat(3*n_nodes,3*n_nodes), Svec(3*n_nodes)
    real (prec)  ::  Smat(3*n_nodes,3*n_nodes)
    real (prec)  :: p, smises                          ! Pressure and Mises stress
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D hyperelastic elements
    !     El props are:

    !     element_properties(1)         Shear Modulus
    !     element_properties(2)         Bulk Modulus

    x = reshape(element_coords,(/3,length_coord_array/3/))
    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    mu1 = element_properties(1)
    K1 = element_properties(2)

    disp = dof_total+dof_increment
    dNdx = 0.d0

    !     --  Loop over integration points
    do kint = 1, n_points
        ! calculate shapefunction derivatives
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        ! calculate Fij
        F = 0.d0
        Finv = 0.d0
        J = 0.d0
        F(1,:) = matmul(disp(1:3*n_nodes-2:3),dNdx(1:n_nodes,1:3))
        F(2,:) = matmul(disp(2:3*n_nodes-1:3),dNdx(1:n_nodes,1:3))
        F(3,:) = matmul(disp(3:3*n_nodes:3),dNdx(1:n_nodes,1:3))
        F(1,1) = F(1,1)+1.d0
        F(2,2) = F(2,2)+1.d0
        F(3,3) = F(3,3)+1.d0
        call invert_small(F, Finv, J)                         ! determinant and inverse
        ! calculate left C-G tensor
        BC = 0.d0
        BCinv = 0.d0
        BCdet = 0.d0
        BC = matmul(F,transpose(F))
        call invert_small(BC, BCinv, BCdet)
        Bkk = BC(1,1) + BC(2,2) + BC(3,3)
        BCv(1) = BC(1,1)
        BCv(2) = BC(2,2)
        BCv(3) = BC(3,3)
        BCv(4) = BC(1,2)
        BCv(5) = BC(1,3)
        BCv(6) = BC(2,3)
        BCinvv(1) = BCinv(1,1)
        BCinvv(2) = BCinv(2,2)
        BCinvv(3) = BCinv(3,3)
        BCinvv(4) = BCinv(1,2)
        BCinvv(5) = BCinv(1,3)
        BCinvv(6) = BCinv(2,3)
        Iv = 0.d0
        Iv(1:3) = 1

        ! calculate dNdy
        dNdy = 0.d0
        dNdy(:,1) = matmul(dNdx(1:n_nodes,:),Finv(:,1))
        dNdy(:,2) = matmul(dNdx(1:n_nodes,:),Finv(:,2))
        dNdy(:,3) = matmul(dNdx(1:n_nodes,:),Finv(:,3))

        !calculate kirchoff stress
        stressk = 0.d0
        stressk(1:3) = (mu1/J**(2.d0/3.d0))*(BCv(1:3)-Bkk/3.d0)+J*K1*(J-1.d0)
        stressk(4:6) = (mu1/J**(2.d0/3.d0))*BCv(4:6)

        ! B is in relation to y
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

        ! construct D
        D = 0.d0
        M1 = 0.d0
        M1(1,1) = 1.d0
        M1(2,2) = 1.d0
        M1(3,3) = 1.d0
        M1(4,4) = .5d0
        M1(5,5) = .5d0
        M1(6,6) = .5d0

        D = (mu1/J**(2.d0/3.d0))*M1+(mu1/(3*J**(2.d0/3.d0)))*(Bkk/3 &
            * spread(Iv,dim=2,ncopies=6)*spread(BCinvv,dim=1,ncopies=6) &
            - spread(Iv,dim=2,ncopies=6)*spread(Iv,dim=1,ncopies=6) &
            - spread(BCv,dim=2,ncopies=6)*spread(BCinvv,dim=1,ncopies=6)) &
            + K1*J*(J-.5d0)*(spread(Iv,dim=2,ncopies=6)*spread(Bcinvv,dim=1,ncopies=6))

        ! fill G
        G = 0.d0
        G(1,1) = 2*BC(1,1)
        G(1,4) = 2*BC(1,2)
        G(1,6) = 2*BC(1,3)
        G(2,2) = 2*BC(2,2)
        G(2,5) = 2*BC(1,3)
        G(2,8) = 2*BC(2,3)
        G(3,3) = 2*BC(3,3)
        G(3,7) = 2*BC(1,3)
        G(3,9) = 2*BC(1,3)
        G(4,1:2) = 2*BC(1,2)
        G(4,4) = 2*BC(2,2)
        G(4,5) = 2*BC(1,1)
        G(4,6) = 2*BC(2,3)
        G(4,8) = 2*BC(1,3)
        G(5,1) = 2*BC(1,3)
        G(5,3) = 2*BC(1,3)
        G(5,4) = 2*BC(2,3)
        G(5,6) = 2*BC(3,3)
        G(5,7) = 2*BC(1,1)
        G(5,9) = 2*BC(1,2)
        G(6,2:3) = 2*BC(2,3)
        G(6,5) = 2*BC(1,3)
        G(6,7) = 2*BC(1,2)
        G(6,8) = 2*BC(3,3)
        G(6,9) = 2*BC(2,2)

        ! construct B*
        Bstar = 0.d0
        Bstar(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        Bstar(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        Bstar(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        Bstar(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        Bstar(5,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        Bstar(6,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        Bstar(7,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        Bstar(8,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        Bstar(9,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)

        ! construct sigma
        sig = 0.d0
        Pmat = 0.d0
        Smat = 0.d0
        S = reshape(matmul(transpose(B),stressk),(/3,length_dof_array/3/))
        ind = 0
        do ind = 1,n_nodes
            Pvec = reshape(spread(transpose(dNdy(ind:ind,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
            Pmat(3*ind-2:3*ind,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
            Svec = reshape(spread(S(1:3,ind),dim=2,ncopies=n_nodes),(/3*n_nodes/))
            Smat(3*ind-2:3*ind,1:3*n_nodes) = spread(Svec,dim=1,ncopies=3)
        end do
        sig = Pmat*transpose(Smat)

        !strain = matmul(B,dof_total)
        !dstrain = matmul(B,dof_increment)
        !stress = matmul(D,strain+dstrain)

        p = sum(stressk(1:3))/3.d0
        sdev = stressk
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stressk(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stressk(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stressk(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stressk(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stressk(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stressk(6)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
 
    end do
  
    return
end subroutine fieldvars_hyperelast_3d

