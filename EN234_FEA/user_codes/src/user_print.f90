subroutine user_print(n_steps)
  use Types
  use ParamIO
  use Globals, only : TIME, DTIME
  use Mesh, only : extract_element_data
  use Mesh, only : extract_node_data
  use Mesh, only : zone,zone_list
  use Printparameters, only : n_user_print_files                  ! No. files specified by the user
  use Printparameters, only : n_user_print_parameters             ! No. user supplied print parameters
  use Printparameters, only : user_print_units                    ! Unit numbers
  use Printparameters, only : user_print_parameters               ! List of user supplied parameters
  use User_Subroutine_Storage, only : length_state_variable_array ! Max no. state variables on any element
  implicit none
  
  integer, intent(in) :: n_steps                                 ! Current step number
  
  integer ::  lmn, lmn_start, lmn_end
  integer ::  status
  integer ::  n_state_vars_per_intpt                                         ! No. state variables per integration point
  integer :: iforce, ninc, idisp, istrain
  real (prec) ::   vol_averaged_strain(6)                                    ! Volume averaged strain in an element
!  real (prec), allocatable ::   vol_averaged_state_variables(:)              ! Volume averaged state variables in an element
!  real (prec) :: J_value                                                     !J integral value
  real (prec) :: vol_averaged_stress(6)                                                   ! volume averaged stress (for hypoelastic)
  real (prec) :: Dof1(6), Dof2(6)
  real (prec) :: xs(11), Forces(11), Disps(11), Strains(11)


!
!  Use this file to process or print time histories of the solution, or to print a non-standard mesh.
!
!  As an example, this subroutine computes the volume averaged infinitesimal strain and the volume average
!  element state variables (if they exist) in an element.   The element is specified by the user.
!  The first six state variables (which are usually the stresses) are printed along with the strains.
!
!

!   allocate(vol_averaged_state_variables(length_state_variable_array), stat=status)

!   if (status/=0) then
!      write(IOW,*) ' Error in subroutine user_print'
!      write(IOW,*) ' Unable to allocate memory for state variables '
!      stop
!   endif

!   lmn = int(user_print_parameters(1))     ! The element number

!   call compute_element_volume_average_3D(lmn,vol_averaged_strain,vol_averaged_state_variables,length_state_variable_array, &
!                                                       n_state_vars_per_intpt)



!    if (TIME<1.d-12) then
!      if (n_state_vars_per_intpt<6) then
!        write(user_print_units(1),'(A)') 'VARIABLES = TIME,e11,e22,e33,e12,e13,e23'
!      else
!         write(user_print_units(1),'(A)') 'VARIABLES = TIME,e11,e22,e33,e12,e13,e23,s11,s22,s33,s12,s13,s23'
!      endif
!    endif

!   if (n_state_vars_per_intpt<6) then
!      write(user_print_units(1),'(7(1x,D12.5))') TIME+DTIME,vol_averaged_strain(1:6)
!   else
!      vol_averaged_state_variables(1:3) = vol_averaged_state_variables(1:3) + vol_averaged_state_variables(7)
!      write(user_print_units(1),'(13(1x,D12.5))') TIME+DTIME,vol_averaged_strain(1:6),vol_averaged_state_variables(1:6)
!   endif

!!! J integral calculation
!    call compute_J_integral(J_value)
 !   write(user_print_units(1),'(A)') 'J-Integral Value:'
 !   write(user_print_units(1),*) J_value

 !! Hypoelastic vol averaged stress/strain
 !   lmn = 1
 !   vol_averaged_strain = 0.d0
 !   vol_averaged_stress = 0.d0
 !   call compute_hypo_stress(lmn,vol_averaged_strain, vol_averaged_stress)
 !   !write(user_print_units(1),'(A)') 'Volume Averaged Strain'
 !   !write(user_print_units(1),*) vol_averaged_strain
 !   !write(user_print_units(1),'(A)') 'Volume Averaged Stress'
 !   !write(user_print_units(1),*) vol_averaged_stress
 !   write(user_print_units(1),*) vol_averaged_strain(1), ',', vol_averaged_strain(2), ',', vol_averaged_stress(1), ';'

 !! Timoshenko Stress Output
  !  You can access the first and last element using
    lmn_start = zone_list(1)%start_element
    lmn_end = zone_list(1)%end_element
    write(user_print_units(1), '(A)') 'Time:'
    write(user_print_units(1), *) TIME
    iforce = user_print_parameters(1)
    idisp = user_print_parameters(2)
    istrain = user_print_parameters(3)
    ninc = 10
    do lmn = lmn_start,lmn_end ! loop through elements
        Dof1 = 0.d0
        Dof2 = 0.d0
        xs = 0.d0
        Forces = 0.d0
        Strains = 0.d0
        call Timoshenko_Print(lmn, Dof1, Dof2, xs, Forces, Disps, Strains, ninc, iforce, idisp, istrain)
        write(user_print_units(1), '(A)') '    Element Number: '
        write(user_print_units(1), *) '    ', lmn
        write(user_print_units(1), '(A)') '    Dof1 - [u1,v1,w1,phi1,psi1,theta1,...]'
        write(user_print_units(1), *) '    ', Dof1
        write(user_print_units(1), '(A)') '    Dof2 - [u1,v1,w1,phi1,psi1,theta1,...]'
        write(user_print_units(1), *) '    ', Dof2
        !clunkiest way to do this, but whatever
        write(user_print_units(1), *) '     x positions (element cs only, matlab readable): [', xs, '];'
        write(user_print_units(1), *) '     Forces (matlab readable): [', Forces, '];'
        write(user_print_units(1), *) '     Displacements (matlab readable): [', Disps, '];'
        write(user_print_units(1), *) '     Strains (matlab readable): [', Strains, '];'

    end do

end subroutine user_print

subroutine compute_element_volume_average_3D(lmn,vol_averaged_strain,vol_averaged_state_vars,length_output_array, &
                                                                                                       n_state_vars_per_intpt)
    use Types
    use ParamIO
    use Mesh, only : extract_element_data
    use Mesh, only : extract_node_data
    use User_Subroutine_Storage
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent ( in )      :: lmn                                          ! Element number
    integer, intent ( in )      :: length_output_array

    real (prec), intent( out )  ::  vol_averaged_strain(6)
    real (prec), intent( out )  ::  vol_averaged_state_vars(length_output_array)

    integer, intent( out )      :: n_state_vars_per_intpt

    ! Local variables

    integer    :: node_identifier                              ! Flag identifying node type
    integer    :: element_identifier                           ! Flag identifying element type (specified in .in file)
    integer    :: n_nodes                                      ! # nodes on the element
    integer    :: n_properties                                 ! # properties for the element
    integer    :: n_state_variables                            ! # state variables for the element


    integer, allocatable    :: node_list(:)                                ! List of nodes on the element (connectivity)

    real( prec ), allocatable   :: element_properties(:)                  ! Element or material properties, stored in order listed in input file
    real( prec ), allocatable   :: initial_state_variables(:)             ! Element state variables at start of step.
    real( prec ), allocatable   :: updated_state_variables(:)             ! State variables at end of time step

    real( prec ), allocatable   :: x(:,:)                                  ! Nodal coords x(i,a) is ith coord of ath node
    real( prec ), allocatable   :: dof_increment(:)                        ! DOF increment, using usual element storage convention
    real( prec ), allocatable   :: dof_total(:)                            ! accumulated DOF, using usual element storage convention

    integer      :: n_points,kint,i
    integer      :: n_coords, n_dof
    integer      :: iof
    integer      :: status

    real (prec)  ::  el_vol
    real (prec), allocatable  ::  B(:,:)               ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  strain(6)                         ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  dstrain(6)                        ! Strain increment vector
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    !
    !  Allocate memory to store element data.
    !  The variables specifying the size of the arrays are stored in the module user_subroutine_storage
    !  They are initialized when the input file is read, and specify the size of the arrays required to store data
    !  for any element in the mesh.  Some elements may require less storage.

    allocate(node_list(length_node_array), stat=status)
    allocate(element_properties(length_property_array), stat=status)
    allocate(initial_state_variables(length_state_variable_array), stat=status)
    allocate(updated_state_variables(length_state_variable_array), stat=status)
    allocate(x(3,length_coord_array/3), stat=status)
    allocate(dof_increment(length_dof_array), stat=status)
    allocate(dof_total(length_dof_array), stat=status)
    allocate(B(6,length_dof_array), stat=status)

    if (status/=0) then
       write(IOW,*) ' Error in subroutine compute_volume_average_3D'
       write(IOW,*) ' Unable to allocate memory for element variables '
       stop
    endif
    !
    ! Extract element and node data from global storage (see module Mesh.f90 for the source code for these subroutines)

    call extract_element_data(lmn,element_identifier,n_nodes,node_list,n_properties,element_properties, &
                                            n_state_variables,initial_state_variables,updated_state_variables)

    do i = 1, n_nodes
        iof = 3*(i-1)+1     ! Points to first DOF for the node in the dof_increment and dof_total arrays
        call extract_node_data(node_list(i),node_identifier,n_coords,x(1:3,i),n_dof, &
                                                 dof_increment(iof:iof+2),dof_total(iof:iof+2))
    end do

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    vol_averaged_strain = 0.d0
    vol_averaged_state_vars = 0.d0
    el_vol = 0.d0
    n_state_vars_per_intpt = n_state_variables/n_points

    if (n_state_vars_per_intpt>size(vol_averaged_state_vars)) then
       write(IOW,*) ' Error detected in subroutine compute_element_volume_average_3d '
       write(IOW,*) ' The element contains ',n_state_vars_per_intpt
       write(IOW,*) ' but the array storing averaged state variables has length ',size(vol_averaged_state_vars)
       stop
    endif
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)

        iof = n_state_vars_per_intpt*(kint-1)+1
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

        strain = matmul(B(1:6,1:3*n_nodes),dof_total(1:3*n_nodes))
        dstrain = matmul(B(1:6,1:3*n_nodes),dof_increment(1:3*n_nodes))

        vol_averaged_strain(1:6) = vol_averaged_strain(1:6) + (strain(1:6)+dstrain(1:6))*w(kint)*determinant

        if (n_state_vars_per_intpt>0) then
           vol_averaged_state_vars(1:n_state_vars_per_intpt) = vol_averaged_state_vars(1:n_state_vars_per_intpt) &
                              + updated_state_variables(iof:iof+n_state_vars_per_intpt-1)*w(kint)*determinant
        endif

        el_vol = el_vol + w(kint)*determinant

    end do

    vol_averaged_strain = vol_averaged_strain/el_vol
    vol_averaged_state_vars = vol_averaged_state_vars/el_vol

    deallocate(node_list)
    deallocate(element_properties)
    deallocate(initial_state_variables)
    deallocate(updated_state_variables)
    deallocate(x)
    deallocate(dof_increment)
    deallocate(dof_total)
    deallocate(B)

    return




end subroutine compute_element_volume_average_3D

subroutine compute_J_integral(J_integral_value)
    use Types
    use ParamIO
    use Mesh, only : extract_element_data
    use Mesh, only : extract_node_data
    use Mesh, only : zone,zone_list
    use User_Subroutine_Storage
    use Element_Utilities, only : N => shape_functions_2D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_2D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_2D
    use Element_Utilities, only : xi => integrationpoints_2D, w => integrationweights_2D
    use Element_Utilities, only : dxdxi => jacobian_2D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    real (prec), intent( out )  ::  J_integral_value


    ! Local variables

    integer    :: node_identifier                              ! Flag identifying node type
    integer    :: element_identifier                           ! Flag identifying element type (specified in .in file)
    integer    :: n_nodes                                      ! # nodes on the element
    integer    :: n_properties                                 ! # properties for the element
    integer    :: n_state_variables                            ! # state variables for the element
    integer    :: n_coords                                     ! No. coords for a node
    integer    :: n_dof                                        ! No. DOFs for a node

    integer      :: status
    integer      :: iof
    integer      :: lmn               ! Element number
    integer      :: lmn_start,lmn_end ! First and last crack tip element
    integer      :: i                 ! Loop counter

!   The arrays below have to be given dimensions large enough to store the data. It doesnt matter if they are too large.

    integer, allocatable    :: node_list(:)                                ! List of nodes on the element (connectivity)

    real( prec ), allocatable   :: element_properties(:)                  ! Element or material properties, stored in order listed in input file
    real( prec ), allocatable   :: initial_state_variables(:)             ! Element state variables at start of step.
    real( prec ), allocatable   :: updated_state_variables(:)             ! State variables at end of time step

    real( prec ), allocatable   :: x(:,:)                                  ! Nodal coords x(i,a) is ith coord of ath node
    real( prec ), allocatable   :: dof_increment(:)                        ! DOF increment, using usual element storage convention
    real( prec ), allocatable   :: dof_total(:)                            ! accumulated DOF, using usual element storage convention

    real (prec), allocatable  ::  B(:,:)                                   ! strain = B*(dof_total+dof_increment)

    ! Local Variables
    integer      :: n_points,kint

    real (prec)  ::  strain(3), dstrain(3)             ! Strain vector contains [e11, e22, 2e12]
    real (prec)  :: totstrain(3)                       ! sum of strain and dstrain
    real (prec)  ::  stress(3)                         ! Stress vector contains [s11, s22, s12]
    real (prec)  ::  D(3,3)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
    real (prec)  :: ptx(2), r                          ! Coords of integration pts
    real (prec)  :: r0                                 ! constant in J-integral
    real (prec)  :: SEdensity, du1dx2                  ! strain energy density and du1/dx2
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    !
    !
    !  The variables specifying the sizes of the arrays listed below are determined while reading the input file
    !  They specify the array dimensions required to store the relevant variables for _any_ element or node in the mesh
    !  The actual data will vary depending on the element or node selected
    !
    allocate(node_list(length_node_array), stat=status)
    allocate(element_properties(length_property_array), stat=status)
    allocate(initial_state_variables(length_state_variable_array), stat=status)
    allocate(updated_state_variables(length_state_variable_array), stat=status)
    allocate(x(2,length_coord_array/2), stat=status)
    allocate(dof_increment(length_dof_array), stat=status)
    allocate(dof_total(length_dof_array), stat=status)
    allocate(B(3,length_dof_array), stat=status)

  !  Write your code to calculate the J integral here

  !  You will need to loop over the crack tip elements, and sum the contribution to the J integral from each element.
  !
  !  You can access the first and last crack tip element using
    lmn_start = zone_list(2)%start_element
    lmn_end = zone_list(2)%end_element

  !  The two subroutines below extract data for elements and nodes (see module Mesh.f90 for the source code for these subroutines)


    !!! begin my code

    J_integral_value = 0.d0
    r0 = .0006d0

    do lmn = lmn_start,lmn_end ! loop through elements
        call extract_element_data(lmn,element_identifier,n_nodes,node_list,n_properties,element_properties, &
                                                n_state_variables,initial_state_variables,updated_state_variables)


        do i = 1, n_nodes
            iof = 2*(i-1)+1     ! Points to first DOF for the node in the dof_increment and dof_total arrays
            call extract_node_data(node_list(i),node_identifier,n_coords,x(1:2,i),n_dof, &
                                                 dof_increment(iof:iof+2),dof_total(iof:iof+2))
        end do



        n_points = 9 ! for all crack tip elements

        call initialize_integration_points(n_points, n_nodes, xi, w)
        !element_residual = 0.d0
        !element_stiffness = 0.d0


        D = 0.d0
        E = element_properties(1)
        xnu = element_properties(2)
        d44 = 0.5D0*E/(1+xnu)
        d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
        d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
        D(1:3,1:3) = d12
        D(1,1) = d11
        D(2,2) = d11
        D(3,3) = d44
        D(1:2,3) = 0.d0
        D(3,1:2) = 0.d0

        !     --  Loop over integration points
        do kint = 1, n_points
            call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
            dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
            call invert_small(dxdxi,dxidx,determinant)
            dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)
            B = 0.d0
            B(1,1:2*n_nodes-1:2) = dNdx(1:n_nodes,1)
            B(2,2:2*n_nodes:2) = dNdx(1:n_nodes,2)
            B(3,1:2*n_nodes-1:2)   = dNdx(1:n_nodes,2)
            B(3,2:2*n_nodes-1:2) = dNdx(1:n_nodes,1)

            strain = matmul(B,dof_total)
            dstrain = matmul(B,dof_increment)
            totstrain = strain+dstrain
            stress = matmul(D,totstrain)

            ptx(1:2) = matmul(x(1:2,1:n_nodes),N(1:n_nodes))
            r = sqrt(ptx(1)**2+ptx(2)**2)
            SEdensity = (totstrain(1)*stress(1) + totstrain(2)*stress(2) + .5*totstrain(3)*stress(3))/2
            du1dx2 = DOT_PRODUCT(B(2,2:2*n_nodes:2),dof_total(1:2*n_nodes-1:2))
            J_integral_value = J_integral_value - (1/(r*r0))*((ptx(1)*stress(1)+ptx(2)*stress(3))*du1dx2 &
                + (ptx(1)*stress(3)+ptx(2)*stress(2))*totstrain(2)-ptx(2)*SEdensity)*w(kint)*determinant
        end do
    end do

    !!! end my code

    deallocate(node_list)
    deallocate(element_properties)
    deallocate(initial_state_variables)
    deallocate(updated_state_variables)
    deallocate(x)
    deallocate(dof_increment)
    deallocate(dof_total)
    deallocate(B)


    return




end subroutine compute_J_integral

subroutine compute_hypo_stress(lmn, vol_averaged_strain, vol_averaged_stress)
    use Types
    use ParamIO
    use Mesh, only : extract_element_data
    use Mesh, only : extract_node_data
    use Mesh, only : zone,zone_list
    use User_Subroutine_Storage
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent ( in )      :: lmn                                          ! Element number

    real (prec), intent( out )  ::  vol_averaged_strain(6)
    real (prec), intent( out )  ::  vol_averaged_stress(6)


    ! Local variables

    integer    :: node_identifier                              ! Flag identifying node type
    integer    :: element_identifier                           ! Flag identifying element type (specified in .in file)
    integer    :: n_nodes                                      ! # nodes on the element
    integer    :: n_properties                                 ! # properties for the element


    integer, allocatable    :: node_list(:)                                ! List of nodes on the element (connectivity)

    real( prec ), allocatable   :: element_properties(:)                  ! Element or material properties, stored in order listed in input file

    real( prec ), allocatable   :: initial_state_variables(:)             ! Element state variables at start of step.
    real( prec ), allocatable   :: updated_state_variables(:)             ! State variables at end of time step
    integer    :: n_state_variables                                         ! # state variables for the element

    real( prec ), allocatable   :: x(:,:)                                  ! Nodal coords x(i,a) is ith coord of ath node
    real( prec ), allocatable   :: dof_increment(:)                        ! DOF increment, using usual element storage convention
    real( prec ), allocatable   :: dof_total(:)                            ! accumulated DOF, using usual element storage convention

    integer      :: n_points,kint,i
    integer      :: n_coords, n_dof
    integer      :: iof
    integer      :: status

    real (prec)  ::  el_vol
    real (prec), allocatable  ::  B(:,:)               ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  strain(6)                         ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  dstrain(6)                        ! Strain increment vector
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  :: s0, e0, Kmat, nmat                       ! Material properties
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  B_aug(6,length_dof_array)         ! array to augment B in the case of B-bar elements
    !
    !  Allocate memory to store element data.
    !  The variables specifying the size of the arrays are stored in the module user_subroutine_storage
    !  They are initialized when the input file is read, and specify the size of the arrays required to store data
    !  for any element in the mesh.  Some elements may require less storage.

    allocate(node_list(length_node_array), stat=status)
    allocate(element_properties(length_property_array), stat=status)
    allocate(initial_state_variables(length_state_variable_array), stat=status)
    allocate(updated_state_variables(length_state_variable_array), stat=status)
    allocate(x(3,length_coord_array/3), stat=status)
    allocate(dof_increment(length_dof_array), stat=status)
    allocate(dof_total(length_dof_array), stat=status)
    allocate(B(6,length_dof_array), stat=status)

    if (status/=0) then
       write(IOW,*) ' Error in subroutine compute_volume_average_3D'
       write(IOW,*) ' Unable to allocate memory for element variables '
       stop
    endif
    !
    ! Extract element and node data from global storage (see module Mesh.f90 for the source code for these subroutines)

    call extract_element_data(lmn,element_identifier,n_nodes,node_list,n_properties,element_properties, &
                                            n_state_variables,initial_state_variables,updated_state_variables)

    do i = 1, n_nodes
        iof = 3*(i-1)+1     ! Points to first DOF for the node in the dof_increment and dof_total arrays
        call extract_node_data(node_list(i),node_identifier,n_coords,x(1:3,i),n_dof, &
                                                 dof_increment(iof:iof+2),dof_total(iof:iof+2))
    end do

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)


    D = 0.d0


    ! need to define volume averages for b-bar -- cannot do inline with B calculation :(
    dNbardx = 0.d0
    el_vol = 0.d0
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        dNbardx = dNbardx + dNdx*w(kint)*determinant
        el_vol = el_vol + w(kint)*determinant
    end do

    dNbardx = (1/el_vol)*dNbardx


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

        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
        s0 = element_properties(1)
        e0 = element_properties(2)
        nmat = element_properties(3)
        Kmat = element_properties(4)

        call calculate_hypoelast_d(strain+dstrain, s0, e0, nmat, Kmat, D, stress)

        vol_averaged_strain(1:6) = vol_averaged_strain(1:6) + (strain(1:6)+dstrain(1:6))*w(kint)*determinant
        vol_averaged_stress(1:6) = vol_averaged_stress(1:6) + stress(1:6)*w(kint)*determinant

    end do

    vol_averaged_strain = vol_averaged_strain/el_vol
    vol_averaged_stress = vol_averaged_stress/el_vol

    deallocate(node_list)
    deallocate(element_properties)
    deallocate(initial_state_variables)
    deallocate(updated_state_variables)
    deallocate(x)
    deallocate(dof_increment)
    deallocate(dof_total)
    deallocate(B)

    return
end subroutine compute_hypo_stress

subroutine Timoshenko_Print(lmn, Dof1, Dof2, xs, Forces, Disps, Strains, ninc, iforce, idisp, istrain)
    use Types
    use ParamIO
    use Mesh, only : extract_element_data
    use Mesh, only : extract_node_data
    use Mesh, only : zone,zone_list
    use User_Subroutine_Storage
    implicit none

    integer, intent( in )  ::  lmn
    integer, intent( in )  ::  ninc, iforce, idisp, istrain !number of increments and which force to output
    real (prec), intent( out )  :: Dof1(6), Dof2(6)
    real (prec), intent( out )  :: xs(ninc), Forces(ninc), Disps(ninc), Strains(ninc) !internal forces and the increments they're at, etc


    ! Local variables

    integer    :: node_identifier                              ! Flag identifying node type
    integer    :: element_identifier                           ! Flag identifying element type (specified in .in file)
    integer    :: n_nodes                                      ! # nodes on the element
    integer    :: n_properties                                 ! # properties for the element
    integer    :: n_state_variables                            ! # state variables for the element
    integer    :: n_coords                                     ! No. coords for a node
    integer    :: n_dof                                        ! No. DOFs for a node

    integer      :: status
    integer      :: iof
    integer      :: lmn_start,lmn_end ! First and last crack tip element
    integer      :: i                 ! Loop counter

!   The arrays below have to be given dimensions large enough to store the data. It doesnt matter if they are too large.

    integer, allocatable    :: node_list(:)                                ! List of nodes on the element (connectivity)

    real( prec ), allocatable   :: element_properties(:)                  ! Element or material properties, stored in order listed in input file
    real( prec ), allocatable   :: initial_state_variables(:)             ! Element state variables at start of step.
    real( prec ), allocatable   :: updated_state_variables(:)             ! State variables at end of time step

    real( prec ), allocatable   :: x(:,:)                                  ! Nodal coords x(i,a) is ith coord of ath node
    real( prec ), allocatable   :: dof_increment(:)                        ! DOF increment, [u1,v1,w1,phi1,psi1,theta1,...]
    real( prec ), allocatable   :: dof_total(:)                            ! accumulated DOF, same convention as increment

    real (prec), allocatable  ::  B(:,:)                                   ! strain = B*(dof_total+dof_increment)

    ! Local Variables
    integer      :: n_points,kint

    real (prec)  ::  intForce(6)                       ! [N, Mz, My, Qy, Qz, Mx]
    real (prec)  ::  strain(6), dstrain(6)             ! generalized strains [e, ky, kz, gy, gz, kx]
    real (prec)  ::  D(6,6)                            ! Material matrix
    real (prec)  ::  N(18), dNdx(18)                   ! array of shape functions and derivatives [N1, N2, Hv1, Hv2, Hw1, Hw2, Ht1, Ht2, Hp1, Hp2, Gv1, Gv2, Gw1, Gw2, Gt1, Gt2, Gp1, Gp2]
    real (prec)  ::  dxidx, det, w, xi                 ! values to calculate
    real (prec)  ::  E, G, A, k, Iz, Iy                ! Material properties
    real (prec)  ::  L                                 ! Calculated Length
    real (prec)  :: ay, by, az, bz                     ! constants in shape functions
    real (prec)  :: F(6,12), dispInt(6)                ! Matrix to interpolate DOFs
    !
    !
    !  The variables specifying the sizes of the arrays listed below are determined while reading the input file
    !  They specify the array dimensions required to store the relevant variables for _any_ element or node in the mesh
    !  The actual data will vary depending on the element or node selected
    !
    allocate(node_list(length_node_array), stat=status)
    allocate(element_properties(length_property_array), stat=status)
    allocate(initial_state_variables(length_state_variable_array), stat=status)
    allocate(updated_state_variables(length_state_variable_array), stat=status)
    allocate(x(3,length_coord_array/3), stat=status)
    allocate(dof_increment(length_dof_array), stat=status)
    allocate(dof_total(length_dof_array), stat=status)
    allocate(B(6,12), stat=status)

  !  The two subroutines below extract data for elements and nodes (see module Mesh.f90 for the source code for these subroutines)


    !!! begin my code



    call extract_element_data(lmn,element_identifier,n_nodes,node_list,n_properties,element_properties, &
                                           n_state_variables,initial_state_variables,updated_state_variables)


    do i = 1, n_nodes
        iof = 6*(i-1)+1     ! Points to first DOF for the node in the dof_increment and dof_total arrays
        call extract_node_data(node_list(i),node_identifier,n_coords,x(1:3,i),n_dof, &
                                                 dof_increment(iof:iof+5),dof_total(iof:iof+5))
    end do

    Dof1 = dof_total(1:6) + dof_increment(1:6)
    Dof2 = dof_total(7:12) + dof_increment(7:12)

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
    w = 1.d0
    dxidx = 1.d0/L !note difference from typical because of luo's convention of xi = 0->1
    det = L/2.d0

    do kint = 0,ninc
    ! calculate shape functions and derivatives
    N = 0.d0
    dNdx = 0.d0
    xi = real(kint)/real(ninc)
    xs(kint+1) = xi*L

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

    !populate F matrix
    F = 0.d0
    F(1,1:7:6) = N(1:2)
    F(2,2:8:6) = N(3:4)
    F(2,6:12:6) = N(7:8)
    F(3,3:9:6) = N(5:6)
    F(3,5:11:6) = N(9:10)
    F(4,4:10:6) = N(1:2)
    F(5,3:9:6) = N(13:14)
    F(5,5:11:6) = N(17:18)
    F(6,2:8:6) = N(11:12)
    F(6,6:12:6) = N(15:16)


    strain = matmul(B,dof_total)
    dstrain = matmul(B,dof_increment)
    intForce = matmul(D,strain+dstrain)
    Forces(kint+1) = intForce(iforce)
    dispInt = matmul(F,dof_total+dof_increment)
    Disps(kint+1) = dispInt(idisp)
    Strains(kint+1) = strain(istrain)+dstrain(istrain)

    end do

    !!! end my code

    deallocate(node_list)
    deallocate(element_properties)
    deallocate(initial_state_variables)
    deallocate(updated_state_variables)
    deallocate(x)
    deallocate(dof_increment)
    deallocate(dof_total)
    deallocate(B)


    return




end subroutine Timoshenko_Print

