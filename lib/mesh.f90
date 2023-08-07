module mesh_m
    implicit none
    private

    type cell_inVTK_t
        integer, private ::n_TYPE
        integer, allocatable :: nodeID(:)
    end type 

    type node_t
    real coordinate(3)
    end type

    type, public :: UnstructuredGrid_inVTK
        type(node_t), allocatable :: node_array(:)
        type(cell_inVTK_t), allocatable :: cell_array(:)

        contains

        procedure :: read_UnstructuredGrid_inVTK
        procedure :: output_UnstructuredGrid_inVTK
        procedure get_numCell, get_numNode, get_nodeCoordinate, get_cellVertices
        procedure set_nodeCoordinate, set_cellVertices
        procedure calculate_volume

    end type


    contains



    subroutine read_UnstructuredGrid_inVTK(self, FNAME, action, cellScalar, cellVector)
        class(UnstructuredGrid_inVTK) self
        character(*), intent(in) :: FNAME
        character(*), intent(in), optional ::  action
        real, allocatable, intent(out), optional :: cellScalar(:), cellVector(:,:)
        integer II,KK, num_node, n_unit, KKMX, IIMX, ios
        integer l, nodeID(8)
        character AAA*7, str*99
        logical dataOnly

        dataOnly = .false.
        if(present(action)) then
            if(action=='dataOnly') dataOnly = .true.
        end if

        print*, 'READ_VTK:', FNAME

        open(newunit=n_unit,FILE=FNAME, status='old', action='read')
            if(dataOnly) then
                do while(index(AAA, 'CELL_DATA')==0)
                    read(n_unit, '(A)') AAA
                end do
            else
                read(n_unit,'()')
                read(n_unit,'()')
                read(n_unit,'()')
                read(n_unit,'()')
                read(n_unit,*) AAA,KKMX
                allocate(self%node_array(KKMX))
                DO KK = 1, KKMX
                    read(n_unit,*) self%node_array(KK)%coordinate(:)
                END DO
                read(n_unit, *) AAA,IIMX
                allocate(self%cell_array(IIMX))
                DO II = 1, IIMX
                    read(n_unit, '(A)') str !一旦1行丸ごと読む
                    read(str, *) num_node
                    read(str, *) num_node, (nodeID(l), l=1,num_node)
                    self%cell_array(II)%nodeID = nodeID(:num_node) + 1
                END DO

                read(n_unit,'()')  !CELL_TYPES
                DO II = 1, IIMX
                    read(n_unit, *) self%cell_array(II)%n_TYPE
                END DO
                if((.not.present(cellScalar)) .and. (.not.present(cellVector))) return
                read(n_unit,'()')
    
                read(n_unit,'()')   !CELL_DATA

            end if 

    end subroutine        

    subroutine calculate_volume(self)
        class(UnstructuredGrid_inVTK) self
        real t1,t2,t3,h1,h2,h3,pr1,pr2,pr3,py1,py2,py3
        real volume_tetra , volume_hexa , volume_prism , volume_pyramid
        t1 = self%node_array(self%cell_array(1)%nodeID(2))%coordinate(1) &
        - self%node_array(self%cell_array(1)%nodeID(1))%coordinate(1)
        t2 = self%node_array(self%cell_array(1)%nodeID(3))%coordinate(2) &
        - self%node_array(self%cell_array(1)%nodeID(1))%coordinate(2)
        t3 = self%node_array(self%cell_array(1)%nodeID(4))%coordinate(3) &
        - self%node_array(self%cell_array(1)%nodeID(1))%coordinate(3)
        volume_tetra = t1*t2*t3/6
        h1 = self%node_array(self%cell_array(2)%nodeID(2))%coordinate(1) &
        - self%node_array(self%cell_array(2)%nodeID(1))%coordinate(1)
        h2 = self%node_array(self%cell_array(2)%nodeID(4))%coordinate(2) &
        - self%node_array(self%cell_array(2)%nodeID(1))%coordinate(2)
        h3 = self%node_array(self%cell_array(2)%nodeID(5))%coordinate(3) &
        - self%node_array(self%cell_array(2)%nodeID(1))%coordinate(3)
        volume_hexa = h1*h2*h3 
        pr1 = self%node_array(self%cell_array(3)%nodeID(3))%coordinate(2) &
        - self%node_array(self%cell_array(3)%nodeID(1))%coordinate(2)
        pr2 = self%node_array(self%cell_array(3)%nodeID(4))%coordinate(3) &
        - self%node_array(self%cell_array(3)%nodeID(1))%coordinate(3)
        pr3 = self%node_array(self%cell_array(3)%nodeID(2))%coordinate(1) &
        - self%node_array(self%cell_array(3)%nodeID(1))%coordinate(1)
        volume_prism = pr1*pr2*pr3/2
        py1 = self%node_array(self%cell_array(4)%nodeID(2))%coordinate(1) &
        - self%node_array(self%cell_array(4)%nodeID(1))%coordinate(1)
        py2 = self%node_array(self%cell_array(4)%nodeID(4))%coordinate(2) &
        - self%node_array(self%cell_array(4)%nodeID(1))%coordinate(2)
        py3 = self%node_array(self%cell_array(4)%nodeID(5))%coordinate(3) &
        - self%node_array(self%cell_array(4)%nodeID(1))%coordinate(3)
        volume_pyramid = py1*py2*py3/3


        print*,"volume_tetra=", volume_tetra
        print*,"volume_hexa=", volume_hexa
        print*,"volume_prism=", volume_prism
        print*,"volume_pyramid=", volume_pyramid

    end subroutine

    subroutine output_UnstructuredGrid_inVTK(self, FNAME, cellScalar, cellVector, scalarName, vectorName)
        class(UnstructuredGrid_inVTK) self
        character(*), intent(in) :: FNAME
        real, intent(in), optional :: cellScalar(:), cellVector(:,:)
        character(*), intent(in), optional :: scalarName, vectorName
        integer II,KK, n_unit, KKMX, IIMX, IITOTAL
        integer, allocatable :: nodeID(:)

        print*, 'OUTPUT_VTK:', FNAME
            
        open(newunit=n_unit, FILE=FNAME, STATUS='replace')
            write(n_unit, '(A)') '# vtk DataFile Version 2.0'
            write(n_unit, '(A)') 'Header'
            write(n_unit, '(A)') 'ASCII'
            write(n_unit, '(A)') 'DATASET UNSTRUCTURED_GRID'

            KKMX = size(self%node_array)
            write(n_unit, '(A,1x,I0,1x,A)') 'POINTS', KKMX, 'float'   
            DO KK = 1, KKMX
                write(n_unit,'(3(e12.5,2X))') self%node_array(KK)%coordinate(:)
            END DO
            write(n_unit,'()')

            IIMX = size(self%cell_array)
            IITOTAL = 0
            do II = 1, IIMX
                IITOTAL = IITOTAL + size(self%cell_array(II)%nodeID)  + 1
            end do
            write(n_unit,'(A,I0,2X,I0)') 'CELLS ', IIMX, IITOTAL
            DO II = 1, IIMX
                nodeID = self%cell_array(II)%nodeID(:) - 1
                write(n_unit, '(*(g0:," "))') size(nodeID), nodeID(:)
            END DO
            write(n_unit,'()')

            write(n_unit,'(A,I0)') 'CELL_TYPES ', IIMX
            DO II = 1, IIMX
                write(n_unit, '(I0)') self%cell_array(II)%n_TYPE
            END DO
    
        end subroutine 

        integer function get_numNode(self)
            class(UnstructuredGrid_inVTK), intent(in) :: self
            get_numNode = size(self%node_array)
        end function
    
        integer function get_numCell(self)
            class(UnstructuredGrid_inVTK), intent(in) :: self
            get_numCell = size(self%cell_array)
        end function
    
        function get_nodeCoordinate(self) result(cdn)
            class(UnstructuredGrid_inVTK), intent(in) :: self
            real, allocatable :: cdn(:,:)
            integer k, num_node
    
            num_node = size(self%node_array)
            allocate(cdn(3, num_node))
    
            do k = 1, num_node
                cdn(:,k) = self%node_array(k)%coordinate(:)
            end do
    
        end function
    
        subroutine get_cellVertices(self, vertices, types)
            class(UnstructuredGrid_inVTK), intent(in) :: self
            integer, allocatable, intent(out) :: vertices(:,:), types(:)
            integer i, num_cell, num_node
    
            num_cell = size(self%cell_array)
            allocate(vertices(6, num_cell))
            allocate(types(num_cell))
    
            do i = 1, num_cell
                num_node = size(self%cell_array(i)%nodeID)
                vertices(1:num_node, i) = self%cell_array(i)%nodeID(:)
                types(i) = self%cell_array(i)%n_TYPE
            end do
    
        end subroutine
    
        subroutine set_nodeCoordinate(self, cdn)
            class(UnstructuredGrid_inVTK) self
            real, intent(in) :: cdn(:,:)
            integer k, num_node
    
            num_node = size(cdn, dim=2)
            allocate(self%node_array(num_node))
    
            do k = 1, num_node
                self%node_array(k)%coordinate(:) = cdn(:,k) 
            end do
    
        end subroutine
    
        subroutine set_cellVertices(self, vertices, types)
            class(UnstructuredGrid_inVTK) self
            integer, intent(in) :: vertices(:,:), types(:)
            integer i, num_cell, num_node
    
            num_cell = size(vertices, dim=2)
            allocate(self%cell_array(num_cell))
    
            do i = 1, num_cell
                select case(types(i))
                    case(10)!tetra
                        num_node = 4
                    case(11)!hexa
                        num_node = 8
                    case(13)!prism
                        num_node = 6
                    case(14)!pyramid
                        num_node = 5
                end select
                self%cell_array(i)%nodeID = vertices(1:num_node, i)
                self%cell_array(i)%n_TYPE = types(i)
            end do
    
        end subroutine   
        
end module 