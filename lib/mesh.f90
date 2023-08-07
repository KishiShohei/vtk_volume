module mesh_m
    implicit none

    type cell_t
        integer, allocatable :: nodeID(:)
        integer n_TYPE
    end type

    type node_t
        real coordinate(3)
    end type

    type UnstructuredGrid_inVTK
        type(cell_t), allocatable :: cells(:)
        type(node_t), allocatable :: nodes(:)

        contains

        procedure read_vtk
        procedure calculate_volume
        procedure output_vtk

    end type

    contains

    subroutine read_vtk(self)
        class(UnstructuredGrid_inVTK) self
        integer n_unit, KKMX, KK, IIMX, II, num_node, l, nodeID(8)
        character(17) AAA
        character(6) BBB
        character(99) str

        open(newunit = n_unit, file="data/old_ver.vtk", status="old")
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,'()')
            read(n_unit,*) AAA,KKMX
            allocate(self%nodes(KKMX))
            do KK=1,KKMX
                read(n_unit,*) self%nodes(KK)%coordinate(:)
            end do

            read(n_unit,*) BBB,IIMX
            allocate(self%cells(IIMX))
            do II = 1, IIMX
                read(n_unit, '(A)') str
                read(str, *) num_node
                read(str, *) num_node, (nodeID(l), l=1,num_node)
                self%cells(II)%nodeID = nodeID(:num_node) + 1
            end do

            read(n_unit, '()')
            do II = 1, IIMX
                read(n_unit, *) self%cells(II)%n_TYPE
            end do

        close(n_unit)

        ! print*, "IIMX=", IIMX

    end subroutine


    subroutine calculate_volume(self)
        class(UnstructuredGrid_inVTK) self
        real V_tetra, V_hexa, V_tr_prism, V_quad_prism, A, B, C, D, E, F, G, H, I, J, K, L

        A=norm2(self%nodes(self%cells(1)%nodeID(2))%coordinate - self%nodes(self%cells(1)%nodeID(1))%coordinate)
        B=norm2(self%nodes(self%cells(1)%nodeID(3))%coordinate - self%nodes(self%cells(1)%nodeID(1))%coordinate)
        C=norm2(self%nodes(self%cells(1)%nodeID(4))%coordinate - self%nodes(self%cells(1)%nodeID(1))%coordinate)
        V_tetra = A*B*C/(2*3)
        print*, "V_tetra =",V_tetra
        D=norm2(self%nodes(self%cells(2)%nodeID(2))%coordinate - self%nodes(self%cells(2)%nodeID(1))%coordinate)
        E=norm2(self%nodes(self%cells(2)%nodeID(4))%coordinate - self%nodes(self%cells(2)%nodeID(1))%coordinate)
        F=norm2(self%nodes(self%cells(2)%nodeID(5))%coordinate - self%nodes(self%cells(2)%nodeID(1))%coordinate)
        V_hexa = D*E*F
        print*, "V_hexa = ",V_hexa
        G=norm2(self%nodes(self%cells(3)%nodeID(2))%coordinate - self%nodes(self%cells(3)%nodeID(1))%coordinate)
        H=norm2(self%nodes(self%cells(3)%nodeID(3))%coordinate - self%nodes(self%cells(3)%nodeID(1))%coordinate)
        I=norm2(self%nodes(self%cells(3)%nodeID(4))%coordinate - self%nodes(self%cells(3)%nodeID(1))%coordinate)
        V_tr_prism = G*H*I/2
        print*, "V_tr_prism = ",V_tr_prism
        J=norm2(self%nodes(self%cells(4)%nodeID(2))%coordinate - self%nodes(self%cells(4)%nodeID(1))%coordinate)
        K=norm2(self%nodes(self%cells(4)%nodeID(4))%coordinate - self%nodes(self%cells(4)%nodeID(1))%coordinate)
        L=norm2(self%nodes(self%cells(4)%nodeID(5))%coordinate - self%nodes(self%cells(4)%nodeID(1))%coordinate)
        V_quad_prism = J*K*L/3
        print*, "V_quad_prism = ", V_quad_prism
    end subroutine

    subroutine output_vtk(self)
        class(UnstructuredGrid_inVTK) self
        integer KK, KKMX, II, IIMX, IITOTAL, n_unit
        integer, allocatable :: nodeID(:)

        open(newunit = n_unit, file ="data/output.vtk", status = "replace")
            write(n_unit, '(A)') '# vtk DataFile Version 2.0'
            write(n_unit, '(A)') 'Header'
            write(n_unit, '(A)') 'ASCII'
            write(n_unit, '(A)') 'DATASET UNSTRUCTURED_GRID'

            KKMX = size(self%nodes)
            write(n_unit, '(*(g0:," "))') 'POINTS', KKMX, 'float'
            do KK = 1, KKMX
                write(n_unit, '(*(g0:," "))') self%nodes(KK)%coordinate(:)
            end do

            IIMX = size(self%cells)
            IITOTAL = 0
            do II = 1, IIMX
                IITOTAL = IITOTAL + size(self%cells(II)%nodeID) + 1
            end do
            write(n_unit, '(*(g0:," "))') 'CELLS', IIMX, IITOTAL
            do II = 1, IIMX
                nodeID = self%cells(II)%nodeID - 1
                write(n_unit, '(*(g0:," "))') size(nodeID), nodeID
            end do

            write(n_unit, '(*(g0:," "))') 'CELL_TYPES', IIMX
            do II = 1, IIMX
                write(n_unit, '(*(g0:," "))') self%cells(II)%n_TYPE
            enddo
        close(n_unit)
    end subroutine

end module