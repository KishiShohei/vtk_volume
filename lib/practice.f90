module mesh_practice
    implicit none

    type cell_t
    !セル情報を書き込む
    !必要な情報:nodeの番号、セルタイプ
        integer, allocatable :: nodeID(:)
        integer n_TYPE
    end type

    type node_t
    !セル一つに対する情報を書き込む
    !必要な情報:座標
        real coordinate(3)
    end type

    type UnstructuredGrid_inVTK
    !クラス情報。手続きを含む。
    ! 今回はread_vtk,calculate_volume,output_vtkの3つとする。
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
        integer KKMX, IIMX, KK, II, num_node, l, nodeID(8)
        character(17) AAA
        character(6) BBB
        character(99) str

        open(newunit = n_unit, file = "old_ver.vtk" status = "old")
            read(n_unit, '()')
            read(n_unit, '()')
            read(n_unit, '()')
            read(n_unit, '()')
            read(n_unit,*) AAA, KKMX
            allocate(self%nodes(KKMX))
            do KK = 1, KKMX
                read(n_unit,*) self%nodes(KK)%coordinate(:)
            enddo

            read(n_unit,*) BBB, IIMX
            allocate(self%cells(IIMX))
            do II = 1, IIMX
                read(n_unit, '(A)') str
                read(str, *) num_node
                read(str, *) num_node, (nodeID(l), l=1,num_node)
                self%cells(II)%nodeID = nodeID(:num_node) + 1
            enddo

            read(n_unit, '()')
            do II = 1, IIMX
                read(n_unit, *) self%cells(II)%n_TYPE
            end do

            close(n_unit)
    end subroutine

    subroutine calculate_volume()
    end subroutine

    subroutine output_vtk()
    end subroutine
end module