module mesh_practice
    implicit none

    type cell_t
    !セル情報を書き込む
    !必要な情報:nodeの番号、セルタイプ
    end type

    type node_t
    !セル一つに対する情報を書き込む
    !必要な情報:座標
    end type

    type UnstructuredGrid_inVTK
    !クラス情報。手続きを含む。
    ! 今回はread_vtk,calculate_volume,output_vtkの3つとする。
    end type

    contains

    subroutine read_vtk()
    end subroutine

    subroutine calculate_volume()
    end subroutine

    subroutine output_vtk()
    end subroutine
end module