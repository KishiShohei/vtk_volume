program main
    use mesh_m
    implicit none
    type(UnstructuredGrid_inVTK) grid
    
    call grid%read_vtk()
    call grid%calculate_volume()
    call grid%output_vtk()
end program