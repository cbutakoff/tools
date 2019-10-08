#cython: language_level=3
import numpy as np
cimport numpy as np
from libc.stdint cimport int64_t, int32_t

def element_list_CGNS2ANSABIN_int64( np.ndarray[np.int64_t, ndim=1] data, int element_id_shift, int nelements, int generate_element_offsets=0 ):
    # data is a np.array of the type [cgns_element_type n1 n2 .. cgns_element_type n1 n2 .. cgns_element_type n1 n2 .. ]
    # generate_element_offsets creates another array of pairs (element_id, element_index), such that result[element_index] points to the element element_id 

    cells_new = np.zeros( nelements+data.shape[0], dtype=np.uint64 )
    element_offsets = None

    if generate_element_offsets==1:
        element_offsets = np.zeros( (nelements,2), dtype=np.uint64 )

	# CGNS element types
	# typedef enum {
	# ElementTypeNull, ElementTypeUserDefined,	/* 0, 1,	*/
	# NODE, BAR_2, BAR_3, 				/* 2, 3, 4, 	*/
	# TRI_3, TRI_6,					/* 5, 6,	*/
	# QUAD_4, QUAD_8, QUAD_9,				/* 7, 8, 9,	*/
	# TETRA_4, TETRA_10, 				/* 10, 11,	*/
	# PYRA_5, PYRA_14, 				/* 12, 13,	*/
	# PENTA_6, PENTA_15, PENTA_18,			/* 14, 15, 16,	*/
	# HEXA_8, HEXA_20, HEXA_27, 			/* 17, 18, 19,	*/
	# MIXED, PYRA_13, NGON_n, NFACE_n			/* 20, 21, 22, 23*/
	# } ElementType_t;

    cgns_nnodes_per_element = np.zeros(18, dtype=np.int8)
    cgns_nnodes_per_element[ [5,7,10,12,14,17] ] = [3,4,4,5,6,8]  #only chosen elements
    		
    cdef int64_t i_old = 0
    cdef int64_t i_new = 0
    cdef int64_t cell_id = 1 + element_id_shift
    cdef int nnodes = 0
    cdef int64_t i = 0
    while i_old<data.shape[0]:
        nnodes = <int>cgns_nnodes_per_element[ data[i_old] ]
        cells_new[ i_new ] = cell_id
        cells_new[ i_new + 1 ] = nnodes

        cells_new[ (i_new + 2):(i_new+2+nnodes) ] = data[ (i_old + 1):(i_old+1+nnodes) ]

        if generate_element_offsets==1:
            element_offsets[i,:] = [cell_id, i_new]

        i_old += nnodes+1
        i_new += nnodes+2
        cell_id += 1
        i += 1

    return cells_new, element_offsets

