#include <stdio.h>
#include <stdlib.h>

#include "grid_io.h"

int writeTimetraces2ascii( int dim0, int dim1, int t_end, double period, 
                           char filename[], double timetraces[dim0][dim1] ) {
//{{{

    size_t
        ii;

    FILE
        *file_pntr;

    // open file in w(rite) mode; might consider using a+ instead
    file_pntr   = fopen( filename, "w" );
    if (file_pntr == NULL) {
        printf( "ERROR: Unable to create file for timetraces.\n" );
        return EXIT_FAILURE;
    } else {
        // NOTE: if return value of printf < 0, then writing failed.
        //       might be good idea to implicetely check this
        //       e.g. if ( (fprintf( file_pntr, "a b c" )) < 0 ) ....
        fprintf( file_pntr, "# T  poynt_z1  poynt_z2  poynt_x1  poynt_x2  poynt_y1  poynt_y2  P_out\n" ); 
        for ( ii=0 ; ii<(t_end/(int)period) ; ++ii )
            fprintf( file_pntr, " %4d  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e\n",
                    (int)timetraces[ii][1], 
                    timetraces[ii][2], timetraces[ii][3],
                    timetraces[ii][4], timetraces[ii][5],
                    timetraces[ii][6], timetraces[ii][7],
                    (timetraces[ii][2]+timetraces[ii][3] + timetraces[ii][4]+timetraces[ii][5] + timetraces[ii][6]+timetraces[ii][7])
                  );
        if ((fclose(file_pntr)) == EOF) {
            printf( "ERROR: could not close file for timetraces.\n" );
        }
    }
    printf( "successfully written timetraces into %s\n", filename );
    return EXIT_SUCCESS;

}//}}}


#ifdef HDF5
int writeMyHDF_v4( int dim0, int dim1, int dim2, char filename[], char dataset[], double array_3D[dim0][dim1][dim2] ) {
    //#{{{

    // hdf related variables
    hid_t       file_id, dataset_id, dataspace_id;      // object identifiers
    hsize_t     dims[3];                                // size used for dimensions
    herr_t      status;                                 // function return value

    // hdf5 related variables used for applying shuffle and compression filter
    hid_t       dcpl;
    hsize_t     chunk[3];
    unsigned int    filter_info;
    int         filter_avail;

    // required for check if hdf5-file already exists
    struct      stat st;

    // assume as default setting, that filters are available
    filter_avail = 1;

    // check if specified hdf5 file already exists
    // if not, create new one; if yes, open and add dataset to it
    if ( stat( filename, &st )==0 ) {
        // open file for read + write access
        file_id = H5Fopen( filename,            // filename
                           H5F_ACC_RDWR,        // allow read & write access (_RDONLY for read only)
                           H5P_DEFAULT);        // file access property list (default one)

        // hdf5 version 1.8.0 introduced H5Lexists to check if link (to group or dataset) exists in hdf5-file
#if H5_VERS_MAJOR>=1 && H5_VERS_MINOR>=8
        if ( H5_VERS_MINOR >= 10 ) {
            printf( "WARNING: hdf5 version 1.10 (or larger is used)\n" );
            printf( "         behavior of H5Lexists was slightly changed in this version\n" );
            printf( "         for details, see https://support.hdfgroup.org/HDF5/doc/RM/RM_H5L.html#Link-Exists\n" );
        }
        if ( H5Lexists( file_id,                // file or group identifier
                        dataset,                // name of link (to group or dataset) to check
                        H5P_DEFAULT )           // link access property list identifiert
                > 0 ) {                         // NOTE: for version 1.8.10, this might be slightly different
            printf( "ERROR: dataset named '%s' already exists in file '%s'\n", dataset, filename );
            printf( "       dataset will NOT be saved (no overwrite by default)\n" );
            status = H5Fclose(file_id);
            if (status < 0) printf( "ERROR: could not close file '%s'\n", filename );
            return EXIT_FAILURE;
        }
#endif
    } else {
        // create a new file using default properties.
        file_id = H5Fcreate( filename,          // filename
                             H5F_ACC_TRUNC,     // how file should be created (removes existing file)
                             H5P_DEFAULT,       // file creating property list
                             H5P_DEFAULT);      // file access property list
    }

    // create simple data space for the dataset
    // (simple = regular N-dimensional array, i.e. data on regular rectangular grid)
    // (complex = e.g.: irregular grids resulting from dynamic refinement of mesh)
    dims[0] = dim0;
    dims[1] = dim1;
    dims[2] = dim2;
    dataspace_id = H5Screate_simple( 3,     // number of dimensions of dataspace 
                                     dims,  // size of array in each dimension
                                     NULL); // allow stretching of data space (NULL=no)

    // check if filters for shuffle and gzip exists, if yes, apply them
    // (check is done, because the filters are optional parts of hdf5 library)
    // check if gzip-filter is available
    if ( !(H5Zfilter_avail( H5Z_FILTER_DEFLATE )) ) {
        printf( "WARNING: gzip filter not available (for hdf5)\n" );
        filter_avail = 0;
    } else {
        status = H5Zget_filter_info( H5Z_FILTER_DEFLATE, &filter_info );
        if (status < 0) printf( "ERROR: could not get hdf5 filter info\n" );
        if ( !(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED) ||
             !(filter_info & H5Z_FILTER_CONFIG_DECODE_ENABLED) ) {
            printf( "WARNING: gzip filter not available for encoding and decoding (for hdf5)\n" );
            filter_avail = 0;
        }
    }
    // check if shuffle-filter is available
    if ( !(H5Zfilter_avail( H5Z_FILTER_SHUFFLE )) ) {
        printf( "WARNING: shuffle filter not available (for hdf5)\n" );
        filter_avail = 0;
    } else {
        status = H5Zget_filter_info( H5Z_FILTER_SHUFFLE, &filter_info );
        if (status < 0) printf( "ERROR: could not get hdf5 filter info\n" );
        if ( !(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED) ||
             !(filter_info & H5Z_FILTER_CONFIG_DECODE_ENABLED) ) {
            printf( "WARNING: shuffle filter not available for encoding and decoding (for hdf5)\n" );
            filter_avail = 0;
        }
    }

    // apply shuffle and gzip filters, if available
    if (filter_avail) {
        // set chunk size to be same as dimension (might not be optimal, but seems same as h5repack)
        chunk[0] = dim0;
        chunk[1] = dim1;
        chunk[2] = dim2;
        // create dataset creation property list
        dcpl = H5Pcreate( H5P_DATASET_CREATE );
        // add shuffle filter and gzip compression filter
        // note that the order of filter is significant: first shuffle!
        // order of filters applied correspond to order in which they are invoked when writin gdata
        status = H5Pset_shuffle( dcpl );
        if (status < 0) printf( "ERROR: could not add shuffle filter\n" );
        status = H5Pset_deflate( dcpl, 9 );
        if (status < 0) printf( "ERROR: could not add gzip filter\n" );
        status = H5Pset_chunk(dcpl,             // dataset creation property list identifier
                              3,                // number of dimensions of each chunk
                              chunk );          // array defining size, in dataset elements, of each chunk
        if (status < 0) printf( "ERROR: could not set chunk size\n" );
        // create the dataset
        dataset_id = H5Dcreate( file_id,        // file identifier (or group identifier)
                                dataset,        // name of dataset (relative to group specified, if speficied)
                                H5T_NATIVE_DOUBLE, // datatype to use when creating dataset
                                dataspace_id,   // dataspace identifier
                                H5P_DEFAULT,    // link creation property list (was dataset creating property list <=v1.6)
                                dcpl,           // dataset creation property list (added in HDF5v1.8)
                                H5P_DEFAULT);   // dataset access property list (added in HDF5v1.8)
    } else {
        // create the dataset
        dataset_id = H5Dcreate( file_id,        // file identifier (or group identifier)
                                dataset,        // name of dataset (relative to group specified, if speficied)
                                H5T_NATIVE_DOUBLE, // datatype to use when creating dataset
                                dataspace_id,   // dataspace identifier
                                H5P_DEFAULT,    // link creation property list (was dataset creating property list <=v1.6)
                                H5P_DEFAULT,    // dataset creation property list (added in HDF5v1.8)
                                H5P_DEFAULT);   // dataset access property list (added in HDF5v1.8)
    }

    // write the dataset
    status = H5Dwrite( dataset_id,          // dataset identifier
                       H5T_NATIVE_DOUBLE,    // informs hdf about format of data in memory of computer
                       H5S_ALL,             // identifier of memory dataspace
                       H5S_ALL,             // file space identifier
                       H5P_DEFAULT,         // data transfer property list
//                       array_2D[0]);        // pointer to data array
                       array_3D);        // pointer to data array
    if (status < 0) 
        printf( "ERROR: could not write dataset '%s' to file '%s'\n", dataset, filename );

    // terminate access and free ressources/identifiers
    // dataset creation property list
    if (filter_avail) {
        status = H5Pclose(dcpl);
        if (status < 0)
            printf( "ERROR: could not close filter\n" );
    }
    // dataset
    status = H5Dclose(dataset_id);
    if (status < 0) printf( "ERROR: could not close dataset '%s'\n", dataset );
    // data space
    status = H5Sclose(dataspace_id);
    if (status < 0) printf( "ERROR: could not close dataspace for dataset '%s'\n", dataset );
    // file 
    status = H5Fclose(file_id);
    if (status < 0) printf( "ERROR: could not close file '%s'\n", filename );
    
    return EXIT_SUCCESS;
}//#}}}
#endif


#ifdef HDF5
int writeConfig2HDF( gridConfiguration *gridCfg, beamConfiguration *beamCfg, char filename[] ) {
    //#{{{

    long        data2write_long[1];
    double      data2write_double[1];

    // hdf related variables
    hid_t       file_id, dataset_id, dataspace_id;      // object identifiers
    hsize_t     dims[1];                                // size used for dimensions
    herr_t      status;                                 // function return value

    // note that shuffle and compression filter is not applied here, as only single values are saved

    // required for check if hdf5-file already exists
    struct      stat st;

    // check if specified hdf5 file already exists
    // if not, create new one; if yes, open and add dataset to it
    if ( stat( filename, &st )==0 ) {
        // open file for read + write access
        file_id = H5Fopen( filename,            // filename
                           H5F_ACC_RDWR,        // allow read & write access (_RDONLY for read only)
                           H5P_DEFAULT);        // file access property list (default one)

        // hdf5 version 1.8.0 introduced H5Lexists to check if link (to group or dataset) exists in hdf5-file
#if H5_VERS_MAJOR>=1 && H5_VERS_MINOR>=8
        if ( H5_VERS_MINOR >= 10 ) {
            printf( "WARNING: hdf5 version 1.10 (or larger is used)\n" );
            printf( "         behavior of H5Lexists was slightly changed in this version\n" );
            printf( "         for details, see https://support.hdfgroup.org/HDF5/doc/RM/RM_H5L.html#Link-Exists\n" );
        }
        if ( H5Lexists( file_id,                // file or group identifier
                        "/config",              // name of link (to group or dataset) to check
                        H5P_DEFAULT )           // link access property list identifiert
                > 0 ) {                         // NOTE: for version 1.8.10, this might be slightly different
            printf( "ERROR: dataset named '/config' already exists in file '%s'\n", filename );
            printf( "       dataset will NOT be saved (no overwrite by default)\n" );
            status = H5Fclose(file_id);
            if (status < 0) printf("ERROR: could not close file '%s'\n", filename);
            return EXIT_FAILURE;
        }
#endif
    } else {
        // create a new file using default properties.
        file_id = H5Fcreate( filename,          // filename
                             H5F_ACC_TRUNC,     // how file should be created (removes existing file)
                             H5P_DEFAULT,       // file creating property list
                             H5P_DEFAULT);      // file access property list
    }

    // create new group for config parameters
    H5Gcreate( file_id, 
               "/config",
               H5P_DEFAULT,
               H5P_DEFAULT,
               H5P_DEFAULT);

    // create simple data space for the dataset
    // (simple = regular N-dimensional array, i.e. data on regular rectangular grid)
    // (complex = e.g.: irregular grids resulting from dynamic refinement of mesh)
    dims[0] = 1;
    dataspace_id = H5Screate_simple( 1,     // number of dimensions of dataspace 
                                     dims,  // size of array in each dimension
                                     NULL); // allow stretching of data space (NULL=no)

    // period
    // create the dataset
    dataset_id = H5Dcreate( file_id,        // file identifier (or group identifier)
                            "/config/period",  // name of dataset (relative to group specified, if speficied)
                            H5T_NATIVE_LONG, // datatype to use when creating dataset
                            dataspace_id,   // dataspace identifier
                            H5P_DEFAULT,    // link creation property list (was dataset creating property list <=v1.6)
                            H5P_DEFAULT,    // dataset creation property list (added in HDF5v1.8)
                            H5P_DEFAULT);   // dataset access property list (added in HDF5v1.8)

    // write the dataset
    data2write_long[0]  = (long)gridCfg->period;
    status = H5Dwrite( dataset_id,          // dataset identifier
                       H5T_NATIVE_LONG,   // informs hdf about format of data in memory of computer
                       H5S_ALL,             // identifier of memory dataspace
                       H5S_ALL,             // file space identifier
                       H5P_DEFAULT,         // data transfer property list
                       data2write_long);     // pointer to data array
    if (status < 0) printf( "ERROR: could not write dataset '/config/period' into file '%s'\n", filename);
    // terminate access and free ressources/identifiers for dataset
    status = H5Dclose(dataset_id);
    if (status < 0) printf("ERROR: could not close dataset '/config/period'\n");

    // d_absorb
    dataset_id = H5Dcreate( file_id, "/config/d_absorb", H5T_NATIVE_LONG,
                            dataspace_id, 
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    data2write_long[0]  = (long)gridCfg->d_absorb;
    status = H5Dwrite( dataset_id, H5T_NATIVE_LONG,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT,
                       data2write_long); 
    if (status < 0) printf( "ERROR: could not write dataset '/config/d_absorb' into file '%s'\n", filename);
    status = H5Dclose(dataset_id);
    if (status < 0) printf("ERROR: could not close dataset '/config/d_absorb'\n");

    // N_x
    dataset_id = H5Dcreate( file_id, "/config/N_x", H5T_NATIVE_LONG,
                            dataspace_id, 
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    data2write_long[0]  = (long)gridCfg->Nx;
    status = H5Dwrite( dataset_id, H5T_NATIVE_LONG,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT,
                       data2write_long); 
    if (status < 0) printf( "ERROR: could not write dataset '/config/N_x' into file '%s'\n", filename);
    status = H5Dclose(dataset_id);
    if (status < 0) printf("ERROR: could not close dataset '/config/N_x'\n");

    // N_y
    dataset_id = H5Dcreate( file_id, "/config/N_y", H5T_NATIVE_LONG,
                            dataspace_id, 
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    data2write_long[0]  = (long)gridCfg->Ny;
    status = H5Dwrite( dataset_id, H5T_NATIVE_LONG,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT,
                       data2write_long); 
    if (status < 0) printf( "ERROR: could not write dataset '/config/N_y' into file '%s'\n", filename);
    status = H5Dclose(dataset_id);
    if (status < 0) printf("ERROR: could not close dataset '/config/N_y'\n");

    // N_z
    dataset_id = H5Dcreate( file_id, "/config/N_z", H5T_NATIVE_LONG,
                            dataspace_id, 
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    data2write_long[0]  = (long)gridCfg->Nz;
    status = H5Dwrite( dataset_id, H5T_NATIVE_LONG,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT,
                       data2write_long); 
    if (status < 0) printf( "ERROR: could not write dataset '/config/N_z' into file '%s'\n", filename);
    status = H5Dclose(dataset_id);
    if (status < 0) printf("ERROR: could not close dataset '/config/N_z'\n");

    // ant_x
    dataset_id = H5Dcreate( file_id, "/config/ant_x", H5T_NATIVE_LONG,
                            dataspace_id, 
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    data2write_long[0]  = (long)beamCfg->ant_x;
    status = H5Dwrite( dataset_id, H5T_NATIVE_LONG,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT,
                       data2write_long); 
    if (status < 0) printf( "ERROR: could not write dataset '/config/ant_x' into file '%s'\n", filename);
    status = H5Dclose(dataset_id);
    if (status < 0) printf("ERROR: could not close dataset '/config/ant_x'\n");

    // ant_y
    dataset_id = H5Dcreate( file_id, "/config/ant_y", H5T_NATIVE_LONG,
                            dataspace_id, 
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    data2write_long[0]  = (long)beamCfg->ant_y;
    status = H5Dwrite( dataset_id, H5T_NATIVE_LONG,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT,
                       data2write_long); 
    if (status < 0) printf( "ERROR: could not write dataset '/config/ant_y' into file '%s'\n", filename);
    status = H5Dclose(dataset_id);
    if (status < 0) printf("ERROR: could not close dataset '/config/ant_y'\n");

    // ant_z
    dataset_id = H5Dcreate( file_id, "/config/ant_z", H5T_NATIVE_LONG,
                            dataspace_id, 
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    data2write_long[0]  = (long)beamCfg->ant_z;
    status = H5Dwrite( dataset_id, H5T_NATIVE_LONG,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT,
                       data2write_long); 
    if (status < 0) printf( "ERROR: could not write dataset '/config/ant_z' into file '%s'\n", filename);
    status = H5Dclose(dataset_id);
    if (status < 0) printf("ERROR: could not close dataset '/config/ant_z'\n");

    // ant_w0x
    dataset_id = H5Dcreate( file_id, "/config/ant_w0x", H5T_NATIVE_DOUBLE,
                            dataspace_id, 
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    data2write_double[0]  = beamCfg->ant_w0x;
    status = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT,
                       data2write_double); 
    if (status < 0) printf( "ERROR: could not write dataset '/config/ant_w0x' into file '%s'\n", filename);
    status = H5Dclose(dataset_id);
    if (status < 0) printf("ERROR: could not close dataset '/config/ant_w0x'\n");

    // ant_w0y
    dataset_id = H5Dcreate( file_id, "/config/ant_w0y", H5T_NATIVE_DOUBLE,
                            dataspace_id, 
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    data2write_double[0]  = beamCfg->ant_w0y;
    status = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT,
                       data2write_double); 
    if (status < 0) printf( "ERROR: could not write dataset '/config/ant_w0y' into file '%s'\n", filename);
    status = H5Dclose(dataset_id);
    if (status < 0) printf("ERROR: could not close dataset '/config/ant_w0y'\n");

    // z2waist
    dataset_id = H5Dcreate( file_id, "/config/z2waist", H5T_NATIVE_DOUBLE,
                            dataspace_id, 
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    data2write_double[0]  = beamCfg->z2waist;
    status = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE,
                       H5S_ALL, H5S_ALL, H5P_DEFAULT,
                       data2write_double); 
    if (status < 0) printf( "ERROR: could not write dataset '/config/z2waist' into file '%s'\n", filename);
    status = H5Dclose(dataset_id);
    if (status < 0) printf("ERROR: could not close dataset '/config/z2waist'\n");

    // terminate access and free ressources/identifiers
    // data space
    status = H5Sclose(dataspace_id);
    if (status < 0) printf("ERROR: could not close dataspace'\n");
    // file 
    status = H5Fclose(file_id);
    if (status < 0) printf( "ERROR: could not close file '%s'\n", filename);
    
    return EXIT_SUCCESS;
}//#}}}
#endif


#ifdef HDF5
int readMyHDF( int dim0, int dim1, int dim2, char filename[], char dataset[], double array_3D[dim0][dim1][dim2]) {
    //#{{{

    // hdf handles
    hid_t           file_id, dset_id;
    herr_t          status;
    //hsize_t         dims[3] = { dim0, dim1, dim2};

    // open file using default properties
    file_id = H5Fopen( filename, H5F_ACC_RDONLY, H5P_DEFAULT);

    // open dataset using default properties
    dset_id = H5Dopen( file_id, dataset, H5P_DEFAULT);

    // Read the data using the default properties.
    status = H5Dread( dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      array_3D );
    if (status < 0)
        printf( "ERROR: could not read dataset '%s' from file '%s'\n", dataset, filename );

    // close the dataset
    status = H5Dclose( dset_id);
    if (status < 0)
        printf( "ERROR: could not close dataset '%s' from file '%s'\n", dataset, filename );

    // close the file
    status = H5Fclose( file_id);
    if (status < 0)
        printf( "ERROR: could not close file '%s'\n", filename );

    return EXIT_SUCCESS;
}//#}}}
#endif

