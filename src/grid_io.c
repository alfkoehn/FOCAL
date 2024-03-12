
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

