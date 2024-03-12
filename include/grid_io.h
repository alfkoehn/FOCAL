#ifndef GRID_IO_H
#define GRID_IO_H

int writeTimetraces2ascii( int dim0, int dim1, int t_end, double period, 
                           char filename[], double timetraces[dim0][dim1] );

#endif  // GRID_IO_H
