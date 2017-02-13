/*  md - Molecular Dynamics Applied to ionic solids.
*   Copyright (C) 2017 Nils Harmening, Marco Manni,
*   Darian Steven Viezzer, Steffanie Kieninger, Henrik Narvaez
*
*   This program is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   This program is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

double _fast_neighbourlist(double **R, int N, double box_length, double r_cutoff)
{
    // assume same size in all N dimension
    int n_cells = (int) box_length / r_cutoff;
    // divide simulation box into small cells of equal size r_c >= r_cutoff
    double r_c = box_length / n_cells;
    
    
    // define head and list 
    int head[n_cells*n_cells*n_cells];
    int cllist[N];
    int x[3];
    for (int c=0; c<n_cells*n_cells*n_cells; c++) head[c] = -1;
    // Loop over all particles to construct header (head) and linked lists (cllist)
    for (int i=0; i<N; i++) {
        /* Vector cell index to which this atom belongs */
        for (int dim=0; dim<3; dim++) x[dim] = R[i][dim]/r_c;
        int cell_index = x[0]*n_cells*n_cells+x[1]*n_cells+x[2];
        // Link to the previous occupant (or EMPTY if it's the 1st)
        cllist[i] = head[cell_index];
        // The last one goes to the header
        head[cell_index] = i;
    }



    // Scan inner cells
    int y[3];
    double rshift[3];
    double dist[N][N];
    for (x[0]=0; x[0]<n_cells; (x[0])++) {
      for (x[1]=0; x[1]<n_cells; (x[1])++) {
        for (x[2]=0; x[2]<n_cells; (x[2])++) {
          int cell_index = x[0]*n_cells*n_cells+x[1]*n_cells+x[2];
          // Scan the neighbor cells (including itself)
          for (y[0]=x[0]-1; y[0]<=x[0]+1; (y[0])++) {
            for (y[1]=x[1]-1; y[1]<=x[1]+1; (y[1])++) {
              for (y[2]=x[2]-1; y[2]<=x[2]+1; (y[2])++) {
                // Periodic boundary condition by shifting coordinates
                for (int dim=0; dim<3; dim++) {
                  if (x[dim] < 0)
                      rshift[dim] = -box_length;
                  else if (y[dim]>=n_cells)
                      rshift[dim] = box_length;
                  else
                      rshift[dim] = 0.0;
                }
                // Calculate the scalar cell index of the neighbor cell */
                int cell_indy = ((y[0]+n_cells)%n_cells)*n_cells
                  +((y[1]+n_cells)%n_cells)*n_cells
                  +((y[2]+n_cells)%n_cells);
                // Scan atom i in cell cell_index
                int i = head[cell_index];
                while (i != -1) {
                    // Scan atom j in neighboring cell
                    int j = head[cell_indy];
                    while (j != -1) {
                        if (i<j)  {
                          double rij=0.0;
                          // Image corrected relative pair position
                          for (int dim=0; dim<3; dim++) rij += R[i][dim]-(R[j][dim]+rshift[dim]);
                          if (rij<=r_cutoff) {
                            dist[i][j]=rij;
                            dist[j][i]=rij;
                          }
                          else {
                            dist[i][j]=0.0;
                            dist[j][i]=0.0;
                          }
                        }
                        j = cllist[j];
                    }
                    i = cllist[i];
                }
              }
            }
          }
        }
      }
    }
    return dist[N][N];
}
