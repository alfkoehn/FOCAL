#ifndef MACROS_BOUNDARY_H
#define MACROS_BOUNDARY_H

/*#define Nx(gridCfg)             gridCfg->Nx
#define Ny(gridCfg)             gridCfg->Ny
#define Nz(gridCfg)             gridCfg->Nz

#define Nx                      Nx(gridCfg)             
#define Ny                      Ny(gridCfg)           
#define Nz                      Nz(gridCfg)
*/

#define EBxG(PML,ii,jj,kk)       PML->EBx[( (ii) * (gridCfg->Ny - 1) + jj) * (gridCfg->Nz - 1) + kk]
#define EByG(PML,ii,jj,kk)       PML->EBy[( (ii) * (gridCfg->Ny - 1) + jj) * (gridCfg->Nz - 1) + kk]
#define EBzG(PML,ii,jj,kk)       PML->EBz[( (ii) * (gridCfg->Ny - 1) + jj) * (gridCfg->Nz - 1) + kk]
#define EBG(PML,ii,jj,kk)        PML->EB[( (ii) * (gridCfg->Ny - 1) + jj) * (gridCfg->Nz - 1) + kk]

#define EBx(ii,jj,kk)           EBxG(PML,ii,jj,kk)
#define EBy(ii,jj,kk)           EByG(PML,ii,jj,kk)
#define EBz(ii,jj,kk)           EBzG(PML,ii,jj,kk)
#define EB(ii,jj,kk)            EBG(PML,ii,jj,kk)

#define bxG(PML,nn)             PML->bx[(nn)]
#define byG(PML,nn)             PML->by[(nn)]
#define bzG(PML,nn)             PML->bz[(nn)]
#define cxG(PML,nn)             PML->cx[(nn)]
#define cyG(PML,nn)             PML->cy[(nn)]
#define czG(PML,nn)             PML->cz[(nn)]

#define bx(nn)                  bxG(PML,nn)             
#define by(nn)                  byG(PML,nn)             
#define bz(nn)                  bzG(PML,nn)             
#define cx(nn)                  cxG(PML,nn)             
#define cy(nn)                  cyG(PML,nn)             
#define cz(nn)                  czG(PML,nn)             

#endif
