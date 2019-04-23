#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

static const double a_rk[] = {0.0, 3.0/4.0, 1.0/3.0};  /*!< Coefficients (left) for 
                                                            SSP-RK3 scheme by Shu-Osher */
static const double b_rk[] = {1.0, 1.0/4.0, 2.0/3.0};  /*!< Coefficients (right) for 
                                                            SSP-RK3 scheme by Shu-Osher */
                                                            
static const double rk4_rks[] = {1.0,1.0/2.0,1.0/2.0,1.0};  /*!< Coefficients for 
                                                            usual RK4 */                                                          

static const double KKK = 1.0/3.0; /*!< Used in the MUSCL scheme */
               
enum CellType {median,     /*!< Median dual cells to be considered */ 
               voronoi     /*!< Voronoi dual cells to be considered */
               }; 
enum Dimension {two,       /*!< 2D problem setup */ 
                axi        /*!< 3D problem setup, but with axial symmetry */
               };

#endif
