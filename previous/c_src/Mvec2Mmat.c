#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "compearth.h"

/*!
 * @brief Converts between two representations of a set of moment tensors
 *        (1) 6 x n and
 *        (2) 3 x 3 x n.
 *
 * @param[in] n       Number of moment tensors.
 * @param[in] Min     If itype is 1 then this is the moment tensor packed
 *                    [M11 M22 M33 M12 M13 M23].  In this case Min is
 *                    an array of dimension [6 x n] with leading dimension 6. \n
 *                    Otherwise this is the moment tensor packed 
 *                    [M11 M12 M13 M21 M22 M23 M31 M32 M33] and Min is
 *                    an array of dimension [9 x n] with elading dimension 9.
 * @param[in] itype   If 1 then convert from 6 x n to 3 x 3 x n. \n
 *                    Otherwise convert from 3 x 3 x n to 6 x n.
 *
 * @param[out] Mout   If itype is 1 then this is the moment tensor packed
 *                    [M11 M12 M13 M21 M22 M23 M31 M32 M33].  In this case
 *                    Mout is an array of dimension [9 x n] with leading 
 *                    dimension 9. \n
 *                    Otherwise this is the moment tensor packed
 *                    [M11 M22 M33 M12 M13 M23] and Mout is an array of
 *                    dimension [6 x n] with leading dimension 6.
 *
 * @author Carl Tape and Converted to C by Ben Baker
 *
 * @date 2016
 *
 * @copyright MIT
 *
 */
void compearth_Mvec2Mmat(const int n,
                         const double *__restrict__ Min,
                         const int itype,
                         double *__restrict__ Mout)
{
    double Mrr, Mtt, Mpp, Mrt, Mrp, Mtp;
    int i;
    // 6 --> 3 x 3
    if (itype == 1)
    {
        for (i=0; i<n; i++)
        {
            Mrr = Min[6*i+0];
            Mtt = Min[6*i+1];
            Mpp = Min[6*i+2];
            Mrt = Min[6*i+3];
            Mrp = Min[6*i+4];
            Mtp = Min[6*i+5];
            Mout[9*i+0] = Mrr;
            Mout[9*i+1] = Mrt;
            Mout[9*i+2] = Mrp;
            Mout[9*i+3] = Mrt;
            Mout[9*i+4] = Mtt;
            Mout[9*i+5] = Mtp;
            Mout[9*i+6] = Mrp;
            Mout[9*i+7] = Mtp;
            Mout[9*i+8] = Mpp;
        }
    }
    // 3 x 3 --> 6
    else
    {
        for (i=0; i<n; i++)
        {
            Mrr = Min[9*i+0];
            Mrt = Min[9*i+1];
            Mrp = Min[9*i+2];
            Mtt = Min[9*i+4];
            Mtp = Min[9*i+5];
            Mpp = Min[9*i+8];
            Mout[6*i+0] = Mrr; //Min[9*i+0]; // (1,1)
            Mout[6*i+1] = Mtt; //Min[9*i+4]; // (2,2)
            Mout[6*i+2] = Mpp; //Min[9*i+8]; // (3,3)
            Mout[6*i+3] = Mrt; //Min[9*i+3]; // (1,2)
            Mout[6*i+4] = Mrp; //Min[9*i+6]; // (1,3)
            Mout[6*i+5] = Mtp; //Min[9*i+7]; // (2,3)
        }
    }
    return;
}
