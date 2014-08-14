/*  ---------------------------------------------------------------  */
/*           The HMM-Based Speech Synthesis System (HTS)             */
/*                       HTS Working Group                           */
/*                                                                   */
/*                  Department of Computer Science                   */
/*                  Nagoya Institute of Technology                   */
/*                               and                                 */
/*   Interdisciplinary Graduate School of Science and Engineering    */
/*                  Tokyo Institute of Technology                    */
/*                                                                   */
/*                     Copyright (c) 2001-2007                       */
/*                       All Rights Reserved.                        */
/*                                                                   */
/*  Permission is hereby granted, free of charge, to use and         */
/*  distribute this software and its documentation without           */
/*  restriction, including without limitation the rights to use,     */
/*  copy, modify, merge, publish, distribute, sublicense, and/or     */
/*  sell copies of this work, and to permit persons to whom this     */
/*  work is furnished to do so, subject to the following conditions: */
/*                                                                   */
/*    1. The source code must retain the above copyright notice,     */
/*       this list of conditions and the following disclaimer.       */
/*                                                                   */
/*    2. Any modifications to the source code must be clearly        */
/*       marked as such.                                             */
/*                                                                   */
/*    3. Redistributions in binary form must reproduce the above     */
/*       copyright notice, this list of conditions and the           */
/*       following disclaimer in the documentation and/or other      */
/*       materials provided with the distribution.  Otherwise, one   */
/*       must contact the HTS working group.                         */
/*                                                                   */
/*  NAGOYA INSTITUTE OF TECHNOLOGY, TOKYO INSTITUTE OF TECHNOLOGY,   */
/*  HTS WORKING GROUP, AND THE CONTRIBUTORS TO THIS WORK DISCLAIM    */
/*  ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL       */
/*  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT   */
/*  SHALL NAGOYA INSTITUTE OF TECHNOLOGY, TOKYO INSTITUTE OF         */
/*  TECHNOLOGY, HTS WORKING GROUP, NOR THE CONTRIBUTORS BE LIABLE    */
/*  FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY        */
/*  DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,  */
/*  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTUOUS   */
/*  ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR          */
/*  PERFORMANCE OF THIS SOFTWARE.                                    */
/*                                                                   */
/*  ---------------------------------------------------------------  */
/*                 File: lbfgs.h: L-BFGS routine                     */
/*  ---------------------------------------------------------------  */

/* 

 This L-BFGS FORTRAN source code is originally distributed in

   http://www.ece.northwestern.edu/~nocedal/lbfgs.html

 We thank Prof. Jorge Nocedal of Northwestern University for permission to 
 redistribute this LBFGS code in the HTS releases.


 L-BFGS: Software for Large-scale Unconstrained Optimization
   
   L-BFGS is a limited-memory quasi-Newton code for unconstrained optimization. 
   The code has been developed at the Optimization Technology Center, a joint 
   venture of Argonne National Laboratory and Northwestern University.

 Condition for Use: 

   This software is freely available for educational or commercial purposes. We 
   expect that all publications describing work using this software quote at 
   least one of the references given below.

 References

   * J. Nocedal, "Updating Quasi-Newton Matrices with Limited Storage," 
     Mathematics of Computation 35, pp. 773-782, 1980.

   * D.C. Liu and J. Nocedal, "On the Limited Memory BFGS Method for Large Scale 
     Optimization," Mathematical Programming B, 45, 3, pp. 503-528, 1989.

*/


/* !HVER!lbfgs:   2.1 [NIT 31/10/07] */

#ifndef _LBFGS_H_
#define _LBFGS_H_

#ifdef __cplusplus
extern "C" {
#endif 

void lbfgs_(int* n, int* m, double* x, double* f, double* g,
            int* diagco, double* diag, int* iprint, double* eps,
            double* xtol, double* w, int* iflag);

#ifdef __cplusplus
}
#endif

#endif  /* _LBFGS_H_ */

/* ------------------------- End of lbfgs.h -------------------------- */

