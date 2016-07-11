//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File:    rand32.h							Author:     Bob Harris
//
//----------
//
// Interface definitions for C support module
//
//----------

//----------
//
// prototypes for routines in this module
//
//----------

void			srand32  (unsigned long seed);
void			ssrand32 (char* seed);
unsigned long	rand32   (void);
long			rrand32  (long lo, long hi);
unsigned long	urand32  (unsigned long lo, unsigned long hi);
float			frand32  (void);
int				crand32  (float p);
int				drand32  (int piValues, unsigned long* piSum);

