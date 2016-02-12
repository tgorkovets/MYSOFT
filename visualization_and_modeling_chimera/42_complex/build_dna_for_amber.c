#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nabcode.h"
extern char NAB_rsbuf[];
static int mytaskid, numtasks;





static MOLECULE_T *m,  *m2;

static INT_T seqlen;

static STRING_T *seq = NULL;

static ATOM_T *a;


int main( argc, argv )
	int	argc;
	char	*argv[];
{
	nabout = stdout; /*default*/

	mytaskid=0; numtasks=1;
static STRING_T *__st0001__ = NULL;
static STRING_T *__st0002__ = NULL;
NAB_strcpy(  &seq, argv[2 - 1] );

seqlen = length( seq );

m2 = fd_helix( STEMP( __st0001__, "lbdna" ),  &seq, STEMP( __st0002__, "dna" ) );









putpdb( argv[3 - 1], m2, NULL );






	exit( 0 );
}
