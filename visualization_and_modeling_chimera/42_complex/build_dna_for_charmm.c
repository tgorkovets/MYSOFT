#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nabcode.h"
extern char NAB_rsbuf[];
static int mytaskid, numtasks;


































INT_T myput_tmp( REAL_T *x, REAL_T *y, REAL_T *z, INT_T *count, INT_T *rescount, STRING_T * *atomname, STRING_T * *resname, STRING_T * *chname, FILE_T * *outfile, INT_T *nres )
{
INT_T rn;

rn =  *rescount;
if( rn > (  *nres / 2 ) )
rn = rn -  *nres / 2;






if( EQ(  *atomname, "O1'" ) )
NAB_strcpy(  & *atomname, "O4'" );

if( EQ(  *atomname, "O1P" ) )
NAB_strcpy(  & *atomname, "OP1" );

if( EQ(  *atomname, "O2P" ) )
NAB_strcpy(  & *atomname, "OP2" );






fprintf(  *outfile, "ATOM  %5d %-4s %3s %s%4d    %8.3lf%8.3lf%8.3lf\n",  *count,  *atomname,  *resname,  *chname, rn,  *x,  *y,  *z );

return( 1 );

}




MOLECULE_T *my_fd_helix( STRING_T * *helix_type, STRING_T * *seq, STRING_T * *acid_type )
{





MOLECULE_T *m;





HASH_T *ade_r = NULL;

HASH_T *ade_phi = NULL;

HASH_T *ade_zz = NULL;


HASH_T *gua_r = NULL;

HASH_T *gua_phi = NULL;

HASH_T *gua_zz = NULL;


HASH_T *thy_r = NULL;

HASH_T *thy_phi = NULL;

HASH_T *thy_zz = NULL;


HASH_T *ura_r = NULL;

HASH_T *ura_phi = NULL;

HASH_T *ura_zz = NULL;


HASH_T *cyt_r = NULL;

HASH_T *cyt_phi = NULL;

HASH_T *cyt_zz = NULL;




REAL_T temp_r, temp_phi, temp_zz;




REAL_T x, y, z, yyr, xrad;




REAL_T current_height, current_rotation;

REAL_T height_increment, rotation_increment;


HASH_T *hxht = NULL;

HASH_T *hxrep = NULL;


STRING_T *tempname = NULL,  *cseq = NULL,  *fullseq = NULL,  *buffer = NULL,  *restype = NULL;

STRING_T *temp = NULL,  *amberhome = NULL;

STRING_T *resout = NULL,  *strname = NULL;


INT_T nres, nresh, i, hxmul, count, chain, begin, end;


FILE_T *infile,  *outfile;





STRING_T *__st0001__ = NULL;
CURHASH_T __cht0001__;
CURHASH_T __cht0002__;
CURHASH_T __cht0003__;
CURHASH_T __cht0004__;
CURHASH_T __cht0005__;
HRF(  &hxht, "arna", 3 ) = 2.810000E+00;HRF(  &hxht, "aprna", 3 ) = 3.000000E+00;HRF(  &hxht, "lbdna", 3 ) = 3.380000E+00;
HRF(  &hxht, "abdna", 3 ) = 3.380000E+00;HRF(  &hxht, "sbdna", 3 ) =  - 3.380000E+00;HRF(  &hxht, "adna", 3 ) = 2.560000E+00;

HRF(  &hxrep, "arna", 3 ) = 3.270000E+01;HRF(  &hxrep, "aprna", 3 ) = 3.000000E+01;HRF(  &hxrep, "lbdna", 3 ) = 3.600000E+01;
HRF(  &hxrep, "abdna", 3 ) = 3.600000E+01;HRF(  &hxrep, "sbdna", 3 ) = 3.600000E+01;HRF(  &hxrep, "adna", 3 ) = 3.270000E+01;

NAB_strcpy(  &temp, wc_complement( seq, STEMP( __st0001__, NAB_strcat(  *acid_type, "amber94.rlb" ) ), acid_type ) );
NAB_strcpy(  &cseq, "" );
for( i = length( temp );i >= 1;i --  ){
NAB_strcpy(  &cseq, NAB_strcat( cseq, substr( temp, i, 1 ) ) );
}
NAB_strcpy(  &fullseq, NAB_strcat(  *seq, cseq ) );

nresh = length(  *seq );

nres = length( fullseq );

if(  !( NAB_strcpy(  &amberhome, getenv( "AMBERHOME" ) ) ) ){
fprintf( stderr, "AMBERHOME not defined.\n" );
exit( 1 );
}

NAB_strcpy(  &temp, NAB_strcat( amberhome, NAB_strcat( "/dat/fd_data/", NAB_strcat(  *helix_type, ".dat" ) ) ) );


infile = fopen( temp, "r" );
if( infile == NULL ){
fprintf( stderr, "Unable to open data file %s; exiting\n", temp );
exit( 1 );
}

outfile = fopen( "nab_tmp.pdb", "w" );




while( NAB_strcpy(  &buffer, NAB_getline( infile ) ) ){

sscanf( buffer, "%s %lf %lf %lf %s", NAB_readstring(  &tempname ),  &temp_r,  &temp_phi,  &temp_zz, NAB_readstring(  &restype ) );

if( ( EQ( restype, "A" ) ) || ( EQ( restype, "a" ) ) ){
HRF(  &ade_r, tempname, 3 ) = temp_r;
HRF(  &ade_phi, tempname, 3 ) = temp_phi;
HRF(  &ade_zz, tempname, 3 ) = temp_zz;
}
else if( ( EQ( restype, "G" ) ) || ( EQ( restype, "g" ) ) ){
HRF(  &gua_r, tempname, 3 ) = temp_r;
HRF(  &gua_phi, tempname, 3 ) = temp_phi;
HRF(  &gua_zz, tempname, 3 ) = temp_zz;
}
else if( ( EQ( restype, "T" ) ) || ( EQ( restype, "t" ) ) ){
HRF(  &thy_r, tempname, 3 ) = temp_r;
HRF(  &thy_phi, tempname, 3 ) = temp_phi;
HRF(  &thy_zz, tempname, 3 ) = temp_zz;
}
else if( ( EQ( restype, "U" ) ) || ( EQ( restype, "u" ) ) ){
HRF(  &ura_r, tempname, 3 ) = temp_r;
HRF(  &ura_phi, tempname, 3 ) = temp_phi;
HRF(  &ura_zz, tempname, 3 ) = temp_zz;
}
else if( ( EQ( restype, "C" ) ) || ( EQ( restype, "c" ) ) ){
HRF(  &cyt_r, tempname, 3 ) = temp_r;
HRF(  &cyt_phi, tempname, 3 ) = temp_phi;
HRF(  &cyt_zz, tempname, 3 ) = temp_zz;
}

}

height_increment = HRF(  &hxht,  *helix_type, 3 );
rotation_increment = HRF(  &hxrep,  *helix_type, 3 );
current_height = 0;
current_rotation = 0;
count = 0;






for( chain = 1;chain <= 2;chain ++  ){

if( chain == 1 ){
begin = 1;
end = nresh;
hxmul =  - 1;
NAB_strcpy(  &strname, "I" );
}

else if( chain == 2 ){
begin = nresh + 1;
end = nres;
hxmul = 1;
NAB_strcpy(  &strname, "J" );
}
for( i = begin;i <= end;i ++  ){
NAB_strcpy(  &restype, substr( fullseq, i, 1 ) );
if( ( EQ( restype, "A" ) ) || ( EQ( restype, "a" ) ) ){
NAB_strcpy(  &resout, "DA" );if( EQ(  *acid_type, "rna" ) )NAB_strcpy(  &resout, "A" );
for( NAB_hfirst( ade_r,  &__cht0001__ );NAB_strcpy(  &tempname, NAB_hnext( ade_r,  &__cht0001__ ) ); ){
count ++ ;
yyr = ( hxmul * ( HRF(  &ade_phi, tempname, 3 ) ) + current_rotation );
xrad = HRF(  &ade_r, tempname, 3 );
x = xrad * ( COS( yyr ) );
y = xrad * ( SIN( yyr ) );
z = hxmul * ( HRF(  &ade_zz, tempname, 3 ) ) + current_height;
myput_tmp(  &x,  &y,  &z,  &count,  &i,  &tempname,  &resout,  &strname,  &outfile,  &nres );
}
}

else if( ( EQ( restype, "G" ) ) || ( EQ( restype, "g" ) ) ){
NAB_strcpy(  &resout, "DG" );if( EQ(  *acid_type, "rna" ) )NAB_strcpy(  &resout, "G" );
for( NAB_hfirst( gua_r,  &__cht0002__ );NAB_strcpy(  &tempname, NAB_hnext( gua_r,  &__cht0002__ ) ); ){
count ++ ;
yyr = ( hxmul * ( HRF(  &gua_phi, tempname, 3 ) ) + current_rotation );
xrad = HRF(  &gua_r, tempname, 3 );
x = xrad * ( COS( yyr ) );
y = xrad * ( SIN( yyr ) );
z = hxmul * ( HRF(  &gua_zz, tempname, 3 ) ) + current_height;
myput_tmp(  &x,  &y,  &z,  &count,  &i,  &tempname,  &resout,  &strname,  &outfile,  &nres );
}

}

else if( ( EQ( restype, "T" ) ) || ( EQ( restype, "t" ) ) ){
NAB_strcpy(  &resout, "DT" );
for( NAB_hfirst( thy_r,  &__cht0003__ );NAB_strcpy(  &tempname, NAB_hnext( thy_r,  &__cht0003__ ) ); ){
count ++ ;
yyr = ( hxmul * ( HRF(  &thy_phi, tempname, 3 ) ) + current_rotation );
xrad = HRF(  &thy_r, tempname, 3 );
x = xrad * ( COS( yyr ) );
y = xrad * ( SIN( yyr ) );
z = hxmul * ( HRF(  &thy_zz, tempname, 3 ) ) + current_height;
myput_tmp(  &x,  &y,  &z,  &count,  &i,  &tempname,  &resout,  &strname,  &outfile,  &nres );
}

}

else if( ( EQ( restype, "U" ) ) || ( EQ( restype, "u" ) ) ){
NAB_strcpy(  &resout, "U" );
for( NAB_hfirst( ura_r,  &__cht0004__ );NAB_strcpy(  &tempname, NAB_hnext( ura_r,  &__cht0004__ ) ); ){
count ++ ;
yyr = ( hxmul * ( HRF(  &ura_phi, tempname, 3 ) ) + current_rotation );
xrad = HRF(  &ura_r, tempname, 3 );
x = xrad * ( COS( yyr ) );
y = xrad * ( SIN( yyr ) );
z = hxmul * ( HRF(  &ura_zz, tempname, 3 ) ) + current_height;
myput_tmp(  &x,  &y,  &z,  &count,  &i,  &tempname,  &resout,  &strname,  &outfile,  &nres );
}

}

else if( ( EQ( restype, "C" ) ) || ( EQ( restype, "c" ) ) ){
NAB_strcpy(  &resout, "DC" );if( EQ(  *acid_type, "rna" ) )NAB_strcpy(  &resout, "C" );
for( NAB_hfirst( cyt_r,  &__cht0005__ );NAB_strcpy(  &tempname, NAB_hnext( cyt_r,  &__cht0005__ ) ); ){
count ++ ;
yyr = ( hxmul * ( HRF(  &cyt_phi, tempname, 3 ) ) + current_rotation );
xrad = HRF(  &cyt_r, tempname, 3 );
x = xrad * ( COS( yyr ) );
y = xrad * ( SIN( yyr ) );
z = hxmul * ( HRF(  &cyt_zz, tempname, 3 ) ) + current_height;
myput_tmp(  &x,  &y,  &z,  &count,  &i,  &tempname,  &resout,  &strname,  &outfile,  &nres );
}

}



current_height += height_increment;
current_rotation += rotation_increment;
}


height_increment =  - height_increment;
rotation_increment =  - rotation_increment;

current_rotation += rotation_increment;
current_height += height_increment;

if( chain == 1 )
fprintf( outfile, "TER\n" );


}

fclose( infile );





fclose( outfile );


m = getpdb( "nab_tmp.pdb", NULL );
unlink( "nab_tmp.pdb" );

return( m );

}






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

m2 = my_fd_helix( STEMP( __st0001__, "lbdna" ),  &seq, STEMP( __st0002__, "dna" ) );









putpdb( argv[3 - 1], m2, NULL );






	exit( 0 );
}
