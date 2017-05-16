/*
 * Library:  ransampl (random number sampling)
 *
 * File:     sampling1.c
 *
 * Contents: Draw samples from a given discrete probability distribution;
 *           specifically: draw representative inhabitants of the nine states
 *
 * Note:     Any modification of this example should be copied to
 *           the manual page source ransampl.pod and to the wiki.
 *
 * Author:   Joachim Wuttke 2013
 * 
 * Licence:  see ../COPYING (FreeBSD)
 * 
 * Homepage: apps.jcns.fz-juelich.de/ransampl
 */
 
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "ransampl.h"

#define RANSAMPL_SHALL_REVEAL_ITS_INTERNALS 1

int main()
{
    const int M=1000000;
    int i, m;

    // Discrete probability distribution example:
    const int n = 9;
    // states of Austria
    const char* names[] = {
        "Wien", "Niederoesterreich", "Oberoesterreich", "Tirol",
        "Kaernten", "Salzburg", "Vorarlberg", "Burgenland", "Steiermark" };
    // inhabitants in millions as of 2011 according to www.statistik.at
    double p[] = { 1.721573, 1.614661, 1.415020, .711161,
                   .558056, .532713, .370833, .285377, .1211506 };

    // Initialize random number generator:
    gsl_rng_env_setup();
    gsl_rng* rng = gsl_rng_alloc( gsl_rng_default );

    // Allocate workspace, and precompute tables:
    printf( "Precomputing tables ...\n" );
    ransampl_ws* ws = ransampl_alloc( n );
    ransampl_set( ws, p );

#ifdef RANSAMPL_SHALL_REVEAL_ITS_INTERNALS
    // Inspect tables:
    printf( "  %-3s  %-3s  %-9s\n", "i", "alias", "prob" );
    for ( int i=0; i<n; ++i )
        printf( "  %3i  %3i  %9.7f\n", i, ws->alias[i], ws->prob[i] );
#endif
    
    // Draw M random samples; accumulate statistic in histogram 'cumul':
    printf( "Drawing %i samples ...\n", M );
    double cumul[n];
    for ( i=0; i<n; ++i )
        cumul[i] = 0;
    for ( m=0; m<M; ++m ) {
        i = ransampl_draw( ws, gsl_rng_uniform(rng), gsl_rng_uniform(rng) );
        cumul[i] += 1;
    }

    // Print given probability and obtained frequency:
    printf( "Result (input->output):\n");
    double sum = 0;
    for ( int i=0; i<n; ++i )
        sum += p[i];
    printf( "  %-18s  %-9s  %-9s  %-9s\n", "state", "N (Mio.)", "rel", "sim" );
    for ( int i=0; i<n; ++i )
        printf( "  %-18s  %9.7f  %9.7f  %9.7f\n",
                names[i], p[i], p[i]/sum, ((double)cumul[i])/M );

    // Free workspace and terminate:
    ransampl_free( ws );
    return 0;
}
