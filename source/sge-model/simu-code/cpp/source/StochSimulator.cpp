/*
__ Francois Bertaux, Imperial College
__ f.bertaux@imperial.ac.uk
__ December 2015
*/

#include "StochSimulator.hpp"

// additionally needed function
void sparmatfill(NRvector<NRsparseCol> &sparmat, MatDoub &fullmat)
{
    Int n,m,nz,nn=fullmat.nrows(),mm=fullmat.ncols();
    if (sparmat.size() != mm) throw("bad sizes");
    for (m=0;m<mm;m++)
    {
        for (nz=n=0;n<nn;n++) if (fullmat[n][m]) nz++;
        sparmat[m].resize(nn,nz);
        for (nz=n=0;n<nn;n++) if (fullmat[n][m])
        {
            sparmat[m].row_ind[nz] = n;
            sparmat[m].val[nz++] = fullmat[n][m];
        }
    }
}


// constructor
StochSimulator::StochSimulator ( ModelParameters* modelParameters , Int seed )
    : mf_ModelParameters(modelParameters), ran(seed),
    mm(mf_ModelParameters->mf_NumReacs), nn(mf_ModelParameters->mf_NumSpecies),
    a(mm,0.), outchg(mm), depend(mm), pr(mm), t(0.),asum(0.),
    dispatch(new rateptr[mm])
{
    Int i,j,k,d;
    describereactions();
    sparmatfill(outchg,outstate);
    MatDoub dep(mm,mm);
    for (i=0;i<mm;i++) for (j=0;j<mm;j++)
    {
        d = 0;
        for (k=0;k<nn;k++) d = d || (instate[k][i] && outstate[k][j]);
        dep[i][j] = d;
    }
    sparmatfill(depend,dep);
    for (i=0;i<mm;i++)
    {
        pr[i] = i;
    }
}


// preparation
void
StochSimulator::prepareForSteps ( CellState* cellState )
{
    t=0;
    Int i;
    asum = 0.;
    for (i=0;i<mm;i++)
    {
        pr[i] = i;
        a[i] = (this->*dispatch[i])(cellState->mf_SpecieCounts);
        asum += a[i];
    }
    if ( a[3]>0 ) { cout << "problem !!!!" << endl ; exit (2) ; } // because no protein degradation
}

// simulation
Doub
StochSimulator::doStep ( CellState* cellState , Doub targetTime )
{
    Int i,n,m,k=0;
    Doub tau,atarg,sum,anew;
    if (asum == 0.) {t = targetTime; return t;}
    tau = -log(ran.doub())/asum;
    if (t+tau>targetTime)
    {
        t=targetTime;
        return t;
    }
    atarg = ran.doub()*asum;
    sum = a[pr[0]];
    while (sum < atarg) sum += a[pr[++k]];
    m = pr[k];
    if ( m == 3) { cout << "problem !!!!" << endl ; exit (2) ; } // because no protein degradation
    if (k > 0) SWAP(pr[k],pr[k-1]);
    if (k == mm-1) asum = sum;
    n = outchg[m].nvals;
    for (i=0;i<n;i++)
    {
        k = outchg[m].row_ind[i];
        cellState->mf_SpecieCounts[k] += outchg[m].val[i];
    }
    n = depend[m].nvals;
    for (i=0;i<n;i++)
    {
        k = depend[m].row_ind[i];
        anew = (this->*dispatch[k])(cellState->mf_SpecieCounts);
        asum += (anew - a[k]);
        a[k] = anew;
    }
    if (t*asum < 0.1)
        for (asum=0.,i=0;i<mm;i++) asum += a[i];
    return (t += tau);
}

void
StochSimulator::simulate ( CellState* cellState , Doub duration )
{
    prepareForSteps (cellState) ;
    while ( t < duration )
    {
        doStep ( cellState , duration ) ;
    }
}

// others function, notably for init structure ( called by constructor )
void
StochSimulator::describereactions ()
{
	instate = MatDoub ( nn , mm , 0. ) ;
	outstate = MatDoub ( nn , mm , 0. ) ;
	instate[0][1] += 1.0 ;
	instate[0][2] += 1.0 ;
	instate[1][3] += 1.0 ;
	outstate[0][0] += 1.0 ;
	outstate[0][1] -= 1.0 ;
	outstate[0][2] += 1.0 ;
	outstate[1][2] += 1.0 ;
	outstate[0][2] -= 1.0 ;
	outstate[1][3] -= 1.0 ;
	dispatch[0] = &StochSimulator::rate0 ;
	dispatch[1] = &StochSimulator::rate1 ;
	dispatch[2] = &StochSimulator::rate2 ;
	dispatch[3] = &StochSimulator::rate3 ;
}

